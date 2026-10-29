#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace gnssins {

// common 模块是全工程的公共底座：
// 1. 定义 WGS-84、GPS/BDS 轨道计算和单位换算所需常量；
// 2. 定义三维向量、三阶矩阵、GPS 时间、经纬高等基础类型；
// 3. 提供 ECEF/LLH/NED 坐标转换和最小二乘、矩阵求逆等基础数学工具。
// 其他 GNSS、INS、Fusion 模块都只依赖这些轻量工具，不再引入第三方库。
constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kD2R = kPi / 180.0;
constexpr double kR2D = 180.0 / kPi;
constexpr double kClight = 299792458.0;
constexpr double kOmegaEarth = 7.2921151467e-5;
constexpr double kMuGps = 3.986005e14;
constexpr double kMuBds = 3.986004418e14;
constexpr double kFRel = -4.442807633e-10;
constexpr double kA = 6378137.0;
constexpr double kF = 1.0 / 298.257223563;
constexpr double kE2 = kF * (2.0 - kF);

struct Vec3 {
    // 三维向量统一采用 x/y/z 命名。
    // 在导航坐标中通常表示 N/E/D，在 ECEF 中表示 X/Y/Z，在 IMU 中表示三轴增量。
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // 按 0/1/2 下标访问 x/y/z，便于矩阵循环中统一处理三轴量。
    double &operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
    // const 版本用于只读表达式，避免复制 Vec3。
    double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }
};

// 三维向量加法，用于坐标、速度、误差状态等三轴量叠加。
inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
// 三维向量减法，用于构造坐标差、速度差和残差。
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
// 三维向量乘标量，用于噪声缩放、积分和坐标增量计算。
inline Vec3 operator*(const Vec3 &a, double s) { return {a.x * s, a.y * s, a.z * s}; }
// 三维向量除标量，主要用于归一化和平均值计算。
inline Vec3 operator/(const Vec3 &a, double s) { return {a.x / s, a.y / s, a.z / s}; }
// 点乘，常用于投影、模长和视线方向计算。
inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
// 叉乘，常用于地球自转、姿态误差和旋转向量相关计算。
inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
// 三维向量二范数。
inline double norm(const Vec3 &a) { return std::sqrt(dot(a, a)); }
// 单位化向量；零向量直接返回零，避免除零。
inline Vec3 normalized(const Vec3 &a) {
    const double n = norm(a);
    return n > 0.0 ? a / n : Vec3{};
}

struct Mat3 {
    double m[3][3] = {};
    // 构造三阶单位矩阵，常用于姿态矩阵和小角度线性化。
    static Mat3 identity() {
        Mat3 r;
        r.m[0][0] = r.m[1][1] = r.m[2][2] = 1.0;
        return r;
    }
};

// 三阶矩阵乘三维向量，用于坐标系转换和杆臂投影。
inline Vec3 matVec(const Mat3 &a, const Vec3 &v) {
    return {
        a.m[0][0] * v.x + a.m[0][1] * v.y + a.m[0][2] * v.z,
        a.m[1][0] * v.x + a.m[1][1] * v.y + a.m[1][2] * v.z,
        a.m[2][0] * v.x + a.m[2][1] * v.y + a.m[2][2] * v.z,
    };
}

// 三阶矩阵乘法，用于姿态矩阵组合和坐标转换矩阵叠加。
inline Mat3 matMul(const Mat3 &a, const Mat3 &b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) r.m[i][j] += a.m[i][k] * b.m[k][j];
        }
    }
    return r;
}

// 三阶矩阵转置。旋转矩阵转置可作为逆变换使用。
inline Mat3 transpose(const Mat3 &a) {
    Mat3 r;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) r.m[i][j] = a.m[j][i];
    return r;
}

struct Llh {
    // 大地坐标：纬度、经度单位为 rad，高程单位为 m。
    double lat = 0.0;
    double lon = 0.0;
    double h = 0.0;
};

struct GTime {
    // GPS 周和周内秒。程序内部统一使用 GPST，避免年月日跨天处理。
    int week = 0;
    double tow = 0.0;
};

// 计算两个 GPS 时间之间的秒差 a-b。
double timediff(const GTime &a, const GTime &b);
// 规范化 GPS 周内秒，处理跨周情况。
GTime normalizeTime(GTime t);
// 将年月日时分秒转换为 GPS 周和周内秒。
GTime epochToGps(int year, int month, int day, int hour, int minute, double sec);
// 去掉字符串首尾空白，读取 RINEX 和 CSV 时使用。
std::string trim(const std::string &s);
// 将字符串转为 double，同时兼容 RINEX 中 D 指数格式。
double parseDouble(std::string s);
// 将字符串转为 int；失败时返回 fallback。
int parseInt(const std::string &s, int fallback = 0);
// ECEF 直角坐标转 WGS-84 经纬高。
Llh ecefToLlh(const Vec3 &ecef);
// WGS-84 经纬高转 ECEF 直角坐标。
Vec3 llhToEcef(const Llh &llh);
// 构造 ECEF 到当地 NED 坐标系的旋转矩阵。
Mat3 ecefToNedMatrix(const Llh &llh);
// 将 ECEF 坐标差转换到参考点处的 NED 增量。
Vec3 ecefToNed(const Vec3 &dr, const Llh &ref);
// 将周内秒归一到 [0,604800) 附近，便于匹配历元。
double wrapGpsTow(double t);
// 平方函数，避免代码中反复写 x*x 影响可读性。
double sqr(double x);

class Matrix {
public:
    // 简单动态矩阵类，只实现本程序需要的矩阵乘法、转置、求逆和线性方程求解。
    // 维度较小，主要用于 RTK 最小二乘和 ESKF 协方差传播/更新。
    Matrix() = default;
    Matrix(int rows, int cols, double value = 0.0);
    // 返回矩阵行数。
    int rows() const { return rows_; }
    // 返回矩阵列数。
    int cols() const { return cols_; }
    // 读写矩阵元素，采用行主序存储。
    double &operator()(int r, int c) { return data_[r * cols_ + c]; }
    // const 版本用于只读访问矩阵元素。
    double operator()(int r, int c) const { return data_[r * cols_ + c]; }
    // 构造 n 维单位矩阵。
    static Matrix identity(int n);
    // 返回当前矩阵的转置。
    Matrix transpose() const;

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<double> data_;
};

// 矩阵加法，用于协方差传播和法方程构造。
Matrix operator+(const Matrix &a, const Matrix &b);
// 矩阵减法，用于协方差更新。
Matrix operator-(const Matrix &a, const Matrix &b);
// 矩阵乘法，用于 Kalman 滤波、最小二乘和坐标线性化。
Matrix operator*(const Matrix &a, const Matrix &b);
// 矩阵乘向量，用于求残差修正量。
std::vector<double> operator*(const Matrix &a, const std::vector<double> &x);
// 高斯-约旦法求逆，小矩阵使用，失败会抛出异常。
Matrix inverse(Matrix a);
// 解线性方程 Ax=b，用于最小二乘法方程。
std::vector<double> solveLinear(Matrix a, std::vector<double> b);
// 普通最小二乘 dx=(H'H)^-1 H'v，可选输出参数协方差。
bool leastSquares(const Matrix &h, const std::vector<double> &v, std::vector<double> &dx, Matrix *q = nullptr);

} // namespace gnssins

