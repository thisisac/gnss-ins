#pragma once

#include "common.h"
#include "gnss.h"

namespace gnssins {

// INS 模块负责 IMU 解码、粗对准和惯导递推。
// 这里的设计目标是：先把姿态/速度/位置跑起来，再交给松组合修正。
struct ImuSample {
    double time = 0.0; // IMU 记录时刻，使用 GPS 周内秒 [s]
    double dt = 0.0;   // 与上一帧的时间间隔，动态数据约 0.005 s
    Vec3 dtheta;       // 陀螺角增量 [rad]
    Vec3 dvel;         // 加速度计速度增量 [m/s]
};

struct Quat {
    // 四元数用于姿态传播，避免欧拉角直接积分产生奇异。
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;
};

// 四元数乘法，用于把本次陀螺增量对应的小旋转叠加到当前姿态上。
inline Quat quatMul(const Quat &a, const Quat &b) {
    return {
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
        a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    };
}

// 四元数共轭；单位四元数的共轭等价于姿态逆变换。
inline Quat quatConj(const Quat &q) { return {q.w, -q.x, -q.y, -q.z}; }

// 四元数归一化，抑制数值积分后模长偏离 1 的问题。
inline Quat quatNormalize(const Quat &q) {
    double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    return n > 0 ? Quat{q.w / n, q.x / n, q.y / n, q.z / n} : Quat{1, 0, 0, 0};
}

// 由旋转轴和旋转角构造四元数，陀螺角增量更新姿态时使用。
inline Quat quatFromAxisAngle(const Vec3 &axis, double angle) {
    double half = angle * 0.5;
    double s = std::sin(half);
    return {std::cos(half), axis.x * s, axis.y * s, axis.z * s};
}

// 由方向余弦矩阵构造四元数，粗对准或外部姿态初始化时使用。
inline Quat quatFromDcm(const Mat3 &C) {
    // 先判断矩阵迹，选择数值稳定的分支，避免在某个分量接近零时放大误差。
    double t = C.m[0][0] + C.m[1][1] + C.m[2][2];
    Quat q;
    if (t > 0) {
        double s = 0.5 / std::sqrt(t + 1);
        q.w = 0.25 / s;
        q.x = (C.m[2][1] - C.m[1][2]) * s;
        q.y = (C.m[0][2] - C.m[2][0]) * s;
        q.z = (C.m[1][0] - C.m[0][1]) * s;
    } else if (C.m[0][0] > C.m[1][1] && C.m[0][0] > C.m[2][2]) {
        double s = 2 * std::sqrt(1 + C.m[0][0] - C.m[1][1] - C.m[2][2]);
        q.w = (C.m[2][1] - C.m[1][2]) / s;
        q.x = 0.25 * s;
        q.y = (C.m[0][1] + C.m[1][0]) / s;
        q.z = (C.m[0][2] + C.m[2][0]) / s;
    } else if (C.m[1][1] > C.m[2][2]) {
        double s = 2 * std::sqrt(1 + C.m[1][1] - C.m[0][0] - C.m[2][2]);
        q.w = (C.m[0][2] - C.m[2][0]) / s;
        q.x = (C.m[0][1] + C.m[1][0]) / s;
        q.y = 0.25 * s;
        q.z = (C.m[1][2] + C.m[2][1]) / s;
    } else {
        double s = 2 * std::sqrt(1 + C.m[2][2] - C.m[0][0] - C.m[1][1]);
        q.w = (C.m[1][0] - C.m[0][1]) / s;
        q.x = (C.m[0][2] + C.m[2][0]) / s;
        q.y = (C.m[1][2] + C.m[2][1]) / s;
        q.z = 0.25 * s;
    }
    return quatNormalize(q);
}

// 将四元数转换为方向余弦矩阵，便于把 IMU 体坐标增量转到导航坐标系。
inline Mat3 quatToDcm(const Quat &q) {
    Mat3 C;
    double ww = q.w * q.w, xx = q.x * q.x, yy = q.y * q.y, zz = q.z * q.z;
    double wx = q.w * q.x, wy = q.w * q.y, wz = q.w * q.z;
    double xy = q.x * q.y, xz = q.x * q.z, yz = q.y * q.z;
    C.m[0][0] = ww + xx - yy - zz;  C.m[0][1] = 2 * (xy - wz);       C.m[0][2] = 2 * (xz + wy);
    C.m[1][0] = 2 * (xy + wz);       C.m[1][1] = ww - xx + yy - zz;  C.m[1][2] = 2 * (yz - wx);
    C.m[2][0] = 2 * (xz - wy);       C.m[2][1] = 2 * (yz + wx);       C.m[2][2] = ww - xx - yy + zz;
    return C;
}

// 四元数转欧拉角，输出 roll/pitch/yaw，主要用于结果保存和报告绘图。
inline Vec3 quatToEuler(const Quat &q) {
    Mat3 C = quatToDcm(q);
    double roll = std::atan2(C.m[2][1], C.m[2][2]);
    double pitch = -std::asin(C.m[2][0]);
    double yaw = std::atan2(C.m[1][0], C.m[0][0]);
    return {roll, pitch, yaw};
}

// 欧拉角转四元数，初始化姿态时使用。
inline Quat eulerToQuat(double roll, double pitch, double yaw) {
    double cr = std::cos(roll * 0.5), sr = std::sin(roll * 0.5);
    double cp = std::cos(pitch * 0.5), sp = std::sin(pitch * 0.5);
    double cy = std::cos(yaw * 0.5), sy = std::sin(yaw * 0.5);
    return {
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy
    };
}

// 构造向量对应的反对称矩阵，使 S(v)*a 等价于 v 与 a 的叉乘。
inline Mat3 skewSymmetric(const Vec3 &v) {
    Mat3 S;
    S.m[0][1] = -v.z; S.m[0][2] =  v.y;
    S.m[1][0] =  v.z; S.m[1][2] = -v.x;
    S.m[2][0] = -v.y; S.m[2][1] =  v.x;
    return S;
}

struct NavState {
    // 组合导航名义状态。ESKF 预测和更新后都维护这一组状态。
    GTime time;
    Llh llh;       // IMU 中心位置
    Vec3 ecef;     // IMU 中心 ECEF 坐标
    Vec3 velNed;   // NED 速度
    Quat qbn;      // body-to-NED quaternion
    Vec3 accBias;  // 加表零偏
    Vec3 gyroBias; // 陀螺零偏
    Vec3 accScale; // 加表比例因子误差
    Vec3 gyroScale;// 陀螺比例因子误差
};

struct ImuModel {
    // ESKF 初始方差和过程噪声参数。
    // 这些参数决定滤波器对初值、IMU 预测模型和随机漂移的信任程度。
    Vec3 initPosStd{100.0, 100.0, 100.0};
    Vec3 initVelStd{1.0, 1.0, 1.0};
    Vec3 initAttStd{5.0 * kD2R, 5.0 * kD2R, 15.0 * kD2R};
    Vec3 initAccBiasStd{2.5e-2, 2.5e-2, 2.5e-2};
    Vec3 initGyroBiasStd{1.0e-1 * kD2R, 1.0e-1 * kD2R, 1.0e-1 * kD2R};
    Vec3 initAccScaleStd{1.0e-3, 1.0e-3, 1.0e-3};
    Vec3 initGyroScaleStd{1.0e-3, 1.0e-3, 1.0e-3};
    Vec3 posNoise{1.0e-2, 1.0e-2, 1.0e-2};
    Vec3 velNoise{1.25e-3, 1.25e-3, 1.25e-3};
    Vec3 attNoise{3.54e-3 * kD2R, 3.54e-3 * kD2R, 3.54e-3 * kD2R};
    Vec3 accBiasNoise{5.62e-5, 5.62e-5, 5.62e-5};
    Vec3 gyroBiasNoise{2.72e-5 * kD2R, 2.72e-5 * kD2R, 2.72e-5 * kD2R};
    Vec3 accScaleNoise{1.0e-7, 1.0e-7, 1.0e-7};
    Vec3 gyroScaleNoise{1.0e-7, 1.0e-7, 1.0e-7};
};

// 解码动态 IMR 文件，输出 IMU 增量序列；textOut 用于保存可检查的文本版数据。
std::vector<ImuSample> decodeImr(const std::string &path, const std::string &textOut);
// coarseAlign 使用起始静态段估计 roll/pitch，yaw 使用外业记录或默认值。
NavState coarseAlign(const std::vector<ImuSample> &imu, const GnssSolution &init, double staticSeconds,
                     double defaultYawDeg);
// mechanize 是纯惯导机械编排函数，在线松组合实际通过 EsKf::predict 调用同类传播思想。
void mechanize(NavState &state, const ImuSample &imu);
// 保存解码后的 IMU 采样，便于核对帧率、时间戳和六轴增量。
void writeInsSamples(const std::string &path, const std::vector<ImuSample> &imu);
// 计算正常重力值，纯惯导机械编排中用于重力补偿。
double normalGravity(double lat, double h);
// 计算地球自转角速度在 NED 坐标系下的投影。
Vec3 earthRateNed(double lat);
// 计算导航坐标系相对地球坐标系的运输角速度。
Vec3 transportRateNed(double lat, double h, const Vec3 &vn);

} // namespace gnssins

