#include "common.h"

namespace gnssins {

// common.cpp
// ---------------------------------------------------------------------------
// 本文件实现全工程共享的基础函数。
//
// 主要分三类：
// 1. 时间工具：GPS 周/周内秒、公历时间转换、周内秒归一化；
// 2. 坐标工具：WGS-84 经纬高 <-> ECEF，ECEF 差值 -> NED；
// 3. 线性代数：小规模矩阵乘法、求逆、线性方程求解和最小二乘。
//
// 程序没有依赖 Eigen 等第三方库，因此 RTK 和 ESKF 所需的矩阵运算都在这里实现。
// ---------------------------------------------------------------------------
// 计算两个 GPS 时间的差值，返回 a-b 的秒数。
double timediff(const GTime &a, const GTime &b) {
    return (a.week - b.week) * 604800.0 + (a.tow - b.tow);
}

// 将 GPS 周内秒归一到一个周内，同时修正 GPS 周号。
GTime normalizeTime(GTime t) {
    while (t.tow < 0.0) {
        t.tow += 604800.0;
        --t.week;
    }
    while (t.tow >= 604800.0) {
        t.tow -= 604800.0;
        ++t.week;
    }
    return t;
}

// 将公历年月日转换为相对 Unix epoch 的天数，用于后续 GPS 时间转换。
static int daysFromCivil(int y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + static_cast<int>(doe) - 719468;
}

// 将年月日时分秒转换为 GPS 周和周内秒。
GTime epochToGps(int year, int month, int day, int hour, int minute, double sec) {
    // 把公历时间转成 GPS 周 / 周内秒，后续所有轨迹统一用这个时间基准。
    const int unixDays = daysFromCivil(year, static_cast<unsigned>(month), static_cast<unsigned>(day));
    const int gpsDays0 = daysFromCivil(1980, 1, 6);
    const double days = static_cast<double>(unixDays - gpsDays0);
    const double full = days * 86400.0 + hour * 3600.0 + minute * 60.0 + sec;
    GTime t;
    t.week = static_cast<int>(std::floor(full / 604800.0));
    t.tow = full - t.week * 604800.0;
    return normalizeTime(t);
}

// 去除字符串首尾空白字符，读取 RINEX/CSV 字段时使用。
std::string trim(const std::string &s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

// 字符串转 double，并兼容 RINEX 中常见的 D/d 科学计数法。
double parseDouble(std::string s) {
    for (char &c : s) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    s = trim(s);
    if (s.empty()) return 0.0;
    return std::stod(s);
}

// 字符串转 int；字段为空时返回 fallback。
int parseInt(const std::string &s, int fallback) {
    const std::string t = trim(s);
    if (t.empty()) return fallback;
    return std::stoi(t);
}

// WGS-84 经纬高转 ECEF 直角坐标。
Vec3 llhToEcef(const Llh &llh) {
    const double sinp = std::sin(llh.lat);
    const double cosp = std::cos(llh.lat);
    const double sinl = std::sin(llh.lon);
    const double cosl = std::cos(llh.lon);
    const double n = kA / std::sqrt(1.0 - kE2 * sinp * sinp);
    return {(n + llh.h) * cosp * cosl, (n + llh.h) * cosp * sinl, (n * (1.0 - kE2) + llh.h) * sinp};
}

// ECEF 直角坐标转 WGS-84 经纬高，纬度通过迭代求解。
Llh ecefToLlh(const Vec3 &r) {
    Llh out;
    out.lon = std::atan2(r.y, r.x);
    const double p = std::sqrt(r.x * r.x + r.y * r.y);
    out.lat = std::atan2(r.z, p * (1.0 - kE2));
    for (int i = 0; i < 8; ++i) {
        const double sinp = std::sin(out.lat);
        const double n = kA / std::sqrt(1.0 - kE2 * sinp * sinp);
        out.h = p / std::cos(out.lat) - n;
        out.lat = std::atan2(r.z, p * (1.0 - kE2 * n / (n + out.h)));
    }
    const double sinp = std::sin(out.lat);
    const double n = kA / std::sqrt(1.0 - kE2 * sinp * sinp);
    out.h = p / std::cos(out.lat) - n;
    return out;
}

// 构造 ECEF 到当地 NED 坐标系的方向余弦矩阵。
Mat3 ecefToNedMatrix(const Llh &llh) {
    const double sphi = std::sin(llh.lat), cphi = std::cos(llh.lat);
    const double slam = std::sin(llh.lon), clam = std::cos(llh.lon);
    Mat3 c;
    c.m[0][0] = -sphi * clam; c.m[0][1] = -sphi * slam; c.m[0][2] = cphi;
    c.m[1][0] = -slam;        c.m[1][1] = clam;         c.m[1][2] = 0.0;
    c.m[2][0] = -cphi * clam; c.m[2][1] = -cphi * slam; c.m[2][2] = -sphi;
    return c;
}

// 将 ECEF 坐标差投影到参考点处的 NED 坐标系。
Vec3 ecefToNed(const Vec3 &dr, const Llh &ref) {
    return matVec(ecefToNedMatrix(ref), dr);
}

// 将两个周内秒之差约束到半周范围，避免跨周比较错误。
double wrapGpsTow(double t) {
    while (t > 302400.0) t -= 604800.0;
    while (t < -302400.0) t += 604800.0;
    return t;
}

// 返回 x 的平方，主要用于方差计算。
double sqr(double x) { return x * x; }

// 创建 rows x cols 矩阵，并用 value 初始化所有元素。
Matrix::Matrix(int rows, int cols, double value) : rows_(rows), cols_(cols), data_(rows * cols, value) {}

// 构造 n 维单位矩阵。
Matrix Matrix::identity(int n) {
    Matrix a(n, n);
    for (int i = 0; i < n; ++i) a(i, i) = 1.0;
    return a;
}

// 返回当前矩阵的转置矩阵。
Matrix Matrix::transpose() const {
    Matrix t(cols_, rows_);
    for (int i = 0; i < rows_; ++i) for (int j = 0; j < cols_; ++j) t(j, i) = (*this)(i, j);
    return t;
}

// 矩阵逐元素相加。
Matrix operator+(const Matrix &a, const Matrix &b) {
    Matrix c(a.rows(), a.cols());
    for (int i = 0; i < a.rows(); ++i) for (int j = 0; j < a.cols(); ++j) c(i, j) = a(i, j) + b(i, j);
    return c;
}

// 矩阵逐元素相减。
Matrix operator-(const Matrix &a, const Matrix &b) {
    Matrix c(a.rows(), a.cols());
    for (int i = 0; i < a.rows(); ++i) for (int j = 0; j < a.cols(); ++j) c(i, j) = a(i, j) - b(i, j);
    return c;
}

// 矩阵乘法，用三重循环实现小规模矩阵相乘。
Matrix operator*(const Matrix &a, const Matrix &b) {
    Matrix c(a.rows(), b.cols());
    for (int i = 0; i < a.rows(); ++i) {
        for (int k = 0; k < a.cols(); ++k) {
            const double aik = a(i, k);
            for (int j = 0; j < b.cols(); ++j) c(i, j) += aik * b(k, j);
        }
    }
    return c;
}

// 矩阵乘向量，用于 Kalman 增益作用到残差和法方程计算。
std::vector<double> operator*(const Matrix &a, const std::vector<double> &x) {
    std::vector<double> y(a.rows(), 0.0);
    for (int i = 0; i < a.rows(); ++i) for (int j = 0; j < a.cols(); ++j) y[i] += a(i, j) * x[j];
    return y;
}

// 矩阵求逆。采用带主元选择的 Gauss-Jordan 消元，适合本工程的小矩阵。
Matrix inverse(Matrix a) {
    // Gauss-Jordan matrix inverse with partial pivoting.
    // Used by Kalman gain computation and least-squares covariance estimation.
    const int n = a.rows();
    Matrix inv = Matrix::identity(n);
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double best = std::abs(a(i, i));
        for (int r = i + 1; r < n; ++r) {
            if (std::abs(a(r, i)) > best) {
                best = std::abs(a(r, i));
                pivot = r;
            }
        }
        if (best < 1e-18) throw std::runtime_error("singular matrix");
        if (pivot != i) {
            for (int c = 0; c < n; ++c) {
                std::swap(a(i, c), a(pivot, c));
                std::swap(inv(i, c), inv(pivot, c));
            }
        }
        const double d = a(i, i);
        for (int c = 0; c < n; ++c) {
            a(i, c) /= d;
            inv(i, c) /= d;
        }
        for (int r = 0; r < n; ++r) {
            if (r == i) continue;
            const double f = a(r, i);
            for (int c = 0; c < n; ++c) {
                a(r, c) -= f * a(i, c);
                inv(r, c) -= f * inv(i, c);
            }
        }
    }
    return inv;
}

// 解线性方程组 A x = b，失败时返回空向量。
std::vector<double> solveLinear(Matrix a, std::vector<double> b) {
    // Solve A x = b by Gaussian elimination.
    // Return an empty vector if the matrix is close to singular.
    const int n = a.rows();
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        double best = std::abs(a(i, i));
        for (int r = i + 1; r < n; ++r) {
            if (std::abs(a(r, i)) > best) {
                best = std::abs(a(r, i));
                pivot = r;
            }
        }
        if (best < 1e-18) return {};
        if (pivot != i) {
            for (int c = i; c < n; ++c) std::swap(a(i, c), a(pivot, c));
            std::swap(b[i], b[pivot]);
        }
        for (int r = i + 1; r < n; ++r) {
            const double f = a(r, i) / a(i, i);
            for (int c = i; c < n; ++c) a(r, c) -= f * a(i, c);
            b[r] -= f * b[i];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int c = i + 1; c < n; ++c) s -= a(i, c) * x[c];
        x[i] = s / a(i, i);
    }
    return x;
}

// 最小二乘求解 H dx ~= v，可选返回参数协方差近似矩阵。
bool leastSquares(const Matrix &h, const std::vector<double> &v, std::vector<double> &dx, Matrix *q) {
    // Least-squares solver for H dx ~= v.
    // It forms the normal equation (H^T H) dx = H^T v.
    // q optionally stores inverse(H^T H), used as covariance-like quality information.
    const Matrix ht = h.transpose();
    Matrix n = ht * h;
    std::vector<double> u = ht * v;
    dx = solveLinear(n, u);
    if (dx.empty()) return false;
    if (q) {
        try {
            *q = inverse(n);
        } catch (...) {
            *q = Matrix(n.rows(), n.cols());
        }
    }
    return true;
}

} // namespace gnssins

