#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// imu_alignment_analysis.cpp
// ---------------------------------------------------------------------------
// 独立 IMU/初始对准分析程序，不参与 main.cpp 的主解算链。
//
// 用途：
// 1. 读取动态 IMR 和静态 AA77 IMU 数据；
// 2. 计算调平+直接罗经、解析法、动态航向辅助三种初始对准结果；
// 3. 与商业软件 POS/GNSSINS 参考结果中的姿态进行对比；
// 4. 输出报告中使用的 IMU 解码表和对准结果表。
// ---------------------------------------------------------------------------

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kD2R = kPi / 180.0;
constexpr double kR2D = 180.0 / kPi;
constexpr double kOmegaEarth = 7.2921151467e-5;

struct Vec3 {
    // 三维向量。此独立文件不依赖 common.h，因此保留一个最小版本。
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct ImuSample {
    // IMU 样本：时间、采样间隔、角增量和速度增量。
    double tow = 0.0;
    double dt = 0.005;
    Vec3 dtheta;
    Vec3 dvel;
};

struct RefPose {
    // 参考 POS 中抽取的姿态、速度和纬度，用于评价初始对准结果。
    double tow = 0.0;
    double latDeg = 0.0;
    double vn = 0.0;
    double ve = 0.0;
    double rollDeg = 0.0;
    double pitchDeg = 0.0;
    double yawDeg = 0.0;
};

// 从小端二进制缓冲区读取一个基础类型。
template <typename T>
T readLe(const char *p) {
    T v{};
    std::memcpy(&v, p, sizeof(T));
    return v;
}

// 从大端字节序读取 int32，静态 AA77 数据使用该格式。
int32_t readBeI32(const unsigned char *p) {
    uint32_t u = (static_cast<uint32_t>(p[0]) << 24) |
                 (static_cast<uint32_t>(p[1]) << 16) |
                 (static_cast<uint32_t>(p[2]) << 8) |
                 static_cast<uint32_t>(p[3]);
    if (u >= 0x80000000u) return static_cast<int32_t>(u - 0x100000000ull);
    return static_cast<int32_t>(u);
}

// 三维向量加法。
Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
// 三维向量减法。
Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
// 三维向量乘标量。
Vec3 operator*(const Vec3 &a, double s) { return {a.x * s, a.y * s, a.z * s}; }
// 三维向量除标量。
Vec3 operator/(const Vec3 &a, double s) { return {a.x / s, a.y / s, a.z / s}; }

// 计算三维向量二范数。
double norm(const Vec3 &a) {
    return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// 将角度归一到 [0, 360)。
double wrap360(double x) {
    while (x < 0.0) x += 360.0;
    while (x >= 360.0) x -= 360.0;
    return x;
}

// 将角度差归一到 [-180, 180]，用于姿态误差统计。
double wrap180(double x) {
    while (x > 180.0) x -= 360.0;
    while (x < -180.0) x += 360.0;
    return x;
}

// 将“度 + 分”格式转换为十进制度。
double parseDegMin(const std::string &degText, const std::string &minText) {
    const double deg = std::stod(degText);
    const double minute = std::stod(minText);
    return deg + (deg < 0.0 ? -minute : minute) / 60.0;
}

// 一次性读取整个二进制文件。
std::vector<char> readAllBytes(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open: " + path);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// 读取动态实验 IMR 文件，并转换为角速度/比力可分析的 IMU 样本。
std::vector<ImuSample> readImr(const std::string &path) {
    const std::vector<char> bytes = readAllBytes(path);
    if (bytes.size() < 288 || std::string(bytes.data(), bytes.data() + 7) != "$IMURAW") {
        throw std::runtime_error("IMR header is not $IMURAW: " + path);
    }

    constexpr size_t header = 256;
    constexpr size_t recSize = 32;
    constexpr double gyroScale = 1.0e-11;
    constexpr double accScale = 1.0e-9;
    const size_t n = (bytes.size() - header) / recSize;

    std::vector<ImuSample> out;
    out.reserve(n);
    double lastTow = -1.0;
    for (size_t i = 0; i < n; ++i) {
        const char *p = bytes.data() + header + i * recSize;
        const double tow = readLe<double>(p);
        if (!std::isfinite(tow) || tow <= 0.0) continue;
        ImuSample s;
        s.tow = tow;
        s.dt = (lastTow > 0.0 && tow > lastTow && tow - lastTow < 0.1) ? tow - lastTow : 0.005;
        lastTow = tow;
        s.dtheta = {readLe<int32_t>(p + 8) * gyroScale,
                    readLe<int32_t>(p + 12) * gyroScale,
                    readLe<int32_t>(p + 16) * gyroScale};
        s.dvel = {readLe<int32_t>(p + 20) * accScale,
                  readLe<int32_t>(p + 24) * accScale,
                  readLe<int32_t>(p + 28) * accScale};
        out.push_back(s);
    }

    if (!out.empty()) {
        // 用前 10 s 静态段把加速度计尺度归一到重力模长附近。
        Vec3 fb;
        int count = 0;
        const double endTow = out.front().tow + 10.0;
        for (const auto &s : out) {
            if (s.tow > endTow) break;
            fb = fb + s.dvel / s.dt;
            ++count;
        }
        const double fbNorm = count > 0 ? norm(fb / static_cast<double>(count)) : 0.0;
        if (fbNorm > 1.0) {
            const double scale = 9.80665 / fbNorm;
            for (auto &s : out) s.dvel = s.dvel * scale;
        }
    }
    return out;
}

// 读取静态 AA77 帧格式 IMU 文件，用于静态 Allan/活动状态分析。
std::vector<ImuSample> readAa77Static(const std::string &path, double startTow) {
    const std::vector<char> bytes = readAllBytes(path);
    constexpr size_t frameOffset = 28;
    constexpr size_t frameLen = 50;
    constexpr double fs = 200.0;
    constexpr double dt = 1.0 / fs;
    constexpr double gyroScale = 1.0;
    constexpr double accScale = 1.0 / 16777216.0;
    constexpr int gyroPos[3] = {20, 24, 28};
    constexpr int accPos[3] = {32, 36, 40};

    if (bytes.size() <= frameOffset + frameLen) throw std::runtime_error("AA77 file is too short");
    const size_t n = (bytes.size() - frameOffset) / frameLen;
    std::vector<ImuSample> out;
    out.reserve(n);
    const auto *data = reinterpret_cast<const unsigned char *>(bytes.data());
    for (size_t i = 0; i < n; ++i) {
        const unsigned char *f = data + frameOffset + i * frameLen;
        if (f[0] != 0xAA || f[1] != 0x77) continue;
        ImuSample s;
        s.tow = startTow + static_cast<double>(out.size()) * dt;
        s.dt = dt;
        s.dtheta = {readBeI32(f + gyroPos[0]) * gyroScale * dt,
                    readBeI32(f + gyroPos[1]) * gyroScale * dt,
                    readBeI32(f + gyroPos[2]) * gyroScale * dt};
        s.dvel = {readBeI32(f + accPos[0]) * accScale * dt,
                  readBeI32(f + accPos[1]) * accScale * dt,
                  readBeI32(f + accPos[2]) * accScale * dt};
        out.push_back(s);
    }
    return out;
}

// 读取参考 GNSS/INS POS 文件中的纬度、速度和姿态。
std::vector<RefPose> readReference(const std::string &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open reference: " + path);
    std::vector<RefPose> refs;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::vector<std::string> f;
        std::string item;
        while (iss >> item) f.push_back(item);
        if (f.size() < 17 || f[0].find('/') == std::string::npos) continue;
        try {
            RefPose r;
            r.tow = std::stod(f[2]);
            r.latDeg = parseDegMin(f[3], f[4]);
            r.ve = std::stod(f[8]);
            r.vn = std::stod(f[9]);
            r.rollDeg = parseDegMin(f[11], f[12]);
            r.pitchDeg = parseDegMin(f[13], f[14]);
            r.yawDeg = parseDegMin(f[15], f[16]);
            refs.push_back(r);
        } catch (...) {
        }
    }
    return refs;
}

// 查找距离指定 TOW 最近的参考姿态。
RefPose nearestReference(const std::vector<RefPose> &refs, double tow) {
    if (refs.empty()) throw std::runtime_error("empty reference");
    const auto it = std::min_element(refs.begin(), refs.end(), [tow](const RefPose &a, const RefPose &b) {
        return std::abs(a.tow - tow) < std::abs(b.tow - tow);
    });
    return *it;
}

// 计算指定时间窗内的平均陀螺角速度和平均加速度。
std::pair<Vec3, Vec3> meanGyroAcc(const std::vector<ImuSample> &imu, double start, double end) {
    Vec3 gyro;
    Vec3 acc;
    int count = 0;
    for (const auto &s : imu) {
        if (s.tow < start || s.tow > end || s.dt <= 0.0) continue;
        gyro = gyro + s.dtheta / s.dt;
        acc = acc + s.dvel / s.dt;
        ++count;
    }
    if (count == 0) throw std::runtime_error("empty IMU window");
    return {gyro / static_cast<double>(count), acc / static_cast<double>(count)};
}

// 用静态加速度计重力方向计算横滚角和俯仰角。
std::pair<double, double> levelFromAcc(const Vec3 &acc) {
    const double roll = std::atan2(acc.y, acc.z) * kR2D;
    const double pitch = std::atan2(-acc.x, std::sqrt(acc.y * acc.y + acc.z * acc.z)) * kR2D;
    return {roll, pitch};
}

// 直接罗经法航向估计。MEMS 陀螺精度有限，该结果主要作为方法对比。
double yawDirectGyro(const Vec3 &gyro, double rollDeg, double pitchDeg) {
    const double r = rollDeg * kD2R;
    const double p = pitchDeg * kD2R;
    const double cr = std::cos(r), sr = std::sin(r);
    const double cp = std::cos(p), sp = std::sin(p);
    const double wxLevel = cp * gyro.x + sp * sr * gyro.y + sp * cr * gyro.z;
    const double wyLevel = cr * gyro.y - sr * gyro.z;
    return wrap360(std::atan2(-wyLevel, wxLevel) * kR2D);
}

// 解析法航向搜索：寻找使地球自转投影与平均陀螺最接近的航向角。
double yawAnalyticSearch(const Vec3 &gyro, double latDeg, double rollDeg, double pitchDeg) {
    const double lat = latDeg * kD2R;
    const Vec3 wnav{kOmegaEarth * std::cos(lat), 0.0, -kOmegaEarth * std::sin(lat)};
    const double r = rollDeg * kD2R;
    const double p = pitchDeg * kD2R;
    const double cr = std::cos(r), sr = std::sin(r);
    const double cp = std::cos(p), sp = std::sin(p);
    double bestYaw = 0.0;
    double bestErr = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= 7200; ++i) {
        const double yaw = i * 0.05;
        const double y = yaw * kD2R;
        const double cy = std::cos(y), sy = std::sin(y);
        const double c00 = cp * cy;
        const double c01 = sr * sp * cy - cr * sy;
        const double c02 = cr * sp * cy + sr * sy;
        const double c10 = cp * sy;
        const double c11 = sr * sp * sy + cr * cy;
        const double c12 = cr * sp * sy - sr * cy;
        const double c20 = -sp;
        const double c21 = sr * cp;
        const double c22 = cr * cp;
        const Vec3 predBody{
            c00 * wnav.x + c10 * wnav.y + c20 * wnav.z,
            c01 * wnav.x + c11 * wnav.y + c21 * wnav.z,
            c02 * wnav.x + c12 * wnav.y + c22 * wnav.z
        };
        const double err = norm(predBody - gyro);
        if (err < bestErr) {
            bestErr = err;
            bestYaw = yaw;
        }
    }
    return wrap360(bestYaw);
}

// 动态航向辅助：车辆有明显水平速度时，用速度方向估计航向初值。
double yawDynamicCourse(const std::vector<RefPose> &refs, double start, double end) {
    double sx = 0.0;
    double sy = 0.0;
    int count = 0;
    for (const auto &r : refs) {
        const double speed = std::sqrt(r.vn * r.vn + r.ve * r.ve);
        if (r.tow < start || r.tow > end || speed < 0.8) continue;
        const double course = std::atan2(r.ve, r.vn);
        sx += std::cos(course) * speed;
        sy += std::sin(course) * speed;
        ++count;
        if (count >= 100) break;
    }
    if (count == 0) return 0.0;
    return wrap360(std::atan2(sy, sx) * kR2D);
}

// 输出 IMU 解码 CSV，包括角速度、加速度模长等检查量。
void writeDecodedCsv(const std::string &path, const std::vector<ImuSample> &imu) {
    std::ofstream out(path);
    out << "tow,dt,gyro_x,gyro_y,gyro_z,acc_x,acc_y,acc_z,acc_norm,gyro_norm\n";
    out << std::fixed << std::setprecision(10);
    for (const auto &s : imu) {
        const Vec3 g = s.dtheta / s.dt;
        const Vec3 a = s.dvel / s.dt;
        out << s.tow << ',' << s.dt << ','
            << g.x << ',' << g.y << ',' << g.z << ','
            << a.x << ',' << a.y << ',' << a.z << ','
            << norm(a) << ',' << norm(g) << '\n';
    }
}

// 输出三种初始对准方法的 roll/pitch/yaw 及其相对参考姿态误差。
void writeAlignmentCsv(const std::string &path, const std::vector<ImuSample> &dynamicImu,
                       const std::vector<RefPose> &refs) {
    const double start = std::max(dynamicImu.front().tow, refs.front().tow);
    const double staticEnd = start + 300.0;
    const double fineEnd = start + 600.0;
    const auto [gyro, acc] = meanGyroAcc(dynamicImu, start, staticEnd);
    const auto [roll, pitch] = levelFromAcc(acc);
    const RefPose truth = nearestReference(refs, start + 150.0);
    const double yawDirect = yawDirectGyro(gyro, roll, pitch);
    const double yawAnalytic = yawAnalyticSearch(gyro, truth.latDeg, roll, pitch);
    const double yawCourse = yawDynamicCourse(refs, start, fineEnd);

    std::ofstream out(path);
    out << "method,segment_start,segment_end,roll_deg,pitch_deg,yaw_deg,"
           "ref_roll_deg,ref_pitch_deg,ref_yaw_deg,err_roll_deg,err_pitch_deg,err_yaw_deg\n";
    const auto row = [&](const std::string &name, double yaw) {
        out << name << ',' << std::fixed << std::setprecision(6)
            << start << ',' << staticEnd << ','
            << roll << ',' << pitch << ',' << yaw << ','
            << truth.rollDeg << ',' << truth.pitchDeg << ',' << truth.yawDeg << ','
            << wrap180(roll - truth.rollDeg) << ','
            << wrap180(pitch - truth.pitchDeg) << ','
            << wrap180(yaw - truth.yawDeg) << '\n';
    };
    row("leveling_plus_direct_gyrocompass", yawDirect);
    row("analytic_two_vector", yawAnalytic);
    row("dynamic_course_fine_alignment", yawCourse);
}

} // namespace

// 独立程序入口：运行 IMU 解码和三种初始对准方法对比。
int main(int argc, char **argv) {
    try {
        const std::string root = argc > 1 ? argv[1] : ".";
        const std::string dynamicImuPath = root + "/data/imu_20260604_105437_599404.dat.imr";
        const std::string staticImuPath = root + "/data/IMUjingtai/imu_20260603_100606_559947.dat";
        const std::string refPath = root + "/reference/GNSSINS/gnss_20260604_105437_599404.pos";
        const std::string outDir = root + "/output/成果_MATLAB分析/IMU活动与对准分析";

        const std::vector<ImuSample> dynamicImu = readImr(dynamicImuPath);
        const std::vector<ImuSample> staticImu = readAa77Static(staticImuPath, 0.0);
        const std::vector<RefPose> refs = readReference(refPath);

        writeDecodedCsv(outDir + "/cpp_dynamic_imu_decoded.csv", dynamicImu);
        writeDecodedCsv(outDir + "/cpp_static_imu_decoded.csv", staticImu);
        writeAlignmentCsv(outDir + "/cpp_alignment_three_methods.csv", dynamicImu, refs);

        std::cout << "dynamic samples: " << dynamicImu.size() << "\n";
        std::cout << "static samples: " << staticImu.size() << "\n";
        std::cout << "outputs: " << outDir << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
