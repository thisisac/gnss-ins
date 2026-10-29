#include "ins.h"
#include <cstring>

namespace gnssins {

// ins.cpp
// ---------------------------------------------------------------------------
// 本文件实现惯性导航相关的基础处理。
//
// 主要内容：
// 1. decodeImr：读取动态实验 IMR 二进制文件，转换为 IMU 角增量和速度增量；
// 2. coarseAlign：利用起始静态段加速度计平均值估计 roll/pitch，
//    航向 yaw 使用外业记录值或默认值；
// 3. mechanize：纯惯导机械编排，用 IMU 增量递推姿态、速度和位置；
// 4. normalGravity / earthRateNed / transportRateNed：机械编排所需的重力、
//    地球自转和导航系运动角速度模型。
//
// 在线松组合中，EsKf::predict 实现了同类机械编排并额外传播协方差。
// ---------------------------------------------------------------------------
// 从二进制缓冲区中按小端格式读取一个基础类型。
template <typename T> static T readLe(const char *p) {
    T v; std::memcpy(&v, p, sizeof(T)); return v;
}

// 解码动态实验 IMR 二进制文件，生成按时间排序的 IMU 增量序列。
std::vector<ImuSample> decodeImr(const std::string &path, const std::string &textOut) {
    // 把二进制 IMR 转成程序内部能递推的 IMU 增量序列。
    // IMR 文件中每条记录包含时间、三轴角增量、三轴速度增量。
    // 解码后得到的 ImuSample 会直接送入机械编排或 ESKF predict。
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open IMR file: " + path);
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (bytes.size() < 288 || std::string(bytes.data(), bytes.data() + 7) != "$IMURAW") {
        throw std::runtime_error("IMR header is not $IMURAW");
    }
    constexpr size_t header = 256, recSize = 32;
    const size_t n = (bytes.size() - header) / recSize;
    std::vector<ImuSample> imu; imu.reserve(n);
    double lastTime = -1;
    constexpr double gyroScale = 1.0e-11;
    constexpr double accScale  = 1.0e-9;
    for (size_t i = 0; i < n; ++i) {
        const char *p = bytes.data() + header + i * recSize;
        double t = readLe<double>(p);
        if (!std::isfinite(t) || t <= 0) continue;
        ImuSample s;
        s.time = t;
        s.dtheta = {readLe<int32_t>(p +  8) * gyroScale, readLe<int32_t>(p + 12) * gyroScale,
                    readLe<int32_t>(p + 16) * gyroScale};
        s.dvel   = {readLe<int32_t>(p + 20) * accScale,  readLe<int32_t>(p + 24) * accScale,
                    readLe<int32_t>(p + 28) * accScale};
        s.dt = (lastTime > 0 && t > lastTime && t - lastTime < 0.1) ? t - lastTime : 0.005;
        lastTime = t;
        imu.push_back(s);
    }
    if (!imu.empty()) {
        // 根据起始静态段的加速度模长做一次比例归一化。
        // 这一步用于把设备原始增量尺度调整到接近 9.80665 m/s^2，
        // 避免后续惯导传播因量纲偏差快速发散。
        const double staticEnd = imu.front().time + 10.0;
        Vec3 fbSum;
        int cnt = 0;
        for (const auto &s : imu) {
            if (s.time > staticEnd) break;
            if (s.dt <= 0.0) continue;
            fbSum = fbSum + s.dvel / s.dt;
            ++cnt;
        }
        const double fbNorm = cnt > 0 ? norm(fbSum / cnt) : 0.0;
        if (fbNorm > 1.0) {
            const double scale = 9.80665 / fbNorm;
            for (auto &s : imu) s.dvel = s.dvel * scale;
        }
    }
    writeInsSamples(textOut, imu);
    return imu;
}

// 将解码后的 IMU 增量保存为文本，便于检查时间戳、角增量和速度增量。
void writeInsSamples(const std::string &path, const std::vector<ImuSample> &imu) {
    std::ofstream out(path);
    out << "# tow dtheta_x dtheta_y dtheta_z dvel_x dvel_y dvel_z\n";
    out << std::fixed << std::setprecision(10);
    for (const auto &s : imu)
        out << s.time << ' ' << s.dtheta.x << ' ' << s.dtheta.y << ' ' << s.dtheta.z << ' '
            << s.dvel.x << ' ' << s.dvel.y << ' ' << s.dvel.z << '\n';
}

// 根据纬度和高程计算正常重力，用于惯导速度更新中的重力补偿。
double normalGravity(double lat, double h) {
    double s2 = std::sin(lat); s2 *= s2;
    return 9.7803267715 * (1 + 0.0052790414 * s2 + 0.0000232718 * s2 * s2)
           + h * (-3.0876910891e-6 + 4.3977311e-9 * s2 - 7.211e-13 * h);
}

// 将地球自转角速度投影到当地 NED 坐标系。
Vec3 earthRateNed(double lat) {
    return {kOmegaEarth * std::cos(lat), 0.0, -kOmegaEarth * std::sin(lat)};
}

// 计算导航系相对于地球系的运输角速度，车辆运动时用于科氏项补偿。
Vec3 transportRateNed(double lat, double h, const Vec3 &vn) {
    double s2 = std::sin(lat); s2 *= s2;
    double Rn = kA / std::sqrt(1 - kE2 * s2);
    double Rm = Rn * (1 - kE2) / (1 - kE2 * s2);
    return {vn.y / (Rn + h), -vn.x / (Rm + h), -vn.y * std::tan(lat) / (Rn + h)};
}

NavState coarseAlign(const std::vector<ImuSample> &imu, const GnssSolution &init,
                     double staticSeconds, double defaultYawDeg) {
    // 静态粗对准：
    // 1. 前 staticSeconds 秒内平均加速度计输出，估计重力在载体系下的方向；
    // 2. 由重力方向计算 roll/pitch；
    // 3. MEMS 陀螺难以可靠测地球自转，因此 yaw 使用外业记录/默认值；
    // 4. GNSS 初始位置速度作为滤波初值的位置速度。
    if (imu.empty()) throw std::runtime_error("empty IMU");
    double end = imu.front().time + staticSeconds;
    Vec3 fbSum;
    Vec3 wbSum;
    int n = 0;
    for (const auto &s : imu) {
        if (s.time > end) break;
        if (s.dt > 0) {
            fbSum = fbSum + s.dvel / s.dt;
            wbSum = wbSum + s.dtheta / s.dt;
            ++n;
        }
    }
    if (n == 0) throw std::runtime_error("no IMU for coarse align");
    Vec3 fb = fbSum / n;

    // 静态粗对准只负责姿态起步，重力能定 roll/pitch，yaw 需要外部辅助。
    // 这里采用 NED/车体系约定下的加速度计调平公式。
    double roll  = std::atan2(fb.y, fb.z);  
    double pitch = std::atan2(-fb.x, std::sqrt(fb.y * fb.y + fb.z * fb.z));
    NavState st;
    st.time = init.time;
    st.time = init.time;
    st.llh = init.llh;
    st.ecef = init.ecef;
    st.velNed = init.velNed;
    st.qbn = eulerToQuat(roll, pitch, defaultYawDeg * kD2R);
    st.accBias = {};

    // 平均陀螺输出只作为初始零偏估计，不作为可靠罗经航向。
    // MEMS 陀螺精度不足时，地球自转角速度容易被噪声淹没。
    st.gyroBias = wbSum / n;
    return st;
}

// 纯惯导机械编排：用单帧 IMU 增量递推姿态、速度和位置。
void mechanize(NavState &s, const ImuSample &imu) {
    // 纯惯导机械编排：只使用 IMU 增量递推姿态、速度和位置。
    // 松组合在线模式中，EsKf::predict 内部实现了同类递推，同时传播协方差。
    double dt = imu.dt > 0 ? imu.dt : 0.005;

    // 先扣除零偏，再做惯导递推。
    Vec3 dtheta = imu.dtheta - s.gyroBias * dt;
    Vec3 dvel   = imu.dvel   - s.accBias  * dt;

    // 圆锥补偿：相邻角增量不可交换会产生二阶误差，使用简化一阶项补偿。
    static Vec3 prevDtheta;
    Vec3 dthetaCorr = dtheta + cross(prevDtheta, dtheta) / 12.0;
    prevDtheta = dtheta;

    // 姿态更新：角增量转四元数小旋转，再与当前姿态相乘。
    double dmag = norm(dthetaCorr);
    Quat dq = dmag > 1e-12 ? quatFromAxisAngle(dthetaCorr / dmag, dmag) : Quat{1, 0, 0, 0};
    s.qbn = quatNormalize(quatMul(s.qbn, dq));

    // 划桨补偿：降低角运动和速度增量耦合带来的速度更新误差。
    static Vec3 prevDvel;
    Vec3 dvelCorr = dvel + cross(prevDtheta, dvel + prevDvel / 6.0) / 2.0;
    prevDvel = dvel;

    // 速度更新：机体系 delta-v 转到导航系，加入重力并扣除地球自转/运输率影响。
    Mat3 Cbn = quatToDcm(s.qbn);
    Vec3 fn = matVec(Cbn, dvelCorr) / dt;

    Vec3 wie_n = earthRateNed(s.llh.lat);
    Vec3 wen_n = transportRateNed(s.llh.lat, s.llh.h, s.velNed);
    Vec3 coriolis = cross(wie_n * 2.0 + wen_n, s.velNed);
    Vec3 gn = {0.0, 0.0, normalGravity(s.llh.lat, s.llh.h)};

    s.velNed = s.velNed + (fn - coriolis + gn) * dt;

    // 位置更新：用子午圈、卯酉圈半径把 NED 速度积分为经纬高变化。
    double s2 = std::sin(s.llh.lat); s2 *= s2;
    double Rn = kA / std::sqrt(1 - kE2 * s2);
    double Rm = Rn * (1 - kE2) / (1 - kE2 * s2);
    s.llh.lat += s.velNed.x / (Rm + s.llh.h) * dt;
    s.llh.lon += s.velNed.y / ((Rn + s.llh.h) * std::cos(s.llh.lat)) * dt;
    s.llh.h  -= s.velNed.z * dt;

    s.ecef = llhToEcef(s.llh);
    s.time.tow += dt;
}

} // namespace gnssins






