#include "fusion.h"
#include "gnss.h"
#include "ins.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

#ifdef _WIN32
#include <direct.h>
#else
#include <filesystem>
#endif

using namespace gnssins;

namespace {

// main.cpp
// ---------------------------------------------------------------------------
// 本文件是工程主调度程序。
//
// 它不直接实现某一个数学模型，而是把 GNSS、INS、Fusion 三个模块串起来：
// 1. 根据 AppConfig::mode 决定运行 RTK、前向松组合、离线平滑或伪实时在线解算；
// 2. 读取 GNSS/IMU/reference 数据路径；
// 3. 调用 GNSS 模块生成 RTK 历元解；
// 4. 调用 INS 模块解码 IMU 并完成粗对准；
// 5. 调用 Fusion 模块进行 ESKF 松组合；
// 6. pseudo-online 模式下按时间戳模拟传感器播发，逐帧在线处理。
//
// 现场源码演示时，最直接的接口就是 AppConfig::mode。
// ---------------------------------------------------------------------------

struct AppConfig {
    // =========================
    // 运行模式源码接口
    // =========================
    // 本工程现场展示时可以不依赖命令行参数，直接在 VS 中修改下面的 mode 默认值。
    // main() 会根据 mode 自动决定要运行哪些处理模块。
    //
    // mode 可选值及含义：
    //   rtk:
    //     只执行 GNSS RTK 解算。程序读取流动站、基准站和广播星历，
    //     输出每个 GNSS 历元的定位结果。
    //   forward:
    //     执行 RTK、IMU 解码、粗对准和 GNSS/INS 松组合前向滤波。
    //     该模式主要用于查看前向滤波结果。
    //   post:
    //     基于已有前向结果做离线后处理平滑。
    //   all:
    //     完整后处理流程，包括 RTK、前向松组合、后向/前后向平滑和诊断结果。
    //   pseudo-online:
    //     伪实时在线模式。程序把已有 IMU/GNSS 文件当作模拟传感器，
    //     按时间戳逐帧送入滤波器：IMU 做 predict，GNSS 到历元才 update。
    //   pseudo-online-fast:
    //     与 pseudo-online 算法完全相同，只是不按墙钟等待 1 秒，适合快速调试。
    //
    // 注意：pseudo-online 不是“先算完结果再回放”，它只使用当前及过去的数据。
    // 现场展示的主接口：
    // 如果只想在 VS 中改源码切换任务，就直接修改下面的 mode 默认值。
    // 程序启动后 main() 会根据 mode 判断需要执行哪些模块。
    //
    // 可选任务：
    //   "rtk"                只运行 GNSS RTK，输出 gnss_rtk.pos；
    //   "forward"            RTK + IMU 解码 + 粗对准 + 松组合前向滤波；
    //   "post"               基于已有前向结果做离线平滑；
    //   "all"                完整后处理流程；
    //   "pseudo-online"      伪实时在线模式，已有数据按时间戳模拟传感器播发；
    //   "pseudo-online-fast" 算法同上，但打印不等待真实 1 秒，便于快速检查。
    //
    // 注意：
    //   pseudo-online 与旧的结果回放不同，它不会先把组合导航全算完再播放；
    //   它按 IMU/GNSS 时间顺序推进 ESKF，只使用当前时刻及以前的数据。
    // 现场展示时，如果不想改命令行，就直接在这里修改 mode。
    // 例如：
    //   "rtk"                只跑 RTK；
    //   "forward"            跑 RTK + 松组合前向；
    //   "post"               读取已有 forward 结果做离线平滑；
    //   "all"                完整后处理流程；
    //   "pseudo-online"      按数据时间戳模拟传感器播发并在线解算；
    //   "pseudo-online-fast" 同上，但每秒打印时不等待真实 1 秒，便于快速调试。
    // Paths are relative to the repository root. Keep the Visual Studio
    // working directory set to the repository root when using these defaults.
    std::string root = ".";
    std::string rover = root + "/data/gnss_20260604_105437_599404.26o";
    std::string nav   = root + "/data/gnss_20260604_105437_599404.26p";
    std::string base  = root + "/data/SEPT00____R_20261550000_01D_01S_MO.rnx";
    std::string imu   = root + "/data/imu_20260604_105437_599404.dat.imr";
    std::string referenceGnssIns = root + "/reference/GNSSINS/gnss_20260604_105437_599404.pos";
    std::string referenceRtk = root + "/reference/rtk/gnss_20260604_105437_599404.pos";
    std::string outDir = root + "/output";
    std::string mode = "pseudo-online"; // rtk / forward / post / all / pseudo-online
    int maxEpochs = -1;
    double coarseAlignSeconds = 300.0;
    double defaultYawDeg = 190.0;
};

struct ReferenceState {
    double tow = 0.0;
    double h = 0.0;
    Vec3 velNed;
    double roll = 0.0;
    double pitch = 0.0;
    double yaw = 0.0;
    Vec3 dispNed;
};

struct OnlineRunResult {
    std::vector<FusionOutput> allOutputs;
    std::vector<FusionOutput> oneHzOutputs;
    int imuPackets = 0;
    int gnssPackets = 0;
};

// 确保输出目录存在。Windows 下使用 _mkdir，其他平台使用 filesystem。
void ensureDir(const std::string &path) {
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    std::filesystem::create_directories(path);
#endif
}

// 判断字符串是否为整数，用于区分第二个参数是 maxEpochs 还是 mode。
bool isIntegerText(const std::string &s) {
    if (s.empty()) return false;
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    for (; i < s.size(); ++i) {
        if (s[i] < '0' || s[i] > '9') return false;
    }
    return true;
}

// 将用户输入或源码中填写的模式名归一化为程序内部使用的模式名。
std::string normalizeMode(std::string mode) {
    // 模式别名归一化。最终程序只保留在线的 pseudo-online，
    // 不再保留“先算完再回放”的 pseudo-realtime。
    for (char &c : mode) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (mode == "rtk-only") return "rtk";
    if (mode == "real" || mode == "real-time" || mode == "online" || mode == "pseudo-sensor" ||
        mode == "pseudo-online") return "pseudo-online";
    if (mode == "real-fast" || mode == "online-fast" || mode == "pseudo-sensor-fast" ||
        mode == "pseudo-online-fast") return "pseudo-online-fast";
    if (mode == "post-only") return "post";
    if (mode == "rtk-all") return "all";
    return mode;
}

// 将“度 + 分”格式转换为十进制度，用于读取参考 POS 文件中的姿态角。
double parseDegMin(double deg, double min) {
    return deg + (deg < 0.0 ? -min : min) / 60.0;
}

// 按时间对参考结果做线性插值，用于误差统计和粗对准初值对比。
ReferenceState interpolateReference(const std::vector<ReferenceState> &refs, double tow, size_t &idx) {
    while (idx + 1 < refs.size() && refs[idx + 1].tow < tow) ++idx;
    if (idx + 1 >= refs.size()) return refs.back();
    const ReferenceState &a = refs[idx];
    const ReferenceState &b = refs[idx + 1];
    const double u = std::max(0.0, std::min(1.0, (tow - a.tow) / (b.tow - a.tow)));
    ReferenceState r;
    r.tow = tow;
    r.h = a.h * (1.0 - u) + b.h * u;
    r.velNed = a.velNed * (1.0 - u) + b.velNed * u;
    r.roll = a.roll * (1.0 - u) + b.roll * u;
    r.pitch = a.pitch * (1.0 - u) + b.pitch * u;
    r.yaw = a.yaw * (1.0 - u) + b.yaw * u;
    r.dispNed = a.dispNed * (1.0 - u) + b.dispNed * u;
    return r;
}

// 读取POSMIND软件 GNSS/INS 参考结果，提取高程、速度和姿态信息。
std::vector<ReferenceState> readReferenceStates(const std::string &path) {
    std::ifstream in(path);
    std::vector<ReferenceState> refs;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        std::vector<std::string> f;
        std::string item;
        while (iss >> item) f.push_back(item);
        if (f.size() < 17 || f[0].find('/') == std::string::npos) continue;
        try {
            ReferenceState r;
            r.tow = std::stod(f[2]);
            r.h = std::stod(f[7]);
            r.velNed = {std::stod(f[9]), std::stod(f[8]), -std::stod(f[10])};
            r.roll = parseDegMin(std::stod(f[11]), std::stod(f[12])) * kD2R;
            r.pitch = parseDegMin(std::stod(f[13]), std::stod(f[14])) * kD2R;
            r.yaw = parseDegMin(std::stod(f[15]), std::stod(f[16])) * kD2R;
            if (!refs.empty()) {
                while (r.yaw - refs.back().yaw > kPi) r.yaw -= 2.0 * kPi;
                while (r.yaw - refs.back().yaw < -kPi) r.yaw += 2.0 * kPi;
            }
            refs.push_back(r);
        } catch (...) {
        }
    }
    for (size_t i = 1; i < refs.size(); ++i) {
        const double dt = refs[i].tow - refs[i - 1].tow;
        refs[i].dispNed = refs[i - 1].dispNed + (refs[i - 1].velNed + refs[i].velNed) * 0.5 * dt;
    }
    return refs;
}

// 读取自研 RTK 输出文件，供 post 模式或组合导航重新使用。
std::vector<GnssSolution> readGnssSolutions(const std::string &path) {
    std::ifstream in(path);
    std::vector<GnssSolution> sols;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        GnssSolution g;
        int week = 0;
        if (!(iss >> week >> g.time.tow >> g.llh.lat >> g.llh.lon >> g.llh.h
                  >> g.baseline.x >> g.baseline.y >> g.baseline.z
                  >> g.velNed.x >> g.velNed.y >> g.velNed.z
                  >> g.q >> g.ns >> g.ratio >> g.stdN >> g.stdE >> g.stdD)) {
            continue;
        }
        g.time.week = week;
        g.llh.lat *= kD2R;
        g.llh.lon *= kD2R;
        g.ecef = llhToEcef(g.llh);
        std::getline(iss, g.mode);
        g.mode = trim(g.mode);
        sols.push_back(g);
    }
    return sols;
}

// 读取组合导航输出文件，供 post 模式进行后处理平滑。
std::vector<FusionOutput> readFusionOutputs(const std::string &path) {
    std::ifstream in(path);
    std::vector<FusionOutput> outs;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        FusionOutput o;
        if (!(iss >> o.state.time.tow
                  >> o.state.ecef.x >> o.state.ecef.y >> o.state.ecef.z
                  >> o.state.llh.lat >> o.state.llh.lon >> o.state.llh.h
                  >> o.state.velNed.x >> o.state.velNed.y >> o.state.velNed.z
                  >> o.gnssQuality >> o.mode)) {
            continue;
        }
        o.state.llh.lat *= kD2R;
        o.state.llh.lon *= kD2R;
        o.state.ecef = llhToEcef(o.state.llh);
        outs.push_back(o);
    }
    return outs;
}

// 保存粗对准结果和参考姿态，方便报告中列出 roll/pitch/yaw 对比。
void writeCoarseAlignmentResult(const std::string &path, double alignSeconds,
                                double imuStartTow, double imuEndTow,
                                const GnssSolution &initGnss, const Vec3 &coarseEuler,
                                const Vec3 &filterEuler, const Vec3 &gyroBias,
                                bool hasInitialRef, const ReferenceState &initialRef) {
    std::ofstream out(path);
    out << "# coarse_alignment_result\n";
    out << std::fixed << std::setprecision(10);
    out << "align_seconds " << alignSeconds << "\n";
    out << "imu_start_tow " << imuStartTow << "\n";
    out << "imu_end_tow " << imuEndTow << "\n";
    out << "init_llh_deg " << initGnss.llh.lat * kR2D << ' ' << initGnss.llh.lon * kR2D << ' ' << initGnss.llh.h << "\n";
    out << "coarse_euler_deg " << coarseEuler.x * kR2D << ' ' << coarseEuler.y * kR2D << ' ' << coarseEuler.z * kR2D << "\n";
    out << "filter_euler_deg " << filterEuler.x * kR2D << ' ' << filterEuler.y * kR2D << ' ' << filterEuler.z * kR2D << "\n";
    out << "gyro_bias " << gyroBias.x << ' ' << gyroBias.y << ' ' << gyroBias.z << "\n";
    out << "has_initial_ref " << (hasInitialRef ? 1 : 0) << "\n";
    if (hasInitialRef) {
        out << "initial_ref " << initialRef.tow << ' '
            << initialRef.h << ' '
            << initialRef.velNed.x << ' ' << initialRef.velNed.y << ' ' << initialRef.velNed.z << ' '
            << initialRef.roll * kR2D << ' ' << initialRef.pitch * kR2D << ' ' << initialRef.yaw * kR2D << "\n";
    }
}

// 生成参考辅助诊断结果，只用于检查数据对齐和可视化
std::vector<FusionOutput> buildReferenceAidedSmooth(const std::vector<ImuSample> &imu,
                                                    const std::vector<GnssSolution> &gnss,
                                                    const std::vector<ReferenceState> &refs) {
    std::vector<FusionOutput> out;
    if (imu.empty() || gnss.empty() || refs.empty()) return out;
    const size_t n = std::min(imu.size(), gnss.size());
    out.reserve(n);
    size_t refIdx = 0;
    for (size_t i = 0; i < n; ++i) {
        const ReferenceState ref = interpolateReference(refs, gnss[i].time.tow, refIdx);
        FusionOutput o;
        o.state.time = gnss[i].time;
        o.state.llh = gnss[i].llh;
        o.state.ecef = gnss[i].ecef;
        o.state.velNed = ref.velNed;
        o.state.qbn = eulerToQuat(ref.roll, ref.pitch, ref.yaw);
        o.gnssQuality = gnss[i].q;
        o.mode = "reference_aided";
        out.push_back(o);
    }
    return out;
}

// 保存伪在线模式每秒播发一次的状态序列，便于做实时轨迹动画或现场展示。
void writeOneHzOnlineStream(const std::string &path, const std::vector<FusionOutput> &outputs) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write online stream: " + path);
    out << "# seq tow lat_deg lon_deg h_m vn ve vd roll pitch yaw q mode\n";
    out << std::fixed << std::setprecision(9);
    for (size_t i = 0; i < outputs.size(); ++i) {
        const FusionOutput &o = outputs[i];
        const Vec3 e = quatToEuler(o.state.qbn);
        out << i << ' '
            << o.state.time.tow << ' '
            << o.state.llh.lat * kR2D << ' '
            << o.state.llh.lon * kR2D << ' '
            << o.state.llh.h << ' '
            << o.state.velNed.x << ' '
            << o.state.velNed.y << ' '
            << o.state.velNed.z << ' '
            << e.x * kR2D << ' '
            << e.y * kR2D << ' '
            << e.z * kR2D << ' '
            << o.gnssQuality << ' '
            << o.mode << '\n';
    }
}

// 在控制台打印一条伪在线输出，模拟在线程序每秒发布当前导航状态。
void printOnlineEpoch(size_t seq, const FusionOutput &o) {
    const Vec3 e = quatToEuler(o.state.qbn);
    std::cout << std::setw(5) << seq
              << " tow=" << std::fixed << std::setprecision(2) << o.state.time.tow
              << " llh=(" << std::setprecision(8) << o.state.llh.lat * kR2D
              << "," << o.state.llh.lon * kR2D
              << "," << std::setprecision(3) << o.state.llh.h << ")"
              << " velNED=(" << o.state.velNed.x
              << "," << o.state.velNed.y
              << "," << o.state.velNed.z << ")"
              << " att=(" << e.x * kR2D
              << "," << e.y * kR2D
              << "," << e.z * kR2D << ")"
              << " q=" << o.gnssQuality
              << " " << o.mode << "\n";
}

// 伪在线处理：按 IMU/GNSS 时间戳顺序喂入滤波器，不提前使用未来数据。
OnlineRunResult runPseudoOnlineFusion(const std::vector<ImuSample> &imu,
                                      const std::vector<GnssSolution> &gnss,
                                      const FusionOptions &opt,
                                      const NavState &init,
                                      double firstTow,
                                      double lastTow,
                                      bool realTimeSleep) {
    // 伪实时核心：
    // 这里不调用 runLooseFusion，因为 runLooseFusion 后面会做后向平滑；
    // 本函数只按时间顺序处理已经到达的 IMU/GNSS 数据，保证不使用未来信息。
    // =========================
    // 伪实时在线处理函数
    // =========================
    // 用已有文件模拟传感器数据流，而不是把已解算结果拿来回放。
    // 本函数只按时间顺序消费 IMU/GNSS 数据：
    //   IMU 到达 -> ESKF predict；
    //   GNSS 历元到达 -> ESKF update；
    //   跨过整数秒 -> 打印当前最新状态。
    // 因为没有后向循环，也没有调用平滑函数，所以不会使用未来数据。
    OnlineRunResult result;
    result.allOutputs.reserve(imu.size());

    EsKf kf(opt, init);
    size_t ig = 0;
    int lastQuality = 0;
    std::string lastMode = "INS";
    double nextConstraintTow = init.time.tow;
    double nextPrintTow = std::ceil(firstTow);
    auto lastWallPrint = std::chrono::steady_clock::now();

    while (ig < gnss.size() && gnss[ig].time.tow < init.time.tow - 1.0e-6) ++ig;

    std::cout << "\nPseudo-online GNSS/INS sensor-stream mode started.\n";
    std::cout << "  source=offline IMU/GNSS files as simulated sensors\n";
    std::cout << "  IMU predict=200Hz by timestamp, GNSS update=1Hz when epoch arrives\n";
    std::cout << "  no future data, no backward/post smoothing\n";
    std::cout << "  wall_clock=" << (realTimeSleep ? "1s/printed epoch" : "fast") << "\n";
    std::cout << "seq,tow,lat_deg,lon_deg,h_m,vn,ve,vd,roll,pitch,yaw,q,mode\n";

    for (const auto &s : imu) {
        // IMU 是高频传感器流。每来一帧 IMU，就做一次 ESKF 预测。
        // Step 1: IMU 数据流到达。
        // 每读取到一个 IMU 历元，就立即做一次 ESKF 预测，不等待未来数据。
        if (s.time < init.time.tow || s.time + 1.0e-9 < firstTow) continue;
        if (s.time > lastTow + 0.01) break;

        // Step 1: IMU 高频预测。
        // predict() 内部完成姿态、速度、位置机械编排，同时传播 21 维误差协方差。
        // Step 2: 惯导预测。
        // predict() 内部完成机械编排和 21 维误差状态协方差传播。
        kf.predict(s);
        ++result.imuPackets;

        if (opt.useVehicleConstraint && s.time + 1.0e-9 >= nextConstraintTow) {
            // NHC 约束按固定时间间隔触发，不依赖未来历元。
            // Step 2: NHC 车辆约束。
            // 车辆正常行驶时车体系侧向和垂向速度接近 0，用该约束抑制横向漂移。
            // Step 3: NHC 非完整性约束更新。
            // 车辆侧向/垂向速度近似为零，用于抑制低速和弱 GNSS 段漂移。
            kf.updateVehicleConstraint();
            nextConstraintTow = s.time + std::max(0.005, opt.vehicleConstraintInterval);
        }

        while (ig < gnss.size() && gnss[ig].time.tow <= s.time + 0.002) {
            // GNSS/RTK 到达当前时间附近时才更新滤波器。
            // 这一步对应“程序收到 GNSS 数据包后做 Kalman 观测更新”。
            if (gnss[ig].time.tow >= firstTow - 1.0e-6 && gnss[ig].time.tow <= lastTow + 0.01) {
                // Step 3: GNSS/RTK 观测更新。
                // 只有当 GNSS 历元已经“到达”当前 IMU 时间，才允许进行 Kalman 更新。
                // 因此这里不会提前使用未来 GNSS 观测。
                // Step 4: GNSS/RTK 观测更新。
                // 只有 GNSS 历元已经到达当前 IMU 时间，才用位置/速度观测修正滤波器。
                kf.updatePosVel(gnss[ig]);
                lastQuality = gnss[ig].q;
                lastMode = gnss[ig].mode;
                ++result.gnssPackets;
            }
            ++ig;
        }

        // Step 4: 保存当前在线状态。
        // allOutputs 保存每个 IMU 历元的高频状态；oneHzOutputs 只保存每秒播发状态。
        // Step 5: 保存当前在线状态。
        // allOutputs 保存高频状态，oneHzOutputs 只保存模拟播发的一秒一次结果。
        FusionOutput out;
        out.state = kf.state();
        out.P = kf.P();
        out.gnssQuality = lastQuality;
        out.mode = lastMode.empty() ? "INS" : lastMode;
        result.allOutputs.push_back(out);

        const double tow = out.state.time.tow;
        if (tow + 1.0e-6 >= nextPrintTow) {
            // 每跨过一个整数秒，就打印当前最新状态，模拟在线系统每秒播发一次结果。
            // Step 5: 1 Hz 结果播发。
            // pseudo-online 会按墙钟等待；pseudo-online-fast 跳过等待，仅用于快速检查。
            // Step 6: 1 Hz 结果播发。
            // pseudo-online 等待真实墙钟，pseudo-online-fast 只快速打印不等待。
            result.oneHzOutputs.push_back(out);
            printOnlineEpoch(result.oneHzOutputs.size() - 1, out);
            std::cout.flush();
            if (realTimeSleep && result.oneHzOutputs.size() > 1) {
                const auto target = lastWallPrint + std::chrono::seconds(1);
                std::this_thread::sleep_until(target);
                lastWallPrint = target;
            } else {
                lastWallPrint = std::chrono::steady_clock::now();
            }
            nextPrintTow += 1.0;
            while (nextPrintTow < tow) nextPrintTow += 1.0;
        }
    }

    std::cout << "Pseudo-online GNSS/INS mode finished: imu_packets=" << result.imuPackets
              << " gnss_updates=" << result.gnssPackets
              << " printed_1Hz=" << result.oneHzOutputs.size() << "\n\n";
    return result;
}

} // namespace

// 主函数：根据 mode 调度 RTK、前向松组合、后处理平滑和伪在线处理。
int main(int argc, char **argv) {
    try {
        AppConfig cfg;
        // 主流程的模式接口就是 cfg.mode。现场可以直接修改 AppConfig 中的默认值，
        // 然后在 VS 中重新生成并运行。
        if (argc > 1) cfg.outDir = argv[1];
        if (argc > 2) {
            if (isIntegerText(argv[2])) cfg.maxEpochs = std::stoi(argv[2]);
            else cfg.mode = normalizeMode(argv[2]);
        }
        if (argc > 3) cfg.mode = normalizeMode(argv[3]);
        ensureDir(cfg.outDir);

        std::cout << "GNSS/INS self-developed C++ processor\n";
        std::cout << "Mode: " << cfg.mode << "\n";
        std::cout << "Modes: rtk | forward | post | all | pseudo-online | pseudo-online-fast\n";

        const bool pseudoOnline = cfg.mode == "pseudo-online" || cfg.mode == "pseudo-online-fast";
        const bool needRtk = cfg.mode == "rtk" || cfg.mode == "forward" || cfg.mode == "post" ||
                             cfg.mode == "all" || pseudoOnline;
        const bool needFusion = cfg.mode == "forward" || cfg.mode == "post" || cfg.mode == "all" ||
                                pseudoOnline;

        RinexObsFile rover;
        RinexObsFile base;
        NavData nav;
        if (needRtk) {
            std::cout << "Reading RINEX observation/navigation files...\n";
            rover = readRinexObs(cfg.rover, cfg.maxEpochs);
            base = readRinexObs(cfg.base, 25000);
            nav = readRinexNav(cfg.nav);
            std::cout << "Rover epochs: " << rover.epochs.size()
                      << ", base epochs: " << base.epochs.size()
                      << ", ephemerides: " << nav.ephs.size() << "\n";
        }

        std::vector<GnssSolution> gnss;
        std::vector<GnssSolution> fusionGnss;
        if (needRtk) {
            std::cout << "Running self-developed GPS/BDS Kalman-filter RTK...\n";
            RtkOptions rtkOpt;
            gnss = runRtkFloat(rover, base, nav, rtkOpt);
            if (gnss.empty()) throw std::runtime_error("GNSS solution is empty");
            writeGnssSolutions(cfg.outDir + "/gnss_rtk.pos", gnss);
            fusionGnss = gnss;
            std::cout << "GNSS solutions: " << gnss.size() << "\n";
            int nFix = 0, nFloat = 0, nFallback = 0;
            for (const auto &s : gnss) {
                if (s.q == 1) ++nFix;
                else if (s.q == 2) ++nFloat;
                else ++nFallback;
            }
            std::cout << "  RTK_FIX=" << nFix << " RTK_FLOAT=" << nFloat << " FALLBACK=" << nFallback << "\n";
            if (cfg.mode == "rtk") {
                std::cout << "RTK-only mode finished. Output folder: " << cfg.outDir << "\n";
                return 0;
            }
        }

        if (needFusion) {
            std::cout << "Decoding IMR and writing text increments...\n";
            std::vector<ImuSample> imu = decodeImr(cfg.imu, cfg.outDir + "/imu_decoded.txt");
            if (imu.empty()) throw std::runtime_error("IMU decoded data is empty");
            if (fusionGnss.empty()) fusionGnss = readGnssSolutions(cfg.outDir + "/gnss_rtk.pos");
            if (fusionGnss.empty()) throw std::runtime_error("GNSS input for fusion is empty");

            double lastUsableGnssTow = fusionGnss.back().time.tow;
            for (auto it = fusionGnss.rbegin(); it != fusionGnss.rend(); ++it) {
                if (it->q > 0 && it->q != 5) {
                    lastUsableGnssTow = it->time.tow;
                    break;
                }
            }
            const double keepUntil = std::max(imu.front().time + cfg.coarseAlignSeconds + 1.0,
                                              lastUsableGnssTow + 0.01);
            imu.erase(std::remove_if(imu.begin(), imu.end(),
                                     [keepUntil](const ImuSample &s) { return s.time > keepUntil; }),
                      imu.end());
            std::cout << "IMU samples: " << imu.size()
                      << ", first tow: " << imu.front().time
                      << ", last tow: " << imu.back().time << "\n";

            std::cout << "Static coarse alignment (first " << cfg.coarseAlignSeconds << " s)...\n";
            const GnssSolution *initGnss = &fusionGnss.front();
            for (const auto &g : fusionGnss) {
                if (g.q == 1) {
                    initGnss = &g;
                    break;
                }
            }
            NavState init = coarseAlign(imu, *initGnss, cfg.coarseAlignSeconds, cfg.defaultYawDeg);
            Vec3 coarseEuler = quatToEuler(init.qbn);
            std::vector<ReferenceState> refs = readReferenceStates(cfg.referenceGnssIns);
            ReferenceState initRef;
            bool haveInitRef = false;
            if (refs.size() >= 2 && initGnss->time.tow >= refs.front().tow && initGnss->time.tow <= refs.back().tow) {
                size_t refIdx = 0;
                initRef = interpolateReference(refs, initGnss->time.tow, refIdx);
                haveInitRef = true;
            }
            Vec3 eulerInit = quatToEuler(init.qbn);
            std::cout << "  coarse roll=" << coarseEuler.x * kR2D << " pitch=" << coarseEuler.y * kR2D
                      << " yaw=" << coarseEuler.z * kR2D << " deg\n";
            std::cout << "  filter init roll=" << eulerInit.x * kR2D << " pitch=" << eulerInit.y * kR2D
                      << " yaw=" << eulerInit.z * kR2D
                      << " deg (coarse alignment + field yaw)\n";
            writeCoarseAlignmentResult(cfg.outDir + "/coarse_alignment_result.txt",
                                        cfg.coarseAlignSeconds, imu.front().time, imu.back().time,
                                        *initGnss, coarseEuler, eulerInit, init.gyroBias,
                                        haveInitRef, initRef);

            double firstUsableGnssTow = fusionGnss.front().time.tow;
            lastUsableGnssTow = fusionGnss.back().time.tow;
            for (const auto &g : fusionGnss) {
                if (g.q == 1) {
                    firstUsableGnssTow = g.time.tow;
                    break;
                }
            }
            for (auto it = fusionGnss.rbegin(); it != fusionGnss.rend(); ++it) {
                if (it->q > 0 && it->q != 5) {
                    lastUsableGnssTow = it->time.tow;
                    break;
                }
            }

            FusionOptions fusionOpt;
            if (pseudoOnline) {
                std::cout << "Running true pseudo-online ESKF sensor-stream processing...\n";
                OnlineRunResult online = runPseudoOnlineFusion(imu, fusionGnss, fusionOpt, init,
                                                               firstUsableGnssTow, lastUsableGnssTow,
                                                               cfg.mode == "pseudo-online");
                writeFusionOutputs(cfg.outDir + "/gnss_ins_pseudo_online.nav", online.allOutputs);
                writeFusionOutputs(cfg.outDir + "/gnss_ins_forward.nav", online.allOutputs);
                writeOneHzOnlineStream(cfg.outDir + "/gnss_ins_pseudo_online_stream.txt", online.oneHzOutputs);
                std::cout << "Pseudo-online mode finished. Output folder: " << cfg.outDir << "\n";
                return 0;
            }

            std::cout << "Running forward ESKF...\n";
            std::vector<FusionOutput> fused = runLooseFusion(imu, fusionGnss, fusionOpt, init);
            fused.erase(std::remove_if(fused.begin(), fused.end(),
                                      [firstUsableGnssTow, lastUsableGnssTow](const FusionOutput &o) {
                                          return o.state.time.tow < firstUsableGnssTow ||
                                                 o.state.time.tow > lastUsableGnssTow + 0.01;
                                      }),
                        fused.end());
            writeFusionOutputs(cfg.outDir + "/gnss_ins_loose.nav", fused);
            writeFusionOutputs(cfg.outDir + "/gnss_ins_forward.nav", fused);
            std::cout << "Fusion outputs: " << fused.size() << "\n";

            if (cfg.mode == "forward") {
                std::cout << "Forward-only mode finished. Output folder: " << cfg.outDir << "\n";
                return 0;
            }

            std::vector<FusionOutput> forwardPost;
            std::vector<FusionOutput> backwardPost;
            if (cfg.mode == "post" || cfg.mode == "all") {
                std::vector<FusionOutput> forwardDisk = cfg.mode == "post" ? readFusionOutputs(cfg.outDir + "/gnss_ins_forward.nav") : fused;
                std::vector<GnssSolution> gnssDisk = cfg.mode == "post" ? readGnssSolutions(cfg.outDir + "/gnss_rtk.pos") : fusionGnss;
                forwardPost = rtsLooseFusionPost(smoothLooseFusionPost(forwardDisk, gnssDisk), imu);
                backwardPost = backwardLooseFusionPost(forwardDisk, gnssDisk);
            }
            if (!forwardPost.empty()) {
                writeFusionOutputs(cfg.outDir + "/gnss_ins_forward_smooth.nav", forwardPost);
                std::cout << "Forward post-smoothed outputs: " << forwardPost.size() << "\n";
            }
            if (!backwardPost.empty()) {
                writeFusionOutputs(cfg.outDir + "/gnss_ins_backward.nav", backwardPost);
                std::cout << "Backward post-smoothed outputs: " << backwardPost.size() << "\n";
            }
            if (forwardPost.size() == backwardPost.size() && !forwardPost.empty()) {
                std::vector<FusionOutput> fbSmooth = forwardPost;
                for (size_t i = 0; i < fbSmooth.size(); ++i) {
                    fbSmooth[i].state.ecef = (forwardPost[i].state.ecef + backwardPost[i].state.ecef) * 0.5;
                    fbSmooth[i].state.llh = ecefToLlh(fbSmooth[i].state.ecef);
                    fbSmooth[i].state.velNed = (forwardPost[i].state.velNed + backwardPost[i].state.velNed) * 0.5;
                    fbSmooth[i].state.qbn = averageAttitude(forwardPost[i].state.qbn, backwardPost[i].state.qbn);
                    if (fbSmooth[i].P.rows() == backwardPost[i].P.rows() &&
                        fbSmooth[i].P.cols() == backwardPost[i].P.cols()) {
                        for (int r = 0; r < fbSmooth[i].P.rows(); ++r) {
                            for (int c = 0; c < fbSmooth[i].P.cols(); ++c) {
                                fbSmooth[i].P(r, c) = 0.50 * (forwardPost[i].P(r, c) + backwardPost[i].P(r, c));
                            }
                        }
                    }
                    fbSmooth[i].mode = forwardPost[i].mode + "_FB";
                }
                writeFusionOutputs(cfg.outDir + "/gnss_ins_smooth.nav", fbSmooth);
                std::cout << "Forward/backward smoothed outputs: " << fbSmooth.size() << "\n";
            }

            std::vector<FusionOutput> refAided = buildReferenceAidedSmooth(imu, fusionGnss, refs);
            if (!refAided.empty()) {
                writeFusionOutputs(cfg.outDir + "/gnss_ins_reference_aided.nav", refAided);
                std::cout << "Reference-aided diagnostic outputs: " << refAided.size() << "\n";
            }
        }

        std::cout << "Finished. Output folder: " << cfg.outDir << "\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
