#pragma once

#include "common.h"

namespace gnssins {

// GNSS 模块负责从原始观测数据得到可供组合导航使用的 GNSS 解：
// 1. 读取流动站、基准站 RINEX 观测文件和广播星历；
// 2. 根据广播星历计算卫星位置和钟差，完成单点定位初值；
// 3. 构造站间单差、站星双差模型，进行 RTK Kalman 滤波与模糊度固定；
// 4. 输出约 1 Hz 的 GnssSolution，作为松组合 ESKF 的位置、速度观测。
enum class System { GPS, BDS, GLO, OTHER };

struct SatId {
    // 卫星编号由系统类型和 PRN 组成，用作 map 的键。
    System sys = System::OTHER;
    int prn = 0;

    bool operator<(const SatId &o) const {
        if (sys != o.sys) return static_cast<int>(sys) < static_cast<int>(o.sys);
        return prn < o.prn;
    }
};

// 将卫星编号格式化为 G01、C19 等字符串，便于输出和调试。
std::string satToString(const SatId &sat);

struct ObsValue {
    // 一个卫星在一个历元上的观测值。
    // 数组 0/1/2 表示不同频点，has* 标志说明该频点观测是否存在。
    double code[3] = {};
    double phase[3] = {};
    double doppler[3] = {};
    double snr[3] = {};
    int lli[3] = {};
    bool hasCode[3] = {};
    bool hasPhase[3] = {};
    bool hasDoppler[3] = {};
};

struct ObsEpoch {
    // 一个观测历元，包含该时刻所有可见卫星的观测量。
    GTime time;
    std::map<SatId, ObsValue> obs;
};

struct RinexObsFile {
    // RINEX 观测文件解析后的结果：近似坐标、观测类型表和所有历元。
    Vec3 approxPos;
    std::map<System, std::vector<std::string>> obsTypes;
    std::vector<ObsEpoch> epochs;
};

struct Eph {
    // 广播星历参数。GPS/BDS 使用开普勒轨道参数，GLO 字段保留给 GLONASS。
    SatId sat;
    GTime toc;
    double af0 = 0.0, af1 = 0.0, af2 = 0.0;
    double iode = 0.0, crs = 0.0, dn = 0.0, m0 = 0.0;
    double cuc = 0.0, e = 0.0, cus = 0.0, sqrtA = 0.0;
    double toe = 0.0, cic = 0.0, omg0 = 0.0, cis = 0.0;
    double i0 = 0.0, crc = 0.0, omg = 0.0, omgDot = 0.0;
    double idot = 0.0;
    int week = 0;
    double toes = 0.0;
    GTime toeTime;
    Vec3 gloPosKm;
    Vec3 gloVelKmS;
    Vec3 gloAccKmS2;
};

struct NavData {
    // 广播星历集合，satelliteState 会按卫星号和时间检索最合适的一组。
    std::vector<Eph> ephs;
};

struct SatState {
    // 卫星位置和钟差结果，供 SPP/RTK 构造观测方程。
    Vec3 pos;
    double clk = 0.0;
    bool valid = false;
};

struct GnssSolution {
    // GNSS/RTK 解算结果，是 GNSS 模块传给组合导航模块的主要接口。
    GTime time;
    Vec3 ecef;      // 天线 ECEF 坐标 [m]
    Vec3 baseline;  // 相对基站的基线向量 [m]
    Vec3 velNed;    // NED 速度 [m/s]
    Llh llh;        // 天线经纬高
    double clock = 0.0;
    int ns = 0;
    int q = 0;      // 1=fixed, 2=float, 5=SPP/fallback, 0=invalid
    double ratio = 0.0;
    double stdN = 5.0;
    double stdE = 5.0;
    double stdD = 10.0;
    std::string mode;
};

struct RtkOptions {
    // RTK 参数集中放在这里，便于说明观测定权、模糊度固定阈值和系统开关。
    double elevationMaskDeg = 10.0;
    double codeSigma = 2.0;
    double phaseSigma = 0.02;
    double processPos = 0.5;
    double processAccel = 0.8;
    double lambdaRatioThreshold = 3.0;
    double fixedPhaseRmsMax = 0.08;
    double fixedPhaseAbsMax = 0.20;
    int minAmbLock = 5;
    bool useGps = true;
    bool useBds = true;
    bool useGlo = false;
    bool useBdsGeo = true;
    bool useGfMwSlip = true;
    bool allowPartialAr = true;
    bool allowSystemAr = true;
    bool allowCandidateOverride = true;
};

// 读取 RINEX 观测文件，解析头文件观测类型和逐历元观测值。
RinexObsFile readRinexObs(const std::string &path, int maxEpochs = -1);
// 读取 RINEX 广播星历文件，为卫星位置和钟差计算提供轨道参数。
NavData readRinexNav(const std::string &path);
// 根据广播星历计算指定卫星在指定时刻的 ECEF 位置和卫星钟差。
bool satelliteState(const NavData &nav, const SatId &sat, const GTime &t, SatState &state);
// 单点定位，主要为 RTK 提供流动站初始坐标，也作为失锁时的保底解。
bool sppSolve(const ObsEpoch &epoch, const NavData &nav, Vec3 initial, GnssSolution &sol);
// 执行双系统双频 RTK Kalman 滤波，输出 float/fix 解序列。
std::vector<GnssSolution> runRtkFloat(const RinexObsFile &rover, const RinexObsFile &base, const NavData &nav,
                                      const RtkOptions &opt);
// 将 GNSS/RTK 结果写成 CSV，供组合导航和可视化程序读取。
void writeGnssSolutions(const std::string &path, const std::vector<GnssSolution> &sols);

} // namespace gnssins
