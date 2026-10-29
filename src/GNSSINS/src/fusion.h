#pragma once

#include "ins.h"

namespace gnssins {

// 松组合模块负责把 RTK 结果和 IMU 递推结果融合起来。
// forward 是按历元顺序跑的在线风格结果，post 是离线平滑增强结果。
struct FusionOptions {
    // 松组合配置：
    // leverArm 用于把 IMU 中心位置改正到 GNSS 天线位置；
    // useVehicleConstraint 控制是否加入 NHC 非完整性约束；
    // gnssRScale/nhcRScale 用于超参数敏感性实验。
    Vec3 leverArm{0.195, -0.185, 0.000};  // measured IMU-to-antenna lever arm [m]
    Vec3 mountAngle{0, 0, 0};             // not used, kept for config
    ImuModel imuModel;
    bool useVehicleConstraint = true;
    double vehicleConstraintInterval = 0.05;
    bool strongGnssClosure = true;
    double gnssRScale = 1.0;  // multiplier for GNSS measurement covariance R
    double nhcRScale = 1.0;   // multiplier for NHC measurement covariance R
};

struct FusionOutput {
    // 每个 IMU 历元输出一次导航状态和协方差，后续可用于画轨迹、误差和 STD 曲线。
    NavState state;
    Matrix P;  // error-state covariance saved for STD analysis
    int gnssQuality = 0;
    std::string mode;
};

struct EsKfState {
    NavState nominal;          // full-state (closed-loop nominal)
    Matrix P;                  // 15x15 error covariance
};

class EsKf {
public:
    explicit EsKf(const FusionOptions &opt, const NavState &init);
    // predict: IMU 到一帧就传播一次，是在线伪实时的 200 Hz 主循环。
    void predict(const ImuSample &imu);
    // updatePosVel: GNSS/RTK 解到达时才调用，构造位置/速度残差完成 Kalman 更新。
    void updatePosVel(const GnssSolution &gnss);
    // updateVehicleConstraint: NHC 约束，抑制车辆侧向/垂向速度误差。
    void updateVehicleConstraint();
    void correctNominal();
    const NavState &state() const { return st_.nominal; }
    const Matrix &P() const { return st_.P; }
    void setState(const NavState &s);
    void setP(const Matrix &p);

private:
    FusionOptions opt_;
    EsKfState st_;
    Vec3 prevDtheta_;
    Vec3 prevDvel_;
};

// RTS backward smoother / post smoother 都是离线后处理函数。
// 注意：pseudo-online 模式不调用这些函数，因此不提前使用未来数据。
void rtsSmooth(std::vector<NavState> &states, std::vector<Matrix> &covs,
               const std::vector<ImuSample> &imu);

std::vector<FusionOutput> runLooseFusion(const std::vector<ImuSample> &imu,
                                         const std::vector<GnssSolution> &gnss,
                                         const FusionOptions &options, const NavState &initial);
std::vector<FusionOutput> rtsLooseFusionPost(const std::vector<FusionOutput> &loose,
                                             const std::vector<ImuSample> &imu);
std::vector<FusionOutput> backwardLooseFusionPost(const std::vector<FusionOutput> &loose,
                                                  const std::vector<GnssSolution> &gnss);
std::vector<FusionOutput> smoothLooseFusionPost(const std::vector<FusionOutput> &loose,
                                                const std::vector<GnssSolution> &gnss);
Quat averageAttitude(const Quat &a, const Quat &b);
void writeFusionOutputs(const std::string &path, const std::vector<FusionOutput> &outputs);

} // namespace gnssins

