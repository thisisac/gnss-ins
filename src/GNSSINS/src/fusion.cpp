#include "fusion.h"

namespace gnssins {

// fusion.cpp
// ---------------------------------------------------------------------------
// 本文件实现 GNSS/INS 松组合算法。
//
// 主要内容：
// 1. EsKf 类：21 维误差状态卡尔曼滤波器；
//    状态顺序为位置误差、速度误差、姿态误差、加表零偏、陀螺零偏、
//    加表比例因子、陀螺比例因子。
// 2. predict：IMU 每来一帧就进行一次名义状态机械编排和协方差预测。
// 3. updatePosVel：GNSS/RTK 历元到达时，用位置/速度残差进行量测更新。
// 4. updateVehicleConstraint：NHC 非完整性约束，抑制车辆侧向和垂向速度。
// 5. smoothLooseFusionPost / backwardLooseFusionPost：离线后处理平滑。
//
// 注意：
// pseudo-online 伪实时模式只直接调用 EsKf::predict / updatePosVel /
// updateVehicleConstraint，不调用离线平滑函数，因此不使用未来数据。
// ---------------------------------------------------------------------------
static constexpr int kDim = 21;
static constexpr int kPos = 0, kVel = 3, kAtt = 6, kAcc = 9, kGyr = 12, kAccScale = 15, kGyrScale = 18;

// 将 NED 坐标差转换为 ECEF 坐标差，常用于把杆臂或平滑改正量加回 ECEF 坐标。
static Vec3 nedToEcef(const Vec3 &dn, const Llh &ref) {
    return matVec(transpose(ecefToNedMatrix(ref)), dn);
}

// 将角度归一到 [-pi, pi]，避免航向跨 180 度时出现跳变。
static double wrapPi(double a) {
    while (a > kPi) a -= 2.0 * kPi;
    while (a < -kPi) a += 2.0 * kPi;
    return a;
}

// 按权重在两个角度之间插值，先处理跨 pi 的角度环绕。
static double blendAngle(double a, double b, double w) {
    return a + wrapPi(b - a) * w;
}

static void applyCourseAttitudeAid(std::vector<FusionOutput> &out, double courseWeight,
                                   double neighborWeight) {
    // 离线姿态辅助平滑：
    // 当车辆存在明显水平速度时，速度方向可近似提供航向观测。
    // 该函数只用于后处理展示，不参与 pseudo-online 在线结果。
    if (out.size() < 2) return;
    for (auto &o : out) {
        Vec3 e = quatToEuler(o.state.qbn);
        const double speedH = std::sqrt(sqr(o.state.velNed.x) + sqr(o.state.velNed.y));
        if (speedH > 0.35) {
            const double courseYaw = std::atan2(o.state.velNed.y, o.state.velNed.x);
            e.z = blendAngle(e.z, courseYaw, courseWeight);
        }
        o.state.qbn = eulerToQuat(e.x, e.y, e.z);
    }

    std::vector<Quat> original;
    original.reserve(out.size());
    for (const auto &o : out) original.push_back(o.state.qbn);
    for (size_t i = 1; i + 1 < out.size(); ++i) {
        Vec3 prev = quatToEuler(original[i - 1]);
        Vec3 now = quatToEuler(original[i]);
        Vec3 next = quatToEuler(original[i + 1]);
        Vec3 e = now;
        e.x = blendAngle(now.x, 0.5 * (prev.x + next.x), neighborWeight);
        e.y = blendAngle(now.y, 0.5 * (prev.y + next.y), neighborWeight);
        e.z = blendAngle(blendAngle(now.z, prev.z, neighborWeight * 0.5),
                         next.z, neighborWeight * 0.5);
        out[i].state.qbn = eulerToQuat(e.x, e.y, e.z);
    }
}

// 两个姿态四元数的简化平均，先转欧拉角再做角度插值。
Quat averageAttitude(const Quat &a, const Quat &b) {
    Vec3 ea = quatToEuler(a);
    Vec3 eb = quatToEuler(b);
    return eulerToQuat(blendAngle(ea.x, eb.x, 0.5),
                       blendAngle(ea.y, eb.y, 0.5),
                       blendAngle(ea.z, eb.z, 0.5));
}

// 判断一个 GNSS/RTK 历元是否适合进入松组合观测更新。
static bool usableGnssUpdate(const GnssSolution &g) {
    if (g.q == 1) {
        if (g.ns < 10) return false;
        if (g.stdN > 0.20 || g.stdE > 0.20 || g.stdD > 0.30) return false;
        if (norm(g.velNed) > 3.0) return false;
        return true;
    }
    if (g.q == 2) {
        if (g.ns < 5) return false;
        if (g.stdN > 6.0 || g.stdE > 6.0 || g.stdD > 12.0) return false;
        if (norm(g.velNed) > 8.0) return false;
        return true;
    }
    return false;
}

// 对协方差矩阵整体缩放，用于离线平滑结果的 STD 表示。
static void scaleCovariance(Matrix &p, double scale) {
    if (p.rows() == 0 || p.cols() == 0) return;
    for (int r = 0; r < p.rows(); ++r) {
        for (int c = 0; c < p.cols(); ++c) {
            p(r, c) *= scale;
        }
    }
}

struct FusionAnchor {
    double tow = 0.0;
    Vec3 posCorrNed;
    Vec3 velCorrNed;
    bool hasVertical = false;
    bool hasVelocity = false;
};

// 从前向松组合结果和 GNSS 结果中提取可用于后处理平滑的时间锚点。
static std::vector<FusionAnchor> collectFusionAnchors(const std::vector<FusionOutput> &loose,
                                                      const std::vector<GnssSolution> &gnss) {
    std::vector<FusionAnchor> anchors;
    anchors.reserve(gnss.size());
    size_t il = 0;
    for (const auto &g : gnss) {
        if (!usableGnssUpdate(g)) continue;
        while (il + 1 < loose.size() &&
               std::abs(loose[il + 1].state.time.tow - g.time.tow) <
                   std::abs(loose[il].state.time.tow - g.time.tow)) {
            ++il;
        }
        const NavState &s = loose[il].state;
        const double dt = std::abs(s.time.tow - g.time.tow);
        if (dt > 0.02) continue;

        FusionAnchor a;
        a.tow = g.time.tow;
        a.posCorrNed = ecefToNed(g.ecef - s.ecef, s.llh);
        a.velCorrNed = g.velNed - s.velNed;
        a.hasVertical = g.q == 1 || g.stdD <= 0.50;
        a.hasVelocity = norm(g.velNed) < (g.q == 1 ? 3.0 : 5.0);
        if (g.q != 1) {
            a.posCorrNed.z = 0.0;
            a.velCorrNed.z = 0.0;
        }
        anchors.push_back(a);
    }
    return anchors;
}

// 在相邻锚点之间插值得到当前时刻的位置/速度改正量。
static FusionAnchor interpolateAnchor(const std::vector<FusionAnchor> &anchors,
                                      size_t ia, double tow, double futureBias = 1.0) {
    FusionAnchor corr = anchors[ia];
    if (ia + 1 < anchors.size()) {
        const FusionAnchor &a = anchors[ia];
        const FusionAnchor &b = anchors[ia + 1];
        const double span = b.tow - a.tow;
        double u = span > 0.0 ? std::max(0.0, std::min(1.0, (tow - a.tow) / span)) : 0.0;
        if (futureBias != 1.0) {
            u = std::pow(u, futureBias);
        }
        corr.tow = tow;
        corr.posCorrNed = a.posCorrNed * (1.0 - u) + b.posCorrNed * u;
        corr.velCorrNed = a.velCorrNed * (1.0 - u) + b.velCorrNed * u;
        corr.hasVertical = a.hasVertical && b.hasVertical;
        corr.hasVelocity = a.hasVelocity && b.hasVelocity;
        if (!corr.hasVertical) corr.posCorrNed.z = 0.0;
        if (!corr.hasVelocity) corr.velCorrNed = {};
    }
    return corr;
}

// 将 21 维误差状态闭环改正到名义导航状态中。
static void applyErrorToState(NavState &s, const std::vector<double> &dx) {
    double s2 = std::sin(s.llh.lat); s2 *= s2;
    double Rn = kA / std::sqrt(1 - kE2 * s2);
    double Rm = Rn * (1 - kE2) / (1 - kE2 * s2);
    s.llh.lat += dx[kPos + 0] / (Rm + s.llh.h);
    s.llh.lon += dx[kPos + 1] / ((Rn + s.llh.h) * std::cos(s.llh.lat));
    s.llh.h -= dx[kPos + 2];
    s.velNed.x += dx[kVel + 0];
    s.velNed.y += dx[kVel + 1];
    s.velNed.z += dx[kVel + 2];

    Vec3 phi = {dx[kAtt + 0], dx[kAtt + 1], dx[kAtt + 2]};
    double phiMag = norm(phi);
    Quat dqCorr = phiMag > 1e-12 ? quatFromAxisAngle(phi / phiMag, -phiMag) : Quat{1, 0, 0, 0};
    s.qbn = quatNormalize(quatMul(dqCorr, s.qbn));

    s.accBias.x += dx[kAcc + 0]; s.accBias.y += dx[kAcc + 1]; s.accBias.z += dx[kAcc + 2];
    s.gyroBias.x += dx[kGyr + 0]; s.gyroBias.y += dx[kGyr + 1]; s.gyroBias.z += dx[kGyr + 2];
    s.accScale.x += dx[kAccScale + 0]; s.accScale.y += dx[kAccScale + 1]; s.accScale.z += dx[kAccScale + 2];
    s.gyroScale.x += dx[kGyrScale + 0]; s.gyroScale.y += dx[kGyrScale + 1]; s.gyroScale.z += dx[kGyrScale + 2];
    s.ecef = llhToEcef(s.llh);
}

EsKf::EsKf(const FusionOptions &opt, const NavState &init) : opt_(opt) {
    // ESKF 初始化：名义状态来自粗对准/GNSS 初值，P 矩阵来自 ImuModel 初始标准差。
    // 后续滤波只估计“小误差状态”，再闭环修正到 nominal 中。
    st_.nominal = init;
    st_.P = Matrix(kDim, kDim);
    const auto &m = opt.imuModel;
    for (int i = 0; i < 3; ++i) {
        st_.P(kPos + i, kPos + i) = sqr(m.initPosStd[i]);
        st_.P(kVel + i, kVel + i) = sqr(m.initVelStd[i]);
        st_.P(kAtt + i, kAtt + i) = sqr(m.initAttStd[i]);
        st_.P(kAcc + i, kAcc + i) = sqr(m.initAccBiasStd[i]);
        st_.P(kGyr + i, kGyr + i) = sqr(m.initGyroBiasStd[i]);
        st_.P(kAccScale + i, kAccScale + i) = sqr(m.initAccScaleStd[i]);
        st_.P(kGyrScale + i, kGyrScale + i) = sqr(m.initGyroScaleStd[i]);
    }
}

// 直接替换名义状态，主要供调试或外部重置使用。
void EsKf::setState(const NavState &s) { st_.nominal = s; }
// 直接替换协方差矩阵，主要供后处理或敏感性实验使用。
void EsKf::setP(const Matrix &p) { st_.P = p; }

// ESKF 预测步：用一帧 IMU 增量推进名义状态，并传播 21 维误差协方差。
void EsKf::predict(const ImuSample &imu) {
    // 预测步对应在线系统的 IMU 数据接收：
    // 每来一帧 IMU，就扣除当前估计的零偏/比例因子，更新名义姿态、速度、位置，
    // 同时用状态转移矩阵 Phi 和过程噪声 Qd 传播误差协方差 P。
    //
    // 该函数可以理解为两条并行链路：
    //   A. 名义状态机械编排：
    //      dtheta/dvel -> 姿态四元数 -> 比力转到导航系 -> 速度 -> 位置；
    //   B. 误差状态协方差传播：
    //      构造连续误差模型 F 的主要块 -> 离散化为 Phi -> P = Phi P Phi^T + Qd。
    // 在线伪实时模式中，这个函数约每 0.005 s 调用一次。
    double dt = imu.dt > 0 ? imu.dt : 0.005;

    // 前向传播：先把惯导名义状态按 IMU 增量推到下一历元。
    NavState &ns = st_.nominal;
    // 1. IMU 误差补偿：
    //    陀螺角增量扣除陀螺零偏和比例因子，加表速度增量同理。
    Vec3 dtheta = {imu.dtheta.x * (1.0 - ns.gyroScale.x) - ns.gyroBias.x * dt,
                   imu.dtheta.y * (1.0 - ns.gyroScale.y) - ns.gyroBias.y * dt,
                   imu.dtheta.z * (1.0 - ns.gyroScale.z) - ns.gyroBias.z * dt};
    Vec3 dvel   = {imu.dvel.x * (1.0 - ns.accScale.x) - ns.accBias.x * dt,
                   imu.dvel.y * (1.0 - ns.accScale.y) - ns.accBias.y * dt,
                   imu.dvel.z * (1.0 - ns.accScale.z) - ns.accBias.z * dt};

    // 2. 圆锥/划桨补偿：
    //    简化的一阶补偿，降低高频姿态和速度增量非交换误差。
    Vec3 dthetaCorr = dtheta + cross(prevDtheta_, dtheta) / 12.0;
    Vec3 dvelCorr   = dvel   + cross(prevDtheta_, dvel + prevDvel_ / 6.0) / 2.0;
    prevDtheta_ = dtheta; prevDvel_ = dvel;

    // 3. 姿态更新：
    //    将角增量转为小旋转四元数，再右乘更新 body-to-NED 姿态。
    double dmag = norm(dthetaCorr);
    Quat dq = dmag > 1e-12 ? quatFromAxisAngle(dthetaCorr / dmag, dmag) : Quat{1, 0, 0, 0};
    ns.qbn = quatNormalize(quatMul(ns.qbn, dq));

    Vec3 fbNow = dvel / dt;
    const double fbNorm = norm(fbNow);
    if (fbNorm > 7.0 && fbNorm < 12.0) {
        // Weak leveling aid: when measured specific force is close to gravity,
        // gently pull roll/pitch toward accelerometer tilt. Yaw is untouched
        // because gravity gives no heading information for a stationary MEMS IMU.
        Vec3 e = quatToEuler(ns.qbn);
        const double rollAcc = std::atan2(-fbNow.y, -fbNow.z);
        const double pitchAcc = std::atan2(fbNow.x, std::sqrt(fbNow.y * fbNow.y + fbNow.z * fbNow.z));
        const double a = 0.002;
        e.x = e.x + a * wrapPi(rollAcc - e.x);
        e.y = e.y + a * wrapPi(pitchAcc - e.y);
        ns.qbn = eulerToQuat(e.x, e.y, e.z);
    }

    // 4. 速度和位置更新：
    //    把机体系速度增量转到导航系，加入重力，扣除地球自转和运输率影响。
    Mat3 Cbn = quatToDcm(ns.qbn);
    Vec3 fn = matVec(Cbn, dvelCorr) / dt;
    Vec3 wie_n = earthRateNed(ns.llh.lat);
    Vec3 wen_n = transportRateNed(ns.llh.lat, ns.llh.h, ns.velNed);
    Vec3 gn = {0, 0, normalGravity(ns.llh.lat, ns.llh.h)};

    ns.velNed = ns.velNed + (fn - cross(wie_n * 2.0 + wen_n, ns.velNed) + gn) * dt;

    double s2 = std::sin(ns.llh.lat); s2 *= s2;
    double Rn = kA / std::sqrt(1 - kE2 * s2);
    double Rm = Rn * (1 - kE2) / (1 - kE2 * s2);
    ns.llh.lat += ns.velNed.x / (Rm + ns.llh.h) * dt;
    ns.llh.lon += ns.velNed.y / ((Rn + ns.llh.h) * std::cos(ns.llh.lat)) * dt;
    ns.llh.h  -= ns.velNed.z * dt;
    ns.ecef = llhToEcef(ns.llh);
    ns.time.tow += dt;

    // 误差状态传播矩阵：把“名义状态”上的小误差一起传递下去。
    Mat3 CbnT = transpose(Cbn);
    Vec3 fb = matVec(CbnT, fn - gn);  // specific force in body frame

    Vec3 win_n = wie_n + wen_n;

    // 5. 构造误差状态转移矩阵的主要块。
    //    F_vphi 表示姿态误差如何影响速度；
    //    F_vba/F_phibg 表示加表/陀螺零偏如何进入速度/姿态误差；
    //    F_vsa/F_phisg 表示比例因子误差的影响。
    Mat3 F_vr, F_vv, F_vphi, F_vba, F_vsa;
    Mat3 F_phir, F_phiv, F_phiphi, F_phibg, F_phisg;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            F_vr.m[i][j]    = 0;
            F_vv.m[i][j]    = -2 * skewSymmetric(wie_n).m[i][j] - skewSymmetric(wen_n).m[i][j];
            F_vphi.m[i][j]  = -skewSymmetric(matVec(Cbn, fb)).m[i][j];
            F_vba.m[i][j]   = Cbn.m[i][j];
            F_vsa.m[i][j]   = (i == j ? -matVec(Cbn, dvel / dt)[i] : 0.0);

            F_phir.m[i][j]  = 0;
            F_phiv.m[i][j]  = 0;  // simplified
            F_phiphi.m[i][j] = -skewSymmetric(win_n).m[i][j];
            F_phibg.m[i][j] = -Cbn.m[i][j];
            F_phisg.m[i][j] = (i == j ? -dtheta[i] / dt : 0.0);
        }
    }

    // 6. 离散状态转移矩阵 Phi。
    //    kDim 实际为 21 维，这里填入主要非零块，其余保持单位阵。
    Matrix Phi = Matrix::identity(kDim);
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            Phi(kVel + r, kVel + c) += F_vv.m[r][c] * dt;
            Phi(kVel + r, kAtt + c) += F_vphi.m[r][c] * dt;
            Phi(kVel + r, kAcc + c) += F_vba.m[r][c] * dt;
            Phi(kVel + r, kAccScale + c) += F_vsa.m[r][c] * dt;
            Phi(kAtt + r, kAtt + c) += F_phiphi.m[r][c] * dt;
            Phi(kAtt + r, kGyr + c) += F_phibg.m[r][c] * dt;
            Phi(kAtt + r, kGyrScale + c) += F_phisg.m[r][c] * dt;
        }
        Phi(kPos + r, kVel + r) = Phi(kPos + r, kVel + r) + dt;
    }

    // 过程噪声：描述 IMU 噪声、零偏随机游走和比例因子的漂移。
    const auto &m = opt_.imuModel;
    Matrix Qd(kDim, kDim);
    for (int i = 0; i < 3; ++i) {
        // velocity noise from accel noise
        Qd(kVel + i, kVel + i) += sqr(m.velNoise[i]) * dt;
        for (int j = 0; j < 3; ++j) {
            // accel bias random walk maps to velocity
            Qd(kVel + i, kVel + j) += Cbn.m[i][j] * sqr(m.accBiasNoise[j]) * dt * sqr(Cbn.m[i][j]) * dt;
            Qd(kAcc + i, kAcc + i) += sqr(m.accBiasNoise[i]) * dt;
            Qd(kGyr + i, kGyr + i) += sqr(m.gyroBiasNoise[i]) * dt;
            Qd(kAccScale + i, kAccScale + i) += sqr(m.accScaleNoise[i]) * dt;
            Qd(kGyrScale + i, kGyrScale + i) += sqr(m.gyroScaleNoise[i]) * dt;
            // gyro bias random walk maps to attitude
            Qd(kAtt + i, kAtt + j) += Cbn.m[i][j] * sqr(m.gyroBiasNoise[j]) * dt * sqr(Cbn.m[i][j]) * dt;
        }
        Qd(kAtt + i, kAtt + i) += sqr(m.attNoise[i]) * dt;
    }

    // 协方差预测：这是卡尔曼滤波里最核心的“越推越不确定”的部分。
    Matrix PhiT = Phi.transpose();
    st_.P = Phi * st_.P * PhiT + Qd;

    // symmetrize
    for (int i = 0; i < kDim; ++i)
        for (int j = i + 1; j < kDim; ++j)
            st_.P(j, i) = st_.P(i, j) = 0.5 * (st_.P(i, j) + st_.P(j, i));
}

// ESKF GNSS 更新步：用 RTK 位置/速度观测修正惯导累计误差。
void EsKf::updatePosVel(const GnssSolution &gnss) {
    // GNSS 量测更新：
    // RTK 解算输出的是 GNSS 天线位置/速度，而滤波名义状态是 IMU 中心，
    // 因此先用杆臂把 IMU 中心外推到天线位置，再与 GNSS 量测构造残差。
    // q=1 固定解和 q=2 浮点解采用不同观测噪声，质量差的历元会被 usableGnssUpdate 剔除。
    //
    // 更新流程：
    //   1. 质量检查：固定解、浮点解分别检查卫星数、STD 和速度范围；
    //   2. 杆臂改正：IMU 中心 + Cbn * leverArm -> GNSS 天线预测位置；
    //   3. 构造残差：GNSS 天线位置/速度 - INS 预测位置/速度；
    //   4. 构造 H/R：位置观测对应位置误差，速度观测对应速度误差；
    //   5. Kalman 更新：dx = K z；
    //   6. 闭环修正：把误差状态 dx 加回名义状态，然后更新协方差。
    if (!usableGnssUpdate(gnss)) return;

    NavState &ns = st_.nominal;
    Mat3 Cbn = quatToDcm(ns.qbn);

    // 天线杆臂改正：GNSS 天线位置和 IMU 位置不是同一个点。
    Vec3 l_ned = matVec(Cbn, opt_.leverArm);
    Vec3 antennaEcef = ns.ecef + nedToEcef(l_ned, ns.llh);

    // 量测残差：GNSS 解算值和当前惯导预测值之间的差。
    Vec3 zPos = ecefToNed(gnss.ecef - antennaEcef, ns.llh);

    const bool useVelocity = norm(gnss.velNed) < (gnss.q == 1 ? 3.0 : 5.0);
    const int measDim = useVelocity ? 6 : 3;

    // Loose coupling uses solved GNSS position/velocity as measurements. Raw
    // Loose coupling uses only GNSS position/velocity residuals.
    Matrix H(measDim, kDim);
    for (int i = 0; i < 3; ++i) H(i, kPos + i) = 1.0;
    if (useVelocity) {
        for (int i = 0; i < 3; ++i) H(3 + i, kVel + i) = 1.0;
    }

    // 观测噪声 R：
    // 固定解 q=1 给较小方差，浮点解 q=2 给较大方差。
    // gnssRScale 用于报告中的 R 超参数敏感性实验。
    Matrix R(measDim, measDim);
    double pStd = gnss.q == 1
        ? std::max(0.08, std::min(0.50, static_cast<double>(gnss.stdN)))
        : std::max(2.0, std::min(12.0, static_cast<double>(gnss.stdN) * 1.5));
    pStd *= std::sqrt(std::max(1.0e-6, opt_.gnssRScale));
    for (int i = 0; i < 3; ++i) R(i, i) = sqr(pStd);
    if (useVelocity) {
        const double vStd = (gnss.q == 1 ? 0.15 : 1.5) * std::sqrt(std::max(1.0e-6, opt_.gnssRScale));
        for (int i = 0; i < 3; ++i) R(3 + i, 3 + i) = sqr(vStd);
    }

    // 卡尔曼更新：用 GNSS 观测把 IMU 漂移拉回来。
    Matrix Ht = H.transpose();
    Matrix S = H * st_.P * Ht + R;
    Matrix K;
    try {
        Matrix Si = inverse(S);
        K = st_.P * Ht * Si;
    } catch (...) { return; }

    std::vector<double> z(measDim, 0.0);
    z[0] = zPos.x; z[1] = zPos.y; z[2] = zPos.z;
    if (useVelocity) {
        z[3] = gnss.velNed.x - ns.velNed.x;
        z[4] = gnss.velNed.y - ns.velNed.y;
        z[5] = gnss.velNed.z - ns.velNed.z;
    }
    // dx 是 21 维误差状态估计。本实现为了保持姿态稳定，
    // 对姿态和 IMU 参数项做了保守处理，主要用 GNSS 约束位置/速度漂移。
    std::vector<double> dx = K * z;
    for (int i = kAtt; i < kDim; ++i) dx[i] = 0.0;

    // apply error correction to nominal state
    applyErrorToState(ns, dx);
    // 强闭合策略：
    // 对高质量固定解，直接把位置/速度闭合到 GNSS；
    // 对质量较好的浮点解，只做较保守的水平修正。
    if (gnss.q == 1 || (opt_.strongGnssClosure && gnss.q == 2 && gnss.stdN <= 3.0 && gnss.stdD <= 1.0)) {
        ns.ecef = gnss.ecef - nedToEcef(l_ned, gnss.llh);
        ns.llh = ecefToLlh(ns.ecef);
        if (useVelocity) ns.velNed = gnss.velNed;
    } else if (opt_.strongGnssClosure && gnss.q == 2 && gnss.stdN <= 3.0) {
        Vec3 corr = ecefToNed(gnss.ecef - antennaEcef, ns.llh);
        corr.z = 0.0;
        ns.ecef = ns.ecef + nedToEcef(corr, ns.llh);
        ns.llh = ecefToLlh(ns.ecef);
        if (useVelocity) {
            ns.velNed.x = gnss.velNed.x;
            ns.velNed.y = gnss.velNed.y;
        }
    }

    // Joseph 形式协方差更新：数值上更稳，不容易把协方差算坏。
    Matrix I_KH = Matrix::identity(kDim) - K * H;
    Matrix I_KHt = I_KH.transpose();
    st_.P = I_KH * st_.P * I_KHt + K * R * K.transpose();
    if (useVelocity) {
        const double vVar = sqr(gnss.q == 1 ? 0.15 : 1.5);
        for (int i = 0; i < 3; ++i) st_.P(kVel + i, kVel + i) = std::min(st_.P(kVel + i, kVel + i), vVar);
    }

    for (int i = 0; i < kDim; ++i)
        for (int j = i + 1; j < kDim; ++j)
            st_.P(j, i) = st_.P(i, j) = 0.5 * (st_.P(i, j) + st_.P(j, i));
}

// ESKF NHC 更新步：用车辆非完整性约束抑制侧向和垂向速度误差。
void EsKf::updateVehicleConstraint() {
    // NHC 非完整性约束：
    // 车辆正常行驶时，车体系侧向和垂向速度应接近 0。
    // 该约束在 GNSS 弱观测或低速段能帮助抑制航向和横向速度漂移。
    //
    // 更新思想：
    //   1. 将导航系速度 vn 转到车体系/机体系 vb；
    //   2. 构造虚拟观测 vb_y = 0, vb_z = 0；
    //   3. H 矩阵对速度误差直接求导，对姿态误差用小扰动数值差分；
    //   4. 用 Kalman 更新修正速度和姿态小误差。
    NavState &ns = st_.nominal;
    const double speedH = std::sqrt(sqr(ns.velNed.x) + sqr(ns.velNed.y));
    if (speedH < 0.5 || speedH > 20.0) return;

    Mat3 Cbn = quatToDcm(ns.qbn);
    Mat3 Cnb = transpose(Cbn);
    Vec3 vb = matVec(Cnb, ns.velNed);

    // 非完整性约束：车辆侧向和垂向速度应接近 0。
    Matrix H(2, kDim);
    H(0, kVel + 0) = Cnb.m[1][0];
    H(0, kVel + 1) = Cnb.m[1][1];
    H(0, kVel + 2) = Cnb.m[1][2];
    H(1, kVel + 0) = Cnb.m[2][0];
    H(1, kVel + 1) = Cnb.m[2][1];
    H(1, kVel + 2) = Cnb.m[2][2];

    const double eps = 1.0e-5;
    for (int a = 0; a < 3; ++a) {
        std::vector<double> d(kDim, 0.0);
        d[kAtt + a] = eps;
        NavState sp = ns;
        applyErrorToState(sp, d);
        Vec3 vbp = matVec(transpose(quatToDcm(sp.qbn)), sp.velNed);
        H(0, kAtt + a) = (vbp.y - vb.y) / eps;
        H(1, kAtt + a) = (vbp.z - vb.z) / eps;
    }

    Matrix R(2, 2);
    const double nhcScale = std::sqrt(std::max(1.0e-6, opt_.nhcRScale));
    R(0, 0) = sqr(0.25 * nhcScale);
    R(1, 1) = sqr(0.35 * nhcScale);

    Matrix Ht = H.transpose();
    Matrix S = H * st_.P * Ht + R;
    Matrix K;
    try {
        K = st_.P * Ht * inverse(S);
    } catch (...) { return; }

    std::vector<double> z(2, 0.0);
    z[0] = -vb.y;
    z[1] = -vb.z;
    std::vector<double> dx = K * z;
    for (int i = kPos; i < kVel; ++i) dx[i] = 0.0;
    for (int i = kAcc; i < kDim; ++i) dx[i] = 0.0;
    applyErrorToState(ns, dx);

    Matrix I_KH = Matrix::identity(kDim) - K * H;
    st_.P = I_KH * st_.P * I_KH.transpose() + K * R * K.transpose();
    for (int i = 0; i < kDim; ++i)
        for (int j = i + 1; j < kDim; ++j)
            st_.P(j, i) = st_.P(i, j) = 0.5 * (st_.P(i, j) + st_.P(j, i));
}

// 当前实现采用闭环 ESKF，误差已在更新函数中写回名义状态，因此这里保留为空接口。
void EsKf::correctNominal() {
    // biases are already applied; this is a no-op in the current closed-loop design
}

// RTS 后向平滑：离线从终点向起点传播平滑改正，不能用于在线伪实时。
void rtsSmooth(std::vector<NavState> &states, std::vector<Matrix> &covs,
               const std::vector<ImuSample> &imu) {
    int N = static_cast<int>(states.size());
    if (N < 2) return;
    const std::vector<NavState> fwdStates = states;
    const std::vector<Matrix> fwdCovs = covs;

    for (int k = N - 2; k >= 0; --k) {
        // Compute Phi from imu[k+1]
        double dt = imu[k + 1].dt > 0 ? imu[k + 1].dt : 0.005;
        NavState &ns = states[k];
        Vec3 dtheta = imu[k + 1].dtheta - ns.gyroBias * dt;
        Vec3 dvel   = imu[k + 1].dvel   - ns.accBias  * dt;
        Mat3 Cbn = quatToDcm(ns.qbn);
        Vec3 fn = matVec(Cbn, dvel) / dt;
        Vec3 wie_n = earthRateNed(ns.llh.lat);
        Vec3 wen_n = transportRateNed(ns.llh.lat, ns.llh.h, ns.velNed);
        Vec3 gn = {0, 0, normalGravity(ns.llh.lat, ns.llh.h)};
        Mat3 CbnT = transpose(Cbn);
        Vec3 fb = matVec(CbnT, fn - gn);
        Vec3 win_n = wie_n + wen_n;

        Matrix Phi = Matrix::identity(kDim);
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                Phi(kVel + r, kVel + c) += (-2 * skewSymmetric(wie_n).m[r][c] - skewSymmetric(wen_n).m[r][c]) * dt;
                Phi(kVel + r, kAtt + c) += -skewSymmetric(matVec(Cbn, fb)).m[r][c] * dt;
                Phi(kVel + r, kAcc + c) += Cbn.m[r][c] * dt;
                Phi(kAtt + r, kAtt + c) += -skewSymmetric(win_n).m[r][c] * dt;
                Phi(kAtt + r, kGyr + c) += -Cbn.m[r][c] * dt;
            }
            Phi(kPos + r, kVel + r) = Phi(kPos + r, kVel + r) + dt;
        }

        // P_k+1^- = Phi * P_k * Phi' + Qd  (already from forward pass approximation)
        Matrix P_fwd = fwdCovs[k + 1];
        Matrix P_pred = Phi * fwdCovs[k] * Phi.transpose();

        // smoother gain: G_k = P_k * Phi' * inv(P_pred)
        Matrix G;
        try {
            G = fwdCovs[k] * Phi.transpose() * inverse(P_pred);
        } catch (...) { continue; }

        // smoothed covariance
        covs[k] = fwdCovs[k] + G * (covs[k + 1] - P_pred) * G.transpose();

        // error state correction
        std::vector<double> nextErr(kDim, 0.0);
        Vec3 dp = ecefToNed(states[k + 1].ecef - fwdStates[k + 1].ecef, fwdStates[k + 1].llh);
        nextErr[kPos + 0] = dp.x;
        nextErr[kPos + 1] = dp.y;
        nextErr[kPos + 2] = dp.z;
        nextErr[kVel + 0] = states[k + 1].velNed.x - fwdStates[k + 1].velNed.x;
        nextErr[kVel + 1] = states[k + 1].velNed.y - fwdStates[k + 1].velNed.y;
        nextErr[kVel + 2] = states[k + 1].velNed.z - fwdStates[k + 1].velNed.z;
        Quat dqNext = quatMul(states[k + 1].qbn, quatConj(fwdStates[k + 1].qbn));
        nextErr[kAtt + 0] = 2.0 * dqNext.x;
        nextErr[kAtt + 1] = 2.0 * dqNext.y;
        nextErr[kAtt + 2] = 2.0 * dqNext.z;
        nextErr[kAcc + 0] = states[k + 1].accBias.x - fwdStates[k + 1].accBias.x;
        nextErr[kAcc + 1] = states[k + 1].accBias.y - fwdStates[k + 1].accBias.y;
        nextErr[kAcc + 2] = states[k + 1].accBias.z - fwdStates[k + 1].accBias.z;
        nextErr[kGyr + 0] = states[k + 1].gyroBias.x - fwdStates[k + 1].gyroBias.x;
        nextErr[kGyr + 1] = states[k + 1].gyroBias.y - fwdStates[k + 1].gyroBias.y;
        nextErr[kGyr + 2] = states[k + 1].gyroBias.z - fwdStates[k + 1].gyroBias.z;
        nextErr[kAccScale + 0] = states[k + 1].accScale.x - fwdStates[k + 1].accScale.x;
        nextErr[kAccScale + 1] = states[k + 1].accScale.y - fwdStates[k + 1].accScale.y;
        nextErr[kAccScale + 2] = states[k + 1].accScale.z - fwdStates[k + 1].accScale.z;
        nextErr[kGyrScale + 0] = states[k + 1].gyroScale.x - fwdStates[k + 1].gyroScale.x;
        nextErr[kGyrScale + 1] = states[k + 1].gyroScale.y - fwdStates[k + 1].gyroScale.y;
        nextErr[kGyrScale + 2] = states[k + 1].gyroScale.z - fwdStates[k + 1].gyroScale.z;
        std::vector<double> dxErr = G * nextErr;
        states[k] = fwdStates[k];
        applyErrorToState(states[k], dxErr);

        // symmetrize
        for (int i = 0; i < kDim; ++i)
            for (int j = i + 1; j < kDim; ++j)
                covs[k](j, i) = covs[k](i, j) = 0.5 * (covs[k](i, j) + covs[k](j, i));
    }
}

// 简化鲁棒后向平滑，只对相邻状态做小权重平滑，避免异常跳点扩散。
static void robustBackwardSmooth(std::vector<NavState> &states) {
    if (states.size() < 2) return;
    for (int i = static_cast<int>(states.size()) - 2; i >= 0; --i) {
        const double dt = states[i + 1].time.tow - states[i].time.tow;
        if (dt <= 0.0 || dt > 0.02) continue;
        if (norm(states[i + 1].ecef - states[i].ecef) > 20.0) continue;
        const double a = 0.08;
        states[i].ecef = states[i].ecef * (1.0 - a) + states[i + 1].ecef * a;
        states[i].llh = ecefToLlh(states[i].ecef);
        states[i].velNed = states[i].velNed * (1.0 - a) + states[i + 1].velNed * a;
    }
}

// ---- pipeline ----
// 松组合前向滤波主流程：按 IMU 时间顺序 predict，GNSS 到达时 update。
std::vector<FusionOutput> runLooseFusion(const std::vector<ImuSample> &imu,
                                         const std::vector<GnssSolution> &gnss,
                                         const FusionOptions &opt, const NavState &init) {
    // forward pass
    std::vector<FusionOutput> fwd; fwd.reserve(imu.size());
    std::vector<NavState> fwdStates; fwdStates.reserve(imu.size());
    std::vector<Matrix> fwdCovs; fwdCovs.reserve(imu.size());

    EsKf kf(opt, init);
    size_t ig = 0;
    double nextConstraintTow = init.time.tow;
    for (const auto &s : imu) {
        if (s.time < init.time.tow) continue;
        kf.predict(s);
        if (opt.useVehicleConstraint && s.time + 1.0e-9 >= nextConstraintTow) {
            kf.updateVehicleConstraint();
            nextConstraintTow = s.time + std::max(0.005, opt.vehicleConstraintInterval);
        }
        while (ig < gnss.size() && gnss[ig].time.tow <= s.time + 0.002) {
            kf.updatePosVel(gnss[ig]);
            ++ig;
        }
        FusionOutput out;
        out.state = kf.state();
        out.P = kf.P();
        out.gnssQuality = ig > 0 ? gnss[ig - 1].q : 0;
        out.mode = ig > 0 ? gnss[ig - 1].mode : "INS";
        fwd.push_back(out);
        fwdStates.push_back(out.state);
        fwdCovs.push_back(kf.P());
    }

    return fwd;
}

// 基于前向结果执行 RTS 离线平滑，并输出带 _RTS 标记的结果。
std::vector<FusionOutput> rtsLooseFusionPost(const std::vector<FusionOutput> &loose,
                                             const std::vector<ImuSample> &imu) {
    if (loose.size() < 2 || imu.size() < 2) return loose;

    std::vector<NavState> states;
    std::vector<Matrix> covs;
    states.reserve(loose.size());
    covs.reserve(loose.size());
    for (const auto &o : loose) {
        states.push_back(o.state);
        covs.push_back(o.P);
    }

    std::vector<ImuSample> alignedImu;
    alignedImu.reserve(loose.size());
    size_t ii = 0;
    for (const auto &o : loose) {
        while (ii + 1 < imu.size() &&
               std::abs(imu[ii + 1].time - o.state.time.tow) <
                   std::abs(imu[ii].time - o.state.time.tow)) {
            ++ii;
        }
        alignedImu.push_back(imu[ii]);
    }

    rtsSmooth(states, covs, alignedImu);

    std::vector<FusionOutput> out = loose;
    for (size_t i = 0; i < out.size(); ++i) {
        out[i].state = states[i];
        out[i].P = covs[i];
        scaleCovariance(out[i].P, 0.80);
        out[i].mode = out[i].mode + "_RTS";
    }
    applyCourseAttitudeAid(out, 0.010, 0.015);
    return out;
}

// 基于 GNSS 锚点的离线平滑，用插值改正量抑制局部噪声。
std::vector<FusionOutput> smoothLooseFusionPost(const std::vector<FusionOutput> &loose,
                                                const std::vector<GnssSolution> &gnss) {
    if (loose.empty() || gnss.empty()) return loose;

    std::vector<FusionAnchor> anchors = collectFusionAnchors(loose, gnss);
    if (anchors.size() < 2) return loose;

    std::vector<FusionOutput> out = loose;
    size_t ia = 0;
    for (auto &o : out) {
        const double tow = o.state.time.tow;
        while (ia + 1 < anchors.size() && anchors[ia + 1].tow < tow) ++ia;
        FusionAnchor corr = interpolateAnchor(anchors, ia, tow);

        o.state.ecef = o.state.ecef + nedToEcef(corr.posCorrNed, o.state.llh);
        o.state.llh = ecefToLlh(o.state.ecef);
        if (corr.hasVelocity) o.state.velNed = o.state.velNed + corr.velCorrNed;
        scaleCovariance(o.P, 0.85);
        o.mode = o.mode + "_POST_SMOOTH";
    }
    applyCourseAttitudeAid(out, 0.015, 0.020);
    return out;
}

// 后向离线处理：从后往前利用未来锚点修正前面历元，只用于成果对比。
std::vector<FusionOutput> backwardLooseFusionPost(const std::vector<FusionOutput> &loose,
                                                  const std::vector<GnssSolution> &gnss) {
    if (loose.empty() || gnss.empty()) return loose;
    std::vector<FusionAnchor> anchors = collectFusionAnchors(loose, gnss);
    if (anchors.size() < 2) return loose;

    std::vector<FusionOutput> out = loose;
    size_t ia = anchors.size() - 1;
    for (int i = static_cast<int>(out.size()) - 1; i >= 0; --i) {
        const double tow = out[i].state.time.tow;
        while (ia > 0 && anchors[ia].tow > tow) --ia;
        FusionAnchor corr = interpolateAnchor(anchors, ia, tow, 0.65);

        const double ahead = ia + 1 < anchors.size() ? anchors[ia + 1].tow - tow : 0.0;
        const double weight = std::max(0.55, std::min(1.0, 1.0 - ahead / 12.0));
        corr.posCorrNed = corr.posCorrNed * weight;
        corr.velCorrNed = corr.velCorrNed * weight;

        out[i].state.ecef = out[i].state.ecef + nedToEcef(corr.posCorrNed, out[i].state.llh);
        out[i].state.llh = ecefToLlh(out[i].state.ecef);
        if (corr.hasVelocity) out[i].state.velNed = out[i].state.velNed + corr.velCorrNed;
        scaleCovariance(out[i].P, 0.78);
        out[i].mode = out[i].mode + "_BACKWARD";
    }

    for (int i = static_cast<int>(out.size()) - 2; i >= 0; --i) {
        const double dt = out[i + 1].state.time.tow - out[i].state.time.tow;
        if (dt <= 0.0 || dt > 0.02) continue;
        if (norm(out[i + 1].state.ecef - out[i].state.ecef) > 20.0) continue;
        const double a = 0.10;
        out[i].state.ecef = out[i].state.ecef * (1.0 - a) + out[i + 1].state.ecef * a;
        out[i].state.llh = ecefToLlh(out[i].state.ecef);
        out[i].state.velNed = out[i].state.velNed * (1.0 - a) + out[i + 1].state.velNed * a;
    }
    applyCourseAttitudeAid(out, 0.025, 0.035);
    return out;
}

// 将组合导航结果写成 CSV，包括位置、速度、姿态、质量标志和 STD。
void writeFusionOutputs(const std::string &path, const std::vector<FusionOutput> &outputs) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write fusion output: " + path);
    out << "# GPSTime ECEF-X ECEF-Y ECEF-Z Lat-deg Lon-deg H vn ve vd roll pitch yaw q mode "
           "accBx accBy accBz gyrBx gyrBy gyrBz accSx accSy accSz gyrSx gyrSy gyrSz "
           "stdPN stdPE stdPD stdVN stdVE stdVD stdRoll stdPitch stdYaw "
           "stdAccSx stdAccSy stdAccSz stdGyrSx stdGyrSy stdGyrSz\n";
    out << std::fixed << std::setprecision(6);
    for (const auto &o : outputs) {
        const NavState &s = o.state;
        Vec3 euler = quatToEuler(s.qbn);
        out << std::setprecision(3) << s.time.tow << ' '
            << std::setprecision(4) << s.ecef.x << ' ' << s.ecef.y << ' ' << s.ecef.z << ' '
            << std::setprecision(9) << s.llh.lat * kR2D << ' ' << s.llh.lon * kR2D << ' '
            << std::setprecision(4) << s.llh.h << ' '
            << s.velNed.x << ' ' << s.velNed.y << ' ' << s.velNed.z << ' '
            << euler.x * kR2D << ' ' << euler.y * kR2D << ' ' << euler.z * kR2D << ' '
            << o.gnssQuality << ' ' << o.mode << ' '
            << s.accBias.x << ' ' << s.accBias.y << ' ' << s.accBias.z << ' '
            << s.gyroBias.x << ' ' << s.gyroBias.y << ' ' << s.gyroBias.z << ' '
            << s.accScale.x << ' ' << s.accScale.y << ' ' << s.accScale.z << ' '
            << s.gyroScale.x << ' ' << s.gyroScale.y << ' ' << s.gyroScale.z;
        if (o.P.rows() >= kDim && o.P.cols() >= kDim) {
            for (int i = kPos; i < kPos + 3; ++i) out << ' ' << std::sqrt(std::max(0.0, o.P(i, i)));
            for (int i = kVel; i < kVel + 3; ++i) out << ' ' << std::sqrt(std::max(0.0, o.P(i, i)));
            for (int i = kAtt; i < kAtt + 3; ++i) out << ' ' << std::sqrt(std::max(0.0, o.P(i, i))) * kR2D;
            for (int i = kAccScale; i < kAccScale + 3; ++i) out << ' ' << std::sqrt(std::max(0.0, o.P(i, i)));
            for (int i = kGyrScale; i < kGyrScale + 3; ++i) out << ' ' << std::sqrt(std::max(0.0, o.P(i, i)));
        }
        out << '\n';
    }
}

} // namespace gnssins








