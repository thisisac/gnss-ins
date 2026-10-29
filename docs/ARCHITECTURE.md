# GNSS/INS Architecture

## 中文概览

本项目是一个“GNSS/RTK 量测 + IMU 机械编排 + 误差状态 Kalman 滤波”的松组合导航程序。主程序在 `main.cpp` 中根据 `AppConfig::mode` 调度各阶段，GNSS、INS 和 Fusion 模块通过明确的数据结构连接。

### 1. GNSS/RTK 模块

`gnss.cpp` 负责读取 RINEX 观测和广播星历，完成卫星状态传播、SPP 初值、站间单差/站星双差建模、Kalman RTK、周跳探测、模糊度浮点处理和候选固定。`runRtkFloat()` 分别执行前向和后向 RTK，再按质量加权组合结果，输出约 1 Hz 的 `GnssSolution`。

### 2. INS 模块

`ins.cpp` 解码 IMR 数据，得到 GPS 周内秒、时间间隔、角增量和速度增量。`coarseAlign()` 用起始静态段估计水平姿态，航向由配置的外部初值提供。`mechanize()` 和 ESKF 预测步骤使用四元数传播姿态，在 NED 坐标系中补偿重力、地球自转和运输角速度。

### 3. Fusion 模块

`fusion.cpp` 使用 15 维误差状态：

```text
[位置误差(3), 速度误差(3), 姿态误差(3),
 加速度计零偏(3), 陀螺零偏(3)]
```

工程还保留比例因子状态字段，用于扩展和参数实验。每帧 IMU 调用 `predict()`；GNSS 到达时调用 `updatePosVel()`；车辆约束开启时调用 `updateVehicleConstraint()`。GNSS 更新支持位置/速度量测定权、RTK 质量分类和 IMU-GNSS 杆臂改正。

### 4. 在线与离线路径

- `forward`：严格按 IMU 时间顺序递推，GNSS 到达后更新；
- `pseudo-online`：使用同一递推逻辑，但把文件数据当作传感器流按时间戳送入；
- `pseudo-online-fast`：在线算法相同，只移除墙钟等待；
- `post`/`all`：在前向结果基础上执行 RTS、后向和前后向平滑。

## English Overview

The project is a loosely coupled navigation pipeline that combines GNSS/RTK measurements, inertial mechanization, and an error-state Kalman filter. `main.cpp` dispatches stages according to `AppConfig::mode`, while GNSS, INS, and Fusion communicate through explicit data structures.

### GNSS/RTK

`gnss.cpp` parses RINEX observations and broadcast ephemerides, propagates satellite states, computes SPP initialization, forms single/double differences, runs Kalman RTK, detects slips, handles float ambiguities, and evaluates fixing candidates. `runRtkFloat()` runs forward and backward passes and combines them into approximately 1 Hz `GnssSolution` records.

### INS

`ins.cpp` decodes IMR records into GPS time-of-week, sample interval, delta angle, and delta velocity. `coarseAlign()` estimates horizontal attitude from an initial static segment while heading is supplied as an external initial value. Quaternion propagation and NED mechanization account for gravity, Earth rotation, and transport rate.

### Fusion

`fusion.cpp` uses a 15-state error model containing position, velocity, attitude, accelerometer-bias, and gyro-bias errors. Each IMU record triggers `predict()`. GNSS records trigger `updatePosVel()`, and optional non-holonomic constraints are applied through `updateVehicleConstraint()`. Measurement weighting, RTK quality, and antenna lever-arm corrections are part of the update model.

### File and pseudo-online paths

`forward` is the ordinary forward filter; `pseudo-online` feeds the same filter one sensor timestamp at a time; `pseudo-online-fast` removes wall-clock sleeps; `post` and `all` add offline RTS/backward/combined smoothing.
