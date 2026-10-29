# GNSS/INS Integrated Navigation

一个基于 C++ 的 GNSS/INS 松组合导航项目，包含自编 GNSS/RTK 解算、IMU 解码、惯导机械编排、ESKF 前向滤波、后向/RTS 平滑、伪实时处理和结果可视化工具。

This repository contains a self-developed C++ GNSS/INS loosely coupled navigation system with GNSS/RTK processing, IMU decoding, inertial mechanization, error-state Kalman filtering, offline smoothing, pseudo-online processing, and analysis tools.

个人实习报告、答辩演示文稿和原始说明文档未随仓库发布；仓库仅保留可公开的源码、工具和脱敏后的示例结果。

Personal internship reports, presentation slides, and original private guides are intentionally excluded. Only source code, tools, and sanitized example outputs are included.

## 1. 项目能力 / Capabilities

### 中文

- 读取 RINEX 观测文件与广播星历，支持 GPS、BDS，并保留 GLONASS 处理接口；
- 使用广播星历计算卫星状态，完成 SPP 初值和双差 RTK 浮点/固定解；
- 解码 IMR 惯性数据，输出可检查的 IMU 增量序列；
- 使用静态初始对准估计横滚角、俯仰角，并用外部航向初值启动组合导航；
- 使用四元数完成姿态传播，使用 NED 导航坐标系进行速度和位置更新；
- 使用 15 维误差状态 Kalman 滤波器融合 GNSS 位置/速度与 IMU 预测；
- 支持非完整性约束（NHC）、GNSS 天线杆臂改正、GNSS 质量加权和强 GNSS 闭合；
- 支持前向滤波、RTS 平滑、后向组合和前后向平滑；
- 支持 `pseudo-online` 与 `pseudo-online-fast`，按时间戳模拟 IMU/GNSS 传感器到达；
- 使用 Python 和 MATLAB 工具生成 Allan 方差、初始对准、误差统计、轨迹和伪实时播放结果。

### English

- Read RINEX observation files and broadcast ephemerides, with GPS/BDS support and a retained GLONASS interface;
- Propagate satellite states and compute SPP initialization plus double-difference RTK float/fixed solutions;
- Decode IMR inertial data and export inspectable IMU increments;
- Estimate roll and pitch from a static segment and initialize heading from an external field value;
- Propagate attitude with quaternions and update velocity/position in the NED navigation frame;
- Fuse GNSS position/velocity measurements with IMU prediction using a 15-state error-state Kalman filter;
- Support non-holonomic constraints, antenna lever-arm correction, quality-dependent GNSS weighting, and strong GNSS closure;
- Run forward filtering, RTS smoothing, backward processing, and forward/backward smoothing;
- Provide `pseudo-online` and `pseudo-online-fast` modes that replay IMU/GNSS arrivals by timestamp;
- Generate Allan-variance, alignment, error-statistics, trajectory, and pseudo-realtime visualizations with Python and MATLAB.

## 2. 程序架构 / Architecture

```text
RINEX rover observations + base observations + broadcast ephemerides
                                  |
                                  v
                       GNSS / RTK module
       readRinexObs -> readRinexNav -> SPP -> DD RTK filter
                                  |
                                  v
                         GnssSolution stream
                                  |
IMR inertial file --> decodeImr --> ImuSample stream
                                  |
                                  v
                         INS mechanization
                 coarse alignment + quaternion propagation
                                  |
                                  v
                         ESKF loose fusion
         IMU predict + GNSS update + NHC + lever-arm correction
                                  |
              +-------------------+-------------------+
              |                                       |
              v                                       v
        forward / pseudo-online              offline smoothing
              |                         RTS / backward / combined
              +-------------------+-------------------+
                                  |
                                  v
                    .nav / .pos / diagnostic outputs
                                  |
                                  v
                   Python and MATLAB analysis tools
```

The implementation is split into the following components:

| 模块 / Module | 主要职责 / Responsibility |
| --- | --- |
| `src/GNSSINS/src/common.*` | WGS-84 constants, vector/matrix types, GPS time, ECEF/LLH/NED conversions, and numerical helpers. |
| `src/GNSSINS/src/gnss.*` | RINEX parsing, broadcast ephemeris propagation, SPP initialization, double-difference RTK, ambiguity handling, and GNSS output. |
| `src/GNSSINS/src/ins.*` | IMR decoding, static coarse alignment, quaternion attitude propagation, and inertial mechanization. |
| `src/GNSSINS/src/fusion.*` | 15-state ESKF prediction/update, GNSS measurement weighting, NHC, lever-arm correction, and smoothing interfaces. |
| `src/GNSSINS/src/main.cpp` | Application configuration, mode dispatch, input/output orchestration, and pseudo-online sensor scheduling. |
| `src/GNSSINS/src/imu_alignment_analysis.cpp` | Auxiliary IMU alignment analysis executable source; it is not part of the default Visual Studio target. |
| `tools/` | Python/MATLAB post-processing and visualization scripts. |
| `results/figures/` | Sanitized example figures from the supplied package. |

## 3. 数据流与关键状态 / Data Flow and States

### 中文

`RinexObsFile` 保存单站 RINEX 观测、观测类型和历元；`NavData` 保存广播星历；`GnssSolution` 是 GNSS 模块与组合导航模块之间的接口，包含 ECEF/LLH、NED 速度、质量标志、标准差和模式。

`ImuSample` 表示一帧 IMR 数据，包括 GPS 周内秒、时间间隔、角增量和速度增量。`NavState` 保存名义导航状态：位置、速度、四元数姿态、加速度计零偏、陀螺零偏及比例因子误差。`EsKfState` 额外保存 15 维误差协方差。

### English

`RinexObsFile` stores station observations, observation-type maps, and epochs. `NavData` stores broadcast ephemerides. `GnssSolution` is the interface between GNSS processing and fusion, carrying ECEF/LLH position, NED velocity, quality flags, standard deviations, and solution mode.

`ImuSample` contains GPS time-of-week, sample interval, delta angle, and delta velocity. `NavState` stores nominal position, velocity, quaternion attitude, accelerometer bias, gyro bias, and scale-factor errors. `EsKfState` adds the 15-state error covariance.

## 4. 安装 / Installation

本项目面向 Windows/MSVC。建议安装 Visual Studio 2022（或支持 C++17 的较新版本），并安装 Desktop development with C++ 工作负载。Python 分析工具需要 Python 3.10+、NumPy 和 Matplotlib；MATLAB 工具需要 MATLAB R2020b 或更高版本。

The project targets Windows/MSVC. Install Visual Studio 2022 (or a newer C++17-capable version) with the Desktop development with C++ workload. Python analysis scripts require Python 3.10+, NumPy, and Matplotlib. MATLAB scripts require MATLAB R2020b or later.

## 5. 构建 / Build

1. 打开 `src/GNSSINS.sln`。
2. 选择 `Debug | x64` 或 `Release | x64`。
3. 将 Visual Studio 的 Working Directory 设置为仓库根目录，使默认相对路径 `data/`、`reference/` 和 `output/` 生效。
4. 生成 `GNSSINS` 项目。默认目标使用 C++17，不依赖第三方 C++ 库。

1. Open `src/GNSSINS.sln`.
2. Select `Debug | x64` or `Release | x64`.
3. Set the Visual Studio working directory to the repository root so that `data/`, `reference/`, and `output/` resolve correctly.
4. Build the `GNSSINS` project. The default target uses C++17 and has no third-party C++ dependency.

## 6. 运行 / Usage

程序参数为：

```text
GNSSINS.exe <output-directory> [max-epochs] [mode]
```

支持模式：

| 模式 / Mode | 作用 / Purpose |
| --- | --- |
| `rtk` | 仅运行 GNSS/RTK，输出 `gnss_rtk.pos`。 |
| `forward` | RTK + IMU 解码 + ESKF 前向松组合。 |
| `post` | 读取前向结果并执行离线平滑。 |
| `all` | 执行完整后处理流程。 |
| `pseudo-online` | 按时间戳在线推进 IMU predict 与 GNSS update，并按墙钟时间输出。 |
| `pseudo-online-fast` | 与 `pseudo-online` 算法相同，但不等待真实秒数。 |

The command-line form is:

```text
GNSSINS.exe <output-directory> [max-epochs] [mode]
```

The included example script runs a short fast pseudo-online demo:

```bat
scripts\run_pseudo_online_fast_50.bat
```

Before running, place authorized RINEX rover/base files, broadcast ephemerides, and IMR data under `data/`, and optional reference POS files under `reference/`. These input files are not included in the repository.

## 7. 输出文件 / Outputs

Typical output files include:

- `gnss_rtk.pos`: GNSS/RTK epoch solutions;
- `imu_decoded.txt`: decoded IMU increments;
- `coarse_alignment_result.txt`: static alignment diagnostics;
- `gnss_ins_forward.nav`: forward ESKF output;
- `gnss_ins_pseudo_online.nav`: pseudo-online output;
- `gnss_ins_forward_smooth.nav`, `gnss_ins_backward.nav`, and `gnss_ins_smooth.nav`: offline post-processing outputs;
- `rtk_debug.txt` and `rtk_debug_bwd.txt`: RTK diagnostic logs.

运行结果默认写入 `output/`，该目录被 `.gitignore` 忽略，以避免提交私有观测数据和本地结果。

Runtime results are written to `output/`, which is ignored by `.gitignore` to avoid committing private observations and local outputs.

## 8. 工具 / Analysis Tools

`tools/` 中的脚本默认以仓库根目录为基准读取 `data/`、`reference/` 和 `output/`。常用工具包括：

- `static_imu_allan_analysis.py`: 静态 IMU Allan 方差；
- `initial_alignment_three_methods.py`: 三种初始对准方法对比；
- `analyze_forward_backward_error_std.py`: 前向/后向误差和 STD 统计；
- `make_pseudo_realtime_player.py`: 生成伪实时轨迹 CSV、PNG 和 HTML 播放器；
- `visualize_final_results_cn.m`: MATLAB 最终结果可视化；
- `pure_ins_analysis_cn.m`: 纯惯导误差分析。

The scripts intentionally do not contain user-specific absolute paths. Supply your own authorized data and reference files using the documented directory layout.

## 9. 复现边界 / Reproducibility

本仓库不包含原始 GNSS/IMU 数据、商业软件参考轨迹和个人课程报告，因此不能保证复现实验中的每个数值。结果会受到数据采样、基站坐标、初始航向、杆臂、滤波噪声和 RTK 质量阈值影响。

The repository does not include raw GNSS/IMU data, commercial reference trajectories, or personal course reports. Exact numerical reproduction is therefore not guaranteed. Results depend on sampling, base coordinates, initial heading, lever arm, filter noise, and RTK quality thresholds.

更多设计说明见 [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) 和 [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)。

## 10. 许可证 / License

本项目采用 MIT License，详见 [`LICENSE`](LICENSE)。发布前请确认输入数据、结果图和第三方依赖均拥有可公开许可。

This project is released under the MIT License; see [`LICENSE`](LICENSE). Confirm that input data, figures, and third-party dependencies are redistributable before publishing derived materials.
