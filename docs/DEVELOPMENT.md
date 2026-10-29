# Development Guide

## 中文

### 工具链

- Visual Studio 2022 或更新版本；
- MSVC C++17；
- Python 3.10+、NumPy、Matplotlib；
- MATLAB（运行 `.m` 分析脚本时需要）。

### 修改流程

1. 在 `src/GNSSINS/src/` 修改 C++ 模块；
2. 使用 `src/GNSSINS.sln` 构建；
3. 将本地输入放在 `data/` 和 `reference/`，结果写入 `output/`；
4. 运行 `tools/` 中相应分析脚本；
5. 提交前检查是否包含绝对路径、个人信息、原始观测文件、二进制程序或本地输出。

### 调试重点

- 检查 IMU 与 GNSS GPS 周内秒是否对齐；
- 检查 `GnssSolution::q`、标准差和卫星数；
- 检查初始航向、杆臂和 NED 轴定义；
- 检查 ESKF 协方差是否保持有限和半正定；
- 对 `pseudo-online` 使用小的 `maxEpochs` 先做快速验证。

## English

Use Visual Studio 2022+, MSVC C++17, Python 3.10+ with NumPy/Matplotlib, and MATLAB for `.m` scripts. Modify C++ files under `src/GNSSINS/src/`, build with `src/GNSSINS.sln`, keep local inputs under `data/` and `reference/`, and write generated files under ignored `output/`.

Before committing, inspect changes for absolute paths, personal information, raw observations, binaries, and local outputs. For debugging, verify GPS time alignment, RTK quality flags, initial heading, lever arm, NED conventions, covariance stability, and small `maxEpochs` pseudo-online runs.
