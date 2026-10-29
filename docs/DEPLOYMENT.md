# Deployment and Data Policy

## 中文

本项目是研究与课程型 Windows 程序，不提供服务器部署或预编译发布包。公开发布时只应提交源码、脚本、公开许可的示例图和文档。

运行部署步骤：

1. 从 Visual Studio 构建 `src/GNSSINS.sln`；
2. 将 `GNSSINS.exe` 与其运行环境放在受控工作目录；
3. 准备有权限使用的 RINEX、广播星历和 IMR 文件；
4. 调整 `main.cpp` 中的相对文件名或通过本地工作副本配置输入；
5. 将输出写入 `output/`，不要把原始数据或含个人信息的结果提交到 Git。

## English

This is a research/course Windows application, not a server service. No prebuilt binary is distributed. Public releases should contain source, scripts, openly licensed example figures, and documentation only.

Build `src/GNSSINS.sln`, run the executable in a controlled working directory, provide authorized RINEX/broadcast/IMR inputs, keep outputs under `output/`, and never commit raw observations or personally identifying results.
