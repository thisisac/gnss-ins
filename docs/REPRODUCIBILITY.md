# Reproducibility Notes

## 中文

仓库不包含原始 GNSS/IMU 数据、商业软件参考轨迹和个人报告。若要复现实验，需要自行准备与程序格式兼容的 RINEX 观测、广播星历、Base/Rover 数据和 IMR 文件，并在本地建立：

```text
data/
  rover.26o
  base.rnx
  navigation.26p
  imu.dat.imr
reference/
  GNSSINS/reference.pos
  rtk/reference.pos
output/
```

由于采样率、时间同步、基站坐标、初始航向、杆臂和滤波器参数都会影响结果，仓库中的示例图不构成逐字节数值基准。

## English

Raw GNSS/IMU data, commercial reference trajectories, and personal reports are not included. To reproduce a run, provide compatible RINEX observations, broadcast ephemerides, Base/Rover data, and IMR input under the documented `data/` and `reference/` layout.

Exact values are not guaranteed because sampling, time synchronization, base coordinates, initial heading, lever arm, and filter parameters all affect the solution. The archived figures are illustrative outputs rather than byte-for-byte regression fixtures.
