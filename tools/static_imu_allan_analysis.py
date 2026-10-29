from pathlib import Path
import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
STATIC_IMU = ROOT / "data" / "IMUjingtai" / "imu_20260603_100606_559947.dat"
OUT_DIR = ROOT / "output" / "成果_MATLAB分析" / "静态IMU_Allan分析"

FRAME_OFFSET = 28
FRAME_LEN = 50
FS = 200.0
DT = 1.0 / FS
MAX_SECONDS = 7200.0

GYRO_POS = (20, 24, 28)
ACC_POS = (32, 36, 40)
ACC_SCALE = 1.0 / (2.0 ** 24)


def int32_be(frames, pos):
    b = frames[:, pos:pos + 4].astype(np.int64)
    u = (b[:, 0] << 24) | (b[:, 1] << 16) | (b[:, 2] << 8) | b[:, 3]
    u = np.where(u >= 2 ** 31, u - 2 ** 32, u)
    return u.astype(np.float64)


def read_static_aa77(path):
    data = path.read_bytes()
    n = (len(data) - FRAME_OFFSET) // FRAME_LEN
    frames = np.frombuffer(data, dtype=np.uint8, count=n * FRAME_LEN, offset=FRAME_OFFSET)
    frames = frames.reshape(n, FRAME_LEN)
    valid = (frames[:, 0] == 0xAA) & (frames[:, 1] == 0x77)
    frames = frames[valid]
    max_n = int(MAX_SECONDS * FS)
    if len(frames) > max_n:
        frames = frames[:max_n]

    gyro_counts = np.column_stack([int32_be(frames, p) for p in GYRO_POS])
    acc_counts = np.column_stack([int32_be(frames, p) for p in ACC_POS])
    acc_ms2 = acc_counts * ACC_SCALE

    return frames, gyro_counts, acc_counts, acc_ms2, valid.sum()


def allan_deviation(x, dt):
    n = x.shape[0]
    max_m = max(1, n // 10)
    m_list = np.unique(np.floor(np.logspace(0, math.log10(max_m), 70)).astype(int))
    tau = m_list.astype(np.float64) * dt
    adev = np.full((len(m_list), x.shape[1]), np.nan)
    for i, m in enumerate(m_list):
        nb = n // m
        if nb < 3:
            continue
        y = x[:nb * m].reshape(nb, m, x.shape[1]).mean(axis=1)
        dy = np.diff(y, axis=0)
        adev[i] = np.sqrt(0.5 * np.mean(dy * dy, axis=0))
    return tau, adev


def safe_font(size):
    candidates = [
        r"C:\Windows\Fonts\msyh.ttc",
        r"C:\Windows\Fonts\simsun.ttc",
        r"C:\Windows\Fonts\arial.ttf",
    ]
    for item in candidates:
        p = Path(item)
        if p.exists():
            try:
                return ImageFont.truetype(str(p), size)
            except OSError:
                pass
    return ImageFont.load_default()


def draw_log_plot(path, tau, adev, title, ylabel):
    w, h = 1200, 820
    margin_l, margin_r, margin_t, margin_b = 110, 40, 70, 105
    img = Image.new("RGB", (w, h), "white")
    d = ImageDraw.Draw(img)
    font = safe_font(26)
    small = safe_font(20)
    tiny = safe_font(16)

    valid = np.isfinite(adev).all(axis=1) & (tau > 0) & np.all(adev > 0, axis=1)
    x = np.log10(tau[valid])
    y = np.log10(adev[valid])
    xmin, xmax = float(np.floor(x.min())), float(np.ceil(x.max()))
    ymin, ymax = float(np.floor(y.min())), float(np.ceil(y.max()))
    if ymin == ymax:
        ymin -= 1
        ymax += 1

    def px(xv):
        return margin_l + (xv - xmin) / (xmax - xmin) * (w - margin_l - margin_r)

    def py(yv):
        return h - margin_b - (yv - ymin) / (ymax - ymin) * (h - margin_t - margin_b)

    d.rectangle([margin_l, margin_t, w - margin_r, h - margin_b], outline=(30, 30, 30), width=2)

    for gx in range(int(xmin), int(xmax) + 1):
        xx = px(gx)
        d.line([xx, margin_t, xx, h - margin_b], fill=(225, 225, 225))
        d.text((xx - 20, h - margin_b + 10), f"1e{gx}", fill=(0, 0, 0), font=tiny)
    for gy in range(int(ymin), int(ymax) + 1):
        yy = py(gy)
        d.line([margin_l, yy, w - margin_r, yy], fill=(225, 225, 225))
        d.text((20, yy - 10), f"1e{gy}", fill=(0, 0, 0), font=tiny)

    colors = [(220, 40, 40), (30, 150, 60), (40, 80, 220)]
    labels = ["X", "Y", "Z"]
    for j in range(3):
        pts = [(px(xv), py(yv)) for xv, yv in zip(x, y[:, j])]
        if len(pts) > 1:
            d.line(pts, fill=colors[j], width=3)
        d.line([w - 220, margin_t + 35 + 28 * j, w - 170, margin_t + 35 + 28 * j], fill=colors[j], width=4)
        d.text((w - 160, margin_t + 22 + 28 * j), labels[j], fill=(0, 0, 0), font=small)

    d.text((w // 2 - 260, 25), title, fill=(0, 0, 0), font=font)
    d.text((w // 2 - 90, h - 55), "平均时间 tau (s)", fill=(0, 0, 0), font=small)
    d.text((15, 25), ylabel, fill=(0, 0, 0), font=small)
    img.save(path)


def write_summary(path, frames, gyro_counts, acc_counts, acc_ms2, tau_g, adev_g, tau_a, adev_a, valid_total):
    duration = len(frames) / FS
    acc_norm = np.linalg.norm(acc_ms2, axis=1)
    gyro_mean = gyro_counts.mean(axis=0)
    gyro_std = gyro_counts.std(axis=0)
    acc_mean = acc_ms2.mean(axis=0)
    acc_std = acc_ms2.std(axis=0)

    def fmt3(v):
        return f"{v[0]:.8e}, {v[1]:.8e}, {v[2]:.8e}"

    lines = []
    lines.append("IMUjingtai 静态 IMU Allan 分析结果")
    lines.append("=" * 36)
    lines.append(f"输入文件: {STATIC_IMU}")
    lines.append("解析格式: AA77 帧头，50 字节定长帧，6 个大端 int32 原始量")
    lines.append(f"文件有效 AA77 帧数: {valid_total}")
    lines.append(f"本次分析使用帧数: {len(frames)}")
    lines.append(f"采样频率: {FS:.3f} Hz")
    lines.append(f"分析时长: {duration:.3f} s")
    lines.append("字段解释: 前 3 个 int32 作为陀螺原始计数；后 3 个 int32 作为加速度计原始计数")
    lines.append("比例说明: 加速度计按 1/2^24 转为 m/s^2，静态范数接近重力；陀螺缺少厂家比例因子，保留原始计数单位。")
    lines.append("")
    lines.append(f"陀螺原始计数均值: {fmt3(gyro_mean)}")
    lines.append(f"陀螺原始计数标准差: {fmt3(gyro_std)}")
    lines.append(f"加速度计均值 (m/s^2): {fmt3(acc_mean)}")
    lines.append(f"加速度计标准差 (m/s^2): {fmt3(acc_std)}")
    lines.append(f"加速度计范数均值 (m/s^2): {acc_norm.mean():.8e}")
    lines.append(f"加速度计范数标准差 (m/s^2): {acc_norm.std():.8e}")
    lines.append("")
    lines.append("陀螺 Allan deviation 最小值，单位 raw count/s:")
    lines.append(f"  X {np.nanmin(adev_g[:, 0]):.8e}")
    lines.append(f"  Y {np.nanmin(adev_g[:, 1]):.8e}")
    lines.append(f"  Z {np.nanmin(adev_g[:, 2]):.8e}")
    lines.append("")
    lines.append("加速度计 Allan deviation 最小值，单位 m/s^2:")
    lines.append(f"  X {np.nanmin(adev_a[:, 0]):.8e}")
    lines.append(f"  Y {np.nanmin(adev_a[:, 1]):.8e}")
    lines.append(f"  Z {np.nanmin(adev_a[:, 2]):.8e}")
    lines.append("")
    lines.append("说明: 该结果真正使用 data/IMUjingtai 中的静态 IMU 数据。由于缺少陀螺厂家比例因子，陀螺 Allan 以原始计数表示，不能直接等同于 rad/s 标定值。")
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frames, gyro_counts, acc_counts, acc_ms2, valid_total = read_static_aa77(STATIC_IMU)
    tau_g, adev_g = allan_deviation(gyro_counts, DT)
    tau_a, adev_a = allan_deviation(acc_ms2, DT)

    np.savetxt(OUT_DIR / "static_gyro_allan_raw_counts.csv",
               np.column_stack([tau_g, adev_g]),
               delimiter=",", header="tau_s,gyro_x_raw,gyro_y_raw,gyro_z_raw", comments="")
    np.savetxt(OUT_DIR / "static_acc_allan_ms2.csv",
               np.column_stack([tau_a, adev_a]),
               delimiter=",", header="tau_s,acc_x_ms2,acc_y_ms2,acc_z_ms2", comments="")

    draw_log_plot(OUT_DIR / "static_allan_gyro_raw.png", tau_g, adev_g,
                  "IMUjingtai 静态陀螺 Allan deviation", "raw count/s")
    draw_log_plot(OUT_DIR / "static_allan_acc.png", tau_a, adev_a,
                  "IMUjingtai 静态加速度计 Allan deviation", "m/s^2")
    write_summary(OUT_DIR / "static_allan_summary.txt",
                  frames, gyro_counts, acc_counts, acc_ms2, tau_g, adev_g, tau_a, adev_a, valid_total)
    print(f"static allan done: {OUT_DIR}")


if __name__ == "__main__":
    main()
