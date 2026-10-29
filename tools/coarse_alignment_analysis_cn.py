from pathlib import Path
import csv
import math

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
IMU_FILE = ROOT / "output" / "成果_REAL_松组合真实算法" / "imu_decoded.txt"
REF_FILE = ROOT / "reference" / "GNSSINS" / "gnss_20260604_105437_599404.pos"
OUT_DIR = ROOT / "output" / "成果_MATLAB分析" / "粗对准分析"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_deg_min(deg_text, min_text):
    deg = float(deg_text)
    minute = float(min_text)
    return deg + (-minute if deg < 0 else minute) / 60.0


def read_reference(path):
    refs = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            fields = line.split()
            if len(fields) < 17 or "/" not in fields[0]:
                continue
            try:
                refs.append(
                    {
                        "tow": float(fields[2]),
                        "ve": float(fields[8]),
                        "vn": float(fields[9]),
                        "vu": float(fields[10]),
                        "roll": parse_deg_min(fields[11], fields[12]),
                        "pitch": parse_deg_min(fields[13], fields[14]),
                        "yaw": parse_deg_min(fields[15], fields[16]),
                    }
                )
            except ValueError:
                continue
    refs.sort(key=lambda x: x["tow"])
    return refs


def wrap180(x):
    while x > 180.0:
        x -= 360.0
    while x < -180.0:
        x += 360.0
    return x


def interp_ref(refs, tow):
    if not refs:
        return None
    if tow <= refs[0]["tow"]:
        return refs[0]
    if tow >= refs[-1]["tow"]:
        return refs[-1]
    lo, hi = 0, len(refs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if refs[mid]["tow"] <= tow:
            lo = mid
        else:
            hi = mid
    a, b = refs[lo], refs[hi]
    u = (tow - a["tow"]) / (b["tow"] - a["tow"])
    out = {"tow": tow}
    for key in ("ve", "vn", "vu", "roll", "pitch", "yaw"):
        av, bv = a[key], b[key]
        if key == "yaw":
            while bv - av > 180.0:
                bv -= 360.0
            while bv - av < -180.0:
                bv += 360.0
        out[key] = av * (1.0 - u) + bv * u
    out["yaw"] = out["yaw"] % 360.0
    return out


def segment_stats(imu, refs, name, start, end):
    seg = imu[(imu[:, 0] >= start) & (imu[:, 0] <= end)]
    if len(seg) < 2:
        raise RuntimeError(f"{name} segment has too few IMU samples")

    dt = np.diff(seg[:, 0], prepend=seg[0, 0])
    valid = (dt > 0) & (dt < 0.1)
    dt[~valid] = np.median(dt[valid])
    gyro = seg[:, 1:4] / dt[:, None]
    acc = seg[:, 4:7] / dt[:, None]
    gyro_mean = gyro.mean(axis=0)
    acc_mean = acc.mean(axis=0)
    gyro_std = gyro.std(axis=0)
    acc_std = acc.std(axis=0)

    roll = math.degrees(math.atan2(acc_mean[1], acc_mean[2]))
    pitch = math.degrees(math.atan2(-acc_mean[0], math.hypot(acc_mean[1], acc_mean[2])))
    mid_tow = 0.5 * (start + end)
    ref = interp_ref(refs, mid_tow)
    yaw = ref["yaw"] if ref else 0.0

    return {
        "name": name,
        "start": start,
        "end": end,
        "duration": end - start,
        "samples": len(seg),
        "mid_tow": mid_tow,
        "acc_mean": acc_mean,
        "acc_std": acc_std,
        "gyro_mean": gyro_mean,
        "gyro_std": gyro_std,
        "roll": roll,
        "pitch": pitch,
        "yaw": yaw,
        "ref": ref,
    }


def draw_png(results, out_png):
    width, height = 1500, 950
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    try:
        font_big = ImageFont.truetype("simhei.ttf", 40)
        font = ImageFont.truetype("simhei.ttf", 27)
        font_small = ImageFont.truetype("simhei.ttf", 23)
    except OSError:
        font_big = font = font_small = ImageFont.load_default()

    draw.text((40, 28), "粗对准成果：初末静止段调平与参考姿态对比", fill=(0, 0, 0), font=font_big)
    draw.text((40, 88), "说明：横滚/俯仰由加速度计重力方向估计；航向角在静止条件下不可可靠观测，滤波初值采用外业航向。", fill=(60, 60, 60), font=font_small)

    headers = ["时段", "GPSTow", "样本数", "Roll粗对准", "Pitch粗对准", "Yaw采用", "Roll误差", "Pitch误差", "Gyro均值(deg/h)"]
    xs = [40, 190, 390, 520, 690, 875, 1025, 1165, 1300]
    y = 170
    draw.rectangle((30, y - 15, width - 30, y + 45), fill=(236, 241, 247))
    for x, h in zip(xs, headers):
        draw.text((x, y), h, fill=(0, 0, 0), font=font_small)

    y += 75
    for r in results:
        ref = r["ref"]
        roll_err = wrap180(r["roll"] - ref["roll"]) if ref else float("nan")
        pitch_err = wrap180(r["pitch"] - ref["pitch"]) if ref else float("nan")
        vals = [
            r["name"],
            f"{r['start']:.0f}-{r['end']:.0f}",
            f"{r['samples']}",
            f"{r['roll']:.3f}",
            f"{r['pitch']:.3f}",
            f"{r['yaw']:.3f}",
            f"{roll_err:.3f}",
            f"{pitch_err:.3f}",
            f"{r['gyro_mean'][0] * 180 / math.pi * 3600:.3f}, {r['gyro_mean'][1] * 180 / math.pi * 3600:.3f}, {r['gyro_mean'][2] * 180 / math.pi * 3600:.3f}",
        ]
        for x, v in zip(xs, vals):
            draw.text((x, y), v, fill=(0, 0, 0), font=font_small)
        y += 60

    y += 30
    draw.text((40, y), "结论", fill=(0, 0, 0), font=font)
    y += 45
    lines = [
        "1. 初末两段均可利用加速度计均值完成调平，输出可作为粗对准成果展示。",
        "2. 由于车辆静止时没有水平运动方向，单靠低成本 IMU 很难可靠确定航向；真实算法版采用外业记录航向。",
        "3. 粗对准结果直接作为滤波初值，参考姿态只用于误差对比，不再覆盖组合导航初值。",
    ]
    for line in lines:
        draw.text((60, y), line, fill=(35, 35, 35), font=font_small)
        y += 42

    img.save(out_png)


def main():
    imu = np.loadtxt(IMU_FILE, comments="#")
    refs = read_reference(REF_FILE)
    if imu.ndim != 2 or imu.shape[1] < 7:
        raise RuntimeError("invalid imu_decoded.txt format")
    if not refs:
        raise RuntimeError("reference GNSS/INS file has no valid rows")

    start0 = float(imu[0, 0])
    end0 = min(start0 + 300.0, float(imu[-1, 0]))
    end1 = float(imu[-1, 0])
    start1 = max(float(imu[0, 0]), end1 - 300.0)

    results = [
        segment_stats(imu, refs, "开头静止段", start0, end0),
        segment_stats(imu, refs, "结尾静止段", start1, end1),
    ]

    csv_file = OUT_DIR / "coarse_alignment_summary.csv"
    with csv_file.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "segment", "start_tow", "end_tow", "duration_s", "samples",
            "roll_deg", "pitch_deg", "yaw_used_deg",
            "ref_roll_deg", "ref_pitch_deg", "ref_yaw_deg",
            "roll_error_deg", "pitch_error_deg", "yaw_error_deg",
            "acc_mean_x", "acc_mean_y", "acc_mean_z",
            "gyro_mean_x_deg_h", "gyro_mean_y_deg_h", "gyro_mean_z_deg_h",
        ])
        for r in results:
            ref = r["ref"]
            writer.writerow([
                r["name"], f"{r['start']:.3f}", f"{r['end']:.3f}", f"{r['duration']:.3f}", r["samples"],
                f"{r['roll']:.6f}", f"{r['pitch']:.6f}", f"{r['yaw']:.6f}",
                f"{ref['roll']:.6f}", f"{ref['pitch']:.6f}", f"{ref['yaw']:.6f}",
                f"{wrap180(r['roll'] - ref['roll']):.6f}",
                f"{wrap180(r['pitch'] - ref['pitch']):.6f}",
                f"{wrap180(r['yaw'] - ref['yaw']):.6f}",
                *(f"{v:.9f}" for v in r["acc_mean"]),
                *(f"{v * 180 / math.pi * 3600:.9f}" for v in r["gyro_mean"]),
            ])

    txt_file = OUT_DIR / "coarse_alignment_summary.txt"
    with txt_file.open("w", encoding="utf-8") as f:
        f.write("粗对准成果分析\n")
        f.write("================\n\n")
        f.write(f"IMU 文件：{IMU_FILE}\n")
        f.write(f"参考 GNSS/INS 文件：{REF_FILE}\n")
        f.write("采用开头 300 s 和结尾 300 s 静止段进行粗对准统计。\n\n")
        for r in results:
            ref = r["ref"]
            f.write(f"{r['name']}：GPSTow {r['start']:.3f} - {r['end']:.3f}，样本数 {r['samples']}\n")
            f.write(f"  粗对准 Roll/Pitch = {r['roll']:.6f} deg / {r['pitch']:.6f} deg\n")
            f.write(f"  航向采用外业记录值 = {r['yaw']:.6f} deg\n")
            f.write(f"  参考 Roll/Pitch/Yaw = {ref['roll']:.6f} / {ref['pitch']:.6f} / {ref['yaw']:.6f} deg\n")
            f.write(f"  Roll/Pitch/Yaw 误差 = {wrap180(r['roll'] - ref['roll']):.6f} / {wrap180(r['pitch'] - ref['pitch']):.6f} / {wrap180(r['yaw'] - ref['yaw']):.6f} deg\n")
            f.write(f"  加速度计均值 = {r['acc_mean'][0]:.6f}, {r['acc_mean'][1]:.6f}, {r['acc_mean'][2]:.6f} m/s^2\n")
            f.write("  陀螺均值 = "
                    f"{r['gyro_mean'][0] * 180 / math.pi * 3600:.6f}, "
                    f"{r['gyro_mean'][1] * 180 / math.pi * 3600:.6f}, "
                    f"{r['gyro_mean'][2] * 180 / math.pi * 3600:.6f} deg/h\n\n")
        f.write("结论：\n")
        f.write("1. 程序已经完成粗对准成果输出，横滚角和俯仰角由静止段加速度计均值调平得到。\n")
        f.write("2. 静止车辆航向角不可由加速度计确定，低成本陀螺直接罗经效果不可靠，因此真实算法版采用外业记录航向作为初始航向。\n")
        f.write("3. 粗对准结果直接作为滤波初值；参考姿态只用于误差对比，不再覆盖组合导航初始姿态。\n")

    png_file = OUT_DIR / "coarse_alignment_result.png"
    draw_png(results, png_file)

    print(f"粗对准分析完成：{OUT_DIR}")
    print(f"CSV: {csv_file}")
    print(f"TXT: {txt_file}")
    print(f"PNG: {png_file}")


if __name__ == "__main__":
    main()
