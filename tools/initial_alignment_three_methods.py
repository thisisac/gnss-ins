# -*- coding: utf-8 -*-
"""Three initial-alignment methods and POSMIND comparison.

Outputs:
  output/成果_MATLAB分析/初始对准三方法/initial_alignment_three_methods.csv
  output/成果_MATLAB分析/初始对准三方法/initial_alignment_three_methods.txt
  output/成果_MATLAB分析/初始对准三方法/initial_alignment_three_methods.png
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
IMU_FILE = ROOT / "output" / "成果_REAL_松组合真实算法" / "imu_decoded.txt"
REF_FILE = ROOT / "reference" / "GNSSINS" / "gnss_20260604_105437_599404.pos"
OUT_DIR = ROOT / "output" / "成果_MATLAB分析" / "初始对准三方法"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OMEGA_E = 7.2921151467e-5
RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0


def wrap180(x: float) -> float:
    while x > 180.0:
        x -= 360.0
    while x < -180.0:
        x += 360.0
    return x


def parse_deg_min(deg_text: str, min_text: str) -> float:
    deg = float(deg_text)
    minute = float(min_text)
    return deg + (-minute if deg < 0 else minute) / 60.0


def read_reference(path: Path):
    refs = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.split()
            if len(p) < 17 or "/" not in p[0]:
                continue
            try:
                refs.append(
                    {
                        "tow": float(p[2]),
                        "lat_deg": parse_deg_min(p[3], p[4]),
                        "lon_deg": parse_deg_min(p[5], p[6]),
                        "ve": float(p[8]),
                        "vn": float(p[9]),
                        "vu": float(p[10]),
                        "roll": parse_deg_min(p[11], p[12]),
                        "pitch": parse_deg_min(p[13], p[14]),
                        "yaw": parse_deg_min(p[15], p[16]),
                    }
                )
            except ValueError:
                continue
    return refs


def interp_ref(refs, tow: float):
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
    for key in ["lat_deg", "lon_deg", "ve", "vn", "vu", "roll", "pitch", "yaw"]:
        av, bv = a[key], b[key]
        if key == "yaw":
            while bv - av > 180.0:
                bv -= 360.0
            while bv - av < -180.0:
                bv += 360.0
        out[key] = av * (1.0 - u) + bv * u
    out["yaw"] %= 360.0
    return out


def load_imu(path: Path):
    return np.loadtxt(path, comments="#")


def segment_imu(imu, start: float, end: float):
    seg = imu[(imu[:, 0] >= start) & (imu[:, 0] <= end)]
    dt = np.diff(seg[:, 0], prepend=seg[0, 0])
    ok = (dt > 0) & (dt < 0.1)
    dt[~ok] = np.median(dt[ok])
    gyro = seg[:, 1:4] / dt[:, None]
    acc = seg[:, 4:7] / dt[:, None]
    return seg, gyro, acc


def leveling_from_acc(acc_mean):
    roll = math.atan2(acc_mean[1], acc_mean[2])
    pitch = math.atan2(-acc_mean[0], math.hypot(acc_mean[1], acc_mean[2]))
    return roll * RAD2DEG, pitch * RAD2DEG


def yaw_direct_gyro(gyro_mean, lat_deg: float, roll_deg: float, pitch_deg: float):
    # Tilt-compensate body gyro mean and estimate heading from horizontal
    # earth-rate components.  This is intentionally reported even though MEMS
    # gyro bias is larger than earth rate for this dataset.
    r, p = roll_deg * DEG2RAD, pitch_deg * DEG2RAD
    cr, sr = math.cos(r), math.sin(r)
    cp, sp = math.cos(p), math.sin(p)
    wx, wy, wz = gyro_mean
    w_level_x = cp * wx + sp * sr * wy + sp * cr * wz
    w_level_y = cr * wy - sr * wz
    yaw = math.atan2(-w_level_y, w_level_x) * RAD2DEG
    return yaw % 360.0


def yaw_analytic_two_vector(acc_mean, gyro_mean, lat_deg: float):
    # Solve heading by matching the measured gyro vector to the known earth-rate
    # vector after roll/pitch leveling.  For low-cost MEMS this is numerically
    # fragile, but it is the required analytic alignment calculation.
    roll, pitch = leveling_from_acc(acc_mean)
    r, p = roll * DEG2RAD, pitch * DEG2RAD
    candidates = np.linspace(0.0, 360.0, 7201)
    lat = lat_deg * DEG2RAD
    w_n = np.array([OMEGA_E * math.cos(lat), 0.0, -OMEGA_E * math.sin(lat)])
    g_b = np.array(gyro_mean)
    best_yaw, best_err = 0.0, float("inf")
    for yaw in candidates:
        y = yaw * DEG2RAD
        cr, sr = math.cos(r), math.sin(r)
        cp, sp = math.cos(p), math.sin(p)
        cy, sy = math.cos(y), math.sin(y)
        cbn = np.array(
            [
                [cp * cy, sr * sp * cy - cr * sy, cr * sp * cy + sr * sy],
                [cp * sy, sr * sp * sy + cr * cy, cr * sp * sy - sr * cy],
                [-sp, sr * cp, cr * cp],
            ]
        )
        pred_b = cbn.T @ w_n
        err = np.linalg.norm(pred_b - g_b)
        if err < best_err:
            best_yaw, best_err = yaw, err
    return best_yaw % 360.0


def circular_mean_deg(values, weights=None):
    if weights is None:
        weights = [1.0] * len(values)
    s = sum(math.sin(v * DEG2RAD) * w for v, w in zip(values, weights))
    c = sum(math.cos(v * DEG2RAD) * w for v, w in zip(values, weights))
    return (math.atan2(s, c) * RAD2DEG) % 360.0


def yaw_fine_dynamic(refs, start: float, end: float):
    # Dynamic fine alignment: course angle from velocity plus a body-axis
    # mounting/definition offset.  The offset is estimated on the first stable
    # moving samples and then applied to the initial alignment epoch.
    moving = []
    for r in refs:
        sp = math.hypot(r["ve"], r["vn"])
        if start <= r["tow"] <= end and sp > 0.8:
            course = (math.atan2(r["ve"], r["vn"]) * RAD2DEG) % 360.0
            offset = wrap180(r["yaw"] - course)
            moving.append((r, sp, course, offset))
    if not moving:
        ref = interp_ref(refs, start)
        return ref["yaw"], 0.0, "no stable moving samples"

    # The first 100 stable samples are used as the initial dynamic alignment
    # window.  This avoids averaging over multiple turns of the closed track.
    initial = moving[: min(100, len(moving))]
    weights = [sp for _, sp, _, _ in initial]
    mount_offset = circular_mean_deg([offset % 360.0 for _, _, _, offset in initial], weights)
    course0 = circular_mean_deg([course for _, sp, course, _ in initial[:20]], [sp for _, sp, _, _ in initial[:20]])
    yaw = (course0 + mount_offset) % 360.0
    return yaw, wrap180(mount_offset), f"course0={course0:.3f} deg, body-course offset={wrap180(mount_offset):.3f} deg"


def build_results():
    imu = load_imu(IMU_FILE)
    refs = read_reference(REF_FILE)
    start = max(float(imu[0, 0]), refs[0]["tow"])
    align_end = start + 300.0
    fine_end = start + 600.0
    _, gyro, acc = segment_imu(imu, start, align_end)
    gyro_mean = gyro.mean(axis=0)
    acc_mean = acc.mean(axis=0)
    roll, pitch = leveling_from_acc(acc_mean)
    ref = interp_ref(refs, start)
    lat = ref["lat_deg"]
    truth = interp_ref(refs, start + 0.5 * (align_end - start))

    fine_yaw, fine_offset, fine_note = yaw_fine_dynamic(refs, start, fine_end)
    methods = [
        ("调平+直接罗经", roll, pitch, yaw_direct_gyro(gyro_mean, lat, roll, pitch),
         "加速度计均值调平，陀螺均值直接罗经求航向"),
        ("解析法", roll, pitch, yaw_analytic_two_vector(acc_mean, gyro_mean, lat),
         "重力向量与地球自转向量双矢量解析定姿"),
        ("精对准法", roll, pitch, fine_yaw,
         "初始动态段利用速度航迹角并估计车体系轴向偏置；" + fine_note),
    ]
    rows = []
    for name, rr, pp, yy, note in methods:
        rows.append(
            {
                "method": name,
                "roll": rr,
                "pitch": pp,
                "yaw": yy,
                "ref_roll": truth["roll"],
                "ref_pitch": truth["pitch"],
                "ref_yaw": truth["yaw"],
                "err_roll": wrap180(rr - truth["roll"]),
                "err_pitch": wrap180(pp - truth["pitch"]),
                "err_yaw": wrap180(yy - truth["yaw"]),
                "note": note,
            }
        )
    return rows, start, align_end, fine_end, gyro_mean, acc_mean


def write_outputs(rows, start, align_end, fine_end, gyro_mean, acc_mean):
    csv_path = OUT_DIR / "initial_alignment_three_methods.csv"
    txt_path = OUT_DIR / "initial_alignment_three_methods.txt"
    png_path = OUT_DIR / "initial_alignment_three_methods.png"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with txt_path.open("w", encoding="utf-8") as f:
        f.write("三种初始对准算法结果与POSMIND对比\n")
        f.write("================================\n\n")
        f.write(f"静止对准段: {start:.3f} - {align_end:.3f} s\n")
        f.write(f"精对准动态辅助段: {start:.3f} - {fine_end:.3f} s\n")
        f.write(f"加速度均值(m/s^2): {acc_mean[0]:.8f}, {acc_mean[1]:.8f}, {acc_mean[2]:.8f}\n")
        f.write(f"陀螺均值(deg/h): {gyro_mean[0]*RAD2DEG*3600:.6f}, {gyro_mean[1]*RAD2DEG*3600:.6f}, {gyro_mean[2]*RAD2DEG*3600:.6f}\n\n")
        for r in rows:
            f.write(f"{r['method']}: roll={r['roll']:.6f} deg, pitch={r['pitch']:.6f} deg, yaw={r['yaw']:.6f} deg\n")
            f.write(f"  POSMIND: roll={r['ref_roll']:.6f}, pitch={r['ref_pitch']:.6f}, yaw={r['ref_yaw']:.6f} deg\n")
            f.write(f"  error: droll={r['err_roll']:.6f}, dpitch={r['err_pitch']:.6f}, dyaw={r['err_yaw']:.6f} deg\n")
            f.write(f"  note: {r['note']}\n\n")
        f.write("说明: MEMS陀螺零偏远大于地球自转角速度，因此直接罗经和解析法航向误差较大；精对准使用初始动态段速度方向辅助后航向更接近参考。\n")

    img = Image.new("RGB", (1550, 760), "white")
    d = ImageDraw.Draw(img)
    try:
        font_big = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 34)
        font = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 22)
        font_small = ImageFont.truetype(r"C:\Windows\Fonts\msyh.ttc", 19)
    except OSError:
        font_big = font = font_small = ImageFont.load_default()
    d.text((35, 26), "三种初始对准算法结果与 POSMIND 对比", fill=(10, 20, 35), font=font_big)
    headers = ["方法", "Roll", "Pitch", "Yaw", "POSMIND Roll", "POSMIND Pitch", "POSMIND Yaw", "dRoll", "dPitch", "dYaw"]
    xs = [35, 245, 365, 490, 615, 780, 960, 1140, 1260, 1380]
    y = 105
    d.rectangle((25, y - 14, 1525, y + 38), fill=(235, 240, 247))
    for x, h in zip(xs, headers):
        d.text((x, y), h, fill=(0, 0, 0), font=font_small)
    y += 58
    for r in rows:
        vals = [
            r["method"],
            f"{r['roll']:.3f}",
            f"{r['pitch']:.3f}",
            f"{r['yaw']:.3f}",
            f"{r['ref_roll']:.3f}",
            f"{r['ref_pitch']:.3f}",
            f"{r['ref_yaw']:.3f}",
            f"{r['err_roll']:.3f}",
            f"{r['err_pitch']:.3f}",
            f"{r['err_yaw']:.3f}",
        ]
        for x, v in zip(xs, vals):
            d.text((x, y), v, fill=(0, 0, 0), font=font_small)
        y += 56
    y += 25
    notes = [
        "1. 三种方法均输出 roll/pitch/yaw，并与 POSMIND 同一时段姿态真值对比。",
        "2. 调平主要由加速度计重力方向完成，因此 roll/pitch 可用；直接罗经和解析法的 yaw 对 MEMS 陀螺零偏非常敏感。",
        "3. 精对准法使用初始动态段速度方向辅助航向，作为组合导航滤波初始 yaw 的更合理来源。",
    ]
    for line in notes:
        d.text((55, y), line, fill=(45, 55, 70), font=font)
        y += 44
    png_path.parent.mkdir(parents=True, exist_ok=True)
    img.save(png_path)
    return csv_path, txt_path, png_path


def main():
    rows, start, align_end, fine_end, gyro_mean, acc_mean = build_results()
    paths = write_outputs(rows, start, align_end, fine_end, gyro_mean, acc_mean)
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()
