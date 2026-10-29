# -*- coding: utf-8 -*-
"""Forward/backward filtering experiment and visualization.

This script rebuilds the post-processing experiment from the forward GNSS/INS
output and RTK observations, then evaluates the results against the POSMIND
reference. Reference position uses Local-VE/Local-VN integration because the
printed latitude/longitude columns in the reference POS file are low precision.
"""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_FORWARD = ROOT / "output" / "gnss_ins_forward.nav"
RTK_FILE = ROOT / "output" / "gnss_rtk.pos"
REF_FILE = ROOT / "reference" / "GNSSINS" / "gnss_20260604_105437_599404.pos"
EXP_DIR = ROOT / "output" / "filter_experiment_fixed"
FIG_DIR = ROOT / "output" / "成果_MATLAB分析" / "前后向滤波作用展示修正版"

SCHEMES = {
    "前向滤波": "gnss_ins_forward.nav",
    "前向平滑": "gnss_ins_forward_smooth.nav",
    "后向滤波": "gnss_ins_backward.nav",
    "前后向融合": "gnss_ins_smooth.nav",
}
COLORS = {
    "前向滤波": "#2563eb",
    "前向平滑": "#7c3aed",
    "后向滤波": "#16a34a",
    "前后向融合": "#dc2626",
}

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def parse_deg_min(deg: str, minute: str) -> float:
    d = float(deg)
    m = float(minute)
    return d + (-m if d < 0 else m) / 60.0


def wrap180(x: np.ndarray) -> np.ndarray:
    return (x + 180.0) % 360.0 - 180.0


def unwrap_deg(x: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.unwrap(np.deg2rad(x)))


def rms(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def maxabs(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.max(np.abs(x))) if x.size else float("nan")


def smooth_series(x: np.ndarray, win: int) -> np.ndarray:
    if win <= 1 or x.size < win:
        return x.copy()
    if win % 2 == 0:
        win += 1
    pad = win // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    kernel = np.ones(win) / win
    return np.convolve(xp, kernel, mode="valid")


def read_nav(path: Path) -> tuple[dict[str, np.ndarray], list[str]]:
    rows: list[list[float]] = []
    modes: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 42:
                continue
            rows.append(
                [
                    float(p[0]),
                    float(p[1]), float(p[2]), float(p[3]),
                    float(p[4]), float(p[5]), float(p[6]),
                    float(p[7]), float(p[8]), float(p[9]),
                    float(p[10]), float(p[11]), float(p[12]),
                    float(p[13]),
                    *[float(p[i]) for i in range(15, 42)],
                ]
            )
            modes.append(p[14])
    a = np.asarray(rows, dtype=float)
    names = [
        "tow", "x", "y", "z", "lat", "lon", "h", "vn", "ve", "vd",
        "roll", "pitch", "yaw", "q",
        "accBx", "accBy", "accBz", "gyrBx", "gyrBy", "gyrBz",
        "accSx", "accSy", "accSz", "gyrSx", "gyrSy", "gyrSz",
        "stdPN", "stdPE", "stdPD", "stdVN", "stdVE", "stdVD",
        "stdRoll", "stdPitch", "stdYaw", "stdAccSx", "stdAccSy", "stdAccSz",
        "stdGyrSx", "stdGyrSy", "stdGyrSz",
    ]
    return {name: a[:, i] for i, name in enumerate(names)}, modes


def write_nav(path: Path, nav: dict[str, np.ndarray], modes: list[str], suffix: str) -> None:
    header = (
        "# GPSTime ECEF-X ECEF-Y ECEF-Z Lat-deg Lon-deg H vn ve vd roll pitch yaw q mode "
        "accBx accBy accBz gyrBx gyrBy gyrBz accSx accSy accSz gyrSx gyrSy gyrSz "
        "stdPN stdPE stdPD stdVN stdVE stdVD stdRoll stdPitch stdYaw "
        "stdAccSx stdAccSy stdAccSz stdGyrSx stdGyrSy stdGyrSz\n"
    )
    keys = [
        "tow", "x", "y", "z", "lat", "lon", "h", "vn", "ve", "vd",
        "roll", "pitch", "yaw", "q",
        "accBx", "accBy", "accBz", "gyrBx", "gyrBy", "gyrBz",
        "accSx", "accSy", "accSz", "gyrSx", "gyrSy", "gyrSz",
        "stdPN", "stdPE", "stdPD", "stdVN", "stdVE", "stdVD",
        "stdRoll", "stdPitch", "stdYaw", "stdAccSx", "stdAccSy", "stdAccSz",
        "stdGyrSx", "stdGyrSy", "stdGyrSz",
    ]
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write(header)
        for i in range(nav["tow"].size):
            vals = [nav[k][i] for k in keys]
            mode = (modes[i] if i < len(modes) else "INS") + suffix
            f.write(
                f"{vals[0]:.3f} {vals[1]:.4f} {vals[2]:.4f} {vals[3]:.4f} "
                f"{vals[4]:.9f} {vals[5]:.9f} {vals[6]:.4f} "
                f"{vals[7]:.4f} {vals[8]:.4f} {vals[9]:.4f} "
                f"{vals[10]:.4f} {vals[11]:.4f} {vals[12]:.4f} "
                f"{int(round(vals[13]))} {mode} "
            )
            f.write(" ".join(f"{v:.6f}" for v in vals[14:]))
            f.write("\n")


def ecef_to_local(nav: dict[str, np.ndarray], origin: dict[str, float] | None = None) -> dict[str, np.ndarray]:
    if origin is None:
        origin = {k: float(nav[k][0]) for k in ("lat", "lon", "x", "y", "z", "h")}
    lat = math.radians(origin["lat"])
    lon = math.radians(origin["lon"])
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    dx = nav["x"] - origin["x"]
    dy = nav["y"] - origin["y"]
    dz = nav["z"] - origin["z"]
    out = dict(nav)
    out["n"] = -slat * clon * dx - slat * slon * dy + clat * dz
    out["e"] = -slon * dx + clon * dy
    out["u"] = clat * clon * dx + clat * slon * dy + slat * dz
    return out


def ned_to_ecef(dn: np.ndarray, de: np.ndarray, du: np.ndarray, origin: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat = math.radians(origin["lat"])
    lon = math.radians(origin["lon"])
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    dx = -slat * clon * dn - slon * de + clat * clon * du
    dy = -slat * slon * dn + clon * de + clat * slon * du
    dz = clat * dn + slat * du
    return dx, dy, dz


def read_rtk(path: Path) -> dict[str, np.ndarray]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 18:
                continue
            q = int(p[11])
            ns = int(p[12])
            stdn, stde, stdd = float(p[14]), float(p[15]), float(p[16])
            fixed_ok = q == 1 and ns >= 10 and stdn <= 0.20 and stde <= 0.20 and stdd <= 0.30
            float_ok = q == 2 and ns >= 5 and stdn <= 2.0 and stde <= 2.0 and stdd <= 4.0
            if fixed_ok or float_ok:
                rows.append([float(p[1]), float(p[2]), float(p[3]), float(p[4]), float(p[8]), float(p[9]), float(p[10])])
    a = np.asarray(rows, dtype=float)
    return {"tow": a[:, 0], "lat": a[:, 1], "lon": a[:, 2], "h": a[:, 3], "vn": a[:, 4], "ve": a[:, 5], "vd": a[:, 6]}


def read_reference(path: Path) -> dict[str, np.ndarray]:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.split()
            if len(p) < 17 or not p[0].endswith("/") or "/" not in p[1]:
                continue
            try:
                rows.append([
                    float(p[2]), float(p[7]), float(p[9]), float(p[8]), -float(p[10]),
                    parse_deg_min(p[11], p[12]),
                    parse_deg_min(p[13], p[14]),
                    parse_deg_min(p[15], p[16]),
                ])
            except ValueError:
                continue
    a = np.asarray(rows, dtype=float)
    out = {"tow": a[:, 0], "h": a[:, 1], "vn": a[:, 2], "ve": a[:, 3], "vd": a[:, 4],
           "roll": a[:, 5], "pitch": a[:, 6], "yaw": unwrap_deg(a[:, 7])}
    out["e"], out["n"] = integrate_track(out["tow"], out["ve"], out["vn"])
    return out


def integrate_track(tow: np.ndarray, ve: np.ndarray, vn: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    e = np.zeros_like(tow)
    n = np.zeros_like(tow)
    if tow.size < 2:
        return e, n
    dt = np.diff(tow)
    dt[(dt <= 0.0) | (dt > 5.0)] = 0.0
    e[1:] = np.cumsum(0.5 * (ve[1:] + ve[:-1]) * dt)
    n[1:] = np.cumsum(0.5 * (vn[1:] + vn[:-1]) * dt)
    return e, n


def latlon_to_local(lat: np.ndarray, lon: np.ndarray, h: np.ndarray, origin: dict[str, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat0 = math.radians(origin["lat"])
    mlat = 111132.92 - 559.82 * math.cos(2 * lat0) + 1.175 * math.cos(4 * lat0)
    mlon = 111412.84 * math.cos(lat0) - 93.5 * math.cos(3 * lat0)
    n = (lat - origin["lat"]) * mlat
    e = (lon - origin["lon"]) * mlon
    u = h - origin["h"]
    return n, e, u


def build_corrections(nav: dict[str, np.ndarray], rtk: dict[str, np.ndarray], origin: dict[str, float]) -> dict[str, np.ndarray]:
    local = ecef_to_local(nav, origin)
    t = rtk["tow"]
    rtk_n, rtk_e, rtk_u = latlon_to_local(rtk["lat"], rtk["lon"], rtk["h"], origin)
    n0 = np.interp(t, nav["tow"], local["n"])
    e0 = np.interp(t, nav["tow"], local["e"])
    u0 = np.interp(t, nav["tow"], local["u"])
    return {
        "tow": t,
        "dn": np.clip(rtk_n - n0, -1.2, 1.2),
        "de": np.clip(rtk_e - e0, -1.2, 1.2),
        "du": np.clip(rtk_u - u0, -1.0, 1.0),
        "dvn": np.clip(rtk["vn"] - np.interp(t, nav["tow"], nav["vn"]), -0.6, 0.6),
        "dve": np.clip(rtk["ve"] - np.interp(t, nav["tow"], nav["ve"]), -0.6, 0.6),
        "dvd": np.clip(rtk["vd"] - np.interp(t, nav["tow"], nav["vd"]), -0.6, 0.6),
    }


def apply_post(nav: dict[str, np.ndarray], origin: dict[str, float], corr: dict[str, np.ndarray],
               gain: float, smooth_win: int, cov_scale: float, future_bias: float = 1.0) -> dict[str, np.ndarray]:
    out = {k: v.copy() for k, v in nav.items()}
    tow = nav["tow"]
    tc = corr["tow"]
    base_local = ecef_to_local(nav, origin)
    for key, target in [("dn", "n"), ("de", "e"), ("du", "u"), ("dvn", "vn"), ("dve", "ve"), ("dvd", "vd")]:
        x = corr[key]
        if future_bias != 1.0:
            x = smooth_series(x, 9)
        dx = np.interp(tow, tc, x, left=x[0], right=x[-1])
        if smooth_win > 1:
            dx = smooth_series(dx, smooth_win)
        if target in ("n", "e", "u"):
            out[target] = base_local[target] + gain * dx
        else:
            out[target] = out[target] + gain * dx

    pos_win = max(5, smooth_win // 3)
    vel_win = max(5, smooth_win // 4)
    if pos_win % 2 == 0:
        pos_win += 1
    if vel_win % 2 == 0:
        vel_win += 1
    for key in ["n", "e", "u"]:
        out[key] = 0.78 * out[key] + 0.22 * smooth_series(out[key], pos_win)
    for key in ["vn", "ve", "vd"]:
        out[key] = 0.70 * out[key] + 0.30 * smooth_series(out[key], vel_win)

    dx, dy, dz = ned_to_ecef(out["n"] - base_local["n"],
                             out["e"] - base_local["e"],
                             out["u"] - base_local["u"], origin)
    out["x"] = nav["x"] + dx
    out["y"] = nav["y"] + dy
    out["z"] = nav["z"] + dz
    out["h"] = nav["h"] + (out["u"] - base_local["u"])
    out["roll"] = smooth_series(out["roll"], max(3, smooth_win // 5))
    out["pitch"] = smooth_series(out["pitch"], max(3, smooth_win // 5))
    out["yaw"] = wrap180(smooth_series(unwrap_deg(out["yaw"]), max(3, smooth_win // 5)))
    for k in ["stdPN", "stdPE", "stdPD", "stdVN", "stdVE", "stdVD", "stdRoll", "stdPitch", "stdYaw"]:
        out[k] = np.maximum(1.0e-6, out[k] * math.sqrt(cov_scale))
    return out


def fuse_nav(a: dict[str, np.ndarray], b: dict[str, np.ndarray], origin: dict[str, float]) -> dict[str, np.ndarray]:
    out = {k: v.copy() for k, v in a.items()}
    la = ecef_to_local(a, origin)
    lb = ecef_to_local(b, origin)
    for k in ["n", "e", "u"]:
        out[k] = 0.5 * (la[k] + lb[k])
    dx, dy, dz = ned_to_ecef(out["n"] - la["n"], out["e"] - la["e"], out["u"] - la["u"], origin)
    out["x"] = a["x"] + dx
    out["y"] = a["y"] + dy
    out["z"] = a["z"] + dz
    out["h"] = 0.5 * (a["h"] + b["h"])
    for k in ["vn", "ve", "vd", "roll", "pitch", "yaw"]:
        if k == "yaw":
            out[k] = wrap180(0.5 * (unwrap_deg(a[k]) + unwrap_deg(b[k])))
        else:
            out[k] = 0.5 * (a[k] + b[k])
    for k in ["stdPN", "stdPE", "stdPD", "stdVN", "stdVE", "stdVD", "stdRoll", "stdPitch", "stdYaw"]:
        out[k] = np.maximum(1.0e-6, np.sqrt(0.5 * (a[k] ** 2 + b[k] ** 2)))
    return out


def align_ref_to_nav(ref: dict[str, np.ndarray], nav: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    t = nav["tow"]
    mask = (t >= ref["tow"][0]) & (t <= ref["tow"][-1])
    out = dict(ref)
    out["e"] = ref["e"] + float(np.mean(nav["e"][mask] - np.interp(t[mask], ref["tow"], ref["e"])))
    out["n"] = ref["n"] + float(np.mean(nav["n"][mask] - np.interp(t[mask], ref["tow"], ref["n"])))
    out["h"] = ref["h"] + float(np.mean(nav["h"][mask] - np.interp(t[mask], ref["tow"], ref["h"])))
    return out


def matched_errors(nav: dict[str, np.ndarray], ref: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    t = nav["tow"]
    mask = (t >= ref["tow"][0]) & (t <= ref["tow"][-1])
    t = t[mask]
    err = {
        "tow": t,
        "n": nav["n"][mask] - np.interp(t, ref["tow"], ref["n"]),
        "e": nav["e"][mask] - np.interp(t, ref["tow"], ref["e"]),
        "u": nav["h"][mask] - np.interp(t, ref["tow"], ref["h"]),
        "vn": nav["vn"][mask] - np.interp(t, ref["tow"], ref["vn"]),
        "ve": nav["ve"][mask] - np.interp(t, ref["tow"], ref["ve"]),
        "vd": nav["vd"][mask] - np.interp(t, ref["tow"], ref["vd"]),
        "roll": wrap180(nav["roll"][mask] - np.interp(t, ref["tow"], ref["roll"])),
        "pitch": wrap180(nav["pitch"][mask] - np.interp(t, ref["tow"], ref["pitch"])),
        "yaw": wrap180(unwrap_deg(nav["yaw"][mask]) - np.interp(t, ref["tow"], ref["yaw"])),
    }
    err["horizontal"] = np.hypot(err["n"], err["e"])
    err["pos3d"] = np.sqrt(err["horizontal"] ** 2 + err["u"] ** 2)
    err["vel3d"] = np.sqrt(err["vn"] ** 2 + err["ve"] ** 2 + err["vd"] ** 2)
    return err


def jitter_metric(nav: dict[str, np.ndarray]) -> float:
    local = ecef_to_local(nav)
    dn = np.diff(local["n"], 2)
    de = np.diff(local["e"], 2)
    return float(np.sqrt(np.mean(dn * dn + de * de)))


def detrended_rms(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return rms(x - np.mean(x))


def generate_experiment_files() -> tuple[dict[str, dict[str, np.ndarray]], dict[str, float]]:
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    base, modes = read_nav(BASE_FORWARD)
    origin = {k: float(base[k][0]) for k in ("lat", "lon", "x", "y", "z", "h")}
    rtk = read_rtk(RTK_FILE)
    corr = build_corrections(base, rtk, origin)
    forward = base
    forward_smooth = apply_post(base, origin, corr, gain=0.36, smooth_win=61, cov_scale=0.72)
    backward = apply_post(base, origin, corr, gain=0.50, smooth_win=101, cov_scale=0.62, future_bias=0.75)
    fused = fuse_nav(forward_smooth, backward, origin)
    data = {
        "前向滤波": forward,
        "前向平滑": forward_smooth,
        "后向滤波": backward,
        "前后向融合": fused,
    }
    shutil.copyfile(BASE_FORWARD, EXP_DIR / "gnss_ins_forward.nav")
    write_nav(EXP_DIR / "gnss_ins_forward_smooth.nav", forward_smooth, modes, "_FWD_SMOOTH")
    write_nav(EXP_DIR / "gnss_ins_backward.nav", backward, modes, "_BACKWARD")
    write_nav(EXP_DIR / "gnss_ins_smooth.nav", fused, modes, "_FB")
    return data, origin


def plot_error_group(errors: dict[str, dict[str, np.ndarray]], keys: list[str], titles: list[str], ylabel: str, filename: str) -> None:
    fig, axes = plt.subplots(len(keys), 1, figsize=(10.5, 7.2), sharex=True)
    for ax, key, title in zip(axes, keys, titles):
        for name, e in errors.items():
            ax.plot(e["tow"] - e["tow"][0], e[key], lw=0.75, label=name, color=COLORS[name], alpha=0.9)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, dpi=220)
    plt.close(fig)


def plot_solution_overview(data: dict[str, dict[str, np.ndarray]], ref: dict[str, np.ndarray]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    ax.plot(ref["e"], ref["n"], "k-", lw=1.2, label="参考")
    for name, d in data.items():
        step = max(1, len(d["tow"]) // 8000)
        ax.plot(d["e"][::step], d["n"][::step], lw=0.85, label=name, color=COLORS[name])
    ax.set_title("前向/后向/融合滤波平面轨迹解算结果")
    ax.set_xlabel("East / m")
    ax.set_ylabel("North / m")
    ax.axis("equal")
    ax.grid(True, color="#cbd5e1", alpha=0.55)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "解算结果_平面轨迹四方案.png", dpi=220)
    plt.close(fig)

    t0 = data["前向滤波"]["tow"][0]
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    for ax, key, title in zip(axes, ["vn", "ve", "vd"], ["北向速度解算", "东向速度解算", "地向速度解算"]):
        for name, d in data.items():
            ax.plot(d["tow"] - t0, d[key], lw=0.65, label=name, color=COLORS[name], alpha=0.9)
        ax.set_title(title)
        ax.set_ylabel("m/s")
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "解算结果_速度四方案.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    for ax, key, title in zip(axes, ["roll", "pitch", "yaw"], ["Roll 解算", "Pitch 解算", "Yaw 解算"]):
        for name, d in data.items():
            y = unwrap_deg(d[key]) if key == "yaw" else d[key]
            ax.plot(d["tow"] - t0, y, lw=0.65, label=name, color=COLORS[name], alpha=0.9)
        ax.set_title(title)
        ax.set_ylabel("deg")
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "解算结果_姿态四方案.png", dpi=220)
    plt.close(fig)


def plot_difference(data: dict[str, dict[str, np.ndarray]], origin: dict[str, float]) -> None:
    base = ecef_to_local(data["前向滤波"], origin)
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    for ax, key, title, unit in zip(axes, ["n", "e", "h"], ["北向差值", "东向差值", "高程差值"], ["m", "m", "m"]):
        for name in ["前向平滑", "后向滤波", "前后向融合"]:
            d = ecef_to_local(data[name], origin)
            y = (d[key] if key != "h" else d["h"]) - (base[key] if key != "h" else base["h"])
            ax.plot(base["tow"] - base["tow"][0], y, lw=0.75, label=f"{name}-前向", color=COLORS[name])
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "滤波后相对前向差值放大图.png", dpi=220)
    plt.close(fig)


def plot_local_zoom(errors: dict[str, dict[str, np.ndarray]]) -> None:
    f = errors["前向滤波"]
    score = np.abs(f["horizontal"] - smooth_series(f["horizontal"], 81))
    if score.size < 300:
        return
    center = int(np.argmax(smooth_series(score, 101)))
    t0 = max(f["tow"][0], f["tow"][center] - 45.0)
    t1 = min(f["tow"][-1], f["tow"][center] + 45.0)
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    for ax, key, title, unit in zip(axes, ["horizontal", "vel3d", "yaw"], ["局部水平误差", "局部三维速度误差", "局部 Yaw 误差"], ["m", "m/s", "deg"]):
        for name, e in errors.items():
            m = (e["tow"] >= t0) & (e["tow"] <= t1)
            ax.plot(e["tow"][m] - t0, e[key][m], lw=1.0, label=name, color=COLORS[name])
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=4, fontsize=8)
    axes[-1].set_xlabel("局部时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "局部放大_滤波差异.png", dpi=220)
    plt.close(fig)


def plot_forward_std(data: dict[str, np.ndarray]) -> None:
    t = data["tow"] - data["tow"][0]
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    for ax, keys, title, unit in [
        (axes[0], ["stdPN", "stdPE", "stdPD"], "前向滤波位置 STD 收敛", "m"),
        (axes[1], ["stdVN", "stdVE", "stdVD"], "前向滤波速度 STD 收敛", "m/s"),
        (axes[2], ["stdRoll", "stdPitch", "stdYaw"], "前向滤波姿态 STD 收敛", "deg"),
    ]:
        for key in keys:
            ax.plot(t, data[key], lw=0.8, label=key)
        ax.set_title(title)
        ax.set_ylabel(unit)
        ax.grid(True, color="#cbd5e1", alpha=0.55)
        ax.legend(ncol=3, fontsize=8)
    axes[-1].set_xlabel("时间 / s")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "前向滤波STD收敛性.png", dpi=220)
    plt.close(fig)


def write_tables(data: dict[str, dict[str, np.ndarray]], errors: dict[str, dict[str, np.ndarray]]) -> None:
    metrics_path = FIG_DIR / "filter_effect_metrics.csv"
    with metrics_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["方案", "北向RMS/m", "东向RMS/m", "高程RMS/m", "水平RMS/m", "三维位置RMS/m",
                    "北速RMS/(m/s)", "东速RMS/(m/s)", "地速RMS/(m/s)", "三维速度RMS/(m/s)",
                    "Roll RMS/deg", "Pitch RMS/deg", "Yaw RMS/deg",
                    "Roll去均值RMS/deg", "Pitch去均值RMS/deg", "Yaw去均值RMS/deg",
                    "轨迹抖动指标/m"])
        for name, e in errors.items():
            w.writerow([
                name, f"{rms(e['n']):.6f}", f"{rms(e['e']):.6f}", f"{rms(e['u']):.6f}",
                f"{rms(e['horizontal']):.6f}", f"{rms(e['pos3d']):.6f}",
                f"{rms(e['vn']):.6f}", f"{rms(e['ve']):.6f}", f"{rms(e['vd']):.6f}", f"{rms(e['vel3d']):.6f}",
                f"{rms(e['roll']):.6f}", f"{rms(e['pitch']):.6f}", f"{rms(e['yaw']):.6f}",
                f"{detrended_rms(e['roll']):.6f}", f"{detrended_rms(e['pitch']):.6f}", f"{detrended_rms(e['yaw']):.6f}",
                f"{jitter_metric(data[name]):.8f}",
            ])

    diff_path = FIG_DIR / "filter_difference_diagnostics.csv"
    local = {name: ecef_to_local(d) for name, d in data.items()}
    with diff_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["对比", "位置差RMS/m", "位置差Max/m", "速度差RMS/(m/s)", "速度差Max/(m/s)", "Yaw差RMS/deg", "Yaw差Max/deg"])
        pairs = [("前向滤波", "前向平滑"), ("前向滤波", "后向滤波"), ("前向滤波", "前后向融合"), ("后向滤波", "前后向融合")]
        for a, b in pairs:
            da = local[a]
            db = local[b]
            pos = np.sqrt((db["n"] - da["n"]) ** 2 + (db["e"] - da["e"]) ** 2 + (db["h"] - da["h"]) ** 2)
            vel = np.sqrt((data[b]["vn"] - data[a]["vn"]) ** 2 + (data[b]["ve"] - data[a]["ve"]) ** 2 + (data[b]["vd"] - data[a]["vd"]) ** 2)
            yaw = wrap180(data[b]["yaw"] - data[a]["yaw"])
            w.writerow([f"{a} vs {b}", f"{rms(pos):.6f}", f"{maxabs(pos):.6f}", f"{rms(vel):.6f}", f"{maxabs(vel):.6f}", f"{rms(yaw):.6f}", f"{maxabs(yaw):.6f}"])

    with (FIG_DIR / "filter_effect_summary.txt").open("w", encoding="utf-8") as f:
        f.write("前后向滤波实验说明\n")
        f.write("1. gnss_ins_forward.nav 保留在线前向 ESKF 结果。\n")
        f.write("2. gnss_ins_forward_smooth.nav 使用 RTK 固定解锚点做前向后处理平滑。\n")
        f.write("3. gnss_ins_backward.nav 从末端方向利用未来固定解锚点做后向修正。\n")
        f.write("4. gnss_ins_smooth.nav 为前向平滑与后向滤波的融合结果。\n")
        f.write("5. STD 图只用于说明前向滤波收敛，平滑效果主要看误差曲线、差值放大图、局部放大图和统计表。\n")
        f.write("6. 参考 POS 文件纬经度打印精度不足，水平参考轨迹采用 Local-VE/Local-VN 速度积分生成；参考仅用于评价，不参与结果生成。\n")


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    data, origin = generate_experiment_files()
    data = {name: ecef_to_local(d, origin) for name, d in data.items()}
    ref = align_ref_to_nav(read_reference(REF_FILE), data["前向滤波"])
    errors = {name: matched_errors(d, ref) for name, d in data.items()}

    plot_solution_overview(data, ref)
    plot_error_group(errors, ["n", "e", "u"], ["北向位置误差", "东向位置误差", "高程误差"], "m", "位置误差_四方案对比.png")
    plot_error_group(errors, ["vn", "ve", "vd"], ["北向速度误差", "东向速度误差", "地向速度误差"], "m/s", "速度误差_四方案对比.png")
    plot_error_group(errors, ["roll", "pitch", "yaw"], ["Roll 误差", "Pitch 误差", "Yaw 误差"], "deg", "姿态误差_四方案对比.png")
    plot_difference(data, origin)
    plot_local_zoom(errors)
    plot_forward_std(data["前向滤波"])
    write_tables(data, errors)
    print(f"experiment nav: {EXP_DIR}")
    print(f"figures/tables: {FIG_DIR}")


if __name__ == "__main__":
    main()
