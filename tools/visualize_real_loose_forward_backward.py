# -*- coding: utf-8 -*-
"""Generate loose GNSS/INS figures with reference overlays and difference panels.

The loose forward, backward and smoothed tracks are often very close at the
hundreds-of-metres scale.  The report therefore uses faceted plots plus
difference series instead of relying on one overplotted red curve.
"""

from __future__ import annotations

import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "output" / "成果_REAL_松组合真实算法"
REF_FILE = ROOT / "reference" / "GNSSINS" / "gnss_20260604_105437_599404.pos"
OUT = ROOT / "output" / "成果_REAL_全部可视化" / "松组合前后向真实算法"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "前向": BASE / "gnss_ins_forward.nav",
    "前向平滑": BASE / "gnss_ins_forward_smooth.nav",
    "后向": BASE / "gnss_ins_backward.nav",
    "前后向平滑": BASE / "gnss_ins_smooth.nav",
}

COLORS = {
    "参考": "#111827",
    "前向": "#9ca3af",
    "前向平滑": "#2563eb",
    "后向": "#16a34a",
    "前后向平滑": "#dc2626",
}


def load_font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    for path in [r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf", r"C:\Windows\Fonts\simsun.ttc"]:
        if Path(path).exists():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


F_TITLE = load_font(30)
F_SUB = load_font(21)
F = load_font(17)
F_SMALL = load_font(14)


def read_nav(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            rows.append(
                {
                    "tow": float(p[0]),
                    "x": float(p[1]),
                    "y": float(p[2]),
                    "z": float(p[3]),
                    "lat": float(p[4]),
                    "lon": float(p[5]),
                    "h": float(p[6]),
                    "vn": float(p[7]),
                    "ve": float(p[8]),
                    "vd": float(p[9]),
                    "roll": float(p[10]),
                    "pitch": float(p[11]),
                    "yaw": float(p[12]),
                }
            )
    return rows


def parse_deg_min(deg: str, minute: str) -> float:
    d = float(deg)
    m = float(minute)
    return d + (-m if d < 0.0 else m) / 60.0


def read_reference(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            p = line.split()
            if len(p) < 17 or "/" not in p[0]:
                continue
            try:
                rows.append(
                    {
                        "tow": float(p[2]),
                        "h": float(p[7]),
                        "ve": float(p[8]),
                        "vn": float(p[9]),
                        "vd": -float(p[10]),
                        "roll": parse_deg_min(p[11], p[12]),
                        "pitch": parse_deg_min(p[13], p[14]),
                        "yaw": parse_deg_min(p[15], p[16]),
                    }
                )
            except ValueError:
                continue
    if not rows:
        return rows
    rows[0]["e"] = 0.0
    rows[0]["n"] = 0.0
    rows[0]["u"] = 0.0
    for i in range(1, len(rows)):
        dt = rows[i]["tow"] - rows[i - 1]["tow"]
        if dt <= 0 or dt > 5:
            dt = 1.0
        rows[i]["e"] = rows[i - 1]["e"] + 0.5 * (rows[i]["ve"] + rows[i - 1]["ve"]) * dt
        rows[i]["n"] = rows[i - 1]["n"] + 0.5 * (rows[i]["vn"] + rows[i - 1]["vn"]) * dt
        rows[i]["u"] = rows[i]["h"] - rows[0]["h"]
    for r in rows:
        r["speed"] = math.sqrt(r["vn"] ** 2 + r["ve"] ** 2 + r["vd"] ** 2)
    return rows


def to_local(rows: list[dict[str, float]], origin: dict[str, float]) -> list[dict[str, float]]:
    lat = math.radians(origin["lat"])
    lon = math.radians(origin["lon"])
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    x0, y0, z0 = origin["x"], origin["y"], origin["z"]
    local = []
    for r in rows:
        dx, dy, dz = r["x"] - x0, r["y"] - y0, r["z"] - z0
        q = dict(r)
        q["n"] = -slat * clon * dx - slat * slon * dy + clat * dz
        q["e"] = -slon * dx + clon * dy
        q["u"] = clat * clon * dx + clat * slon * dy + slat * dz
        q["speed"] = math.sqrt(r["vn"] ** 2 + r["ve"] ** 2 + r["vd"] ** 2)
        local.append(q)
    return local


def align_reference(ref: list[dict[str, float]], smooth: list[dict[str, float]]) -> list[dict[str, float]]:
    if not ref:
        return []
    offsets = []
    ismooth = 0
    for r in ref:
        while ismooth + 1 < len(smooth) and abs(smooth[ismooth + 1]["tow"] - r["tow"]) < abs(smooth[ismooth]["tow"] - r["tow"]):
            ismooth += 1
        s = smooth[ismooth]
        if abs(s["tow"] - r["tow"]) <= 0.55:
            offsets.append((s["e"] - r["e"], s["n"] - r["n"], s["h"] - r["h"]))
    if offsets:
        de = sum(v[0] for v in offsets) / len(offsets)
        dn = sum(v[1] for v in offsets) / len(offsets)
        dh = sum(v[2] for v in offsets) / len(offsets)
    else:
        target = smooth[0]["tow"]
        idx = min(range(len(ref)), key=lambda i: abs(ref[i]["tow"] - target))
        de = smooth[0]["e"] - ref[idx]["e"]
        dn = smooth[0]["n"] - ref[idx]["n"]
        dh = smooth[0]["h"] - ref[idx]["h"]
    aligned = []
    for r in ref:
        q = dict(r)
        q["e"] = r["e"] + de
        q["n"] = r["n"] + dn
        q["h"] = r["h"] + dh
        q["u"] = q["h"] - smooth[0]["h"]
        aligned.append(q)
    return aligned


def thin(rows: list[dict[str, float]], target: int = 2500) -> list[dict[str, float]]:
    step = max(1, len(rows) // target)
    return rows[::step]


def make_canvas(title: str, size: tuple[int, int] = (1400, 920)) -> tuple[Image.Image, ImageDraw.ImageDraw]:
    im = Image.new("RGB", size, "white")
    d = ImageDraw.Draw(im)
    d.text((44, 24), title, fill="#111827", font=F_TITLE)
    return im, d


def draw_box(d: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x1, y1, x2, y2 = box
    d.rectangle(box, outline="#374151", width=2)
    for i in range(1, 4):
        x = x1 + (x2 - x1) * i / 4
        y = y1 + (y2 - y1) * i / 4
        d.line((x, y1, x, y2), fill="#e5e7eb")
        d.line((x1, y, x2, y), fill="#e5e7eb")


def scaler(points: list[tuple[float, float]], box: tuple[int, int, int, int], equal: bool = False):
    x1, y1, x2, y2 = box
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if equal:
        span = max(xmax - xmin, ymax - ymin, 1e-6)
        cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
        xmin, xmax = cx - span * 0.58, cx + span * 0.58
        ymin, ymax = cy - span * 0.58, cy + span * 0.58
    else:
        padx = max(1e-6, 0.08 * (xmax - xmin))
        pady = max(1e-6, 0.08 * (ymax - ymin))
        xmin, xmax = xmin - padx, xmax + padx
        ymin, ymax = ymin - pady, ymax + pady

    def pxy(x: float, y: float) -> tuple[int, int]:
        px = x1 + (x - xmin) / (xmax - xmin) * (x2 - x1)
        py = y2 - (y - ymin) / (ymax - ymin) * (y2 - y1)
        return int(px), int(py)

    return pxy


def draw_legend(d: ImageDraw.ImageDraw, x: int, y: int, names: list[str]) -> None:
    for name in names:
        d.line((x, y, x + 42, y), fill=COLORS[name], width=4)
        d.text((x + 52, y - 12), name, fill="#111827", font=F)
        y += 29


def angle_diff(a: float, b: float) -> float:
    d = a - b
    while d > 180.0:
        d -= 360.0
    while d < -180.0:
        d += 360.0
    return d


def plot_track_facets(data: dict[str, list[dict[str, float]]], ref: list[dict[str, float]], out: Path) -> None:
    im, d = make_canvas("松组合轨迹分面对比：每幅均叠加参考轨迹", (1500, 1040))
    boxes = {
        "前向": (78, 118, 700, 480),
        "前向平滑": (800, 118, 1422, 480),
        "后向": (78, 590, 700, 952),
        "前后向平滑": (800, 590, 1422, 952),
    }
    all_pts = [(r["e"], r["n"]) for rows in data.values() for r in thin(rows)] + [(r["e"], r["n"]) for r in ref]
    for name, box in boxes.items():
        draw_box(d, box)
        d.text((box[0], box[1] - 34), name + " / 黑线为参考", fill="#111827", font=F_SUB)
        pxy = scaler(all_pts, box, equal=True)
        ref_pts = [pxy(r["e"], r["n"]) for r in ref]
        own_pts = [pxy(r["e"], r["n"]) for r in thin(data[name])]
        if len(ref_pts) > 1:
            d.line(ref_pts, fill=COLORS["参考"], width=3)
        d.line(own_pts, fill=COLORS[name], width=3)
    draw_legend(d, 1050, 980, ["参考", "前向", "前向平滑", "后向", "前后向平滑"])
    im.save(out)


def plot_track_diff(data: dict[str, list[dict[str, float]]], out: Path) -> None:
    im, d = make_canvas("松组合相对最终平滑的三维差值", (1400, 900))
    smooth = data["前后向平滑"]
    names = ["前向", "前向平滑", "后向"]
    comps = [("e", "dE / m"), ("n", "dN / m"), ("u", "dU / m")]
    for panel, (key, ylabel) in enumerate(comps):
        box = (95, 105 + panel * 250, 1305, 290 + panel * 250)
        draw_box(d, box)
        series = {}
        for name in names:
            rows = data[name]
            step = max(1, len(rows) // 3200)
            points = []
            for i in range(0, min(len(rows), len(smooth)), step):
                points.append((rows[i]["tow"] - rows[0]["tow"], rows[i][key] - smooth[i][key]))
            series[name] = points
        xs = [x for pts in series.values() for x, _ in pts]
        ys = [y for pts in series.values() for _, y in pts]
        pxy = scaler(list(zip(xs, ys)), box)
        y0 = pxy(xs[0], 0.0)[1] if min(ys) < 0.0 < max(ys) else None
        if y0:
            d.line((box[0], y0, box[2], y0), fill="#9ca3af")
        d.text((28, (box[1] + box[3]) // 2 - 12), ylabel, fill="#111827", font=F)
        for name, pts in series.items():
            d.line([pxy(x, y) for x, y in pts], fill=COLORS[name], width=2)
    draw_legend(d, 1080, 38, names)
    d.text((620, 858), "时间 / s", fill="#111827", font=F)
    im.save(out)


def plot_scalar_facets(data: dict[str, list[dict[str, float]]], ref: list[dict[str, float]],
                       key: str, title: str, ylabel: str, out: Path) -> None:
    im, d = make_canvas(title, (1500, 1040))
    boxes = {
        "前向": (78, 118, 700, 480),
        "前向平滑": (800, 118, 1422, 480),
        "后向": (78, 590, 700, 952),
        "前后向平滑": (800, 590, 1422, 952),
    }
    ref_points = [(r["tow"] - ref[0]["tow"], r[key]) for r in ref]
    for name, box in boxes.items():
        draw_box(d, box)
        d.text((box[0], box[1] - 34), name + " / 黑线为参考", fill="#111827", font=F_SUB)
        rows = thin(data[name])
        own_points = [(r["tow"] - data[name][0]["tow"], r[key]) for r in rows]
        pts_for_scale = ref_points + own_points
        pxy = scaler(pts_for_scale, box)
        d.line([pxy(x, y) for x, y in ref_points], fill=COLORS["参考"], width=2)
        d.line([pxy(x, y) for x, y in own_points], fill=COLORS[name], width=2)
        d.text((box[0] + 8, box[1] + 8), ylabel, fill="#374151", font=F_SMALL)
    draw_legend(d, 1050, 980, ["参考", "前向", "前向平滑", "后向", "前后向平滑"])
    im.save(out)


def plot_yaw_speed(data: dict[str, list[dict[str, float]]], ref: list[dict[str, float]], out: Path) -> None:
    im, d = make_canvas("松组合姿态与速度分面：每幅均叠加参考", (1500, 1120))
    names = ["前向", "前向平滑", "后向", "前后向平滑"]
    panels = [("yaw", "航向角 / deg"), ("speed", "速度 / m/s")]
    for row, name in enumerate(names):
        for col, (key, ylabel) in enumerate(panels):
            x0 = 82 + col * 710
            y0 = 112 + row * 240
            box = (x0, y0, x0 + 620, y0 + 175)
            draw_box(d, box)
            d.text((box[0], box[1] - 28), f"{name} - {ylabel}", fill="#111827", font=F)
            ref_points = [(r["tow"] - ref[0]["tow"], r[key]) for r in ref]
            rows = thin(data[name])
            own_points = [(r["tow"] - data[name][0]["tow"], r[key]) for r in rows]
            pxy = scaler(ref_points + own_points, box)
            d.line([pxy(x, y) for x, y in ref_points], fill=COLORS["参考"], width=2)
            d.line([pxy(x, y) for x, y in own_points], fill=COLORS[name], width=2)
            if key == "yaw":
                d.text((box[0] + 8, box[3] - 22), "跨 ±180 deg 会出现跳变", fill="#6b7280", font=F_SMALL)
    draw_legend(d, 1140, 972, ["参考", "前向", "前向平滑", "后向", "前后向平滑"])
    im.save(out)


def plot_enu_facets(data: dict[str, list[dict[str, float]]], ref: list[dict[str, float]], out: Path) -> None:
    im, d = make_canvas("松组合三维局部坐标序列分面：每幅均叠加参考", (1600, 1220))
    panels = [("e", "东向 / m"), ("n", "北向 / m"), ("u", "天向 / m")]
    names = ["前向", "前向平滑", "后向", "前后向平滑"]
    for row, name in enumerate(names):
        for col, (key, ylabel) in enumerate(panels):
            x0 = 78 + col * 505
            y0 = 112 + row * 260
            box = (x0, y0, x0 + 455, y0 + 190)
            draw_box(d, box)
            d.text((box[0], box[1] - 28), f"{name} - {ylabel}", fill="#111827", font=F)
            ref_points = [(r["tow"] - ref[0]["tow"], r[key]) for r in ref]
            rows = thin(data[name])
            own_points = [(r["tow"] - data[name][0]["tow"], r[key]) for r in rows]
            pxy = scaler(ref_points + own_points, box)
            d.line([pxy(x, y) for x, y in ref_points], fill=COLORS["参考"], width=2)
            d.line([pxy(x, y) for x, y in own_points], fill=COLORS[name], width=2)
    draw_legend(d, 1230, 1060, ["参考", "前向", "前向平滑", "后向", "前后向平滑"])
    im.save(out)


def main() -> None:
    raw = {name: read_nav(path) for name, path in FILES.items()}
    origin = raw["前后向平滑"][0]
    data = {name: to_local(rows, origin) for name, rows in raw.items()}
    ref = align_reference(read_reference(REF_FILE), data["前后向平滑"])

    plot_track_facets(data, ref, OUT / "10-2_真实算法松组合轨迹总览.png")
    plot_track_diff(data, OUT / "10-3_真实算法松组合轨迹差异放大.png")
    plot_scalar_facets(data, ref, "h", "松组合高程序列分面对比", "高程 / m", OUT / "10-4_真实算法松组合高程序列.png")
    plot_scalar_facets(data, ref, "speed", "松组合速度序列分面对比", "速度 / m/s", OUT / "10-5_真实算法松组合速度序列.png")
    plot_yaw_speed(data, ref, OUT / "10-6_真实算法前后向姿态速度.png")
    plot_enu_facets(data, ref, OUT / "10-7_真实算法三维轨迹平面投影.png")
    plot_track_facets(data, ref, OUT / "附加_真实算法松组合分面检查.png")
    print(OUT)


if __name__ == "__main__":
    main()
