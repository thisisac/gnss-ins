# -*- coding: utf-8 -*-
"""Build a self-contained pseudo-realtime GNSS/INS trajectory player.

Data source rule:
Only the C++ pseudo-online NAV output is used for the realtime display.  The
script deliberately does not read the commercial/reference trajectory files or
replace pseudo-online with forward post-processing output.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "成果_伪实时处理_展示图"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def find_pseudo_online_nav() -> Path:
    preferred = [
        ROOT / "output" / "伪实时_online测试_fast" / "gnss_ins_pseudo_online.nav",
        ROOT / "output" / "伪实时_online测试_after_remove" / "gnss_ins_pseudo_online.nav",
        ROOT / "output" / "gnss_ins_pseudo_online.nav",
    ]
    for path in preferred:
        if path.exists():
            return path
    candidates = list((ROOT / "output").rglob("gnss_ins_pseudo_online.nav"))
    if not candidates:
        raise FileNotFoundError(
            "没有找到 C++ 伪实时输出 gnss_ins_pseudo_online.nav。请先运行："
            "GNSSINS.exe <输出目录> 3000 pseudo-online-fast"
        )
    return max(candidates, key=lambda p: p.stat().st_size)


def read_nav(path: Path) -> dict[str, np.ndarray]:
    rows: list[list[float]] = []
    modes: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 15:
                continue
            rows.append([float(p[i]) for i in range(14)])
            modes.append(p[14])
    a = np.asarray(rows, dtype=float)
    keys = ["tow", "x", "y", "z", "lat", "lon", "h", "vn", "ve", "vd", "roll", "pitch", "yaw", "q"]
    nav = {k: a[:, i] for i, k in enumerate(keys)}
    nav["mode"] = np.asarray(modes, dtype=object)
    nav["speed"] = np.sqrt(nav["vn"] ** 2 + nav["ve"] ** 2 + nav["vd"] ** 2)
    return nav


def ecef_to_local(nav: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    lat = math.radians(float(nav["lat"][0]))
    lon = math.radians(float(nav["lon"][0]))
    slat, clat = math.sin(lat), math.cos(lat)
    slon, clon = math.sin(lon), math.cos(lon)
    dx = nav["x"] - float(nav["x"][0])
    dy = nav["y"] - float(nav["y"][0])
    dz = nav["z"] - float(nav["z"][0])
    out = dict(nav)
    out["north"] = -slat * clon * dx - slat * slon * dy + clat * dz
    out["east"] = -slon * dx + clon * dy
    out["up"] = clat * clon * dx + clat * slon * dy + slat * dz
    return out


def pick_rate(nav: dict[str, np.ndarray], hz: float = 5.0) -> np.ndarray:
    t = nav["tow"]
    keep = [0]
    next_t = t[0] + 1.0 / hz
    for i in range(1, len(t)):
        if t[i] + 1e-9 >= next_t:
            keep.append(i)
            next_t += 1.0 / hz
    return np.asarray(keep, dtype=int)


def decimate_for_html(idx: np.ndarray, max_points: int = 9000) -> np.ndarray:
    if len(idx) <= max_points:
        return idx
    return idx[np.linspace(0, len(idx) - 1, max_points, dtype=int)]


def write_stream(nav: dict[str, np.ndarray], idx: np.ndarray) -> Path:
    path = OUT / "cpp_pseudo_online_stream_5hz.csv"
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "seq", "tow", "elapsed_s", "east_m", "north_m", "height_m",
            "vn_mps", "ve_mps", "vd_mps", "speed_mps",
            "roll_deg", "pitch_deg", "yaw_deg", "quality", "mode",
        ])
        t0 = float(nav["tow"][idx[0]])
        for seq, i in enumerate(idx):
            w.writerow([
                seq,
                f"{nav['tow'][i]:.3f}",
                f"{nav['tow'][i] - t0:.3f}",
                f"{nav['east'][i]:.4f}",
                f"{nav['north'][i]:.4f}",
                f"{nav['h'][i]:.4f}",
                f"{nav['vn'][i]:.4f}",
                f"{nav['ve'][i]:.4f}",
                f"{nav['vd'][i]:.4f}",
                f"{nav['speed'][i]:.4f}",
                f"{nav['roll'][i]:.4f}",
                f"{nav['pitch'][i]:.4f}",
                f"{nav['yaw'][i]:.4f}",
                int(nav["q"][i]),
                str(nav["mode"][i]),
            ])
    return path


def make_snapshot(nav: dict[str, np.ndarray], idx: np.ndarray) -> Path:
    marks = [0.25, 0.50, 0.75, 1.00]
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2))
    e = nav["east"][idx]
    n = nav["north"][idx]
    for ax, frac in zip(axes.ravel(), marks):
        k = max(2, int(len(idx) * frac))
        sel = idx[:k]
        ax.plot(e, n, color="#d1d5db", lw=0.8, label="自编结果完整范围")
        ax.plot(nav["east"][sel], nav["north"][sel], color="#2563eb", lw=1.5, label="已接收历元")
        ax.scatter(nav["east"][sel[-1]], nav["north"][sel[-1]], s=36, color="#dc2626", zorder=3, label="当前历元")
        ax.set_title(f"接收进度 {int(frac * 100)}%")
        ax.set_xlabel("East / m")
        ax.set_ylabel("North / m")
        ax.axis("equal")
        ax.grid(True, color="#cbd5e1", alpha=0.55)
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9)
    fig.suptitle("自编 GNSS/INS 伪实时轨迹按历元推进", fontsize=15)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    path = OUT / "伪实时轨迹推进快照_不含参考轨迹.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_quality_plot(nav: dict[str, np.ndarray], idx: np.ndarray) -> Path:
    t = nav["tow"][idx] - nav["tow"][idx[0]]
    fig, axes = plt.subplots(3, 1, figsize=(10.5, 7.2), sharex=True)
    axes[0].plot(t, nav["speed"][idx], color="#2563eb", lw=1.0)
    axes[0].set_ylabel("speed / m/s")
    axes[0].set_title("伪实时速度输出")
    axes[1].plot(t, nav["h"][idx], color="#16a34a", lw=1.0)
    axes[1].set_ylabel("height / m")
    axes[1].set_title("伪实时高程输出")
    axes[2].step(t, nav["q"][idx], where="post", color="#dc2626", lw=1.0)
    axes[2].set_ylabel("quality")
    axes[2].set_xlabel("elapsed time / s")
    axes[2].set_title("GNSS/INS 逐历元质量标记")
    for ax in axes:
        ax.grid(True, color="#cbd5e1", alpha=0.55)
    fig.tight_layout()
    path = OUT / "伪实时速度高程质量时间序列.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_html(nav: dict[str, np.ndarray], idx: np.ndarray) -> Path:
    points = []
    t0 = float(nav["tow"][idx[0]])
    for i in idx:
        points.append({
            "t": round(float(nav["tow"][i] - t0), 3),
            "tow": round(float(nav["tow"][i]), 3),
            "e": round(float(nav["east"][i]), 4),
            "n": round(float(nav["north"][i]), 4),
            "h": round(float(nav["h"][i]), 4),
            "speed": round(float(nav["speed"][i]), 4),
            "roll": round(float(nav["roll"][i]), 4),
            "pitch": round(float(nav["pitch"][i]), 4),
            "yaw": round(float(nav["yaw"][i]), 4),
            "q": int(nav["q"][i]),
            "mode": str(nav["mode"][i]),
        })
    js_points = json.dumps(points, ensure_ascii=False)
    html = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>GNSS/INS C++ 伪实时轨迹播放器</title>
  <style>
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; font-family: "Microsoft YaHei", "Segoe UI", Arial, sans-serif; color: #111827; background: #f8fafc; }}
    .app {{ display: grid; grid-template-columns: 1fr 340px; min-height: 100vh; }}
    main {{ padding: 18px; }}
    aside {{ border-left: 1px solid #dbe3ef; background: #ffffff; padding: 18px; }}
    h1 {{ margin: 0 0 10px; font-size: 22px; }}
    .note {{ margin: 0 0 14px; color: #475569; font-size: 13px; }}
    canvas {{ width: 100%; height: calc(100vh - 132px); min-height: 520px; background: white; border: 1px solid #cbd5e1; border-radius: 6px; }}
    .controls {{ display: flex; gap: 8px; align-items: center; margin-top: 10px; }}
    button {{ border: 1px solid #94a3b8; background: #fff; border-radius: 6px; padding: 8px 12px; cursor: pointer; }}
    button.primary {{ background: #2563eb; border-color: #2563eb; color: white; }}
    input[type="range"] {{ flex: 1; }}
    .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 10px; }}
    .metric {{ border: 1px solid #e2e8f0; border-radius: 6px; padding: 10px; background: #f8fafc; }}
    .metric span {{ display: block; color: #64748b; font-size: 12px; }}
    .metric strong {{ display: block; margin-top: 4px; font-size: 18px; }}
    .status {{ margin: 12px 0; padding: 10px; border-radius: 6px; background: #eff6ff; color: #1e3a8a; font-size: 13px; }}
    .bar {{ height: 9px; background: #e2e8f0; border-radius: 999px; overflow: hidden; margin: 12px 0; }}
    .bar > div {{ height: 100%; width: 0%; background: #2563eb; }}
    .small {{ color: #64748b; font-size: 12px; line-height: 1.6; }}
    @media (max-width: 900px) {{ .app {{ grid-template-columns: 1fr; }} aside {{ border-left: 0; border-top: 1px solid #dbe3ef; }} canvas {{ height: 58vh; }} }}
  </style>
</head>
<body>
<div class="app">
  <main>
    <h1>GNSS/INS C++ 伪实时轨迹播放器</h1>
    <p class="note">数据源：自编 C++ pseudo-online 模式输出的 gnss_ins_pseudo_online.nav。页面只播放伪实时导航结果，不使用商业参考轨迹。</p>
    <canvas id="scene"></canvas>
    <div class="controls">
      <button id="play" class="primary">播放</button>
      <button id="reset">重置</button>
      <label>速度
        <select id="speed">
          <option value="1">1x</option>
          <option value="5" selected>5x</option>
          <option value="10">10x</option>
          <option value="25">25x</option>
          <option value="50">50x</option>
        </select>
      </label>
      <input id="scrub" type="range" min="0" max="{len(points)-1}" value="0" />
    </div>
  </main>
  <aside>
    <h1>实时状态</h1>
    <div class="bar"><div id="progress"></div></div>
    <div class="grid">
      <div class="metric"><span>历元序号</span><strong id="seq">0</strong></div>
      <div class="metric"><span>TOW / s</span><strong id="tow">0</strong></div>
      <div class="metric"><span>East / m</span><strong id="east">0</strong></div>
      <div class="metric"><span>North / m</span><strong id="north">0</strong></div>
      <div class="metric"><span>Height / m</span><strong id="height">0</strong></div>
      <div class="metric"><span>Speed / m/s</span><strong id="spd">0</strong></div>
      <div class="metric"><span>Roll / deg</span><strong id="roll">0</strong></div>
      <div class="metric"><span>Pitch / deg</span><strong id="pitch">0</strong></div>
      <div class="metric"><span>Yaw / deg</span><strong id="yaw">0</strong></div>
      <div class="metric"><span>Quality</span><strong id="quality">0</strong></div>
    </div>
    <div class="status" id="mode">mode</div>
    <p class="small">现场展示口径：该页面展示 C++ 伪实时程序逐历元输出后的实时显示过程，只使用伪实时输出点。拖动进度条可模拟任意时刻已接收数据。</p>
  </aside>
</div>
<script>
const pts = {js_points};
const canvas = document.getElementById('scene');
const ctx = canvas.getContext('2d');
const playBtn = document.getElementById('play');
const resetBtn = document.getElementById('reset');
const scrub = document.getElementById('scrub');
const speedSel = document.getElementById('speed');
let i = 0, playing = false, last = performance.now(), acc = 0;
const minE = Math.min(...pts.map(p => p.e)), maxE = Math.max(...pts.map(p => p.e));
const minN = Math.min(...pts.map(p => p.n)), maxN = Math.max(...pts.map(p => p.n));
function resize() {{
  const r = canvas.getBoundingClientRect();
  canvas.width = Math.max(600, Math.floor(r.width * devicePixelRatio));
  canvas.height = Math.max(420, Math.floor(r.height * devicePixelRatio));
  draw();
}}
function xy(p) {{
  const pad = 46 * devicePixelRatio;
  const w = canvas.width, h = canvas.height;
  const span = Math.max(maxE - minE, maxN - minN) || 1;
  const cx = (minE + maxE) / 2, cy = (minN + maxN) / 2;
  const x0 = cx - span * 0.58, y0 = cy - span * 0.58;
  const x = pad + (p.e - x0) / (span * 1.16) * (w - 2 * pad);
  const y = h - pad - (p.n - y0) / (span * 1.16) * (h - 2 * pad);
  return [x, y];
}}
function line(points, color, width) {{
  if (points.length < 2) return;
  ctx.beginPath();
  points.forEach((p, k) => {{
    const [x, y] = xy(p);
    if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }});
  ctx.strokeStyle = color;
  ctx.lineWidth = width * devicePixelRatio;
  ctx.stroke();
}}
function draw() {{
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.strokeStyle = '#e2e8f0';
  ctx.lineWidth = devicePixelRatio;
  for (let k = 0; k <= 10; k++) {{
    const x = canvas.width * k / 10, y = canvas.height * k / 10;
    ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, canvas.height); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(canvas.width, y); ctx.stroke();
  }}
  line(pts, '#d1d5db', 1.2);
  line(pts.slice(0, i + 1), '#2563eb', 2.2);
  const p = pts[i], [x, y] = xy(p);
  ctx.beginPath(); ctx.arc(x, y, 6 * devicePixelRatio, 0, Math.PI * 2); ctx.fillStyle = '#dc2626'; ctx.fill();
  ctx.fillStyle = '#111827'; ctx.font = `${{14 * devicePixelRatio}}px Microsoft YaHei`;
  ctx.fillText('完整范围', 24 * devicePixelRatio, 28 * devicePixelRatio);
  ctx.fillStyle = '#2563eb'; ctx.fillText('已接收轨迹', 24 * devicePixelRatio, 50 * devicePixelRatio);
  updatePanel();
}}
function updatePanel() {{
  const p = pts[i], pct = i / (pts.length - 1) * 100;
  document.getElementById('seq').textContent = `${{i + 1}} / ${{pts.length}}`;
  document.getElementById('tow').textContent = p.tow.toFixed(1);
  document.getElementById('east').textContent = p.e.toFixed(2);
  document.getElementById('north').textContent = p.n.toFixed(2);
  document.getElementById('height').textContent = p.h.toFixed(2);
  document.getElementById('spd').textContent = p.speed.toFixed(2);
  document.getElementById('roll').textContent = p.roll.toFixed(2);
  document.getElementById('pitch').textContent = p.pitch.toFixed(2);
  document.getElementById('yaw').textContent = p.yaw.toFixed(2);
  document.getElementById('quality').textContent = p.q;
  document.getElementById('mode').textContent = p.mode;
  document.getElementById('progress').style.width = `${{pct}}%`;
  scrub.value = i;
}}
function tick(now) {{
  const dt = (now - last) / 1000; last = now;
  if (playing) {{
    acc += dt * Number(speedSel.value) * 5;
    while (acc >= 1 && i < pts.length - 1) {{ i++; acc--; }}
    if (i >= pts.length - 1) {{ playing = false; playBtn.textContent = '播放'; }}
    draw();
  }}
  requestAnimationFrame(tick);
}}
playBtn.onclick = () => {{ playing = !playing; playBtn.textContent = playing ? '暂停' : '播放'; }};
resetBtn.onclick = () => {{ i = 0; playing = false; playBtn.textContent = '播放'; draw(); }};
scrub.oninput = () => {{ i = Number(scrub.value); draw(); }};
window.addEventListener('resize', resize);
resize(); requestAnimationFrame(tick);
</script>
</body>
</html>
"""
    path = OUT / "pseudo_realtime_trajectory_player.html"
    path.write_text(html, encoding="utf-8")
    return path


def write_summary(stream: Path, html: Path, snapshot: Path, quality: Path, n: int, duration: float) -> Path:
    path = OUT / "C++伪实时轨迹播放器说明.md"
    path.write_text(
        "\n".join([
            "# 伪实时轨迹播放器成果说明",
            "",
            "本成果只使用自编 C++ 程序 `pseudo-online` 模式输出的 `gnss_ins_pseudo_online.nav`，不读取、不绘制商业参考轨迹，也不使用前向后处理结果替代。",
            "",
            "## 输出文件",
            "",
            f"- `{html.name}`：可交互伪实时轨迹播放器，可播放、暂停、调速、拖动进度。",
            f"- `{stream.name}`：由 C++ 伪实时输出抽取的 5 Hz 显示数据流。",
            f"- `{snapshot.name}`：25%、50%、75%、100%四个接收进度的轨迹快照。",
            f"- `{quality.name}`：伪实时速度、高程和质量标记时间序列。",
            "",
            "## 数据规模",
            "",
            f"- 播放点数：{n}",
            f"- 覆盖时长：{duration:.1f} s",
            "",
            "## 现场讲解",
            "",
            "该播放器模拟 C++ 伪实时导航状态按历元逐点到达后的实时显示。页面中的灰色线只是 C++ 伪实时结果的完整范围预览，蓝色线表示当前已经接收并显示的历元，红点表示当前最新历元。",
        ]),
        encoding="utf-8",
    )
    return path


def main() -> None:
    src_nav = find_pseudo_online_nav()
    nav = ecef_to_local(read_nav(src_nav))
    idx_5hz = pick_rate(nav, 5.0)
    stream = write_stream(nav, idx_5hz)
    idx_html = decimate_for_html(idx_5hz)
    snapshot = make_snapshot(nav, idx_html)
    quality = make_quality_plot(nav, idx_html)
    html = make_html(nav, idx_html)
    summary = write_summary(stream, html, snapshot, quality, len(idx_html), float(nav["tow"][idx_html[-1]] - nav["tow"][idx_html[0]]))
    print(f"source: {src_nav}")
    print(f"stream: {stream}")
    print(f"html: {html}")
    print(f"summary: {summary}")


if __name__ == "__main__":
    main()
