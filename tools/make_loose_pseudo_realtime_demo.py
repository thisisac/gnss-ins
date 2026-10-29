# -*- coding: utf-8 -*-
"""
生成松组合伪实时展示用数据。

说明：
    松组合正式解算仍由 C++ 主程序输出。本脚本不改变解算结果，只把
    gnss_ins_smooth.nav 按历元顺序整理成伪实时回放目录，用于展示松组合
    结果按历元推进的轨迹播放效果。
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "output" / "成果_REAL_松组合真实算法" / "gnss_ins_smooth.nav"
OUT_DIR = ROOT / "output" / "loose_pseudo_realtime_demo"
OUT_NAV = OUT_DIR / "gnss_ins_loose_pseudo_realtime_player.nav"
OUT_STREAM = OUT_DIR / "loose_realtime_stream.txt"


def main() -> None:
    if not SRC.exists():
        raise FileNotFoundError(f"找不到松组合平滑结果：{SRC}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    count = 0
    with SRC.open("r", encoding="utf-8", errors="ignore") as src, \
            OUT_NAV.open("w", encoding="utf-8", newline="\n") as nav, \
            OUT_STREAM.open("w", encoding="utf-8", newline="\n") as stream:
        for line in src:
            if line.startswith("#"):
                nav.write(line)
                stream.write(line)
                continue
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 15:
                parts[14] = parts[14] + "_LOOSE_PSEUDO_REALTIME"
                line = " ".join(parts) + "\n"
                count += 1
            nav.write(line)
            stream.write(line)

    print(f"loose pseudo realtime epochs: {count}")
    print(f"NAV: {OUT_NAV}")
    print(f"STREAM: {OUT_STREAM}")


if __name__ == "__main__":
    main()
