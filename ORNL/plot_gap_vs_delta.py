"""Gap knob: physical gap E1 - E0 vs dimerisation delta, for several N (exact diagonalisation).

    python plot_gap_vs_delta.py                          # N = 8 10 12 (about a minute), J2 from the config (0.2411)
    python plot_gap_vs_delta.py --N 8 12 16 --set J2=0   # plain Heisenberg
    python plot_gap_vs_delta.py --points 81 --out figs/gap_vs_delta.pdf

H = sum_i J1 (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}. delta > 0 strengthens the end bonds
(dimers, gap opens to ~1.8 J1 at delta -> 1); delta < 0 weakens them (nearly free edge spins, gap -> 0).
Dots mark the dimerisations that hit the study's target gaps (cfg.gaps).
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np

from ornl import Chain, StudyConfig
from ornl.config import add_cli


def main():
    ap = add_cli(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter))
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12])
    ap.add_argument("--points", type=int, default=61, help="delta values across cfg.delta_range")
    ap.add_argument("--out", default="figs/gap_vs_delta.png")
    a = ap.parse_args()
    cfg = StudyConfig.from_cli(a.config, a.set)
    deltas = np.linspace(*cfg.delta_range, a.points)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    for N in a.N:
        gaps = [Chain(N, d, cfg.J1, cfg.J2).gap for d in deltas]
        line, = ax1.plot(deltas, gaps, label=f"N = {N}")
        ax2.semilogy(deltas, gaps, color=line.get_color())
        k = int(np.argmin(gaps))                                  # below this delta the gap is ~0 with level crossings
        g_up, d_up = np.array(gaps[k:]), deltas[k:]
        tgt = [g for g in cfg.gaps if g_up.min() <= g <= g_up.max()]
        ds = np.interp(np.log(tgt), np.log(np.maximum.accumulate(g_up)), d_up)
        ax2.plot(ds, tgt, "o", ms=3, color=line.get_color())
        print(f"N={N}: gap {min(gaps):.4f} .. {max(gaps):.4f}; at delta=0: {Chain(N, 0.0, cfg.J1, cfg.J2).gap:.4f}")
    for ax in (ax1, ax2):
        ax.set_xlabel("dimerisation δ")
        ax.set_ylabel("gap E1 − E0  (J1)")
        ax.axvline(0, color="0.7", lw=0.8)
    ax1.legend()
    ax1.set_title(f"J2 = {cfg.J2}")
    ax2.set_title("log scale; dots = study's target gaps")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print("saved", a.out)


if __name__ == "__main__":
    main()
