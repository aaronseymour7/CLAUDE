"""Gap x overlap CNOT maps at eps = 1e-2 and fixed N, dyadic SBC filter, measured by simulation.

    python sbpipe_study/gapov_report/make_gridmap.py
Reads results/gapov.jsonl and results/gapov_grid.jsonl. Filter-only CNOTs (trial circuit excluded).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_figures as M  # noqa: E402
from make_heat import _lab  # noqa: E402

plt = M.plt
FIGD = os.path.join(HERE, "figs")
GAPS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
GAMS = [0.03, 0.1, 0.3, 0.5, 0.7, 0.9]
NS = (12, 14)
EPS = 1e-2


def load():
    from derive_gapov_filter import load_cases
    P = {}
    for r in load_cases(("gapov.jsonl", "gapov_grid.jsonl")):
        if r["eps"] == EPS and r["N"] in NS and r.get("meas_dyadic"):
            P[(r["N"], r["gap_target"], r["gamma_target"])] = r
    return P


def _filter_cnot(r, fam):
    m = r.get(f"meas_{fam}")
    if not m:
        return None, None
    fc = m["cnot_filter"] if "cnot_filter" in m else m["cnot"] * m["p"] - r["trial_cnot"]
    return fc, m["p"]


def one(P, kind, fname, title):
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.8), sharex=True, sharey=True)
    ye = np.arange(len(GAMS), dtype=float)          # equal-height rows, one per overlap value
    yedge = np.arange(len(GAMS) + 1) - 0.5
    xedge = np.concatenate([[GAPS[0] - 0.025], np.array(GAPS[:-1]) + 0.025, [GAPS[-1] + 0.025]])
    Zs = {}
    for fam in ("dyadic", "optimized"):
        for N in NS:
            Z = np.full((len(GAMS), len(GAPS)), np.nan)
            for (n, gap, gm), r in P.items():
                if n != N or gap not in GAPS or gm not in GAMS:
                    continue
                fc, p = _filter_cnot(r, fam)
                if fc:
                    Z[GAMS.index(gm), GAPS.index(gap)] = np.log10(fc / p if kind == "succ" else fc)
            Zs[(fam, N)] = Z
    allz = np.concatenate([z.ravel() for z in Zs.values()])
    lo, hi = np.nanmin(allz), np.nanmax(allz)
    for row, fam in enumerate(("dyadic", "optimized")):
        for col, N in enumerate(NS):
            a, Z = axs[row, col], Zs[(fam, N)]
            pc = a.pcolormesh(xedge, yedge, Z, cmap="viridis", vmin=lo, vmax=hi, shading="flat")
            if np.isfinite(Z).all():
                a.contour(GAPS, ye, Z, levels=np.arange(np.floor(lo * 4) / 4, hi + 0.01, 0.25), colors="white",
                          linewidths=0.6, alpha=0.7)
            for i, gm in enumerate(ye):
                for j, gap in enumerate(GAPS):
                    if np.isfinite(Z[i, j]):
                        a.text(gap, gm, _lab(10 ** Z[i, j]), ha="center", va="center", fontsize=5.6,
                               color="white" if Z[i, j] < (lo + hi) / 2 else "black")
                    else:
                        a.text(gap, gm, "n/a", ha="center", va="center", fontsize=5.5, color=M.INK2)
            a.set_title(f"{fam} SBC, $N={N}$", loc="left", fontsize=7.5)
            a.set_xticks(GAPS)
            a.set_yticks(ye); a.set_yticklabels([f"{g:g}" for g in GAMS])
            a.grid(False)
    cb = fig.colorbar(pc, ax=axs, orientation="vertical", fraction=0.03, pad=0.02)
    cb.set_label(r"$\log_{10}$ CNOTs", fontsize=7)
    for a in axs[1]:
        a.set_xlabel(r"gap $\Delta$ (held fixed across $N$ by tuning $J_2$)")
    for a in axs[:, 0]:
        a.set_ylabel(r"trial overlap $\gamma$")
    fig.suptitle(title, fontsize=8.5, x=0.01, ha="left")
    fig.savefig(os.path.join(FIGD, fname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGD, fname + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    P = load()
    one(P, "succ", "fig_gridmap",
        r"Filter CNOTs per success, measured by simulation, $\varepsilon=10^{-2}$ (contours every factor 1.8)")
    one(P, "shot", "fig_gridmap_shot",
        r"Filter CNOTs per shot, measured by simulation, $\varepsilon=10^{-2}$ (contours every factor 1.8)")
    print({N: int(sum(1 for k in P if k[0] == N)) for N in NS})


if __name__ == "__main__":
    main()
