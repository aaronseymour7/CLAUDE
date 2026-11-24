"""Overlay f(E)^2 for the cost functions in qpf_logcost_filters.py, with the ideal step.
Usage: python j1j2/qpf_logcost_fig.py RESULTS_DIR OUT_PREFIX"""
import os
import sys

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, pre = sys.argv[1], sys.argv[2]
D = np.load(os.path.join(d, "logcost_designs.npy"), allow_pickle=True).item()
SHOW = [("v4 (authors)", "#2a78d6", "-", 2.2), ("log v4", "#9ec5f4", "--", 1.6),
        ("log-mean, floor 0.0001", "#eda100", "-", 1.6), ("log-mean, floor 1e-12", "#eb6834", "-", 1.6),
        ("minimax (log max)", "#1baf7a", "-", 2.2)]
for c in (1.0, 1.5):
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for j, g in enumerate((0.1, 0.3)):
        E = np.linspace(0, 1, 8001)
        for name, col, ls, lw in SHOW:
            t = D[f"{g}|{c}|{name}"]
            f2 = np.prod(np.cos(np.outer(E, t)) ** 2, axis=1)
            mx = f2[E >= g].max()
            axs[0, j].plot(E, f2, ls, color=col, lw=lw, label=f"{name}: max above Δ {mx:.0e}")
            axs[1, j].semilogy(E, np.maximum(f2, 1e-12), ls, color=col, lw=lw * 0.7)
        axs[0, j].plot([0, g, g, 1], [1, 1, 0, 0], color="#0b0b0b", lw=1.2, ls=":", label="ideal step")
        for ax in axs[:, j]:
            ax.axvline(g, color="#52514e", lw=0.8, ls="--"); ax.grid(alpha=0.25)
        axs[0, j].set_title(f"scaled gap Δ = {g}, T = {c:g}·π/Δ", fontsize=11)
        axs[0, j].legend(fontsize=7.5, frameon=False, loc="upper right")
        axs[1, j].set_ylim(1e-8, 2); axs[1, j].set_xlabel("scaled energy (E − E0)/W")
    axs[0, 0].set_ylabel("kept weight f(E)²  (linear)")
    axs[1, 0].set_ylabel("kept weight f(E)²  (log)")
    fig.tight_layout(); fig.savefig(f"{pre}_c{c:g}.png", dpi=150)
