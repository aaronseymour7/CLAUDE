"""Filter amplitude f(E) = prod cos(t_i E) on the scaled band (the 'step'): authors' v4 vs the minimax option."""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from j1j2.route_costs import pick
from qpf import design_plus as dp

out = sys.argv[1]
C = ["#2a78d6", "#eb6834", "#1baf7a"]
gaps = (0.1, 0.3)
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
for j, g in enumerate(gaps):
    E = np.linspace(0, 1, 8001)
    sets = [("authors v4, T = π/Δ", pick(g, "v4")["times"]),
            ("minimax, T = π/Δ", dp.design(g, "minimax", 1.0)["times"]),
            ("minimax, T = 1.5π/Δ", dp.design(g, "minimax", 1.5)["times"])]
    for col, (lab, t) in zip(C, sets):
        f = np.prod(np.cos(np.outer(E, t)), axis=1)
        lab2 = f"{lab} (m={len(t)}, max f² above Δ = {dp.max_leak(g, t):.0e})"
        axs[0, j].plot(E, f ** 2, color=col, lw=1.8, label=lab2)
        axs[1, j].semilogy(E, np.maximum(f ** 2, 1e-12), color=col, lw=1.4)
    axs[0, j].plot([0, g, g, 1], [1, 1, 0, 0], color="#52514e", lw=1, ls=":", label="ideal step")
    for ax in axs[:, j]:
        ax.axvline(g, color="#52514e", lw=0.8, ls="--"); ax.grid(alpha=0.25)
    axs[0, j].set_title(f"scaled gap Δ = {g}", fontsize=11)
    axs[0, j].legend(fontsize=8, frameon=False, loc="upper right")
    axs[1, j].set_ylim(1e-8, 2); axs[1, j].set_xlabel("scaled energy (E − E0)/W")
axs[0, 0].set_ylabel("kept weight f(E)²  (linear)")
axs[1, 0].set_ylabel("kept weight f(E)²  (log)")
fig.tight_layout(); fig.savefig(out, dpi=150)
