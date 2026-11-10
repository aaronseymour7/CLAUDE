"""Cost with and without the SBC S^2 stage vs the gap, on the original grid (triplet gap 0.005-1.7 J1).

    python plot_s2_crossover.py      # reads results/map_N*.json and results/su2_D_triplet_grid_s2stage/
                                     # -> figs/s2_crossover.png, results/s2_crossover.csv

Same Hamiltonians and same staggered-Rz trials in both; "with" adds the S^2 stage (its CNOTs and success
probability included) and then filters at the singlet gap. Lower panels: the share of the S^2 stage in the CNOTs
of one attempt, i.e. the price of resolving the spin contamination.
"""
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GAMMAS = (0.5, 0.1, 0.01)
EPS, PR = "0.01", "rus"
COL = {8: "#7fb3d5", 10: "#2471a3", 12: "#154360"}


def load(N):
    a = json.load(open(f"results/map_N{N}.json"))
    d = json.load(open(f"results/su2_D_triplet_grid_s2stage/map_N{N}.json"))
    return a, d


rows = []
fig, axs = plt.subplots(2, len(GAMMAS), figsize=(4.2 * len(GAMMAS), 6.4), sharex=True,
                        gridspec_kw=dict(height_ratios=(2.2, 1)))
for N in (8, 10, 12):
    if not os.path.exists(f"results/su2_D_triplet_grid_s2stage/map_N{N}.json"):
        continue
    a, d = load(N)
    for col, g in enumerate(GAMMAS):
        j = a["gammas"].index(g)
        x, ya, yd, share = [], [], [], []
        for k in sorted(a["rows"], key=float):
            ra, rd = a["rows"][k], d["rows"].get(k)
            if rd is None:
                continue
            ba, bd = ra["best"][j][EPS][PR], rd["best"][j][EPS][PR]
            if not (ba and bd):
                continue
            x.append(ra["gap"] / ra["W"]); ya.append(ba["cost"]); yd.append(bd["cost"])
            share.append(bd["s2_cnot"] / bd["cnot"])
            rows.append(dict(N=N, gamma=g, gap_triplet=ra["gap"], gap_singlet=rd["gap_singlet"], W=ra["W"],
                             cost_without=ba["cost"], cost_with=bd["cost"], ratio=bd["cost"] / ba["cost"],
                             s2_share_of_attempt=share[-1], s2_steps=bd["s2_n"], p_without=ba["p"], p_with=bd["p"]))
        axs[0, col].loglog(x, ya, "o-", ms=3.5, color=COL[N], label=f"N={N} without S²")
        axs[0, col].loglog(x, yd, "s--", ms=3.5, color=COL[N], mfc="none", label=f"N={N} with S²")
        axs[1, col].semilogx(x, share, "s--", ms=3.5, color=COL[N], mfc="none")
for col, g in enumerate(GAMMAS):
    axs[0, col].set_title(f"γ = {g}, ε = 1e-2 (RUS)")
    axs[0, col].grid(alpha=0.3, which="both")
    axs[1, col].grid(alpha=0.3, which="both")
    axs[1, col].set_ylim(0, 1)
    axs[1, col].set_xlabel("triplet gap / W  (the gap without the S² stage)")
axs[0, 0].set_ylabel("CNOT per success")
axs[1, 0].set_ylabel("S² stage share\nof one attempt")
axs[0, 0].legend(fontsize=7, ncol=1)
fig.tight_layout()
os.makedirs("figs", exist_ok=True)
fig.savefig("figs/s2_crossover.png", dpi=150)
with open("results/s2_crossover.csv", "w", newline="") as f:
    w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
# crossover: smallest triplet gap/W where with >= without, per (N, gamma)
for N in (8, 10, 12):
    for g in GAMMAS:
        s = sorted((r["gap_triplet"] / r["W"], r["ratio"]) for r in rows if r["N"] == N and r["gamma"] == g)
        cross = next((x for x, q in s if q >= 1), None)
        print(f"N={N} gamma={g}: with/without from {s[0][1]:.3f} (gap/W {s[0][0]:.4f}) to {s[-1][1]:.2f}; "
              f"first gap/W where S^2 stops paying: {cross}")
