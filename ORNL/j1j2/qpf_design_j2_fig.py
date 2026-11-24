"""CNOT per success vs J2 at fixed N for the three designs (from qpf_design_eps_scan.py --J2s ... --trial mps).
Usage: python j1j2/qpf_design_j2_fig.py RESULTS_DIR N OUT_PNG"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, N, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
runs = sorted((json.load(open(f)) for f in glob.glob(os.path.join(d, f"eps_scan_mps_N{N}_J2*.json"))),
              key=lambda r: r["J2"])
COL = {"v4": "#2a78d6", "v5 same T no cap": "#eb6834", "v5 1.5T no cap": "#1baf7a"}
J2 = np.array([r["J2"] for r in runs])
fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
for eps, ls in ((1e-2, "-"), (1e-3, "--")):
    for name, col in COL.items():
        y = [next((x.get("cx_per_success") for x in r["rows"] if x["design"] == name and x["eps"] == eps), None)
             for r in runs]
        y = np.array([np.nan if v is None else v for v in y], float)
        axs[0].semilogy(J2, y, ls, marker="o", color=col, lw=1.8, label=f"{name}, ε = {eps:g}")
        R = [next((x.get("R") for x in r["rows"] if x["design"] == name and x["eps"] == eps), None) for r in runs]
        if eps == 1e-2:
            for xx, yy, rr in zip(J2, y, R):
                if rr and np.isfinite(yy):
                    axs[0].annotate(f"R{rr}", (xx, yy), textcoords="offset points", xytext=(4, 4), fontsize=7, color=col)
axs[0].set_ylabel("CNOT per success (RUS)"); axs[0].set_title(f"N = {N}: cost per success (labels: blocks R at ε = 1e-2)",
                                                               fontsize=10)
axs[0].legend(fontsize=7, frameon=False)
axs[1].plot(J2, [r["gamma"] for r in runs], "o-", color="#52514e")
axs[1].set_ylabel("trial overlap γ = |⟨E0|trial⟩|²"); axs[1].set_title("trial quality", fontsize=10)
axs[2].plot(J2, [r["gap_scaled"] for r in runs], "o-", color="#52514e")
axs[2].set_ylabel("scaled gap Δ/W (global)"); axs[2].set_title("gap the filter is designed for", fontsize=10)
for ax in axs:
    ax.set_xlabel("J2"); ax.grid(alpha=0.25, which="both")
fig.suptitle("Uniform open J1-J2 chain, DMRG MPS → mps_to_circuit (L = 1) trial, energy filter only, 2nd-order Trotter, "
             "analytic 6-CNOT model (trial circuit CNOTs included)", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(out, dpi=150)
