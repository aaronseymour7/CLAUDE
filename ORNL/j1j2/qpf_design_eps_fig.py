"""Plots for qpf_design_eps_scan.py: CNOT per block, Trotter steps and CNOT per success vs eps, per design and N.
Usage: python j1j2/qpf_design_eps_fig.py RESULTS_DIR OUT_PNG"""
import glob
import json
import os
import sys

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
runs = sorted((json.load(open(f)) for f in glob.glob(os.path.join(d, "eps_scan_N*.json"))), key=lambda r: r["N"])
COL = {"v4": "#2a78d6", "v5 same T no cap": "#eb6834", "v5 1.5T no cap": "#1baf7a"}
MK = {8: "o", 10: "s", 12: "^", 14: "D"}
LS = {8: ":", 10: "--", 12: "-", 14: "-."}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
for r in runs:
    N = r["N"]
    for name, col in COL.items():
        rows = [x for x in r["rows"] if x["design"] == name and x.get("cx_block") is not None]
        e = np.array([x["eps"] for x in rows])
        for ax, key in zip(axs, ("cx_block", "steps", "cx_per_success")):
            ax.loglog(e, [x[key] for x in rows], LS[N], marker=MK[N], ms=5, color=col, lw=1.6)
for ax, yl in zip(axs, ("CNOT per block (one pass of all pulses)", "Trotter steps, total (R blocks)",
                        "CNOT per success (RUS)")):
    ax.axvline(1e-2, color="#52514e", lw=0.8, ls="--")
    ax.set_xlabel("target error ε = 1 − F"); ax.set_ylabel(yl); ax.grid(alpha=0.25, which="both")
    ax.invert_xaxis()
h = [plt.Line2D([], [], color=c, lw=2, label=n) for n, c in COL.items()]
h += [plt.Line2D([], [], color="#52514e", ls=LS[r["N"]], marker=MK[r["N"]], label=f"N = {r['N']} "
      f"(γ = {r['gamma']:.2f}, Δ/W = {r['gap_scaled']:.3f}, {r['cx_per_step']} CNOT/step)") for r in runs]
fig.legend(handles=h, loc="lower center", ncol=3, fontsize=8, frameon=False)
fig.suptitle(f"Uniform J1-J2 chain, J2 = {runs[0]['J2']}, Néel trial, energy filter only, 2nd-order Trotter, "
             "analytic 6-CNOT model", fontsize=10)
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(out, dpi=150)
