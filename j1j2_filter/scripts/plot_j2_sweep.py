"""Plot the J2 sweep (results/data_sweep) for N = 4,6,8,10: certified vs empirical design cost, gap and trial overlap.

usage: python scripts/plot_j2_sweep.py   -> results/figures_sweep/*.png and results/tables/j2_sweep.csv"""
import csv
import glob
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figures_sweep"
OUT.mkdir(parents=True, exist_ok=True)
NS = (4, 6, 8, 10)
COL = dict(zip(NS, ("tab:blue", "tab:orange", "tab:green", "tab:red")))

rows = []
for f in glob.glob(str(ROOT / "results" / "data_sweep" / "dmrg_N*_eps0.01.json")):
    d = json.load(open(f))
    inp = d["inputs"]
    r = dict(N=d["N"], J2=d["J2"], gap_scaled=inp["gap_scaled"], gap_raw=inp["gap_raw"], gamma=inp["gamma"],
             alpha=inp["alpha"], W=inp["W_cert"], pauli_terms=inp["n_pauli_terms"], no_filter=not d["designs"])
    for w in ("guaranteed", "empirical"):
        g = d["designs"].get(w)
        if g:
            r[f"{w}_n"] = g["n"]
            r[f"{w}_T_phys"] = g["T_phys"]
            r[f"{w}_rzz_filter"] = g["counts_rzz_a2a"]["filter"]["rzz"]
            r[f"{w}_1mF"] = g["eps_actual"]
    rows.append(r)
rows.sort(key=lambda r: (r["N"], r["J2"]))
keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("N", "J2"), k))
with open(ROOT / "results" / "tables" / "j2_sweep.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, keys)
    w.writeheader()
    w.writerows(rows)


def series(N, key):
    pts = [(r["J2"], r[key]) for r in rows if r["N"] == N and key in r]
    return [p[0] for p in pts], [p[1] for p in pts]


def panel(ax, base, ylabel, log=True):
    for N in NS:
        ax.plot(*series(N, f"guaranteed_{base}"), "o-", color=COL[N], label=f"N={N} certified")
        ax.plot(*series(N, f"empirical_{base}"), "s--", color=COL[N], mfc="none", label=f"N={N} empirical")
    ax.set_xlabel("$J_2/J_1$")
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")
    ax.grid(alpha=0.3)


for base, ylabel, name in (("n", "Trotter steps $n$", "steps_vs_J2"),
                           ("rzz_filter", "filter RZZ (level 3)", "rzz_vs_J2"),
                           ("T_phys", "evolution time $T$ ($1/J_1$)", "T_vs_J2")):
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    panel(ax, base, ylabel, log=base != "T_phys")
    ax.axvline(0.5, color="gray", lw=0.8, ls=":")
    ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax.set_title(r"First-order filter, $\varepsilon=10^{-2}$ ($J_2=0.5$: trial exact, no filter)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
for N in NS:
    axs[0].plot(*series(N, "gap_raw"), "o-", color=COL[N], label=f"N={N}")
    axs[1].plot(*series(N, "gamma"), "o-", color=COL[N], label=f"N={N}")
axs[0].set_ylabel("gap $E_1-E_0$ ($J_1$)")
axs[1].set_ylabel(r"trial overlap $\gamma$")
for a in axs:
    a.set_xlabel("$J_2/J_1$")
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "gap_gamma_vs_J2.png", dpi=200)
print("wrote", OUT, len(rows), "points")
