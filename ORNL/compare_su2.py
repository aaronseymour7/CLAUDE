"""Three ways to handle total spin, on the same Hamiltonians (dimerisation set by the singlet gap):

    A  staggered-Rz knob, no S^2 stage       (the current study: the filter must resolve the triplet gap)
    B  staggered-Rz knob + SBC S^2 stage     (Stetcu et al.: project J = 0, then filter at the singlet gap)
    C  singlet (SU(2) brick) knob            (trial is already a singlet; filter at the singlet gap)

    python run_maps.py --config configs/su2_A_staggered.json           # and _B_, _C_ (see configs/)
    python compare_su2.py                                               # -> results/su2_compare.csv, figs/su2_*.png

Costs are CNOTs per success (filter + knob layer + S^2 stage; the exact trial base costs nothing).
"""
import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

VARIANTS = dict(A="su2_A_staggered", B="su2_B_staggered_s2stage", C="su2_C_singlet_knob")
LABEL = dict(A="A: Rz knob, no S² stage", B="B: Rz knob + S² stage", C="C: singlet knob")
COLOR = dict(A="#c0392b", B="#2471a3", C="#1e8449")


def load(res="results"):
    out = {}
    for v, d in VARIANTS.items():
        for f in glob.glob(os.path.join(res, d, "map_N*.json")):
            m = json.load(open(f))
            out[(v, m["N"])] = m
    return out


def points(maps):
    rows = []
    for (v, N), m in maps.items():
        for key, r in m["rows"].items():
            for j, gam in enumerate(m["gammas"]):
                if r["gamma"][j] is None:
                    continue
                for e in m["eps"]:
                    b = r["best"][j][f"{e:g}"]
                    for pr in ("rus", "aa"):
                        x = b[pr]
                        rows.append(dict(variant=v, N=N, gap_singlet=r["gap_singlet"], gap_triplet=r["gap_triplet"],
                                         W=r["W"], delta=r["delta"], gamma=gam, eps=e, protocol=pr,
                                         cost=None if x is None else x["cost"], cnot=None if x is None else x["cnot"],
                                         p=None if x is None else x["p"], K=None if x is None else x["K"],
                                         r=None if x is None else x["r"], n=None if x is None else x["n"],
                                         s2_n=None if x is None else x.get("s2_n"),
                                         s2_cnot=None if x is None else x.get("s2_cnot"),
                                         s2_p=None if x is None else x.get("s2_p")))
    return rows


def table(rows, path):
    keys = list(rows[0])
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)


def lookup(rows):
    return {(r["variant"], r["N"], round(r["delta"], 8), r["gamma"], r["eps"], r["protocol"]): r for r in rows}


def summary(rows, eps=1e-2, protocol="rus"):
    L = lookup(rows)
    print(f"cost per success, eps = {eps:g}, {protocol}; ratios to A at the same Hamiltonian and gamma")
    for N in sorted({r["N"] for r in rows}):
        for v in ("B", "C"):
            rat = []
            for (vv, NN, d, g, e, pr), r in L.items():
                if vv != v or NN != N or e != eps or pr != protocol or r["cost"] is None:
                    continue
                a = L.get(("A", N, d, g, e, pr))
                if a and a["cost"]:
                    rat.append(r["cost"] / a["cost"])
            if rat:
                q = np.quantile(rat, [0, 0.5, 1])
                print(f"  N={N:2d} {v}/A: min {q[0]:.3g}  median {q[1]:.3g}  max {q[2]:.3g}  ({len(rat)} points)")


def fig_cuts(rows, N, eps=1e-2, protocol="rus", gammas=(0.5, 0.1, 0.01), out="figs/su2_cuts.png"):
    fig, axs = plt.subplots(1, len(gammas), figsize=(4 * len(gammas), 3.4), sharey=True)
    for ax, g in zip(axs, gammas):
        for v in VARIANTS:
            pts = sorted((r["gap_singlet"] / r["W"], r["cost"]) for r in rows
                         if r["variant"] == v and r["N"] == N and r["gamma"] == g and r["eps"] == eps
                         and r["protocol"] == protocol and r["cost"] is not None)
            if pts:
                x, y = zip(*pts)
                ax.loglog(x, y, "o-", ms=4, color=COLOR[v], label=LABEL[v])
        ax.set_title(f"N = {N}, γ = {g}")
        ax.set_xlabel("singlet gap / W")
        ax.grid(alpha=0.3, which="both")
    axs[0].set_ylabel(f"CNOT per success ({protocol}), ε = {eps:g}")
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fit_exponents(rows, eps=1e-2, protocol="rus", gamma_max=0.95):
    """ln C = c + a ln(W / gap_singlet) + b ln(1 / gamma), per variant and N."""
    out = []
    for v in VARIANTS:
        for N in sorted({r["N"] for r in rows}):
            s = [r for r in rows if r["variant"] == v and r["N"] == N and r["eps"] == eps and r["protocol"] == protocol
                 and r["cost"] is not None and r["gamma"] <= gamma_max and 1 - r["gamma"] > 10 * eps]
            if len(s) < 5:
                continue
            X = np.column_stack([np.ones(len(s)), [np.log(r["W"] / r["gap_singlet"]) for r in s],
                                 [np.log(1 / r["gamma"]) for r in s]])
            y = np.log([r["cost"] for r in s])
            c, *_ = np.linalg.lstsq(X, y, rcond=None)
            out.append(dict(variant=v, N=N, eps=eps, protocol=protocol, a_gap=c[1], b_gamma=c[2], C0=np.exp(c[0]),
                            n=len(s)))
    return out


if __name__ == "__main__":
    rows = points(load())
    table(rows, "results/su2_compare.csv")
    for e in (1e-2, 1e-3):
        for pr in ("rus", "aa"):
            summary(rows, e, pr)
    fits = [f for e in (1e-2, 1e-3) for pr in ("rus", "aa") for f in fit_exponents(rows, e, pr)]
    table(fits, "results/su2_fits.csv")
    print("\nexponents vs singlet gap: ln C = c + a ln(W/Delta_singlet) + b ln(1/gamma)")
    for f in fits:
        print(f"  {f['variant']} N={f['N']:2d} eps={f['eps']:g} {f['protocol']}: a = {f['a_gap']:.2f}, "
              f"b = {f['b_gamma']:.2f} ({f['n']} pts)")
    for N in sorted({r["N"] for r in rows}):
        fig_cuts(rows, N, out=f"figs/su2_cuts_N{N}.png")
