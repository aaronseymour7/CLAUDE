"""N-by-epsilon CNOT maps for the dyadic SBC filter at fixed gap (0.40) and overlap (0.1, 0.5).

    python sbpipe_study/gapov_report/make_heat.py
Reads results/gapov.jsonl and results/gapov_heat.jsonl. Filter-only CNOTs (trial circuit excluded).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_figures as M  # noqa: E402
from derive_gapov_filter import counts  # noqa: E402
from sbpipe import j1j2_chain  # noqa: E402

plt = M.plt
FIGD = os.path.join(HERE, "figs")
EPS = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
GAMS = (0.1, 0.5)


FILTERS = ("dyadic", "optimized")


def rows():
    """one dict per (N, eps, gamma) point at gap 0.40, dyadic from gapov / gapov_heat, optimized (with retries)
    attached by derive_gapov_filter.load_cases."""
    from derive_gapov_filter import load_cases
    P = {}
    for r in load_cases(("gapov.jsonl", "gapov_heat.jsonl")):
        if abs(r["gap_target"] - 0.40) > 1e-9 or r["gamma_target"] not in GAMS or r["eps"] not in EPS:
            continue
        P[(r["N"], r["gamma_target"], r["eps"])] = r
    return list(P.values())


def grid(R, Ns):
    """filter-only CNOTs per shot and per success on the (N, eps) grid, keyed (filter, meas|cert, shot|succ)."""
    H = {}
    out = {(f, w, k): {g: np.full((len(Ns), len(EPS)), np.nan) for g in GAMS}
           for f in FILTERS for w in ("meas", "cert") for k in ("shot", "succ")}
    for r in R:
        i, j = Ns.index(r["N"]), EPS.index(r["eps"])
        key = (r["N"], r["J2"])
        if key not in H:
            H[key] = j1j2_chain(r["N"], 1.0, r["J2"])
        g = r["gamma_target"]
        for f in FILTERS:
            md, cd = r.get(f"meas_{f}"), r.get(f"cert_{f}")
            if md:
                c = md["cnot_filter"] if "cnot_filter" in md else counts(H[key], md["steps"], md["m"], 1)[0]
                out[(f, "meas", "shot")][g][i, j] = c
                out[(f, "meas", "succ")][g][i, j] = c / md["p"]
            if cd:
                c, _ = counts(H[key], cd["steps"], cd["m"], 1)
                p = (c + r["trial_cnot"]) / cd["cnot"]
                out[(f, "cert", "shot")][g][i, j] = c
                out[(f, "cert", "succ")][g][i, j] = c / p
    return out


def _lab(v):
    e = int(np.floor(np.log10(v)))
    return f"${v / 10 ** e:.1f}\\times10^{{{e}}}$"


def _edges(Ns):
    xe = np.log10(EPS)
    xedge = np.concatenate([[xe[0] + 0.25], (xe[:-1] + xe[1:]) / 2, [xe[-1] - 0.25]])
    yedge = np.concatenate([[Ns[0] - 1], (np.array(Ns[:-1]) + np.array(Ns[1:])) / 2, [Ns[-1] + 1]])
    return xe, xedge, yedge


def _axes(a, xe, xedge, Ns):
    a.set_xlim(xedge[0], xedge[-1])
    a.set_xticks(xe)
    a.set_xticklabels([f"$10^{{{int(v)}}}$" if v == int(v) else "" for v in xe])
    a.set_yticks(Ns)
    a.grid(False)


def heat_fig(G, Ns, w, kind, fname, title, cols=None):
    """cols: list of (filter, meas|cert) panels; default both filters at the given w."""
    cols = cols or [(f, w) for f in FILTERS]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True, sharey=True)
    xe, xedge, yedge = _edges(Ns)
    for col, (f, w) in enumerate(cols):
        vals = np.concatenate([G[(f, w, kind)][g].ravel() for g in GAMS])
        lo, hi = np.log10(np.nanmin(vals)), np.log10(np.nanmax(vals))
        if len({c[1] for c in cols}) == 1:          # same quantity in both columns: one shared scale
            vals = np.concatenate([G[(ff, ww, kind)][g].ravel() for ff, ww in cols for g in GAMS])
            lo, hi = np.log10(np.nanmin(vals)), np.log10(np.nanmax(vals))
        levels = np.arange(np.floor(lo * 4) / 4, hi + 0.01, 0.25)
        for row, g in enumerate(GAMS):
            a = axs[row, col]
            Z = np.log10(G[(f, w, kind)][g])
            pc = a.pcolormesh(xedge, yedge, Z, cmap="viridis", vmin=lo, vmax=hi, shading="flat")
            if np.isfinite(Z).all():
                a.contour(xe, Ns, Z, levels=levels, colors="white", linewidths=0.6, alpha=0.7)
            for i, N in enumerate(Ns):
                for j, e in enumerate(xe):
                    if np.isfinite(Z[i, j]):
                        a.text(e, N, _lab(10 ** Z[i, j]), ha="center", va="center", fontsize=5.6,
                               color="white" if Z[i, j] < (lo + hi) / 2 else "black")
                    else:
                        a.text(e, N, "pending", ha="center", va="center", fontsize=5, color=M.INK2)
            wn = "measured" if w == "meas" else "certified"
            a.set_title(f"{f} SBC, {wn}, $\\gamma={g}$", loc="left", fontsize=7.5)
            _axes(a, xe, xedge, Ns)
        if len({c[1] for c in cols}) > 1 or col == len(cols) - 1:
            cb = fig.colorbar(pc, ax=axs[:, col] if len({c[1] for c in cols}) > 1 else axs,
                              orientation="horizontal", fraction=0.04, pad=0.1, aspect=40)
            cb.set_label(r"$\log_{10}$ CNOTs (contours every 0.25)", fontsize=7)
    for a in axs[1]:
        a.set_xlabel(r"target infidelity $\varepsilon$")
    for a in axs[:, 0]:
        a.set_ylabel("qubits $N$")
    fig.suptitle(title, fontsize=8.5, x=0.01, ha="left")
    fig.savefig(os.path.join(FIGD, fname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGD, fname + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def ratio_fig(G, Ns, fname):
    """optimized / dyadic CNOTs per success, measured and certified."""
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.4), sharex=True, sharey=True)
    xe, xedge, yedge = _edges(Ns)
    for col, w in enumerate(("meas", "cert")):
        for row, g in enumerate(GAMS):
            a = axs[row, col]
            Z = np.log2(G[("optimized", w, "succ")][g] / G[("dyadic", w, "succ")][g])
            pc = a.pcolormesh(xedge, yedge, Z, cmap="RdBu_r", vmin=-2, vmax=2, shading="flat")
            for i, N in enumerate(Ns):
                for j, e in enumerate(xe):
                    if np.isfinite(Z[i, j]):
                        a.text(e, N, f"{2 ** Z[i, j]:.2f}", ha="center", va="center", fontsize=6,
                               color="white" if abs(Z[i, j]) > 1.2 else "black")
            a.set_title(f"{'measured' if w == 'meas' else 'certified'}, $\\gamma={g}$", loc="left", fontsize=7.5)
            _axes(a, xe, xedge, Ns)
    cb = fig.colorbar(pc, ax=axs, orientation="horizontal", fraction=0.04, pad=0.1, aspect=50,
                      ticks=[-2, -1, 0, 1, 2])
    cb.ax.set_xticklabels(["1/4", "1/2", "1", "2", "4"])
    cb.set_label("optimized / dyadic, CNOTs per success (blue: optimized cheaper)", fontsize=7)
    for a in axs[1]:
        a.set_xlabel(r"target infidelity $\varepsilon$")
    for a in axs[:, 0]:
        a.set_ylabel("qubits $N$")
    fig.suptitle(r"Optimized versus dyadic SBC filter, gap fixed at $\Delta=0.40$", fontsize=8.5, x=0.01, ha="left")
    fig.savefig(os.path.join(FIGD, fname + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIGD, fname + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    R = rows()
    Ns = sorted({r["N"] for r in R})
    G = grid(R, Ns)
    for w, wn in (("meas", "measured by simulation"), ("cert", "certified by Theorem 1")):
        for kind, kn in (("succ", "per success"), ("shot", "per shot")):
            heat_fig(G, Ns, w, kind, f"fig_heat_{w}_{kind}",
                     f"Filter CNOTs {kn}, {wn}; gap fixed at $\\Delta=0.40$ for every $N$")
    ratio_fig(G, Ns, "fig_heat_ratio")
    for kind, kn in (("succ", "per success"), ("shot", "per shot")):
        heat_fig(G, Ns, None, kind, f"fig_heat_dyadic_{kind}",
                 f"Dyadic SBC filter, filter CNOTs {kn}; gap fixed at $\\Delta=0.40$ for every $N$",
                 cols=[("dyadic", "meas"), ("dyadic", "cert")])
    print(len(R), "points;", {"/".join(k): int(sum(np.isfinite(v[g]).sum() for g in GAMS)) for k, v in G.items()})
