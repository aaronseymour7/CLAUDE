"""Measured-only variants of the gap/overlap figures and fit table, for paper/paper_measured.tex.

    python sbpipe_study/make_figures_measured.py

Reads results/gapov_filter.json and paper/generated/tab_gapov_fit.tex (from make_figures.py) and writes
paper/figs/fig_filter_{A,B,mech}_meas.pdf and paper/generated/tab_gapov_fit_meas.tex.
"""
import os

import numpy as np

import make_figures as M

plt = M.plt
D = M.load("gapov_filter.json")
NS = sorted({r["N"] for r in D})
NCOL = M.GO_NCOL
PANELS = (("cnot", False, "CNOTs per shot"), ("rz", False, "$R_z$ rotations per shot"),
          ("time", False, "evolution time per shot ($1/J_1$)"), ("cnot", True, "CNOTs per success"))


def L_of(r):
    return np.log((1 - r["gamma"]) / (r["gamma"] * r["eps"]))


def val(r, meth, met, per_success=False):
    m = r[meth]
    return m[met] / (m["p"] if per_success else 1.0)


def study_fig(tag):
    eps = 1e-3
    if tag == "A":
        sel = lambda r: abs(r["gap_target"] - 0.40) < 1e-9 and r["eps"] == eps
        xk, xl = "gamma", r"trial overlap $\gamma$"
    else:
        sel = lambda r: r["gamma_target"] == 0.10 and r["eps"] == eps
        xk, xl = "gap", r"gap $\Delta$"
    fig, axs = plt.subplots(3, 4, figsize=(7.2, 6.3), sharex=True)
    for row, fam in enumerate(("dyadic", "optimized", "qetu")):
        for a, (met, ps, lab) in zip(axs[row], PANELS):
            for N in NS:
                R = sorted([r for r in D if r["N"] == N and sel(r)], key=lambda r: r[xk])
                if R:
                    a.plot([r[xk] for r in R], [val(r, "meas_" + fam, met, ps) for r in R], "o-", color=NCOL[N],
                           ms=2.8, lw=1.3)
            a.set_xscale("log"); a.set_yscale("log")
            a.yaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
            if row == 0:
                a.set_title(lab, loc="left", fontsize=7.5)
            if row == 2:
                a.set_xlabel(xl)
            if tag == "B":
                a.set_xticks([0.3, 0.4, 0.5]); a.set_xticklabels(["0.3", "0.4", "0.5"])
                a.xaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
        axs[row, 0].set_ylabel({"dyadic": "dyadic SBC", "optimized": "optimized SBC", "qetu": "QETU"}[fam])
    h = [plt.Line2D([], [], color=NCOL[N], lw=1.5, label=f"$N={N}$") for N in NS]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=6.8, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(M.FIG, f"fig_filter_{tag}_meas.pdf"))
    plt.close(fig)


def mechanism_fig():
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    sty = (("meas_dyadic", "dyadic SBC", M.SCHED_COL["dyadic"], "o"),
           ("meas_optimized", "optimized SBC", M.SCHED_COL["optimized"], "s"),
           ("meas_qetu", "QETU", M.QETU_COL, "^"))
    x = np.array([L_of(r) / r["gap"] for r in D])
    e = np.array([r["eps"] for r in D]) * np.exp(np.random.default_rng(1).uniform(-0.15, 0.15, len(D)))
    for meth, name, col, mk in sty:
        axs[0].plot(x, [r[meth]["time"] for r in D], mk, mfc=col, mec=col, ms=2.4, alpha=0.7, label=name)
        axs[2].plot(x, [r[meth]["steps"] for r in D], mk, mfc=col, mec=col, ms=2.4, alpha=0.7)
        axs[1].plot(e, [r[meth]["dt"] for r in D], mk, mfc=col, mec=col, ms=2.4, alpha=0.6)
    for a in axs:
        a.set_xscale("log"); a.set_yscale("log")
        a.xaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
        a.yaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
    axs[0].set_xlabel(r"$L/\Delta$"); axs[0].set_ylabel(r"evolution time ($1/J_1$)")
    axs[1].set_xlabel(r"target $\varepsilon$"); axs[1].set_ylabel(r"Trotter step $dt$ ($1/J_1$)")
    axs[1].invert_xaxis()
    axs[2].set_xlabel(r"$L/\Delta$"); axs[2].set_ylabel("Trotter steps per shot")
    axs[0].legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(M.FIG, "fig_filter_mech_meas.pdf"))
    plt.close(fig)


def fit_table():
    rows = [l for l in open(os.path.join(M.GEN, "tab_gapov_fit.tex")) if "measured" in l]
    with open(os.path.join(M.GEN, "tab_gapov_fit_meas.tex"), "w") as fh:
        fh.write("".join(rows))


if __name__ == "__main__":
    study_fig("A"); study_fig("B"); mechanism_fig(); fit_table()
