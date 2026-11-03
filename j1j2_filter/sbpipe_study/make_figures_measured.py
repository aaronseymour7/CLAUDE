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


def motivation_fig():
    """Chains versus 4-leg cylinders at the same N (J2 = 0, DMRG): residual weight 1 - gamma of the optimal MPS
    against bond dimension, and the gap against N. Sources: results/inputs.json, results/ext_inputs2d.json."""
    ch = {r["N"]: r for r in M.load("inputs.json") if r["J2"] == 0.0}
    cy = {r["N"]: r for r in M.load("ext_inputs2d.json") if r["J2"] == 0.0}
    col = {12: "#7fb0e8", 16: "#4f8fd8", 24: "#2a6bc0", 32: "#14509e", 48: "#0b2f66"}
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6))
    a = axs[0]
    for N, c in col.items():
        for R, ls, mk in ((cy, "-", "o"), (ch, "--", "s")):
            if N not in R:
                continue
            g = {int(k): 1 - v for k, v in R[N]["gamma"].items() if int(k) >= 2 and 1 - v > 1e-12}
            x = sorted(g)
            a.plot(x, [g[k] for k in x], ls=ls, marker=mk, color=c, ms=3.5, lw=1.5, mfc=c if ls == "-" else "white")
    a.set_xscale("log", base=2); a.set_yscale("log")
    a.set_xlabel(r"bond dimension $\chi$"); a.set_ylabel(r"residual weight $1-\gamma$")
    a.yaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
    h = [plt.Line2D([], [], color=c, lw=1.5, label=f"$N={N}$") for N, c in col.items()]
    h += [plt.Line2D([], [], color=M.INK2, ls="-", marker="o", ms=3.5, label=r"cylinder $4\times N/4$"),
          plt.Line2D([], [], color=M.INK2, ls="--", marker="s", mfc="white", ms=3.5, label="chain")]
    a.legend(handles=h, fontsize=6.2, loc="lower left", ncol=2)
    a.set_title("(a) MPS trial quality", loc="left", fontsize=7.5)
    b = axs[1]
    for R, ls, mk, lab in ((cy, "-", "o", r"cylinder $4\times N/4$"), (ch, "--", "s", "chain")):
        N = sorted(R)
        b.plot(N, [R[n]["gap"] for n in N], ls=ls, marker=mk, color=M.INK2, ms=3.5, lw=1.5,
               mfc=M.INK2 if ls == "-" else "white", label=lab)
    b.set_xscale("log", base=2); b.set_yscale("log")
    b.set_xticks([8, 12, 16, 24, 32, 48]); b.set_xticklabels(["8", "12", "16", "24", "32", "48"])
    b.xaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
    b.yaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
    b.set_yticks([0.1, 0.2, 0.3, 0.5]); b.set_yticklabels(["0.1", "0.2", "0.3", "0.5"])
    b.set_xlabel(r"sites $N$"); b.set_ylabel(r"gap $E_1-E_0$ ($J_1$)")
    b.legend(fontsize=6.5, loc="lower left")
    b.set_title("(b) gap", loc="left", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(M.FIG, "fig_motivation_meas.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    study_fig("A"); study_fig("B"); mechanism_fig(); fit_table(); motivation_fig()
