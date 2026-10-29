"""Filter-only analysis of the gap/overlap study: CNOTs, Rz rotations, evolution time, Trotter steps.

    python sbpipe_study/derive_gapov_filter.py
    python sbpipe_study/gapov_report/make_report_filter.py && cd sbpipe_study/gapov_report && pdflatex report_filter
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_figures as M  # noqa: E402

plt = M.plt
FIGD = os.path.join(HERE, "figs")
os.makedirs(FIGD, exist_ok=True)
D = M.load("gapov_filter.json")
NS = sorted({r["N"] for r in D})
NCOL = M.GO_NCOL
METHODS = (("meas_dyadic", "dyadic, measured"), ("meas_qetu", "QETU, measured"),
           ("cert_dyadic", "dyadic, certified"), ("cert_qetu", "QETU, certified"))
sci = M.sci


def L_of(r):
    return np.log((1 - r["gamma"]) / (r["gamma"] * r["eps"]))


def val(r, meth, met, per_success=False):
    m = r[meth]
    return m[met] / (m["p"] if per_success else 1.0)


def fit(meth, met, per_success=False, cols=("N", "gap", "L", "eps")):
    X, y = [], []
    for r in D:
        f = dict(N=np.log(r["N"]), gap=np.log(r["gap"]), L=np.log(L_of(r)), eps=np.log(r["eps"]),
                 gamma=np.log(r["gamma"]))
        X.append([1] + [f[c] for c in cols])
        y.append(np.log(val(r, meth, met, per_success)))
    X, y = np.array(X), np.array(y)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    return c, float(np.exp(np.sqrt(np.mean(res ** 2)))), float(np.exp(np.abs(res).max()))


# ----------------------------------------------------------------------------- figures
PANELS = (("cnot", False, "CNOTs per shot"), ("rz", False, "$R_z$ rotations per shot"),
          ("time", False, "evolution time per shot ($1/J_1$)"), ("cnot", True, "CNOTs per success"))


def study_fig(tag):
    eps = 1e-3
    if tag == "A":
        sel = lambda r: abs(r["gap_target"] - 0.40) < 1e-9 and r["eps"] == eps
        xk, xl = "gamma", r"trial overlap $\gamma$"
    else:
        sel = lambda r: r["gamma_target"] == 0.10 and r["eps"] == eps
        xk, xl = "gap", r"gap $\Delta$"
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.4), sharex=True)
    for row, fam in enumerate(("dyadic", "qetu")):
        for a, (met, ps, lab) in zip(axs[row], PANELS):
            for N in NS:
                R = sorted([r for r in D if r["N"] == N and sel(r)], key=lambda r: r[xk])
                if not R:
                    continue
                x = [r[xk] for r in R]
                a.plot(x, [val(r, "meas_" + fam, met, ps) for r in R], "o-", color=NCOL[N], ms=2.8, lw=1.3)
                a.plot(x, [val(r, "cert_" + fam, met, ps) for r in R], "D--", color=NCOL[N], ms=2.2, lw=0.8)
            a.set_xscale("log"); a.set_yscale("log")
            a.yaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
            if row == 0:
                a.set_title(lab, loc="left", fontsize=7.5)
            if row == 1:
                a.set_xlabel(xl)
            if tag == "B":
                a.set_xticks([0.3, 0.4, 0.5]); a.set_xticklabels(["0.3", "0.4", "0.5"])
                a.xaxis.set_minor_formatter(M.matplotlib.ticker.NullFormatter())
        axs[row, 0].set_ylabel("dyadic SBC" if fam == "dyadic" else "QETU")
    h = [plt.Line2D([], [], color=NCOL[N], lw=1.5, label=f"$N={N}$") for N in NS]
    h += [plt.Line2D([], [], color=M.INK2, marker="o", ms=2.8, lw=1.3, label="measured"),
          plt.Line2D([], [], color=M.INK2, marker="D", ls="--", ms=2.2, lw=0.8, label="certified")]
    fig.legend(handles=h, loc="lower center", ncol=6, fontsize=6.8, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(FIGD, f"fig_filter_{tag}.pdf"))
    fig.savefig(os.path.join(FIGD, f"fig_filter_{tag}.png"), dpi=200)
    plt.close(fig)


def mechanism_fig():
    """(a) evolution time collapses on L/gap; (b) Trotter step against eps; (c) steps against L/gap."""
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    sty = {"meas_dyadic": (M.SCHED_COL["dyadic"], "o"), "meas_qetu": (M.QETU_COL, "^"),
           "cert_dyadic": (M.SCHED_COL["dyadic"], "D"), "cert_qetu": (M.QETU_COL, "s")}
    for meth, name in METHODS:
        col, mk = sty[meth]
        fill = col if meth.startswith("meas") else "none"
        x = np.array([L_of(r) / r["gap"] for r in D])
        axs[0].plot(x, [r[meth]["time"] for r in D], mk, mfc=fill, mec=col, ms=2.4, alpha=0.7, label=name)
        axs[2].plot(x, [r[meth]["steps"] for r in D], mk, mfc=fill, mec=col, ms=2.4, alpha=0.7)
        e = np.array([r["eps"] for r in D]) * np.exp(np.random.default_rng(1).uniform(-0.15, 0.15, len(D)))
        axs[1].plot(e, [r[meth]["dt"] for r in D], mk, mfc=fill, mec=col, ms=2.4, alpha=0.6)
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
    fig.savefig(os.path.join(FIGD, "fig_filter_mech.pdf"))
    fig.savefig(os.path.join(FIGD, "fig_filter_mech.png"), dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------------------- tables
def tab_fit():
    L = []
    rows = (("time", False, "evolution time"), ("dt", False, "Trotter step $dt$"), ("steps", False, "Trotter steps"),
            ("cnot", False, "CNOTs"), ("rz", False, "$R_z$"))
    for meth, name in METHODS:
        L.append(f"\\multicolumn{{8}}{{l}}{{\\emph{{{name}, per shot}}}} \\\\")
        for met, ps, lab in rows:
            c, rms, mx = fit(meth, met, ps)
            L.append(f"\\quad {lab} & ${sci(np.exp(c[0]), 2)}$ & {c[1]:.2f} & {-c[2]:.2f} & {c[3]:.2f} & "
                     f"{-c[4]:.2f} & {rms:.2f} & {mx:.1f} \\\\")
    return "\n".join(L)


def tab_succ():
    L = []
    for meth, name in METHODS:
        for met, lab in (("cnot", "CNOTs"), ("rz", "$R_z$"), ("time", "time")):
            c, rms, mx = fit(meth, met, True, cols=("N", "gap", "L", "eps", "gamma"))
            L.append(f"{name} & {lab} & ${sci(np.exp(c[0]), 2)}$ & {c[1]:.2f} & {-c[2]:.2f} & {c[3]:.2f} & "
                     f"{-c[4]:.2f} & {-c[5]:.2f} & {rms:.2f} \\\\")
        L.append("\\addlinespace")
    return "\n".join(L[:-1])


def ratios_A():
    """per-shot growth when gamma drops 0.9 -> 0.03 at gap 0.40 (range over N, eps)."""
    out = {}
    for meth, _ in METHODS:
        for met in ("cnot", "time"):
            v = []
            for e in (1e-2, 1e-3, 1e-4):
                for N in NS:
                    R = {r["gamma_target"]: r for r in D if r["N"] == N and abs(r["gap_target"] - 0.4) < 1e-9
                         and r["eps"] == e}
                    if 0.03 in R and 0.9 in R:
                        v.append(R[0.03][meth][met] / R[0.9][meth][met])
            out[(meth, met)] = (min(v), max(v))
    return out


def ratios_B():
    """per-shot growth when the gap closes from its largest to its smallest value at gamma 0.1, N = 12, 14."""
    out = {}
    for meth, _ in METHODS:
        for met in ("cnot", "time"):
            v = []
            for e in (1e-2, 1e-3, 1e-4):
                for N in (12, 14):
                    R = sorted([r for r in D if r["N"] == N and r["gamma_target"] == 0.1 and r["eps"] == e],
                               key=lambda r: r["gap"])
                    v.append(R[0][meth][met] / R[-1][meth][met])
            out[(meth, met)] = (min(v), max(v))
    return out


def ratio_range(a, b, met):
    v = [r[a][met] / r[b][met] for r in D]
    return min(v), float(np.median(v)), max(v)


def rng(t):
    return f"{t[0]:.1f}--{t[1]:.1f}$\\times$"


study_fig("A"); study_fig("B"); mechanism_fig()
RA, RB = ratios_A(), ratios_B()
SUB = dict(FIT=tab_fit(), SUCC=tab_succ(), NCASES=str(len(D)))
for meth, _ in METHODS:
    t = "".join(w[0] for w in meth.split("_")).upper()
    for met in ("cnot", "time"):
        SUB[f"A{t}{met.upper()}"] = rng(RA[(meth, met)])
        SUB[f"B{t}{met.upper()}"] = rng(RB[(meth, met)])
for fam in ("dyadic", "qetu"):
    for met in ("time", "steps"):
        lo, med, hi = ratio_range("cert_" + fam, "meas_" + fam, met)
        SUB[f"CM{fam.upper()}{met.upper()}"] = f"{lo:.1f}--{hi:.0f}$\\times$ (median {med:.1f}$\\times$)"
lo, med, hi = ratio_range("meas_qetu", "meas_dyadic", "time")
SUB["QDTIME"] = f"{med:.2f}"
lo, med, hi = ratio_range("meas_qetu", "meas_dyadic", "steps")
SUB["QDSTEPS"] = f"{med:.2f}"
rz = [r["meas_dyadic"]["rz"] / r["meas_dyadic"]["cnot"] for r in D]
SUB["RZCNOT"] = f"{np.median(rz):.2f}"
cps = [r["meas_dyadic"]["cnot"] / r["meas_dyadic"]["steps"] / r["N"] for r in D]
SUB["CPS"] = f"{np.median(cps):.0f}"
for e in (1e-2, 1e-3, 1e-4):
    v = [r["meas_dyadic"]["dt"] for r in D if r["eps"] == e]
    SUB[f"DTD{int(-np.log10(e))}"] = f"{np.median(v):.2f}"
    v = [r["cert_dyadic"]["dt"] for r in D if r["eps"] == e]
    SUB[f"DTC{int(-np.log10(e))}"] = f"{np.median(v):.3f}"
v = [r["meas_qetu"]["dt"] for r in D]
SUB["DTQ"] = f"{min(v):.2f}--{max(v):.2f}"
txt = open(os.path.join(HERE, "report_filter_template.tex")).read()
for k in sorted(SUB, key=len, reverse=True):
    txt = txt.replace(f"@{k}@", SUB[k])
open(os.path.join(HERE, "report_filter.tex"), "w").write(txt)
print({k: v for k, v in SUB.items() if len(v) < 60})
