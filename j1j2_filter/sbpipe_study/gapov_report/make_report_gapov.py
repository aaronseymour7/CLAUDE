"""Standalone LaTeX report on the gap-vs-overlap study, read only from results/gapov*.json(l).

    python sbpipe_study/gapov_report/make_report_gapov.py && cd sbpipe_study/gapov_report && pdflatex report_gapov
"""
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_figures as M  # noqa: E402

FIGD = os.path.join(HERE, "figs")
os.makedirs(FIGD, exist_ok=True)
T = M.load("gapov_tune.json")
G = M._gapov()
M.gapov()
for f in ("fig_gapov_A", "fig_gapov_B", "fig_gapov_fit"):
    shutil.copy(os.path.join(M.FIG, f + ".pdf"), FIGD)
sci = M.sci
NS = sorted({r["N"] for r in G})
GAPS = sorted({p["gap"] for p in T["points"]})
GAMMAS = sorted({t["gamma"] for t in T["theta"]})


def tab_tune():
    J = {(p["N"], p["gap"]): p["J2"] for p in T["points"]}
    L = []
    for g in GAPS:
        L.append(f"{g:.2f} & " + " & ".join(f"{J[(N, g)]:.4f}" if (N, g) in J else "--" for N in NS) + " \\\\")
    return "\n".join(L)


def tab_theta():
    J = {(p["N"], p["gap"]): p["J2"] for p in T["points"]}
    L = []
    for gm in GAMMAS:
        cells = []
        for N in NS:
            th = [t["theta"] for t in T["theta"] if t["N"] == N and abs(t["J2"] - J.get((N, 0.40), -1)) < 1e-9
                  and t["gamma"] == gm]
            cells.append(f"{th[0]:.3f}" if th else "--")
        L.append(f"{gm:g} & " + " & ".join(cells) + " \\\\")
    mps = []
    for N in NS:
        q = [t["gamma_mps"] for t in T["theta"] if t["N"] == N and abs(t["J2"] - J.get((N, 0.40), -1)) < 1e-9]
        mps.append(f"{q[0]:.4f}" if q else "--")
    L.append("\\midrule\n$\\gamma$ of $\\chi=16$ MPS ($\\theta=0$) & " + " & ".join(mps) + " \\\\")
    return "\n".join(L)


def tab_fit():
    L = []
    for key, name in (("meas_dyadic", "dyadic, measured"), ("meas_qetu", "QETU, measured"),
                      ("cert_dyadic", "dyadic, certified"), ("cert_qetu", "QETU, certified")):
        c, rms, mx, X, y = M._go_fit(G, key)
        L.append(f"{name} & ${sci(np.exp(c[0]), 2)}$ & {c[1]:.2f} & {-c[2]:.2f} & {c[3]:.2f} & {-c[4]:.2f} & "
                 f"{rms:.2f}$\\times$ & {mx:.1f}$\\times$ & {len(y)} \\\\")
    return "\n".join(L)


def tab_cases():
    """measured / certified CNOTs per success for the dyadic filter at gap 0.40, eps 1e-3."""
    L = []
    for gm in GAMMAS:
        cells = []
        for N in NS:
            r = [r for r in G if r["N"] == N and abs(r["gap_target"] - 0.40) < 1e-9 and r["gamma_target"] == gm
                 and r["eps"] == 1e-3]
            if r and r[0]["meas_dyadic"] and r[0]["cert_dyadic"]:
                cells.append(f"${sci(r[0]['meas_dyadic']['cnot'])}$ / ${sci(r[0]['cert_dyadic']['cnot'])}$")
            else:
                cells.append("--")
        L.append(f"{gm:g} & " + " & ".join(cells) + " \\\\")
    return "\n".join(L)


def ratios():
    out = {}
    for key in ("dyadic", "qetu"):
        r = [x[f"cert_{key}"]["cnot"] / x[f"meas_{key}"]["cnot"] for x in G if x[f"cert_{key}"] and x[f"meas_{key}"]]
        out[key] = (min(r), float(np.median(r)), max(r))
    s = [x["meas_qetu"]["cnot"] / x["meas_dyadic"]["cnot"] for x in G if x["meas_qetu"] and x["meas_dyadic"]]
    out["sv"] = (min(s), float(np.median(s)), max(s), float(np.mean(np.array(s) > 1)))
    return out


R = ratios()
SUB = dict(
    TUNE=tab_tune(), THETA=tab_theta(), FIT=tab_fit(), CASES=tab_cases(),
    NCOLS="c" * len(NS), NHEAD=" & ".join(f"$N={N}$" for N in NS), NCASES=str(len(G)),
    NMIN=str(min(NS)), NMAX=str(max(NS)),
    CRD=f"{R['dyadic'][0]:.0f}--{R['dyadic'][2]:.0f}$\\times$ (median {R['dyadic'][1]:.0f}$\\times$)",
    CRQ=f"{R['qetu'][0]:.0f}--{R['qetu'][2]:.0f}$\\times$ (median {R['qetu'][1]:.0f}$\\times$)",
    SV=f"{R['sv'][1]:.2f}", SVFRAC=f"{100 * R['sv'][3]:.0f}",
)
for k, v in [(k, M._go_fit(G, k)) for k in ("meas_dyadic", "meas_qetu", "cert_dyadic", "cert_qetu")]:
    c = v[0]
    t = "".join(w[0] for w in k.split("_")).upper()
    SUB[f"{t}A"] = sci(np.exp(c[0]), 2)
    SUB[f"{t}N"], SUB[f"{t}GAP"], SUB[f"{t}L"], SUB[f"{t}GAM"] = (f"{c[1]:.2f}", f"{-c[2]:.2f}", f"{c[3]:.2f}",
                                                                   f"{-c[4]:.2f}")
    SUB[f"{t}RMS"], SUB[f"{t}MAX"] = f"{v[1]:.2f}", f"{v[2]:.1f}"
txt = open(os.path.join(HERE, "report_gapov_template.tex")).read()
for k in sorted(SUB, key=len, reverse=True):
    txt = txt.replace(f"@{k}@", SUB[k])
open(os.path.join(HERE, "report_gapov.tex"), "w").write(txt)
print("wrote report_gapov.tex;", {k: v for k, v in SUB.items() if len(v) < 40})
