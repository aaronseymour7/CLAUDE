"""Standalone LaTeX report on the 4 x Lx cylinder data, read only from results/*.json.

    python sbpipe_study/report2d/make_report2d.py && cd sbpipe_study/report2d && pdflatex report2d
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import make_figures as M  # noqa: E402  (style, loaders, frontier helper)

plt = M.plt
OUT = HERE
FIGD = os.path.join(HERE, "figs")
os.makedirs(FIGD, exist_ok=True)
INP = sorted(M.load("ext_inputs2d.json"), key=lambda r: (r["J2"], r["Lx"]))
REG = [r for f in ("ext_regime_cyl.json", "ext_regime_cyl_12.json") for r in (M.load(f) or [])]
FR = [r for r in M.load("emp_frontier.json") if r["kind"] == "cyl"]
ER = [r for r in M.load("emp_regime.json") if r["kind"] == "cyl"]
CACHE = {}
for f in os.listdir(os.path.join(M.RES, "cache2d")):
    if f.endswith(".json"):
        d = json.load(open(os.path.join(M.RES, "cache2d", f)))
        CACHE[(d["Lx"], d["J2"])] = d
CHIS = [2, 4, 8, 16, 32, 64]
import sbpipe  # noqa: E402
import run_extended as RX  # noqa: E402
WEYL = {(r["Lx"], r["J2"]): sbpipe.weyl_bounds(RX.model("cyl", r["Lx"], r["J2"]))[1] for r in INP}
sci = M.sci


def tab_inputs():
    L = []
    for r in INP:
        ed = ""
        if "ed_E0" in r:
            ed = f"${sci(abs(r['dE0']), 0)}$ / ${sci(abs(r['dE1']), 0)}$"
        L.append(f"$4\\times{r['Lx']}$ & {r['N']} & {r['J2']:g} & {r['E0']:.6f} & {r['E1']:.6f} & {r['gap']:.4f} & "
                 f"${sci(r['E0_err'])}$ & ${sci(r['E1_err'])}$ & {r.get('weyl_hi', WEYL[(r['Lx'], r['J2'])]):.1f} & "
                 f"{r['dmrg_time'] / 3600:.2f} & {ed or '--'} \\\\")
    return "\n".join(L)


def tab_conv():
    L = []
    for r in INP:
        h = {c: e for c, e, _ in r["ref_history"]}
        E = h[256]
        cells = [f"${sci(abs(h[c] - E), 0)}$" if c != 256 else "ref" for c in (16, 32, 64, 128, 256)]
        L.append(f"$4\\times{r['Lx']}$ & {r['J2']:g} & " + " & ".join(cells) + " \\\\")
    return "\n".join(L)


def tab_overlap():
    L = []
    for r in INP:
        a = CACHE[(r["Lx"], r["J2"])]["ansatz"]
        cells = []
        for c in CHIS:
            inf = 1 - a[str(c)]["gamma"]
            cells.append(f"${sci(inf, 1)}$" if inf > 1e-9 else "$<10^{-9}$")
        L.append(f"$4\\times{r['Lx']}$ & {r['J2']:g} & " + " & ".join(cells) + " \\\\")
    return "\n".join(L)


def tab_cost():
    L = []
    for Lx in (3, 4, 6, 8, 12):
        a = CACHE[(Lx, 0.0)]["ansatz"]
        L.append(f"$4\\times{Lx}$ & " + " & ".join(f"${sci(a[str(c)]['cost']['cnot'])}$" for c in CHIS) + " \\\\")
    return "\n".join(L)


def tab_regime():
    from collections import defaultdict
    g = defaultdict(list)
    for r in REG:
        g[(r["N"], r["J2"], r["eps"])].append(r)
    L = []
    for J2 in (0.0, 0.3):
        for N in (12, 16, 24, 32, 48):
            cells = []
            for eps in (1e-2, 1e-3, 1e-4):
                v = g[(N, J2, eps)]
                dr = [r for r in v if r["schedule"] == "auto" and 1 - r["gamma"] <= eps]
                fl = [r for r in v if r["plan"] and r["plan"]["m"] > 0]
                db = min(dr, key=lambda r: r["trial_cost"]["cnot"]) if dr else None
                fb = min(fl, key=lambda r: r["plan"]["res_expected_cnot"]) if fl else None
                ef = [r for r in ER if r["N"] == N and r["J2"] == J2 and r["eps"] == eps and not r.get("mps_alone")
                      and r.get("cnot")]
                eb = min(ef, key=lambda r: r["cnot"]) if ef else None
                cells.append(f"${sci(db['trial_cost']['cnot'])}$ ({db['chi']})" if db else "--")
                cells.append(f"${sci(fb['plan']['res_expected_cnot'])}$ ({fb['chi']})" if fb else "--")
                cells.append(f"${sci(eb['cnot'])}$ ({eb['chi']})" if eb else "--")
            L.append(f"$4\\times{N // 4}$ & {J2:g} & " + " & ".join(cells) + " \\\\")
    return "\n".join(L)


def tab_sim():
    L = []
    for r in sorted([r for r in ER if r.get("sim")], key=lambda r: (r["N"], r["chi"], -r["eps"])):
        s = r["sim"]
        L.append(f"$4\\times{r['N'] // 4}$ & {r['J2']:g} & {r['chi']} & {r['gamma']:.4f} & ${sci(r['eps'], 0)}$ & "
                 f"{r['d_model']} & {s['d']} & ${sci(s['infidelity'])}$ & {s['p_succ']:.3f} \\\\")
    return "\n".join(L), sum(1 for r in ER if r.get("sim")), sum(1 for r in ER if r.get("sim") and
                                                                    r["sim"]["infidelity"] <= r["eps"] and
                                                                    r["sim"]["d"] == r["d_model"])


def crossover_rows():
    rows = []
    for J2 in (0.0, 0.3):
        for N in (12, 16, 24, 32, 48):
            g, m, f, _, _ = M._frontiers(FR, REG, "cyl", N, J2)
            ch = np.isfinite(m) & np.isfinite(f) & (f < m)
            iv, start = [], None
            for i in range(len(g)):
                if ch[i] and start is None:
                    start = i
                if start is not None and (not ch[i] or i == len(g) - 1):
                    end = i - 1 if not ch[i] else i
                    if g[start] / g[end] > 1.5:
                        sel = slice(start, end + 1)
                        iv.append((g[start], g[end], float(np.max(m[sel] / f[sel]))))
                    start = None
            none = ~np.isfinite(m)
            lim = g[none].max() if none.any() else None
            rows.append((N, J2, iv, lim))
    L = []
    for N, J2, iv, lim in rows:
        w = "; ".join(f"${sci(a, 0)}$--${sci(b, 0)}$ (${x:.1f}\\times$)" for a, b, x in iv) or "none"
        L.append(f"$4\\times{N // 4}$ & {J2:g} & {w} & " + (f"${sci(lim, 0)}$" if lim else "--") + " \\\\")
    return "\n".join(L)


# ----------------------------------------------------------------------------- figures
def fig_overlap():
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.5))
    cols = {3: "#9aa0a6", 4: "#1baf7a", 6: "#eb6834", 8: "#2a78d6", 12: "#0b0b0b"}
    for Lx, col in cols.items():
        for J2, ls in ((0.0, "-"), (0.3, "--")):
            a = CACHE[(Lx, J2)]["ansatz"]
            inf = [max(1 - a[str(c)]["gamma"], 1e-9) for c in CHIS]
            ax[0].plot(CHIS, inf, "o" + ls, color=col, ms=3, lw=1.2,
                       label=f"$4\\times{Lx}$" if J2 == 0 else None)
        r = [x for x in INP if x["Lx"] == Lx and x["J2"] == 0.0][0]
        h = {c: e for c, e, _ in r["ref_history"]}
        ax[1].plot([16, 32, 64, 128], [max(abs(h[c] - h[256]) / r["N"], 1e-12) for c in (16, 32, 64, 128)],
                   "o-", color=col, ms=3, lw=1.2)
    for a in ax[:2]:
        a.set_xscale("log", base=2); a.set_yscale("log"); a.set_xlabel(r"bond dimension $\chi$")
    ax[0].set_ylabel(r"$1-\gamma$"); ax[0].set_title(r"(a) trial infidelity (solid $J_2=0$, dashed 0.3)",
                                                    loc="left", fontsize=7.5)
    ax[0].legend(fontsize=6, loc="lower left")
    ax[1].set_ylabel(r"$|E_0(\chi)-E_0(256)|/N$"); ax[1].set_title("(b) DMRG convergence, $J_2=0$", loc="left",
                                                                     fontsize=7.5)
    for J2, mk in ((0.0, "o"), (0.3, "s")):
        v = [x for x in INP if x["J2"] == J2]
        ax[2].plot([x["N"] for x in v], [x["gap"] for x in v], mk + "-", color=M.QETU_COL if J2 == 0 else M.SCHED_COL["dyadic"],
                   ms=3.5, label=f"$J_2={J2:g}$")
    ax[2].set_xlabel("$N$"); ax[2].set_ylabel("gap $E_1-E_0$"); ax[2].set_title("(c) triplet gap", loc="left", fontsize=7.5)
    ax[2].legend(fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGD, "fig_inputs.pdf"))
    plt.close(fig)


def fig_crossover():
    fig, axs = plt.subplots(2, 5, figsize=(7.4, 3.6), sharey=True, sharex=True)
    MPS_C, EMP_C, CERT_C = M.SCHED_COL["optimized"], M.QETU_COL, M.SCHED_COL["dyadic"]
    for row, J2 in enumerate((0.0, 0.3)):
        for a, N in zip(axs[row], (12, 16, 24, 32, 48)):
            g, m, f, cert, mps = M._frontiers(FR, REG, "cyl", N, J2)
            a.plot(g, m, color=MPS_C, lw=1.4, drawstyle="steps-pre")
            for inf, c, chi in mps:
                if 1e-7 <= inf <= 1e-1:
                    a.plot(inf, c, "o", color=MPS_C, ms=2.8)
            a.plot(g, f, "-", color=EMP_C, lw=1.4)
            a.plot([e for e, _ in cert], [c for _, c in cert], "D--", color=CERT_C, ms=2.8, lw=0.8)
            ch = np.isfinite(m) & np.isfinite(f) & (f < m)
            a.fill_between(g, 1e1, 1e8, where=ch, color=EMP_C, alpha=0.18, lw=0, step="pre")
            a.fill_between(g, 1e1, 1e8, where=~np.isfinite(m), step="pre", facecolor=M.GRID, alpha=0.7, lw=0,
                           hatch="////", edgecolor="white")
            a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1.5e-1, 7e-8); a.set_ylim(1e2, 3e7)
            a.set_xticks([1e-1, 1e-4, 1e-7])
            a.set_title(f"$4\\times{N // 4}$, $J_2={J2:g}$", loc="left", fontsize=7)
            if row == 1:
                a.set_xlabel("$1-F$")
        axs[row][0].set_ylabel("CNOTs")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGD, "fig_crossover2d.pdf"))
    plt.close(fig)


def main():
    fig_overlap(); fig_crossover()
    sim_tab, n_sim, n_ok = tab_sim()
    tex = open(os.path.join(HERE, "report2d_template.tex")).read()
    for k, v in dict(TAB_INPUTS=tab_inputs(), TAB_CONV=tab_conv(), TAB_OVERLAP=tab_overlap(), TAB_COST=tab_cost(),
                     TAB_REGIME=tab_regime(), TAB_SIM=sim_tab, NSIM=str(n_sim), NOK=str(n_ok),
                     TAB_CROSS=crossover_rows()).items():
        tex = tex.replace(f"@@{k}@@", v)
    open(os.path.join(HERE, "report2d.tex"), "w").write(tex)


if __name__ == "__main__":
    main()
