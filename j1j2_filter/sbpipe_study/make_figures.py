"""Figures and LaTeX tables/macros for the paper, read only from results/*.json.

    python sbpipe_study/make_figures.py
"""
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "paper", "figs")
GEN = os.path.join(HERE, "paper", "generated")
os.makedirs(FIG, exist_ok=True)
os.makedirs(GEN, exist_ok=True)

# reference categorical slots 1-3 (validated all-pairs, light surface); text in neutral ink
SCHED_COL = {"optimized": "#2a78d6", "dyadic": "#eb6834", "uniform": "#1baf7a"}
QETU_COL = "#7a5af0"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
J2_NAME = {0.0: r"$J_2=0$", 0.2411: r"$J_2=0.2411$", 0.4: r"$J_2=0.4$"}
J2S = [0.0, 0.2411, 0.4]

plt.rcParams.update({
    "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 8.5, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.6, "lines.markersize": 4.5, "savefig.dpi": 300, "font.family": "serif",
    "mathtext.fontset": "cm", "legend.frameon": False,
})


def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


MACROS = []


def macro(name, value):
    MACROS.append(f"\\newcommand{{\\{name}}}{{{value}}}")


def sci(x, d=1):
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    if round(m, d) >= 10:
        m, e = m / 10, e + 1
    return f"{m:.{d}f}\\times10^{{{e}}}" if e not in (0,) else f"{x:.{d}f}"


def num(x):
    return f"{x:,.0f}".replace(",", "{,}")


def fit_power(N, y):
    N, y = np.asarray(N, float), np.asarray(y, float)
    ok = np.isfinite(y) & (y > 0)
    p = np.polyfit(np.log(N[ok]), np.log(y[ok]), 1)
    return p[0]


# ----------------------------------------------------------------------------- inputs table
def inputs_table():
    rows = load("inputs.json")
    if rows is None:
        return
    lines = []
    worst_E0 = worst_E1 = worst_g = 0.0
    for r in rows:
        if "ed_E0" not in r:
            continue
        worst_E0 = max(worst_E0, abs(r["dE0"]))
        worst_E1 = max(worst_E1, abs(r["dE1"]))
        for c, g in r["ed_gamma"].items():
            worst_g = max(worst_g, abs(g - r["gamma"][c]))
        worst_W = max(abs(r["weyl_hi"] - r["ed_Etop"]) for r in rows if "ed_E0" in r)
    macro("dmrgEzeroErr", sci(worst_E0))
    macro("dmrgEoneErr", sci(worst_E1))
    macro("dmrgGammaErr", sci(worst_g))
    macro("weylErr", sci(max(worst_W, 1e-16)))
    for r in rows:
        if r["N"] not in (8, 12, 16, 24, 32, 48):
            continue
        g = r["gamma"]
        lines.append(f"{r['N']} & {r['J2']:g} & {r['E0']:.6f} & {r['gap']:.4f} & "
                     f"{float(g['1']):.3f} & {float(g['2']):.3f} & {float(g['4']):.3f} & {float(g['8']):.4f} \\\\")
    with open(os.path.join(GEN, "tab_inputs.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------- validation
def validation():
    V = load("validate.json")
    if V is None:
        return
    sims = [v for v in V if "sim" in v]
    n_ok = sum(bool(v["sim"]["bound_ok"]) for v in sims)
    macro("nSim", len(sims))
    macro("nSimOK", n_ok)
    macro("nPlans", len(V))
    ratios = [v["plan"]["infidelity_bound"] / max(v["sim"]["infidelity_trotter"], 1e-300) for v in sims]
    macro("boundRatioMin", f"{min(ratios):.1f}")
    macro("boundRatioMax", sci(max(ratios)))
    macro("boundRatioMed", f"{np.median(ratios):.0f}")
    lr = [v["plan"]["leak"] / max(v["sim"]["infidelity_exact"], 1e-300) for v in sims]
    macro("leakRatioMin", f"{min(lr):.1f}")
    macro("leakRatioMax", sci(max(lr)))
    macro("leakRatioMed", f"{np.median(lr):.1f}")
    tr = [v["sim"]["trotter_dist_bound"] / max(v["sim"]["trotter_state_dist"], 1e-300) for v in sims]
    macro("trotRatioMin", f"{min(tr):.0f}")
    macro("trotRatioMax", sci(max(tr)))
    macro("trotRatioMed", f"{np.median(tr):.0f}")
    pr = [v["sim"]["p_succ_trotter"] / v["plan"]["p_g"] for v in sims]
    macro("psuccRatioMin", f"{min(pr):.2f}")
    macro("psuccRatioMax", f"{max(pr):.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.6))
    off = {"uniform": -0.45, "dyadic": 0.0, "optimized": 0.45}
    for v in sims:
        s, o = v["schedule"], v["order"]
        kw = dict(color=SCHED_COL[s], marker="o" if o == 2 else "s",
                  mfc=SCHED_COL[s] if o == 2 else "white", mew=1.0, ls="none", ms=4)
        x = v["N"] + off[s] + (0.12 if o == 1 else -0.12)
        ax[0].plot(x, v["sim"]["infidelity_trotter"] / v["plan"]["infidelity_bound"], **kw)
        ax[1].plot(v["sim"]["trotter_dist_bound"], max(v["sim"]["trotter_state_dist"], 1e-9), **kw)
    ax[0].axhline(1.0, color=INK2, lw=0.8, ls="--")
    ax[0].text(7.2, 1.3, "bound", color=INK2, fontsize=7)
    ax[0].set_yscale("log"); ax[0].set_ylim(1e-6, 3)
    ax[0].set_xticks([8, 10, 12, 14]); ax[0].set_xlim(7, 15)
    ax[0].set_xlabel("$N$"); ax[0].set_ylabel(r"measured $1-F$ / certified bound")
    lim = np.array([1e-9, 1])
    ax[1].plot(lim, lim, color=INK2, lw=0.8, ls="--")
    ax[1].text(2.2e-2, 6e-3, "bound = measured", rotation=28, color=INK2, fontsize=7)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlim(3e-3, 0.2); ax[1].set_ylim(1e-9, 0.2)
    ax[1].set_xlabel(r"bound on Trotter distance $\epsilon_T/\sqrt{p_g}$")
    ax[1].set_ylabel("measured Trotter distance")
    h = [plt.Line2D([], [], color=SCHED_COL[s], marker="o", ls="none", label=s) for s in SCHED_COL]
    h += [plt.Line2D([], [], color=INK2, marker="o", ls="none", label="2nd order"),
          plt.Line2D([], [], color=INK2, marker="s", mfc="white", ls="none", label="1st order")]
    ax[1].legend(handles=h, loc="lower right", ncol=1)
    ax[0].set_title(r"(a) total infidelity ($\epsilon=10^{-2},10^{-3}$)", loc="left")
    ax[1].set_title("(b) Trotter part", loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_validation.pdf"))
    plt.close(fig)

    # filter shapes for one case
    c = [v for v in V if v["N"] == 12 and v["J2"] == 0.2411 and v["eps"] == 1e-2 and v["order"] == 2]
    if c:
        fig, ax = plt.subplots(figsize=(3.4, 2.4))
        delta = c[0]["plan"]["delta"]
        E = np.linspace(0, 1, 40001)
        for v in sorted(c, key=lambda v: v["schedule"]):
            t, ph = np.array(v["times"]), np.array(v["phases"])
            F = np.prod(np.cos(np.outer(E, t) + ph), axis=1)
            ax.plot(E, np.abs(F) + 1e-12, color=SCHED_COL[v["schedule"]], lw=0.9,
                    label=f"{v['schedule']} ($m={v['plan']['m']}$)")
        ax.axvspan(0, delta, color=GRID, alpha=0.8, lw=0)
        ax.text(delta * 1.15, 3e-6, r"$\Delta$", color=INK2)
        ax.set_yscale("log"); ax.set_ylim(1e-6, 2)
        ax.set_xscale("log"); ax.set_xlim(1e-3, 1)
        ax.set_xlabel(r"scaled energy $(E-s)/W$"); ax.set_ylabel(r"$|F(E)|$")
        ax.legend(loc="lower left", fontsize=6.5)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "fig_filters.pdf"))
        plt.close(fig)

    # schedule table at N=12, J2=0.2411
    lines = []
    for eps in (1e-2, 1e-3):
        for s in ("uniform", "dyadic", "optimized"):
            for o in (1, 2):
                vv = [v for v in V if v["N"] == 12 and v["J2"] == 0.2411 and v["eps"] == eps
                      and v["schedule"] == s and v["order"] == o]
                if not vv:
                    continue
                v = vv[0]
                p = v["plan"]
                meas = (f"${sci(v['sim']['infidelity_trotter'])}$" if "sim" in v else "--")
                lines.append(f"$10^{{{int(round(np.log10(eps)))}}}$ & {s} & {o} & {p['m']} & {p['T_phys']:.1f} & "
                             f"{p['p_g']:.2f} & {num(p['res_steps'])} & ${sci(p['res_filter_cnot'])}$ & "
                             f"${sci(p['res_expected_cnot'])}$ & ${sci(p['res_attempt_rotations'])}$ & "
                             f"${sci(p['res_synth_delta'], 0)}$ & ${sci(p['res_expected_t'])}$ & {meas} \\\\")
        lines.append("\\midrule")
    with open(os.path.join(GEN, "tab_schedules.tex"), "w") as fh:
        fh.write("\n".join(lines[:-1]) + "\n")


# ----------------------------------------------------------------------------- scaling
def scaling():
    S = load("scaling.json")
    if S is None:
        return
    fails = [r for r in S if r["plan"] is None]
    macro("nScaleFail", len(fails))
    macro("nScale", len(S))
    S = [r for r in S if r["plan"] is not None]
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.5), sharey=False)
    eps, chi = 1e-2, 2
    fits = {}
    for k, J2 in enumerate(J2S):
        a = ax[k]
        for s in ("optimized", "dyadic"):
            for o in (2, 1):
                rows = sorted([r for r in S if r["J2"] == J2 and r["eps"] == eps and r["chi"] == chi
                               and r["schedule"] == s and r["order"] == o], key=lambda r: r["N"])
                N = [r["N"] for r in rows]
                y = [r["plan"]["res_expected_cnot"] for r in rows]
                a.plot(N, y, color=SCHED_COL[s], ls="-" if o == 2 else "--", marker="o" if o == 2 else "s",
                       mfc=SCHED_COL[s] if o == 2 else "white", mew=1.0,
                       label=f"{s}, order {o}")
                big = [i for i, n in enumerate(N) if n >= 16]
                fits[(J2, s, o, "cnot")] = fit_power([N[i] for i in big], [y[i] for i in big])
                fits[(J2, s, o, "steps")] = fit_power([N[i] for i in big],
                                                      [rows[i]["plan"]["res_steps"] for i in big])
        Q = load("ext_qetu_scaling.json") or []
        rows = sorted([r for r in Q if r["J2"] == J2 and r["eps"] == eps and r["order"] == 2 and r["plan"]],
                      key=lambda r: r["N"])
        if rows:
            a.plot([r["N"] for r in rows], [r["plan"]["res_expected_cnot"] for r in rows], color=QETU_COL,
                   ls=":", marker="^", ms=4, label="QETU, order 2")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xticks([8, 16, 32, 48]); a.set_xticklabels(["8", "16", "32", "48"])
        a.minorticks_off()
        a.set_xlabel("$N$")
        a.set_title(J2_NAME[J2], loc="left")
    ax[0].set_ylabel("expected CNOTs per success")
    ax[0].legend(loc="upper left", fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_scaling.pdf"))
    plt.close(fig)
    for (J2, s, o, q), p in fits.items():
        tag = {0.0: "A", 0.2411: "B", 0.4: "C"}[J2] + {"optimized": "Opt", "dyadic": "Dya"}[s] + \
            ("One" if o == 1 else "Two") + ("Cnot" if q == "cnot" else "Steps")
        macro(f"fit{tag}", f"{p:.1f}")

    # table: the cheaper schedule (second order) per case, with first-order steps of that schedule
    def pick(N, J2, e, o, sched=None):
        c = [r for r in S if r["N"] == N and r["J2"] == J2 and r["eps"] == e and r["chi"] == chi
             and r["order"] == o and (sched is None or r["schedule"] == sched)]
        return min(c, key=lambda r: r["plan"]["res_expected_cnot"]) if c else None

    lines = []
    for N in (16, 24, 32, 48):
        for J2 in J2S:
            b = pick(N, J2, eps, 2)
            if b is None:
                continue
            b1 = pick(N, J2, eps, 1, b["schedule"])
            p1, p2 = (b1 or {}).get("plan"), b["plan"]
            lines.append(f"{N} & {J2:g} & {b['gamma']:.3f} & {p2['delta']:.4f} & {b['schedule'][:3]}. & "
                         f"{p2['T_phys']:.1f} & {p2['p_g']:.2f} & {num(p1['res_steps']) if p1 else '--'} & "
                         f"{num(p2['res_steps'])} & ${sci(p2['res_attempt_cnot'])}$ & "
                         f"${sci(p2['res_expected_cnot'])}$ & ${sci(p2['res_attempt_rotations'])}$ & "
                         f"${sci(p2['res_synth_delta'], 0)}$ & ${sci(p2['res_expected_t'])}$ \\\\")
    with open(os.path.join(GEN, "tab_scaling.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    b = pick(48, 0.2411, 1e-2, 2)
    macro("cnotFortyEight", f"${sci(b['plan']['res_expected_cnot'])}$")
    macro("tFortyEight", f"${sci(b['plan']['res_expected_t'])}$")
    macro("schedFortyEight", b["schedule"])
    # order-2 vs order-1 (same schedule) and optimized vs dyadic (second order) gains
    g12, gdo = [], []
    for r in S:
        if r["order"] != 2 or r["plan"]["m"] == 0:
            continue
        m1 = pick(r["N"], r["J2"], r["eps"], 1, r["schedule"])
        if m1:
            g12.append(m1["plan"]["res_expected_cnot"] / r["plan"]["res_expected_cnot"])
        if r["schedule"] == "optimized":
            md = pick(r["N"], r["J2"], r["eps"], 2, "dyadic")
            if md:
                gdo.append(md["plan"]["res_expected_cnot"] / r["plan"]["res_expected_cnot"])
    macro("gainOrderMin", f"{min(g12):.0f}"); macro("gainOrderMax", f"{max(g12):.0f}")
    macro("gainOrderMed", f"{np.median(g12):.0f}")
    macro("gainOptMin", f"{min(gdo):.1f}"); macro("gainOptMax", f"{max(gdo):.1f}")
    macro("gainOptMed", f"{np.median(gdo):.1f}")
    macro("nDyaWins", sum(g < 1 for g in gdo))
    macro("nOptCmp", len(gdo))


# ----------------------------------------------------------------------------- QETU baseline
def qetu():
    Q, S = load("ext_qetu_scaling.json"), load("scaling.json")
    if Q is None or S is None:
        return
    best = {}
    for r in S:
        if r["plan"] and r["order"] == 2:
            best[(r["N"], r["J2"], r["eps"], r["schedule"])] = r["plan"]["res_expected_cnot"]
    rd, ro, rb, tr = [], [], [], []
    for r in Q:
        if r["order"] != 2 or not r["plan"]:
            continue
        k = (r["N"], r["J2"], r["eps"])
        c = r["plan"]["res_expected_cnot"]
        d, o = best.get(k + ("dyadic",)), best.get(k + ("optimized",))
        rd.append(d / c)
        if o:
            ro.append(o / c)
        rb.append(min(x for x in (d, o) if x) / c)
    macro("qetuFail", sum(1 for r in Q if not r["plan"]))
    macro("qetuN", sum(1 for r in Q if r["order"] == 2))
    macro("qetuDyaMed", f"{np.median(rd):.1f}")
    macro("qetuDyaMin", f"{min(rd):.1f}"); macro("qetuDyaMax", f"{max(rd):.1f}")
    macro("qetuOptMed", f"{np.median(ro):.2f}")
    macro("qetuOptMin", f"{min(ro):.2f}"); macro("qetuOptMax", f"{max(ro):.1f}")
    macro("qetuOptN", len(ro))
    macro("qetuBestMed", f"{np.median(rb):.2f}")
    q48 = [r for r in Q if r["N"] == 48 and r["J2"] == 0.2411 and r["eps"] == 1e-2 and r["order"] == 2][0]
    macro("qetuFortyEight", f"${sci(q48['plan']['res_expected_cnot'])}$")
    V = load("ext_qetu_validate_chain.json") or []
    sims = [v for v in V if "sim" in v]
    macro("qetuNSim", len(sims))
    macro("qetuNSimOK", sum(bool(v["sim"]["bound_ok"]) for v in sims))
    if sims:
        macro("qetuBoundMed", f"{np.median([v['plan']['infidelity_bound'] / v['sim']['infidelity_trotter'] for v in sims]):.1f}")
        macro("qetuLeakMed", f"{np.median([v['plan']['leak'] / v['sim']['infidelity_exact'] for v in sims]):.1f}")


# ----------------------------------------------------------------------------- what drives the cost
def drivers():
    D = load("ext_drivers.json")
    if D is None:
        return
    from collections import defaultdict
    ov = [r for r in D if r["study"] == "overlap" and r["plan"]]
    gp = [r for r in D if r["study"] == "gap" and r["plan"]]
    cost = lambda r: r["plan"]["res_expected_cnot"]
    # exponents: gap (within fixed N, J2, eps, schedule, gamma), overlap (within fixed gap), N (common gap)
    grp = defaultdict(list)
    for r in gp:
        grp[(r["N"], r["J2"], r["eps"], r["schedule"], r["gamma"])].append(r)
    ag = [-np.polyfit(np.log([r["plan"]["delta"] for r in v]), np.log([cost(r) for r in v]), 1)[0]
          for v in grp.values()]
    grp = defaultdict(list)
    for r in ov:
        if r["gapmode"] == "true":
            grp[(r["N"], r["J2"], r["eps"], r["schedule"])].append(r)
    cg, craw = [], []
    for (N, J2, eps, s), v in grp.items():
        g = np.array([r["gamma"] for r in v]); c = np.array([cost(r) for r in v])
        L = np.log((1 - g) / (g * eps))
        cg.append(-np.polyfit(np.log(g), np.log(c / L ** 1.5), 1)[0])
        craw.append(-np.polyfit(np.log(g), np.log(c), 1)[0])
    grp = defaultdict(list)
    for r in ov:
        if r["gapmode"] == "common" and r["N"] >= 16:
            grp[(r["J2"], r["eps"], r["schedule"], r["gamma"])].append(r)
    nN = [np.polyfit(np.log([r["N"] for r in v]), np.log([cost(r) for r in v]), 1)[0] for v in grp.values()]
    macro("drvGapExp", f"{np.median(ag):.2f}")
    macro("drvGapExpMin", f"{min(ag):.2f}"); macro("drvGapExpMax", f"{max(ag):.2f}")
    macro("drvGamExp", f"{np.median(cg):.2f}")
    macro("drvGamExpMin", f"{min(cg):.2f}"); macro("drvGamExpMax", f"{max(cg):.2f}")
    macro("drvGamRaw", f"{np.median(craw):.2f}")
    macro("drvNExp", f"{np.median(nN):.1f}")
    macro("drvNExpMin", f"{min(nN):.1f}"); macro("drvNExpMax", f"{max(nN):.1f}")
    V = load("ext_drivers_validate.json")
    if V:
        ok = [r for r in V if r["sim"] and r["sim"]["bound_ok"]]
        ratio = [r["sim"]["infidelity_trotter"] / r["plan"]["infidelity_bound"] for r in V if r["sim"]]
        macro("drvValN", str(len(V))); macro("drvValOk", str(len(ok)))
        macro("drvValMed", f"{np.median(ratio):.2f}"); macro("drvValMax", f"{max(ratio):.2f}")
        macro("drvValNmax", str(max(r["N"] for r in V)))
    macro("drvNPlans", len(D))
    sel = lambda r, N: r["N"] == N and r["J2"] == 0.2411 and r["eps"] == 1e-2 and r["schedule"] == "qetu"
    lo = [r for r in ov if sel(r, 48) and r["gapmode"] == "true"]
    macro("drvGamSpread", sci(max(map(cost, lo)) / min(map(cost, lo)), 0))
    # figure
    fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.7))
    cols = {8: "#9aa0a6", 16: "#1baf7a", 24: "#eb6834", 32: "#2a78d6", 48: "#0b0b0b"}
    for N, col in cols.items():
        v = sorted([r for r in ov if sel(r, N) and r["gapmode"] == "common"], key=lambda r: r["gamma"])
        ax[0].plot([r["gamma"] for r in v], [cost(r) for r in v], "o-", color=col, ms=3.5, label=f"$N={N}$")
    g = np.logspace(-2, np.log10(0.97), 100)
    v48 = sorted([r for r in ov if sel(r, 48) and r["gapmode"] == "common"], key=lambda r: r["gamma"])
    m = lambda g: g ** (-1.25) * np.log((1 - g) / (g * 1e-2)) ** 1.5
    ax[0].plot(g, m(g) * cost(v48[3]) / m(v48[3]["gamma"]), color=INK2, ls="--", lw=0.8)
    ax[0].text(0.011, 3.5e5, r"$\gamma^{-5/4}\ln^{3/2}\frac{1-\gamma}{\gamma\epsilon}$", color=INK2, fontsize=7)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel(r"trial overlap $\gamma$"); ax[0].set_ylabel("filter CNOTs per success")
    ax[0].set_title("(a) overlap, common gap", loc="left"); ax[0].legend(fontsize=6.5, loc="upper right")
    mk = {0.0: "o", 0.2411: "s", 0.4: "^", 0.5: "D"}
    for N, col in ((16, cols[16]), (32, cols[32]), (48, cols[48])):
        for J2 in (0.0, 0.2411, 0.4, 0.5):
            v = sorted([r for r in gp if r["N"] == N and r["J2"] == J2 and r["eps"] == 1e-2 and r["gamma"] == 0.5
                        and r["schedule"] == "qetu"], key=lambda r: r["f"])
            ax[1].plot([r["plan"]["delta"] for r in v], [cost(r) for r in v], marker=mk[J2], color=col,
                       ms=3.5, lw=0.8, mfc="white" if J2 in (0.4, 0.5) else col)
    x = np.array([5e-4, 3e-2])
    ax[1].plot(x, 2e5 * (x / 1e-2) ** -1.5, color=INK2, ls="--", lw=0.8)
    ax[1].text(1.2e-2, 1.6e5, r"$\Delta^{-3/2}$", color=INK2, fontsize=7)
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel(r"scaled gap $\Delta$"); ax[1].set_title(r"(b) gap, $\gamma=0.5$", loc="left")
    h = [plt.Line2D([], [], color=INK2, marker=mk[j], ls="none", mfc="white" if j in (0.4, 0.5) else INK2,
                    label=f"$J_2={j:g}$") for j in mk]
    h += [plt.Line2D([], [], color=cols[n], lw=1.5, label=f"$N={n}$") for n in (16, 32, 48)]
    ax[1].legend(handles=h, fontsize=6, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_drivers.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------- what the filter actually needs
def empirical():
    T = load("emp_trotter.json")
    if T is None:
        return
    from collections import defaultdict
    T = [r for r in T if r.get("d_min")]
    syn = [r for r in T if str(r["trial"]).startswith("g")]
    act, cer = (lambda r: r["cnot_per_success"]), (lambda r: r["cert_cnot_per_success"])
    macro("empN", len(T))
    macro("empKone", sum(r["k_min"] == 1 for r in T))
    macro("empDegRatio", f"{np.median([r['d_cert'] / r['d_min'] for r in T]):.2f}")
    macro("empDegRatioMax", f"{max(r['d_cert'] / r['d_min'] for r in T):.1f}")
    rat = [cer(r) / act(r) for r in T]
    macro("empCertRatio", f"{np.median(rat):.0f}")
    macro("empCertRatioMin", f"{min(rat):.0f}"); macro("empCertRatioMax", sci(max(rat), 0))
    macro("empCertK", f"{np.median([r['cert_k_per_query'] for r in T]):.0f}")
    # overlap exponent at fixed gap (matched family), with and without the logarithm
    g = defaultdict(list)
    for r in syn:
        if r["family"] == "matched":
            g[(r["N"], r["eps"])].append(r)
    ex, exl, exc = [], [], []
    for (N, eps), v in g.items():
        G = np.array([r["gamma"] for r in v])
        L = np.log((1 - G) / (G * eps))
        ex.append(-np.polyfit(np.log(G), np.log([act(r) for r in v]), 1)[0])
        exl.append(-np.polyfit(np.log(G), np.log([act(r) for r in v] / L), 1)[0])
        exc.append(-np.polyfit(np.log(G), np.log([cer(r) for r in v]), 1)[0])
    macro("empGamRaw", f"{np.median(ex):.2f}"); macro("empGamLog", f"{np.median(exl):.2f}")
    macro("empGamCert", f"{np.median(exc):.2f}")
    # gap exponent at fixed N and overlap (gap scan)
    g = defaultdict(list)
    for r in syn:
        if r["family"] == "gapscan":
            g[(r["N"], r["trial"], r["eps"])].append(r)
    eg, ed, ec = [], [], []
    for v in g.values():
        if len(v) < 3:
            continue
        x = np.log([r["delta"] for r in v])
        eg.append(-np.polyfit(x, np.log([act(r) for r in v]), 1)[0])
        ed.append(-np.polyfit(x, np.log([r["d_min"] for r in v]), 1)[0])
        ec.append(-np.polyfit(x, np.log([cer(r) for r in v]), 1)[0])
    macro("empGapExp", f"{np.median(eg):.2f}"); macro("empGapDeg", f"{np.median(ed):.2f}")
    macro("empGapCert", f"{np.median(ec):.2f}")
    # N at fixed physical gap (matched) and at fixed J2
    for fam, nm in (("matched", "Match"), ("fixedJ2", "Fixed")):
        g = defaultdict(list)
        for r in syn:
            if r["family"] == fam:
                g[(r["trial"], r["eps"])].append(r)
        fit = lambda f: np.median([np.polyfit(np.log([r["N"] for r in v]), np.log([f(r) for r in v]), 1)[0]
                                   for v in g.values()])
        macro(f"empN{nm}", f"{fit(act):.1f}"); macro(f"empN{nm}Cert", f"{fit(cer):.1f}")
        macro(f"empN{nm}Deg", f"{fit(lambda r: r['d_min']):.1f}")
    M = load("emp_match.json")
    macro("empMatchJ", ", ".join(f"{float(M[str(N)]):.3f}" for N in (8, 10, 12, 14, 16)))
    # figure: (a) overlap at fixed gap, (b) gap at fixed overlap, (c) N at fixed gap vs fixed J2
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.5))
    cols = {8: "#9aa0a6", 10: "#1baf7a", 12: "#eb6834", 14: "#2a78d6", 16: "#0b0b0b"}
    for N in (8, 12, 16):
        v = sorted([r for r in syn if r["family"] == "matched" and r["N"] == N and r["eps"] == 1e-3],
                   key=lambda r: r["gamma"])
        ax[0].plot([r["gamma"] for r in v], [act(r) for r in v], "o-", color=cols[N], ms=3.5, label=f"$N={N}$")
        ax[0].plot([r["gamma"] for r in v], [cer(r) for r in v], "o--", color=cols[N], ms=3, mfc="white", lw=0.9)
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel(r"trial overlap $\gamma$"); ax[0].set_ylabel("filter CNOTs per success")
    ax[0].set_title(r"(a) overlap, gap $=0.40$", loc="left"); ax[0].legend(fontsize=6.5, loc="lower left")
    from sbpipe.design import eta_target
    fam_col = {"matched": QETU_COL, "fixedJ2": SCHED_COL["dyadic"], "gapscan": SCHED_COL["uniform"]}
    fam_lab = {"matched": "gap 0.40", "fixedJ2": r"$J_2=0.2411$", "gapscan": "gap scan"}
    xs = []
    for fam in ("fixedJ2", "matched", "gapscan"):
        v = [r for r in T if r["family"] == fam and r["gamma"] < 1 - r["eps"]]
        x = [np.arccosh(1 / eta_target(r["eps"], r["gamma"])) for r in v]
        xs += x
        ax[1].plot(x, [r["d_min"] * r["delta"] for r in v], "o", color=fam_col[fam], ms=2.5, alpha=0.7,
                   mew=0, label=fam_lab[fam])
    x = np.array([min(xs), max(xs)])
    ax[1].plot(x, 2 / np.pi * x, color=INK2, ls="--", lw=0.8)
    ax[1].text(0.6 * x[1], 2 / np.pi * 0.6 * x[1] - 1.2, r"$d\,\Delta=\frac{2}{\pi}\cosh^{-1}\frac{1}{\eta}$",
               color=INK2, fontsize=7)
    ax[1].set_xlabel(r"$\cosh^{-1}(1/\eta)$,  $\eta^2=\gamma\epsilon/(1-\gamma)$")
    ax[1].set_ylabel(r"degree $\times$ scaled gap, $d_{\min}\Delta$")
    ax[1].set_title("(b) degree, all cases", loc="left"); ax[1].legend(fontsize=6.5, loc="upper left")
    for fam, col, lab in (("matched", QETU_COL, "gap 0.40"), ("fixedJ2", SCHED_COL["dyadic"], r"$J_2=0.2411$")):
        v = sorted([r for r in syn if r["family"] == fam and r["eps"] == 1e-3 and r["trial"] == "g0.3"],
                   key=lambda r: r["N"])
        ax[2].plot([r["N"] for r in v], [act(r) for r in v], "o-", color=col, ms=3.5, label=lab)
        ax[2].plot([r["N"] for r in v], [cer(r) for r in v], "o--", color=col, ms=3, mfc="white", lw=0.9)
    ax[2].set_xscale("log"); ax[2].set_yscale("log")
    ax[2].set_xticks([8, 10, 12, 14, 16]); ax[2].set_xticklabels(["8", "10", "12", "14", "16"])
    ax[2].minorticks_off()
    ax[2].set_xlabel(r"$N$"); ax[2].set_title(r"(c) size, $\gamma=0.3$", loc="left")
    ax[2].legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_empirical.pdf"))
    fig.savefig(os.path.join(FIG, "fig_empirical.png"), dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------------------- regime: MPS alone vs filter
def regime():
    R = []
    for f in ("ext_regime_chain.json", "ext_regime_cyl.json", "ext_regime_cyl_12.json"):
        R += load(f) or []
    if not R:
        return
    from collections import defaultdict
    g = defaultdict(list)
    for r in R:
        g[(r["kind"], r["N"], r["J2"], r["eps"])].append(r)
    lines, n_direct, n_none = [], 0, 0
    inputs = {}
    for f in sorted(os.listdir(os.path.join(RES, "cache2d"))):
        if f.endswith(".json"):
            d = json.load(open(os.path.join(RES, "cache2d", f)))
            inputs[(d["N"], d["J2"])] = d
    for kind in ("chain", "cyl"):
        for N in sorted({k[1] for k in g if k[0] == kind}):
            for J2 in sorted({k[2] for k in g if k[0] == kind and k[1] == N}):
                cells = []
                for eps in (1e-2, 1e-3, 1e-4):
                    v = g[(kind, N, J2, eps)]
                    dr = [r for r in v if r["schedule"] == "auto" and 1 - r["gamma"] <= eps]
                    fl = [r for r in v if r["plan"] and r["plan"]["m"] > 0]
                    db = min(dr, key=lambda r: r["trial_cost"]["cnot"]) if dr else None
                    fb = min(fl, key=lambda r: r["plan"]["res_expected_cnot"]) if fl else None
                    n_direct += db is not None
                    n_none += db is None
                    cells.append(f"{db['chi']}" if db else "--")
                    cells.append(f"${sci(db['trial_cost']['cnot'])}$" if db else "--")
                    cells.append(f"${sci(fb['plan']['res_expected_cnot'])}$ ({fb['chi']})" if fb else "--")
                name = f"chain {N}" if kind == "chain" else f"$4\\times{N // 4}$"
                if kind == "cyl":
                    e = inputs[(N, J2)]["est"]
                else:
                    e = json.load(open(os.path.join(RES, "cache", f"N{N}_J2{J2:g}.json")))["est"]
                gap = f"{e['E1'] - e['E0']:.3f}"
                lines.append(f"{name} & {J2:g} & {gap} & " + " & ".join(cells) + " \\\\")
    with open(os.path.join(GEN, "tab_regime.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    macro("regDirect", n_direct); macro("regNone", n_none); macro("regTotal", n_direct + n_none)
    c = inputs[(48, 0.0)]["ansatz"]
    macro("cylGamSixteen", f"{c['16']['gamma']:.2f}"); macro("cylGamSixtyfour", f"{c['64']['gamma']:.4f}")


# ----------------------------------------------------------------------------- infidelity on the (N, CNOT) plane
def fid_plane():
    R = []
    for f in ("ext_regime_chain.json", "ext_regime_cyl.json", "ext_regime_cyl_12.json"):
        R += load(f) or []
    if not R:
        return
    blues = {2: "#9ec2f0", 4: "#6aa0e6", 8: "#2a78d6", 16: "#1d5fb3", 32: "#174a8c", 64: "#0b1f3d"}
    oranges = {1e-2: "#f2a27e", 1e-3: "#eb6834", 1e-4: "#a3401a"}
    fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex="row")
    for row, (kind, J2, name) in enumerate((("chain", 0.2411, r"chain, $J_2=0.2411$"),
                                             ("cyl", 0.0, r"$4\times L_x$ cylinder, $J_2=0$"))):
        v = [r for r in R if r["kind"] == kind and r["J2"] == J2]
        Ns = sorted({r["N"] for r in v})
        a0, a1 = ax[row]
        for chi in sorted({r["chi"] for r in v}):
            pts = sorted({(r["N"], r["gamma"], r["trial_cost"]["cnot"]) for r in v if r["chi"] == chi})
            N = [p[0] for p in pts]
            inf = np.array([max(1 - p[1], 1e-6) for p in pts])
            a0.plot(N, inf, "-", color=blues[chi], label=f"MPS $\\chi={chi}$", marker="o", ms=4)
            low = inf <= 1e-6                        # at or below the precision of the reference state
            a0.plot(np.array(N)[low], inf[low], "o", ms=5, mfc="white", mec=blues[chi], mew=1.2)
            a1.plot(N, [p[2] / p[1] for p in pts], "o-", color=blues[chi], ms=4)
        for eps in (1e-2, 1e-3, 1e-4):
            pts = []
            for n in Ns:
                fl = [r for r in v if r["N"] == n and r["eps"] == eps and r["plan"] and r["plan"]["m"] > 0]
                if fl:
                    pts.append((n, min(r["plan"]["res_expected_cnot"] for r in fl)))
            N = [p[0] for p in pts]
            lab = f"filter $\\epsilon=10^{{{int(np.log10(eps))}}}$"
            a0.plot(N, [eps] * len(N), "D--", color=oranges[eps], ms=4, lw=1.2, label=lab)
            a1.plot(N, [p[1] / (1 - eps) for p in pts], "D--", color=oranges[eps], ms=4, lw=1.2)
        for a in (a0, a1):
            a.set_yscale("log"); a.set_xticks(Ns); a.set_xlabel(r"$N$")
        a0.set_ylabel(r"infidelity $1-F$")
        a1.set_ylabel(r"CNOTs / $F$")
        a0.set_title(f"({'ac'[row]}) {name}: infidelity", loc="left")
        a1.set_title(f"({'bd'[row]}) {name}: CNOTs per unit fidelity", loc="left")
        a0.set_ylim(5e-7, 2)
        a0.axhspan(5e-7, 1.5e-6, color=GRID, alpha=0.6, lw=0, zorder=0)
    h0, l0 = ax[1][0].get_legend_handles_labels()
    fig.legend(h0, l0, loc="lower center", ncol=5, fontsize=6.8, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_fidplane.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, "fig_fidplane.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

# ----------------------------------------------------------------------------- cheapest route to a given accuracy
def _frontiers(F, R, kind, N, J2):
    """MPS alone, empirical filter and certified filter: cheapest CNOTs to reach infidelity <= eps, on a
    fine grid of eps (descending). NaN where no route with chi <= 64 reaches eps."""
    v = [r for r in F if r["kind"] == kind and r["N"] == N and r["J2"] == J2]
    fgrid = np.array(sorted({p[0] for r in v for p in r["points"]}, reverse=True))
    f_raw = np.array([min([p[2] for r in v for p in r["points"] if p[0] == e and p[1] > 0], default=np.nan)
                      for e in fgrid])
    f_cost = np.array([np.nanmin(f_raw[i:]) if np.isfinite(f_raw[i:]).any() else np.nan
                       for i in range(len(fgrid))])                    # a tighter plan also meets a looser target
    grid = np.logspace(-1, -7, 601)
    f = np.exp(np.interp(np.log(grid[::-1]), np.log(fgrid[::-1]), np.log(f_cost[::-1])))[::-1]
    mps = sorted([(max(1 - r["gamma"], 0.0), r["trial_cnot"], r["chi"]) for r in v])
    m = np.array([min([c for inf, c, _ in mps if inf <= e], default=np.nan) for e in grid])
    c = [r for r in R if r["kind"] == kind and r["N"] == N and r["J2"] == J2 and r["plan"] and r["plan"]["m"] > 0]
    cert = sorted({(e, min(r["plan"]["res_expected_cnot"] for r in c if r["eps"] == e)) for e in {r["eps"] for r in c}},
                  reverse=True)
    return grid, m, f, cert, mps


def crossover():
    F = load("emp_frontier.json")
    if F is None:
        return
    R = []
    for f in ("ext_regime_chain.json", "ext_regime_cyl.json", "ext_regime_cyl_12.json"):
        R += load(f) or []
    systems = [("chain", 16, 0.2411, "chain $N=16$"), ("chain", 32, 0.2411, "chain $N=32$"),
               ("chain", 48, 0.2411, "chain $N=48$"), ("cyl", 16, 0.0, "cylinder $4\\times4$"),
               ("cyl", 24, 0.0, "cylinder $4\\times6$"), ("cyl", 32, 0.0, "cylinder $4\\times8$"),
               ("cyl", 48, 0.0, "cylinder $4\\times12$")]
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.2), sharey=True)
    axs = axs.ravel()
    MPS_C, EMP_C, CERT_C = SCHED_COL["optimized"], QETU_COL, SCHED_COL["dyadic"]
    rows = []
    for a, (kind, N, J2, name) in zip(axs, systems):
        grid, m, f, cert, mps = _frontiers(F, R, kind, N, J2)
        a.plot(grid, m, color=MPS_C, lw=1.6, drawstyle="steps-pre")
        for inf, c, chi in mps:
            if 1e-7 <= inf <= 1e-1:
                a.plot(inf, c, "o", color=MPS_C, ms=3.5)
                a.annotate(f"{chi}", (inf, c), xytext=(-2, 4), textcoords="offset points", fontsize=5.5,
                           color=INK2, ha="right")
        a.plot(grid, f, "-", color=EMP_C, lw=1.6)
        a.plot([e for e, _ in cert], [c for _, c in cert], "D--", color=CERT_C, ms=3.5, lw=0.9)
        both = np.isfinite(m) & np.isfinite(f)
        cheaper = both & (f < m)
        none = ~np.isfinite(m)
        a.fill_between(grid, 1e1, 1e8, where=cheaper, color=EMP_C, alpha=0.18, lw=0, step="pre")
        a.fill_between(grid, 1e1, 1e8, where=none, step="pre", facecolor=GRID, alpha=0.7, lw=0, hatch="////", edgecolor="white")
        x = grid[np.argmax(cheaper)] if cheaper.any() else None
        rows.append((name, x, float(grid[cheaper].min()) if cheaper.any() else None,
                     float(grid[none].max()) if none.any() else None))
        a.set_xscale("log"); a.set_yscale("log"); a.invert_xaxis()
        a.set_xlim(1.5e-1, 7e-8); a.set_ylim(1e2, 3e7)
        a.set_title(name, loc="left", fontsize=7.5)
        a.set_xticks([1e-1, 1e-3, 1e-5, 1e-7])
    for a in axs[[0, 4]]:
        a.set_ylabel("CNOTs to reach $1-F$")
    for a in axs[3:7]:
        a.set_xlabel(r"infidelity $1-F$")
    axs[7].axis("off")
    h = [plt.Line2D([], [], color=MPS_C, marker="o", ms=3.5, label="MPS alone (number: $\\chi$)"),
         plt.Line2D([], [], color=EMP_C, lw=1.6, label="MPS + filter, empirical"),
         plt.Line2D([], [], color=CERT_C, marker="D", ls="--", ms=3.5, lw=0.9, label="MPS + filter, certified"),
         plt.Rectangle((0, 0), 1, 1, color=EMP_C, alpha=0.18, lw=0, label="empirical filter cheaper"),
         plt.Rectangle((0, 0), 1, 1, facecolor=GRID, alpha=0.7, hatch="////", edgecolor="white",
                       label="no MPS with $\\chi\\leq64$ reaches it")]
    axs[7].legend(handles=h, loc="center", fontsize=6.8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_crossover.pdf"))
    fig.savefig(os.path.join(FIG, "fig_crossover.png"), dpi=200)
    plt.close(fig)
    for name, x, lo, none in rows:
        print("crossover", name, "first", x, "last", lo, "mps limit", none)


# ----------------------------------------------------------------------------- wider sensitivity
def sens_wide():
    C = load("ext_sensitivity.json")
    if C is None:
        return
    scheds = ("dyadic", "optimized", "qetu")
    lines = []
    for kind, fs in (("gap", (1.0, 1.1, 1.25, 1.5, 2.0, 3.0)), ("gamma", (1.02, 1.05))):
        for f in fs:
            cells = []
            for s in scheds:
                rs = [r for r in C if r["kind"] == kind and r["f"] == f and r["schedule"] == s and r.get("sim")]
                miss = [r for r in rs if r["sim"]["infidelity_trotter"] > r["plan"]["eps"]]
                worst = max(r["sim"]["infidelity_trotter"] for r in rs)
                cells.append(f"{len(miss)}/{len(rs)} & ${sci(worst)}$")
            name = "gap" if kind == "gap" else r"$\gamma$"
            lines.append(f"{name} & {f:g} & " + " & ".join(cells) + " \\\\")
    with open(os.path.join(GEN, "tab_sens_wide.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    tot = lambda kind, f: [r for r in C if r["kind"] == kind and r["f"] == f and r.get("sim")]
    miss = lambda rs: sum(r["sim"]["infidelity_trotter"] > r["plan"]["eps"] for r in rs)
    for f, tag in ((1.0, "One"), (1.1, "OneTen"), (1.25, "OneQuarter"), (1.5, "OneHalf")):
        macro(f"sensGap{tag}Miss", miss(tot("gap", f)))
    macro("sensGapN", len(tot("gap", 1.0)))
    macro("sensGapOneTenWorst", f"${sci(max(r['sim']['infidelity_trotter'] for r in tot('gap', 1.1)))}$")
    dy = [r for r in tot("gap", 1.25) if r["schedule"] == "dyadic"]
    macro("sensDyaQuarterMiss", miss(dy))
    g = [r for r in tot("gamma", 1.02) + tot("gamma", 1.05) if r["sim"]["infidelity_trotter"] > r["plan"]["eps"]]
    macro("sensGammaMiss", len(g))
    macro("sensGammaN", len(tot("gamma", 1.02) + tot("gamma", 1.05)))
    macro("sensGammaNoFilter", sum(r["plan"]["m"] == 0 for r in g))


# ----------------------------------------------------------------------------- sensitivity
def sensitivity():
    C = load("sensitivity.json")
    if C is None:
        return
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4))
    eps = C[0]["eps"]
    for k, (kind, xl) in enumerate([("gap", r"assumed gap / true gap"),
                                    ("E0err", r"$\delta E_0$ / gap"),
                                    ("gamma", r"assumed $\gamma$ / true $\gamma$")]):
        rows = sorted([c for c in C if c["kind"] == kind], key=lambda c: c["f"])
        f = [c["f"] for c in rows]
        meas = [c["sim"]["infidelity_trotter"] for c in rows]
        cost = [c["plan"]["res_expected_cnot"] for c in rows]
        ref = [c for c in rows if (c["f"] == 1.0 if kind != "E0err" else c["f"] == min(f))][0]
        a = ax[k]
        a.plot(f, meas, color=SCHED_COL["optimized"], marker="o", label=r"measured $1-F$")
        a.plot(f, [c / ref["plan"]["res_expected_cnot"] * eps for c in cost], color=SCHED_COL["dyadic"],
               marker="s", mfc="white", mew=1.0, ls="--", label=r"cost / cost$_{\rm ref}$ $\times\,\epsilon$")
        a.axhline(eps, color=INK2, lw=0.8, ls=":")
        a.text(f[0], eps * 1.25, r"target $\epsilon$", color=INK2, fontsize=7)
        a.set_yscale("log")
        if kind == "E0err":
            a.set_xscale("log")
        a.set_xlabel(xl)
        bad = [c for c in rows if c["sim"]["infidelity_trotter"] > eps]
        if kind == "gap":
            macro("gapFailFactor", f"{min(c['f'] for c in bad):g}" if bad else "none")
            ok = [c for c in rows if c["f"] < 1]
            macro("gapMarginCost", f"{max(c['plan']['res_expected_cnot'] for c in ok) / ref['plan']['res_expected_cnot']:.1f}")
            macro("gapMarginF", f"{min(c['f'] for c in ok):g}")
            worst = max(rows, key=lambda c: c["sim"]["infidelity_trotter"])
            macro("gapWorstInf", f"${sci(worst['sim']['infidelity_trotter'])}$")
            macro("gapWorstF", f"{worst['f']:g}")
            macro("nGapFail", len(bad))
        if kind == "gamma":
            macro("nGammaFail", len(bad))
            if bad:
                macro("gammaFailF", f"{min(c['f'] for c in bad):g}")
            low = min(rows, key=lambda c: c["f"])
            macro("gammaLowCost", f"{low['plan']['res_expected_cnot'] / ref['plan']['res_expected_cnot']:.1f}")
        if kind == "E0err":
            hi = max(rows, key=lambda c: c["f"])
            macro("eZeroHiCost", f"{hi['plan']['res_expected_cnot'] / ref['plan']['res_expected_cnot']:.1f}")
            macro("eZeroHiF", f"{hi['f']:g}")
            macro("nEzeroFail", len(bad))
    ax[0].set_ylabel(r"infidelity, scaled cost")
    ax[0].legend(loc="lower right", fontsize=6.5)
    for a, t in zip(ax, ["(a) gap", r"(b) $E_0$ error bar", r"(c) overlap"]):
        a.set_title(t, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_sensitivity.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------- ansatz
def ansatz():
    A = load("ansatz.json")
    if A is None:
        return
    macro("nAnsatzFail", sum(r["plan"] is None for r in A))
    A = [r for r in A if r["plan"] is not None]
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.4), sharey=True)
    eps = 1e-2
    Ns = sorted({a["N"] for a in A})
    shades = ["#9ec5f0", "#5f9fe5", "#2a78d6", "#1a4f93"][-len(Ns):]   # one hue, light -> dark = N
    best = {}
    for k, J2 in enumerate(J2S):
        a = ax[k]
        for N, col in zip(Ns, shades):
            rows = sorted([r for r in A if r["N"] == N and r["J2"] == J2 and r["eps"] == eps],
                          key=lambda r: r["chi"])
            chi = [r["chi"] for r in rows]
            tot = [r["plan"]["res_expected_cnot"] for r in rows]
            a.plot(chi, tot, color=col, marker="o", label=f"$N={N}$")
            best[(N, J2)] = rows[int(np.argmin(tot))]["chi"]
            # direct preparation: the trial alone suffices if 1-gamma <= eps
            for r in rows:
                if 1 - r["gamma"] <= eps:
                    a.plot(r["chi"], r["trial_cost"]["cnot"], marker="*", color=col, ms=8, ls="none")
        a.set_xscale("log", base=2); a.set_yscale("log")
        a.set_xticks([1, 2, 4, 8]); a.set_xticklabels(["1", "2", "4", "8"]); a.minorticks_off()
        a.set_xlabel(r"trial bond dimension $\chi$")
        a.set_title(J2_NAME[J2], loc="left")
    ax[0].set_ylabel("expected CNOTs per success")
    ax[0].legend(loc="upper right", fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_ansatz.pdf"))
    plt.close(fig)
    lines = []
    for N in Ns:
        for J2 in J2S:
            rows = sorted([r for r in A if r["N"] == N and r["J2"] == J2 and r["eps"] == eps],
                          key=lambda r: r["chi"])
            by = {r["chi"]: r for r in rows}
            cells = " & ".join((f"{by[c]['gamma']:.3f} & ${sci(by[c]['plan']['res_expected_cnot'])}$"
                                if c in by else "-- & --") for c in (1, 2, 4, 8))
            lines.append(f"{N} & {J2:g} & {cells} \\\\")
    with open(os.path.join(GEN, "tab_ansatz.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    macro("bestChiList", ", ".join(sorted({str(v) for v in best.values()})))
    none = [r for r in A if r["plan"]["m"] == 0]
    macro("nNoFilter", len(none))
    macro("nNoFilterChiEight", sum(r["chi"] == 8 for r in none))
    macro("nChiEight", sum(r["chi"] == 8 for r in A))
    r4 = [r for r in A if r["N"] == 48 and r["J2"] == 0.2411 and r["chi"] == 4 and r["eps"] == eps][0]
    r8 = [r for r in A if r["N"] == 48 and r["J2"] == 0.2411 and r["chi"] == 8 and r["eps"] == eps][0]
    macro("cnotFortyEightChiFour", f"${sci(r4['plan']['res_expected_cnot'])}$")
    macro("tFortyEightChiFour", f"${sci(r4['plan']['res_expected_t'])}$")
    macro("cnotFortyEightChiEight", num(r8['plan']['res_expected_cnot']))
    macro("gammaFortyEightChiEight", f"{r8['gamma']:.4f}")


# ----------------------------------------------------------------------------- measured-only crossover
def _actual_rows(A, sys_):
    R = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
    eps = np.array(sorted({r["eps"] for r in R}, reverse=True))
    out = dict(eps=eps)
    for m in ("mps", "dyadic", "optimized", "qetu"):
        c, chi = [], []
        for e in eps:
            Q = [r for r in R if r["eps"] == e]
            if m == "mps":
                q = [(r["trial_cnot"], r["chi"]) for r in Q if r["mps_alone"]]
            else:
                q = [(r[m]["cnot"], r["chi"]) for r in Q if not r["mps_alone"] and r.get(m)]
            b = min(q) if q else (np.nan, None)
            c.append(b[0]); chi.append(b[1])
        out[m] = (np.array(c, float), chi)
    return out


def actual_crossover():
    A = load("emp_actual.json")
    if A is None:
        return
    systems = [("chain", 16, 0.2411, "chain $N=16$, $J_2=0.2411$"), ("cyl", 3, 0.0, "cylinder $4\\times3$, $J_2=0$"),
               ("cyl", 3, 0.3, "cylinder $4\\times3$, $J_2=0.3$"), ("cyl", 4, 0.0, "cylinder $4\\times4$, $J_2=0$"),
               ("cyl", 4, 0.3, "cylinder $4\\times4$, $J_2=0.3$")]
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey=True)
    axs = axs.ravel()
    MPS_C = INK2
    wins, ratios = 0, []
    for a, (kind, sz, J2, name) in zip(axs, systems):
        o = _actual_rows(A, (kind, sz, J2))
        e = o["eps"]
        m, mchi = o["mps"]
        sbc = np.fmin(o["dyadic"][0], o["optimized"][0])
        a.plot(e, m, "o-", color=MPS_C, ms=3.5, lw=1.2, drawstyle="steps-mid")
        for x, y, c in zip(e, m, mchi):
            if c:
                a.annotate(f"{c}", (x, y), xytext=(0, 4), textcoords="offset points", fontsize=5.5, color=MPS_C, ha="center")
        a.plot(e, o["dyadic"][0], "o-", color=SCHED_COL["dyadic"], ms=3.5)
        for x, y, c in zip(e, *o["dyadic"]):
            if c:
                a.annotate(f"{c}", (x, y), xytext=(0, -8), textcoords="offset points", fontsize=5.5,
                           color=SCHED_COL["dyadic"], ha="center")
        ok = np.isfinite(o["optimized"][0])
        a.plot(e[ok], o["optimized"][0][ok], "s", color=SCHED_COL["optimized"], ms=3.5)
        a.plot(e, o["qetu"][0], "^--", color=QETU_COL, ms=3, lw=0.9)
        cheaper = np.isfinite(sbc) & (~np.isfinite(m) | (sbc < m))
        none = ~np.isfinite(m)
        wins += int(cheaper.sum())
        ratios += list((m / sbc)[cheaper & ~none])
        lo, hi = np.log10(e), None
        edges = 10 ** np.concatenate([[lo[0] + 0.25], (lo[:-1] + lo[1:]) / 2, [lo[-1] - 0.25]])
        for i in range(len(e)):
            if cheaper[i]:
                a.axvspan(edges[i], edges[i + 1], color=SCHED_COL["dyadic"], alpha=0.12, lw=0)
            if none[i]:
                a.axvspan(edges[i], edges[i + 1], facecolor=GRID, alpha=0.7, hatch="////", edgecolor="white", lw=0)
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlim(edges[0], edges[-1]); a.set_ylim(1e2, 2e5)
        a.set_title(name, loc="left", fontsize=7.5)
        a.set_xticks([1e-2, 1e-3, 1e-4, 1e-5])
    for a in axs[[0, 3]]:
        a.set_ylabel("CNOTs per success")
    for a in axs[2:5]:
        a.set_xlabel(r"target infidelity $\varepsilon$")
    axs[5].axis("off")
    h = [plt.Line2D([], [], color=MPS_C, marker="o", ms=3.5, lw=1.2, label="MPS alone (number: $\\chi$)"),
         plt.Line2D([], [], color=SCHED_COL["dyadic"], marker="o", ms=3.5, label="MPS + dyadic SBC filter\n(number: trial $\\chi$)"),
         plt.Line2D([], [], color=SCHED_COL["optimized"], marker="s", ls="", ms=3.5, label="MPS + optimized SBC filter"),
         plt.Line2D([], [], color=QETU_COL, marker="^", ls="--", ms=3, lw=0.9, label="MPS + QETU (reference)"),
         plt.Rectangle((0, 0), 1, 1, color=SCHED_COL["dyadic"], alpha=0.12, lw=0, label="SBC filter cheaper than MPS"),
         plt.Rectangle((0, 0), 1, 1, facecolor=GRID, alpha=0.7, hatch="////", edgecolor="white",
                       label="no MPS with $\\chi\\leq64$ reaches $\\varepsilon$")]
    axs[5].legend(handles=h, loc="center", fontsize=6.6)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_actual.pdf"))
    fig.savefig(os.path.join(FIG, "fig_actual.png"), dpi=200)
    plt.close(fig)
    macro("actWins", str(wins)); macro("actRatioMax", f"{max(ratios):.1f}")
    cyl = []
    for kind, sz, J2, _ in systems:
        if kind != "cyl":
            continue
        o = _actual_rows(A, (kind, sz, J2))
        m = o["mps"][0]
        sbc = np.fmin(o["dyadic"][0], o["optimized"][0])
        cyl += [m[i] / sbc[i] for i, e in enumerate(o["eps"]) if e <= 3.1e-3 and np.isfinite(m[i])]
    macro("actRatioLo", f"{min(cyl):.1f}")
    print("actual: SBC cheaper in", wins, "of", 5 * 7, "targets; ratios", np.round(sorted(ratios), 2))


def _actual1d_data():
    p = os.path.join(RES, "emp_actual1d.jsonl")
    if not os.path.exists(p):
        return None
    A = [json.loads(l) for l in open(p)]
    return A + [r for r in (load("emp_actual.json") or []) if r["kind"] == "chain"]


def actual1d():
    """Measured crossover on nine chains (N = 12, 14, 16 x J2 = 0, 0.2411, 0.4): fig_actual1d and act1* macros."""
    A = _actual1d_data()
    if A is None:
        return
    Ns, J2s = (12, 14, 16), (0.0, 0.2411, 0.4)
    fig, axs = plt.subplots(3, 3, figsize=(7.2, 6.4), sharex=True, sharey=True)
    MPS_C = INK2
    wins, tot, ratios = 0, 0, []
    for i, J2 in enumerate(J2s):
        for j, N in enumerate(Ns):
            a = axs[i, j]
            o = _actual_rows(A, ("chain", N, J2))
            e = o["eps"]
            m, mchi = o["mps"]
            sbc = np.fmin(o["dyadic"][0], o["optimized"][0])
            a.plot(e, m, "o-", color=MPS_C, ms=3.2, lw=1.2, drawstyle="steps-mid")
            for x, y, c in zip(e, m, mchi):
                if c:
                    a.annotate(f"{c}", (x, y), xytext=(0, 4), textcoords="offset points", fontsize=5.5, color=MPS_C, ha="center")
            a.plot(e, o["dyadic"][0], "o-", color=SCHED_COL["dyadic"], ms=3.2)
            for x, y, c in zip(e, *o["dyadic"]):
                if c:
                    a.annotate(f"{c}", (x, y), xytext=(0, -8), textcoords="offset points", fontsize=5.5,
                               color=SCHED_COL["dyadic"], ha="center")
            ok = np.isfinite(o["optimized"][0])
            a.plot(e[ok], o["optimized"][0][ok], "s", color=SCHED_COL["optimized"], ms=3.2)
            a.plot(e, o["qetu"][0], "^--", color=QETU_COL, ms=2.8, lw=0.9)
            cheaper = np.isfinite(sbc) & (sbc < m)
            wins += int(cheaper.sum()); tot += len(e)
            ratios += list((m / sbc)[cheaper])
            lo = np.log10(e)
            edges = 10 ** np.concatenate([[lo[0] + 0.25], (lo[:-1] + lo[1:]) / 2, [lo[-1] - 0.25]])
            for k in range(len(e)):
                if cheaper[k]:
                    a.axvspan(edges[k], edges[k + 1], color=SCHED_COL["dyadic"], alpha=0.12, lw=0)
            a.set_xscale("log"); a.set_yscale("log")
            a.set_xlim(edges[0], edges[-1]); a.set_ylim(1e2, 1e5)
            a.set_title(f"$N={N}$, $J_2={J2:g}$", loc="left", fontsize=7.5)
            a.set_xticks([1e-2, 1e-3, 1e-4, 1e-5])
            if j == 0:
                a.set_ylabel("CNOTs per success")
            if i == 2:
                a.set_xlabel(r"target infidelity $\varepsilon$")
    h = [plt.Line2D([], [], color=MPS_C, marker="o", ms=3.2, lw=1.2, label="MPS alone (number: $\\chi$)"),
         plt.Line2D([], [], color=SCHED_COL["dyadic"], marker="o", ms=3.2, label="MPS + dyadic SBC (number: trial $\\chi$)"),
         plt.Line2D([], [], color=SCHED_COL["optimized"], marker="s", ls="", ms=3.2, label="MPS + optimized SBC"),
         plt.Line2D([], [], color=QETU_COL, marker="^", ls="--", ms=2.8, lw=0.9, label="MPS + QETU"),
         plt.Rectangle((0, 0), 1, 1, color=SCHED_COL["dyadic"], alpha=0.12, lw=0, label="SBC cheaper than MPS alone")]
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=6.6, frameon=False)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(FIG, "fig_actual1d.pdf"))
    fig.savefig(os.path.join(FIG, "fig_actual1d.png"), dpi=200)
    plt.close(fig)
    # every filtered route against the cheapest MPS alone at the same (system, eps), by how far the trial is from eps
    pts = []
    for r in A:
        if r["mps_alone"]:
            continue
        m = min(q["trial_cnot"] for q in A if q["mps_alone"] and
                (q["size"], q["J2"], q["eps"]) == (r["size"], r["J2"], r["eps"]))
        c = min(r[k]["cnot"] for k in ("dyadic", "optimized") if r.get(k))
        pts.append(((1 - r["gamma"]) / r["eps"], m / c))
    pts = np.array(pts)
    win = pts[:, 1] > 1
    far = pts[:, 0] > pts[win, 0].max()
    macro("actOneWins", str(wins)); macro("actOneTot", str(tot)); macro("actOneLoss", str(tot - wins))
    macro("actOneRatioMax", f"{max(ratios):.1f}")
    macro("actOneFarMax", f"{pts[win, 0].max():.1f}")
    macro("actOneNearWins", str(int(win.sum()))); macro("actOneNear", str(int((~far).sum())))
    macro("actOneFar", str(int(far.sum()))); macro("actOneNearPlusFar", str(len(far)))
    macro("actOneFarMed", f"{1 / np.median(pts[far, 1]):.0f}")
    print("actual1d: SBC cheaper in", wins, "of", tot, "; ratios", np.round(sorted(ratios), 2),
          "; win (1-g)/eps <=", pts[win, 0].max(), "; far cases", far.sum(), "median MPS advantage", 1 / np.median(pts[far, 1]))


# ----------------------------------------------------------------------------- gap vs overlap
def _gapov():
    """gap/overlap cases with the optimized SBC results attached (derive_gapov_filter.load_cases)."""
    if not os.path.exists(os.path.join(RES, "gapov.jsonl")):
        return None
    from derive_gapov_filter import load_cases
    return load_cases()


GO_NCOL = {8: "#9ec5f0", 10: "#5a9be0", 12: "#2a6bc0", 14: "#0d3a80"}


def _go_cost(r, key):
    """Expected CNOTs per success: measured from simulation, certified from the theorem."""
    m = r[key]
    return m["cnot"] if m else np.nan


def _go_fit(G, key):
    """log-linear fit of the expected CNOTs per success C against N, gap, L = ln[(1-g)/(g eps)] and gamma."""
    X, y = [], []
    for r in G:
        m = r[key]
        if not m:
            continue
        g, e = r["gamma"], r["eps"]
        X.append([1, np.log(r["N"]), np.log(r["gap"]), np.log(np.log((1 - g) / (g * e))), np.log(g)])
        y.append(np.log(m["cnot"]))
    X, y = np.array(X), np.array(y)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    return c, float(np.exp(np.sqrt(np.mean(res ** 2)))), float(np.exp(np.abs(res).max())), X, y


def gapov():
    G = _gapov()
    if G is None:
        return
    Ns = sorted({r["N"] for r in G})
    EPS = (1e-2, 1e-3, 1e-4)
    # --- study A: gap fixed, overlap varied; study B: overlap fixed, gap varied
    for tag, gfix, xkey, xlab, sel in (
            ("A", 0.40, "gamma", r"trial overlap $\gamma$", lambda r, e: abs(r["gap_target"] - 0.40) < 1e-9 and r["eps"] == e),
            ("B", 0.10, "gap", r"gap $\Delta$",
             lambda r, e: r["gamma_target"] == 0.10 and r["eps"] == e)):
        fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=True)
        for a, e in zip(axs, EPS):
            for N in Ns:
                R = sorted([r for r in G if r["N"] == N and sel(r, e)], key=lambda r: r[xkey])
                if not R:
                    continue
                x = [r[xkey] for r in R]
                a.plot(x, [_go_cost(r, "meas_dyadic") for r in R], "o-", color=GO_NCOL[N], ms=3.2, label=f"$N={N}$")
                a.plot(x, [_go_cost(r, "cert_dyadic") for r in R], "D--", color=GO_NCOL[N], ms=2.6, lw=0.9)
            a.set_xscale("log"); a.set_yscale("log")
            a.set_title(rf"$\varepsilon=10^{{{int(np.log10(e))}}}$", loc="left")
            a.set_xlabel(xlab)
            if tag == "B":
                a.set_xticks([0.3, 0.4, 0.5]); a.set_xticklabels(["0.3", "0.4", "0.5"]); a.minorticks_off()
        axs[0].set_ylabel("CNOTs per success")
        h = [plt.Line2D([], [], color=GO_NCOL[N], lw=1.6, label=f"$N={N}$") for N in Ns]
        h += [plt.Line2D([], [], color=INK2, marker="o", ms=3.2, label="measured (simulation)"),
              plt.Line2D([], [], color=INK2, marker="D", ls="--", ms=2.6, lw=0.9, label="certified (Theorem 1)")]
        fig.legend(handles=h, fontsize=6.6, loc="center right")
        sub = (r"gap fixed at $\Delta=0.40$" if tag == "A" else r"overlap fixed at $\gamma=0.10$")
        fig.suptitle(f"Dyadic SBC filter, {sub}", fontsize=8.5, x=0.01, ha="left")
        fig.tight_layout(rect=(0, 0, 0.8, 1))
        fig.savefig(os.path.join(FIG, f"fig_gapov_{tag}.pdf"))
        fig.savefig(os.path.join(FIG, f"fig_gapov_{tag}.png"), dpi=200)
        plt.close(fig)
    # --- fits and collapse
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.4))
    lines, TAB = [], []
    keys = (("meas_dyadic", "dyadic, measured"), ("meas_optimized", "optimized, measured"),
            ("meas_qetu", "QETU, measured"), ("cert_dyadic", "dyadic, certified"),
            ("cert_optimized", "optimized, certified"), ("cert_qetu", "QETU, certified"))
    for a, (key, name) in zip(axs.ravel(), keys):
        c, rms, mx, X, y = _go_fit(G, key)
        col = SCHED_COL["dyadic"] if "dyadic" in key else SCHED_COL["optimized"] if "optim" in key else QETU_COL
        a.plot(np.exp(X @ c), np.exp(y), "o", color=col, ms=2, alpha=0.6)
        lo, hi = np.exp(min(y.min(), (X @ c).min())), np.exp(max(y.max(), (X @ c).max()))
        a.plot([lo, hi], [lo, hi], color=INK2, lw=0.8)
        a.set_xscale("log"); a.set_yscale("log")
        a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        a.set_title(f"{name}\nrms {rms:.2f}$\\times$, max {mx:.1f}$\\times$", loc="left", fontsize=7)
        a.set_xlabel("fit")
        tagm = "".join(w[0] for w in key.split("_")).upper()   # MD, MQ, CD, CQ
        macro(f"go{tagm}a", f"{np.exp(c[0]):.3g}"); macro(f"go{tagm}N", f"{c[1]:.2f}")
        macro(f"go{tagm}gap", f"{-c[2]:.2f}"); macro(f"go{tagm}L", f"{c[3]:.2f}")
        macro(f"go{tagm}rms", f"{rms:.2f}"); macro(f"go{tagm}max", f"{mx:.1f}")
        macro(f"go{tagm}gam", f"{-c[4]:.2f}")
        lines.append((name, c, rms, mx))
        TAB.append(f"{name} & {len(y)} & ${sci(np.exp(c[0]), 1)}$ & {c[1]:.2f} & {-c[2]:.2f} & {c[3]:.2f} & {-c[4]:.2f} & "
                   f"{rms:.2f} & {mx:.1f} \\\\")
    for a in axs[:, 0]:
        a.set_ylabel("CNOTs per success")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_gapov_fit.pdf"))
    fig.savefig(os.path.join(FIG, "fig_gapov_fit.png"), dpi=200)
    plt.close(fig)
    with open(os.path.join(GEN, "tab_gapov_fit.tex"), "w") as fh:
        fh.write("\n".join(TAB) + "\n")
    rat = [r["cert_dyadic"]["cnot"] / r["meas_dyadic"]["cnot"] for r in G if r["cert_dyadic"] and r["meas_dyadic"]]
    macro("goCertRatioMed", f"{np.median(rat):.0f}"); macro("goCertRatioMin", f"{min(rat):.0f}")
    macro("goCertRatioMax", f"{max(rat):.0f}"); macro("goNcases", str(len(G)))
    sv = [r["meas_qetu"]["cnot"] / r["meas_dyadic"]["cnot"] for r in G if r["meas_qetu"] and r["meas_dyadic"]]
    macro("goQetuOverDyMed", f"{np.median(sv):.2f}")
    for name, c, rms, mx in lines:
        print("gapov fit", name, np.round(c, 3), "rms", round(rms, 3), "max", round(mx, 2))
    print("cert/meas dyadic", np.round(np.percentile(rat, [0, 50, 100]), 1), "qetu/dyadic meas", np.round(np.percentile(sv, [0, 50, 100]), 2))


# ----------------------------------------------------------------------------- filter-only (gap vs overlap)
def gapov_filter():
    """Macros for the filter-only analysis (results/gapov_filter.json from derive_gapov_filter.py) and copies
    of the gap/overlap report figures into the paper's figure folder."""
    import shutil
    rep = os.path.join(HERE, "gapov_report", "figs")
    for f in ("fig_filter_A", "fig_filter_B", "fig_filter_mech", "fig_gridmap", "fig_gridmap_shot", "fig_heat_dyadic_succ",
              "fig_heat_ratio", "fig_heat_meas_succ"):
        if os.path.exists(os.path.join(rep, f + ".pdf")):
            shutil.copy(os.path.join(rep, f + ".pdf"), FIG)
    D = load("gapov_filter.json")
    if D is None:
        return
    Ns = sorted({r["N"] for r in D})

    def rng(meth, met, sel, a, b):
        v = []
        for e in (1e-2, 1e-3, 1e-4):
            for N in Ns:
                R = [r for r in D if r["N"] == N and r["eps"] == e and sel(r)]
                A = [r for r in R if a(r)]
                B = [r for r in R if b(r)]
                if A and B:
                    v.append(A[0][meth][met] / B[0][meth][met])
        return min(v), max(v)

    selA = lambda r: abs(r["gap_target"] - 0.4) < 1e-9
    for meth, t in (("meas_dyadic", "MD"), ("meas_qetu", "MQ"), ("cert_dyadic", "CD"), ("cert_qetu", "CQ")):
        for met in ("cnot", "time"):
            lo, hi = rng(meth, met, selA, lambda r: r["gamma_target"] == 0.03, lambda r: r["gamma_target"] == 0.9)
            macro(f"fA{t}{met}", f"{lo:.1f}--{hi:.1f}")
            v = []
            for e in (1e-2, 1e-3, 1e-4):
                for N in (12, 14):
                    R = sorted([r for r in D if r["N"] == N and r["gamma_target"] == 0.1 and r["eps"] == e],
                               key=lambda r: r["gap"])
                    v.append(R[0][meth][met] / R[-1][meth][met])
            macro(f"fB{t}{met}", f"{min(v):.1f}--{max(v):.1f}")
    for fam, t in (("dyadic", "D"), ("qetu", "Q")):
        for met in ("time", "steps"):
            v = [r["cert_" + fam][met] / r["meas_" + fam][met] for r in D]
            macro(f"fCM{t}{met}", f"{np.median(v):.1f}")
            macro(f"fCM{t}{met}Lo", f"{min(v):.1f}"); macro(f"fCM{t}{met}Hi", f"{max(v):.0f}")
    macro("fQDtime", f"{np.median([r['meas_qetu']['time'] / r['meas_dyadic']['time'] for r in D]):.2f}")
    macro("fQDsteps", f"{np.median([r['meas_qetu']['steps'] / r['meas_dyadic']['steps'] for r in D]):.2f}")
    macro("fRzCnot", f"{np.median([r['meas_dyadic']['rz'] / r['meas_dyadic']['cnot'] for r in D]):.2f}")
    macro("fCps", f"{np.median([r['meas_dyadic']['cnot'] / r['meas_dyadic']['steps'] / r['N'] for r in D]):.0f}")
    for e, t in ((1e-2, "Two"), (1e-3, "Three"), (1e-4, "Four")):
        macro(f"fDtD{t}", f"{np.median([r['meas_dyadic']['dt'] for r in D if r['eps'] == e]):.2f}")
        macro(f"fDtC{t}", f"{np.median([r['cert_dyadic']['dt'] for r in D if r['eps'] == e]):.3f}")
    # per-shot fits X ~ a N^al gap^-be L^la eps^-e
    def fit(meth, met):
        X, y = [], []
        for r in D:
            g, e = r["gamma"], r["eps"]
            X.append([1, np.log(r["N"]), np.log(r["gap"]), np.log(np.log((1 - g) / (g * e))), np.log(e)])
            y.append(np.log(r[meth][met]))
        c, *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
        return c
    for meth, t in (("meas_dyadic", "MD"), ("meas_qetu", "MQ"), ("cert_dyadic", "CD"), ("cert_qetu", "CQ")):
        for met, m in (("time", "T"), ("steps", "S"), ("cnot", "C")):
            c = fit(meth, met)
            macro(f"fit{t}{m}N", f"{c[1]:.2f}"); macro(f"fit{t}{m}gap", f"{-c[2]:.2f}")
            macro(f"fit{t}{m}L", f"{c[3]:.2f}"); macro(f"fit{t}{m}eps", f"{-c[4]:.2f}")
    # optimized design retries, and dyadic vs optimized vs QETU, measured and certified
    from derive_gapov_filter import _jsonl
    fix = _jsonl("gapov_opt_fix.jsonl")
    macro("optRetried", str(sum(1 for r in fix if "meas_optimized" in r)))
    macro("optCertRetried", str(sum(1 for r in fix if "cert_optimized" in r)))
    macro("optPoints", str(len(_jsonl("gapov_opt.jsonl"))))
    G = _gapov()
    for w, t in (("meas", "M"), ("cert", "C")):
        v = [r[f"{w}_optimized"]["cnot"] / r[f"{w}_dyadic"]["cnot"] for r in G if r[f"{w}_optimized"]]
        macro(f"od{t}Med", f"{np.median(v):.2f}"); macro(f"od{t}Lo", f"{min(v):.2f}")
        macro(f"od{t}Hi", f"{max(v):.2f}"); macro(f"od{t}Wins", str(sum(x < 1 for x in v)))
    v = [r["meas_optimized"]["p"] / r["gamma"] for r in G]
    macro("optPgLo", f"{min(v):.2f}"); macro("optPgMed", f"{np.median(v):.2f}")
    for e, t in ((1e-2, "Two"), (1e-3, "Three"), (1e-4, "Four")):
        R = [r for r in G if r["eps"] == e]
        best = [min(("dyadic", "optimized", "qetu"), key=lambda f: r[f"meas_{f}"]["cnot"]) for r in R]
        for f, ft in (("dyadic", "Dy"), ("optimized", "Op"), ("qetu", "Qe")):
            macro(f"best{ft}{t}", str(best.count(f)))
    Do = load("gapov_filter.json")
    for met in ("time", "steps", "cnot"):
        v = [r["meas_optimized"][met] / r["meas_dyadic"][met] for r in Do]
        macro(f"odShot{met.capitalize()}", f"{np.median(v):.2f}")
    macro("fDtOFour", f"{np.median([r['meas_optimized']['dt'] for r in Do if r['eps'] == 1e-4]):.2f}")
    macro("fDtOTwo", f"{np.median([r['meas_optimized']['dt'] for r in Do if r['eps'] == 1e-2]):.2f}")
    for meth, t in (("meas_optimized", "MO"), ("cert_optimized", "CO")):
        for met, m in (("time", "T"), ("steps", "S"), ("cnot", "C")):
            X, y = [], []
            for r in Do:
                g, e = r["gamma"], r["eps"]
                X.append([1, np.log(r["N"]), np.log(r["gap"]), np.log(np.log((1 - g) / (g * e))), np.log(e)])
                y.append(np.log(r[meth][met]))
            c, *_ = np.linalg.lstsq(np.array(X), np.array(y), rcond=None)
            macro(f"fit{t}{m}N", f"{c[1]:.2f}"); macro(f"fit{t}{m}gap", f"{-c[2]:.2f}")
            macro(f"fit{t}{m}L", f"{c[3]:.2f}"); macro(f"fit{t}{m}eps", f"{-c[4]:.2f}")
    v = [r["cert_optimized"]["time"] / r["meas_optimized"]["time"] for r in Do]
    macro("fCMOtime", f"{np.median(v):.1f}")
    v = [r["cert_optimized"]["steps"] / r["meas_optimized"]["steps"] for r in Do]
    macro("fCMOsteps", f"{np.median(v):.0f}")
    v = [r["cert_optimized"]["cnot"] / r["meas_optimized"]["cnot"] for r in G]
    rat = [r["cert_optimized"]["cnot"] / r["meas_optimized"]["cnot"] for r in G]
    macro("goCertRatioOptMed", f"{np.median(rat):.0f}"); macro("goCertRatioOptMin", f"{min(rat):.0f}")
    macro("goCertRatioOptMax", f"{max(rat):.0f}")
    # dyadic vs QETU, measured, by eps
    for e, t in ((1e-2, "Two"), (1e-3, "Three"), (1e-4, "Four")):
        R = [r for r in G if r["eps"] == e]
        macro(f"goQetuWins{t}", str(sum(r["meas_qetu"]["cnot"] < r["meas_dyadic"]["cnot"] for r in R)))
        macro(f"goCases{t}", str(len(R)))


# ----------------------------------------------------------------------------- oracle
def oracle():
    O = load("oracle.json")
    if O is None:
        return
    r = [o["ratio"] for o in O]
    macro("oracleMin", f"{min(r):.1f}"); macro("oracleMax", f"{max(r):.0f}")
    macro("oracleMed", f"{np.median(r):.0f}")
    lines = []
    for o in sorted(O, key=lambda o: (o["eps"], o["schedule"], o["J2"], o["N"]), reverse=False):
        pass
    for eps in (1e-2, 1e-3):
        for s in ("dyadic", "optimized"):
            cells = []
            for J2 in J2S:
                for N in (8, 10, 12):
                    q = [o for o in O if o["N"] == N and o["J2"] == J2 and o["eps"] == eps and o["schedule"] == s]
                    cells.append(f"{q[0]['ratio']:.0f}" if q else "--")
            lines.append(f"${sci(eps, 0)}$ & {s} & " + " & ".join(cells) + " \\\\")
    with open(os.path.join(GEN, "tab_oracle.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------- measured vs certified, robustness
CYLS = [("cyl", 3, 0.0), ("cyl", 3, 0.3), ("cyl", 4, 0.0), ("cyl", 4, 0.3)]
CHAINS9 = [("chain", N, J2) for N in (12, 14, 16) for J2 in (0.0, 0.2411, 0.4)]


def _qsd(n):
    return 0 if n <= 1 else round((23 / 48) * 4 ** n - 1.5 * 2 ** n + 4 / 3)


def _iso(m, n):
    """Lower bound on CNOTs for an isometry from m to n qubits (Iten et al., PRA 93, 032318 (2016))."""
    return max(0, int(np.ceil((2 ** (n + m + 1) - 4 ** m - 3 * n - 1) / 4)))


def _iso_lb(n):
    return _iso(n - 1, n)


def _pb(n):
    """CNOTs of generic n-qubit state preparation (Plesch and Brukner, PRA 83, 032302 (2011))."""
    if n % 2 == 0:
        return round(23 / 24 * 2 ** n - 1.5 * 2 ** (n // 2) + 4 / 3)
    return round(115 / 96 * 2 ** n - 2 ** ((n + 1) // 2) + 4 / 3)


def _bonds(N, chi):
    """Qubits carried across each cut of an open MPS with bond dimension <= chi."""
    return [0] + [int(np.ceil(np.log2(min(chi, 2 ** k, 2 ** (N - k))))) for k in range(1, N)] + [0]


def _mps_cnot(N, chi, model="qsd", f=1.0):
    """CNOTs of a chi-MPS circuit. qsd: the uniform staircase bound used in the study, (N - b) site unitaries
    on b + 1 qubits; ach: each site unitary sized by its actual bonds, capped by generic state preparation
    (both constructive); lb: the isometry lower bound per site with actual bonds, capped by the lower bound
    for generic state preparation (no compiler of these circuits can do better)."""
    b = int(np.ceil(np.log2(max(chi, 1))))
    if b == 0:
        return 0.0
    if model == "qsd":
        return f * (N - b) * _qsd(b + 1)
    q = _bonds(N, chi)
    if model == "ach":
        return f * min(sum(_qsd(q[k] + 1) for k in range(N)), _pb(N))
    return f * min(sum(_iso(q[k + 1], q[k] + 1) for k in range(N)), _iso(0, N))


def _generic(N, model):
    """Exact (eps = 0) state preparation as a baseline; none under the study's uniform model."""
    return {"qsd": np.inf, "ach": _pb(N), "lb": _iso(0, N)}[model]


def _crossover_counts(A, systems, f=1.0, model="qsd"):
    """Measured crossover with the MPS circuits re-costed (model, factor f): (wins, targets, max finite ratio)."""
    wins, tot, rat = 0, 0, []
    for sys_ in systems:
        R = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
        for e in sorted({r["eps"] for r in R}):
            Q = [r for r in R if r["eps"] == e]
            N = Q[0].get("N", Q[0]["size"])
            m = min([_mps_cnot(N, r["chi"], model, f) for r in Q if r["mps_alone"]] + [f * _generic(N, model)])
            best = np.inf
            for r in Q:
                if r["mps_alone"]:
                    continue
                for k in ("dyadic", "optimized"):
                    d = r.get(k)
                    if d:
                        filt = d["cnot"] * d["p"] - r["trial_cnot"]
                        best = min(best, (filt + _mps_cnot(N, r["chi"], model, f)) / d["p"])
            tot += 1
            if best < m:
                wins += 1
                if np.isfinite(m):
                    rat.append(m / best)
    return wins, tot, (max(rat) if rat else np.nan)


def _cert_crossover(systems):
    """Certified crossover from results/cert_crossover.jsonl: wins of the certified SBC ('auto') and QETU
    routes against the cheapest MPS alone, per target."""
    p = os.path.join(RES, "cert_crossover.jsonl")
    if not os.path.exists(p):
        return None
    C = [json.loads(l) for l in open(p)]
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + load("emp_actual.json")
    out = dict(sbc=0, qetu=0, tot=0, ratios=[], missing=0, best_sbc_over_mps=[])
    for sys_ in systems:
        R = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
        for e in sorted({r["eps"] for r in R}):
            m = min([r["trial_cnot"] for r in R if r["eps"] == e and r["mps_alone"]] or [np.inf])
            c = {s: [q["plan"]["res_expected_cnot"] for q in C if (q["kind"], q["size"], q["J2"]) == sys_
                     and q["eps"] == e and q["schedule"] == s and q.get("plan")] for s in ("auto", "qetu")}
            out["tot"] += 1
            bs = min(c["auto"] or [np.inf]); bq = min(c["qetu"] or [np.inf])
            out["sbc"] += bs < m; out["qetu"] += bq < m
            if np.isfinite(m) and np.isfinite(bs):
                out["best_sbc_over_mps"].append(bs / m)
            if bs < m and np.isfinite(m):
                out["ratios"].append(m / bs)
    return out


def meas_cert():
    """Table: every comparison measured (oracle-tuned) and certified (implementable)."""
    from derive_gapov_filter import load_cases
    rows = load_cases()
    lines = []
    for e, tag in ((1e-2, "Two"), (1e-3, "Three"), (1e-4, "Four")):
        rs = [r for r in rows if r["eps"] == e]
        mb = {"d": 0, "o": 0, "q": 0}; cb = {"d": 0, "o": 0, "q": 0}
        for r in rs:
            mc = {"d": r["meas_dyadic"]["cnot"], "q": (r["meas_qetu"] or {}).get("cnot", np.inf),
                  "o": (r["meas_optimized"] or {}).get("cnot", np.inf)}
            cc = {"d": r["cert_dyadic"]["cnot"], "q": r["cert_qetu"]["cnot"],
                  "o": (r["cert_optimized"] or {}).get("cnot", np.inf)}
            mb[min(mc, key=mc.get)] += 1; cb[min(cc, key=cc.get)] += 1
        macro(f"cbestDy{tag}", str(cb["d"])); macro(f"cbestOp{tag}", str(cb["o"])); macro(f"cbestQe{tag}", str(cb["q"]))
        lines.append(f"Cheapest filter, $\\epsilon={sci(e, 0)}$ (dy/opt/QETU) & {mb['d']}/{mb['o']}/{mb['q']} & "
                     f"{cb['d']}/{cb['o']}/{cb['q']} \\\\")
    macro("cQetuBestTot", str(sum(
        min({"d": r["cert_dyadic"]["cnot"], "q": r["cert_qetu"]["cnot"],
             "o": (r["cert_optimized"] or {}).get("cnot", np.inf)}.items(), key=lambda x: x[1])[0] == "q"
        for r in rows)))
    qd = [r["cert_qetu"]["cnot"] / r["cert_dyadic"]["cnot"] for r in rows]
    macro("cQetuOverDyMed", f"{np.median(qd):.2f}")
    macro("cQetuWinsDy", str(sum(x < 1 for x in qd)))
    A2 = load("emp_actual.json")
    A1 = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + [r for r in A2 if r["kind"] == "chain"]
    m1 = _crossover_counts(A1, CHAINS9); m2 = _crossover_counts(A2, CYLS)
    c1 = _cert_crossover(CHAINS9); c2 = _cert_crossover(CYLS)
    def cc(c):
        return "--" if c is None else f"{c['sbc']} of {c['tot']} (QETU {c['qetu']})"
    lines.insert(0, f"Optimized vs dyadic SBC, median ratio per success & {MEAS_OD[0]} & {MEAS_OD[1]} \\\\") if MEAS_OD else None
    lines.append(f"MPS + filter beats MPS alone, 9 chains & {m1[0]} of {m1[1]} & {cc(c1)} \\\\")
    lines.append(f"MPS + filter beats MPS alone, 4 cylinders & {m2[0]} of {m2[1]} & {cc(c2)} \\\\")
    with open(os.path.join(GEN, "tab_meascert.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    for c, t in ((c1, "One"), (c2, "Two")):
        if c:
            macro(f"cert{t}Wins", str(c["sbc"])); macro(f"cert{t}QWins", str(c["qetu"])); macro(f"cert{t}Tot", str(c["tot"]))
            r = np.array(c["best_sbc_over_mps"])
            macro(f"cert{t}OverMed", f"{np.median(r):.1f}" if len(r) else "--")
            macro(f"cert{t}OverMin", f"{np.min(r):.1f}" if len(r) else "--")
            macro(f"cert{t}RatioMax", f"{max(c['ratios']):.1f}" if c["ratios"] else "--")
    print("meas/cert:", m1, m2, c1 and {k: v for k, v in c1.items() if k not in ('ratios', 'best_sbc_over_mps')},
          c2 and {k: v for k, v in c2.items() if k not in ('ratios', 'best_sbc_over_mps')})
    # robustness of the measured crossovers to the MPS compiler
    rl = []
    for model, name in (("qsd", "QSD bound (used)"), ("lb", "isometry lower bound")):
        cells1, cells2 = [], []
        for f in (1, 1 / 2, 1 / 4, 1 / 8, 1 / 16):
            w2 = _crossover_counts(A2, CYLS, f, model); w1 = _crossover_counts(A1, CHAINS9, f, model)
            cells2.append(f"{w2[0]}"); cells1.append(f"{w1[0]}")
            if model == "lb" and f == 1:
                macro("robLbTwo", str(w2[0])); macro("robLbOne", str(w1[0]))
            if model == "qsd":
                macro("robQ" + {1: "One", 0.5: "Half", 0.25: "Quarter", 0.125: "Eighth", 0.0625: "Sixteenth"}[f] + "Two", str(w2[0]))
                macro("robQ" + {1: "One", 0.5: "Half", 0.25: "Quarter", 0.125: "Eighth", 0.0625: "Sixteenth"}[f] + "One", str(w1[0]))
        rl.append(f"{name}, cylinders (of 28) & " + " & ".join(cells2) + " \\\\")
        rl.append(f"{name}, chains (of 63) & " + " & ".join(cells1) + " \\\\")
    with open(os.path.join(GEN, "tab_robust.tex"), "w") as fh:
        fh.write("\n".join(rl) + "\n")


def _cert_counts(systems, model="qsd"):
    """Certified crossover with every MPS circuit (trial and MPS alone) re-costed under `model`:
    (SBC wins, QETU wins, targets, certified best filtered / best alternative per finite target)."""
    C = [json.loads(l) for l in open(os.path.join(RES, "cert_crossover.jsonl"))]
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + load("emp_actual.json")
    ws = wq = tot = 0
    over = []
    for sys_ in systems:
        R = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
        N = R[0]["N"]
        for e in sorted({r["eps"] for r in R}):
            m = min([_mps_cnot(N, r["chi"], model) for r in R if r["eps"] == e and r["mps_alone"]]
                    + [_generic(N, model)])
            b = {}
            for q in C:
                if (q["kind"], q["size"], q["J2"]) == sys_ and q["eps"] == e and q.get("plan"):
                    P = q["plan"]
                    c = (P["res_filter_cnot"] + _mps_cnot(N, q["chi"], model)) / P["res_p_succ"]
                    b[q["schedule"]] = min(b.get(q["schedule"], np.inf), c)
            tot += 1
            ws += b.get("auto", np.inf) < m
            wq += b.get("qetu", np.inf) < m
            if np.isfinite(m):
                over.append(min(b.values()) / m)
    return ws, wq, tot, over


def baselines():
    """Table: measured and certified crossovers under three MPS cost models, with exact generic state
    preparation as an extra baseline in the two bond-aware models."""
    A2 = load("emp_actual.json")
    A1 = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + [r for r in A2 if r["kind"] == "chain"]
    groups = (("One", A1, CHAINS9), ("Three", A2, CYLS[:2]), ("Four", A2, CYLS[2:]))
    names = {"qsd": "uniform QSD bound (study)", "ach": "bond-aware QSD, or generic prep.",
             "lb": "lower bounds (isometry, generic state)"}
    rows = []
    for model in ("qsd", "ach", "lb"):
        cells = []
        for tag, A, S in groups:
            w = _crossover_counts(A, S, 1.0, model)
            c = _cert_counts(S, model)
            cells.append(f"{w[0]}/{w[1]}" + (f" ({w[2]:.1f}$\\times$)" if np.isfinite(w[2]) else ""))
            cells.append(f"{c[0]}/{c[2]}")
            M = model.capitalize()
            macro(f"bl{M}{tag}", str(w[0])); macro(f"bl{M}{tag}Max", f"{w[2]:.1f}" if np.isfinite(w[2]) else "--")
            macro(f"bl{M}{tag}Cert", str(c[0])); macro(f"bl{M}{tag}CertQ", str(c[1]))
            macro(f"bl{M}{tag}CertMed", f"{np.median(c[3]):.1f}")
            print("baseline", model, tag, w, c[:3], f"cert/alt median {np.median(c[3]):.2f} min {np.min(c[3]):.2f}")
        rows.append(names[model] + " & " + " & ".join(cells) + " \\\\")
    with open(os.path.join(GEN, "tab_baselines.tex"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    for N in (12, 14, 16):
        macro("pb" + {12: "Twelve", 14: "Fourteen", 16: "Sixteen"}[N], f"{_pb(N):,}".replace(",", "\\,"))
        macro("gl" + {12: "Twelve", 14: "Fourteen", 16: "Sixteen"}[N], f"{_iso(0, N):,}".replace(",", "\\,"))


def qetu_tau():
    """QETU query-length sweep (results/qetu_tau.jsonl): cost ratio to tau = pi/2, cheapest measured filter in the
    controlled study with QETU at its best tau, and the measured crossovers with it."""
    p = os.path.join(RES, "qetu_tau.jsonl")
    if not os.path.exists(p):
        return
    Q = [json.loads(l) for l in open(p)]
    from derive_gapov_filter import load_cases
    rows = {(r["N"], round(r["J2"], 7), r["gamma_target"], r["eps"]): r for r in load_cases()}
    best = lambda q: min([x["cnot"] for x in (q["tau_pi2"], q["3pi/4"], q["pi/(1+delta)"]) if x] or [np.inf])
    G = [q for q in Q if q["src"] == "gapov"]
    rpi = [q["pi/(1+delta)"]["cnot"] / q["tau_pi2"]["cnot"] for q in Q if q["pi/(1+delta)"] and q["tau_pi2"]]
    r34 = [q["3pi/4"]["cnot"] / q["tau_pi2"]["cnot"] for q in Q if q["3pi/4"] and q["tau_pi2"]]
    macro("tauPiMed", f"{np.median(rpi):.2f}"); macro("tauPiMin", f"{np.min(rpi):.2f}"); macro("tauPiMax", f"{np.max(rpi):.2f}")
    macro("tauTqMed", f"{np.median(r34):.2f}")
    macro("tauNcases", str(len(Q))); macro("tauNgo", str(len(G)))
    lines = []
    for e, tag in ((1e-2, "Two"), (1e-3, "Three"), (1e-4, "Four")):
        old = {"d": 0, "o": 0, "q": 0}; new = {"d": 0, "o": 0, "q": 0}
        for q in G:
            if q["eps"] != e:
                continue
            r = rows[(q["N"], round(q["J2"], 7), q["gamma_target"], q["eps"])]
            mc = {"d": r["meas_dyadic"]["cnot"], "o": (r["meas_optimized"] or {}).get("cnot", np.inf),
                  "q": (r["meas_qetu"] or {}).get("cnot", np.inf)}
            old[min(mc, key=mc.get)] += 1
            mc["q"] = min(mc["q"], best(q))
            new[min(mc, key=mc.get)] += 1
        for k, n in (("d", "Dy"), ("o", "Op"), ("q", "Qe")):
            macro(f"tauBest{n}{tag}", str(new[k]))
        lines.append(f"$\\epsilon={sci(e, 0)}$ & {old['d']}/{old['o']}/{old['q']} & {new['d']}/{new['o']}/{new['q']} \\\\")
    with open(os.path.join(GEN, "tab_qetutau.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    # measured crossovers: is QETU at its best tau ever the cheapest route?
    A2 = load("emp_actual.json")
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + A2
    qa = {(q["kind"], q["size"], round(q["J2"], 9), q["chi"], q["eps"]): best(q) for q in Q if q["src"] == "actual"}
    out = {}
    for name, S in (("One", CHAINS9), ("Two", CYLS)):
        qwin = qbest = 0
        for sys_ in S:
            R = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
            for e in sorted({r["eps"] for r in R}):
                Rq = [r for r in R if r["eps"] == e]
                m = min([r["trial_cnot"] for r in Rq if r["mps_alone"]] or [np.inf])
                sb = min([r[k]["cnot"] for r in Rq if not r["mps_alone"] for k in ("dyadic", "optimized") if r.get(k)] or [np.inf])
                qb = min([qa.get((r["kind"], r["size"], round(r["J2"], 9), r["chi"], e), np.inf) for r in Rq
                          if not r["mps_alone"]] + [(r.get("qetu") or {}).get("cnot", np.inf) for r in Rq if not r["mps_alone"]])
                qwin += qb < m
                qbest += qb < min(m, sb)
        macro(f"tauX{name}Wins", str(qwin)); macro(f"tauX{name}Best", str(qbest))
        out[name] = (qwin, qbest)
    print("qetu tau:", len(Q), "pi/(1+d) med", np.median(rpi), "3pi/4 med", np.median(r34), lines, out)


def aposteriori():
    """Planned filters with an a-posteriori Trotter step (results/aposteriori.jsonl): reliability and cost of the
    energy-convergence and size-transfer rules against the certified step and the oracle step."""
    p = os.path.join(RES, "aposteriori.jsonl")
    if not os.path.exists(p):
        return
    import run_aposteriori as AP
    rows = [json.loads(l) for l in open(p)]
    for r in rows:
        r["curve"] = AP.dedupe(r["curve"])
    R = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}
    cost = lambda r, q: (q["filter_cnot"] + r["trial_cnot"]) / q["p"]
    lines, agg = [], {k: [] for k in ("e_ok", "e_cert", "e_or", "shots", "t_ok", "t_cert", "t_or", "or_cert", "e_ckcnot")}
    groups = ((("chain", 12), "chain, $N=12$"), (("chain", 14), "chain, $N=14$"), (("chain", 16), "chain, $N=16$"),
              (("cyl", 3), "$4\\times3$"), (("cyl", 4), "$4\\times4$"))
    for sched, sname in (("auto", "SBC"), ("qetu", "QETU")):
        for g, gname in groups:
            G = [r for r in rows if (r["kind"], r["size"]) == g and r["schedule"] == sched]
            if not G:
                continue
            st = {k: [] for k in agg}
            for r in G:
                c = r["curve"]
                o = AP.coarsest_pass(r)
                if o is not None:
                    st["or_cert"].append(cost(r, c[o]) / r["cert_cnot"])
                tol = r["eps"] * r["gap"] / 2
                i = next((i for i in range(1, len(c)) if abs(c[i]["E"] - c[i - 1]["E"]) <= tol), None)
                if i is not None:
                    st["e_ok"].append(c[i]["inf"] <= r["eps"])
                    st["e_cert"].append(cost(r, c[i]) / r["cert_cnot"])
                    if o is not None:
                        st["e_or"].append(cost(r, c[i]) / cost(r, c[o]))
                    sh = sum(2 * np.sqrt(q["var"]).sum() ** 2 / tol ** 2 / q["p"] for q in c[:i + 1])
                    st["shots"].append(sh)
                    st["e_ckcnot"].append(sum(2 * np.sqrt(q["var"]).sum() ** 2 / tol ** 2 / q["p"]
                                              * (q["filter_cnot"] + r["trial_cnot"]) for q in c[:i + 1]) / cost(r, c[i]))
                else:
                    st["e_ok"].append(False)
                dt = AP.transfer_dt(R, r)
                if dt is not None:
                    j = [q for q in c if q["dt"] <= dt + 1e-12]
                    if j:
                        q = max(j, key=lambda q: q["dt"])
                        st["t_ok"].append(q["inf"] <= r["eps"])
                        st["t_cert"].append(cost(r, q) / r["cert_cnot"])
                        if o is not None:
                            st["t_or"].append(cost(r, q) / cost(r, c[o]))
            for k in agg:
                agg[k] += [(sched, v) for v in st[k]]
            med = lambda k: f"{np.median(st[k]):.2f}" if st[k] else "--"
            tr = f"{sum(st['t_ok'])}/{len(st['t_ok'])}" if st["t_ok"] else "--"
            lines.append(f"{sname}, {gname} & {len(G)} & {med('or_cert')} & {sum(st['e_ok'])}/{len(st['e_ok'])} & "
                         f"{med('e_cert')} & {('$' + sci(np.median(st['shots']), 0) + '$') if st['shots'] else '--'} & {tr} & {med('t_cert')} \\\\")
    with open(os.path.join(GEN, "tab_apost.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    v = lambda k, s=None: [x for t, x in agg[k] if s is None or t == s]
    macro("apN", str(len(rows)))
    macro("apEok", str(sum(v("e_ok")))); macro("apEtot", str(len(v("e_ok"))))
    macro("apTok", str(sum(v("t_ok")))); macro("apTtot", str(len(v("t_ok"))))
    for s_, n in (("auto", "S"), ("qetu", "Q")):
        macro(f"apOrCert{n}", f"{np.median(v('or_cert', s_)):.2f}")
        macro(f"apEcert{n}", f"{np.median(v('e_cert', s_)):.2f}")
        macro(f"apEor{n}", f"{np.median(v('e_or', s_)):.1f}")
        if v("t_cert", s_):
            macro(f"apTcert{n}", f"{np.median(v('t_cert', s_)):.2f}")
            macro(f"apTor{n}", f"{np.median(v('t_or', s_)):.1f}")
    sh = v("shots")
    macro("apShotsMin", sci(np.min(sh), 0)); macro("apShotsMax", sci(np.max(sh), 0)); macro("apShotsMed", sci(np.median(sh), 0))
    macro("apCkMed", sci(np.median(v("e_ckcnot")), 0))
    # end to end: the planned filter with the a-posteriori step against the MPS alone (and generic preparation)
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + load("emp_actual.json")
    for model, M in (("qsd", "Q"), ("ach", "A")):
        for name, S in (("One", CHAINS9), ("Three", CYLS[:2]), ("Four", CYLS[2:])):
            w = {"e": 0, "t": 0}; tot = 0; ntr = 0; rat = []; ratt = []
            for sys_ in S:
                Rs = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
                N = Rs[0]["N"]
                for e in sorted({r["eps"] for r in Rs}):
                    m = min([_mps_cnot(N, r["chi"], model) for r in Rs if r["eps"] == e and r["mps_alone"]]
                            + [_generic(N, model)])
                    be = bt = np.inf
                    for sch in ("auto", "qetu"):
                        r = R.get(sys_ + (e, sch))
                        if r is None:
                            continue
                        c = r["curve"]
                        tr = _mps_cnot(N, r["chi"], model)
                        tol = r["eps"] * r["gap"] / 2
                        i = next((i for i in range(1, len(c)) if abs(c[i]["E"] - c[i - 1]["E"]) <= tol), None)
                        if i is not None and c[i]["inf"] <= e:
                            be = min(be, (c[i]["filter_cnot"] + tr) / c[i]["p"])
                        dt = AP.transfer_dt(R, r)
                        j = [q for q in c if dt is not None and q["dt"] <= dt + 1e-12]
                        if j:
                            q = max(j, key=lambda q: q["dt"])
                            if q["inf"] <= e:
                                bt = min(bt, (q["filter_cnot"] + tr) / q["p"])
                    tot += 1
                    w["e"] += be < m
                    if np.isfinite(bt):
                        ntr += 1
                        w["t"] += bt < m
                    if np.isfinite(m) and np.isfinite(be):
                        rat.append(m / be)
                    if np.isfinite(m) and np.isfinite(bt):
                        ratt.append(m / bt)
            macro(f"apX{M}{name}E", str(w["e"])); macro(f"apX{M}{name}T", str(w["t"]))
            macro(f"apX{M}{name}Tot", str(tot)); macro(f"apX{M}{name}Ttot", str(ntr))
            macro(f"apX{M}{name}Max", f"{max(rat):.1f}" if rat and max(rat) > 1 else "--")
            macro(f"apX{M}{name}TMax", f"{max(ratt):.1f}" if ratt and max(ratt) > 1 else "--")
            print("apost crossover", model, name, w, tot, ntr, f"max {max(rat) if rat else np.nan:.2f}", f"tmax {max(ratt) if ratt else np.nan:.2f}")
    print("aposteriori:", lines)


def ext_transfer():
    """Larger cylinders (4x6, 4x8, 4x12; DMRG inputs, planned filters from ext_regime_cyl*.json) with the Trotter
    step transferred from 4x4 (half its coarsest passing step, same J2, eps, schedule). Estimate, unvalidated at these
    sizes: S <= T/dt + m Strang steps (SBC) or m ceil(tau/dt) (QETU), success probability p_g = gamma F_pass^2."""
    import math
    import run_aposteriori as AP
    import run_extended as RX
    p = os.path.join(RES, "aposteriori.jsonl")
    if not os.path.exists(p):
        return
    rows = [json.loads(l) for l in open(p)]
    for r in rows:
        r["curve"] = AP.dedupe(r["curve"])
    Dd = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}
    E = load("ext_regime_cyl.json") + load("ext_regime_cyl_12.json")
    lines, wins, wins_lo, advs = [], {"qsd": 0, "ach": 0}, {"qsd": 0, "ach": 0}, []
    tot = only = 0
    for s in (6, 8, 12):
        for J2 in (0.0, 0.3):
            Bn = len(RX.model("cyl", s, J2).without_identity()) // 3
            N = 4 * s
            cells = []
            for e in (1e-2, 1e-3, 1e-4):
                Q = [r for r in E if r["size"] == s and r["J2"] == J2 and r["eps"] == e and r["plan"]]
                tot += 1
                res = {}
                # step factor relative to the 4x4 coarsest passing step: 1/2 is the rule; 1 (no safety factor) and
                # (1/2) sqrt(16/N) (a further N^(-1/2) shrink, the energy-shift scaling of Sec. apost) bracket it
                facs = {"mid": 0.5, "hi": 1.0, "lo": 0.5 * math.sqrt(16 / N)}
                for model in ("qsd", "ach"):
                    alone = min([_mps_cnot(N, r["chi"], model) for r in Q if r["plan"]["m"] == 0] or [np.inf])
                    best = {k: np.inf for k in facs}
                    for r in Q:
                        P = r["plan"]
                        if P["m"] == 0:
                            continue
                        sm = Dd[("cyl", 4, J2, e, r["schedule"])]
                        for k, fac in facs.items():
                            dt = sm["curve"][AP.coarsest_pass(sm)]["dt"] * fac
                            if r["schedule"] == "qetu":
                                S = P["m"] * math.ceil(math.pi / 2 / P["W"] / dt - 1e-12)
                            else:
                                S = math.ceil(P["T_phys"] / dt) + P["m"]
                            best[k] = min(best[k], (S * (12 * Bn - 12) + 6 * P["m"] + _mps_cnot(N, r["chi"], model)) / P["p_g"])
                    wins[model] += best["mid"] < alone
                    wins_lo[model] += best["lo"] < alone
                    res[model] = (alone, best)
                only += not np.isfinite(res["qsd"][0])
                a, b = res["ach"]
                if np.isfinite(a):
                    adv = f"{a / b['mid']:.2f} ({a / b['lo']:.2f}--{a / b['hi']:.2f})"
                    advs.append((a / b["lo"], a / b["hi"]))
                else:
                    adv = "only"
                cells.append(f"${sci(b['mid'], 1)}$ & " + adv)
            lines.append(f"$4\\times{s}$, $J_2={J2:g}$ & " + " & ".join(cells) + " \\\\")
    with open(os.path.join(GEN, "tab_exttransfer.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    macro("extTrQlo", str(wins_lo["qsd"])); macro("extTrAlo", str(wins_lo["ach"]))
    close = sum(1 for lo, hi in advs if lo < 1 < hi)
    macro("extTrStraddle", str(close))
    macro("extTrQ", str(wins["qsd"])); macro("extTrA", str(wins["ach"])); macro("extTrTot", str(tot)); macro("extTrOnly", str(only))
    print("ext transfer:", wins, wins_lo, tot, only, "straddle", close)


def _iso_c(a, b):
    """Isometry lower bound (Iten et al.) for dimensions a -> b, not only powers of two; equals _iso there."""
    return max(0, int(np.ceil((2 * a * b - a * a - 3 * np.log2(b) - 1) / 4)))


def _mps_lbd(N, chi):
    """Lower bound for a chi-MPS circuit with the actual (any integer) bond dimensions D_k = min(chi, 2^k, 2^(N-k)):
    the per-site isometry bound D_(k+1) -> 2 D_k, capped by the generic-state bound. Most generous to the MPS."""
    D = [1] + [min(chi, 2 ** k, 2 ** (N - k)) for k in range(1, N)] + [1]
    return float(min(sum(_iso_c(D[k + 1], 2 * D[k]) for k in range(N)), _iso(0, N)))


def chi_mid():
    """MPS-alone route with intermediate bond dimensions (results/chi_mid.json), against the measured filtered route
    and the planned filter with the size-transfer step. Models: qsd and ach (cost set by ceil(log2 chi), so an
    intermediate chi costs as much as the next power of two) and lbd (the dimension-resolved lower bound)."""
    p = os.path.join(RES, "chi_mid.json")
    if not os.path.exists(p):
        return
    import run_aposteriori as AP
    import run_empirical as RE
    INF = json.load(open(p))
    A2 = load("emp_actual.json")
    A1 = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))] + [r for r in A2 if r["kind"] == "chain"]
    rows = [json.loads(l) for l in open(os.path.join(RES, "aposteriori.jsonl"))]
    for r in rows:
        r["curve"] = AP.dedupe(r["curve"])
    R = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}
    grid = set(RE.ACT_CHIS)

    def cost(model, N, chi):
        return _mps_lbd(N, chi) if model == "lbd" else _mps_cnot(N, chi, model)

    def gen(model, N):
        return _generic(N, "lb" if model == "lbd" else model)

    out, lines = {}, []
    for model, mname in (("qsd", "uniform QSD bound"), ("ach", "bond-aware QSD"), ("lbd", "dimension lower bound")):
        cells = {}
        for allchi in (False, True):
            for tag, A, S in (("One", A1, CHAINS9), ("Three", A2, CYLS[:2]), ("Four", A2, CYLS[2:])):
                wm = wt = tot = ntr = 0
                for sys_ in S:
                    Rs = [r for r in A if (r["kind"], r["size"], r["J2"]) == sys_]
                    N = Rs[0].get("N", Rs[0]["size"])
                    inf = {int(c): v for c, v in INF[RE._act_key(*sys_)].items()}
                    for e in sorted({r["eps"] for r in Rs}):
                        ch = [c for c in inf if inf[c] <= e and (allchi or c in grid)]
                        m = min([cost(model, N, c) for c in ch] + [gen(model, N)])
                        trial = lambda c: _mps_cnot(N, c, "lb" if model == "lbd" else model)
                        bm = np.inf
                        for r in Rs:
                            if r["eps"] != e or r["mps_alone"]:
                                continue
                            for k in ("dyadic", "optimized"):
                                d = r.get(k)
                                if d:
                                    bm = min(bm, (d["cnot"] * d["p"] - r["trial_cnot"] + trial(r["chi"])) / d["p"])
                        bt = np.inf
                        for sch in ("auto", "qetu"):
                            r = R.get(sys_ + (e, sch))
                            dt = AP.transfer_dt(R, r) if r else None
                            j = [q for q in r["curve"] if q["dt"] <= dt + 1e-12] if dt is not None else []
                            if j:
                                q = max(j, key=lambda q: q["dt"])
                                if q["inf"] <= e:
                                    bt = min(bt, (q["filter_cnot"] + trial(r["chi"])) / q["p"])
                        tot += 1
                        wm += bm < m
                        if np.isfinite(bt):
                            ntr += 1
                            wt += bt < m
                cells[(allchi, tag)] = (wm, tot, wt, ntr)
                M = {"qsd": "Q", "ach": "A", "lbd": "L"}[model]
                macro(f"cm{M}{tag}{'All' if allchi else 'Grid'}", str(wm))
                macro(f"cm{M}{tag}{'All' if allchi else 'Grid'}T", str(wt))
                if model == "qsd" and not allchi:
                    macro(f"cm{tag}Tot", str(tot)); macro(f"cm{tag}Ttot", str(ntr))
        for allchi in (False, True):
            lab = mname + (", all $\\chi\\le64$" if allchi else ", $\\chi=2^k$")
            lines.append(lab + " & " + " & ".join(f"{cells[(allchi, t)][0]} & "
                                                  + (f"{cells[(allchi, t)][2]}" if cells[(allchi, t)][3] else "--")
                                                  for t in ("One", "Three", "Four")) + " \\\\")
        out[model] = cells
    with open(os.path.join(GEN, "tab_chimid.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("chi mid:", out)


def cross20():
    """Implementable crossover on the 4x5 cylinder (results/cross20.jsonl, cross20_mps.json): planned filter with the
    transferred step against the MPS alone, for the 14 targets of the 4x4 crossover."""
    p = os.path.join(RES, "cross20.jsonl")
    if not os.path.exists(p):
        return
    rows = [json.loads(l) for l in open(p)]
    M = json.load(open(os.path.join(RES, "cross20_mps.json")))
    grid = {2, 4, 8, 16, 32, 64}
    ok = [r["curve"][1]["inf"] <= r["eps"] for r in rows if len(r["curve"]) > 1]
    macro("xfOk", str(sum(ok))); macro("xfTot", str(len(ok)))
    macro("xfPrev", str(sum(r["curve"][0]["inf"] <= r["eps"] for r in rows if r["curve"])))
    trial = lambda model, c: _mps_cnot(20, c, "lb" if model == "lbd" else model)
    cost = lambda model, c: _mps_lbd(20, c) if model == "lbd" else _mps_cnot(20, c, model)
    gen = lambda model: _generic(20, "lb" if model == "lbd" else model)
    out = {}
    for model, Mn in (("qsd", "Q"), ("ach", "A"), ("lbd", "L")):
        for tag, allchi, cap in (("Grid", False, 64), ("All", True, 64), ("Big", True, 256)):
            w = tot = 0
            for J2 in (0.0, 0.3):
                inf = {int(c): v for c, v in M[f"{J2:g}"].items()}
                for e in sorted({r["eps"] for r in rows}):
                    ch = [c for c in inf if inf[c] <= e and c <= cap and (allchi or c in grid)]
                    m = min([cost(model, c) for c in ch] + [gen(model)])
                    bt = np.inf
                    for r in rows:
                        if r["J2"] == J2 and r["eps"] == e and len(r["curve"]) > 1 and r["curve"][1]["inf"] <= e:
                            q = r["curve"][1]
                            bt = min(bt, (q["filter_cnot"] + trial(model, r["chi"])) / q["p"])
                    tot += 1
                    w += bt < m
            macro(f"xf{Mn}{tag}", str(w))
            out[(model, tag)] = (w, tot)
    macro("xfTargets", str(tot))
    print("cross20:", out, sum(ok), len(ok))


def margin():
    """Gap margin from DMRG convergence (results/margin.json, margin_behave.jsonl): error bars equal to the last change
    along the bond-dimension ladder, tested against exact energies and, on the systems with a full spectrum, by the
    exact filter function."""
    p = os.path.join(RES, "margin.json")
    if not os.path.exists(p):
        return
    R = json.load(open(p))
    E0x = {}
    for l in open(os.path.join(RES, "aposteriori.jsonl")):
        r = json.loads(l)
        E0x[(r["kind"], r["size"], r["J2"])] = r["E0"]
    for J2 in (0.0, 0.3):
        E0x[("cyl", 5, J2)] = json.load(open(os.path.join(RES, "cache_t20", f"cyl5_J2{J2:g}.json")))["E01"][0]
    cov = {0: 0, 1: 0}; n = 0; over = []
    for r in R:
        e0 = E0x[(r["kind"], r["size"], r["J2"])]; e1 = e0 + r["gap_exact"]
        for i in range(1, len(r["chis"])):
            d0 = abs(r["E0"][i] - r["E0"][i - 1]); d1 = abs(r["E1"][i] - r["E1"][i - 1])
            n += 1
            over.append((r["E1"][i] - r["E0"][i]) / r["gap_exact"])
            for c in cov:
                cov[c] += (r["E1"][i] - c * d1 <= e1 + 1e-9) and (r["E0"][i] - c * d0 <= e0 + 1e-9)
    macro("mgN", str(n)); macro("mgSys", str(len(R))); macro("mgCovRaw", str(cov[0])); macro("mgCov", str(cov[1]))
    macro("mgOverMax", f"{100 * (max(over) - 1):.0f}")
    B = [json.loads(l) for l in open(os.path.join(RES, "margin_behave.jsonl"))]
    D = {(r["kind"], r["size"], r["J2"], r["chi_dmrg"], r["c"], r["chi"], r["eps"], r["schedule"]): r for r in B}
    ok = [r["inf"] <= r["eps"] for r in B if "inf" in r]
    macro("mbPlans", str(len(ok))); macro("mbPass", str(sum(ok)))
    def ratios(lo, hi):
        x = [r["cert_cnot"] / D[k[:4] + (0,) + k[5:]]["cert_cnot"] for k, r in D.items()
             if r["c"] == 1 and lo <= r["chi_dmrg"] <= hi and "cert_cnot" in r]
        return x
    hi = ratios(64, 10 ** 6); mid = ratios(32, 32)
    macro("mbHiMed", f"{np.median(hi):.2f}"); macro("mbHiMax", f"{max(hi):.1f}"); macro("mbMidMed", f"{np.median(mid):.1f}")
    lo = [r for r in B if r["c"] == 1 and r["chi_dmrg"] <= 16]
    macro("mbLoInf", str(sum("error" in r for r in lo))); macro("mbLoN", str(len(lo)))
    macro("mbSys", str(len({(r["kind"], r["size"], r["J2"]) for r in B})))
    print("margin:", cov, n, len(ok), sum(ok), np.median(hi), max(hi), np.median(mid))


def transfer20():
    """Size transfer one step beyond any earlier simulation: chain 16 -> 20 and cylinder 4x4 -> 4x5, simulated
    (results/transfer20.jsonl). Also the cross-J2 test of the same rule on results/aposteriori.jsonl."""
    import run_aposteriori as AP
    rows = [json.loads(l) for l in open(os.path.join(RES, "aposteriori.jsonl"))]
    for r in rows:
        r["curve"] = AP.dedupe(r["curve"])
    D = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}
    ok = tot = beyond = changed = 0
    for (k, s_, j, e, sc), r in D.items():
        for (k2, s2, j2, e2, sc2), r2 in D.items():
            if (k2, s2, e2, sc2) != (k, s_, e, sc) or j2 == j:
                continue
            dt = r2["curve"][AP.coarsest_pass(r2)]["dt"] / 2
            changed += r["curve"][AP.coarsest_pass(r)]["dt"] != r2["curve"][AP.coarsest_pass(r2)]["dt"]
            c = [q for q in r["curve"] if abs(q["dt"] - dt) < 1e-12]
            if c:
                tot += 1
                ok += c[0]["inf"] <= e
            else:
                beyond += 1
    macro("jtOk", str(ok)); macro("jtTot", str(tot)); macro("jtChanged", str(changed)); macro("jtBeyond", str(beyond))
    # step versus size among the simulated sizes
    sh = tt = 0
    for k, a, b in (("chain", 12, 14), ("chain", 14, 16), ("cyl", 3, 4)):
        for (kk, s_, j, e, sc), r in D.items():
            if kk == k and s_ == a:
                r2 = D[(k, b, j, e, sc)]
                tt += 1
                sh += r2["curve"][AP.coarsest_pass(r2)]["dt"] < r["curve"][AP.coarsest_pass(r)]["dt"]
    macro("szShrink", str(sh)); macro("szTot", str(tt))
    print("cross-J2 transfer:", ok, tot, "beyond", beyond, "step changed", changed, "| size shrink", sh, tt)
    p = os.path.join(RES, "transfer20.jsonl")
    if not os.path.exists(p):
        return
    T = [json.loads(l) for l in open(p)]
    T = [t for t in T if t["curve"]]
    lines = []
    allok = allprev = alln = 0
    cc = []
    for kind, size, name in (("chain", 20, "chain, $N=20$ (from 16)"), ("cyl", 5, "$4\\times5$ (from $4\\times4$)")):
        for sched, sn in (("auto", "SBC"), ("qetu", "QETU")):
            G = [t for t in T if t["kind"] == kind and t["schedule"] == sched]
            if not G:
                continue
            prev = sum(t["curve"][0]["inf"] <= t["eps"] for t in G)
            tr = [next((q for q in t["curve"] if abs(q["dt"] - t["dt_prev"] / 2) < 1e-12), None) for t in G]
            okc = sum(q is not None and q["inf"] <= t["eps"] for q, t in zip(tr, G))
            rc = [q["cnot"] / t["cert_cnot"] for q, t in zip(tr, G) if q is not None]
            worst = max(q["inf"] / t["eps"] for q, t in zip(tr, G) if q is not None)
            cc += rc
            allok += okc; allprev += prev; alln += len(G)
            lines.append(f"{sn}, {name} & {len(G)} & {prev}/{len(G)} & {okc}/{len(G)} & {worst:.2f} & {np.median(rc):.2f} \\\\")
    with open(os.path.join(GEN, "tab_transfer20.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    macro("twOk", str(allok)); macro("twTot", str(alln)); macro("twPrev", str(allprev))
    macro("twCert", f"{np.median(cc):.2f}")
    print("transfer20:", lines)


MEAS_OD = None


if __name__ == "__main__":
    for f in (inputs_table, validation, scaling, qetu, drivers, empirical, regime, fid_plane, crossover, actual_crossover, actual1d, gapov, gapov_filter, sens_wide, sensitivity, ansatz, oracle, meas_cert, baselines, qetu_tau, aposteriori, ext_transfer, chi_mid, transfer20, cross20, margin):
        f()
    with open(os.path.join(GEN, "macros.tex"), "w") as fh:
        fh.write("% generated by sbpipe_study/make_figures.py -- do not edit\n" + "\n".join(MACROS) + "\n")
    print("\n".join(MACROS))
