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
    macro("boundRatioMin", f"{min(ratios):.0f}")
    macro("boundRatioMax", sci(max(ratios)))
    macro("boundRatioMed", f"{np.median(ratios):.0f}")
    lr = [v["plan"]["leak"] / max(v["sim"]["infidelity_exact"], 1e-300) for v in sims]
    macro("leakRatioMin", f"{min(lr):.1f}")
    macro("leakRatioMax", f"{max(lr):.0f}")
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
        ax[0].plot(x, v["sim"]["infidelity_trotter"] / v["eps"], **kw)
        ax[1].plot(v["sim"]["trotter_dist_bound"], max(v["sim"]["trotter_state_dist"], 1e-9), **kw)
    ax[0].axhline(1.0, color=INK2, lw=0.8, ls="--")
    ax[0].text(7.4, 1.25, "target / bound", color=INK2, fontsize=7)
    ax[0].set_yscale("log"); ax[0].set_ylim(1e-6, 3)
    ax[0].set_xticks([8, 10, 12, 14]); ax[0].set_xlim(7, 15)
    ax[0].set_xlabel("$N$"); ax[0].set_ylabel(r"measured $1-F$ / target $\epsilon$")
    lim = np.array([1e-9, 1])
    ax[1].plot(lim, lim, color=INK2, lw=0.8, ls="--")
    ax[1].text(4e-3, 6e-3, "bound = measured", rotation=28, color=INK2, fontsize=7)
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


if __name__ == "__main__":
    for f in (inputs_table, validation, scaling, sensitivity, ansatz, oracle):
        f()
    with open(os.path.join(GEN, "macros.tex"), "w") as fh:
        fh.write("% generated by sbpipe_study/make_figures.py -- do not edit\n" + "\n".join(MACROS) + "\n")
    print("\n".join(MACROS))
