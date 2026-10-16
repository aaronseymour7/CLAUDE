"""build_results_tex.py -- build results/RESULTS.tex (a self-contained LaTeX report) and its figures from results/data/*.json,
results/dmrg_validation.json and results/bandwidth_comparison.json.  Same data as RESULTS.md; every quantity is defined in Section 1.

usage:  python scripts/build_results_tex.py        (writes results/RESULTS.tex and results/figures_tex/*.png)
"""
import glob, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullFormatter

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures_tex"
FIG.mkdir(exist_ok=True)

R = {}
for f in glob.glob(str(RES / "data" / "dmrg_N*_eps*.json")):
    r = json.load(open(f))
    R[(r["N"], round(r["J2"], 4), r["eps"])] = r
VAL = json.load(open(RES / "dmrg_validation.json")) if (RES / "dmrg_validation.json").exists() else {}
BWC = json.load(open(RES / "bandwidth_comparison.json")) if (RES / "bandwidth_comparison.json").exists() else None
J2S = (0.0, 0.2411, 0.4)
J2NAME = {0.0: "0", 0.2411: "0.2411", 0.4: "0.4"}
COL = {0.0: "#1f4e79", 0.2411: "#2e8b57", 0.4: "#c0392b"}
MK = {0.0: "o", 0.2411: "s", 0.4: "^"}
EPS = 0.01


# ------------------------------------------------------------------------------------------------ helpers
def one(N, J2, eps=EPS):
    return R.get((N, round(J2, 4), eps))


def des(N, J2, which="guaranteed", eps=EPS):
    r = one(N, J2, eps)
    d = None if r is None else r["designs"].get(which)
    return d if d and "n" in d else None


def num(x, d=0):
    """number with thin-space thousands separators."""
    if x is None:
        return "--"
    s = f"{x:,.{d}f}"
    return s.replace(",", "\\,")


def sci(x, d=1):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    x = float(x)
    if x == 0:
        return "$0$"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    if e == 0:
        return f"${m:.{d}f}$"
    return f"${m:.{d}f}\\times10^{{{e}}}$"


def fx(x, d=3):
    return "--" if x is None else f"{x:.{d}f}"


def tab(header, rows, colspec, caption, label, small=True, long=True):
    """booktabs longtable (allows page breaks)."""
    out = []
    out.append("\\begin{longtable}{" + colspec + "}" if long else "\\begin{tabular}{" + colspec + "}")
    out.append("\\caption{" + caption + "}\\label{" + label + "}\\\\" if long else "")
    out.append("\\toprule")
    out.append(" & ".join(header) + " \\\\")
    out.append("\\midrule")
    if long:
        out.append("\\endfirsthead")
        out.append("\\toprule")
        out.append(" & ".join(header) + " \\\\")
        out.append("\\midrule")
        out.append("\\endhead")
        out.append("\\bottomrule")
        out.append("\\endfoot")
    for r in rows:
        out.append(" & ".join(str(c) for c in r) + " \\\\")
    if long:
        out.append("\\end{longtable}")
    else:
        out.append("\\bottomrule\n\\end{tabular}")
    body = "\n".join(x for x in out if x != "")
    return ("{\\small\n" if small else "{\n") + body + "\n}\n"


def figure(name, caption, label, width="0.95\\linewidth"):
    return (f"\\begin{{figure}}[htbp]\n\\centering\n\\includegraphics[width={width}]{{figures_tex/{name}}}\n"
            f"\\caption{{{caption}}}\\label{{{label}}}\n\\end{{figure}}\n")


plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25})


def plain_x(ax, ticks):
    ax.set_xticks(ticks); ax.xaxis.set_major_formatter(ScalarFormatter()); ax.xaxis.set_minor_formatter(NullFormatter())


def legend_j2(ax, **kw):
    ax.legend(frameon=False, fontsize=8, **kw)


# ------------------------------------------------------------------------------------------------ figures
# F1: DMRG accuracy
fig, ax = plt.subplots(figsize=(6.2, 3.4))
for J2 in J2S:
    xs, e0, e1 = [], [], []
    for N in (4, 6, 8, 10, 12, 14):
        r = one(N, J2)
        e = (r or {}).get("inputs", {}).get("dmrg_vs_ed")
        if e:
            xs.append(N); e0.append(max(e["E0_diff"], 1e-17)); e1.append(max(e["E1_diff"], 1e-17))
    for k, v in VAL.items():
        if k.startswith("ed_") and abs(v["J2"] - J2) < 1e-9:
            xs.append(v["N"]); e0.append(float(v["dE0"])); e1.append(float(v["dE1"]))
    o = np.argsort(xs)
    xs, e0, e1 = np.array(xs)[o], np.array(e0)[o], np.array(e1)[o]
    ax.semilogy(xs, e0, MK[J2] + "-", color=COL[J2], ms=4, label=f"$|E_0-E_0^{{\\rm ED}}|$, $J_2={J2NAME[J2]}$")
    ax.semilogy(xs, e1, MK[J2] + "--", color=COL[J2], ms=4, mfc="white", label=f"$|E_1-E_1^{{\\rm ED}}|$, $J_2={J2NAME[J2]}$")
    cs = sorted((v["N"], max(float(v["dE0"]), float(v["dE1"]))) for k, v in VAL.items() if k.startswith("consist_") and abs(v["J2"] - J2) < 1e-9)
    if cs:
        ax.semilogy([c[0] for c in cs], [c[1] for c in cs], MK[J2], color=COL[J2], ms=7, mfc="none", mew=1.3)
ax.axhline(1e-13, color="k", lw=0.8, ls=":")
ax.set_xlabel("number of spins $N$"); ax.set_ylabel("DMRG error (energy units, $J_1=1$)")
ax.set_title("DMRG accuracy: against exact diagonalisation ($N\\leq20$), run-to-run (large open markers, $N\\geq24$)", fontsize=8)
ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper left")
ax.set_ylim(1e-17, 1e-8)
fig.tight_layout(); fig.savefig(FIG / "dmrg_accuracy.png", dpi=200); plt.close(fig)

# F2: trial state weight
fig, ax = plt.subplots(figsize=(6.2, 3.2))
for J2 in J2S:
    pts = [(N, one(N, J2)["inputs"]["gamma"]) for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32) if one(N, J2)]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], MK[J2] + "-", color=COL[J2], ms=4, label=f"$J_2={J2NAME[J2]}$")
    bad = [N for N in (28, 32) if one(N, J2) and des(N, J2) is None]
    ax.plot(bad, [one(N, J2)["inputs"]["gamma"] for N in bad], "x", color="k", ms=8, mew=1.6)
ax.plot([], [], "x", color="k", ms=8, mew=1.6, label="no certified design at $\\varepsilon=10^{-2}$")
ax.set_xlabel("number of spins $N$"); ax.set_ylabel("ground-state weight $\\gamma$ of the trial state")
ax.legend(frameon=False, fontsize=8); ax.set_ylim(0.25, 1.02)
fig.tight_layout(); fig.savefig(FIG / "gamma_vs_N.png", dpi=200); plt.close(fig)

# F3: steps guaranteed vs empirical
fig, ax = plt.subplots(figsize=(6.2, 3.4))
for J2 in J2S:
    g = [(N, des(N, J2)["n"]) for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32) if des(N, J2)]
    e = [(N, des(N, J2, "empirical")["n"]) for N in (4, 6, 8, 10, 12, 14) if des(N, J2, "empirical")]
    ax.semilogy([p[0] for p in g], [p[1] for p in g], MK[J2] + "-", color=COL[J2], ms=4, label=f"guaranteed, $J_2={J2NAME[J2]}$")
    ax.semilogy([p[0] for p in e], [p[1] for p in e], MK[J2] + "--", color=COL[J2], ms=4, mfc="white", label=f"empirical (oracle), $J_2={J2NAME[J2]}$")
ax.set_xlabel("number of spins $N$"); ax.set_ylabel("Trotter steps $n$ for $\\varepsilon=10^{-2}$")
ax.legend(frameon=False, fontsize=6.5, ncol=2)
fig.tight_layout(); fig.savefig(FIG / "steps_vs_N.png", dpi=200); plt.close(fig)

# F4: RZZ and non-Clifford rotations (Rz + RZZ), log-log, both compilations; fits on level 3 RZZ
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))
FITTXT = {}
for key, ax, lab in (("cx", axes[0], "total RZZ"), ("rz", axes[1], "total non-Clifford rotations")):
    for J2 in J2S:
        p3, pd = [], []
        for N in (6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
            d = des(N, J2)
            if d and one(N, J2)["design"]["source"] in ("floor", "none"):
                if key == "cx":
                    p3.append((N, d["counts_rzz_a2a"]["total"]["rzz"])); pd.append((N, d["counts_rzz_direct"]["filter"]["rzz"]))
                else:
                    p3.append((N, d["counts_rzz_a2a"]["total"]["nc_total"])); pd.append((N, d["counts_rzz_direct"]["filter"]["nc_total"]))
        if not p3:
            continue
        ax.loglog([p[0] for p in p3], [p[1] for p in p3], MK[J2], color=COL[J2], ms=4)
        ax.loglog([p[0] for p in pd], [p[1] for p in pd], MK[J2], color=COL[J2], ms=5, mfc="white")
        for pts_, ls in ((p3, "-"), (pd, ":")):
            sl, c = np.polyfit(np.log([p[0] for p in pts_]), np.log([p[1] for p in pts_]), 1)
            xx = np.array([pts_[0][0], pts_[-1][0]])
            ax.loglog(xx, np.exp(c) * xx ** sl, ls, color=COL[J2], lw=0.9, alpha=0.7)
            FITTXT[(key, J2, ls)] = (sl, pts_[0][0], pts_[-1][0])
    ax.set_xlabel("number of spins $N$"); ax.set_ylabel(lab); plain_x(ax, [6, 10, 16, 24, 32])
for J2 in J2S:
    axes[0].plot([], [], MK[J2], color=COL[J2], ms=4, label=f"$J_2={J2NAME[J2]}$")
axes[0].plot([], [], "o", color="k", ms=4, label="level 3 (re-synthesised)")
axes[0].plot([], [], "o", color="k", ms=5, mfc="white", label="direct")
axes[0].legend(frameon=False, fontsize=6.5)
fig.tight_layout(); fig.savefig(FIG / "cost_vs_N.png", dpi=200); plt.close(fig)

# F5: gap and physical time
fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2))
for J2 in J2S:
    gp = [(N, one(N, J2)["inputs"]["gap_raw"]) for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32) if one(N, J2)]
    axes[0].loglog([p[0] for p in gp], [p[1] for p in gp], MK[J2] + "-", color=COL[J2], ms=4, label=f"$J_2={J2NAME[J2]}$")
    tp = [(N, des(N, J2)["T_phys"]) for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32) if des(N, J2)]
    axes[1].plot([p[0] for p in tp], [p[1] for p in tp], MK[J2] + "-", color=COL[J2], ms=4, label=f"$J_2={J2NAME[J2]}$")
axes[0].set_xlabel("number of spins $N$"); axes[0].set_ylabel("energy gap $g=E_1-E_0$ ($J_1=1$)"); plain_x(axes[0], [4, 8, 16, 32])
axes[1].set_xlabel("number of spins $N$"); axes[1].set_ylabel("physical evolution time $T_{\\rm phys}$ ($1/J_1$)")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "gap_and_time.png", dpi=200); plt.close(fig)

# F6: measured infidelity of the guaranteed designs against the proven bound
fig, ax = plt.subplots(figsize=(5.6, 3.3))
for J2 in J2S:
    pts = [(N, des(N, J2)["eps_actual"]) for N in (4, 6, 8, 10, 12, 14) if des(N, J2) and des(N, J2).get("simulated") and des(N, J2).get("eps_actual")]
    ax.semilogy([p[0] for p in pts], [p[1] for p in pts], MK[J2] + "-", color=COL[J2], ms=5, label=f"measured, $J_2={J2NAME[J2]}$")
ax.axhline(EPS, color="k", lw=1.0, ls="--", label="proven bound $\\approx\\varepsilon=10^{-2}$")
ax.set_xlabel("number of spins $N$"); ax.set_ylabel("infidelity $1-F$ of the guaranteed design")
ax.legend(frameon=False, fontsize=8, loc="center right"); ax.set_ylim(1e-4, 3e-2)
fig.tight_layout(); fig.savefig(FIG / "bound_vs_measured.png", dpi=200); plt.close(fig)

# F7: epsilon scan
fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.0), sharey=False)
for ax, (N, J2) in zip(axes, ((6, 0.0), (8, 0.0), (6, 0.4))):
    es = sorted({e for (n_, j_, e) in R if n_ == N and abs(j_ - J2) < 1e-9}, reverse=True)
    for which, ls in (("guaranteed", "-"), ("empirical", "--")):
        pts = [(e, des(N, J2, which, e)["n"]) for e in es if des(N, J2, which, e)]
        ax.loglog([p[0] for p in pts], [p[1] for p in pts], "o" + ls, color="#1f4e79" if which == "guaranteed" else "#c0392b", ms=4, label=which)
    ax.set_title(f"$N={N}$, $J_2={J2NAME[J2]}$", fontsize=9); ax.set_xlabel("target $\\varepsilon$"); ax.invert_xaxis()
axes[0].set_ylabel("Trotter steps $n$"); axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "eps_scan.png", dpi=200); plt.close(fig)

# F8: per step cost
fig, ax = plt.subplots(figsize=(5.4, 3.2))
for J2 in J2S:
    pts = [(N, one(N, J2)["per_step"]["rzz"]["a2a"]["one"]["rzz"], one(N, J2)["per_step"]["rzz"]["line"]["avg3"]["rzz"]) for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32) if one(N, J2)]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], MK[J2] + "-", color=COL[J2], ms=4, label=f"all-to-all, $J_2={J2NAME[J2]}$")
    ax.plot([p[0] for p in pts], [p[2] for p in pts], MK[J2] + ":", color=COL[J2], ms=4, mfc="white", label=f"line, $J_2={J2NAME[J2]}$")
ax.set_xlabel("number of spins $N$"); ax.set_ylabel("RZZ per Trotter step of one pulse")
ax.legend(frameon=False, fontsize=6.5, ncol=2)
fig.tight_layout(); fig.savefig(FIG / "rzz_per_step.png", dpi=200); plt.close(fig)


# ------------------------------------------------------------------------------------------------ tables
T = {}

# inputs
rows = []
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
        r = one(N, J2)
        if not r:
            continue
        i = r["inputs"]
        rows.append([N, J2NAME[J2], fx(i["E0"], 6), fx(i["gap_raw"], 5), fx(i["W"], 3), fx(i["gap_scaled"], 5), fx(i["gamma"], 4), sci(1 - i["gamma"], 2), fx(i["alpha"], 4), i["n_pauli_terms"], i["trial_cost_rzz"]["rzz"] and int(i["trial_cost_rzz"]["rzz"])])
T["inputs"] = tab(["$N$", "$J_2$", "$E_0$", "$g=E_1-E_0$", "$W$", "$\\Delta=g/W$", "$\\gamma$", "$1-\\gamma$", "$\\alpha$", "Pauli terms", "trial RZZ"],
                  rows, "rrrrrrrrrrr", "Inputs to the filter for every case (DMRG values; $J_1=1$). Definitions in Section~\\ref{sec:glossary}.", "tab:inputs")

# DMRG vs ED
rows = []
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14):
        r = one(N, J2)
        e = r and r["inputs"].get("dmrg_vs_ed")
        if e:
            mb = max(x["max_bond"] for x in r["inputs"]["dmrg_runs"])
            rows.append([N, J2NAME[J2], sci(max(e["E0_diff"], 1e-17)), sci(max(e["E1_diff"], 1e-17)), sci(max(e["E_top_diff"], 1e-17)), sci(max(abs(1 - r["inputs"]["dmrg_ed_overlap"]), 1e-17)), mb, "yes"])
    for k, v in sorted(VAL.items(), key=lambda kv: kv[1].get("N", 0)):
        if k.startswith("ed_") and abs(v["J2"] - J2) < 1e-9:
            rows.append([v["N"], J2NAME[J2], sci(float(v["dE0"])), sci(float(v["dE1"])), sci(float(v["dEtop"])), sci(float(v["one_minus_overlap"])), max(v["max_bond"]), "yes" if v["converged"] else "NO"])
T["ed"] = tab(["$N$", "$J_2$", "$|E_0-E_0^{\\rm ED}|$", "$|E_1-E_1^{\\rm ED}|$", "$|E_{\\rm top}-E_{\\rm top}^{\\rm ED}|$", "$1-|\\langle\\psi_0|\\psi_0^{\\rm ED}\\rangle|^2$", "max bond", "conv."],
              rows, "rrrrrrrr", "DMRG against exact diagonalisation (ED). All energies in units of $J_1$.", "tab:ed")
rows = []
for k, v in sorted(VAL.items(), key=lambda kv: (kv[1].get("J2", 0), kv[1].get("N", 0))):
    if k.startswith("consist_"):
        rows.append([v["N"], J2NAME[v["J2"]], sci(float(v["dE0"])), sci(float(v["dE1"])), sci(float(v["dEtop"])), sci(float(v["one_minus_overlap"])), sci(float(v["var0"])), sci(float(v["var1"])), "yes" if v["converged"] else "NO"])
T["cons"] = tab(["$N$", "$J_2$", "$|\\Delta E_0|$", "$|\\Delta E_1|$", "$|\\Delta E_{\\rm top}|$", "$1-|\\langle\\psi_0|\\psi_0'\\rangle|^2$", "var.\\ of $\\psi_0$", "var.\\ of $\\psi_1$", "conv."],
                rows, "rrrrrrrrr", "No exact reference: two independent DMRG runs from different random starts. $\\Delta E$ is the difference between the runs; the variance is $\\langle H^2\\rangle-\\langle H\\rangle^2$ of each state.", "tab:cons")
rows = []
for k, v in sorted(VAL.items(), key=lambda kv: (kv[1].get("J2", 0), kv[1].get("N", 0))):
    if k.startswith("gamma_"):
        rows.append([v["N"], J2NAME[v["J2"]], fx(v["gamma_ed_dense"], 8), fx(v["gamma_mps_contraction"], 8), sci(abs(float(v["diff"])))])
T["gam"] = tab(["$N$", "$J_2$", "$\\gamma$ (dense, exact GS)", "$\\gamma$ (MPS contraction)", "difference"], rows, "rrrrr",
               "The trial-state weight $\\gamma$ computed by MPS contraction (used for $N>14$) against the dense overlap with the exact ground state (testing only).", "tab:gam")

# guaranteed designs small N
rows = []
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14):
        d = des(N, J2)
        if not d:
            continue
        r = one(N, J2)
        rows.append([N, J2NAME[J2], d["n_pulses_used"], num(d["n"]), fx(d["T_scaled"], 1), fx(d["T_phys"], 2), fx(d["eta"], 3), fx(d["p_succ_lb"], 3), sci(d["eps_bound"], 2), sci(d.get("eps_actual"), 2) if d.get("simulated") else "--"])
T["gd"] = tab(["$N$", "$J_2$", "$m$", "$n$", "$T$", "$T_{\\rm phys}$", "$\\eta$", "$p_g$", "bound on $1-F$", "measured $1-F$"], rows, "rrrrrrrrrr",
              "Guaranteed designs for $\\varepsilon=10^{-2}$, $N\\leq14$. $m$ = pulses, $n$ = Trotter steps, $T=\\sum t_i$ in units of $1/W$, $T_{\\rm phys}=T/W$ in units of $1/J_1$, $\\eta$ = certified suppression, $p_g$ = proven lower bound on the success probability, ``bound'' = proven upper bound on the infidelity $1-F$, ``measured'' = simulated against the exact ground state.", "tab:gd")
rows = []
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14):
        d = des(N, J2)
        if not d:
            continue
        r = one(N, J2); t = d["counts_rzz_a2a"]["total"]; f_ = d["counts_rzz_a2a"]["filter"]; ps = d["p_succ"] if d.get("simulated") else d["p_succ_lb"]
        dd = d["counts_rzz_direct"]["filter"]
        rows.append([N, J2NAME[J2], int(r["inputs"]["trial_cost_rzz"]["rzz"]), num(f_["rzz"]), num(t["rzz"]), num(dd["rzz"]), num(t["nc_total"]), num(dd["nc_total"]), num(t["depth"]), num(t["rzz"] / ps)])
T["gc"] = tab(["$N$", "$J_2$", "RZZ trial", "RZZ filter", "RZZ total", "RZZ filter (direct)", "non-Cliff.\\ rot. (level 3)", "non-Cliff.\\ rot. (direct)", "depth", "RZZ / accept"], rows, "rrrrrrrrrr",
              "Gate counts of the guaranteed designs ($\\varepsilon=10^{-2}$, all-to-all connectivity). Totals include the trial circuit. ``RZZ / accept'' = RZZ total divided by the success probability: the expected RZZ per accepted run. Non-Clifford rotations = $R_z$ plus RZZ gates whose angle is not a multiple of $\\pi/2$.", "tab:gc")

# guaranteed vs empirical
rows = []
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14):
        g, e = des(N, J2), des(N, J2, "empirical")
        if g and e:
            rows.append([N, J2NAME[J2], num(g["n"]), num(e["n"]), f"{g['n'] / e['n']:.0f}$\\times$", num(g["counts_rzz_a2a"]["total"]["rzz"]), num(e["counts_rzz_a2a"]["total"]["rzz"]), sci(e.get("eps_actual"), 2)])
T["ge"] = tab(["$N$", "$J_2$", "$n$ guaranteed", "$n$ empirical", "ratio", "RZZ guaranteed", "RZZ empirical", "measured $1-F$ (emp.)"], rows, "rrrrrrrr",
              "Guaranteed versus empirical (oracle) designs, $\\varepsilon=10^{-2}$. The empirical design is the smallest $n$ whose \\emph{measured} infidelity is below $\\varepsilon$; it uses the exact ground state.", "tab:ge")

# per step
rows = []
for J2 in J2S:
    for N in (4, 8, 12, 16, 20, 24, 28, 32):
        r = one(N, J2)
        if not r:
            continue
        a = r["per_step"]["rzz"]["a2a"]["one"]; l = r["per_step"]["rzz"]["line"]["avg3"]
        dr_ = r["per_step"]["rzz"]["direct"]
        rows.append([N, J2NAME[J2], r["inputs"]["n_pauli_terms"], int(a["rzz"]), int(a["nc_total"]), int(dr_["rzz"]), int(dr_["nc_total"]), int(a["depth"]), f"{l['rzz']:.0f}"])
T["ps"] = tab(["$N$", "$J_2$", "Pauli terms", "RZZ (level 3)", "non-Cliff.\\ rot. (level 3)", "RZZ (direct)", "non-Cliff.\\ rot. (direct)", "depth (level 3)", "RZZ (line, level 3)"], rows, "rrrrrrrrr",
              "Transpiled cost of \\emph{one} Trotter step of one pulse (every term of $H$, controlled by the ancilla). ``Line'': ancilla at the end of a linear chain of qubits (routing included; average of three transpiler seeds).", "tab:ps")

# CX -> RZZ comparison (legacy CX-basis counts are kept in the data files)
rows = []
for J2 in J2S:
    for N in (6, 10, 14, 16, 20, 24, 32):
        d = des(N, J2)
        if not d:
            continue
        c3 = d["counts_a2a"]["total"]; z3 = d["counts_rzz_a2a"]["total"]; cd = d["counts_direct"]["filter"]; zd = d["counts_rzz_direct"]["filter"]
        rows.append([N, J2NAME[J2], num(c3["cx"]), num(z3["rzz"]), f"{z3['rzz'] / c3['cx']:.2f}", num(cd["cx"]), num(zd["rzz"]), f"{zd['rzz'] / cd['cx']:.2f}", num(c3["rz_nonclifford"]), num(z3["nc_total"])])
T["cxrzz"] = tab(["$N$", "$J_2$", "CX (level 3)", "RZZ (level 3)", "ratio", "CX (direct)", "RZZ (direct)", "ratio", "non-Cliff.\\ $R_z$ (CX, level 3)", "non-Cliff.\\ rot. (RZZ, level 3)"], rows, "rrrrrrrrrr",
                 "Previous CX basis versus the RZZ basis, guaranteed designs at $\\varepsilon=10^{-2}$. The direct compilation maps each CX onto one Clifford RZZ($\\pi/2$) (ratio 1.00); level 3 needs only two RZZ per ancilla-controlled Pauli rotation, one of them a non-Clifford RZZ($\\theta$).", "tab:cxrzz")

# large N
rows = []
for J2 in J2S:
    for N in (14, 16, 18, 20, 24, 28, 32):
        r = one(N, J2)
        if not r:
            continue
        d = des(N, J2)
        i = r["inputs"]
        if not d:
            rows.append([N, J2NAME[J2], fx(i["gamma"], 3), fx(i["gap_scaled"], 4), "--", "not achievable", "--", "--", "--", "--", "--"])
            continue
        t = d["counts_rzz_a2a"]["total"]
        flag = "" if r["design"]["source"] in ("floor", "none") else "$^\\ddagger$"
        rows.append([N, J2NAME[J2], fx(i["gamma"], 3), fx(i["gap_scaled"], 4), d["n_pulses_used"], num(d["n"]) + flag, fx(d["T_phys"], 1), sci(t["rzz"], 2), sci(d["counts_rzz_direct"]["filter"]["rzz"], 2), sci(d["counts_rzz_direct"]["filter"]["nc_total"], 2), fx(d["p_succ_lb"], 3)])
T["large"] = tab(["$N$", "$J_2$", "$\\gamma$", "$\\Delta$", "$m$", "$n$", "$T_{\\rm phys}$", "RZZ total (level 3)", "RZZ (direct)", "non-Cliff.\\ rot. (direct)", "$p_g$"], rows, "rrrrrrrrrrr",
                 "Guaranteed designs for $N\\geq14$, $\\varepsilon=10^{-2}$ (DMRG/MPS only; circuits counted, not simulated). ``Not achievable'': no filter with at most 8 pulses meets the target for this trial state (Fig.~\\ref{fig:gamma}). $^\\ddagger$: fallback design.", "tab:large")
rows = []
for J2 in J2S:
    for N in (14, 16, 18, 20, 24, 28, 32):
        r = one(N, J2, 0.1)
        d = des(N, J2, "guaranteed", 0.1)
        if r and d:
            t = d["counts_rzz_a2a"]["total"]
            rows.append([N, J2NAME[J2], d["n_pulses_used"], num(d["n"]), fx(d["T_phys"], 1), sci(t["rzz"], 2), sci(d["counts_rzz_direct"]["filter"]["rzz"], 2), sci(d["counts_rzz_direct"]["filter"]["nc_total"], 2), fx(d["p_succ_lb"], 3)])
T["large1"] = tab(["$N$", "$J_2$", "$m$", "$n$", "$T_{\\rm phys}$", "RZZ total (level 3)", "RZZ (direct)", "non-Cliff.\\ rot. (direct)", "$p_g$"], rows, "rrrrrrrrr",
                  "Guaranteed designs for $N\\geq14$ at the looser target $\\varepsilon=10^{-1}$.", "tab:large1")

# eps scan
rows = []
for (N, J2) in ((6, 0.0), (8, 0.0), (6, 0.4)):
    for e in sorted({e for (n_, j_, e) in R if n_ == N and abs(j_ - J2) < 1e-9}, reverse=True):
        for which in ("guaranteed", "empirical"):
            d = des(N, J2, which, e)
            if d:
                rows.append([N, J2NAME[J2], f"{e:g}", which, num(d["n"]), d["n_pulses_used"], fx(d["T_phys"], 2), sci(d["counts_rzz_a2a"]["total"]["rzz"], 2), sci(d.get("eps_actual"), 2) if d.get("simulated") else "--"])
T["eps"] = tab(["$N$", "$J_2$", "$\\varepsilon$", "design", "$n$", "$m$", "$T_{\\rm phys}$", "RZZ total", "measured $1-F$"], rows, "rrrlrrrrr",
               "Dependence on the target $\\varepsilon$ for three small cases.", "tab:eps")

# bandwidth
rows = []
for J2 in J2S:
    for N in (4, 8, 12, 16, 20, 24, 28, 32):
        r = one(N, J2)
        if r:
            i = r["inputs"]
            rows.append([N, J2NAME[J2], fx(i["W"], 5), sci(abs(i["W"] - i["W_exact_dmrg"]), 1), fx(3 * (i["W"] + i["E0"]) - i["E0"], 4), f"{(3 * (i['W'] + i['E0']) - i['E0']) / i['W']:.2f}"])
T["bw"] = tab(["$N$", "$J_2$", "$W$ (Weyl)", "$|W-(E_{\\rm top}-E_0)_{\\rm DMRG}|$", "$W$ (Pauli-$\\ell_1$, earlier)", "$\\ell_1$ / Weyl"], rows, "rrrrrr",
              "The certified bandwidth $W=\\sum_{\\rm bonds}\\lambda_{\\max}-E_0$ against the DMRG value of $E_{\\rm top}-E_0$ and against the looser Pauli bound used earlier.", "tab:bw")

# validation + timing
vals = [(k, r["validation_rzz"]) for k, r in R.items() if r.get("validation_rzz") and "l3" in r["validation_rzz"]]
vz3 = max(abs(v["l3"]["est_rzz"] - v["l3"]["full_rzz"]) for _, v in vals); vzd = max(abs(v["direct"]["est_rzz"] - v["direct"]["full_rzz"]) for _, v in vals)
vnd = max(abs(v["direct"]["est_nc"] - v["direct"]["full_nc"]) for _, v in vals); vo3 = 100 * max((v["l3"]["est_nc"] - v["l3"]["full_nc"]) / v["l3"]["full_nc"] for _, v in vals)
dpe = [abs(v["l3"]["est_depth"] - v["l3"]["full_depth"]) / v["l3"]["full_depth"] * 100 for _, v in vals]
timing = {}
for (N, J2, e), r in R.items():
    timing.setdefault(N, []).append((r["inputs"]["build_s"], r["design"].get("design_s", 0.0), r["wall_s"]))
rows = []
for N in sorted(timing):
    a = np.array(timing[N])
    rows.append([N, len(a), f"{np.median(a[:, 0]):.0f}", f"{np.median(a[:, 1]):.0f}", f"{np.median(a[:, 2]):.0f}", f"{a[:, 2].max():.0f}"])
T["time"] = tab(["$N$", "points", "context build (median s)", "filter design (median s)", "total per point (median s)", "total per point (max s)"], rows, "rrrrrr",
                "Classical compute time per point (single core). The context build includes the three DMRG runs, the trial circuit and its overlap, $\\alpha$ and the per-step transpile.", "tab:time")

# not achievable
na = []
for (N, J2, e), r in sorted(R.items()):
    d = r["designs"].get("guaranteed")
    if d is not None and "n" not in d:
        na.append((N, J2, e, r["inputs"]["gamma"]))
nat = ""
if na:
    nat = "\\begin{itemize}\n" + "\n".join(f"\\item $N={N}$, $J_2={J2NAME[round(J2, 4)]}$, $\\varepsilon={e:g}$: $\\gamma={g:.3f}$" for N, J2, e, g in na) + "\n\\end{itemize}\n"


# summary numbers
def rng(vals_):
    return min(vals_), max(vals_)


smallg = [des(N, J2) for J2 in J2S for N in (4, 6, 8, 10, 12, 14) if des(N, J2)]
bigg = [(N, J2, des(N, J2)) for J2 in J2S for N in (16, 18, 20, 24, 28, 32) if des(N, J2) and one(N, J2)["design"]["source"] in ("floor", "none")]
sm_n = rng([d["n"] for d in smallg]); sm_cx = rng([d["counts_rzz_a2a"]["total"]["rzz"] for d in smallg])
bg_n = rng([d["n"] for _, _, d in bigg]); bg_cx = rng([d["counts_rzz_a2a"]["total"]["rzz"] for _, _, d in bigg]); bg_rz = rng([d["counts_rzz_direct"]["filter"]["nc_total"] for _, _, d in bigg]); bg_nc3 = rng([d["counts_rzz_a2a"]["total"]["nc_total"] for _, _, d in bigg]); bg_cxd = rng([d["counts_rzz_direct"]["filter"]["rzz"] for _, _, d in bigg])
bg_T = rng([d["T_phys"] for _, _, d in bigg])
ed_err = max(max(float(v["dE0"]), float(v["dE1"]), float(v["dEtop"])) for k, v in VAL.items() if k.startswith("ed_"))
for J2 in J2S:
    for N in (4, 6, 8, 10, 12, 14):
        r = one(N, J2); e = r and r["inputs"].get("dmrg_vs_ed")
        if e:
            ed_err = max(ed_err, e["E0_diff"], e["E1_diff"], e["E_top_diff"])
cons_err = max(max(float(v["dE0"]), float(v["dE1"]), float(v["dEtop"])) for k, v in VAL.items() if k.startswith("consist_"))
ratios = [g["n"] / e["n"] for J2 in J2S for N in (4, 6, 8, 10, 12, 14) for g, e in [(des(N, J2), des(N, J2, "empirical"))] if g and e]
nviol = sum(1 for r in R.values() for d in r["designs"].values() if d and d.get("simulated") and d.get("eps_actual") is not None and d["eps_actual"] > d["eps_bound"] + 1e-9)
fits_txt = ", ".join(f"$N^{{{FITTXT[('cx', J2, '-')][0]:.1f}}}$ ($J_2={J2NAME[J2]}$, $N={FITTXT[('cx', J2, '-')][1]}$--${FITTXT[('cx', J2, '-')][2]}$)" for J2 in J2S if ('cx', J2, '-') in FITTXT)
cmp_txt = (f"Switching from the earlier Pauli-$\\ell_1$ bandwidth to the Weyl bandwidth changed the guaranteed step count by a factor {BWC['n_ratio_min']:.2f}--{BWC['n_ratio_max']:.2f} "
           f"(median {BWC['n_ratio_median']:.2f}) and the physical time by {BWC['T_ratio_min']:.2f}--{BWC['T_ratio_max']:.2f}, over {BWC['count']} matching designs.") if BWC else ""

# ------------------------------------------------------------------------------------------------ document
GLOSS = [
    ("$N$", "number of spin-$\\tfrac12$ sites in the open chain."),
    ("$J_1,J_2$", "nearest- and next-nearest-neighbour exchange couplings; $J_1=1$ sets the energy unit. $J_2=0.2411$ is the point where the chain's ground state changes from critical to dimerised."),
    ("$H$", "Hamiltonian $H=J_1\\sum_i\\mathbf S_i\\!\\cdot\\!\\mathbf S_{i+1}+J_2\\sum_i\\mathbf S_i\\!\\cdot\\!\\mathbf S_{i+2}$, written as a sum of Pauli strings (``Pauli terms'')."),
    ("$E_0,\\,E_1,\\,E_{\\rm top}$", "lowest, second-lowest (first excited state of any symmetry) and highest energy of $H$. $g=E_1-E_0$ is the energy gap."),
    ("DMRG", "density-matrix renormalisation group: a classical algorithm that finds low-energy states of 1D systems as matrix-product states (MPS). ``max bond'' is the largest MPS bond dimension reached."),
    ("$W$", "bandwidth, a rigorous upper bound on $E_{\\rm top}-E_0$ obtained from Weyl's inequality (Section~\\ref{sec:bandwidth}); here it equals $E_{\\rm top}-E_0$."),
    ("$H_s,\\ \\Delta$", "scaled Hamiltonian $H_s=(H-E_0)/W$, whose ground state sits at $0$ and whose other levels lie in $[\\Delta,1]$ with scaled gap $\\Delta=g/W$."),
    ("trial state, $\\gamma$", "the starting state prepared by a shallow circuit compiled from the DMRG state; $\\gamma$ is its weight on the true ground state ($1-\\gamma$ is the part the filter must remove)."),
    ("filter $F(E)$", "$F(E)=\\prod_{i=1}^m\\cos(Et_i+\\phi_i)$ with $m$ pulses of times $t_i$ and phases $\\phi_i$; applied to the state it multiplies the weight of the energy-$E$ component by $F(E)$."),
    ("pulse, ancilla", "one pulse is a Hadamard, a rotation, a controlled evolution $e^{-it_iH_s\\otimes Z}$ and a Hadamard on one extra (ancilla) qubit; reading the ancilla as $0$ applies $\\cos(H_st_i+\\phi_i)$ to the system. The run is kept only if all ancilla readings are $0$ (post-selection)."),
    ("$T,\\ T_{\\rm phys}$", "total evolution time $T=\\sum_it_i$ in units of $1/W$; $T_{\\rm phys}=T/W$ in units of $1/J_1$."),
    ("$\\eta$", "certified ratio of the largest value of $|F|$ on $[\\Delta,1]$ to $|F(0)|$: how strongly the filter suppresses excited states relative to the ground state."),
    ("$p_g$, $P_{\\rm succ}$", "$p_g=\\gamma F(0)^2$ is a proven lower bound on the success (post-selection) probability $P_{\\rm succ}$; the expected number of repetitions is at most $1/p_g$."),
    ("$\\varepsilon$, $1-F$", "target error: $1-F$ is the infidelity ($1$ minus the squared overlap) of the final state with the exact ground state; the design must prove $1-F\\leq\\varepsilon$."),
    ("Trotter steps $n$", "the controlled evolution is approximated by $n$ (in total, over all pulses) first-order product-formula steps; $k_i$ steps in pulse $i$, $t_i=k_i\\,{\\rm d}t$."),
    ("$\\alpha$", "commutator constant $\\alpha=\\sum_g\\|[H_g,\\sum_{g'>g}H_{g'}]\\|$ of the scaled Hamiltonian: it sets the Trotter error $\\alpha T^2/(2n)$."),
    ("guaranteed design", "smallest $n$ for which the \\emph{proven} bound on $1-F$ is $\\leq\\varepsilon$ (needs no ground state)."),
    ("empirical (oracle) design", "smallest $n$ for which the \\emph{measured} $1-F$ is $\\leq\\varepsilon$; requires simulation and the exact ground state, so it is a reference, only available for $N\\leq14$."),
    ("RZZ, $R_z$, non-Clifford rotations, depth", "two-qubit gates ${\\rm RZZ}(\\theta)=e^{-i\\theta ZZ/2}$ (CX is not in the gate set; ${\\rm RZZ}(\\pi/2)$ is Clifford), single-qubit $z$-rotations, and circuit depth after transpiling to the basis $\\{{\\rm RZZ},R_z,\\sqrt{X},X\\}$. ``Non-Clifford rotations'' = $R_z$ plus RZZ gates whose angle is not a multiple of $\\pi/2$: the expensive gates on fault-tolerant hardware. Clifford single-qubit gates are not tracked."),
    ("direct / level 3", "two compilations of the same circuit. \\emph{Direct} (qiskit optimisation level 1): every ancilla-controlled Pauli rotation is a ladder of four Clifford RZZ($\\pi/2$) around one $R_z$, giving $4\\times$(Pauli terms) RZZ and (Pauli terms $+1$) non-Clifford rotations per Trotter step. \\emph{Level 3}: blocks are re-synthesised, giving two RZZ per Pauli term (one of them a non-Clifford RZZ($\\theta$)), i.e. half the RZZ count, plus some extra non-Clifford $R_z$."),
    ("all-to-all / line", "qubit connectivity assumed in the transpile: every qubit couples to every other, or the qubits sit on a line (extra RZZ from routing)."),
    ("conditional certificate", "the proven bound holds if the DMRG inputs ($E_0$, $E_1$, $\\gamma$) are exact; they cannot be certified by DMRG itself (see Section~\\ref{sec:limits})."),
]
gl_rows = [[a, b] for a, b in GLOSS]
glossary = ("\\begin{longtable}{p{0.2\\linewidth}p{0.74\\linewidth}}\n\\toprule\nTerm & Meaning \\\\\n\\midrule\n\\endhead\n" +
            "\n".join(f"{a} & {b} \\\\" for a, b in gl_rows) + "\n\\bottomrule\n\\end{longtable}\n")

doc = r"""\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{booktabs,longtable,array,graphicx,float}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\setlength{\LTcapwidth}{\linewidth}
\title{Resource Costs of the Cosine-Filter Pipeline for the Open $J_1$--$J_2$ Chain\\[2pt]\large Results, with every quantity defined}
\author{}
\date{}
\begin{document}
\maketitle

\section{What this report shows}
A shallow circuit built from a DMRG state leaves some excited-state weight. A \emph{cosine filter} (the Stetcu--Baroni--Carlson projection) removes it using one extra qubit, controlled time evolution, and post-selection. This report gives the \textbf{gate cost of that filter} for the spin-$\tfrac12$ $J_1$--$J_2$ chain with $N$ up to """ + str(max(N for (N, _, _) in R)) + r""", at the target error $\varepsilon=10^{-2}$ (and $10^{-1}$), together with the accuracy of the classical DMRG inputs and a check of every proven bound. Every symbol is defined in Section~\ref{sec:glossary}. The error guarantees are \emph{conditional certificates}: they hold if the DMRG energies and the trial-state weight are exact, which DMRG cannot itself certify (Section~\ref{sec:limits}); the derivation of every chosen value is in \texttt{CERTIFICATION.tex}.

\section{Headline results}
\begin{itemize}
\item \textbf{DMRG inputs are accurate.} Energies agree with exact diagonalisation to $""" + f"{ed_err:.1e}" + r"""$ ($N\leq20$); independent DMRG runs agree to $""" + f"{cons_err:.1e}" + r"""$ for $N=24$--$32$.
\item \textbf{Small chains ($N\leq14$, $\varepsilon=10^{-2}$).} Guaranteed designs need """ + num(sm_n[0]) + "--" + num(sm_n[1]) + r""" Trotter steps and """ + num(sm_cx[0]) + "--" + num(sm_cx[1]) + r""" RZZ gates; the measured error stays below the proven bound in all """ + str(sum(1 for r in R.values() for d in r["designs"].values() if d and d.get("simulated"))) + r""" simulated designs (violations: """ + str(nviol) + r"""). The guaranteed design needs """ + f"{min(ratios):.0f}--{max(ratios):.0f}" + r"""$\times$ the steps of the oracle design.
\item \textbf{Larger chains ($N=16$--$32$, $\varepsilon=10^{-2}$).} """ + num(bg_n[0]) + "--" + num(bg_n[1]) + r""" steps, """ + sci(bg_cx[0], 1) + "--" + sci(bg_cx[1], 1) + r""" RZZ with block re-synthesis (""" + sci(bg_cxd[0], 1) + "--" + sci(bg_cxd[1], 1) + r""" RZZ with direct compilation) and """ + sci(bg_rz[0], 1) + "--" + sci(bg_rz[1], 1) + r""" non-Clifford rotations (direct compilation; """ + sci(bg_nc3[0], 1) + "--" + sci(bg_nc3[1], 1) + r""" with re-synthesis), with physical evolution time """ + f"{bg_T[0]:.1f}--{bg_T[1]:.1f}" + r""" (units of $1/J_1$). Total RZZ grows like """ + fits_txt + r""".
\item \textbf{The bandwidth does not change the cost.} """ + cmp_txt + r"""
\item \textbf{Limit.} For $N=28,32$ at $J_2=0$ and $0.2411$ the shallow trial circuit has too little ground-state weight, and no certified filter with at most 8 pulses reaches $\varepsilon=10^{-2}$ (Fig.~\ref{fig:gamma}).
\end{itemize}

\section{Glossary: every term used}\label{sec:glossary}
""" + glossary + r"""

\section{How the numbers are produced}\label{sec:method}
For each chain $(N,J_2)$ the pipeline does the following; all steps use only DMRG and circuits (no ground-state knowledge) except where marked.
\begin{enumerate}
\item \textbf{DMRG} finds $E_0$ and the ground-state MPS; a second DMRG run finds $E_1$ as the lowest state orthogonal to it, whatever its symmetry (the symmetry of the first excitation is not assumed: exact diagonalization finds a triplet for every $J_2\in[0,1]$ on the grid, open chains $N=8$--$16$).
\item \textbf{Scale} $H\to H_s=(H-E_0)/W$ so that the spectrum lies in $\{0\}\cup[\Delta,1]$ (Section~\ref{sec:bandwidth}).
\item \textbf{Trial circuit}: the MPS is compiled to a shallow brickwork circuit ($L=1$ layer); $\gamma$ is its overlap with the DMRG ground state (computed from an MPS simulation of the circuit for $N>14$).
\item \textbf{Design the filter}: choose $m\leq8$ pulses, times $t_i$ and phases $\phi_i$ that maximise $F(0)^2$ subject to a \emph{certified} bound $\eta\leq\eta_{\rm target}(\varepsilon,\gamma)$.
\item \textbf{Choose $n$}: the smallest number of Trotter steps for which the proven bound $1-F\leq\big(\sqrt{\bar\ell}+\alpha T^2/(2n\sqrt{p_g})\big)^2\leq\varepsilon$ holds, where $\bar\ell$ is the leakage of the exact filter. (\emph{Empirical} design: simulate and take the smallest $n$ that works.)
\item \textbf{Count}: transpile one Trotter step once, multiply by $n$, add the trial circuit; checked against a full transpile (Section~\ref{sec:checks}).
\end{enumerate}

\section{Case inputs}\label{sec:inputs}
Table~\ref{tab:inputs} lists, for each chain, the DMRG energies and the quantities derived from them. The gap shrinks like $1/N$, so $W/g$ and hence the filter time grow with $N$ (Fig.~\ref{fig:gap}). The trial-state weight $\gamma$ falls with $N$ for fixed circuit depth, most strongly at small $J_2$ (Fig.~\ref{fig:gamma}); $1-\gamma$ is the infidelity of the trial \emph{circuit}, not an error of DMRG.
""" + T["inputs"] + figure("gamma_vs_N.png", "Ground-state weight $\\gamma$ of the $L=1$ trial circuit against $N$. Crosses mark cases where no certified filter reaches $\\varepsilon=10^{-2}$.", "fig:gamma", "0.75\\linewidth") + figure("gap_and_time.png", "(a) Energy gap $g=E_1-E_0$ against $N$. (b) Physical evolution time $T_{\\rm phys}$ of the guaranteed design ($\\varepsilon=10^{-2}$): it grows because the filter must resolve the shrinking gap ($T_{\\rm phys}=x\\pi/g$). The jumps reflect the discrete choice of $x=T\\Delta/\\pi\\in[0.4,2.5]$ made by the design search, not the physics.", "fig:gap") + r"""

\section{Accuracy of the DMRG inputs}\label{sec:dmrg}
The filter design depends on $E_0$, $E_1$ (through $\Delta$) and $\gamma$, all computed by DMRG. Settings: SVD cutoff $10^{-14}$, energy tolerance $10^{-12}$, inner eigensolver tolerance $10^{-10}$, bond dimension up to 256. The excited state is found by adding the rank-one penalty $\lambda|\psi_0\rangle\langle\psi_0|$ to each local two-site problem, which returns the lowest state of any symmetry orthogonal to the ground state without building the large projector operator. Table~\ref{tab:ed} compares with exact diagonalisation (ED, used for testing only; $N\leq20$); Table~\ref{tab:cons} compares two independent runs where no exact answer exists ($N=24$--$32$); Table~\ref{tab:gam} checks the MPS computation of $\gamma$ used above $N=14$. Fig.~\ref{fig:dmrg} summarises.
""" + figure("dmrg_accuracy.png", "DMRG energy errors. Small markers: against exact diagonalisation. Large open markers ($N\\geq24$): largest run-to-run difference of $E_0,E_1$. Dotted line: $10^{-13}$.", "fig:dmrg") + T["ed"] + T["cons"] + T["gam"] + r"""
Run-to-run agreement shows that the solves are converged and reproducible; it does not prove they sit at the true minimum. The weight $\gamma$ is the least accurate input (first order in the state error): it agrees with the exact overlap to $""" + f"{max(abs(float(v['diff'])) for k, v in VAL.items() if k.startswith('gamma_')):.1e}" + r"""$.

\section{The certified bandwidth}\label{sec:bandwidth}
The scaling $H_s=(H-E_0)/W$ requires $W\geq E_{\rm top}-E_0$. DMRG on $-H$ is variational from the wrong side ($E_{\rm top}^{\rm DMRG}\leq E_{\rm top}$), so it cannot supply a safe $W$. Instead we use \textbf{Weyl's inequality}: the largest eigenvalue of a sum is at most the sum of the largest eigenvalues. A bond $J\,\mathbf S\!\cdot\!\mathbf S$ has eigenvalues $J/4$ (triplet) and $-3J/4$ (singlet), so for $J_1,J_2\geq0$
\[
E_{\rm top}\;\leq\;\tfrac{N-1}{4}J_1+\tfrac{N-2}{4}J_2,\qquad W=\tfrac{N-1}{4}J_1+\tfrac{N-2}{4}J_2-E_0 .
\]
The fully polarised state attains the bound, so $W=E_{\rm top}-E_0$ exactly (Table~\ref{tab:bw}). The Pauli-term bound used earlier ($\sum|c_\beta|-E_0$) is $1.6$--$1.9\times$ larger. \textbf{$W$ cancels in the cost:} $\alpha_{\rm scaled}=\alpha/W^2$ and $T=x\pi W/g$, so $n\propto\alpha T^2\propto\alpha/g^2$ and $T_{\rm phys}=x\pi/g$ do not depend on $W$; a looser $W$ only widens the suppression window. """ + cmp_txt + r"""
""" + T["bw"] + r"""

\section{Small chains: guaranteed designs}\label{sec:small}
Table~\ref{tab:gd} gives the design (pulses, steps, times) and the proven and measured infidelity for $N\leq14$ at $\varepsilon=10^{-2}$; Table~\ref{tab:gc} the gate counts. The proven bound is below $\varepsilon$ by construction, and the measured infidelity is smaller by a factor of about 6--40 (Fig.~\ref{fig:bound}), because the Trotter part of the bound is a worst case over all states.
""" + T["gd"] + T["gc"] + figure("bound_vs_measured.png", "Measured infidelity $1-F$ of the guaranteed designs ($N\\leq14$, simulated against the exact ground state) compared with the proven bound (dashed), which equals the target $\\varepsilon=10^{-2}$. The measured values lie well below the bound: it is never violated.", "fig:bound", "0.55\\linewidth") + r"""

\section{Guaranteed versus empirical designs}\label{sec:emp}
The empirical design is the smallest $n$ whose \emph{measured} infidelity meets $\varepsilon$; it needs the exact ground state and a simulation, so it shows what an ideal user could save. Table~\ref{tab:ge} and Fig.~\ref{fig:steps} compare. The proven design needs a large multiple of the empirical step count (and of the RZZ count, which is proportional to $n$) because its Trotter error bound is pessimistic.
""" + T["ge"] + figure("steps_vs_N.png", "Trotter steps $n$ against $N$ for the guaranteed (solid) and empirical (dashed, $N\\leq14$) designs, $\\varepsilon=10^{-2}$.", "fig:steps") + r"""

\section{Cost of one Trotter step}\label{sec:step}
All circuit costs are $n\times$ (cost of one step) plus the trial circuit. One step applies every term of $H$ once, controlled by the ancilla (Table~\ref{tab:ps}, Fig.~\ref{fig:step}). RZZ per step grows linearly in $N$ (one block per bond); $J_2\neq0$ roughly doubles it because the number of bonds doubles.
""" + T["ps"] + T["cxrzz"] + figure("rzz_per_step.png", "RZZ per Trotter step against $N$, all-to-all connectivity (solid) and line connectivity (dotted).", "fig:step", "0.7\\linewidth") + r"""

\section{Larger chains, $N$ up to """ + str(max(N for (N, _, _) in R)) + r"""}\label{sec:large}
Above $N=14$ nothing is simulated as a state vector: DMRG gives the energies and MPS, the trial circuit is simulated as an MPS to obtain $\gamma$, the filter is designed from $(\gamma,\Delta,\alpha)$, and circuits are counted. Tables~\ref{tab:large} and~\ref{tab:large1} give the guaranteed designs at $\varepsilon=10^{-2}$ and $10^{-1}$; Fig.~\ref{fig:cost} shows how the total RZZ and non-Clifford rotation counts grow (fits: """ + fits_txt + r"""). For $N=28,32$ at $J_2=0$ and $0.2411$ the trial circuit has $\gamma=0.31$--$0.51$ and the required suppression ($\eta_{\rm target}\approx0.07$--$0.10$) cannot be certified with at most 8 pulses; the looser target $\varepsilon=10^{-1}$ is feasible. This is a limit of the one-layer trial circuit, not of the filter.
""" + T["large"] + T["large1"] + figure("cost_vs_N.png", "(a) Total RZZ and (b) total non-Clifford rotations of the guaranteed designs ($\\varepsilon=10^{-2}$), log--log, with power-law fits (lines).", "fig:cost") + r"""

\section{Dependence on the target $\varepsilon$}\label{sec:eps}
Tightening $\varepsilon$ needs both a stronger filter (smaller $\eta$) and fewer Trotter errors (more steps). Table~\ref{tab:eps} and Fig.~\ref{fig:eps} show three small cases.
""" + T["eps"] + figure("eps_scan.png", "Trotter steps against target $\\varepsilon$ for guaranteed (solid) and empirical (dashed) designs.", "fig:eps") + r"""

\section{Checks and classical compute cost}\label{sec:checks}
\textbf{Gate-count method.} The per-step transpile times $n$ was compared with a full transpile of the real circuit in the RZZ basis at $n=48$ (4 pulses of 12 steps, all-to-all) for """ + str(len(vals)) + r""" cases: the RZZ count differed by at most """ + f"{vz3:.0f}" + r""" (level 3) and """ + f"{vzd:.0f}" + r""" (direct), the direct non-Clifford count by at most """ + f"{vnd:.0f}" + r""", and the level-3 non-Clifford estimate was at most """ + f"{vo3:.1f}" + r"""\% above the full transpile; the level-3 depth estimate differed by at most """ + f"{max(dpe):.1f}" + r"""\%. \textbf{Bound violations:} """ + str(nviol) + r""" among the simulated designs. \textbf{Classical cost:} Table~\ref{tab:time}.
""" + T["time"] + (r"""
\textbf{Points without a certified design:}
""" + nat if na else "") + r"""

\section{Limits and assumptions}\label{sec:limits}
\begin{itemize}
\item \textbf{Conditional certificates.} The proven bounds assume the DMRG shift $E_0$, the gap and the weight $\gamma$ are exact. Interlacing shows the DMRG gap can overestimate the true gap only by the solver's own optimisation error, and DMRG's $E_0$ errs on the safe side, but these errors are not computable from DMRG; validated accuracy is $\sim10^{-13}$ (Section~\ref{sec:dmrg}). The derivation, including how the bound changes for given error bounds, is in \texttt{CERTIFICATION.tex}.
\item \textbf{Costs are ideal gate counts.} Ancilla measurement and reset time, classical feed-forward, native controlled-Pauli-rotation implementation on real hardware, and error-correction or mitigation overhead are not included; RZZ counts come from a heuristic compiler and are not minimal; the expected-repetition figure $1/p_g$ assumes independent attempts. No $T$-count is reported (it depends on synthesis constants); the metric is the number of non-Clifford rotations ($R_z$ and RZZ($\\theta$) with angle not a multiple of $\\pi/2$).
\item \textbf{Numerical certificate.} The bound on $\eta$ is computed in double precision without interval arithmetic (rounding $\lesssim10^{-14}$, far below the pruning tolerance).
\item \textbf{Scope.} One model family, three values of $J_2$; a non-degenerate singlet ground state; first-order Trotter formula.
\end{itemize}

\section{Reproducing}
\begin{verbatim}
export OMP_NUM_THREADS=1
for i in 0 1 2; do python scripts/run_results.py --worker $i --of 3 & done; wait   # resumable
python scripts/validate_dmrg.py ed 16 18 20; python scripts/validate_dmrg.py consist 24 28 32
python scripts/compare_bandwidth.py OLD_DATA_DIR      # optional
python scripts/audit_results.py [--alpha]             # independent re-derivation of every stored result
python scripts/build_results.py                       # RESULTS.md
python scripts/build_results_tex.py                   # RESULTS.tex and figures_tex/
\end{verbatim}
\end{document}
"""
(RES / "RESULTS.tex").write_text(doc)
print("wrote", RES / "RESULTS.tex", len(doc), "chars;", len(list(FIG.glob('*.png'))), "figures")
