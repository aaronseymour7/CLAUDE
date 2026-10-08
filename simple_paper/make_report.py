"""make_report.py -- turn results/*.json into figures, tables and statistics (results/tables.json, stats.json)."""
import json, glob, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = pathlib.Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figs"
M = {}
for f in sorted(glob.glob(str(RES / "main_N*_J2*.json"))):
    d = json.load(open(f)); M[(d["inputs"]["N"], d["inputs"]["J2"])] = d
keys = sorted(M, key=lambda k: (k[1], k[0]))
load = lambda n: json.load(open(RES / n)) if (RES / n).exists() else None
TR, SE, DE, BA = load("extra_trials.json"), load("extra_sens.json"), load("extra_design.json"), load("extra_baseline.json")
T, S = {}, {}
sci = lambda x: f"{x:.1e}".replace("e-0", "e-").replace("e+0", "e+")
cnt = lambda x: f"{int(round(x)):,}"
mr = lambda xs: (float(min(xs)), float(np.median(xs)), float(max(xs)))


def case(k, eps):
    for c in M[k]["cases"]:
        if abs(c["eps"] - eps) < 1e-12:
            return c


def best_direct(N, J2, eps):
    """cheapest direct preparation (mps-to-circuit exact / approximate, generic) with infidelity <= eps."""
    if not BA:
        return None, None
    for r in BA:
        if r["N"] == N and r["J2"] == J2:
            cands = [(e["cx"], f"sequential MPS, $\\chi={e['chi']}$") for e in r["exact"] if "cx" in e and e["infid"] <= eps]
            cands += [(a["cx"], f"approx. circuit, $L={a['L']}$") for a in r["approx"] if a["infid"] <= eps]
            if r["generic"]["cx"]:
                cands.append((r["generic"]["cx"], "generic"))
            return min(cands)
    return None, None


# ------------------------------------------------------------------ front end: DMRG inputs vs exact diagonalisation
rows = ["| $N$ | $J_2$ | $\\chi_{\\mathrm{ref}}$ | $|E_0^{\\mathrm{DMRG}}-E_0|$ | $|E_1^{\\mathrm{DMRG}}-E_1|$ | $|\\Delta^{\\mathrm{DMRG}}-\\Delta|/\\Delta$ | $\\gamma$ (DMRG-ref.) | $\\gamma$ (exact) | trial CX | $\\alpha$ |", "|---|---|---|---|---|---|---|---|---|---|"]
S["input_err"] = {}
for k in keys:
    i = M[k]["inputs"]
    e0, e1 = abs(i["E0_dmrg"] - i["E0_ed"]), abs(i["E1_dmrg"] - i["E1_ed"])
    dg = abs(i["gap_dmrg"] - i["gap_true"]) / i["gap_true"]
    S["input_err"][f"N{k[0]}_J{k[1]}"] = dict(E0=e0, E1=e1, gap_rel=dg, gamma=abs(i["gamma_est"] - i["gamma_true"]))
    rows.append(f"| {k[0]} | {k[1]} | {i['chi_ref']} | {sci(e0)} | {sci(e1)} | {sci(dg)} | {i['gamma_est']:.4f} | {i['gamma_true']:.4f} | {i['trial_cx']} | {i['alpha']:.3f} |")
T["inputs"] = "\n".join(rows)
S["max_E_err"] = max(max(v["E0"], v["E1"]) for v in S["input_err"].values())
S["max_gap_rel"] = max(v["gap_rel"] for v in S["input_err"].values())
S["max_gamma_err"] = max(v["gamma"] for v in S["input_err"].values())

# ------------------------------------------------------------------ certificate accuracy and step-count rules (eps = 1e-2)
rows = ["| $N$ | $J_2$ | $\\gamma$ | $\\Delta$ | $n_{\\mathrm{bound}}$ | $n_{\\mathrm{hyb}}$ | $n_{\\mathrm{meas}}$ | bound $\\bar D$ | $D$ at $n_{\\mathrm{bound}}$ | leakage: bound / meas. | Trotter: bound / meas. |", "|---|---|---|---|---|---|---|---|---|---|---|"]
for k in keys:
    i, c = M[k]["inputs"], case(k, 0.01)
    if c is None:
        continue
    if "n_bound" not in c or c.get("n_bound") is None:
        rows.append(f"| {k[0]} | {k[1]} | {i['gamma_est']:.3f} | {i['gap_dmrg']:.3f} | – | – | – | – | {c['D_trial']:.3f}$^\\dagger$ | – | – |")
        continue
    a = c.get("at_n_bound")
    if a:
        rows.append(f"| {k[0]} | {k[1]} | {i['gamma_est']:.3f} | {i['gap_dmrg']:.3f} | {cnt(c['n_bound'])} | {cnt(c['n_hyb']) if c.get('n_hyb') else '–'} | {c['n_meas'] if c.get('n_meas') else '–'} | {a['bound']:.3f} | {a['meas']:.3f} | {a['bound_leak']/a['meas_leak']:.1f}$\\times$ | {a['bound_trot']/a['meas_trot']:.0f}$\\times$ |")
    else:
        rows.append(f"| {k[0]} | {k[1]} | {i['gamma_est']:.3f} | {i['gap_dmrg']:.3f} | {cnt(c['n_bound'])} | {c['n_hyb'] if c.get('n_hyb') else '–'} | {c['n_meas'] if c.get('n_meas') else '–'} | {c['bound_point']['bound']:.3f} | not simulated | – | – |")
T["main"] = "\n".join(rows)

pts = []
for k in keys:
    for c in M[k]["cases"]:
        if c.get("at_n_bound"):
            pts.append(c["at_n_bound"])
        pts += c.get("sweep", [])
S["validity"] = dict(points=len(pts), **{n: sum(bool(p[n]) for p in pts) for n in ("bound_ok", "leak_ok", "trot_ok", "triangle_ok")},
                     energy_ok=sum(p["energy_err"] <= p["energy_err_bound"] for p in pts),
                     psucc_ok=sum(p["psucc_exact"] >= p["psucc_lb"] - 1e-12 for p in pts))
fil = [(k, case(k, 0.01)) for k in keys if case(k, 0.01) and case(k, 0.01).get("at_n_bound")]
a_ = lambda f: [f(c["at_n_bound"]) for _, c in fil]
S["total_ratio"] = mr(a_(lambda a: a["bound"] / a["meas"]))
S["leak_ratio"] = mr(a_(lambda a: a["bound_leak"] / a["meas_leak"]))
S["trot_ratio"] = mr(a_(lambda a: a["bound_trot"] / a["meas_trot"]))
S["energy_ratio"] = mr(a_(lambda a: a["energy_err_bound"] / a["energy_err"]))
S["psucc_rel_dev"] = mr(a_(lambda a: abs(a["psucc"] / a["psucc_lb"] - 1)))
S["meas_leak_range"] = (min(a_(lambda a: a["meas_leak"])), max(a_(lambda a: a["meas_leak"])))
S["meas_trot_max"] = max(a_(lambda a: a["meas_trot"]))
allc = [(k, c) for k in keys for c in M[k]["cases"] if c.get("n_bound") and c.get("n_meas")]
S["n_ratio_bound_meas"] = mr([c["n_bound"] / c["n_meas"] for _, c in allc])
S["n_ratio_bound_hyb"] = mr([c["n_bound"] / c["n_hyb"] for _, c in allc if c.get("n_hyb")])
S["n_ratio_hyb_meas"] = mr([c["n_hyb"] / c["n_meas"] for _, c in allc if c.get("n_hyb")])
hy = [(k, c) for k, c in allc if c.get("hybrid")]
S["hybrid_ok"] = dict(cases=len(hy), met=sum(c["hybrid"]["meas"] <= np.sqrt(c["eps"]) + 1e-12 for _, c in hy),
                      margin=mr([np.sqrt(c["eps"]) / c["hybrid"]["meas"] for _, c in hy]) if hy else None)
rich = [(s["rich_trot"] / s["meas_trot"]) for k in keys for c in M[k]["cases"] for s in c.get("sweep", []) if "rich_trot" in s and c.get("n_meas") and s["n"] >= c["n_meas"]]
S["rich"] = mr(rich) if rich else None
S["rich_count"] = len(rich)
sl = {}
for k in keys:
    c = case(k, 0.01)
    if c and c.get("sweep"):
        sw = [s for s in c["sweep"] if s["bound_trot"] < 0.9]
        if len(sw) >= 3:
            n = np.log([s["n"] for s in sw])
            sl[f"N{k[0]}_J{k[1]}"] = (float(np.polyfit(n, np.log([s["bound_trot"] for s in sw]), 1)[0]), float(np.polyfit(n, np.log([s["meas_trot"] for s in sw]), 1)[0]))
S["slopes"] = sl

# hybrid table
rows = ["| $N$ | $J_2$ | $\\varepsilon$ | $n_{\\mathrm{bound}}$ | $n_{\\mathrm{hyb}}$ | $n_{\\mathrm{meas}}$ | $D$ at $n_{\\mathrm{hyb}}$ | $\\sqrt\\varepsilon$ | steps saved vs. certified |", "|---|---|---|---|---|---|---|---|---|"]
for k in keys:
    for c in M[k]["cases"]:
        if c.get("hybrid") and c["eps"] in (0.01, 0.001):
            rows.append(f"| {k[0]} | {k[1]} | {c['eps']} | {cnt(c['n_bound'])} | {c['n_hyb']} | {c['n_meas']} | {c['hybrid']['meas']:.3f} | {np.sqrt(c['eps']):.3f} | {c['n_bound']/c['n_hyb']:.0f}$\\times$ |")
T["hybrid"] = "\n".join(rows)

# eps scaling (N <= 8)
rows = ["| $N$ | $J_2$ | $n_{\\mathrm{bound}}$ ($\\varepsilon=10^{-1}$) | ($10^{-2}$) | ($10^{-3}$) | ratio $10^{-3}/10^{-2}$ | $n_{\\mathrm{hyb}}$ ($10^{-2}$ / $10^{-3}$) | $n_{\\mathrm{meas}}$ ($10^{-2}$ / $10^{-3}$) |", "|---|---|---|---|---|---|---|---|"]
for k in keys:
    cs = [case(k, e) for e in (0.1, 0.01, 0.001)]
    if cs[2] is None or not cs[2].get("n_bound") or not cs[1].get("n_bound"):
        continue
    g = lambda c, f: (c.get(f) or "–")
    rows.append(f"| {k[0]} | {k[1]} | {cnt(cs[0]['n_bound']) if cs[0].get('n_bound') else '–'} | {cnt(cs[1]['n_bound'])} | {cnt(cs[2]['n_bound'])} | {cs[2]['n_bound']/cs[1]['n_bound']:.1f} | {g(cs[1],'n_hyb')} / {g(cs[2],'n_hyb')} | {g(cs[1],'n_meas')} / {g(cs[2],'n_meas')} |")
T["eps"] = "\n".join(rows)

# ------------------------------------------------------------------ cost (eps = 1e-2)
rows = ["| $N$ | $J_2$ | trial CX | CX / step | filter CX: certified | filter CX: hybrid | filter CX: measured | total, certified | total, hybrid | best direct MPS prep. |", "|---|---|---|---|---|---|---|---|---|---|"]
S["cost"] = {}
for k in keys:
    i, c = M[k]["inputs"], case(k, 0.01)
    if not c or not c.get("n_bound"):
        continue
    bd, bl = best_direct(k[0], k[1], 0.01)
    S["cost"][f"N{k[0]}_J{k[1]}"] = dict(bound=c["cx_total_bound"], hyb=c.get("cx_total_hyb"), meas=c.get("cx_total_meas"), direct=bd, step=i["cx_step_a2a"], trial=i["trial_cx"],
                                         filt_bound=c["cx_filter_bound"], filt_hyb=c.get("cx_filter_hyb"), filt_meas=c.get("cx_filter_meas"))
    rows.append(f"| {k[0]} | {k[1]} | {i['trial_cx']} | {i['cx_step_a2a']} | {cnt(c['cx_filter_bound'])} | {cnt(c['cx_filter_hyb']) if c.get('cx_filter_hyb') else '–'} | {cnt(c['cx_filter_meas']) if c.get('cx_filter_meas') else '–'} | {cnt(c['cx_total_bound'])} | {cnt(c['cx_total_hyb']) if c.get('cx_total_hyb') else '–'} | {cnt(bd) if bd else '–'} ({bl}) |")
T["cost"] = "\n".join(rows)
cx = S["cost"]
for J2 in (0.0, 0.4):
    sel = [(k[0], v) for k, v in ((k, cx.get(f"N{k[0]}_J{k[1]}")) for k in keys if k[1] == J2) if v]
    sel = [(n, v) for n, v in sel if n >= 6]
    if len(sel) >= 3:
        lg = np.log([n for n, _ in sel])
        fit = lambda f: float(np.polyfit(lg, np.log([v[f] for _, v in sel]), 1)[0])
        S[f"fit_J{J2}"] = dict(step=fit("step"), filt_bound=fit("filt_bound"), filt_hyb=fit("filt_hyb") if all(v["filt_hyb"] for _, v in sel) else None,
                              filt_meas=fit("filt_meas") if all(v["filt_meas"] for _, v in sel) else None)
        S[f"fit_J{J2}"]["n_bound"] = float(np.polyfit(lg, np.log([case((n, J2), 0.01)["n_bound"] for n, _ in sel]), 1)[0])
        S[f"fit_J{J2}"]["n_meas"] = float(np.polyfit(lg, np.log([case((n, J2), 0.01)["n_meas"] for n, _ in sel]), 1)[0])
        S[f"fit_J{J2}"]["N_range"] = (sel[0][0], sel[-1][0])
S["cx_bound_over_hyb"] = mr([v["filt_bound"] / v["filt_hyb"] for v in cx.values() if v["filt_hyb"]])
S["cx_hyb_over_meas"] = mr([v["filt_hyb"] / v["filt_meas"] for v in cx.values() if v["filt_hyb"] and v["filt_meas"]])
S["total_over_direct"] = {n: dict(bound=v["bound"] / v["direct"], hyb=(v["hyb"] / v["direct"] if v["hyb"] else None), meas=(v["meas"] / v["direct"] if v["meas"] else None)) for n, v in cx.items() if v["direct"]}
S["line_over_a2a"] = mr([M[k]["inputs"]["cx_step_line"] / M[k]["inputs"]["cx_step_a2a"] for k in keys])

# ------------------------------------------------------------------ extras tables
if TR:
    rows = ["| $J_2$ | trial state | $\\gamma$ | trial CX | $n_{\\mathrm{bound}}$ | total CX, certified | total CX, hybrid |", "|---|---|---|---|---|---|---|"]
    for r in TR:
        nm = f"DMRG $\\chi={r['param']}$" if r["kind"] == "dmrg" else f"approx. circuit $L={r['param']}$"
        if r.get("note"):
            rows.append(f"| {r['J2']} | {nm} | {r['gamma_est']:.4f} | {r['trial_cx']} | no filter needed | {r['trial_cx']} | {r['trial_cx']} |")
        else:
            rows.append(f"| {r['J2']} | {nm} | {r['gamma_est']:.4f} | {r['trial_cx']} | {cnt(r['n_bound'])} | {cnt(r['cx_total_bound'])} | {cnt(r['cx_total_hyb']) if r.get('cx_total_hyb') else '–'} |")
    T["trials"] = "\n".join(rows)
if SE:
    rows = ["| $N$ | input error | $n_{\\mathrm{bound}}$ | claimed $\\bar D$ | $\\bar D$ with true inputs | measured $D$ | claim violated? | target $\\sqrt\\varepsilon$ missed? |", "|---|---|---|---|---|---|---|---|"]
    for r in SE:
        if r["kind"] == "gap":
            lab = f"gap {r['param']*100:+.0f}\\%"
        else:
            lab = f"$\\gamma$ over-estimated ({r['gamma_used']:.3f} vs {r['gamma_true']:.3f})"
        rows.append(f"| {r['N']} | {lab} | {cnt(r['n_bound'])} | {r['bound_claimed']:.3f} | {r['bound_with_truth']:.3f} | {('%.3f' % r['meas']) if r['meas'] is not None else '–'} | {'yes' if r['claim_violated'] else 'no'} | {'yes' if r['target_missed'] else 'no'} |")
    T["sens"] = "\n".join(rows)
if DE:
    rows = ["| leakage share $\\varepsilon_\\ell/\\varepsilon$ | $m=4$ | $m=6$ | $m=8$ |", "|---|---|---|---|"]
    for fr in sorted({r["frac"] for r in DE}):
        cells = []
        for m in (4, 6, 8):
            r = [x for x in DE if x["frac"] == fr and x["m"] == m][0]
            cells.append("infeasible" if "error" in r else f"{cnt(r['n_bound'])} ({r['m_used']} pulses)")
        rows.append(f"| {fr} | " + " | ".join(cells) + " |")
    T["design"] = "\n".join(rows)
    vals = [r["n_bound"] for r in DE if "n_bound" in r]
    S["design_range"] = (min(vals), max(vals))
if BA:
    rows = ["| $N$ | $J_2$ | MPS $\\chi=2$ | MPS $\\chi=4$ | MPS $\\chi=8$ | approx. $L=1$ | approx. $L=3$ | generic |", "|---|---|---|---|---|---|---|---|"]
    for r in BA:
        ex = {e["chi"]: e for e in r["exact"] if "cx" in e}
        ap = {a["L"]: a for a in r["approx"]}
        f = lambda e: f"{e['cx']} ({sci(max(e['infid'], 1e-16))})" if e else "–"
        rows.append(f"| {r['N']} | {r['J2']} | {f(ex.get(2))} | {f(ex.get(4))} | {f(ex.get(8))} | {f(ap.get(1))} | {f(ap.get(3))} | {r['generic']['cx'] or '–'} |")
    T["baseline"] = "\n".join(rows)

(RES / "tables.json").write_text(json.dumps(T, indent=1))
(RES / "stats.json").write_text(json.dumps(S, indent=1, default=float))

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
C = {"bound": "#1f4e79", "hyb": "#e08a1e", "meas": "#c0392b", "trot": "#2e8b57", "direct": "#444444"}

sel = [(k, case(k, 0.01)) for k in keys if k[1] == 0.0 and k[0] in (4, 6, 8, 10) and case(k, 0.01) and case(k, 0.01).get("sweep")]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
for (k, c), ls in zip(sel, ["-", "--", ":", "-."]):
    sw = c["sweep"]; n = np.array([s["n"] for s in sw])
    axes[0].loglog(n, [s["bound_truth"] for s in sw], ls, color=C["bound"], lw=1.4)
    axes[0].loglog(n, [s["meas"] for s in sw], ls, color=C["meas"], lw=1.4)
    axes[1].loglog(n, [s["bound_trot"] for s in sw], ls, color=C["bound"], lw=1.4, label=f"$N={k[0]}$")
    axes[1].loglog(n, [s["meas_trot"] for s in sw], ls, color=C["trot"], lw=1.4)
axes[0].axhline(0.1, color="k", lw=0.8, ls=":")
axes[0].set_xlabel("Trotter steps $n$"); axes[0].set_ylabel("trace distance to $|E_0\\rangle$")
axes[0].set_title("(a) total: certificate (blue), measured (red)", fontsize=9)
axes[1].set_xlabel("Trotter steps $n$"); axes[1].set_ylabel("Trotter part")
axes[1].set_title("(b) Trotter part: certificate (blue), measured (green)", fontsize=9)
axes[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "fig_sweep.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.0, 4.4))
for J2, mk, mf in [(0.0, "o", True), (0.4, "s", False)]:
    ks_ = [k for k in keys if k[1] == J2 and cx.get(f"N{k[0]}_J{k[1]}")]
    N_ = [k[0] for k in ks_]
    fc = lambda c: C[c] if mf else "white"
    ax.semilogy(N_, [cx[f"N{k[0]}_J{k[1]}"]["bound"] for k in ks_], mk + "-", color=C["bound"], mfc=fc("bound"), label=f"filter, certified $n$, $J_2={J2}$")
    hy_ = [(k[0], cx[f"N{k[0]}_J{k[1]}"]["hyb"]) for k in ks_ if cx[f"N{k[0]}_J{k[1]}"]["hyb"]]
    ax.semilogy([a for a, _ in hy_], [b for _, b in hy_], mk + "-", color=C["hyb"], mfc=fc("hyb"), label=f"filter, hybrid $n$, $J_2={J2}$")
    me_ = [(k[0], cx[f"N{k[0]}_J{k[1]}"]["meas"]) for k in ks_ if cx[f"N{k[0]}_J{k[1]}"]["meas"]]
    ax.semilogy([a for a, _ in me_], [b for _, b in me_], mk + "--", color=C["meas"], mfc=fc("meas"), label=f"filter, measured $n$ (oracle), $J_2={J2}$")
    dd_ = [(k[0], cx[f"N{k[0]}_J{k[1]}"]["direct"]) for k in ks_ if cx[f"N{k[0]}_J{k[1]}"]["direct"]]
    ax.semilogy([a for a, _ in dd_], [b for _, b in dd_], mk + ":", color=C["direct"], mfc=fc("direct"), label=f"best direct MPS preparation, $J_2={J2}$")
ax.set_xlabel("$N$"); ax.set_ylabel("total CX (trial + filter), $\\varepsilon=10^{-2}$"); ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout(); fig.savefig(FIG / "fig_cx.png", dpi=200); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 1.9)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(2, 31)
boxes = [(1, "classical front end", "DMRG: $E_0,E_1,E_{\\mathrm{top}}$\nMPS, compiled to a circuit\n$\\gamma$, $\\Delta$, $W$"),
         (21, "scale", "$H_s=(H-E_0)/W$\nspec $\\subset\\{0\\}\\cup[\\Delta,1]$\n$\\alpha$ (commutators)"),
         (41, "filter", "$F(E)=\\prod_i\\cos(Et_i{+}\\phi_i)$\ncertified $\\eta$, $F(0)^2$\ntotal time $T$"),
         (61, "choose $n$", "certified: $\\bar D(n)\\leq\\sqrt{\\varepsilon}$\nhybrid: $\\sqrt{\\bar\\ell}$ + n vs 2n\nTrotter estimate"),
         (81, "run / check", "Trotterised filter\n$\\langle H\\rangle-E_0$, $P_{\\mathrm{succ}}$\ncost: CX")]
for x, title, body in boxes:
    ax.add_patch(FancyBboxPatch((x, 4), 17.5, 24, boxstyle="round,pad=0.4,rounding_size=1.2", fc="#eef3f8", ec="#1f4e79", lw=1))
    ax.text(x + 8.75, 25.2, title, ha="center", va="center", fontsize=8, weight="bold", color="#1f4e79")
    ax.text(x + 8.75, 14.5, body, ha="center", va="center", fontsize=6.6)
for x in (18.9, 38.9, 58.9, 78.9):
    ax.annotate("", xy=(x + 2.4, 16), xytext=(x - 0.2, 16), arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.2))
fig.savefig(FIG / "fig_workflow.png", dpi=220, bbox_inches="tight"); plt.close(fig)
print(json.dumps({k: v for k, v in S.items() if k not in ("input_err", "cost", "total_over_direct", "slopes")}, indent=1, default=float))
