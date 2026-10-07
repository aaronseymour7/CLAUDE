"""make_report.py -- turn results/*.json into figures and markdown tables (results/tables.md)."""
import json, glob, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
RES, FIG = HERE / "results", HERE / "figs"
cases = [json.loads(open(f).read()) for f in sorted(glob.glob(str(RES / "case_*.json")))]
cases = [c for c in cases if "error" not in c]
main = [c for c in cases if c["chi"] == 2]
key = lambda c: (c["eps"] != 0.01, c["N"], c["J2"])
main.sort(key=key)

def fmt(x, d=3):
    return "–" if x is None else (f"{x:.{d}g}")

# ------------------------------------------------------------------ tables
T = {}
rows = ["| $N$ | $J_2$ | $\\varepsilon$ | $\\gamma$ | $\\Delta$ | $\\alpha$ | $T$ | $\\eta$ | $n_{\\mathrm{bound}}$ | $n_{\\mathrm{emp}}$ | bound $\\bar D$ | measured $D$ | $\\bar D/D$ |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for c in main:
    if "at_n_bound" not in c:
        rows.append(f"| {c['N']} | {c['J2']} | {c['eps']} | {c['gamma']:.4f} | {c['gap_scaled']:.4f} | {c['alpha']:.3f} | – | – | – | – | – | {c['D_trial']:.3f} (trial) | – |")
        continue
    a = c["at_n_bound"]
    rows.append(f"| {c['N']} | {c['J2']} | {c['eps']} | {c['gamma']:.4f} | {c['gap_scaled']:.4f} | {c['alpha']:.3f} | {a['T']:.1f} | {a['eta']:.3f} | {c['n_bound']} | {c['n_emp']} | {a['bound']:.4f} | {a['meas']:.4f} | {a['bound']/a['meas']:.1f} |")
T["main"] = "\n".join(rows)

rows = ["| $N$ | $J_2$ | $\\bar D_{\\mathrm{leak}}$ | $D_{\\mathrm{leak}}$ | ratio | $\\bar D_{\\mathrm{Trot}}$ | $D_{\\mathrm{Trot}}$ | ratio | $D\\le D_{\\mathrm{leak}}+D_{\\mathrm{Trot}}$ | $P_{\\mathrm{succ}}$ | $\\gamma F(0)^2$ |",
        "|---|---|---|---|---|---|---|---|---|---|---|"]
for c in main:
    if c["eps"] != 0.01 or "at_n_bound" not in c:
        continue
    a = c["at_n_bound"]
    rows.append(f"| {c['N']} | {c['J2']} | {a['bound_leak']:.4f} | {a['meas_leak']:.4f} | {a['bound_leak']/a['meas_leak']:.2f} | {a['bound_trot']:.4f} | {a['meas_trot']:.5f} | {a['bound_trot']/a['meas_trot']:.0f} | {'yes' if a['triangle_ok'] else 'NO'} | {a['psucc']:.3f} | {a['psucc_lb']:.3f} |")
T["components"] = "\n".join(rows)

rows = ["| $N$ | $J_2$ | $n$ | $\\bar D$ | $D$ | $W\\bar D^2$ | $\\langle H\\rangle-E_0$ | Richardson $\\hat D_{\\mathrm{Trot}}$ | $D_{\\mathrm{Trot}}$ |",
        "|---|---|---|---|---|---|---|---|---|"]
for c in main:
    if c["eps"] != 0.01 or "at_n_bound" not in c:
        continue
    a = c["at_n_bound"]
    r = a.get("rich_trot")
    rows.append(f"| {c['N']} | {c['J2']} | {a['n']} | {a['bound']:.4f} | {a['meas']:.4f} | {a['energy_err_bound']:.4f} | {a['energy_err']:.5f} | {fmt(r,3) if r is not None else '–'} | {a['meas_trot']:.5f} |")
T["energy"] = "\n".join(rows)

chis = sorted([c for c in cases if c["N"] == 8 and c["J2"] == 0.0 and c["eps"] == 0.01], key=lambda c: c["chi"])
rows = ["| $\\chi$ | $\\gamma$ | $\\sqrt{1-\\gamma}$ (trial) | $T$ | $\\eta$ | $n_{\\mathrm{bound}}$ | $n_{\\mathrm{emp}}$ | bound $\\bar D$ | measured $D$ |", "|---|---|---|---|---|---|---|---|---|"]
for c in chis:
    if "at_n_bound" not in c:
        rows.append(f"| {c['chi']} | {c['gamma']:.4f} | {c['D_trial']:.3f} | – | – | no filter needed | – | – | {c['D_trial']:.3f} |")
        continue
    a = c["at_n_bound"]
    rows.append(f"| {c['chi']} | {c['gamma']:.4f} | {c['D_trial']:.3f} | {a['T']:.1f} | {a['eta']:.3f} | {c['n_bound']} | {c['n_emp']} | {a['bound']:.4f} | {a['meas']:.4f} |")
T["chi"] = "\n".join(rows)

# validity statistics over every simulated point
pts = []
for c in cases:
    if "at_n_bound" in c:
        pts.append(c["at_n_bound"])
    pts += c.get("sweep", [])
T["stats"] = dict(points=len(pts), bound_ok=sum(p["bound_ok"] for p in pts), leak_ok=sum(p["leak_ok"] for p in pts),
                  trot_ok=sum(p["trot_ok"] for p in pts), tri_ok=sum(p["triangle_ok"] for p in pts),
                  energy_ok=sum(p["energy_err"] <= p["energy_err_bound"] for p in pts),
                  psucc_ok=sum(p["psucc_exact"] >= p["psucc_lb"] - 1e-12 for p in pts))
(RES / "tables.json").write_text(json.dumps(T, indent=1))
print(T["stats"])

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
C = {"bound": "#1f4e79", "meas": "#c0392b", "leak": "#7f8c8d", "trot": "#2e8b57"}

sel = [c for c in main if c["eps"] == 0.01 and "sweep" in c and c["J2"] == 0.0]
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
for c, ls in zip(sel, ["-", "--", ":", "-."]):
    sw = c["sweep"]
    n = np.array([s["n"] for s in sw])
    lab = f"$N={c['N']}$"
    axes[0].loglog(n, [s["bound"] for s in sw], ls, color=C["bound"], lw=1.4)
    axes[0].loglog(n, [s["meas"] for s in sw], ls, color=C["meas"], lw=1.4)
    axes[1].loglog(n, [s["bound_trot"] for s in sw], ls, color=C["bound"], lw=1.4, label=lab)
    axes[1].loglog(n, [s["meas_trot"] for s in sw], ls, color=C["trot"], lw=1.4)
axes[0].axhline(np.sqrt(0.01), color="k", lw=0.8, ls=":")
axes[0].text(axes[0].get_xlim()[0] * 1.2, 0.108, r"$\sqrt{\varepsilon}$", fontsize=8)
axes[0].set_xlabel("Trotter steps $n$"); axes[0].set_ylabel("trace distance to $|E_0\\rangle$")
axes[0].set_title("(a) total: bound (blue) vs measured (red)", fontsize=9)
axes[1].set_xlabel("Trotter steps $n$"); axes[1].set_ylabel("Trotter part")
axes[1].set_title("(b) Trotter part: bound (blue) vs measured (green)", fontsize=9)
axes[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "fig_sweep.png", dpi=200); plt.close(fig)

# n_bound vs n_emp
fig, ax = plt.subplots(figsize=(3.6, 3.0))
for J2, mk in [(0.0, "o"), (0.4, "s")]:
    cc = [c for c in main if c["eps"] == 0.01 and c["J2"] == J2 and c.get("n_bound")]
    ax.semilogy([c["N"] for c in cc], [c["n_bound"] for c in cc], mk + "-", color=C["bound"], mfc="white" if J2 else C["bound"], label=f"bound, $J_2={J2}$")
    ax.semilogy([c["N"] for c in cc], [c["n_emp"] for c in cc], mk + "--", color=C["meas"], mfc="white" if J2 else C["meas"], label=f"measured, $J_2={J2}$")
ax.set_xlabel("$N$"); ax.set_ylabel("Trotter steps for $D\\leq 0.1$"); ax.legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(FIG / "fig_steps.png", dpi=200); plt.close(fig)
print("figures written")

# ------------------------------------------------------------------ extra statistics for the text
S = {}
fil = [c for c in main if "at_n_bound" in c]
r = lambda xs: (float(min(xs)), float(np.median(xs)), float(max(xs)))
S["n_ratio"] = r([c["n_bound"] / c["n_emp"] for c in fil])
S["total_ratio"] = r([c["at_n_bound"]["bound"] / c["at_n_bound"]["meas"] for c in fil])
S["leak_ratio"] = r([c["at_n_bound"]["bound_leak"] / c["at_n_bound"]["meas_leak"] for c in fil])
S["trot_ratio"] = r([c["at_n_bound"]["bound_trot"] / c["at_n_bound"]["meas_trot"] for c in fil])
S["energy_ratio"] = r([c["at_n_bound"]["energy_err_bound"] / c["at_n_bound"]["energy_err"] for c in fil])
S["psucc_rel_dev"] = r([abs(c["at_n_bound"]["psucc"] / c["at_n_bound"]["psucc_lb"] - 1) for c in fil])
rich = [(s["rich_trot"] / s["meas_trot"], s["n"], c["n_emp"]) for c in main for s in c.get("sweep", []) if "rich_trot" in s and c.get("n_emp") and s["n"] >= c["n_emp"]]
S["rich_ratio_n_ge_nemp"] = r([x[0] for x in rich]) if rich else None
S["rich_count"] = len(rich)
S["sweep_bound_slope"] = {}
for c in main:
    if "sweep" in c and c["eps"] == 0.01:
        sw = [s for s in c["sweep"] if s["bound_trot"] < 0.9]
        if len(sw) >= 3:
            n = np.log([s["n"] for s in sw]); b = np.log([s["bound_trot"] for s in sw]); m_ = np.log([s["meas_trot"] for s in sw])
            S["sweep_bound_slope"][f"N{c['N']}_J2{c['J2']}"] = (float(np.polyfit(n, b, 1)[0]), float(np.polyfit(n, m_, 1)[0]))
(RES / "stats.json").write_text(json.dumps(S, indent=1))
print(json.dumps(S, indent=1))

# ------------------------------------------------------------------ workflow figure
from matplotlib.patches import FancyBboxPatch
fig, ax = plt.subplots(figsize=(7.4, 1.9)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(2, 31)
boxes = [(1, "inputs", "$H=\\sum_g H_g$\nenergies $E_0,E_1,E_{\\mathrm{top}}$\ntrial state $|\\psi\\rangle$"),
         (21, "scaled quantities", "$H_s=(H-E_0)/W$,  $\\Delta$\n$\\gamma=|\\langle E_0|\\psi\\rangle|^2$\n$\\alpha$ (commutators)"),
         (41, "filter", "$F(E)=\\prod_i\\cos(Et_i{+}\\phi_i)$\ncertified $\\eta$, $F(0)^2$\ntotal time $T$"),
         (61, "bound", "$\\bar D=\\sqrt{\\bar\\ell}+\\epsilon_T/\\sqrt{p_g}$\n$\\epsilon_T=\\alpha T^2/2n$\nchoose $n$: $\\bar D\\leq\\sqrt{\\varepsilon}$"),
         (81, "run / check", "Trotterised filter\n$D$, $\\langle H\\rangle-E_0$,\n$P_{\\mathrm{succ}}$")]
for x, title, body in boxes:
    ax.add_patch(FancyBboxPatch((x, 4), 17.5, 24, boxstyle="round,pad=0.4,rounding_size=1.2", fc="#eef3f8", ec="#1f4e79", lw=1))
    ax.text(x + 8.75, 25.2, title, ha="center", va="center", fontsize=8.5, weight="bold", color="#1f4e79")
    ax.text(x + 8.75, 14.5, body, ha="center", va="center", fontsize=6.9)
for x in (18.9, 38.9, 58.9, 78.9):
    ax.annotate("", xy=(x + 2.4, 16), xytext=(x - 0.2, 16), arrowprops=dict(arrowstyle="->", color="#1f4e79", lw=1.2))
fig.savefig(FIG / "fig_workflow.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ resource scaling (needs results/resources.json from resources.py)
rf = RES / "resources.json"
if rf.exists():
    R = sorted(json.loads(rf.read_text()), key=lambda r: (r["J2"], r["N"]))
    rows = ["| $N$ | $J_2$ | CX / step (all-to-all) | CX / step (line) | $n_{\\mathrm{bound}}$ | $n_{\\mathrm{emp}}$ | CX bound (all-to-all) | CX measured-$n$ (all-to-all) | CX bound (line) | exact generic state prep |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for r in R:
        lc = r["cx_step_line_amortised"]
        rows.append(f"| {r['N']} | {r['J2']} | {r['cx_step_a2a']} | {lc:.0f} | {r['n_bound']} | {r['n_emp']} | {r['cx_bound']:,} | {r['cx_emp']:,} | {int(round(lc*r['n_bound'])):,} | {r['generic_prep_cx']:,} |")
    T["resources"] = "\n".join(rows)
    RS = {}
    for J2 in (0.0, 0.4):
        sel_ = [r for r in R if r["J2"] == J2 and r["N"] >= 6]
        N_ = np.log([r["N"] for r in sel_])
        RS[f"J2_{J2}"] = dict(
            cx_step=float(np.polyfit(N_, np.log([r["cx_step_a2a"] for r in sel_]), 1)[0]),
            n_bound=float(np.polyfit(N_, np.log([r["n_bound"] for r in sel_]), 1)[0]),
            n_emp=float(np.polyfit(N_, np.log([r["n_emp"] for r in sel_]), 1)[0]),
            cx_bound=float(np.polyfit(N_, np.log([r["cx_bound"] for r in sel_]), 1)[0]),
            cx_emp=float(np.polyfit(N_, np.log([r["cx_emp"] for r in sel_]), 1)[0]))
    RS["generic_growth_per_2qubits"] = float(np.mean([R_["generic_prep_cx"] for R_ in R if R_["N"] == 12 and R_["J2"] == 0.0]) / np.mean([R_["generic_prep_cx"] for R_ in R if R_["N"] == 10 and R_["J2"] == 0.0]))
    RS["cx_ratio_bound_over_emp"] = r_ = (min(x["cx_bound"] / x["cx_emp"] for x in R), max(x["cx_bound"] / x["cx_emp"] for x in R))
    RS["emp_over_generic"] = {f"N{x['N']}_J2{x['J2']}": x["cx_emp"] / x["generic_prep_cx"] for x in R}
    RS["bound_over_generic"] = {f"N{x['N']}_J2{x['J2']}": x["cx_bound"] / x["generic_prep_cx"] for x in R}
    S["resources"] = RS
    (RES / "stats.json").write_text(json.dumps(S, indent=1))
    (RES / "tables.json").write_text(json.dumps(T, indent=1))
    fig, ax = plt.subplots(figsize=(3.9, 3.1))
    for J2, mk in [(0.0, "o"), (0.4, "s")]:
        rr = [r for r in R if r["J2"] == J2]
        mf = "white" if J2 else None
        ax.semilogy([r["N"] for r in rr], [r["cx_bound"] for r in rr], mk + "-", color=C["bound"], mfc=mf or C["bound"], label=f"filter, bound $n$, $J_2={J2}$")
        ax.semilogy([r["N"] for r in rr], [r["cx_emp"] for r in rr], mk + "--", color=C["meas"], mfc=mf or C["meas"], label=f"filter, measured $n$, $J_2={J2}$")
    g_ = sorted({(r["N"], r["generic_prep_cx"]) for r in R})
    ax.semilogy([a for a, _ in g_], [b for _, b in g_], "k:", marker="^", label="exact generic state prep")
    ax.set_xlabel("$N$"); ax.set_ylabel("CX count"); ax.legend(frameon=False, fontsize=6)
    fig.tight_layout(); fig.savefig(FIG / "fig_cx.png", dpi=200); plt.close(fig)
    print(json.dumps(RS, indent=1))
