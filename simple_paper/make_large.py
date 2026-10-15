"""make_large.py -- tables/statistics for the DMRG/MPS-only regime (N = 14-32), read from ../j1j2_filter/results.
Writes results/large.json, merged into the paper by build_paper.py."""
import json, glob, pathlib
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "j1j2_filter" / "results" / "data"
VAL = json.load(open(HERE.parent / "j1j2_filter" / "results" / "dmrg_validation.json"))
sci = lambda x: f"{float(x):.1e}".replace("e-0", "e-").replace("e+0", "e+")
cnt = lambda x: f"{int(round(x)):,}"
R = {}
for f in glob.glob(str(DATA / "dmrg_N*_eps*.json")):
    r = json.load(open(f)); R[(r["N"], round(r["J2"], 4), r["eps"])] = r
J2S = (0.0, 0.2411, 0.4)
Ns = (14, 16, 18, 20, 24, 28, 32)
out, S = {}, {}

rows = ["| $N$ | $J_2$ | $\\gamma$ | $\\Delta$ | $\\alpha$ | pulses | $n$ (guaranteed) | $T$ (phys.) | CX total | non-Clifford Rz | $P_{\\mathrm{succ}}\\ge$ |", "|---|---|---|---|---|---|---|---|---|---|---|"]
missing, pts = [], []
for J2 in J2S:
    for N in Ns:
        r = R.get((N, J2, 0.01))
        if r is None:
            continue
        d = r["designs"]["guaranteed"]
        i = r["inputs"]
        if "n" not in d:
            rows.append(f"| {N} | {J2:g} | {i['gamma']:.3f} | {i['gap_scaled']:.4f} | {i['alpha']:.4f} | – | not achievable | – | – | – | – |")
            missing.append((N, J2, i["gamma"])); continue
        t = d["counts_a2a"]["total"]
        fl = r["design"]["source"] in ("floor", "none")
        rows.append(f"| {N} | {J2:g} | {i['gamma']:.3f} | {i['gap_scaled']:.4f} | {i['alpha']:.4f} | {d['n_pulses_used']} | {cnt(d['n'])}{'' if fl else '$^\\ddagger$'} | {d['T_phys']:.1f} | {t['cx']:.2e} | {t['rz_nonclifford']:.2e} | {d['p_succ_lb']:.3f} |")
        if fl and N >= 16:
            pts.append((N, J2, d["n"], t["cx"], t["rz_nonclifford"], d["T_phys"]))
out["large"] = "\n".join(rows)
S["large_missing"] = [(N, J2, round(g, 3)) for N, J2, g in missing]
S["large_n"] = (min(p[2] for p in pts), max(p[2] for p in pts))
S["large_cx"] = (min(p[3] for p in pts), max(p[3] for p in pts))
S["large_rz"] = (min(p[4] for p in pts), max(p[4] for p in pts))
S["large_T"] = (min(p[5] for p in pts), max(p[5] for p in pts))
S["large_Nmax"] = max(p[0] for p in pts)
fits = {}
for J2 in J2S:
    q = [p for p in pts if p[1] == J2] + [(N, J2, R[(N, J2, 0.01)]["designs"]["guaranteed"]["n"], R[(N, J2, 0.01)]["designs"]["guaranteed"]["counts_a2a"]["total"]["cx"], 0, 0)
                                          for N in (6, 8, 10, 12, 14) if (N, J2, 0.01) in R and "n" in R[(N, J2, 0.01)]["designs"]["guaranteed"] and R[(N, J2, 0.01)]["design"]["source"] == "floor"]
    q.sort()
    if len(q) >= 3:
        lg = np.log([p[0] for p in q])
        fits[J2] = dict(cx=float(np.polyfit(lg, np.log([p[3] for p in q]), 1)[0]), n=float(np.polyfit(lg, np.log([p[2] for p in q]), 1)[0]), Nrange=(q[0][0], q[-1][0]))
S["large_fits"] = {str(k): v for k, v in fits.items()}

# DMRG accuracy beyond the main table: ED rows (N = 16-20) and run-to-run consistency (N = 24-32)
rows = ["| $N$ | $J_2$ | $|E_0-E_0^{\\mathrm{ED}}|$ | $|E_1-E_1^{\\mathrm{ED}}|$ | $|E_{\\mathrm{top}}-E_{\\mathrm{top}}^{\\mathrm{ED}}|$ | $1-|\\langle\\psi_0|\\psi_0^{\\mathrm{ED}}\\rangle|^2$ |", "|---|---|---|---|---|---|"]
ed = [v for k, v in VAL.items() if k.startswith("ed_")]
for v in sorted(ed, key=lambda v: (v["J2"], v["N"])):
    rows.append(f"| {v['N']} | {v['J2']:g} | {sci(v['dE0'])} | {sci(v['dE1'])} | {sci(v['dEtop'])} | {sci(v['one_minus_overlap'])} |")
out["dmrg_ed"] = "\n".join(rows)
S["ed_max"] = max(max(float(v["dE0"]), float(v["dE1"]), float(v["dEtop"])) for v in ed)
rows = ["| $N$ | $J_2$ | run-to-run $|\\Delta E_0|$ | $|\\Delta E_1|$ | $|\\Delta E_{\\mathrm{top}}|$ | $1-|\\langle\\psi_0|\\psi_0'\\rangle|^2$ | variance of $\\psi_0$ | variance of $\\psi_1$ |", "|---|---|---|---|---|---|---|---|"]
co = [v for k, v in VAL.items() if k.startswith("consist_")]
for v in sorted(co, key=lambda v: (v["J2"], v["N"])):
    rows.append(f"| {v['N']} | {v['J2']:g} | {sci(v['dE0'])} | {sci(v['dE1'])} | {sci(v['dEtop'])} | {sci(v['one_minus_overlap'])} | {sci(v['var0'])} | {sci(v['var1'])} |")
out["dmrg_cons"] = "\n".join(rows)
S["cons_max"] = max(max(float(v["dE0"]), float(v["dE1"]), float(v["dEtop"])) for v in co)
S["cons_var_max"] = max(max(float(v["var0"]), float(v["var1"])) for v in co)
g = [v for k, v in VAL.items() if k.startswith("gamma_")]
S["gamma_diff_max"] = max(abs(float(v["diff"])) for v in g)
out["stats"] = S
(HERE / "results" / "large.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(S, indent=1, default=float))
