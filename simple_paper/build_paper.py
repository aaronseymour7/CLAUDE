"""build_paper.py -- fill paper_template.md with tables/numbers from results/ and run pandoc (html, docx, tex)."""
import json, pathlib, subprocess, re
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
T = json.loads((HERE / "results" / "tables.json").read_text())
S = json.loads((HERE / "results" / "stats.json").read_text())
P = json.loads((HERE / "results" / "pinsker.json").read_text())
txt = (HERE / "paper_template.md").read_text()

V = dict(T)
def rng(name, key, d=1, scale=1.0):
    lo, med, hi = S[key]
    V[f"{name}_lo"], V[f"{name}_med"], V[f"{name}_hi"] = (f"{x*scale:.{d}f}" for x in (lo, med, hi))
for name, key, d in [("tot", "total_ratio", 1), ("leak", "leak_ratio", 1), ("trot", "trot_ratio", 0), ("en", "energy_ratio", 0),
                     ("nbm", "n_ratio_bound_meas", 0), ("nbh", "n_ratio_bound_hyb", 0), ("nhm", "n_ratio_hyb_meas", 1),
                     ("hmarg", "hybrid_ok" , 1) if False else ("cxbh", "cx_bound_over_hyb", 0), ("cxhm", "cx_hyb_over_meas", 1),
                     ("line", "line_over_a2a", 1), ("rich", "rich", 3)]:
    rng(name, key, d)
V["hmarg_lo"], V["hmarg_med"], V["hmarg_hi"] = (f"{x:.1f}" for x in S["hybrid_ok"]["margin"])
V["hyb_met"], V["hyb_cases"] = f"{int(S['hybrid_ok']['met'])}", f"{S['hybrid_ok']['cases']}"
v = S["validity"]
assert all(v[k] == v['points'] for k in ('bound_ok','leak_ok','trot_ok','triangle_ok','energy_ok','psucc_ok')), v
V.update(npts=str(v["points"]), max_E=f"{S['max_E_err']:.0e}".replace("e-0", "e-"), max_gap=f"{S['max_gap_rel']:.0e}".replace("e-0", "e-"),
         max_gam=f"{S['max_gamma_err']:.0e}".replace("e-0", "e-"), psucc_max=f"{S['psucc_rel_dev'][2]*100:.2f}",
         leak_meas_lo=f"{S['meas_leak_range'][0]:.3f}", leak_meas_hi=f"{S['meas_leak_range'][1]:.3f}", trot_meas_max=f"{S['meas_trot_max']:.3f}",
         rich_n=str(S["rich_count"]), des_lo=f"{S['design_range'][0]:,}", des_hi=f"{S['design_range'][1]:,}")
sl = list(S["slopes"].values())
V["slope_b_lo"], V["slope_b_hi"] = f"{min(a for a,_ in sl):.2f}", f"{max(a for a,_ in sl):.2f}"
V["slope_m_lo"], V["slope_m_hi"] = f"{min(b for _,b in sl):.2f}", f"{max(b for _,b in sl):.2f}"
for J, tag in ((0.0, "0"), (0.4, "4")):
    f = S[f"fit_J{J}"]
    for k in ("step", "filt_bound", "filt_hyb", "filt_meas", "n_bound", "n_meas"):
        V[f"fit{tag}_{k}"] = f"{f[k]:.1f}"
V["fit_Nlo"], V["fit_Nhi"] = S["fit_J0.0"]["N_range"]
tod = S["total_over_direct"]
for k in ("bound", "hyb", "meas"):
    xs = [x[k] for x in tod.values() if x[k]]
    V[f"dir_{k}_lo"], V[f"dir_{k}_hi"] = f"{min(xs):,.0f}", f"{max(xs):,.0f}"
V["pk_lo"], V["pk_hi"] = f"{min(p['P_cap'] for p in P):.2f}", f"{max(p['P_cap'] for p in P):.2f}"
V["pP_lo"], V["pP_hi"] = f"{min(p['P'] for p in P):.2f}", f"{max(p['P'] for p in P):.2f}"
V["pn"] = str(len(P))
for k, val in V.items():
    txt = txt.replace("{{" + k + "}}", str(val))
left = re.findall(r"\{\{[A-Za-z_0-9]+\}\}", txt)
assert not left, left
(HERE / "paper.md").write_text(txt)
common = ["pandoc", "paper.md", "--from=markdown+tex_math_dollars+pipe_tables", "--resource-path=."]
subprocess.run(common + ["-o", "paper.docx"], cwd=HERE, check=True)
subprocess.run(common + ["-o", "paper.tex", "--standalone"], cwd=HERE, check=True)
(HERE / "paper.css").write_text("body{max-width:54em;margin:2em auto;padding:0 1em;font:16px/1.5 Georgia,serif}table{border-collapse:collapse;font-size:.78em;margin:1em 0}"
                               "td,th{border-bottom:1px solid #ccc;padding:2px 6px}img{max-width:100%}h1.title{font-size:1.6em}h2{margin-top:1.6em}")
subprocess.run(common + ["-o", "paper.html", "--standalone", "--mathml", "--embed-resources", "--css=paper.css"], cwd=HERE, check=True)
print("built paper.md, paper.docx, paper.tex, paper.html")
