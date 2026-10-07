"""build_paper.py -- fill paper_template.md with the numbers/tables from results/ and run pandoc (html, docx, tex)."""
import json, pathlib, subprocess, re
HERE = pathlib.Path(__file__).resolve().parent
T = json.loads((HERE / "results" / "tables.json").read_text())
S = json.loads((HERE / "results" / "stats.json").read_text())
txt = (HERE / "paper_template.md").read_text()

def rng(key, d=1, scale=1.0):
    lo, med, hi = S[key]
    f = lambda v: f"{v*scale:.{d}f}"
    return f(lo), f(med), f(hi)

vals = {k: T[k] for k in ("main", "components", "energy", "chi")}
for key, d in [("total_ratio", 1), ("n_ratio", 0), ("leak_ratio", 1), ("trot_ratio", 0), ("energy_ratio", 0)]:
    lo, med, hi = rng(key, d)
    vals.update({f"{key}_lo": lo, f"{key}_med": med, f"{key}_hi": hi})
lo, med, hi = rng("rich_ratio_n_ge_nemp", 2)
vals.update(rich_lo=lo, rich_med=med, rich_hi=hi)
R = S["resources"]
vals.update(resources=T["resources"],
    cxstep0=f"{R['J2_0.0']['cx_step']:.1f}", cxstep4=f"{R['J2_0.4']['cx_step']:.1f}",
    nb0=f"{R['J2_0.0']['n_bound']:.1f}", ne0=f"{R['J2_0.0']['n_emp']:.1f}",
    cxb0=f"{R['J2_0.0']['cx_bound']:.1f}", cxb4=f"{R['J2_0.4']['cx_bound']:.1f}",
    cxe0=f"{R['J2_0.0']['cx_emp']:.1f}", cxe4=f"{R['J2_0.4']['cx_emp']:.1f}",
    cxr_lo=f"{R['cx_ratio_bound_over_emp'][0]:.0f}", cxr_hi=f"{R['cx_ratio_bound_over_emp'][1]:.0f}",
    gen_growth=f"{R['generic_growth_per_2qubits']:.1f}$\\times$",
    eg4=f"{R['emp_over_generic']['N4_J20.0']:.0f}", eg12=f"{R['emp_over_generic']['N12_J20.0']:.1f}",
    bg_lo=f"{min(R['bound_over_generic'].values()):.0f}", bg_hi=f"{max(R['bound_over_generic'].values()):.0f}")
vals["psucc_dev_max"] = f"{S['psucc_rel_dev'][2]*100:.1f}"
for k, v in vals.items():
    txt = txt.replace("{{" + k + "}}", v)
left = re.findall(r"\{\{[a-z_]+\}\}", txt)
assert not left, left
(HERE / "paper.md").write_text(txt)

common = ["pandoc", "paper.md", "--from=markdown+tex_math_dollars+pipe_tables", "--resource-path=."]
subprocess.run(common + ["-o", "paper.docx"], cwd=HERE, check=True)
subprocess.run(common + ["-o", "paper.tex", "--standalone"], cwd=HERE, check=True)
css = HERE / "paper.css"
css.write_text("body{max-width:52em;margin:2em auto;padding:0 1em;font:16px/1.5 Georgia,serif}table{border-collapse:collapse;font-size:.8em;margin:1em 0}"
               "td,th{border-bottom:1px solid #ccc;padding:2px 6px}img{max-width:100%}h1.title{font-size:1.6em}h2{margin-top:1.6em}")
subprocess.run(common + ["-o", "paper.html", "--standalone", "--mathml", "--embed-resources", "--css=paper.css"], cwd=HERE, check=True)
print("built paper.md, paper.docx, paper.tex, paper.html")
