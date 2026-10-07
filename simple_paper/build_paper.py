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
