"""Section A fits (ln C = c + a ln(W/Delta) + b ln(1/gamma) + alpha ln N) for the authors' v4 and the improved option,
J2 -> H, staggered knob, analytic model, pooled N. Usage: python j1j2/qpf_design_compare_fits.py ANALYSIS_DIR"""
import sys, os
import numpy as np
import pandas as pd

d = sys.argv[1]
w = pd.read_csv(os.path.join(d, "design_compare_wide.csv"))
pts = pd.read_csv(os.path.join(d, "design_compare_points.csv"))
W = pts.groupby(["N", "gap_target"])["W"].first().rename("W").reset_index()
g = pts.groupby(["N", "gap_target", "knob", "gamma_target"])["gamma"].first().rename("gamma").reset_index()
w = w.merge(W, on=["N", "gap_target"]).merge(g, on=["N", "gap_target", "knob", "gamma_target"])
rows = []
for (route, knob, proto, eps), s in w.groupby(["route", "knob", "protocol", "eps"]):
    if route != "J2 -> H":
        continue
    for col in ("QPF (t_i) v4", "ms best", "minimax any", "all options"):
        if col not in s:
            continue
        q = s.dropna(subset=[col])
        X = np.c_[np.ones(len(q)), np.log(q.W / q.gap_target), np.log(1 / q.gamma), np.log(q.N)]
        y = np.log(q[col])
        co, *_ = np.linalg.lstsq(X, y, rcond=None)
        r2 = 1 - np.sum((y - X @ co) ** 2) / np.sum((y - y.mean()) ** 2)
        rows.append(dict(route=route, knob=knob, protocol=proto, eps=eps, schedule=col, n=len(q), a_gap=co[1],
                         b_overlap=co[2], alpha_N=co[3], R2=r2))
f = pd.DataFrame(rows)
f.to_csv(os.path.join(d, "design_compare_fits.csv"), index=False)
print(f.round(3).to_string())
