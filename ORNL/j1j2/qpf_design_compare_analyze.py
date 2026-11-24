"""Analyse j1j2/qpf_design_compare.py runs against the manuscript's published-QPF points.

    python j1j2/qpf_design_compare_analyze.py --raw DIR --manuscript MANUSCRIPT/qpf_maps/qpf_points.csv --out DIR

Writes design_compare_points.csv (one row per point and schedule, analytic gate model), design_compare_summary.csv
(ratios to the authors' v4 and to the manuscript's 'best' published schedule), and figures.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KEY = ["N", "gap_target", "knob", "gamma_target", "eps", "route", "protocol"]
MM_CAP = ["minimax c=1", "minimax c=1.25", "minimax c=1.5", "minimax c=2"]
PH = ["minimax+ph p0=0.9 c=1", "minimax+ph p0=0.9 c=1.25", "minimax+ph p0=0.9 c=1.5", "minimax+ph p0=0.7 c=1",
      "minimax+ph p0=0.7 c=1.25"]
LM = ["logmean1e-4 c=1", "logmean1e-4 c=1.25", "logmean1e-4 c=1.5", "logmean1e-4 c=2"]
V4C = ["QPF (t_i) v4", "v4 c=1.25", "v4 c=1.5", "v4 c=2"]


def load(raw):
    rows, designs = [], []
    for f in sorted(glob.glob(os.path.join(raw, "cmp_N*_gap*.json")) + glob.glob(os.path.join(raw, "cmpph_N*_gap*.json")) + glob.glob(os.path.join(raw, "cmplog_N*_gap*.json"))
                    + glob.glob(os.path.join(raw, "cmpv4nc_N*_gap*.json"))):
        r = json.load(open(f))
        for x in r["results"]:
            if not x.get("reached"):
                continue
            rows.append(dict(N=r["N"], gap_target=r["gap_target"], W=r["W"], **{k: x.get(k) for k in (
                "knob", "gamma_target", "gamma", "eps", "route", "gap_used", "schedule", "protocol", "cx_per_success",
                "rz_per_success", "cx_attempt", "p", "infid", "R", "m_E", "n_E", "k_J", "time_E_scaled",
                "trotter_steps_E", "cx_E_block", "cx_J", "k")}))
        for k, d in r["designs"].items():
            designs.append(dict(N=r["N"], gap_target=r["gap_target"], key=k, **d))
    # one row per point and schedule (a point appears in both cmp_ and cmpph_ files with different schedules)
    return pd.DataFrame(rows), pd.DataFrame(designs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--manuscript", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    df, des = load(a.raw)
    df.to_csv(os.path.join(a.out, "design_compare_points.csv"), index=False)
    des.to_csv(os.path.join(a.out, "design_compare_designs.csv"), index=False)

    ms = pd.read_csv(a.manuscript, low_memory=False)
    ms = ms[ms.gate_model == "analytic"]
    piv = df.pivot_table(index=KEY, columns="schedule", values="cx_per_success", aggfunc="first")
    msp = ms.pivot_table(index=KEY, columns="schedule", values="cx_per_success", aggfunc="first")
    piv = piv.join(msp.rename(columns={"QPF (t_i)": "ms v4", "QPF (t_i, δ_i)": "ms v3", "best": "ms best"}), how="left")
    piv["minimax best c"] = piv[MM_CAP].min(axis=1)
    piv["minimax any"] = piv[MM_CAP + ["minimax c=1 nocap"]].min(axis=1)
    piv["which c"] = piv[MM_CAP + ["minimax c=1 nocap"]].idxmin(axis=1)
    ph = [c for c in PH if c in piv]
    if ph:
        piv["minimax+ph best"] = piv[ph].min(axis=1)
        piv["all options"] = piv[MM_CAP + ["minimax c=1 nocap"] + ph].min(axis=1)
    if all(c in piv for c in LM):
        piv["logmean best c"] = piv[LM].min(axis=1)
    if all(c in piv for c in V4C):
        piv["v4 best c"] = piv[V4C].min(axis=1)
    piv = piv.reset_index()
    piv.to_csv(os.path.join(a.out, "design_compare_wide.csv"), index=False)

    out = []
    for (route, proto, eps), g in piv.groupby(["route", "protocol", "eps"]):
        rep = (g["QPF (t_i) v4"] / g["ms v4"]).dropna()
        extra = [c for c in PH + ["minimax+ph best", "all options"] + LM + ["logmean best c"] + V4C[1:] + ["v4 best c", "v4 nocap c=1"]
                 if c in g]
        for col in MM_CAP + ["minimax c=1 nocap", "minimax best c", "minimax any"] + extra:
            for base in ("QPF (t_i) v4", "ms best"):
                r = (g[col] / g[base]).replace([np.inf, -np.inf], np.nan).dropna()
                if not len(r):
                    continue
                out.append(dict(route=route, protocol=proto, eps=eps, option=col, baseline=base, n=len(r),
                                median=r.median(), q10=r.quantile(0.1), q90=r.quantile(0.9), min=r.min(), max=r.max(),
                                frac_cheaper=(r < 1).mean(), reproduce_v4_median=rep.median(),
                                reproduce_v4_maxdev=(rep - 1).abs().max() if len(rep) else np.nan))
    s = pd.DataFrame(out)
    s.to_csv(os.path.join(a.out, "design_compare_summary.csv"), index=False)
    pd.set_option("display.width", 250)
    show = ["minimax c=1", "minimax best c", "minimax any", "logmean1e-4 c=1", "logmean best c", "v4 best c", "v4 nocap c=1", "minimax c=1 nocap"]
    print(s[s.option.isin(show)][["route", "protocol", "eps", "option", "baseline", "n", "median", "q10", "q90", "frac_cheaper"]].round(3).to_string())
    print(piv.groupby(["route", "protocol", "eps"])["which c"].value_counts().unstack().fillna(0).astype(int))


if __name__ == "__main__":
    main()
