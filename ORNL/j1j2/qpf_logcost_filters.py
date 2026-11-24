"""Filter-level comparison of cost functions for the QPF energy pulses (all prod_i cos(t_i E), zero phases):
  v4           the authors' |1 - 1/sqrt(sum_n f^2)|, their SLSQP (qpf.design via route_costs.pick)
  log v4       log(sum_n f^2) with the authors' SLSQP, same start, bounds and constraints (same minimiser as v4)
  log-mean     mean_E log(f^2 + floor), floor 1e-4 / 1e-8 / 1e-12 (qpf.design_plus, objective "logmean")
  minimax      log max_E f^2 (qpf.design_plus, objective "minimax")
at T = pi / gap and 1.5 pi / gap, with the authors' cap t_i <= T / 3.
Usage: python j1j2/qpf_logcost_filters.py OUT_DIR"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from j1j2.route_costs import pick, M_RANGE  # noqa: E402
from qpf import design, design_plus as dp  # noqa: E402


def authors_log_v4(gap, m, T):
    """The authors' SLSQP setup (halving start, bounds, |sum t - T| equality) on log(sum_n f^2)."""
    E = design.baroni_spectrum(gap)
    fun = lambda t: np.log(np.prod(np.cos(E[:, None] * t) ** 2, axis=1).sum())
    cons = [{"type": "eq", "fun": lambda t: abs(np.sum(t) - T)}]
    r = minimize(fun, design.halving_guess(T, m), method="SLSQP", bounds=[(0, T / 3)] * m, constraints=cons,
                 options=dict(maxiter=5000, ftol=1e-12), tol=1e-13)
    return r.x


def pick_log_v4(gap, T):
    best = None
    for m in M_RANGE:
        t = authors_log_v4(gap, m, T)
        s = dp.mean_leak(gap, t)
        if best is None or s < best[0] - 1e-15:
            best = (s, np.sort(t))
    return best[1]


def run(gaps=(0.1, 0.2, 0.3), cs=(1.0, 1.5)):
    rows, designs = [], {}
    for g in gaps:
        for c in cs:
            T = c * np.pi / g
            sets = {}
            if c == 1.0:
                sets["v4 (authors)"] = pick(g, "v4")["times"]
            else:                                            # authors' SLSQP at the longer T
                best = None
                for m in M_RANGE:
                    q = design.optimize(g, m, T=T, method="v4")
                    if best is None or q["cost"] < best["cost"] - 1e-9:
                        best = q
                sets["v4 (authors)"] = np.sort(best["times"])
            sets["log v4"] = pick_log_v4(g, T)
            for fl in (1e-4, 1e-8, 1e-12):
                sets[f"log-mean, floor {fl:g}"] = dp.design(g, "logmean", c, floor=fl, n_random=10)["times"]
            sets["minimax (log max)"] = dp.design(g, "minimax", c)["times"]
            for k, t in sets.items():
                t = np.asarray(t)
                Ef = np.linspace(g, 1, 20001)
                f2 = np.prod(np.cos(np.outer(Ef, t)) ** 2, axis=1)
                rows.append(dict(gap=g, c=c, cost=k, m=int(np.sum(t > 1e-3 * T)), max_leak=f2.max(),
                                 f2_at_gap=f2[0], mean_leak=dp.mean_leak(g, t),
                                 geo_mean_leak=float(np.exp(np.mean(np.log(f2 + 1e-300)))),
                                 frac_above_1e3=float(np.mean(f2 > 1e-3)), times=np.round(t, 4).tolist()))
                designs[(g, c, k)] = t
    return pd.DataFrame(rows), designs


if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    df, designs = run()
    df.to_csv(os.path.join(out, "logcost_filters.csv"), index=False)
    pd.set_option("display.width", 200)
    print(df.drop(columns="times").to_string(float_format=lambda x: f"{x:.2e}"))
    np.save(os.path.join(out, "logcost_designs.npy"), {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in designs.items()})
