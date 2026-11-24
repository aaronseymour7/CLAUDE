"""Published QPF with the authors' v4 energy pulses vs the improved design option (qpf.design_plus), on the
Section A grid (manuscript/qpf_maps): same chain, gaps, trials, J^2 ladder, Trotter search, R <= 4, RUS and AA.

    python j1j2/qpf_design_compare.py --N 8 --out DIR

Only the energy-pulse times change. Everything else is j1j2.qpf_maps.run_point unchanged; this script swaps its
schedule list and its design picker, and skips the Qiskit transpile of each pulse (the analytic 6-CNOT gate model,
the manuscript's primary one, does not need it). Rows with gate_model == "transpile" are therefore empty and are
dropped by the analysis. AA reflection counts come from Qiskit's mcx, which changes between
versions: use Qiskit 2.6.2 to match the manuscript (2.0.3 gives about 1.7x more reflection CNOTs).

Schedules (all zero phases, prod_i cos(t_i H)):
  QPF (t_i) v4             the authors' design, as in the manuscript (qpf.design.optimize via route_costs.pick)
  minimax c=X              qpf.design_plus.design(objective="minimax", c=X), authors' cap t_i <= T / 3
  minimax c=1 nocap        same at the authors' T = pi / gap, without the T / 3 cap
  --set phases:
  minimax+ph p0=P c=X      qpf.design_plus.design_phases: minimax of f(E)^2 / f(0)^2 with phases, f(0)^2 >= P
  --set log:
  logmean1e-4 c=X          qpf.design_plus.design(objective="logmean", floor=1e-4): mean_E log(f^2 + 1e-4)
  v4 c=X                   the authors' v4 design and SLSQP at the longer T = c pi / gap (time knob alone)
  --set v4nocap:
  v4 nocap c=1             the v4 objective (multistart, qpf.design_plus objective "v4") at T = pi / gap, no T/3 cap
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import qpf_maps as qm  # noqa: E402
from j1j2.route_costs import pick as authors_pick  # noqa: E402
from qpf import design_plus  # noqa: E402

SCHEDS = {
    "QPF (t_i) v4": "v4",
    "minimax c=1": "mm|1.0|cap",
    "minimax c=1.25": "mm|1.25|cap",
    "minimax c=1.5": "mm|1.5|cap",
    "minimax c=2": "mm|2.0|cap",
    "minimax c=1 nocap": "mm|1.0|nocap",
}
SCHEDS_PHASES = {
    "minimax+ph p0=0.9 c=1": "ph|1.0|0.9",
    "minimax+ph p0=0.9 c=1.25": "ph|1.25|0.9",
    "minimax+ph p0=0.9 c=1.5": "ph|1.5|0.9",
    "minimax+ph p0=0.7 c=1": "ph|1.0|0.7",
    "minimax+ph p0=0.7 c=1.25": "ph|1.25|0.7",
}
SCHEDS_LOG = {
    "logmean1e-4 c=1": "lm|1.0|1e-4",
    "logmean1e-4 c=1.25": "lm|1.25|1e-4",
    "logmean1e-4 c=1.5": "lm|1.5|1e-4",
    "logmean1e-4 c=2": "lm|2.0|1e-4",
    "v4 c=1.25": "v4c|1.25",
    "v4 c=1.5": "v4c|1.5",
    "v4 c=2": "v4c|2.0",
}
SCHEDS_V4NOCAP = {"v4 nocap c=1": "v4nc|1.0"}
_CACHE = {}


def authors_pick_T(gap, c):
    """The authors' v4 design (their SLSQP, cap T/3, m = 4..12 by design cost) at T = c pi / gap."""
    import numpy as np
    from qpf import design
    best = None
    for m in range(4, 13):
        q = design.optimize(gap, m, T=c * np.pi / gap, method="v4")
        if best is None or q["cost"] < best["cost"] - 1e-9:
            best = dict(q, m=m)
    o = np.argsort(best["times"])
    return dict(times=best["times"][o], phases=best["phases"][o], m=best["m"], cost=best["cost"], f0=1.0)


def pick(gap, method):
    key = (round(gap, 12), method)
    if key not in _CACHE:
        if method == "v4":
            _CACHE[key] = authors_pick(gap, "v4")
        elif method.startswith("lm|"):
            _, c, fl = method.split("|")
            d = design_plus.design(gap, "logmean", float(c), floor=float(fl), n_random=10)
            _CACHE[key] = dict(times=d["times"], phases=d["phases"], m=d["m"], cost=d["logmean_leak"], f0=1.0)
        elif method.startswith("v4nc|"):
            d = design_plus.design(gap, "v4", float(method.split("|")[1]), cap_frac=1.0)
            _CACHE[key] = dict(times=d["times"], phases=d["phases"], m=d["m"], cost=d["v4_cost"], f0=1.0)
        elif method.startswith("v4c|"):
            _CACHE[key] = authors_pick_T(gap, float(method.split("|")[1]))
        elif method.startswith("ph|"):
            _, c, p0 = method.split("|")
            d = design_plus.design_phases(gap, float(c), float(p0))
            _CACHE[key] = dict(times=d["times"], phases=d["phases"], m=d["m"], cost=d["max_leak"], f0=d["f0"])
        else:
            _, c, cap = method.split("|")
            d = design_plus.design(gap, "minimax", float(c), cap_frac=1 / 3 if cap == "cap" else 1.0)
            _CACHE[key] = dict(times=d["times"], phases=d["phases"], m=d["m"], cost=d["max_leak"], f0=1.0)
    return _CACHE[key]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--gaps", type=float, nargs="+", default=list(qm.GAPS))
    ap.add_argument("--grid-dir", required=True, help="directory with the manuscript's grid_N*.json (delta per gap)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", choices=("times", "phases", "log", "v4nocap"), default="times")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    qm.SCHEDS = {"times": SCHEDS, "phases": SCHEDS_PHASES, "log": SCHEDS_LOG, "v4nocap": SCHEDS_V4NOCAP}[a.set]
    prefix = {"times": "cmp", "phases": "cmpph", "log": "cmplog", "v4nocap": "cmpv4nc"}[a.set]
    qm.pick = pick
    qm.robust_counts = lambda *args: (0, 0)                # analytic model only
    deltas = json.load(open(os.path.join(a.grid_dir, f"grid_N{a.N}.json")))["deltas"]
    idx = qm.sector(a.N)
    S2 = qm.matrix(qm.as_sector(qm.trotter.total_spin_sq(a.N)), idx)
    t0 = time.time()
    for gt in a.gaps:
        dl = deltas.get(f"{gt:g}")
        path = os.path.join(a.out, f"{prefix}_N{a.N}_gap{gt:g}.json")
        if dl is None or os.path.exists(path):
            continue
        r = qm.run_point(a.N, dl, list(qm.GAMMAS), list(qm.KNOBS), idx, S2, log=lambda s: print(s, flush=True))
        r["results"] = [x for x in r["results"] if x["gate_model"] == "analytic"]
        r.update(gap_target=gt, J2=qm.J2, eps_list=list(qm.EPS), K_J=qm.K_J, N_E=qm.N_E, R_MAX=qm.R_MAX)
        json.dump(qm._clean(r), open(path + ".tmp", "w"))
        os.replace(path + ".tmp", path)
    print(f"N={a.N} done [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
