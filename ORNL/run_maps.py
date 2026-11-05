"""Measured filter-cost maps over (gap, overlap) for several N and target infidelities.

    python run_maps.py                       # N = 8 10 12 14 16, eps = 1e-2 1e-3 (N = 16 takes a few hours)
    python run_maps.py --N 8 10 --eps 1e-2   # a subset

One JSON per N in results/, written after every gap row, so an interrupted run resumes where it stopped.
Each row holds, per target overlap, per eps and per protocol (rus / aa), the cheapest measured filter.
"""
import argparse
import json
import os
import time

import numpy as np

import filtercost as fc

GAPS = np.geomspace(0.005, 1.7, 14)                    # physical gap E1 - E0, units of J1
GAMMAS = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]


def run(N, eps_list, out_dir="results"):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"map_N{N}.json")
    data = json.load(open(path)) if os.path.exists(path) else \
        dict(N=N, J1=fc.J1, J2=fc.J2, gaps=GAPS.tolist(), gammas=GAMMAS, eps=eps_list, rows={})
    assert data["gammas"] == GAMMAS and data["eps"] == eps_list, f"{path} was made with another grid; move it away"
    t0 = time.time()
    for g in GAPS:
        key = f"{g:.6g}"
        if key in data["rows"]:
            continue
        d = fc.delta_for_gap(N, g)
        ch = fc.Chain(N, d)
        rows, eps_knob = fc.measure_point(ch, GAMMAS, eps_list)
        data["rows"][key] = dict(
            gap=ch.gap, delta=d, W=ch.W, gap_scaled=ch.gap_s, E0=ch.E0, E1=ch.E1, Emax=ch.Emax,
            eps_knob=eps_knob,
            gamma=[None if e is None else ch.overlap(e) for e in eps_knob],
            best=[{f"{e:g}": v for e, v in r.items()} for r in rows])
        json.dump(data, open(path + ".tmp", "w"))
        os.replace(path + ".tmp", path)
        print(f"N={N:2d} gap {ch.gap:.4f} (delta {d:+.4f}): {time.time() - t0:.0f} s", flush=True)
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12, 14, 16])
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-2, 1e-3])
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    for N in a.N:
        run(N, a.eps, a.out)
