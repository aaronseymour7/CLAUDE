"""Measured filter-cost maps over (gap, overlap) for every N and eps in the config.

    python run_maps.py                                     # defaults in ornl/config.py (N = 16 takes ~2 h)
    python run_maps.py --set Ns=[8,10] --set eps=[0.01]    # a subset
    python run_maps.py --config my_study.json              # any subset of StudyConfig fields

One JSON per N in cfg.out_dir, written after every gap value, so an interrupted run resumes where it stopped.
"""
import argparse
import json
import os
import time

from ornl import StudyConfig, chain_for_gap, measure_point
from ornl.config import add_cli


def run(N, cfg):
    os.makedirs(cfg.out_dir, exist_ok=True)
    path = os.path.join(cfg.out_dir, f"map_N{N}.json")
    grid = dict(N=N, J1=cfg.J1, J2=cfg.J2, gaps=list(cfg.gaps), gammas=list(cfg.gammas), eps=list(cfg.eps),
                trial=cfg.trial)
    data = json.load(open(path)) if os.path.exists(path) else dict(grid, rows={})
    for k in ("J1", "J2", "gammas", "eps", "trial"):
        if data.get(k, "exact" if k == "trial" else None) != grid[k]:
            raise SystemExit(f"{path} was made with another {k}; use another out_dir or move it away")
    t0 = time.time()
    for g in cfg.gaps:
        key = f"{g:.6g}"
        if key in data["rows"]:
            continue
        ch = chain_for_gap(N, g, cfg)
        rows, eps_knob = measure_point(ch, cfg.gammas, cfg.eps, cfg)
        data["rows"][key] = dict(
            gap=ch.gap, delta=ch.delta, W=ch.W, gamma0=ch.gamma0, prep_cnot=ch.prep_cnot, gap_scaled=ch.gap_s, E0=ch.E0, E1=ch.E1, Emax=ch.Emax,
            eps_knob=eps_knob, gamma=[None if e is None else ch.overlap(e) for e in eps_knob],
            best=[{f"{e:g}": v for e, v in r.items()} for r in rows])
        json.dump(data, open(path + ".tmp", "w"))
        os.replace(path + ".tmp", path)
        print(f"N={N:2d} gap {ch.gap:.4f} (delta {ch.delta:+.4f}): {time.time() - t0:.0f} s", flush=True)
    return data


if __name__ == "__main__":
    a = add_cli(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)).parse_args()
    cfg = StudyConfig.from_cli(a.config, a.set)
    for N in cfg.Ns:
        run(N, cfg)
