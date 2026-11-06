"""MPS runs of the same study (needs quimb), for chains beyond exact diagonalisation.

Three modes:

    python run_mps.py search  --N 20 24 --gaps 0.3 1.0 --gammas 0.1 0.5
        measured cheapest filter per (gap, gamma, eps, protocol), same rules as run_maps.py
        -> <out_dir>/mps_N{N}.json, written after every point, resumes where it stopped

    python run_mps.py filter  --N 20 24 32 --gaps 0.3 --gammas 0.1 --K 5 --r 2 --n 1
        one fixed filter: p, 1 - F, CNOTs, RUS and AA cost, bond dimension, time (fast scaling cuts)

    python run_mps.py check   --N 8 10 12 --gaps 0.3 1.0 --gammas 0.1 0.5
        reruns the exact engine's chosen filters on the MPS and prints both (validation)

Reference state and truncation come from the config: --set mps_reference='"dmrg"' (or "ed", "auto"),
--set mps_max_bond=256, --set mps_cutoff=1e-12. Other config fields (J2, eps, n_grid, ...) apply as usual.
"""
import argparse
import json
import os
import time

from ornl import StudyConfig, chain_for_gap, costs, measure_point
from ornl import mps
from ornl.config import add_cli


def _row(mch):
    return dict(N=mch.N, gap=mch.gap, delta=mch.delta, W=mch.W, gap_scaled=mch.gap_s, E0=mch.E0, E1=mch.E1,
                Emax=mch.Emax, reference=mch.reference)


def search(a, cfg):
    os.makedirs(cfg.out_dir, exist_ok=True)
    for N in a.N:
        path = os.path.join(cfg.out_dir, f"mps_N{N}.json")
        data = json.load(open(path)) if os.path.exists(path) else dict(N=N, J1=cfg.J1, J2=cfg.J2, rows={})
        for g in a.gaps:
            mch = None
            for gam in a.gammas:
                key = f"{g:.6g}/{gam:.6g}"
                if key in data["rows"]:
                    continue
                t0 = time.time()
                mch = mch or mps.mps_chain(N, g, cfg)
                rows, ek = mps.measure_point(mch, [gam], a.eps, cfg, verbose=a.verbose)
                data["rows"][key] = dict(_row(mch), gamma=gam, eps_knob=ek[0], seconds=time.time() - t0,
                                         max_bond=cfg.mps_max_bond, cutoff=cfg.mps_cutoff,
                                         best={f"{e:g}": v for e, v in rows[0].items()})
                json.dump(data, open(path + ".tmp", "w"), indent=1)
                os.replace(path + ".tmp", path)
                b = rows[0][a.eps[0]]
                print(f"N={N} gap {mch.gap:.4f} gamma {gam}: " +
                      " ".join(f"{pr} K{v['K']} r{v['r']} n{v['n']} cost {v['cost']:.3g}" for pr, v in b.items() if v)
                      + f"  ({time.time() - t0:.0f} s)", flush=True)


def fixed(a, cfg):
    print(f"{'N':>3} {'gap':>6} {'gamma':>6} {'p':>9} {'1-F':>9} {'CNOT':>9} {'RUS':>10} {'AA':>10} {'chi':>4} {'s':>6}")
    for N in a.N:
        for g in a.gaps:
            mch = mps.mps_chain(N, g, cfg)
            for gam in a.gammas:
                t0 = time.time()
                ek = mch.eps_for_overlap(gam)
                v, p, inf = mch.filter(mch.trial(ek), a.K, a.r, a.n)
                c = costs.cnot_per_attempt(N, a.K, a.r, a.n, cfg.bond_cnots)
                refl = costs.reflection_cnots(N, a.K, a.r, cfg.mcz_cnots_per_qubit)
                print(f"{N:3d} {mch.gap:6.3f} {gam:6.3f} {p:9.6f} {inf:9.2e} {c:9d} {costs.cost(c, p, 'rus'):10.4g} "
                      f"{costs.cost(c, p, 'aa', refl):10.4g} {v.max_bond():4d} {time.time() - t0:6.0f}", flush=True)


def check(a, cfg):
    print(f"{'N':>3} {'gap':>5} {'gam':>5} {'eps':>6} {'K':>2} {'r':>2} {'n':>2} {'p exact':>10} {'p MPS':>10} "
          f"{'1-F exact':>10} {'1-F MPS':>10}")
    worst = [0.0, 0.0]
    for N in a.N:
        for g in a.gaps:
            ch = chain_for_gap(N, g, cfg)
            mch = mps.mps_chain(N, g, cfg, delta=ch.delta)
            rows, ek = measure_point(ch, a.gammas, a.eps, cfg)
            for gam, row, e in zip(a.gammas, rows, ek):
                if e is None:
                    continue
                psi = mch.trial(e)
                for eps in a.eps:
                    b = row[eps]["rus"]
                    _, p, inf = mch.filter(psi, b["K"], b["r"], b["n"])
                    worst = [max(worst[0], abs(p - b["p"])), max(worst[1], abs(inf - b["infid"]))]
                    print(f"{N:3d} {g:5.2f} {gam:5.2f} {eps:6.0e} {b['K']:2d} {b['r']:2d} {b['n']:2d} {b['p']:10.7f} "
                          f"{p:10.7f} {b['infid']:10.3e} {inf:10.3e}", flush=True)
    print(f"largest differences: p {worst[0]:.1e}, 1-F {worst[1]:.1e} (reference: {cfg.mps_reference})")


if __name__ == "__main__":
    ap = add_cli(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter))
    ap.add_argument("mode", choices=["search", "filter", "check"])
    ap.add_argument("--N", type=int, nargs="+", default=None, help="default: cfg.mps_Ns (check: 8 10 12)")
    ap.add_argument("--gaps", type=float, nargs="+", default=[0.3, 1.0])
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.1, 0.5])
    ap.add_argument("--eps", type=float, nargs="+", default=None, help="default: cfg.eps")
    ap.add_argument("--K", type=int, default=5)
    ap.add_argument("--r", type=int, default=2)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    cfg = StudyConfig.from_cli(a.config, a.set)
    a.N = a.N or ([8, 10, 12] if a.mode == "check" else list(cfg.mps_Ns))
    a.eps = a.eps or list(cfg.eps)
    {"search": search, "filter": fixed, "check": check}[a.mode](a, cfg)
