"""Resources vs target error eps for three QPF energy-pulse designs on the uniform J1-J2 chain.

    python j1j2/qpf_design_eps_scan.py --N 8 10 12 --J2 0.3 --out DIR

Chain: open, uniform (delta = 0), J1 = 1, J2 given, S^z = 0 sector, exact spectrum.
Trial: Neel state (S^z = 0 sector), or with --trial mps the mps_to_circuit statevector of the DMRG MPS
(j1j2/mps_trial_gen.py, L = 1 layer), which leaks slightly out of S^z = 0, so that case runs in the full 2^N space.
Route: energy filter only ("H only"), designed for the global gap E1 - E0 scaled by W.
Designs (zero phases, prod_i cos(t_i H)):
  v4                    the authors' design (qpf.design, T = pi / gap, cap T/3, m by design cost)
  v5 same T no cap      qpf.design_plus minimax, T = pi / gap, no cap
  v5 1.5T no cap        qpf.design_plus minimax, T = 1.5 pi / gap, no cap
Circuit: the manuscript's pipeline (j1j2.qpf_maps): second-order Trotter, one dt set by the longest pulse
(n_E steps on it), block repeated R <= R_MAX times, analytic 6-CNOT gate model. For each eps the cheapest (n_E, R)
in CNOT per success (RUS) with 1 - F <= eps there and one Trotter step finer is kept.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from j1j2 import qpf_maps as qm  # noqa: E402
from j1j2.route_costs import pick as authors_pick  # noqa: E402
from qpf import design_plus  # noqa: E402

EPS = (3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4)


def designs(gap):
    a = authors_pick(gap, "v4")
    b = design_plus.design(gap, "minimax", 1.0, cap_frac=1.0)
    c = design_plus.design(gap, "minimax", 1.5, cap_frac=1.0)
    return {"v4": a["times"], "v5 same T no cap": b["times"], "v5 1.5T no cap": c["times"]}


def neel(N, idx):
    s = sum(1 << i for i in range(1, N, 2))
    v = np.zeros(len(idx), complex)
    v[np.searchsorted(idx, s)] = 1.0
    return v


def mps_vec(N, J2, idx, mps_dir, layers=1):
    z = np.load(os.path.join(mps_dir, f"mps_trial_N{N}_J2{J2:g}_L{layers}.npz"))
    v = z["psi"][idx]
    return v / np.linalg.norm(v), int(z["cx"])


def run(N, J2, trial="neel", mps_dir=None):
    qm.J2 = J2
    qm.EPS = EPS
    qm.R_MAX = 8
    qm.robust_counts = lambda *a: (0, 0)
    idx = qm.sector(N) if trial == "neel" else np.arange(1 << N)
    S2 = qm.matrix(qm.as_sector(qm.trotter.total_spin_sq(N)), idx)
    ch = qm.Chain(N, 0.0, idx, S2).spectrum()
    W = ch.W
    gap = ch.gap_global / W
    v0, trial_cx = (neel(N, idx), 0) if trial == "neel" else mps_vec(N, J2, idx, mps_dir)
    gamma = float(abs(np.vdot(ch.g0, v0)) ** 2)
    Hs = qm.SectorSwapSum(N, ch.Hop.pairs, ch.Hop.weights / W, (ch.Hop.const - ch.E0) / W)
    B = len(ch.Hop.pairs)
    out = dict(N=N, J2=J2, trial=trial, trial_cx=trial_cx, W=W, gap=ch.gap_global, gap_singlet=ch.gap_singlet, gap_scaled=gap, gamma=gamma,
               bonds=B, cx_per_step=6 * (2 * B - 2), designs={}, rows=[])
    for name, t in designs(gap).items():
        t = np.sort(np.asarray(t)[np.asarray(t) > 1e-6])
        E = qm.Stage(Hs, ("H", N, J2), t, np.zeros(len(t)), qm.N_E)
        tab = qm.evaluate_energy(E, [(v0[:, None], np.ones(1))], ch.g0, idx)
        out["designs"][name] = dict(times=t.tolist(), m=len(t), T=float(t.sum()),
                                    max_leak=design_plus.max_leak(gap, t))
        for eps in EPS:
            best = None
            for (gi, gj, R), (inf, P) in tab.items():
                if not qm.valid(tab, gi, gj, R, 0, eps, 1, False):
                    continue
                c = E.counts(gj)
                cx_att = R * c["cx_an"] + trial_cx
                cand = dict(eps=eps, design=name, n_E=E.grid[gj]["n_longest"], ns=E.grid[gj]["ns"], R=R,
                            steps_block=c["steps"], steps=R * c["steps"], cx_block=c["cx_an"], cx_attempt=cx_att,
                            p=float(P[0]), infid=float(inf[0]), cx_per_success=cx_att / float(P[0]),
                            time_scaled=R * c["time"])
                if best is None or cand["cx_per_success"] < best["cx_per_success"]:
                    best = cand
            out["rows"].append(best or dict(eps=eps, design=name, reached=False))
        print(f"N={N} {name}: m={len(t)} T={t.sum():.2f} max f2={out['designs'][name]['max_leak']:.1e}", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12])
    ap.add_argument("--J2", type=float, default=0.3)
    ap.add_argument("--J2s", type=float, nargs="+", help="several J2 values (overrides --J2)")
    ap.add_argument("--trial", choices=("neel", "mps"), default="neel")
    ap.add_argument("--mps-dir")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for J2 in (a.J2s or [a.J2]):
        for N in a.N:
            r = run(N, J2, a.trial, a.mps_dir)
            tag = "" if a.trial == "neel" else "_mps"
            json.dump(qm._clean(r), open(os.path.join(a.out, f"eps_scan{tag}_N{N}_J2{J2:g}.json"), "w"), indent=1)
