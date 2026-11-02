"""Intermediate bond dimensions for the MPS-alone route.

The crossover gives the MPS-alone route the grid chi = 2, 4, ..., 64. Here we add chi = 3, 5, 6, 10, 12, 20,
24, 40, 48 for every crossover system (chains N = 12, 14, 16; cylinders 4x3, 4x4), with the same variational
MPS (trial_mps) and the exact ground state of results/cache_emp, and record the infidelity 1 - |<g|psi>|^2.
-> results/chi_mid.json   {system key: {chi: infidelity}}

    python sbpipe_study/run_chi_mid.py
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import run_empirical as RE
import run_extended as RX

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "chi_mid.json")
CHIS = [3, 5, 6, 10, 12, 20, 24, 40, 48]
SYSTEMS = [("chain", N, J2) for N in (12, 14, 16) for J2 in (0.0, 0.2411, 0.4)] + \
          [("cyl", L, J2) for L in (3, 4) for J2 in (0.0, 0.3)]


def _sys(s):
    kind, size, J2 = s
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    H = RX.model(kind, size, J2)
    z = np.load(os.path.join(RE.CACHEE, RE._act_key(kind, size, J2) + ".npz"))
    g = z["ground"]
    mpo = pauli_to_mpo(H)
    out = {}
    for chi in CHIS + RE.ACT_CHIS:
        if chi in RE.ACT_CHIS:                       # the grid trials as stored, for reference
            v = z[f"trial{chi}"]
        else:
            _, psi = trial_mps(mpo, chi)
            v = np.asarray(psi.to_dense()).reshape(-1)
        v = v / np.linalg.norm(v)
        out[chi] = max(0.0, 1.0 - float(abs(np.vdot(g, v)) ** 2))
    return RE._act_key(kind, size, J2), out


def main():
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    todo = [s for s in SYSTEMS if RE._act_key(*s) not in res]
    with ProcessPoolExecutor(2) as ex:
        for k, o in ex.map(_sys, todo):
            res[k] = o
            print(k, {c: f"{v:.1e}" for c, v in sorted(o.items())}, flush=True)
            json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
