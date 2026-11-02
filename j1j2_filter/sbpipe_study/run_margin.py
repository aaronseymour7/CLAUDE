"""Gap margin from DMRG convergence data, tested against exact gaps.

The certified plans are planned from a DMRG gap Delta(chi) = E1(chi) - E0(chi) (E1 = lowest S^z = 1 state), and
an overestimated gap is the input error that breaks them. A rule that uses only DMRG data:

    Delta_plan = Delta(chi) - c * |Delta(chi) - Delta(chi/2)|        (c = 1, 2, 4)

i.e. a margin set by the last change of the gap along the bond-dimension ladder. For every system with an exact
gap (chains N = 12, 14, 16 and cylinders 4x3, 4x4, 4x5, all J2) we run the DMRG ladder chi = 4 ... 128 for the
ground state and the S^z = 1 state, and record at each final chi whether Delta_plan <= Delta_exact (the margin
covers the error) and f = Delta_plan / Delta_exact (how much gap the plan gives away; f < 1 is safe).
-> results/margin.json

    python sbpipe_study/run_margin.py
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
import time
from concurrent.futures import ProcessPoolExecutor

import run_extended as RX
import run_transfer20 as RT

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "margin.json")
LADDER = (4, 8, 16, 32, 64, 128)


def systems():
    S = []
    for r in map(json.loads, open(os.path.join(RES, "aposteriori.jsonl"))):
        k = (r["kind"], r["size"], r["J2"])
        if k not in [s[:3] for s in S]:
            S.append(k + (r["gap"],))
    for J2 in (0.0, 0.3):
        E = json.load(open(os.path.join(RT.CACHE, RT._key("cyl", 5, J2) + ".json")))["E01"]
        S.append(("cyl", 5, J2, E[1] - E[0]))
    return S


def _one(s):
    from sbpipe.frontends import dmrg_quimb as DQ
    kind, size, J2, gap = s
    t0 = time.time()
    H = RX.model(kind, size, J2)
    mpo = DQ.pauli_to_mpo(H)
    g = DQ.ground_state(mpo, bond_dims=LADDER)
    e = DQ.sector_ground_state(mpo, DQ.total_sz(H.num_qubits), 1.0, 2.0, bond_dims=LADDER)
    return dict(kind=kind, size=size, J2=J2, N=H.num_qubits, gap_exact=gap,
                E0=[h[1] for h in g["history"]], E1=[h[1] for h in e["history"]], chis=list(LADDER),
                time=time.time() - t0)


def main():
    res = json.load(open(OUT)) if os.path.exists(OUT) else []
    have = {(r["kind"], r["size"], r["J2"]) for r in res}
    todo = [s for s in systems() if s[:3] not in have]
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 2))) as ex:
        for r in ex.map(_one, todo):
            res.append(r)
            print(r["kind"], r["size"], r["J2"], f"{r['time']:.0f}s", flush=True)
            json.dump(res, open(OUT, "w"), indent=1)


if __name__ == "__main__":
    main()
