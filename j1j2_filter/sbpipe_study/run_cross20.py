"""Implementable crossover on the 4x5 cylinder (20 qubits): planned filter with the size-transfer step against the
MPS alone, for every target of the 4x4 crossover (J2 = 0, 0.3; eps = 1e-2 ... 1e-5).

The 4x4 measured crossover needs the full spectrum (every S^z block diagonalised), which is out of reach at 20
qubits. The implementable route needs none of it: the filter is planned from DMRG energies (results/cache_t20)
and the trial overlap, and run at half the coarsest passing step of the 4x4 curve (same J2, eps, schedule).

    mps   : variational MPS chi = 2 ... 256 (powers of two and intermediates), infidelity to the exact ground state
            -> results/cross20_mps.json; trials chi = 2 ... 64 added to results/cache_t20/cyl5_J2*.npz
            (gamma for chi not in the DMRG inputs is the exact overlap; where both exist they differ by < 3e-5)
    cases : per (J2, eps, schedule) the trial chi (2 ... 64, 1 - gamma > eps, trial cheaper than the cheapest
            MPS alone under the study's uniform cost) whose certified plan is cheapest, then the planned circuit
            simulated at the transferred step (run_transfer20._case) -> results/cross20.jsonl

    python sbpipe_study/run_cross20.py mps
    python sbpipe_study/run_cross20.py
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

import run_aposteriori as AP
import run_extended as RX
import run_transfer20 as RT

RES = os.path.join(HERE, "results")
OUTM = os.path.join(RES, "cross20_mps.json")
OUT = os.path.join(RES, "cross20.jsonl")
J2S = (0.0, 0.3)
EPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)
TRIALS = (2, 4, 8, 16, 32, 64)
CHIS = sorted(set(TRIALS) | {3, 5, 6, 10, 12, 20, 24, 40, 48, 128, 256})


def _mps(J2):
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    key = RT._key("cyl", 5, J2)
    H = RX.model("cyl", 5, J2)
    mpo = pauli_to_mpo(H)
    z = dict(np.load(os.path.join(RT.CACHE, key + ".npz")))
    g = z["ground"]
    info = json.load(open(os.path.join(RT.CACHE, key + ".json")))
    inf = {}
    for chi in CHIS:
        _, psi = trial_mps(mpo, chi)
        v = np.asarray(psi.to_dense()).reshape(-1)
        v = v / np.linalg.norm(v)
        ov = float(abs(np.vdot(g, v)) ** 2)
        inf[chi] = max(0.0, 1.0 - ov)
        if chi in TRIALS and f"trial{chi}" not in z:
            z[f"trial{chi}"] = v
        if chi in TRIALS:
            info["gamma_exact"].setdefault(str(chi), ov)
            info["gamma"].setdefault(str(chi), ov)
        print(J2, chi, f"{inf[chi]:.2e}", flush=True)
    np.savez(os.path.join(RT.CACHE, key + ".npz"), **z)
    json.dump(info, open(os.path.join(RT.CACHE, key + ".json"), "w"))
    return J2, inf


def run_mps():
    res = {}
    with ProcessPoolExecutor(2) as ex:
        for J2, inf in ex.map(_mps, J2S):
            res[f"{J2:g}"] = inf
    json.dump(res, open(OUTM, "w"), indent=1)


def _plan(args):
    J2, eps, sched, chi = args
    from sbpipe import Ansatz, EnergyEstimates
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    from sbpipe.pipeline import Config, refine
    H = RX.model("cyl", 5, J2)
    info = json.load(open(os.path.join(RT.CACHE, RT._key("cyl", 5, J2) + ".json")))
    cost = mps_circuit_cost(H.num_qubits, chi)
    ans = Ansatz(gamma=info["gamma"][str(chi)], statevector=None, cost=cost, name=f"chi={chi}", meta=dict(chi=chi))
    try:
        plan = refine(H, EnergyEstimates(**info["est"]), ans, Config(eps=eps, schedule=sched, order=2))
        return args, plan.resources["expected_cnot"]
    except Exception as e:                                   # no plan reaches eps from this trial
        return args, None


def cases():
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    M = json.load(open(OUTM))
    P = RT._prev()
    todo, meta = [], {}
    for J2 in J2S:
        inf = {int(c): v for c, v in M[f"{J2:g}"].items()}
        for eps in EPS:
            alone = [mps_circuit_cost(20, c)["cnot"] for c in TRIALS if inf[c] <= eps]
            cap = min(alone) if alone else np.inf
            for s in ("auto", "qetu"):
                ch = [c for c in TRIALS if inf[c] > eps and mps_circuit_cost(20, c)["cnot"] < cap]
                r4 = P.get(("cyl", 4, J2, eps, s))
                dt = r4["curve"][AP.coarsest_pass(r4)]["dt"] if r4 else None
                meta[(J2, eps, s)] = (ch, dt)
                todo += [(J2, eps, s, c) for c in ch]
    return todo, meta


def main():
    todo, meta = cases()
    plans = {}
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 3))) as ex:
        for a, c in ex.map(_plan, todo):
            plans[a] = c
            print("plan", a, c, flush=True)
    done = {}
    if os.path.exists(OUT):
        for r in map(json.loads, open(OUT)):
            done[(r["J2"], r["eps"], r["schedule"], r["chi"])] = r
    T20 = {(r["J2"], r["eps"], r["schedule"], r["chi"]): r for r in map(json.loads, open(RT.OUT)) if r["kind"] == "cyl"}
    jobs = []
    for (J2, eps, s), (ch, dt) in meta.items():
        ok = [(plans[(J2, eps, s, c)], c) for c in ch if plans.get((J2, eps, s, c)) is not None]
        if not ok or dt is None:
            print("no case", J2, eps, s, ch, dt, flush=True)
            continue
        chi = min(ok)[1]
        k = (J2, eps, s, chi)
        if k in done:
            continue
        if k in T20:                                          # already simulated in the transfer check
            done[k] = dict(T20[k], cert_choice=dict(zip(map(str, ch), [plans.get((J2, eps, s, c)) for c in ch])))
            continue
        jobs.append(("cyl", 5, J2, eps, s, chi, dt))
    with open(OUT, "w") as fh:
        fh.write("".join(json.dumps(r) + "\n" for r in done.values()))
    jobs.sort(key=lambda c: -c[3])
    print(len(jobs), "circuits", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 3))) as ex, open(OUT, "a") as fh:
        for r in ex.map(RT._case, jobs):
            fh.write(json.dumps(r) + "\n")
            fh.flush()


if __name__ == "__main__":
    run_mps() if sys.argv[1:] == ["mps"] else main()
