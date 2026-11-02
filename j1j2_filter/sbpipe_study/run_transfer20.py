"""Size-transfer check at 20 qubits: chain N = 20 (from N = 16) and cylinder 4x5 (from 4x4).

The larger-cylinder estimates take the Trotter step from the 4x4 a-posteriori curves (half the coarsest
passing step) without simulating the larger circuit. Here we simulate that rule one size further than any
earlier run. For each J2, eps in {1e-2, 1e-3, 1e-4} and schedule (auto = SBC, qetu) we take the trial chi
used at the smaller size, plan the filter at 20 qubits from DMRG inputs (refine, as in run_aposteriori), and
run the Trotterised planned circuit exactly (statevector, fastsim.Strang) at

    dt_prev       the coarsest passing step of the smaller system (does the step shrink with N?)
    dt_prev / 2   the transferred step (the rule used in the paper)
    dt_prev / 4, dt_prev / 8 only if the previous point failed,

recording infidelity to the exact ground state, success probability and CNOT counts.
-> results/transfer20.jsonl

    python sbpipe_study/run_transfer20.py inputs     # DMRG / exact ground states, cached in results/cache_t20
    python sbpipe_study/run_transfer20.py            # the circuits
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict

import numpy as np

from sbpipe import Ansatz, EnergyEstimates
from sbpipe.pipeline import Config, refine
from sbpipe.resources import filter_counts

import fastsim as FS
import run_aposteriori as AP
import run_extended as RX

RES = os.path.join(HERE, "results")
CACHE = os.path.join(RES, "cache_t20")
OUT = os.path.join(RES, "transfer20.jsonl")
EPS = (1e-2, 1e-3, 1e-4)
SYSTEMS = [("chain", 20, J2, ("chain", 16)) for J2 in (0.0, 0.2411, 0.4)] + \
          [("cyl", 5, J2, ("cyl", 4)) for J2 in (0.0, 0.3)]
os.makedirs(CACHE, exist_ok=True)


def _prev():
    rows = [json.loads(l) for l in open(AP.OUT)]
    for r in rows:
        r["curve"] = AP.dedupe(r["curve"])
    return {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}


def _key(kind, size, J2):
    return f"{kind}{size}_J2{J2:g}"


def _inputs(args):
    kind, size, J2, chis = args
    import scipy.sparse.linalg as sla
    from sbpipe.frontends import dmrg_quimb as DQ
    t0 = time.time()
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    mpo = DQ.pauli_to_mpo(H)
    if kind == "chain":
        rec = json.load(open(os.path.join(RX.CACHE, RX.RS._key(size, J2) + ".json")))
        est = rec["est"]
        gam = {int(c): v["gamma"] for c, v in rec["ansatz"].items()}
    else:
        DQ.ansatz_from_dmrg.__defaults__ = (0,)     # no dense statevector inside; built below
        e, ans, _ = DQ.dmrg_inputs(H, chi_trial=list(chis), bond_dims=(16, 32, 64, 128, 256))
        est = asdict(e)
        gam = {a.meta["chi"]: a.gamma for a in ans}
    arrays = {}
    for chi in chis:
        _, psi = DQ.trial_mps(mpo, chi)
        v = np.asarray(psi.to_dense()).reshape(-1)
        arrays[f"trial{chi}"] = v / np.linalg.norm(v)
    Hsp = H.to_sparse().tocsr()
    w, V = sla.eigsh(Hsp, k=2, which="SA", tol=1e-12, v0=arrays[f"trial{max(chis)}"].real)
    o = np.argsort(w)
    arrays["ground"] = V[:, o[0]]
    arrays["E01"] = w[o]
    np.savez(os.path.join(CACHE, _key(kind, size, J2) + ".npz"), **arrays)
    json.dump(dict(est=est, gamma=gam, E01=w[o].tolist(),
                   gamma_exact={c: float(abs(np.vdot(arrays["ground"], arrays[f"trial{c}"])) ** 2) for c in chis},
                   time=time.time() - t0),
              open(os.path.join(CACHE, _key(kind, size, J2) + ".json"), "w"))
    return _key(kind, size, J2), time.time() - t0


def run_inputs():
    P = _prev()
    todo = []
    for kind, size, J2, pk in SYSTEMS:
        if os.path.exists(os.path.join(CACHE, _key(kind, size, J2) + ".json")):
            continue
        chis = sorted({P[pk + (J2, e, s)]["chi"] for e in EPS for s in ("auto", "qetu")})
        todo.append((kind, size, J2, chis))
    with ProcessPoolExecutor(2) as ex:                  # memory: 2^20 sparse H per worker
        for k, dt in ex.map(_inputs, todo):
            print("[inputs]", k, f"{dt:.0f}s", flush=True)


def _case(args):
    kind, size, J2, eps, sched, chi, dt_prev = args
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost


    t0 = time.time()
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    info = json.load(open(os.path.join(CACHE, _key(kind, size, J2) + ".json")))
    z = np.load(os.path.join(CACHE, _key(kind, size, J2) + ".npz"))
    psi, g = z[f"trial{chi}"], z["ground"]
    est = EnergyEstimates(**info["est"])
    cost = mps_circuit_cost(N, chi)
    ans = Ansatz(gamma=info["gamma"][str(chi)], statevector=None, cost=cost, name=f"chi={chi}", meta=dict(chi=chi))
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    des, Hs, W = plan.design, plan.H_scaled, plan.window["W"]
    S = FS.Strang(Hs)
    anc = 2 if des.schedule == "qetu" else 1
    out = dict(kind=kind, size=size, N=N, J2=J2, eps=eps, schedule=sched, chi=chi, dt_prev=dt_prev,
               design=des.schedule, m=des.m, W=W, gamma_exact=info["gamma_exact"][str(chi)],
               trial_cnot=cost["cnot"], cert_steps=int(np.sum(plan.steps)),
               cert_cnot=plan.resources["expected_cnot"], cert_p=plan.resources["p_succ"], curve=[])
    if des.m == 0:
        out["time"] = time.time() - t0
        return out
    t = np.abs(des.times)
    for f in (1, 2, 4, 8):
        dtp = dt_prev / f
        k = np.maximum(1, np.ceil(t / (dtp * W) - 1e-12)).astype(int)
        if des.schedule == "qetu":
            dm = des.meta
            v = FS.qetu_trotter(S, psi, dm["degree"], dm["tau"], dm["x_s"], k)
        else:
            v = FS.filter_trotter(S, psi, des.times, des.phases, k)
        p = float(np.vdot(v, v).real)
        inf = max(0.0, 1.0 - abs(np.vdot(g, v / math.sqrt(p))) ** 2)
        cn = filter_counts(Hs, k, 2, anc_rot_per_pulse=anc)["cnot"]
        out["curve"].append(dict(dt=dtp, steps=int(k.sum()), filter_cnot=cn, inf=inf, p=p,
                                 cnot=(cn + cost["cnot"]) / p))
        print(kind, size, J2, eps, sched, f"dt={dtp:g} inf={inf:.2e}", flush=True)
        if f >= 2 and inf <= eps:
            break
    out["time"] = time.time() - t0
    return out


def cases():
    P = _prev()
    C = []
    for kind, size, J2, pk in SYSTEMS:
        for eps in EPS:
            for s in ("auto", "qetu"):
                r = P[pk + (J2, eps, s)]
                i = AP.coarsest_pass(r)
                C.append((kind, size, J2, eps, s, r["chi"], r["curve"][i]["dt"]))
    return C


def main():
    done = set()
    if os.path.exists(OUT):
        done = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]) for r in map(json.loads, open(OUT))}
    todo = [c for c in cases() if c[:5] not in done
            and os.path.exists(os.path.join(CACHE, _key(*c[:3]) + ".json"))]
    todo.sort(key=lambda c: (-c[3], c[0]))
    print(len(todo), "cases", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 4))) as ex, open(OUT, "a") as fh:
        for r in ex.map(_case, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()


if __name__ == "__main__":
    run_inputs() if sys.argv[1:] == ["inputs"] else main()
