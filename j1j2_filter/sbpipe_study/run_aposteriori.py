"""A-posteriori Trotter step for the planned (Theorem 1) filters.

The certified plan sizes the filter from the estimated gap and overlap (no oracle) but takes its Trotter
step from a worst-case commutator bound, 60-75 times finer than what the circuit needs. Here the filter is
the planned one and only the step is chosen without the bound. For every chain and cylinder target of the
crossover we take the trial chi whose certified plan is cheapest ('auto' = SBC, and 'qetu'), and run the
Trotterised planned circuit on a common grid of physical steps dt = 4, 2, 1, ..., 2^-9 (units 1/J1;
k_i = ceil(t_i / dt) Strang steps per pulse or per QETU query), recording at each step

    infidelity to the exact ground state (the oracle check),
    success probability p and post-selected energy E (what a device measures),
    the variances of the XX, YY and ZZ parts of H in the post-selected state (shot cost of measuring E).

Two step rules are then evaluated offline from these curves (make_figures.py):
    energy convergence : halve dt until |E(dt) - E(2 dt)| <= eps * gap / 2, use dt
    size transfer      : use half the coarsest passing dt of the next-smaller system (same J2, eps)
-> results/aposteriori.jsonl

    python sbpipe_study/run_aposteriori.py
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

import numpy as np

from sbpipe.pauli import PauliSum
from sbpipe.pipeline import Config, refine
from sbpipe.resources import filter_counts
from sbpipe.simulate import filter_trotter, qetu_trotter

import run_empirical as RE
import run_extended as RX

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "aposteriori.jsonl")
DTS = [4.0 / 2 ** n for n in range(12)]


def cases():
    C = [json.loads(l) for l in open(os.path.join(RES, "cert_crossover.jsonl"))]
    best = {}
    for q in C:
        if not q.get("plan"):
            continue
        k = (q["kind"], q["size"], q["J2"], q["eps"], q["schedule"])
        c = q["plan"]["res_expected_cnot"]
        if k not in best or c < best[k][1]:
            best[k] = (q["chi"], c)
    return sorted((k[0], k[1], k[2], k[3], v[0], k[4]) for k, v in best.items())


def _groups(H):
    """Sparse XX-, YY- and ZZ-parts of H (the three bases a device measures the Heisenberg energy in)."""
    h = H.without_identity()
    out = []
    for P in "XYZ":
        idx = [i for i, l in enumerate(h.labels) if set(l) <= {"I", P}]
        out.append(PauliSum([h.labels[i] for i in idx], h.coeffs[idx]).to_sparse().tocsr())
    assert sum(len([1 for l in h.labels if set(l) <= {"I", P}]) for P in "XYZ") == len(h)
    return out


def _case(args, only=None):
    kind, size, J2, eps, chi, sched = args
    t0 = time.time()
    rec, est, ans, _ = RX.get_case(kind, size, J2, chi)
    ans.statevector = None
    H = RX.model(kind, size, J2)
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    des, Hs, W = plan.design, plan.H_scaled, plan.window["W"]
    z = np.load(os.path.join(RE.CACHEE, RE._act_key(kind, size, J2) + ".npz"))
    E, psi, g = z["E"], z[f"trial{chi}"], z["ground"]
    gap = float(E[np.argmax(E > E[0] + 1e-8)] - E[0])
    Hsp = H.to_sparse().tocsr()
    G = _groups(H)
    trial = float(z[f"cost{chi}"])
    anc = 2 if des.schedule == "qetu" else 1
    out = dict(kind=kind, size=size, N=H.num_qubits, J2=J2, eps=eps, chi=chi, schedule=sched,
               E0=float(E[0]), gap=gap, W=W, design=des.schedule, m=des.m, T_phys=des.T / W,
               trial_cnot=trial, cert_steps=int(np.sum(plan.steps)),
               cert_cnot=plan.resources["expected_cnot"], cert_p=plan.resources["p_succ"],
               curve=[])
    t = np.abs(des.times)
    passed = 0
    for dtp in (DTS if only is None else only):
        dt = dtp * W                                    # scaled units
        k = np.maximum(1, np.ceil(t / dt - 1e-12)).astype(int)
        if des.schedule == "qetu":
            dm = des.meta
            v = qetu_trotter(Hs, psi, dm["degree"], dm["tau"], dm["x_s"], k)
        else:
            v = filter_trotter(Hs, psi, des.times, des.phases, k, 2)
        p = float(np.vdot(v, v).real)
        u = v / math.sqrt(p)
        inf = max(0.0, 1.0 - abs(np.vdot(g, u)) ** 2)
        en = float(np.vdot(u, Hsp @ u).real)
        var = []
        for O in G:
            Ou = O @ u
            var.append(max(0.0, float(np.vdot(Ou, Ou).real - np.vdot(u, Ou).real ** 2)))
        cn = filter_counts(Hs, k, 2, anc_rot_per_pulse=anc)["cnot"]
        pt = dict(dt=dtp, steps=int(k.sum()), filter_cnot=cn, inf=inf, p=p, E=en, var=var)
        c = out["curve"]
        if c and c[-1]["steps"] == pt["steps"]:      # same circuit as the coarser dt: keep the finer label
            c[-1] = pt
            continue
        c.append(pt)
        passed = passed + 1 if inf <= eps else 0
        if only is not None:
            continue
        conv = len(c) >= 2 and abs(c[-1]["E"] - c[-2]["E"]) <= eps * gap / 2
        if passed >= 2 and conv:
            break
        if int(k.sum()) > 40 * out["cert_steps"]:
            break
    out["time"] = time.time() - t0
    return out


PREV = {("chain", 14): ("chain", 12), ("chain", 16): ("chain", 14), ("cyl", 4): ("cyl", 3)}


def dedupe(curve):
    c = []
    for q in curve:
        if c and q["steps"] == c[-1]["steps"]:
            c[-1] = q
        else:
            c.append(q)
    return c


def coarsest_pass(r):
    """Index of the coarsest point from which every finer point passes, or None."""
    c, best = r["curve"], None
    for i in range(len(c) - 1, -1, -1):
        if c[i]["inf"] > r["eps"]:
            break
        best = i
    return best


def transfer_dt(R, r):
    """Size transfer: half the coarsest passing step of the next-smaller system (same J2, eps, schedule)."""
    pk = PREV.get((r["kind"], r["size"]))
    if pk is None:
        return None
    s = R.get(pk + (r["J2"], r["eps"], r["schedule"]))
    if s is None or coarsest_pass(s) is None:
        return None
    return s["curve"][coarsest_pass(s)]["dt"] / 2


def _fill(args):
    key, dt = args
    return key, _case(key, only=[dt])["curve"]


def fill():
    """Simulate the transferred step wherever the stored curve stopped before reaching it."""
    rows = [json.loads(l) for l in open(OUT)]
    for r in rows:
        r["curve"] = dedupe(r["curve"])
    R = {(r["kind"], r["size"], r["J2"], r["eps"], r["schedule"]): r for r in rows}
    todo = []
    for r in rows:
        dt = transfer_dt(R, r)
        if dt is not None and dt < min(q["dt"] for q in r["curve"]) - 1e-12:
            todo.append(((r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"]), dt))
    print(len(todo), "transfer points", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 4))) as ex:
        for key, c in ex.map(_fill, todo):
            r = R[key[:4] + (key[5],)]
            r["curve"] = dedupe(r["curve"] + c)
    with open(OUT, "w") as fh:
        fh.write("".join(json.dumps(r) + "\n" for r in rows))


def main():
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l)
            done.add((r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"]))
    todo = [c for c in cases() if c not in done]
    todo.sort(key=lambda c: (c[1] if c[0] == "chain" else 4 * c[1]))
    print(len(todo), "cases", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 4))) as ex, open(OUT, "a") as fh:
        for r in ex.map(_case, todo):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            c = r["curve"][-1]
            print(r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"], len(r["curve"]),
                  f"dt={c['dt']:.3g} inf={c['inf']:.2e}", f"{r['time']:.0f}s", flush=True)


if __name__ == "__main__":
    fill() if sys.argv[1:] == ["fill"] else main()
