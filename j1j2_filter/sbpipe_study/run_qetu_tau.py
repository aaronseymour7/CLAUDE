"""QETU query length sweep. The study runs QETU with tau = pi/2 in scaled units (the choice of Dong, Lin and
Tong); any tau <= pi/(1 + delta) keeps every excited state inside |cos(tau E)| <= x_s = cos(tau delta), so
the even Chebyshev filter still applies, and a longer query means fewer queries, each a longer (coarser)
Trotterised evolution. For every gap/overlap case (results/gapov.jsonl) and every measured crossover case
(results/emp_actual*.json[l]) we redo the measured QETU search at tau = 3 pi/4 and tau = pi/(1 + delta):
smallest even degree whose exact filter reaches eps, then the fewest Strang steps per query (or a slightly
higher degree) whose Trotterised circuit still does.  -> results/qetu_tau.jsonl

    python sbpipe_study/run_qetu_tau.py
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

import sbpipe
from sbpipe import j1j2_chain
from sbpipe import design as D
from sbpipe.resources import filter_counts
from sbpipe.simulate import qetu_trotter, state_metrics

import run_empirical as RE
import run_gapov as GO

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "qetu_tau.jsonl")
TAUS = ("3pi/4", "pi/(1+delta)")


def _tau(name, delta):
    return 0.75 * math.pi if name == "3pi/4" else math.pi / (1 + delta)


def _min_degree(Es, w, tau, x_s, eps, d_max=20000):
    d = 0
    while d < d_max:
        d += 2
        inf, _ = RE._infid(w, D.qetu_values(d, tau, x_s, Es))
        if inf <= eps:
            return d
    return None


def _trotter_search(Hs, psi, g, dmin, tau, x_s, eps, trial_cnot):
    best = None
    for dd in range(dmin, dmin + 40, 2):
        for k in (1, 2, 3, 4):
            if best is not None and dd * k >= best["d"] * best["k"]:
                break
            v = qetu_trotter(Hs, psi, dd, tau, x_s, [k] * dd)
            inf = state_metrics(v, g)["infidelity"]
            if inf <= eps:
                p = float(np.vdot(v, v).real)
                c = filter_counts(Hs, [k] * dd, 2, anc_rot_per_pulse=2)["cnot"]
                best = dict(d_exact=dmin, d=dd, k=k, inf=inf, p=p, cnot=(c + trial_cnot) / p)
                break
        if best is not None and best["k"] == 1:
            break
    return best


def _spectrum_weights(H, g, psi, E0):
    """Exact spectrum and overlaps of psi, per S^z block (N <= 14 dense; enough for the chains here)."""
    import scipy.linalg
    N = H.num_qubits
    Hsp = H.to_sparse().tocsr()
    idx = np.arange(2 ** N)
    pop = np.array([bin(i).count("1") for i in idx])
    E, W = [], []
    for n_up in range(N + 1):
        b = idx[pop == n_up]
        Hb = Hsp[b][:, b].toarray().real
        e, U = scipy.linalg.eigh(Hb, check_finite=False, driver="evr")
        E.append(e)
        W.append(np.abs(U.T @ psi[b]) ** 2)
    E, W = np.concatenate(E), np.concatenate(W)
    o = np.argsort(E)
    return E[o], W[o]


def _gapov_case(r):
    t0 = time.time()
    N, J2, gamma, eps = r["N"], r["J2"], r["gamma_target"], r["eps"]
    s, psi = GO._trial(N, J2, gamma)
    g = s["ground"]
    H = j1j2_chain(N, 1.0, J2)
    E0, W, delta = float(s["E0"]), r["W"], r["delta"]
    Hs = H.shifted_scaled(E0, W)
    E, w = _spectrum_weights(H, g, psi, E0)
    Es = (E - E0) / W
    out = dict(src="gapov", N=N, J2=J2, gamma_target=gamma, gap_target=r["gap_target"], study=r["study"],
               eps=eps, delta=delta, trial_cnot=r["trial_cnot"],
               tau_pi2=r.get("meas_qetu"))
    for name in TAUS:
        tau = _tau(name, delta)
        x_s = math.cos(tau * delta)
        d = _min_degree(Es, w, tau, x_s, eps)
        out[name] = _trotter_search(Hs, psi, g, d, tau, x_s, eps, r["trial_cnot"]) if d else None
    out["time"] = time.time() - t0
    return out


def _act_case(r):
    import run_extended as RX
    t0 = time.time()
    kind, size, J2, chi, eps = r["kind"], r["size"], r["J2"], r["chi"], r["eps"]
    z = np.load(os.path.join(RE.CACHEE, RE._act_key(kind, size, J2) + ".npz"))
    E, w, psi, g = z["E"], z[f"w{chi}"], z[f"trial{chi}"], z["ground"]
    H = RX.model(kind, size, J2)
    lo, hi = sbpipe.weyl_bounds(H)
    Wr = hi - E[0]
    Es = (E - E[0]) / Wr
    delta = r["delta"]
    Hs = H.shifted_scaled(E[0], Wr)
    out = dict(src="actual", kind=kind, size=size, N=r["N"], J2=J2, chi=chi, eps=eps, delta=delta,
               trial_cnot=r["trial_cnot"], tau_pi2=r.get("qetu"))
    for name in TAUS:
        tau = _tau(name, delta)
        x_s = math.cos(tau * delta)
        d = _min_degree(Es, w, tau, x_s, eps)
        out[name] = _trotter_search(Hs, psi, g, d, tau, x_s, eps, r["trial_cnot"]) if d else None
    out["time"] = time.time() - t0
    return out


def _run(arg):
    kind, r = arg
    return _gapov_case(r) if kind == "gapov" else _act_case(r)


def _key(o):
    if o["src"] == "gapov":
        return ("gapov", o["N"], round(o["J2"], 9), o["gamma_target"], o["eps"])
    return ("actual", o["kind"], o["size"], round(o["J2"], 9), o["chi"], o["eps"])


def todo():
    G = [json.loads(l) for l in open(os.path.join(RES, "gapov.jsonl"))]
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))]
    A += json.load(open(os.path.join(RES, "emp_actual.json")))
    jobs = [("gapov", dict(r, src="gapov")) for r in G if r.get("meas_qetu")]
    jobs += [("actual", dict(r, src="actual")) for r in A if r.get("qetu")]
    done = set()
    if os.path.exists(OUT):
        done = {_key(json.loads(l)) for l in open(OUT)}
    jobs = [j for j in jobs if _key(j[1]) not in done]
    jobs.sort(key=lambda j: j[1]["N"])
    return jobs


def main():
    jobs = todo()
    print(len(jobs), "cases", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 4))) as ex, open(OUT, "a") as fh:
        for o in ex.map(_run, jobs):
            fh.write(json.dumps(o) + "\n")
            fh.flush()
            b = {n: (o[n]["cnot"] if o[n] else None) for n in TAUS}
            a = o["tau_pi2"]["cnot"] if o.get("tau_pi2") else None
            print(_key(o), f"pi/2 {a}", b, f"{o['time']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
