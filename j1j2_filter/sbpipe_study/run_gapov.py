"""Gap versus trial overlap, measured and certified, on the J1-J2 chain (N = 8-16, exact state vectors).

The gap is held fixed over N by tuning J2 (bisection on ED), and the overlap is held fixed over N by
tuning a real trial circuit: the variational MPS (chi = 16) followed by one layer of single-qubit
rotations exp(-i theta Y / 2) on every other qubit, theta chosen by bisection so that |<E0|psi>|^2 = gamma.
For each (N, gap, gamma, eps):

  certified : the pipeline's plan (dyadic and QETU, second order) with exact inputs
  measured  : the shortest filter that actually reaches eps (exact filter on the state vector), then
              the coarsest Trotter step at which the full Trotterised circuit still reaches eps

    python sbpipe_study/run_gapov.py tune      # J2(N, gap) and theta(N, gap, gamma)
    python sbpipe_study/run_gapov.py run       # all cases
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.sparse.linalg as sla
from scipy.sparse.linalg import eigsh

import sbpipe
from sbpipe import j1j2_chain
from sbpipe import design as D
from sbpipe.inputs import Ansatz, EnergyEstimates
from sbpipe.pipeline import Config, refine
from sbpipe.resources import filter_counts
from sbpipe.simulate import filter_trotter, qetu_trotter, state_metrics

import run_study as RS
import run_empirical as RE

dump, load, JOBS = RS.dump, RS.load, RS.JOBS
CACHEG = os.path.join(RS.RES, "cache_gapov")
os.makedirs(CACHEG, exist_ok=True)

N_GO = [8, 10, 12, 14, 16]
GAPS = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
GAMMAS = [0.03, 0.1, 0.3, 0.5, 0.7, 0.9]
EPSS = [1e-2, 1e-3, 1e-4]
GAP_A = (0.40, 0.50)            # study A: gap fixed, gamma varied
GAMMA_B = (0.1, 0.5)            # study B: gamma fixed, gap varied
CHI_BASE = 16


# ----------------------------------------------------------------------------- tuning
def _gap_J2(args):
    N, target = args
    lo, hi = 0.0, 0.75
    g = lambda J2: RE._gap(N, J2)
    if g(lo) >= target:                             # the gap is flat below J2 ~ 0.3 and rises above
        return N, target, (0.0 if abs(g(lo) - target) < 2e-3 else None)
    if g(hi) <= target:
        return N, target, None
    lo = 0.25 if g(0.25) < target else 0.0
    for _ in range(32):
        mid = 0.5 * (lo + hi)
        if g(mid) < target:
            lo = mid
        else:
            hi = mid
    return N, target, round(0.5 * (lo + hi), 7)


def _rotate(psi, N, theta):
    """exp(-i theta Y/2) on the even qubits of a state vector (real rotation)."""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    v = psi.reshape([2] * N)
    for q in range(0, N, 2):                       # one sublattice: a global rotation would leave the
        v = np.moveaxis(v, q, 0)                    # singlet ground state invariant
        a, b = v[0].copy(), v[1].copy()
        v = np.stack([c * a - s * b, s * a + c * b])
        v = np.moveaxis(v, 0, q)
    return v.reshape(-1)


def _system(N, J2):
    """Ground state, gap, chi = 16 MPS trial vector; cached."""
    p = os.path.join(CACHEG, f"N{N}_J2{J2:.7g}.npz")
    if os.path.exists(p):
        z = np.load(p)
        return {k: z[k] for k in z.files}
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    H = j1j2_chain(N, 1.0, J2)
    w, V = eigsh(H.to_sparse(), k=6, which="SA", tol=1e-12)
    o = np.argsort(w)
    w, V = w[o], V[:, o]
    i1 = int(np.argmax(w > w[0] + 1e-8))
    _, mps = trial_mps(pauli_to_mpo(H), CHI_BASE)
    t = np.asarray(mps.to_dense()).reshape(-1)
    t = t / np.linalg.norm(t)
    out = dict(E0=w[0], E1=w[i1], ground=V[:, 0].astype(complex), trial=t.astype(complex))
    np.savez(p, **out)
    return out


def _theta_for(args):
    N, J2, gamma = args
    s = _system(N, J2)
    g, t = s["ground"], s["trial"]
    f = lambda th: abs(np.vdot(g, _rotate(t, N, th))) ** 2
    if f(0.0) < gamma:
        return N, J2, gamma, None
    lo, hi = 0.0, math.pi / 2
    if f(hi) > gamma:
        hi = math.pi
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if f(mid) > gamma:
            lo = mid
        else:
            hi = mid
    th = 0.5 * (lo + hi)
    return N, J2, gamma, dict(theta=th, gamma=f(th), gamma_mps=f(0.0))


def run_tune():
    with ProcessPoolExecutor(JOBS) as ex:
        J = list(ex.map(_gap_J2, [(N, g) for N in N_GO for g in GAPS]))
    pts = [(N, g, J2) for N, g, J2 in J if J2 is not None]
    for N, g, J2 in J:
        print(f"[tune] N={N} gap={g} J2={J2}", flush=True)
    with ProcessPoolExecutor(JOBS) as ex:
        list(ex.map(_system, [N for N, g, J2 in pts], [J2 for N, g, J2 in pts]))
        T = list(ex.map(_theta_for, [(N, J2, gm) for N, g, J2 in pts for gm in GAMMAS]))
    out = dict(points=[dict(N=N, gap=g, J2=J2) for N, g, J2 in pts],
               theta=[dict(N=N, J2=J2, gamma=gm, theta=r["theta"], gamma_true=r["gamma"], gamma_mps=r["gamma_mps"])
                      for N, J2, gm, r in T if r])
    for N, J2, gm, r in T:
        print(f"[tune] N={N} J2={J2:.4f} gamma={gm} -> {r}", flush=True)
    dump("gapov_tune.json", out)


# ----------------------------------------------------------------------------- one case
def _trial(N, J2, gamma):
    s = _system(N, J2)
    th = [r for r in load("gapov_tune.json")["theta"] if r["N"] == N and abs(r["J2"] - J2) < 1e-9
          and r["gamma"] == gamma][0]["theta"]
    return s, _rotate(s["trial"], N, th)


def _qetu_measured(Hsp, Hs, psi, g, delta, eps, trial_cnot):
    """Smallest even degree whose exact Chebyshev filter reaches eps (one recurrence pass), then the
    fewest Strang steps per query (or a slightly higher degree at one step) that still reach eps."""
    tau = math.pi / 2
    x_s = math.cos(tau * delta)
    X = lambda v: 0.5 * (sla.expm_multiply(-1j * tau * Hsp, v) + sla.expm_multiply(1j * tau * Hsp, v))
    v0 = psi.astype(complex)
    v1 = X(v0) / x_s
    d, dmin = 1, None
    while d < 5000:
        v0, v1 = v1, 2.0 * X(v1) / x_s - v0
        d += 1
        if d % 2 == 0:
            n = np.vdot(v1, v1).real
            if 1 - abs(np.vdot(g, v1)) ** 2 / n <= eps:
                dmin = d
                break
    if dmin is None:
        return None
    best = None
    for dd in range(dmin, dmin + 40, 2):
        for k in (1, 2, 3):
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


def _dyadic_measured(Hsp, Hs, psi, g, delta, eps, trial_cnot, n_cand=2):
    """Shortest dyadic ladders (K, r) whose exact filter reaches eps on this state, then the coarsest
    common Strang step at which the Trotterised circuit still does."""
    tau = math.pi / (1 + delta)
    cands = []
    for K in range(1, 16):
        v = psi.astype(complex)
        for r in range(1, 200):
            for k in range(K):
                t = tau * 2 ** k
                v = 0.5 * (sla.expm_multiply(-1j * t * Hsp, v) + sla.expm_multiply(1j * t * Hsp, v))
            n = np.vdot(v, v).real
            if 1 - abs(np.vdot(g, v)) ** 2 / n <= eps:
                cands.append((r * tau * (2 ** K - 1), K, r))
                break
            if cands and r * tau * (2 ** K - 1) > 2 * min(c[0] for c in cands):
                break
    cands.sort()
    best = None
    for T, K, r in cands[:n_cand]:
        times = np.tile(tau * 2.0 ** np.arange(K), r)
        des = D.FilterDesign("dyadic", times, np.zeros(len(times)), delta, 0.0, 0.0, 1.0, {})
        tr = RE._sbc_trotter(None, Hs, psi, g, des, eps)
        if tr:
            dt, k, inf, p = tr
            c = (filter_counts(Hs, k, 2)["cnot"] + trial_cnot) / p
            if best is None or c < best["cnot"]:
                best = dict(T=T, K=K, r=r, m=len(times), steps=int(k.sum()), dt=dt, inf=inf, p=p, cnot=c)
    return best


def _case(args):
    study, N, gap, J2, gamma, eps = args
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    t0 = time.time()
    s, psi = _trial(N, J2, gamma)
    g = s["ground"]
    gam = abs(np.vdot(g, psi)) ** 2
    H = j1j2_chain(N, 1.0, J2)
    lo, hi = sbpipe.weyl_bounds(H)
    E0, E1 = float(s["E0"]), float(s["E1"])
    W = hi - E0
    delta = (E1 - E0) / W
    Hs = H.shifted_scaled(E0, W)
    Hsp = Hs.to_sparse().tocsr()
    mc = mps_circuit_cost(N, CHI_BASE)
    trial_cost = dict(cnot=mc["cnot"], rotations=mc["rotations"] + N // 2)
    out = dict(study=study, N=N, gap=E1 - E0, gap_target=gap, J2=J2, gamma=gam, gamma_target=gamma, eps=eps,
               W=W, delta=delta, trial_cnot=trial_cost["cnot"])
    # certified, exact inputs (tiny error bars), dyadic and QETU
    est = EnergyEstimates(E0=E0, E0_err=1e-9, E1=E1, E1_err=1e-9, E_max=None, source="ED")
    ans = Ansatz(gamma=gam, statevector=None, cost=trial_cost, name="mps16+Y", meta={})
    for sched in ("dyadic", "qetu"):
        try:
            pl = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
            sm = pl.summary()
            out[f"cert_{sched}"] = dict(cnot=sm["res_expected_cnot"], m=sm["m"], T=sm["T_scaled"],
                                        steps=sm["res_steps"])
        except RuntimeError:
            out[f"cert_{sched}"] = None
    out["meas_qetu"] = _qetu_measured(Hsp, Hs, psi, g, delta, eps, trial_cost["cnot"])
    out["meas_dyadic"] = _dyadic_measured(Hsp, Hs, psi, g, delta, eps, trial_cost["cnot"])
    out["time"] = time.time() - t0
    return out


def cases():
    T = load("gapov_tune.json")
    pts = {(p["N"], p["gap"]): p["J2"] for p in T["points"]}
    ok = {(r["N"], round(r["J2"], 7), r["gamma"]) for r in T["theta"]}
    c = []
    for (N, gap), J2 in pts.items():
        for eps in EPSS:
            if gap in GAP_A:
                c += [("A", N, gap, J2, gm, eps) for gm in GAMMAS if (N, round(J2, 7), gm) in ok]
            c += [("B", N, gap, J2, gm, eps) for gm in GAMMA_B
                  if (N, round(J2, 7), gm) in ok and not (gap in GAP_A)]
    nmax = int(os.environ.get("GO_NMAX", 14))         # N=16 runs ~1.5 h per case; off by default
    c = [a for a in c if a[1] <= nmax]
    c.sort(key=lambda a: (a[1], a[5] < 1e-3))          # small N first so the scaling fills in early
    return c


def run_run():
    import json
    path = os.path.join(RS.RES, "gapov.jsonl")
    out = [json.loads(l) for l in open(path)] if os.path.exists(path) else []
    done = {(r["study"], r["N"], r["gap_target"], r["gamma_target"], r["eps"]) for r in out}
    todo = [a for a in cases() if (a[0], a[1], a[2], a[4], a[5]) not in done]
    with ProcessPoolExecutor(JOBS) as ex, open(path, "a") as fh:
        for r in ex.map(_case, todo):
            fh.write(json.dumps(r, default=RS._f) + "\n")
            fh.flush()
            q, d = r["meas_qetu"], r["meas_dyadic"]
            cq, cd = r["cert_qetu"], r["cert_dyadic"]
            print(f"[gapov] {r['study']} N={r['N']} gap={r['gap']:.3f} gamma={r['gamma']:.3f} eps={r['eps']:g} "
                  f"meas qetu={q['cnot'] if q else float('nan'):.2e} dyadic={d['cnot'] if d else float('nan'):.2e} "
                  f"cert qetu={cq['cnot'] if cq else float('nan'):.2e} dyadic={cd['cnot'] if cd else float('nan'):.2e} "
                  f"({r['time']:.0f}s)", flush=True)
            out.append(r)
    dump("gapov.json", out)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "tune,run"
    for k in what.split(","):
        t = time.time()
        dict(tune=run_tune, run=run_run)[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)
