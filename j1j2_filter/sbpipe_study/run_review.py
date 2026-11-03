"""Runs added after the internal review of the measured draft (paper/paper_measured.tex).

    python sbpipe_study/run_review.py trunc    # second overlap mechanism: bond-truncated ground states
    python sbpipe_study/run_review.py mono     # is "passes" monotone in the Trotter step?
    python sbpipe_study/run_review.py bign     # dyadic / optimized / QETU at N = 16, 18, 20 (gap 0.40)

All measured searches use the same rules as run_gapov.py (shortest exact filter, then the coarsest common
Strang step that still reaches eps), with the fast bond-SWAP Strang formula of fastsim.py.
-> results/review_trunc.jsonl, results/review_mono.jsonl, results/review_bign.jsonl
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import itertools
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import scipy.sparse.linalg as sla
from scipy.sparse.linalg import eigsh

import sbpipe
from sbpipe import j1j2_chain
from sbpipe import design as D
from sbpipe.resources import filter_counts

from fastsim import Strang, filter_trotter as ftrot, qetu_trotter as qtrot

RES = os.path.join(HERE, "results")
CACHE = os.path.join(RES, "cache_review")
os.makedirs(CACHE, exist_ok=True)
JOBS = int(os.environ.get("JOBS", "4"))
EPSS = [1e-2, 1e-3, 1e-4]
GAMMAS = [0.03, 0.1, 0.3, 0.5, 0.7, 0.9]


def _f(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


# ----------------------------------------------------------------------------- systems
def gap_J2(N, target):
    """J2 at which the singlet-triplet gap of the open chain equals target (bisection on Lanczos)."""
    tune = os.path.join(RES, "gapov_tune.json")
    if os.path.exists(tune):
        for p in json.load(open(tune))["points"]:
            if p["N"] == N and abs(p["gap"] - target) < 1e-9:
                return p["J2"]
    p = os.path.join(CACHE, f"J2_N{N}_g{target}.json")
    if os.path.exists(p):
        return json.load(open(p))["J2"]

    def gap(J2):
        w = eigsh(j1j2_chain(N, 1.0, J2).to_sparse(), k=4, which="SA", tol=1e-10, return_eigenvectors=False)
        w = np.sort(w)
        return float(w[np.argmax(w > w[0] + 1e-8)] - w[0])
    lo, hi = 0.3, 0.75
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if gap(mid) < target:
            lo = mid
        else:
            hi = mid
    J2 = round(0.5 * (lo + hi), 7)
    json.dump(dict(N=N, gap=target, J2=J2), open(p, "w"))
    return J2


def system(N, J2, chi=16):
    """Ground state, E0, E1 and the variational chi = 16 MPS trial vector (cached)."""
    p = os.path.join(CACHE, f"sys_N{N}_J2{J2:.7g}.npz")
    if os.path.exists(p):
        z = np.load(p)
        return {k: z[k] for k in z.files}
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    H = j1j2_chain(N, 1.0, J2)
    w, V = eigsh(H.to_sparse(), k=6, which="SA", tol=1e-12)
    o = np.argsort(w)
    w, V = w[o], V[:, o]
    i1 = int(np.argmax(w > w[0] + 1e-8))
    _, mps = trial_mps(pauli_to_mpo(H), chi)
    t = np.asarray(mps.to_dense()).reshape(-1)
    out = dict(E0=w[0], E1=w[i1], ground=V[:, 0].astype(complex), trial=(t / np.linalg.norm(t)).astype(complex))
    np.savez(p, **out)
    return out


def rotate(psi, N, theta):
    """exp(-i theta Y/2) on the even qubits (the first overlap mechanism, as in run_gapov.py)."""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    v = psi.reshape([2] * N)
    for q in range(0, N, 2):
        v = np.moveaxis(v, q, 0)
        a, b = v[0].copy(), v[1].copy()
        v = np.stack([c * a - s * b, s * a + c * b])
        v = np.moveaxis(v, 0, q)
    return v.reshape(-1)


def theta_for(g, t, N, gamma):
    f = lambda th: abs(np.vdot(g, rotate(t, N, th))) ** 2
    lo, hi = 0.0, math.pi / 2
    if f(hi) > gamma:
        hi = math.pi
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > gamma else (lo, mid)
    return 0.5 * (lo + hi)


def truncate(psi, N, cuts, chi_max=16):
    """Left-to-right SVD of psi with Schmidt rank 1 on the bonds in `cuts` (bond i joins sites i, i+1)
    and at most chi_max elsewhere: an MPS of bond dimension <= 16 whose error is an entanglement
    truncation, the second overlap mechanism."""
    v = psi.reshape(1, -1)
    As = []
    for i in range(N - 1):
        Dl = v.shape[0]
        U, S, Vh = np.linalg.svd(v.reshape(Dl * 2, -1), full_matrices=False)
        k = 1 if i in cuts else min(len(S), chi_max)
        As.append(U[:, :k])
        v = S[:k, None] * Vh[:k]
    out = v
    for A in reversed(As):
        out = (A @ out.reshape(A.shape[1], -1)).reshape(-1)
    out = out.reshape(-1)
    return out / np.linalg.norm(out)


# ----------------------------------------------------------------------------- measured searches
import run_cluster as _RC


class Case(_RC.Case):
    """run_cluster.Case (matrix-free, about 10x faster than sparse expm) plus the helpers used here."""

    def __init__(self, H, E0, E1, g, psi, trial_cnot):
        super().__init__(H, E0, E1, g, psi, trial_cnot)
        self.H = H

    def qetu_scan(self, d, tau, kmax=6):
        x_s = math.cos(tau * self.delta)
        out = []
        for k in range(1, kmax + 1):
            v0 = self.psi.copy()
            v1 = self.Xt(v0, tau, k) / x_s
            for _ in range(2, d + 1):
                v0, v1 = v1, 2.0 * self.Xt(v1, tau, k) / x_s - v0
            out.append((k, self.inf(v1)[0]))
        return out

    def all(self, eps, scheds=("dyadic", "optimized", "qetu")):
        out = dict(E0=self.E0, gap=self.E1 - self.E0, W=self.W, delta=self.delta, gamma=self.gamma,
                   trial_cnot=self.trial_cnot, eps=eps)
        for f in scheds:
            if f == "qetu":
                for name in ("qetu_pi2", "qetu_3pi4", "qetu_max"):
                    out[name] = self.run(name, eps)
            else:
                out[f] = self.run(f, eps)
        return out


def trial_cost(N, chi=16, extra_rot=0):
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    return mps_circuit_cost(N, chi)["cnot"]


def _append(path, rows):
    with open(path, "a") as fh:
        for r in rows:
            fh.write(json.dumps(r, default=_f) + "\n")


def _done(path, keyf):
    if not os.path.exists(path):
        return set()
    return {keyf(json.loads(l)) for l in open(path)}


# ----------------------------------------------------------------------------- 1. truncation mechanism
TR_N = [8, 10, 12, 14]
TR_GAPS = [0.40, 0.50]


def residual_energy(g_E, w_E, E0, gap):
    """Mean excitation energy of the residual, in units of the gap."""
    m = w_E[1:].sum()
    return float(np.sum(w_E[1:] * (g_E[1:] - E0)) / m / gap) if m > 0 else 0.0


def trunc_trials(N, J2):
    """Bond subsets whose rank-1 truncation of the ground state lands closest to each target overlap."""
    p = os.path.join(CACHE, f"trunc_N{N}_J2{J2:.7g}.json")
    if os.path.exists(p):
        return json.load(open(p))
    s = system(N, J2)
    g = s["ground"]
    found = []
    bonds = range(N - 1)
    for size in range(1, 6):
        for cuts in itertools.combinations(bonds, size):
            v = truncate(g, N, set(cuts))
            found.append((abs(np.vdot(g, v)) ** 2, list(cuts)))
    out = []
    for gm in GAMMAS:
        gam, cuts = min(found, key=lambda x: abs(math.log(x[0] / gm)))
        if abs(math.log(gam / gm)) < math.log(1.5) and all(o["cuts"] != cuts for o in out):
            out.append(dict(gamma_target=gm, gamma=gam, cuts=cuts))
    json.dump(out, open(p, "w"))
    return out


def _trunc_case(args):
    N, gap, gm, eps = args
    t0 = time.time()
    J2 = gap_J2(N, gap)
    s = system(N, J2)
    sel = [r for r in trunc_trials(N, J2) if r["gamma_target"] == gm][0]
    psi = truncate(s["ground"], N, set(sel["cuts"])).astype(complex)
    H = j1j2_chain(N, 1.0, J2)
    c = Case(H, float(s["E0"]), float(s["E1"]), s["ground"], psi, trial_cost(N))
    out = dict(mech="trunc", N=N, gap_target=gap, J2=J2, gamma_target=gm, cuts=sel["cuts"])
    out.update(c.all(eps, ("dyadic", "optimized", "qetu")))
    out["resid_E"] = _resid(c, psi)
    out["time"] = time.time() - t0
    return out


def _resid(c, psi, k=40):
    """Mean excitation energy of the residual (psi minus its ground component) in units of the gap,
    <r|H|r>/<r|r> - E0, from the sparse Hamiltonian."""
    r = psi - np.vdot(c.g, psi) * c.g
    n = np.vdot(r, r).real
    if n < 1e-14:
        return 0.0
    Hsp = c.H.to_sparse()
    return float((np.vdot(r, Hsp @ r).real / n - c.E0) / (c.E1 - c.E0))


def _rot_resid(args):
    """Residual energy of the rotation trials, for comparison with the truncation trials."""
    N, gap, gm = args
    J2 = gap_J2(N, gap)
    s = system(N, J2)
    th = theta_for(s["ground"], s["trial"], N, gm)
    psi = rotate(s["trial"], N, th)
    H = j1j2_chain(N, 1.0, J2)
    c = Case(H, float(s["E0"]), float(s["E1"]), s["ground"], psi, 0)
    return dict(mech="rot", N=N, gap_target=gap, gamma_target=gm, gamma=c.gamma, resid_E=_resid(c, psi))


def run_trunc():
    path = os.path.join(RES, "review_trunc.jsonl")
    for N in TR_N:                                   # tune serially (cached), then fan out
        for gap in TR_GAPS:
            trunc_trials(N, gap_J2(N, gap))
    done = _done(path, lambda r: (r["N"], r["gap_target"], r["gamma_target"], r["eps"]))
    todo = []
    for N in TR_N:
        for gap in TR_GAPS:
            got = {r["gamma_target"] for r in trunc_trials(N, gap_J2(N, gap))}
            todo += [(N, gap, gm, e) for gm in GAMMAS if gm in got for e in EPSS if (N, gap, gm, e) not in done]
    todo.sort(key=lambda a: (a[0], -a[3]))
    print(len(todo), "truncation cases", flush=True)
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_trunc_case, todo):
            _append(path, [r])
            d = r.get("dyadic") or {}
            print(f"[trunc] N={r['N']} gap={r['gap_target']} gamma={r['gamma']:.3f} eps={r['eps']:g} "
                  f"dy={d.get('cnot', float('nan')):.3e} rE={r['resid_E']:.2f} ({r['time']:.0f}s)", flush=True)
    rp = os.path.join(RES, "review_rot_resid.json")
    with ProcessPoolExecutor(JOBS) as ex:
        R = list(ex.map(_rot_resid, [(N, gap, gm) for N in TR_N for gap in TR_GAPS for gm in GAMMAS]))
    json.dump(R, open(rp, "w"), default=_f)


# ----------------------------------------------------------------------------- 2. monotonicity
MONO_N = [8, 10, 12]


def _mono_case(args):
    """Dense scan of the common Strang step for the measured dyadic and optimized filters, and of
    k = 1..6 steps per query for QETU (tau = pi/2): does infidelity <= eps hold on a single interval
    dt <= dt*?"""
    N, gap, gm, eps = args
    t0 = time.time()
    J2 = gap_J2(N, gap)
    s = system(N, J2)
    psi = rotate(s["trial"], N, theta_for(s["ground"], s["trial"], N, gm)).astype(complex)
    H = j1j2_chain(N, 1.0, J2)
    c = Case(H, float(s["E0"]), float(s["E1"]), s["ground"], psi, trial_cost(N))
    out = dict(N=N, gap_target=gap, gamma_target=gm, eps=eps, W=c.W)
    for name, f in (("dyadic", c.dyadic), ("optimized", c.optimized)):
        b = f(eps)
        if b is None:
            out[name] = None
            continue
        times, phases = np.asarray(b["times"]), np.asarray(b["phases"])
        tmax = float(np.abs(times).max())
        grid = np.geomspace(tmax, b["dt"] / 8, 60)
        scan = []
        seen = set()
        for dt in grid:
            k = tuple(np.maximum(1, np.ceil(np.abs(times) / dt - 1e-12)).astype(int))
            if k in seen:
                continue
            seen.add(k)
            _, inf, p = c.sbc_run(times, phases, dt)
            scan.append((float(dt), int(sum(k)), inf))
        out[name] = dict(dt_bisect=b["dt"], scan=scan)
    q = c.qetu(eps, math.pi / 2)
    if q:
        out["qetu"] = dict(d=q["d"], scan=c.qetu_scan(q["d"], q["tau"]))
    out["time"] = time.time() - t0
    return out


def run_mono():
    path = os.path.join(RES, "review_mono.jsonl")
    todo = []
    for N in MONO_N:
        for e in EPSS:
            todo += [(N, 0.40, gm, e) for gm in (0.03, 0.3, 0.9)]
            todo += [(N, gp, 0.1, e) for gp in (0.30, 0.45, 0.55)]
    done = _done(path, lambda r: (r["N"], r["gap_target"], r["gamma_target"], r["eps"]))
    todo = [a for a in todo if a not in done]
    print(len(todo), "monotonicity cases", flush=True)
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_mono_case, todo):
            _append(path, [r])
            print(f"[mono] N={r['N']} gap={r['gap_target']} gamma={r['gamma_target']} eps={r['eps']:g} "
                  f"({r['time']:.0f}s)", flush=True)


# ----------------------------------------------------------------------------- 3. larger N
BIG_N = [16, 18, 20]


def _big_case(args):
    N, gap, gm, eps = args
    t0 = time.time()
    J2 = gap_J2(N, gap)
    s = system(N, J2)
    psi = rotate(s["trial"], N, theta_for(s["ground"], s["trial"], N, gm)).astype(complex)
    H = j1j2_chain(N, 1.0, J2)
    c = Case(H, float(s["E0"]), float(s["E1"]), s["ground"], psi, trial_cost(N))
    out = dict(N=N, gap_target=gap, J2=J2, gamma_target=gm)
    out.update(c.all(eps))
    out["time"] = time.time() - t0
    return out


def run_bign():
    path = os.path.join(RES, "review_bign.jsonl")
    for N in BIG_N:
        J2 = gap_J2(N, 0.40)
        system(N, J2)
        print(f"[bign] N={N} J2={J2}", flush=True)
    done = _done(path, lambda r: (r["N"], r["gap_target"], r["gamma_target"], r["eps"]))
    todo = [(N, 0.40, gm, e) for N in BIG_N for gm in (0.1, 0.5) for e in EPSS if (N, 0.40, gm, e) not in done]
    print(len(todo), "large-N cases", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("BIG_JOBS", 2))) as ex:
        for r in ex.map(_big_case, todo):
            _append(path, [r])
            print(f"[bign] N={r['N']} gamma={r['gamma']:.3f} eps={r['eps']:g} "
                  + " ".join(f"{k}={(r[k] or {}).get('cnot', float('nan')):.3e}"
                             for k in ("dyadic", "optimized", "qetu_pi2", "qetu_max"))
                  + f" ({r['time']:.0f}s)", flush=True)


if __name__ == "__main__":
    for k in (sys.argv[1] if len(sys.argv) > 1 else "mono").split(","):
        t = time.time()
        dict(trunc=run_trunc, mono=run_mono, bign=run_bign, sens2=run_sens2, eckart=run_eckart)[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)


# ----------------------------------------------------------------------------- 4. wider wrong-input study
SENS_FILTERS = ("dyadic", "optimized", "qetu")


def _sens2_case(args):
    """Plan from perturbed inputs (certified Trotter step, so only the design is tested), simulate on the
    TRUE trial and ground state. kinds:
      E0     : assumed E0 = E0 + f * gap (E1 kept), so f < 0 also widens the apparent gap
      Emax   : assumed top of spectrum E0 + f * (E_top - E0) (default: the Weyl bound)
      gamma  : assumed gamma = f * gamma, no margin (as in ext_sensitivity.json)
      gammaM : the same with a safety margin: design for gamma_assumed / 1.1"""
    kind, f, N, J2, sched = args
    from dataclasses import replace
    import run_extended as RX
    from sbpipe import Config, refine, validate
    rec, est, ans, g = RX.get_case("chain", N, J2, 2)
    H = j1j2_chain(N, 1.0, J2)
    e, a = est, ans
    if kind == "E0":
        e = replace(est, E0=est.E0 + f * est.gap, E0_err=0.0)
    elif kind == "Emax":
        top = float(eigsh(H.to_sparse(), k=1, which="LA", return_eigenvectors=False, tol=1e-10)[0])
        e = replace(est, E_max=est.E0 + f * (top - est.E0))
    elif kind in ("gamma", "gammaM"):
        gm = min(1.0, f * ans.gamma)
        a = replace(ans, gamma=gm / 1.1 if kind == "gammaM" else gm)
    try:
        plan = refine(H, e, a, Config(eps=1e-2, schedule=sched, order=2))
    except (RuntimeError, ValueError, AssertionError) as err:
        return dict(kind=kind, f=f, N=N, J2=J2, schedule=sched, plan=None, error=str(err))
    sim = validate(plan, H, ans, g)
    return dict(kind=kind, f=f, N=N, J2=J2, schedule=sched, plan=plan.summary(), sim=sim)


def run_sens2():
    import run_study as RS
    path = os.path.join(RES, "review_sens2.json")
    cases = [("E0", f, N, J2, s) for f in (-0.25, -0.1, 0.1, 0.25) for N in RS.N_VAL for J2 in RS.J2S
             for s in SENS_FILTERS]
    cases += [("Emax", f, N, J2, s) for f in (0.9, 0.95, 1.0, 1.25) for N in RS.N_VAL for J2 in RS.J2S
              for s in SENS_FILTERS]
    cases += [("gammaM", f, N, J2, s) for f in (1.0, 1.02, 1.05, 1.1) for N in RS.N_VAL for J2 in RS.J2S
              for s in SENS_FILTERS]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_sens2_case, cases))
    json.dump(out, open(path, "w"), default=_f)
    print(len(out), "wrong-input cases")


# ----------------------------------------------------------------------------- 5. certified (Eckart) step rule
DTS = [4.0 / 2 ** n for n in range(12)]                # physical step grid, units 1/J1


def _eckart_case(args):
    """DMRG-planned filter (refine) with the Trotter step chosen by a certificate a device can check:
    the Eckart bound 1 - F <= (<H> - E0) / (E1 - E0) on the post-selected state, with DMRG E0 and E1.
    Halve dt from 4/J1 until the bound is <= eps. If no step certifies the plan designed at eps (the
    filter's own residual sits too high in the spectrum), redesign at eps/2, eps/4, ... For the plan at
    eps the oracle-coarsest passing step is recorded too, for comparison."""
    kind, size, J2, eps, chi, sched = args
    import run_extended as RX
    import run_aposteriori as AP
    from sbpipe import Config, refine
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    t0 = time.time()
    rec, est, ans, _ = RX.get_case(kind, size, J2, chi)
    ans.statevector = None
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    hh = _RC.Heis(H)
    w, V = _RC.lowest(hh, k=6)
    g = V[:, 0].astype(complex)
    E0x, E1x = w[0], w[int(np.argmax(w > w[0] + 1e-8))]
    _, mps = trial_mps(pauli_to_mpo(H), chi)
    psi = np.asarray(mps.to_dense()).reshape(-1).astype(complex)
    psi /= np.linalg.norm(psi)
    G = AP._groups(H)
    E0, gap = est.E0, est.E1 - est.E0                   # what the experiment knows (DMRG)
    out = dict(kind=kind, size=size, N=N, J2=J2, eps=eps, chi=chi, schedule=sched, E0_dmrg=E0, gap_dmrg=gap,
               E0_ed=E0x, gap_ed=E1x - E0x, gamma=abs(np.vdot(g, psi)) ** 2, attempts=[])
    for frac in (1.0, 0.5, 0.25, 0.125, 0.0625):
        try:
            plan = refine(H, est, ans, Config(eps=eps * frac, schedule=sched, order=2))
        except RuntimeError:
            continue
        des, Hs, W = plan.design, plan.H_scaled, plan.window["W"]
        if des.m == 0:
            break
        h = _RC.Heis(Hs)
        t = np.abs(des.times)
        curve = []
        for dtp in DTS:
            dt = dtp * W
            k = np.maximum(1, np.ceil(t / dt - 1e-12)).astype(int)
            if curve and curve[-1]["steps"] == int(k.sum()):
                continue
            if des.schedule == "qetu":
                dm = des.meta
                x_s = dm["x_s"]
                Xt = lambda v, kk: 0.5 * (h.strang(v, dm["tau"], kk) + h.strang(v, -dm["tau"], kk))
                it = iter(k)
                v0 = psi.copy()
                v1 = Xt(v0, next(it)) / x_s
                for _ in range(2, dm["degree"] + 1):
                    v0, v1 = v1, 2.0 * Xt(v1, next(it)) / x_s - v0
                v = v1 / np.cosh(dm["degree"] * np.arccosh(1 / x_s))
            else:
                v = psi.copy()
                for ti, ph, kk in zip(des.times, des.phases, k):
                    a, b = h.strang(v, ti, int(kk)), h.strang(v, -ti, int(kk))
                    v = 0.5 * (np.exp(-1j * ph) * a + np.exp(1j * ph) * b)
            nv = float(np.vdot(v, v).real)
            u = v / math.sqrt(nv)
            inf = max(0.0, 1.0 - abs(np.vdot(g, u)) ** 2)
            en = float(np.vdot(u, hh.matvec(u)).real)
            var = []
            for O in G:
                Ou = O @ u
                var.append(max(0.0, float(np.vdot(Ou, Ou).real - np.vdot(u, Ou).real ** 2)))
            anc = 2 if des.schedule == "qetu" else 1
            fc = filter_counts(Hs, k, 2, anc_rot_per_pulse=anc)["cnot"]
            curve.append(dict(dt=dtp, steps=int(k.sum()), filter_cnot=fc, inf=inf, p=nv, E=en, var=var,
                              bound=(en - E0) / gap))
            if (en - E0) / gap <= eps or int(k.sum()) > 4000:
                break
        out["attempts"].append(dict(frac=frac, schedule=des.schedule, m=des.m, curve=curve,
                                    trial_cnot=float(ans.cost["cnot"])))
        if curve and curve[-1]["bound"] <= eps:
            break
    out["time"] = time.time() - t0
    return out


def run_eckart():
    import run_aposteriori as AP
    path = os.path.join(RES, "review_eckart.jsonl")
    done = _done(path, lambda r: (r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"]))
    todo = [c for c in AP.cases() if tuple(c) not in done]
    todo.sort(key=lambda c: (c[1] if c[0] == "chain" else 4 * c[1]))
    print(len(todo), "Eckart cases", flush=True)
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_eckart_case, todo):
            _append(path, [r])
            a = r["attempts"][-1] if r["attempts"] else None
            c = a["curve"][-1] if a and a["curve"] else {}
            print(f"[eckart] {r['kind']}{r['size']} J2={r['J2']} eps={r['eps']:g} {r['schedule']} "
                  f"frac={a and a['frac']} bound={c.get('bound', float('nan')):.2e} inf={c.get('inf', float('nan')):.2e} "
                  f"({r['time']:.0f}s)", flush=True)
