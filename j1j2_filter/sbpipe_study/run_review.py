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
class Case:
    def __init__(self, H, E0, E1, g, psi, trial_cnot):
        self.H, self.g, self.psi, self.trial_cnot = H, g, psi, trial_cnot
        lo, hi = sbpipe.weyl_bounds(H)
        self.E0, self.E1, self.W = E0, E1, hi - E0
        self.delta = (E1 - E0) / self.W
        self.Hs = H.shifted_scaled(E0, self.W)
        self.Hsp = self.Hs.to_sparse().tocsr()
        self.S = Strang(self.Hs)
        self.gamma = abs(np.vdot(g, psi)) ** 2

    def inf(self, v):
        p = float(np.vdot(v, v).real)
        return max(0.0, 1 - abs(np.vdot(self.g, v)) ** 2 / p), p

    def X(self, v, t):
        return 0.5 * (sla.expm_multiply(-1j * t * self.Hsp, v) + sla.expm_multiply(1j * t * self.Hsp, v))

    def sbc_run(self, times, phases, dt):
        k = np.maximum(1, np.ceil(np.abs(times) / dt - 1e-12)).astype(int)
        v = ftrot(self.S, self.psi, times, phases, k)
        inf, p = self.inf(v)
        return k, inf, p

    def sbc_step(self, times, phases, eps, log=None):
        """Coarsest common step (scaled units) at which the Trotterised SBC filter reaches eps."""
        t = np.abs(times)
        best, dt = None, float(t.max())
        for _ in range(14):
            k, inf, p = self.sbc_run(times, phases, dt)
            if log is not None:
                log.append((dt, inf))
            if inf <= eps:
                best = (dt, k, inf, p)
                break
            dt /= 2
        if best is None:
            return None
        lo, hi = best[0], 2 * best[0]
        if lo < t.max():
            for _ in range(5):
                mid = 0.5 * (lo + hi)
                k, inf, p = self.sbc_run(times, phases, mid)
                if log is not None:
                    log.append((mid, inf))
                if inf <= eps:
                    lo, best = mid, (mid, k, inf, p)
                else:
                    hi = mid
        return best

    def dyadic(self, eps, n_cand=2):
        tau = math.pi / (1 + self.delta)
        cands = []
        for K in range(1, 16):
            v = self.psi.copy()
            for r in range(1, 200):
                for k in range(K):
                    v = self.X(v, tau * 2 ** k)
                if self.inf(v)[0] <= eps:
                    cands.append((r * tau * (2 ** K - 1), K, r))
                    break
                if cands and r * tau * (2 ** K - 1) > 2 * min(c[0] for c in cands):
                    break
        cands.sort()
        best = None
        for T, K, r in cands[:n_cand]:
            times = np.tile(tau * 2.0 ** np.arange(K), r)
            tr = self.sbc_step(times, np.zeros(len(times)), eps)
            if tr:
                dt, k, inf, p = tr
                fc = filter_counts(self.Hs, k, 2)["cnot"]
                c = (fc + self.trial_cnot) / p
                if best is None or c < best["cnot"]:
                    best = dict(T=T, K=K, r=r, m=len(times), steps=int(k.sum()), dt=dt, dt_phys=dt / self.W,
                                inf=inf, p=p, cnot_filter=fc, cnot=c, times=times.tolist())
        return best

    def optimized(self, eps):
        best = None
        for seed in (0, 1):
            for frac in (1.0, 0.5):
                des = D.design_filter("optimized", self.delta, 0.0, eps * frac, self.gamma, seed=seed)
                if des is None:
                    continue
                tr = self.sbc_step(des.times, des.phases, eps)
                if tr:
                    dt, k, inf, p = tr
                    fc = filter_counts(self.Hs, k, 2)["cnot"]
                    c = (fc + self.trial_cnot) / p
                    if best is None or c < best["cnot"]:
                        best = dict(T=des.T, m=des.m, steps=int(k.sum()), dt=dt, dt_phys=dt / self.W, inf=inf,
                                    p=p, cnot_filter=fc, cnot=c, seed=seed, frac=frac,
                                    times=np.asarray(des.times).tolist(), phases=np.asarray(des.phases).tolist())
            if best is not None and best["p"] >= 0.5 * self.gamma:
                break
        return best

    def qetu(self, eps, tau):
        x_s = math.cos(tau * self.delta)
        v0 = self.psi.copy()
        v1 = self.X(v0, tau) / x_s
        d, dmin = 1, None
        while d < 4000:
            v0, v1 = v1, 2.0 * self.X(v1, tau) / x_s - v0
            d += 1
            if d % 2 == 0 and self.inf(v1)[0] <= eps:
                dmin = d
                break
        if dmin is None:
            return None
        best = None
        for dd in range(dmin, dmin + 40, 2):
            for k in (1, 2, 3, 4):
                if best is not None and dd * k >= best["d"] * best["k"]:
                    break
                v = qtrot(self.S, self.psi, dd, tau, x_s, [k] * dd)
                inf, p = self.inf(v)
                if inf <= eps:
                    fc = filter_counts(self.Hs, [k] * dd, 2, anc_rot_per_pulse=2)["cnot"]
                    best = dict(d_exact=dmin, d=dd, k=k, tau=tau, inf=inf, p=p, cnot_filter=fc,
                                cnot=(fc + self.trial_cnot) / p, dt_phys=tau / k / self.W)
                    break
            if best is not None and best["k"] == 1:
                break
        return best

    def all(self, eps, scheds=("dyadic", "optimized", "qetu")):
        out = dict(E0=self.E0, gap=self.E1 - self.E0, W=self.W, delta=self.delta, gamma=self.gamma,
                   trial_cnot=self.trial_cnot, eps=eps)
        if "dyadic" in scheds:
            out["dyadic"] = self.dyadic(eps)
        if "optimized" in scheds:
            out["optimized"] = self.optimized(eps)
        if "qetu" in scheds:
            for name, tau in (("qetu_pi2", math.pi / 2), ("qetu_3pi4", 0.75 * math.pi),
                              ("qetu_max", math.pi / (1 + self.delta))):
                out[name] = self.qetu(eps, tau)
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
        if abs(math.log(gam / gm)) < math.log(1.5):
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
        times = np.asarray(b["times"])
        phases = np.asarray(b.get("phases", np.zeros(len(times))))
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
        x_s = math.cos(q["tau"] * c.delta)
        scan = []
        for k in range(1, 7):
            v = qtrot(c.S, psi, q["d"], q["tau"], x_s, [k] * q["d"])
            scan.append((k, c.inf(v)[0]))
        out["qetu"] = dict(d=q["d"], scan=scan)
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
        dict(trunc=run_trunc, mono=run_mono, bign=run_bign)[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)
