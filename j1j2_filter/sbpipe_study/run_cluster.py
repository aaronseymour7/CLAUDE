"""Large state-vector studies for a cluster (SLURM array jobs; one case per task).

Studies
  chain : dyadic SBC, optimized SBC and QETU (three query lengths) at N = 16-24, gap held at 0.40 by
          tuning J2, chi = 16 MPS trial rotated to gamma = 0.1 and 0.5, eps = 1e-2, 1e-3, 1e-4.
          Answers: how does the SBC/QETU crossover target move with N?
  cyl   : filtered route against the MPS alone on 4 x 5 and 4 x 6 cylinders (J2 = 0, 0.3): snake MPS
          trials chi = 2-64 and two non-snake trials (products of exact 4-site ring and 8-site
          two-column block ground states), eps = 1e-2 ... 1e-5.

Stages (prep must finish before case; prep is one task per system, case one task per index):
    python run_cluster.py count STUDY prep|case     # number of tasks
    python run_cluster.py prep STUDY IDX            # J2 tuning, ground state, E0, E1, trials -> cache
    python run_cluster.py case STUDY IDX            # one measured search -> results/cluster/STUDY/IDX.json
    python run_cluster.py merge STUDY               # -> results/cluster_STUDY.jsonl

Rules are those of run_gapov.py / run_empirical.py: shortest exact filter that reaches eps, then the coarsest
common Strang step (bisection) at which the Trotterised circuit still does; cost = (trial + filter CNOTs)/p.
Memory: ~3 GB at N = 20, ~10 GB at N = 24. Needs numpy, scipy, quimb and this repository (sbpipe).
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

import numpy as np

import sbpipe
from sbpipe import j1j2_chain
from sbpipe.models import j1j2_square
from sbpipe import design as D
from sbpipe.resources import filter_counts

from bigsim import Heis, rotate, lowest

RES = os.path.join(HERE, "results")
CACHE = os.environ.get("CLUSTER_CACHE", os.path.join(RES, "cache_cluster"))
OUTD = os.path.join(RES, "cluster")

CHAIN_N = [int(x) for x in os.environ.get("CHAIN_N", "16,18,20,22,24").split(",")]
CHAIN_GAMMAS = [0.1, 0.5]
EPSS = [1e-2, 1e-3, 1e-4]
CHAIN_FILTERS = ["dyadic", "optimized", "qetu_pi2", "qetu_3pi4", "qetu_max"]
J2_KNOWN = {16: 0.5553835}                        # gap 0.40 (results/gapov_tune.json)

CYL = [(Lx, J2) for Lx in map(int, os.environ.get("CYL_LX", "5,6").split(",")) for J2 in (0.0, 0.3)]
CYL_CHIS = [2, 4, 8, 16, 32, 64]
CYL_TRIALS = [f"mps{c}" for c in CYL_CHIS] + ["ring", "block"]
CYL_EPS = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]
CYL_FILTERS = ["dyadic", "optimized", "qetu_max"]
LY = 4


def dump(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o))
    os.replace(tmp, path)


def mps_cnot(N, chi):
    """Staircase circuit of a chi = 2^b MPS, each unitary by the quantum Shannon decomposition."""
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    return mps_circuit_cost(N, chi)["cnot"]


def prep_cnot(n):
    """Generic n-qubit state preparation (Plesch-Brukner): (23/24) 2^n - (3/2) 2^(n/2) + 4/3, n even."""
    return round((23 / 24) * 2 ** n - 1.5 * 2 ** (n / 2) + 4 / 3)


# ----------------------------------------------------------------------------- prep
def prep_tasks(study):
    return [("chain", N) for N in CHAIN_N] if study == "chain" else [("cyl", Lx, J2) for Lx, J2 in CYL]


def _gap(N, J2):
    w, _ = lowest(Heis(j1j2_chain(N, 1.0, J2)), k=5, tol=1e-9)
    return float(w[np.argmax(w > w[0] + 1e-8)] - w[0])


def prep_chain(N):
    p = os.path.join(CACHE, f"chain_N{N}.npz")
    if os.path.exists(p):
        return
    J2 = J2_KNOWN.get(N)
    if J2 is None:                                    # the gap rises with J2 here; bracket from N <= 16
        lo, hi = 0.50, 0.62
        for _ in range(16):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if _gap(N, mid) < 0.40 else (lo, mid)
        J2 = 0.5 * (lo + hi)
    H = j1j2_chain(N, 1.0, J2)
    w, V = lowest(Heis(H), k=6)
    i1 = int(np.argmax(w > w[0] + 1e-8))
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    _, mps = trial_mps(pauli_to_mpo(H), 16)
    t = np.asarray(mps.to_dense()).reshape(-1)
    np.savez(p, J2=J2, E0=w[0], E1=w[i1], ground=V[:, 0], trial=t / np.linalg.norm(t))


def _block_trial(Lx, J2, cols):
    """Product of exact ground states of consecutive column blocks (cols columns each, last block may be
    shorter), each block an open-in-x, periodic-in-y J1-J2 cylinder of its own."""
    parts, cost, x = [], 0, 0
    while x < Lx:
        c = min(cols, Lx - x)
        Hb = j1j2_square(c, LY, 1.0, J2)
        w, V = lowest(Heis(Hb), k=1) if c * LY > 8 else np.linalg.eigh(Hb.to_sparse().toarray())
        parts.append(np.asarray(V[:, 0]).reshape(-1))
        cost += prep_cnot(c * LY)
        x += c
    v = parts[0]
    for q in parts[1:]:
        v = np.kron(v, q)
    return v / np.linalg.norm(v), cost


def prep_cyl(Lx, J2):
    p = os.path.join(CACHE, f"cyl_{Lx}x{LY}_J2{J2:g}.npz")
    if os.path.exists(p):
        return
    H = j1j2_square(Lx, LY, 1.0, J2)
    N = H.num_qubits
    w, V = lowest(Heis(H), k=6)
    i1 = int(np.argmax(w > w[0] + 1e-8))
    g = V[:, 0]
    out = dict(E0=w[0], E1=w[i1], ground=g)
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    mpo = pauli_to_mpo(H)
    for chi in CYL_CHIS:
        _, mps = trial_mps(mpo, chi)
        t = np.asarray(mps.to_dense()).reshape(-1)
        t = t / np.linalg.norm(t)
        out[f"trial_mps{chi}"], out[f"cost_mps{chi}"] = t, mps_cnot(N, chi)
    for name, cols in (("ring", 1), ("block", 2)):
        t, c = _block_trial(Lx, J2, cols)
        out[f"trial_{name}"], out[f"cost_{name}"] = t, c
    out["gammas"] = np.array([abs(np.vdot(g, out[f"trial_{n}"])) ** 2 for n in CYL_TRIALS])
    out["generic_cnot"] = prep_cnot(N)
    np.savez(p, **out)


# ----------------------------------------------------------------------------- measured searches
class Case:
    def __init__(self, H, E0, E1, g, psi, trial_cnot):
        lo, hi = sbpipe.weyl_bounds(H)
        self.E0, self.E1, self.W = E0, E1, hi - E0
        self.delta = (E1 - E0) / self.W
        self.Hs = H.shifted_scaled(E0, self.W)
        self.h = Heis(self.Hs)
        self.g, self.psi, self.trial_cnot = g.astype(complex), psi.astype(complex), trial_cnot
        self.gamma = abs(np.vdot(self.g, self.psi)) ** 2
        self.log = []

    def inf(self, v):
        p = float(np.vdot(v, v).real)
        return max(0.0, 1 - abs(np.vdot(self.g, v)) ** 2 / p), p

    def X(self, v, t):                                # exact cos(t Hs) v
        return 0.5 * (self.h.expm(v, t, -0.02, 1.02) + self.h.expm(v, -t, -0.02, 1.02))

    def Xt(self, v, t, k):                            # Trotterised controlled pulse, ancilla-0 block
        return 0.5 * (self.h.strang(v, t, k) + self.h.strang(v, -t, k))

    def sbc_run(self, times, phases, dt):
        k = np.maximum(1, np.ceil(np.abs(times) / dt - 1e-12)).astype(int)
        v = self.psi.copy()
        for t, ph, kk in zip(times, phases, k):
            a, b = self.h.strang(v, t, int(kk)), self.h.strang(v, -t, int(kk))
            v = 0.5 * (np.exp(-1j * ph) * a + np.exp(1j * ph) * b)
        inf, p = self.inf(v)
        self.log.append((float(dt), int(k.sum()), inf))
        return k, inf, p

    def sbc_step(self, times, phases, eps):
        t = np.abs(times)
        best, dt = None, float(t.max())
        for _ in range(14):
            k, inf, p = self.sbc_run(times, phases, dt)
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
                if inf <= eps:
                    lo, best = mid, (mid, k, inf, p)
                else:
                    hi = mid
        return best

    def _sbc_result(self, times, phases, eps, **extra):
        tr = self.sbc_step(times, phases, eps)
        if not tr:
            return None
        dt, k, inf, p = tr
        fc = filter_counts(self.Hs, k, 2)["cnot"]
        return dict(extra, times=np.asarray(times).tolist(), phases=np.asarray(phases).tolist(), m=len(times),
                    T=float(np.abs(times).sum()), steps=int(k.sum()), dt=dt,
                    dt_phys=dt / self.W, inf=inf, p=p, cnot_filter=fc, cnot=(fc + self.trial_cnot) / p)

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
            if cands and tau * (2 ** K - 1) > 2 * min(c[0] for c in cands):
                break
        best = None
        for T, K, r in sorted(cands)[:n_cand]:
            times = np.tile(tau * 2.0 ** np.arange(K), r)
            b = self._sbc_result(times, np.zeros(len(times)), eps, K=K, r=r)
            if b and (best is None or b["cnot"] < best["cnot"]):
                best = b
        return best

    def optimized(self, eps):
        best = None
        for seed in (0, 1):
            for frac in (1.0, 0.5):
                des = D.design_filter("optimized", self.delta, 0.0, eps * frac, self.gamma, seed=seed)
                if des is None:
                    continue
                b = self._sbc_result(np.asarray(des.times), np.asarray(des.phases), eps, seed=seed, frac=frac)
                if b and (best is None or b["cnot"] < best["cnot"]):
                    best = b
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
                v0 = self.psi.copy()
                v1 = self.Xt(v0, tau, k) / x_s
                for _ in range(2, dd + 1):
                    v0, v1 = v1, 2.0 * self.Xt(v1, tau, k) / x_s - v0
                v = v1 / np.cosh(dd * np.arccosh(1 / x_s))     # p(x) = T_d(x/x_s) / T_d(1/x_s)
                inf, p = self.inf(v)
                self.log.append((tau / k, dd * k, inf))
                if inf <= eps:
                    fc = filter_counts(self.Hs, [k] * dd, 2, anc_rot_per_pulse=2)["cnot"]
                    best = dict(d_exact=dmin, d=dd, k=k, tau=tau, inf=inf, p=p, cnot_filter=fc,
                                cnot=(fc + self.trial_cnot) / p, dt_phys=tau / k / self.W, steps=dd * k)
                    break
            if best is not None and best["k"] == 1:
                break
        return best

    def run(self, filt, eps):
        if filt == "dyadic":
            return self.dyadic(eps)
        if filt == "optimized":
            return self.optimized(eps)
        tau = {"qetu_pi2": math.pi / 2, "qetu_3pi4": 0.75 * math.pi,
               "qetu_max": math.pi / (1 + self.delta)}[filt]
        return self.qetu(eps, tau)


# ----------------------------------------------------------------------------- cases
def case_tasks(study):
    if study == "chain":
        return [(N, gm, e, f) for N in CHAIN_N for gm in CHAIN_GAMMAS for e in EPSS for f in CHAIN_FILTERS]
    out = []
    for Lx, J2 in CYL:
        for tr in CYL_TRIALS:
            for e in CYL_EPS:
                for f in CYL_FILTERS:
                    if f == "optimized" and e not in EPSS:
                        continue
                    out.append((Lx, J2, tr, e, f))
    return out


def run_case(study, idx):
    a = case_tasks(study)[idx]
    path = os.path.join(OUTD, study, f"{idx}.json")
    if os.path.exists(path):
        return
    t0 = time.time()
    if study == "chain":
        N, gm, eps, filt = a
        z = np.load(os.path.join(CACHE, f"chain_N{N}.npz"))
        J2 = float(z["J2"])
        g, t = z["ground"].astype(complex), z["trial"].astype(complex)
        f = lambda th: abs(np.vdot(g, rotate(t, N, th))) ** 2
        lo, hi = 0.0, math.pi / 2 if f(math.pi / 2) < gm else math.pi
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) > gm else (lo, mid)
        psi = rotate(t, N, 0.5 * (lo + hi))
        H = j1j2_chain(N, 1.0, J2)
        c = Case(H, float(z["E0"]), float(z["E1"]), g, psi, mps_cnot(N, 16))
        rec = dict(study=study, idx=idx, N=N, J2=J2, gamma_target=gm, eps=eps, filter=filt)
    else:
        Lx, J2, tr, eps, filt = a
        z = np.load(os.path.join(CACHE, f"cyl_{Lx}x{LY}_J2{J2:g}.npz"))
        H = j1j2_square(Lx, LY, 1.0, J2)
        g, psi = z["ground"].astype(complex), z[f"trial_{tr}"].astype(complex)
        c = Case(H, float(z["E0"]), float(z["E1"]), g, psi, float(z[f"cost_{tr}"]))
        rec = dict(study=study, idx=idx, Lx=Lx, Ly=LY, N=Lx * LY, J2=J2, trial=tr, eps=eps, filter=filt,
                   generic_cnot=int(z["generic_cnot"]))
    rec.update(E0=c.E0, gap=c.E1 - c.E0, W=c.W, delta=c.delta, gamma=c.gamma, trial_cnot=c.trial_cnot)
    if 1 - c.gamma <= eps:
        rec.update(mps_alone=True, result=None)
    else:
        rec.update(mps_alone=False, result=c.run(filt, eps))
    rec.update(trotter_log=c.log, time=time.time() - t0)
    dump(path, rec)


def merge(study):
    d = os.path.join(OUTD, study)
    rows = [json.load(open(os.path.join(d, f))) for f in sorted(os.listdir(d), key=lambda s: int(s.split(".")[0]))
            if f.endswith(".json")]
    with open(os.path.join(RES, f"cluster_{study}.jsonl"), "w") as fh:
        fh.write("".join(json.dumps(r) + "\n" for r in rows))
    print(len(rows), "of", len(case_tasks(study)), "cases merged")


if __name__ == "__main__":
    os.makedirs(CACHE, exist_ok=True)
    cmd, study = sys.argv[1], sys.argv[2]
    if cmd == "count":
        print(len(prep_tasks(study) if sys.argv[3] == "prep" else case_tasks(study)))
    elif cmd == "prep":
        a = prep_tasks(study)[int(sys.argv[3])]
        prep_chain(a[1]) if a[0] == "chain" else prep_cyl(a[1], a[2])
    elif cmd == "case":
        run_case(study, int(sys.argv[3]))
    elif cmd == "merge":
        merge(study)
