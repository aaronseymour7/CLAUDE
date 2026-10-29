"""Empirical filter behaviour on the exact spectrum (N <= 16): how the gap and the trial overlap set
what the filter actually needs, as opposed to what the a-priori certificate asks for.

For each chain the trial state is decomposed on the exact eigenbasis (block-diagonal in S^z, dense
diagonalisation per block). Any filter p(H_s) then acts on the weights w_k analytically, so the
smallest filter that reaches a target infidelity can be found exactly. Trotter steps are found by
statevector simulation (N <= 12), as in run_study's oracle.

Families
  fixedJ2  : J2 = 0.2411 for every N (the gap shrinks with N)
  matched  : J2 chosen per N so that the true gap E1 - E0 equals GAP_TARGET (ED, triplet checked)
Trial states: the variationally optimal MPS at chi = 1, 2, 4, and synthetic states
sqrt(gamma)|E0> + sqrt(1-gamma)|xi>, xi the excited part of the chi = 2 MPS (fixed spectral shape).

    python sbpipe_study/run_empirical.py match     # J2(N) with a fixed gap
    python sbpipe_study/run_empirical.py spectra   # eigen-decomposition of every trial state (cached)
    python sbpipe_study/run_empirical.py filters   # smallest QETU / dyadic filter that reaches eps
    python sbpipe_study/run_empirical.py trotter   # smallest Trotter step count at that filter, and the certified plan
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import gc
import json
import math

import scipy.linalg
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.sparse.linalg import eigsh

import sbpipe
from sbpipe import j1j2_chain
from sbpipe import design as D

import run_study as RS

RES = RS.RES
CACHEE = os.path.join(RES, "cache_emp")
os.makedirs(CACHEE, exist_ok=True)
JOBS = RS.JOBS
dump, load = RS.dump, RS.load

N_EMP = [8, 10, 12, 14, 16]
GAP_TARGET = 0.40
J2_FIXED = 0.2411
CHIS = [1, 2, 4]
GAMMAS = [0.03, 0.1, 0.3, 0.5, 0.7, 0.9]
EPSS = [1e-2, 1e-3, 1e-4]


def _gap(N, J2):
    w = eigsh(j1j2_chain(N, 1.0, J2).to_sparse(), k=4, which="SA", tol=1e-11,
              return_eigenvectors=False)
    w = np.sort(w)
    return float(w[np.argmax(w > w[0] + 1e-8)] - w[0])


def _match(N):
    lo, hi = 0.0, 0.65
    if _gap(N, lo) > GAP_TARGET:            # already above the target at J2 = 0: bisect on [0, 0.3]
        lo, hi = 0.25, 0.0
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if (_gap(N, mid) < GAP_TARGET) == (lo < hi):
            lo = mid
        else:
            hi = mid
    return N, round(0.5 * (lo + hi), 6)


def run_match():
    with ProcessPoolExecutor(JOBS) as ex:
        out = {N: J2 for N, J2 in ex.map(_match, N_EMP)}
    for N, J2 in out.items():
        print(f"[match] N={N} J2={J2} gap={_gap(N, J2):.5f}", flush=True)
    dump("emp_match.json", out)


N_SCAN = [10, 12, 14]
GAP_SCAN = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55]


def _scan(args):
    """J2 in [0.25, 0.7] with gap = target (the gap is flat below J2 ~ 0.3 and rises above it)."""
    N, target = args
    lo, hi = 0.25, 0.7
    if not _gap(N, lo) < target < _gap(N, hi):
        return N, target, None
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if _gap(N, mid) < target:
            lo = mid
        else:
            hi = mid
    return N, target, round(0.5 * (lo + hi), 6)


def run_gapscan():
    """Physical gap scan at fixed N: J2 chosen for each target gap, then spectra as for the families."""
    with ProcessPoolExecutor(JOBS) as ex:
        out = [dict(N=N, gap=g, J2=J2) for N, g, J2 in ex.map(_scan, [(N, g) for N in N_SCAN for g in GAP_SCAN])]
    for r in out:
        print(f"[gapscan] N={r['N']} gap={r['gap']} J2={r['J2']}", flush=True)
    dump("emp_gapscan.json", out)


def families():
    m = load("emp_match.json")
    out = [("fixedJ2", N, J2_FIXED) for N in N_EMP] + [("matched", N, float(m[str(N)])) for N in N_EMP]
    out += [("gapscan", r["N"], r["J2"]) for r in (load("emp_gapscan.json") if os.path.exists(os.path.join(RES, "emp_gapscan.json")) else [])
            if r["J2"] is not None]
    return out


# ----------------------------------------------------------------------------- spectra
def _key(N, J2):
    return f"N{N}_J2{J2:.6g}"


def _spectrum(args):
    """Exact eigen-decomposition of H in every S^z block, projected on the trial states."""
    fam, N, J2 = args
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    p = os.path.join(CACHEE, _key(N, J2) + ".npz")
    if os.path.exists(p):
        return fam, N, J2, 0.0
    t0 = time.time()
    H = j1j2_chain(N, 1.0, J2)
    Hsp = H.to_sparse().tocsr()
    mpo = pauli_to_mpo(H)
    trials = {}
    for chi in CHIS:
        _, psi = trial_mps(mpo, chi)
        v = np.asarray(psi.to_dense()).reshape(-1)
        trials[chi] = v / np.linalg.norm(v)
    idx = np.arange(2 ** N)
    pop = np.array([bin(i).count("1") for i in idx])
    E, W = [], {c: [] for c in CHIS}
    vecs_low = []                                   # keep the ground state
    for n_up in range(N + 1):
        b = idx[pop == n_up]
        Hb = Hsp[b][:, b].toarray()
        if np.abs(Hb.imag).max() < 1e-12 if np.iscomplexobj(Hb) else True:
            Hb = np.ascontiguousarray(Hb.real)
        e, U = scipy.linalg.eigh(Hb, overwrite_a=True, check_finite=False, driver="evr")
        del Hb
        E.append(e)
        for c in CHIS:
            W[c].append(np.abs(U.T.conj() @ trials[c][b]) ** 2)
        vecs_low.append((e[0], b, U[:, 0].copy()))
        del U
        gc.collect()
    E = np.concatenate(E)
    o = np.argsort(E)
    out = dict(E=E[o])
    for c in CHIS:
        out[f"w{c}"] = np.concatenate(W[c])[o]
    e0, b, u = min(vecs_low, key=lambda t: t[0])
    g = np.zeros(2 ** N, complex)
    g[b] = u
    out["ground"] = g
    for c in CHIS:
        out[f"trial{c}"] = trials[c]
    np.savez(p, **out)
    return fam, N, J2, time.time() - t0


def run_spectra():
    with ProcessPoolExecutor(JOBS) as ex:
        for fam, N, J2, dt in ex.map(_spectrum, families()):
            print(f"[spectra] {fam} N={N} J2={J2:g} {dt:.0f}s", flush=True)


def spectral_weights(N, J2, trial):
    """(E, w) for trial = chi (int) or ('synthetic', gamma): weights on the exact eigenstates."""
    z = np.load(os.path.join(CACHEE, _key(N, J2) + ".npz"))
    E = z["E"]
    if isinstance(trial, tuple):
        gamma = trial[1]
        w2 = z["w2"].copy()
        exc = w2.copy()
        exc[0] = 0.0
        exc /= exc.sum()                            # excited part of the chi = 2 MPS, normalised
        w = (1 - gamma) * exc
        w[0] = gamma
        return E, w
    return E, z[f"w{trial}"]


# ----------------------------------------------------------------------------- filters on the exact spectrum
def _infid(w, f):
    a = w * f ** 2
    return 1.0 - a[0] / a.sum(), a.sum()


def _qetu_vals(d, x_s, Es):
    return D.qetu_values(d, math.pi / 2, x_s, Es)


def _min_qetu(Es, w, delta, eps):
    """Smallest even degree with the Chebyshev filter designed for the TRUE scaled gap reaching eps."""
    x_s = math.cos(math.pi / 2 * delta)
    d = 0
    while d < 200000:
        d += 2
        inf, p = _infid(w, _qetu_vals(d, x_s, Es))
        if inf <= eps:
            return d, p
    return None, None


def _min_dyadic(Es, w, delta, eps):
    """Smallest total time over (K, r) of the dyadic ladder reaching eps."""
    tau = math.pi / (1 + delta)
    best = None
    for K in range(1, 30):
        blk = np.ones_like(Es)
        for k in range(K):
            blk = blk * np.cos(2 ** k * Es * tau)
        for r in range(1, 400):
            inf, p = _infid(w, blk ** r)
            if inf <= eps:
                T = r * tau * (2 ** K - 1)
                if best is None or T < best[0]:
                    best = (T, K, r, p)
                break
            if best is not None and r * tau * (2 ** K - 1) > best[0]:
                break
    return best


def _filter_case(args):
    fam, N, J2, trial, eps = args
    E, w = spectral_weights(N, J2, trial)
    lo, hi = sbpipe.weyl_bounds(j1j2_chain(N, 1.0, J2))
    Wr = hi - E[0]
    Es = (E - E[0]) / Wr
    gap = float(E[np.argmax(E > E[0] + 1e-8)] - E[0])
    delta = gap / Wr
    gamma = float(w[0])
    # where the excited weight sits: weight-averaged excitation energy and weight within 2x the gap
    exc = w[1:] / max(w[1:].sum(), 1e-300)
    e_mean = float(np.sum(exc * (E[1:] - E[0])))
    near = float(np.sum(exc[(E[1:] - E[0]) <= 2 * gap]))
    out = dict(family=fam, N=N, J2=J2, trial=trial if isinstance(trial, int) else f"g{trial[1]}",
               eps=eps, gamma=gamma, gap=gap, W=Wr, delta=delta, e_mean_exc=e_mean, weight_near_gap=near)
    if 1 - gamma <= eps:
        return dict(out, d_min=0, d_cert=0)
    d, p = _min_qetu(Es, w, delta, eps)
    eta_t = D.eta_target(eps, gamma)
    dc = D.qetu(delta, 0.0, eta_t)
    dy = _min_dyadic(Es, w, delta, eps)
    dyc = D.dyadic(delta, 0.0, eta_t)
    return dict(out, d_min=d, p_succ_min=p, d_cert=dc.m if dc else None,
                T_min_qetu=d * math.pi / 2 if d else None, T_cert_qetu=dc.T if dc else None,
                T_min_dyadic=dy[0] if dy else None, T_cert_dyadic=dyc.T if dyc else None)


def run_filters():
    cases = []
    for fam, N, J2 in families():
        for eps in EPSS:
            cases += [(fam, N, J2, chi, eps) for chi in CHIS]
            cases += [(fam, N, J2, ("synthetic", g), eps) for g in GAMMAS]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_filter_case, cases, chunksize=4))
    dump("emp_filters.json", out)


# ----------------------------------------------------------------------------- Trotter steps actually needed (N <= 12)
def _trotter_case(args):
    fam, N, J2, trial, eps = args
    from sbpipe.pauli import PauliSum
    from sbpipe.simulate import qetu_trotter, state_metrics
    from sbpipe.resources import filter_counts
    from sbpipe.trotter import beta2, steps_for_budget
    r = _filter_case(args)
    if not r.get("d_min"):
        return dict(r, k_min=0)
    z = np.load(os.path.join(CACHEE, _key(N, J2) + ".npz"))
    g = z["ground"]
    if isinstance(trial, tuple):
        t = z["trial2"]
        xi = t - np.vdot(g, t) * g
        xi /= np.linalg.norm(xi)
        ph = np.vdot(g, t) / abs(np.vdot(g, t))
        psi = np.sqrt(trial[1]) * g * ph + np.sqrt(1 - trial[1]) * xi
    else:
        psi = z[f"trial{trial}"]
    H = j1j2_chain(N, 1.0, J2)
    E0 = z["E"][0]
    Hs = H.shifted_scaled(E0, r["W"])
    d, x_s = r["d_min"], math.cos(math.pi / 2 * r["delta"])
    # the target is met by the exact filter with slack (inf_exact <= eps); find the smallest common k
    # per query with the Trotterised filter still at <= eps
    def inf_k(k):
        v = qetu_trotter(Hs, psi, d, math.pi / 2, x_s, [k] * d)
        return state_metrics(v, g)["infidelity"], float(np.vdot(v, v).real)
    k = 1
    inf, p = inf_k(k)
    while inf > eps and k < 64:                     # k is small in practice: plain upward search
        k += 1
        inf, p = inf_k(k)
    # a slightly higher degree at one step per query can be cheaper than d_min at k > 1
    d0 = d
    for dd in range(d0 + 2, int(d0 * k) + 1, 2):
        if dd * 1 >= d0 * k:
            break
        def inf_k(kk, dd=dd):
            v = qetu_trotter(Hs, psi, dd, math.pi / 2, x_s, [kk] * dd)
            return state_metrics(v, g)["infidelity"], float(np.vdot(v, v).real)
        i1, p1 = inf_k(1)
        if i1 <= eps:
            d, k, inf, p = dd, 1, i1, p1
            break
    cnot = filter_counts(Hs, [k] * d, 2, anc_rot_per_pulse=2)["cnot"]
    # what the certificate asks for on the same input (true gap, true overlap, QETU, order 2)
    from sbpipe.pipeline import Config, _plan_for, trotter_constant
    # the pipeline's own choice: best split of eps between leakage and Trotter/synthesis
    cfg = Config(eps=eps, schedule="qetu", order=2)
    cert, const = {}, trotter_constant(Hs, 2)
    for f in cfg.leak_fracs:
        des = D.design_filter("qetu", r["delta"], 0.0, f * eps, r["gamma"])
        pc = _plan_for(H, Hs, None, des, r["gamma"], cfg, const, {}) if des else None
        if pc and (not cert or pc[2]["expected_cnot"] < cert["cert_cnot_per_success"]):
            st, ce, res = pc
            cert = dict(cert_d=des.m, cert_steps=int(np.sum(st)), cert_k_per_query=float(np.mean(st)),
                        cert_cnot_per_success=res["expected_cnot"], cert_p_succ=ce["p_succ_bound"],
                        cert_leak_frac=f)
    return dict(r, d_trot=d, k_min=k, steps_min=k * d, inf_trotter=inf, p_succ_trotter=p,
                cnot_filter=cnot, cnot_per_success=cnot / p, **cert)


def run_trotter():
    cases = []
    for fam, N, J2 in families():
        for eps in EPSS:
            cases += [(fam, N, J2, chi, eps) for chi in CHIS]
            cases += [(fam, N, J2, ("synthetic", g), eps) for g in GAMMAS]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_trotter_case, cases):
            print(f"[trotter] {r['family']} N={r['N']} {r['trial']} eps={r['eps']:g} d={r.get('d_min')} "
                  f"k={r.get('k_min')} cnot/succ={r.get('cnot_per_success', float('nan')):.2e}", flush=True)
            out.append(r)
    dump("emp_trotter.json", out)


# ----------------------------------------------------------------------------- the empirical filter on the regime systems
def _emp_regime_case(args):
    """Filter cost from the two empirical laws (QETU degree from the true-gap design at the full target,
    one Strang step per query), checked by state-vector simulation where a trial vector exists."""
    kind, size, J2, chi, eps = args
    import run_extended as RX
    from sbpipe.simulate import qetu_trotter, state_metrics
    from sbpipe.resources import filter_counts
    rec, est, ans, g = RX.get_case(kind, size, J2, chi)
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    if kind == "chain" and ans.statevector is None and N <= 16:
        p = os.path.join(CACHEE, _key(N, J2) + ".npz")
        if os.path.exists(p):
            z = np.load(p)
            if f"trial{chi}" in z:
                ans.statevector, g = z[f"trial{chi}"], z["ground"]
    out = dict(kind=kind, size=size, N=N, J2=J2, chi=chi, eps=eps, gamma=ans.gamma,
               trial_cnot=ans.cost["cnot"])
    if 1 - ans.gamma <= eps:
        return dict(out, d=0, cnot=ans.cost["cnot"], mps_alone=True)
    lo, hi = sbpipe.weyl_bounds(H)
    win = est.window(hi)
    Hs = H.shifted_scaled(win["shift"], win["W"])
    des = D.qetu(win["delta"], win["w0"], D.eta_target(eps, ans.gamma))
    if des is None:
        return dict(out, d=None)
    d = des.m
    sim = None
    if ans.statevector is not None and g is not None:
        x_s = math.cos(math.pi / 2 * win["delta"])
        for _ in range(40):
            v = qetu_trotter(Hs, ans.statevector, d, math.pi / 2, x_s, [1] * d)
            m = state_metrics(v, g)
            sim = dict(d=d, infidelity=m["infidelity"], p_succ=float(np.vdot(v, v).real))
            if m["infidelity"] <= eps:
                break
            d += 2
    p = sim["p_succ"] if sim else ans.gamma
    cnot = (filter_counts(Hs, [1] * d, 2, anc_rot_per_pulse=2)["cnot"] + ans.cost["cnot"]) / p
    return dict(out, d=d, d_model=des.m, sim=sim, cnot=cnot, mps_alone=False)


def run_emp_regime():
    import run_extended as RX
    cases = [("chain", N, 0.2411, chi, eps) for N in (16, 32, 48) for chi in (2, 4, 8, 16, 32)
             for eps in (1e-2, 1e-3, 1e-4)]
    cases += [("cyl", Lx, J2, chi, eps) for Lx in RX.LX for J2 in RX.J2S_2D for chi in RX.CHIS_2D
              for eps in (1e-2, 1e-3, 1e-4)]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_emp_regime_case, cases):
            s = r.get("sim")
            print(f"[emp_regime] {r['kind']} N={r['N']} chi={r['chi']} eps={r['eps']:g} d={r.get('d')} "
                  f"cnot={r.get('cnot') or float('nan'):.2e}"
                  + (f" sim inf={s['infidelity']:.1e} (d_model={r['d_model']})" if s else ""), flush=True)
            out.append(r)
    dump("emp_regime.json", out)


# ----------------------------------------------------------------------------- cost to reach a given infidelity
EPS_GRID = np.logspace(-1, -7, 61)


def _frontier_case(args):
    """Empirical filter cost (degree law, one step per query) on a fine eps grid for one trial."""
    kind, size, J2, chi = args
    import run_extended as RX
    from sbpipe.resources import filter_counts
    rec, est, ans, _ = RX.get_case(kind, size, J2, chi)
    H = RX.model(kind, size, J2)
    lo, hi = sbpipe.weyl_bounds(H)
    win = est.window(hi)
    Hs = H.shifted_scaled(win["shift"], win["W"])
    per_query = filter_counts(Hs, [1, 1], 2, anc_rot_per_pulse=2)["cnot"] - \
        filter_counts(Hs, [1], 2, anc_rot_per_pulse=2)["cnot"]
    base = filter_counts(Hs, [1], 2, anc_rot_per_pulse=2)["cnot"] - per_query
    pts = []
    for eps in EPS_GRID:
        if 1 - ans.gamma <= eps:
            pts.append((float(eps), 0, float(ans.cost["cnot"])))
            continue
        des = D.qetu(win["delta"], win["w0"], D.eta_target(eps, ans.gamma))
        if des is None:
            continue
        pts.append((float(eps), des.m, (base + per_query * des.m + ans.cost["cnot"]) / ans.gamma))
    return dict(kind=kind, size=size, N=H.num_qubits, J2=J2, chi=chi, gamma=ans.gamma,
                trial_cnot=ans.cost["cnot"], gap=est.E1 - est.E0, delta=win["delta"], points=pts)


def run_frontier():
    import run_extended as RX
    cases = [("chain", N, 0.2411, chi) for N in (16, 32, 48) for chi in (1, 2, 4, 8, 16, 32)]
    cases += [("cyl", Lx, J2, chi) for Lx in RX.LX for J2 in RX.J2S_2D for chi in RX.CHIS_2D]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_frontier_case, cases):
            out.append(r)
    dump("emp_frontier.json", out)


# ----------------------------------------------------------------------------- the Stetcu-Baroni-Carlson filter, empirically
SBC_SCHEDS = ("dyadic", "optimized")


def _sbc_trotter(H, Hs, psi, g, des, eps):
    """Largest common step dt (scaled units, k_i = ceil(t_i / dt)) at which the Trotterised SBC filter
    still reaches eps; returns (dt, steps, infidelity, p_succ)."""
    from sbpipe.simulate import filter_trotter, state_metrics
    t = np.abs(des.times)

    def run(dt):
        k = np.maximum(1, np.ceil(t / dt - 1e-12)).astype(int)
        v = filter_trotter(Hs, psi, des.times, des.phases, k, 2)
        m = state_metrics(v, g)
        return k, m["infidelity"], float(np.vdot(v, v).real)

    best = None
    dt = float(t.max())                             # one step for the longest pulse, then halve
    for _ in range(14):
        k, inf, p = run(dt)
        if inf <= eps:
            best = (dt, k, inf, p)
            break
        dt /= 2
    if best is None:
        return None
    lo, hi = best[0], 2 * best[0]                   # refine between the passing dt and the failing 2 dt
    if lo < t.max():
        for _ in range(5):
            mid = 0.5 * (lo + hi)
            k, inf, p = run(mid)
            if inf <= eps:
                lo, best = mid, (mid, k, inf, p)
            else:
                hi = mid
    return best


def _sbc_case(args):
    """SBC filter designed by the pipeline at the full target with the true gap; its real leakage on the
    exact spectrum and the Trotter steps it really needs, against the certified plan."""
    fam, N, J2, trial, eps, sched = args
    from sbpipe.resources import filter_counts
    from sbpipe.pipeline import Config, _plan_for, trotter_constant
    r = _filter_case((fam, N, J2, trial, eps))
    base = {k: r[k] for k in ("family", "N", "J2", "trial", "eps", "gamma", "gap", "W", "delta")}
    base["schedule"] = sched
    if not r.get("d_min"):
        return dict(base, skip=True)
    E, w = spectral_weights(N, J2, trial)
    Es = (E - E[0]) / r["W"]
    z = np.load(os.path.join(CACHEE, _key(N, J2) + ".npz"))
    g = z["ground"]
    if isinstance(trial, tuple):
        t2 = z["trial2"]
        xi = t2 - np.vdot(g, t2) * g
        xi /= np.linalg.norm(xi)
        ph = np.vdot(g, t2) / abs(np.vdot(g, t2))
        psi = np.sqrt(trial[1]) * g * ph + np.sqrt(1 - trial[1]) * xi
    else:
        psi = z[f"trial{trial}"]
    H = j1j2_chain(N, 1.0, J2)
    Hs = H.shifted_scaled(z["E"][0], r["W"])
    cfg = Config(eps=eps, schedule=sched, order=2)
    const = trotter_constant(Hs, 2)
    # certified: the pipeline's best split of eps
    cert = None
    fracs = cfg.leak_fracs if sched != "optimized" else (0.25,)   # the optimizer takes ~40 s per design
    for f in fracs:
        des = D.design_filter(sched, r["delta"], 0.0, f * eps, r["gamma"])
        pc = _plan_for(H, Hs, None, des, r["gamma"], cfg, const, {}) if des else None
        if pc and (cert is None or pc[2]["expected_cnot"] < cert["cert_cnot_per_success"]):
            cert = dict(cert_T=des.T, cert_m=des.m, cert_steps=int(np.sum(pc[0])),
                        cert_cnot_per_success=pc[2]["expected_cnot"])
    # empirical: the design for the full eps (all of it to leakage), real leakage, real Trotter steps
    des = D.design_filter(sched, r["delta"], 0.0, eps, r["gamma"])
    if des is None:
        return dict(base, **(cert or {}), design=None)
    f = des.values(Es)
    a = w * f ** 2
    leak_exact = 1 - a[0] / a.sum()
    tr = _sbc_trotter(H, Hs, psi, g, des, eps)
    out = dict(base, **(cert or {}), emp_T=des.T, emp_m=des.m, leak_exact=float(leak_exact),
               times=list(map(float, des.times)), phases=list(map(float, des.phases)))
    if tr:
        dt, k, inf, p = tr
        cnot = filter_counts(Hs, k, 2, anc_rot_per_pulse=1)["cnot"]
        out.update(dt=dt, steps=int(np.sum(k)), inf_trotter=inf, p_succ=p, cnot_filter=cnot,
                   cnot_per_success=cnot / p)
    return out


def run_sbc():
    cases = []
    for fam, N, J2 in families():
        if fam == "gapscan":
            continue
        for eps in EPSS:
            cases += [(fam, N, J2, ("synthetic", gm), eps, "dyadic") for gm in (0.03, 0.3, 0.9)]
            cases += [(fam, N, J2, chi, eps, "dyadic") for chi in (1, 2)]
            if N in (8, 12, 16) and eps in (1e-2, 1e-4):
                cases += [(fam, N, J2, ("synthetic", gm), eps, "optimized") for gm in (0.03, 0.3, 0.9)]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_sbc_case, cases):
            if not r.get("skip"):
                print(f"[sbc] {r['family']} N={r['N']} {r['trial']} eps={r['eps']:g} {r['schedule']} "
                      f"T={r.get('emp_T', float('nan')):.1f} leak={r.get('leak_exact', float('nan')):.1e} "
                      f"dt={r.get('dt', float('nan')):.2f} steps={r.get('steps')} cert_steps={r.get('cert_steps')} "
                      f"ratio={r.get('cert_cnot_per_success', float('nan')) / r.get('cnot_per_success', float('nan')):.1f}",
                      flush=True)
            out.append(r)
    dump("emp_sbc.json", out)


# ----------------------------------------------------------------------------- SBC on the regime systems
def _sbc_regime_case(args):
    """SBC (dyadic, or optimized) designed for the full target with the pipeline's gap, one Strang step
    per pulse; simulated where a trial vector exists. With grid=True only the cost model on EPS_GRID."""
    kind, size, J2, chi, eps_list, sched, simulate, dt_model = args
    import run_extended as RX
    from sbpipe.simulate import filter_trotter, state_metrics
    from sbpipe.resources import filter_counts
    rec, est, ans, g = RX.get_case(kind, size, J2, chi)
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    if kind == "chain" and ans.statevector is None and N <= 16:
        p = os.path.join(CACHEE, _key(N, J2) + ".npz")
        if os.path.exists(p):
            z = np.load(p)
            if f"trial{chi}" in z:
                ans.statevector, g = z[f"trial{chi}"], z["ground"]
    lo, hi = sbpipe.weyl_bounds(H)
    win = est.window(hi)
    Hs = H.shifted_scaled(win["shift"], win["W"])
    per_pulse = filter_counts(Hs, [1, 1], 2)["cnot"] - filter_counts(Hs, [1], 2)["cnot"]
    pts = []
    for eps in eps_list:
        if 1 - ans.gamma <= eps:
            pts.append(dict(eps=float(eps), m=0, cnot=float(ans.cost["cnot"])))
            continue
        des = D.design_filter(sched, win["delta"], win["w0"], eps, ans.gamma)
        if des is None:
            continue
        dt = dt_model[0] * eps ** dt_model[1]          # empirical Strang step law (scaled units)
        k = np.maximum(1, np.ceil(np.abs(des.times) / dt - 1e-12)).astype(int)
        cnot_f = filter_counts(Hs, k, 2)["cnot"]
        pt = dict(eps=float(eps), m=des.m, T=des.T, steps_model=int(k.sum()),
                  cnot=(cnot_f + ans.cost["cnot"]) / ans.gamma)
        if simulate and ans.statevector is not None and g is not None:
            tr = _sbc_trotter(H, Hs, ans.statevector, g, des, eps)
            v = filter_trotter(Hs, ans.statevector, des.times, des.phases, k, 2)
            m = state_metrics(v, g)
            pt.update(model_inf=m["infidelity"], model_p=float(np.vdot(v, v).real))
            if tr:
                dt, ks, inf, p = tr
                pt.update(sim_dt=dt, sim_steps=int(ks.sum()), sim_inf=inf, sim_p=p,
                          cnot_sim=(filter_counts(Hs, ks, 2)["cnot"] + ans.cost["cnot"]) / p)
        pts.append(pt)
    return dict(kind=kind, size=size, N=N, J2=J2, chi=chi, schedule=sched, gamma=ans.gamma,
                trial_cnot=ans.cost["cnot"], points=pts)


def run_sbc_regime():
    import run_extended as RX
    systems = [("chain", N, 0.2411, (1, 2, 4, 8, 16, 32)) for N in (16, 32, 48)]
    systems += [("cyl", Lx, J2, tuple(RX.CHIS_2D)) for Lx in RX.LX for J2 in RX.J2S_2D]
    # smallest step that sufficed for the dyadic filter on the exact chains (N = 8-16), per target:
    # 5.8, 3.3, 1.6 at eps = 1e-2, 1e-3, 1e-4  ->  dt = 20 eps^0.27 (a lower envelope)
    dt_model = (20.0, 0.27)
    cases = [(k, s, J2, chi, list(EPS_GRID), "dyadic", False, dt_model) for k, s, J2, chis in systems for chi in chis]
    small = [(k, s, J2, chis) for k, s, J2, chis in systems if (k == "chain" and s == 16) or (k == "cyl" and s <= 4)]
    cases += [(k, s, J2, chi, [1e-2, 1e-3, 1e-4], sc, True, dt_model) for k, s, J2, chis in small for chi in chis
              for sc in ("dyadic", "optimized")]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_sbc_regime_case, cases):
            for p in r["points"]:
                if "model_inf" in p:
                    print(f"[sbc_regime] {r['kind']} N={r['N']} J2={r['J2']} chi={r['chi']} {r['schedule']} "
                          f"eps={p['eps']:g} m={p['m']} model_inf={p['model_inf']:.1e} "
                          f"ok={p['model_inf'] <= p['eps']} sim_dt={p.get('sim_dt')}", flush=True)
            out.append(r)
    dump("emp_sbc_regime.json", out)


# ----------------------------------------------------------------------------- actual filter behaviour, no models
ACT_SYSTEMS = [("chain", 16, 0.2411), ("cyl", 3, 0.0), ("cyl", 3, 0.3), ("cyl", 4, 0.0), ("cyl", 4, 0.3)]
ACT_CHIS = [2, 4, 8, 16, 32, 64]
EPS_ACT = [1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]


def _act_key(kind, size, J2):
    return f"act_{kind}{size}_J2{J2:g}"


def _act_spectrum(args):
    """Exact eigen-decomposition (S^z blocks) of every MPS trial of one system; cached."""
    kind, size, J2 = args
    import run_extended as RX
    from sbpipe.frontends.dmrg_quimb import pauli_to_mpo, trial_mps
    p = os.path.join(CACHEE, _act_key(kind, size, J2) + ".npz")
    if os.path.exists(p):
        return args, 0.0
    t0 = time.time()
    H = RX.model(kind, size, J2)
    N = H.num_qubits
    trials, costs = {}, {}
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    for chi in ACT_CHIS:
        try:
            rec, est, ans, g = RX.get_case(kind, size, J2, chi)
        except KeyError:
            ans = None
        if ans is None or ans.statevector is None:  # chain: build the same variational MPS here
            _, psi = trial_mps(pauli_to_mpo(H), chi)
            v = np.asarray(psi.to_dense()).reshape(-1)
        else:
            v = np.asarray(ans.statevector).reshape(-1)
        trials[chi] = v / np.linalg.norm(v)
        costs[chi] = mps_circuit_cost(N, chi)["cnot"]
    Hsp = H.to_sparse().tocsr()
    idx = np.arange(2 ** N)
    pop = np.array([bin(i).count("1") for i in idx])
    E, W, low = [], {c: [] for c in ACT_CHIS}, []
    for n_up in range(N + 1):
        b = idx[pop == n_up]
        Hb = Hsp[b][:, b].toarray()
        if np.iscomplexobj(Hb) and np.abs(Hb.imag).max() < 1e-12:
            Hb = np.ascontiguousarray(Hb.real)
        e, U = scipy.linalg.eigh(Hb, overwrite_a=True, check_finite=False, driver="evr")
        del Hb
        E.append(e)
        for c in ACT_CHIS:
            W[c].append(np.abs(U.T.conj() @ trials[c][b]) ** 2)
        low.append((e[0], b, U[:, 0].copy()))
        del U
        gc.collect()
    E = np.concatenate(E)
    o = np.argsort(E)
    out = dict(E=E[o])
    for c in ACT_CHIS:
        out[f"w{c}"] = np.concatenate(W[c])[o]
        out[f"trial{c}"] = trials[c]
        out[f"cost{c}"] = costs[c]
    e0, b, u = min(low, key=lambda t: t[0])
    gvec = np.zeros(2 ** N, complex)
    gvec[b] = u
    out["ground"] = gvec
    np.savez(p, **out)
    return args, time.time() - t0


def _dyadic_candidates(Es, w, delta, eps, n=2):
    """Shortest dyadic ladders (K, r) whose exact filter reaches eps on this trial (exact spectrum)."""
    tau = math.pi / (1 + delta)
    out = []
    for K in range(1, 24):
        blk = np.ones_like(Es)
        for k in range(K):
            blk = blk * np.cos(2 ** k * Es * tau)
        for r in range(1, 400):
            inf, _ = _infid(w, blk ** r)
            if inf <= eps:
                out.append((r * tau * (2 ** K - 1), K, r))
                break
    out.sort()
    return out[:n]


def _act_case(args):
    kind, size, J2, chi, eps = args
    import run_extended as RX
    from sbpipe.resources import filter_counts
    from sbpipe.simulate import qetu_trotter, state_metrics
    z = np.load(os.path.join(CACHEE, _act_key(kind, size, J2) + ".npz"))
    E, w, psi, g = z["E"], z[f"w{chi}"], z[f"trial{chi}"], z["ground"]
    trial_cnot = float(z[f"cost{chi}"])
    H = RX.model(kind, size, J2)
    lo, hi = sbpipe.weyl_bounds(H)
    Wr = hi - E[0]
    Es = (E - E[0]) / Wr
    gap = float(E[np.argmax(E > E[0] + 1e-8)] - E[0])
    delta = gap / Wr
    gamma = float(w[0])
    Hs = H.shifted_scaled(E[0], Wr)
    out = dict(kind=kind, size=size, N=H.num_qubits, J2=J2, chi=chi, eps=eps, gamma=gamma, gap=gap,
               delta=delta, trial_cnot=trial_cnot, mps_alone=bool(1 - gamma <= eps))
    if out["mps_alone"]:
        return out
    res = {}
    # dyadic: the shortest exact-passing ladders, each with the coarsest Trotter step that still passes
    best = None
    for T, K, r in _dyadic_candidates(Es, w, delta, eps):
        tau = math.pi / (1 + delta)
        times = np.tile(tau * 2.0 ** np.arange(K), r)
        des = D.FilterDesign("dyadic", times, np.zeros(len(times)), delta, 0.0, 0.0, 1.0, {})
        tr = _sbc_trotter(H, Hs, psi, g, des, eps)
        if tr:
            dt, k, inf, p = tr
            c = (filter_counts(Hs, k, 2)["cnot"] + trial_cnot) / p
            if best is None or c < best["cnot"]:
                best = dict(cnot=c, T=T, K=K, r=r, steps=int(k.sum()), inf=inf, p=p, dt=dt)
    res["dyadic"] = best
    # optimized (Stetcu-Baroni-Carlson recipe), designed for the full target with the true gap
    if eps in (1e-2, 1e-3, 1e-4):
        des = D.design_filter("optimized", delta, 0.0, eps, gamma)
        if des is not None:
            tr = _sbc_trotter(H, Hs, psi, g, des, eps)
            if tr:
                dt, k, inf, p = tr
                res["optimized"] = dict(cnot=(filter_counts(Hs, k, 2)["cnot"] + trial_cnot) / p, T=des.T,
                                        m=des.m, steps=int(k.sum()), inf=inf, p=p, dt=dt)
    # QETU: smallest exact-passing degree, then the fewest steps per query (or a slightly higher degree)
    d, _ = _min_qetu(Es, w, delta, eps)
    if d:
        x_s = math.cos(math.pi / 2 * delta)
        q = None
        for dd in range(d, int(1.5 * d) + 3, 2):
            for kk in (1, 2, 3, 4):
                v = qetu_trotter(Hs, psi, dd, math.pi / 2, x_s, [kk] * dd)
                inf = state_metrics(v, g)["infidelity"]
                if inf <= eps:
                    p = float(np.vdot(v, v).real)
                    c = (filter_counts(Hs, [kk] * dd, 2, anc_rot_per_pulse=2)["cnot"] + trial_cnot) / p
                    if q is None or c < q["cnot"]:
                        q = dict(cnot=c, d=dd, k=kk, inf=inf, p=p)
                    break
            if q is not None and dd * 1 > q["cnot"]:
                break
            if q is not None:
                break
        res["qetu"] = q
    return dict(out, **res)


def run_actual():
    with ProcessPoolExecutor(2) as ex:                  # N = 16 diagonalisations need ~4 GB each
        for a, dt in ex.map(_act_spectrum, ACT_SYSTEMS):
            print(f"[actual] spectrum {a} {dt:.0f}s", flush=True)
    cases = [(k, s, J2, chi, eps) for k, s, J2 in ACT_SYSTEMS for chi in ACT_CHIS for eps in EPS_ACT]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_act_case, cases):
            msg = "MPS alone" if r["mps_alone"] else " ".join(
                f"{s}={r[s]['cnot']:.2e}" for s in ("dyadic", "optimized", "qetu") if r.get(s))
            print(f"[actual] {r['kind']}{r['size']} J2={r['J2']} chi={r['chi']} eps={r['eps']:g} "
                  f"1-gamma={1 - r['gamma']:.1e} mps={r['trial_cnot']:.2e} {msg}", flush=True)
            out.append(r)
    dump("emp_actual.json", out)


ACT1D_SYSTEMS = [("chain", N, J2) for J2 in (0.0, 0.2411, 0.4) for N in (12, 14, 16)
                 if (N, J2) != (16, 0.2411)]


def run_actual1d():
    """The measured crossover of run_actual on more chains (N = 12-16, J2 = 0, 0.2411, 0.4); one json line per
    case in results/emp_actual1d.jsonl, resumable."""
    import json
    with ProcessPoolExecutor(2) as ex:
        for a, dt in ex.map(_act_spectrum, ACT1D_SYSTEMS):
            print(f"[actual1d] spectrum {a} {dt:.0f}s", flush=True)
    path = os.path.join(RES, "emp_actual1d.jsonl")
    done = set()
    if os.path.exists(path):
        done = {(r["size"], r["J2"], r["chi"], r["eps"]) for r in map(json.loads, open(path))}
    cases = [(k, s, J2, chi, eps) for k, s, J2 in ACT1D_SYSTEMS for chi in ACT_CHIS for eps in EPS_ACT
             if (s, J2, chi, eps) not in done]
    with ProcessPoolExecutor(JOBS) as ex, open(path, "a") as fh:
        for r in ex.map(_act_case, cases):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()
            msg = "MPS alone" if r["mps_alone"] else " ".join(
                f"{s}={r[s]['cnot']:.2e}" for s in ("dyadic", "optimized", "qetu") if r.get(s))
            print(f"[actual1d] {r['kind']}{r['size']} J2={r['J2']} chi={r['chi']} eps={r['eps']:g} "
                  f"1-gamma={1 - r['gamma']:.1e} mps={r['trial_cnot']:.2e} {msg}", flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "match,spectra,filters,trotter"
    steps = dict(match=run_match, gapscan=run_gapscan, spectra=run_spectra, filters=run_filters, trotter=run_trotter, emp_regime=run_emp_regime, frontier=run_frontier, sbc=run_sbc, sbc_regime=run_sbc_regime, actual=run_actual, actual1d=run_actual1d)
    for k in what.split(","):
        t = time.time()
        steps[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)
