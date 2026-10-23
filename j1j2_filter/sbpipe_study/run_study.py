"""Application study: DMRG inputs -> sbpipe -> certified Stetcu-Baroni filter for the open J1-J2 chain.

    python sbpipe_study/run_study.py inputs      # DMRG front end (+ ED cross-check for N <= 14), cached
    python sbpipe_study/run_study.py validate    # A: simulate planned circuits, check the certificate
    python sbpipe_study/run_study.py scaling     # B: resources vs N (DMRG only, no simulation)
    python sbpipe_study/run_study.py sensitivity # C: wrong / conservative energy estimates
    python sbpipe_study/run_study.py ansatz      # D: trial bond dimension vs total cost
    python sbpipe_study/run_study.py oracle      # E: smallest step count that actually meets eps
    python sbpipe_study/run_study.py all

Every number in the paper is read from results/*.json written here.
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace

import numpy as np

import sbpipe
from sbpipe import Ansatz, Config, EnergyEstimates, j1j2_chain, refine, validate
from sbpipe.simulate import filter_trotter, state_metrics

RES = os.path.join(HERE, "results")
CACHE = os.path.join(RES, "cache")
os.makedirs(CACHE, exist_ok=True)
JOBS = int(os.environ.get("JOBS", "4"))


def _install_design_cache():
    """Filter designs depend only on (schedule, delta, w0, eps_leak, gamma); the optimized ones take
    minutes at N = 48, and both Trotter orders and several studies reuse them. Cache them on disk."""
    import hashlib
    import pickle

    import sbpipe.design as SD
    if getattr(SD.design_filter, "_cached", False):
        return
    orig = SD.design_filter
    d = os.path.join(CACHE, "designs")
    os.makedirs(d, exist_ok=True)

    def cached(schedule, delta, w0, eps_leak, gamma, **kw):
        key = repr((schedule, float(delta), float(w0), float(eps_leak), float(gamma), sorted(kw.items())))
        p = os.path.join(d, hashlib.sha1(key.encode()).hexdigest() + ".pkl")
        if os.path.exists(p):
            with open(p, "rb") as fh:
                return pickle.load(fh)
        out = orig(schedule, delta, w0, eps_leak, gamma, **kw)
        tmp = f"{p}.{os.getpid()}.tmp"
        with open(tmp, "wb") as fh:
            pickle.dump(out, fh)
        os.replace(tmp, p)
        return out

    cached._cached = True
    SD.design_filter = cached


_install_design_cache()

J2S = [0.0, 0.2411, 0.4]
N_VAL = [8, 10, 12, 14]
N_SCALE = [8, 12, 16, 24, 32, 48]
N_INPUTS = [8, 10, 12, 14, 16, 20, 24, 32, 40, 48]      # DMRG inputs computed (and tabulated)
CHIS = [1, 2, 4, 8]
ED_MAX_N = 14


def _f(x):
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def dump(name, obj):
    with open(os.path.join(RES, name), "w") as fh:
        json.dump(obj, fh, indent=1, default=_f)


def load(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


# ----------------------------------------------------------------------------- inputs (cached)
def _key(N, J2):
    return f"N{N}_J2{J2:g}"


def compute_inputs(args):
    N, J2 = args
    from scipy.sparse.linalg import eigsh

    from sbpipe.frontends.dmrg_quimb import dmrg_inputs
    t0 = time.time()
    H = j1j2_chain(N, 1.0, J2)
    bd = (16, 32, 64, 128, 256) if N <= 24 else (16, 32, 64, 128)
    est, ans_list, ex = dmrg_inputs(H, chi_trial=CHIS, bond_dims=bd)
    rec = dict(N=N, J2=J2, est=asdict(est), dmrg_time=time.time() - t0,
               ref_history=ex["ref_history"], E1_history=ex["E1_history"],
               ansatz={a.meta["chi"]: dict(gamma=a.gamma, cost=a.cost, meta=a.meta) for a in ans_list})
    arrays = {}
    if N <= ED_MAX_N:
        Hsp = H.to_sparse()
        w, v = eigsh(Hsp, k=8, which="SA", tol=1e-13, v0=np.ones(2 ** N))
        o = np.argsort(w)
        w, v = w[o], v[:, o]
        g = v[:, 0]
        # lowest level of any symmetry above E0 (first gap) and its S^z content
        E1_ed = float(w[np.argmax(w > w[0] + 1e-8)])
        lo, hi = sbpipe.weyl_bounds(H)
        top = float(eigsh(Hsp, k=1, which="LA", return_eigenvectors=False, tol=1e-12)[0])
        rec["ed"] = dict(E0=float(w[0]), E1=E1_ed, levels=w.tolist(), E_top=top, weyl_hi=hi,
                         gamma={a.meta["chi"]: float(abs(np.vdot(g, a.statevector)) ** 2) for a in ans_list})
        arrays["ground"] = g
        for a in ans_list:
            arrays[f"trial_chi{a.meta['chi']}"] = a.statevector
        np.savez(os.path.join(CACHE, _key(N, J2) + ".npz"), **arrays)
    with open(os.path.join(CACHE, _key(N, J2) + ".json"), "w") as fh:
        json.dump(rec, fh, default=_f)
    return _key(N, J2), time.time() - t0


def get_inputs(N, J2, chi=2):
    with open(os.path.join(CACHE, _key(N, J2) + ".json")) as fh:
        rec = json.load(fh)
    est = EnergyEstimates(**rec["est"])
    a = rec["ansatz"][str(chi)]
    sv, g = None, None
    p = os.path.join(CACHE, _key(N, J2) + ".npz")
    if os.path.exists(p):
        z = np.load(p)
        sv, g = z[f"trial_chi{chi}"], z["ground"]
    ans = Ansatz(gamma=a["gamma"], statevector=sv, cost=a["cost"], name=f"chi={chi}", meta=a["meta"])
    return rec, est, ans, g


def run_inputs():
    todo = [(N, J2) for J2 in J2S for N in sorted(N_INPUTS, reverse=True)]
    todo = [t for t in todo if not os.path.exists(os.path.join(CACHE, _key(*t) + ".json"))]
    with ProcessPoolExecutor(JOBS) as ex:
        for k, dt in ex.map(compute_inputs, todo):
            print(f"[inputs] {k} {dt:.0f}s", flush=True)
    rows = []
    for J2 in J2S:
        for N in sorted(N_INPUTS):
            rec, est, ans, _ = get_inputs(N, J2)
            r = dict(N=N, J2=J2, E0=est.E0, E1=est.E1, gap=est.gap, E0_err=est.E0_err,
                     E1_err=est.E1_err, dmrg_time=rec["dmrg_time"],
                     gamma={c: rec["ansatz"][c]["gamma"] for c in rec["ansatz"]},
                     max_bond_ref=rec["ref_history"][-1][2])
            if "ed" in rec:
                ed = rec["ed"]
                r.update(ed_E0=ed["E0"], ed_E1=ed["E1"], dE0=est.E0 - ed["E0"], dE1=est.E1 - ed["E1"],
                         ed_gamma=ed["gamma"], ed_Etop=ed["E_top"], weyl_hi=ed["weyl_hi"])
            rows.append(r)
    dump("inputs.json", rows)


# ----------------------------------------------------------------------------- A: validation
def _val_case(args):
    N, J2, eps, sched, order = args
    rec, est, ans, g = get_inputs(N, J2)
    H = j1j2_chain(N, 1.0, J2)
    t0 = time.time()
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=order))
    out = dict(N=N, J2=J2, eps=eps, schedule=sched, order=order, plan=plan.summary(),
               plan_time=time.time() - t0, candidates=plan.candidates,
               times=plan.design.times, phases=plan.design.phases, steps=plan.steps)
    rot = plan.resources["filter_rotations"] * (2 ** N)
    if rot <= 1e10:
        t1 = time.time()
        out["sim"] = validate(plan, H, ans, g)
        out["sim_time"] = time.time() - t1
    return out


def run_validate():
    cases = [(N, J2, eps, s, o) for N in N_VAL for J2 in J2S for eps in (1e-2, 1e-3)
             for s in ("uniform", "dyadic", "optimized") for o in (1, 2)]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_val_case, cases):
            s = r.get("sim", {})
            print(f"[A] N={r['N']} J2={r['J2']} eps={r['eps']:g} {r['schedule']:9s} o{r['order']} "
                  f"steps={r['plan']['res_steps']} bound={r['plan']['infidelity_bound']:.2e} "
                  f"meas={s.get('infidelity_trotter', float('nan')):.2e} ok={s.get('bound_ok')}",
                  flush=True)
            out.append(r)
    dump("validate.json", out)


# ----------------------------------------------------------------------------- B: scaling
def _scale_case(args):
    N, J2, eps, sched, order, chi = args
    rec, est, ans, _ = get_inputs(N, J2, chi)
    ans.statevector = None
    base = dict(N=N, J2=J2, eps=eps, schedule=sched, order=order, chi=chi, gamma=ans.gamma)
    try:
        plan = refine(j1j2_chain(N, 1.0, J2), est, ans, Config(eps=eps, schedule=sched, order=order))
    except RuntimeError as e:
        return dict(base, plan=None, error=str(e))
    return dict(base, plan=plan.summary())


def run_scaling():
    cases = [(N, J2, eps, s, o, chi) for N in N_SCALE for J2 in J2S for eps in (1e-2, 1e-3)
             for s in ("dyadic", "optimized") for o in (1, 2) for chi in (2,)]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_scale_case, cases))
    dump("scaling.json", out)


# ----------------------------------------------------------------------------- C: sensitivity
def _sens_case(args):
    kind, f, N, J2, eps = args
    rec, est, ans, g = get_inputs(N, J2)
    H = j1j2_chain(N, 1.0, J2)
    gap = est.gap
    if kind == "gap":          # assumed E1 = E0 + f * gap (f > 1: over-estimated gap)
        e = replace(est, E1=est.E0 + f * gap, E1_err=0.0)
    elif kind == "E0err":      # error bar on E0 as a fraction of the gap
        e = replace(est, E0_err=f * gap)
    else:                      # gamma used in the design = f * true gamma
        e = est
    a = replace(ans, gamma=min(1.0, f * ans.gamma)) if kind == "gamma" else ans
    plan = refine(H, e, a, Config(eps=eps, schedule="optimized", order=2))
    sim = validate(plan, H, ans, g)            # always simulate the TRUE trial state
    return dict(kind=kind, f=f, N=N, J2=J2, eps=eps, plan=plan.summary(), sim=sim)


def run_sensitivity():
    N, J2, eps = 12, 0.2411, 1e-2
    cases = [("gap", f, N, J2, eps) for f in (0.5, 0.7, 0.85, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0)]
    cases += [("E0err", f, N, J2, eps) for f in (1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1)]
    cases += [("gamma", f, N, J2, eps) for f in (0.25, 0.5, 0.8, 1.0, 1.02)]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_sens_case, cases))
    dump("sensitivity.json", out)


# ----------------------------------------------------------------------------- D: ansatz bond dimension
def _ansatz_case(args):
    N, J2, eps, chi = args
    rec, est, ans, _ = get_inputs(N, J2, chi)
    ans.statevector = None
    try:
        plan = refine(j1j2_chain(N, 1.0, J2), est, ans, Config(eps=eps, schedule="auto", order=2))
    except (RuntimeError, ValueError) as e:      # e.g. gamma ~ 1e-6 for a product state
        return dict(N=N, J2=J2, eps=eps, chi=chi, gamma=ans.gamma, trial_cost=ans.cost,
                    plan=None, error=str(e))
    return dict(N=N, J2=J2, eps=eps, chi=chi, gamma=ans.gamma, trial_cost=ans.cost,
                plan=plan.summary())


def run_ansatz():
    cases = [(N, J2, eps, chi) for N in (16, 32, 48) for J2 in J2S for eps in (1e-2,)
             for chi in CHIS]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_ansatz_case, cases))
    dump("ansatz.json", out)


# ----------------------------------------------------------------------------- E: oracle step count
def _oracle_case(args):
    N, J2, eps, sched = args
    rec, est, ans, g = get_inputs(N, J2)
    H = j1j2_chain(N, 1.0, J2)
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    d = plan.design
    Hsp = H.to_sparse()

    def meas(s):
        k = np.maximum(1, np.ceil(s * plan.steps)).astype(int)
        v = filter_trotter(plan.H_scaled, ans.statevector, d.times, d.phases, k, 2)
        return state_metrics(v, g)["infidelity"], int(k.sum())

    lo, hi = 0.0, 1.0
    inf_hi, n_hi = meas(1.0)
    assert inf_hi <= eps
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        inf, n = meas(mid)
        if inf <= eps:
            hi, n_hi = mid, n
        else:
            lo = mid
        if n_hi <= len(d.times):
            break
    return dict(N=N, J2=J2, eps=eps, schedule=sched, n_cert=int(plan.steps.sum()), n_oracle=n_hi,
                m=d.m, ratio=plan.steps.sum() / n_hi)


def run_oracle():
    cases = [(N, J2, eps, s) for N in (8, 10, 12) for J2 in J2S for eps in (1e-2, 1e-3)
             for s in ("dyadic", "optimized")]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_oracle_case, cases))
    dump("oracle.json", out)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = dict(inputs=run_inputs, validate=run_validate, scaling=run_scaling,
                 sensitivity=run_sensitivity, ansatz=run_ansatz, oracle=run_oracle)
    for k in (steps if what == "all" else what.split(",")):
        t = time.time()
        steps[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)
