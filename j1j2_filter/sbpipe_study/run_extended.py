"""Extended study for the second referee round (items 2 and 3 of the review):

  * a QETU baseline (Chebyshev polynomial filter, Dong-Lin-Tong) costed with the same Trotter,
    synthesis and gate model as the Stetcu-Baroni-Carlson schedules;
  * tighter targets (eps down to 1e-4) and larger trial bond dimension on the chain;
  * a two-dimensional model (J1-J2 on 4-leg square cylinders), where MPS circuits scale worse.

    python sbpipe_study/run_extended.py inputs2d     # DMRG (+ ED for N <= 16) on Ly = 4 cylinders, cached
    python sbpipe_study/run_extended.py chi_chain    # chain trial states with chi = 16, 32, cached
    python sbpipe_study/run_extended.py qetu_val     # simulate QETU plans, chain and cylinders (N <= 16)
    python sbpipe_study/run_extended.py qetu_scale   # QETU vs SBC schedules, chain, chi = 2 (Table IV setting)
    python sbpipe_study/run_extended.py regime       # cheapest route to eps vs chi: MPS alone, SBC, QETU

Results go to results/ext_*.json; run_study.py's results are not touched.
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict

import numpy as np

import sbpipe
from sbpipe import Ansatz, Config, EnergyEstimates, j1j2_chain, refine, validate
from sbpipe.models import j1j2_square

import run_study as RS        # shares the design cache, the chain inputs and the json helpers

RES, CACHE = RS.RES, RS.CACHE
CACHE2D = os.path.join(RES, "cache2d")
CACHECHI = os.path.join(RES, "cache_chi")
for d in (CACHE2D, CACHECHI):
    os.makedirs(d, exist_ok=True)
JOBS = RS.JOBS
dump, load, _f = RS.dump, RS.load, RS._f

LY = 4
LX = [3, 4, 6, 8, 12]                  # N = 12, 16, 24, 32, 48
J2S_2D = [0.0, 0.3]                    # ED: the lowest excitation is the triplet (not so at J2 = 0.5)
CHIS_2D = [2, 4, 8, 16, 32, 64]
CHIS_CHAIN_EXTRA = [16, 32]
ED_MAX_N_2D = 16


def model(kind, size, J2):
    return j1j2_chain(size, 1.0, J2) if kind == "chain" else j1j2_square(size, LY, 1.0, J2)


def _k2(Lx, J2):
    return f"Lx{Lx}_Ly{LY}_J2{J2:g}"


# ----------------------------------------------------------------------------- inputs: cylinders
def compute_inputs2d(args):
    Lx, J2 = args
    from scipy.sparse.linalg import eigsh

    from sbpipe.frontends.dmrg_quimb import dmrg_inputs
    t0 = time.time()
    H = j1j2_square(Lx, LY, 1.0, J2)
    N = H.num_qubits
    bd = (16, 32, 64, 128, 256)
    est, ans_list, ex = dmrg_inputs(H, chi_trial=CHIS_2D, bond_dims=bd)
    rec = dict(Lx=Lx, Ly=LY, N=N, J2=J2, est=asdict(est), dmrg_time=time.time() - t0,
               ref_history=ex["ref_history"], E1_history=ex["E1_history"],
               ansatz={a.meta["chi"]: dict(gamma=a.gamma, cost=a.cost, meta=a.meta) for a in ans_list})
    if N <= ED_MAX_N_2D:
        Hsp = H.to_sparse()
        w, v = eigsh(Hsp, k=6, which="SA", tol=1e-13, v0=np.ones(2 ** N))
        o = np.argsort(w)
        w, v = w[o], v[:, o]
        g = v[:, 0]
        top = float(eigsh(Hsp, k=1, which="LA", return_eigenvectors=False, tol=1e-12)[0])
        rec["ed"] = dict(E0=float(w[0]), E1=float(w[np.argmax(w > w[0] + 1e-8)]), levels=w.tolist(),
                         E_top=top, weyl_hi=sbpipe.weyl_bounds(H)[1],
                         gamma={a.meta["chi"]: float(abs(np.vdot(g, a.statevector)) ** 2) for a in ans_list})
        arrays = {"ground": g}
        for a in ans_list:
            arrays[f"trial_chi{a.meta['chi']}"] = a.statevector
        np.savez(os.path.join(CACHE2D, _k2(Lx, J2) + ".npz"), **arrays)
    with open(os.path.join(CACHE2D, _k2(Lx, J2) + ".json"), "w") as fh:
        json.dump(rec, fh, default=_f)
    return _k2(Lx, J2), time.time() - t0


def run_inputs2d():
    todo = [(Lx, J2) for Lx in sorted(LX) for J2 in J2S_2D]
    todo = [t for t in todo if not os.path.exists(os.path.join(CACHE2D, _k2(*t) + ".json"))]
    with ProcessPoolExecutor(JOBS) as ex:
        for k, dt in ex.map(compute_inputs2d, todo):
            print(f"[inputs2d] {k} {dt:.0f}s", flush=True)
    rows = []
    for J2 in J2S_2D:
        for Lx in LX:
            rec = json.load(open(os.path.join(CACHE2D, _k2(Lx, J2) + ".json")))
            est = rec["est"]
            r = dict(Lx=Lx, Ly=LY, N=rec["N"], J2=J2, E0=est["E0"], E1=est["E1"],
                     gap=est["E1"] - est["E0"], E0_err=est["E0_err"], E1_err=est["E1_err"],
                     dmrg_time=rec["dmrg_time"], gamma={c: v["gamma"] for c, v in rec["ansatz"].items()},
                     ref_history=rec["ref_history"])
            if "ed" in rec:
                ed = rec["ed"]
                r.update(ed_E0=ed["E0"], ed_E1=ed["E1"], dE0=est["E0"] - ed["E0"], dE1=est["E1"] - ed["E1"],
                         ed_gamma=ed["gamma"], ed_Etop=ed["E_top"], weyl_hi=ed["weyl_hi"])
            rows.append(r)
    dump("ext_inputs2d.json", rows)


# ----------------------------------------------------------------------------- inputs: chain, chi 16 and 32
def compute_chi_chain(args):
    N, J2 = args
    from sbpipe.frontends.dmrg_quimb import ansatz_from_dmrg, ground_state, pauli_to_mpo
    t0 = time.time()
    H = j1j2_chain(N, 1.0, J2)
    mpo = pauli_to_mpo(H)
    bd = (16, 32, 64, 128, 256) if N <= 24 else (16, 32, 64, 128)
    g = ground_state(mpo, bond_dims=bd)
    out = dict(N=N, J2=J2, E0_rerun=g["E0"], ansatz={})
    for chi in [8] + CHIS_CHAIN_EXTRA:          # chi = 8 again, as a check against the cached run
        a = ansatz_from_dmrg(mpo, g["psi"], chi, N, dense_max_n=0)
        out["ansatz"][chi] = dict(gamma=a.gamma, cost=a.cost, meta=a.meta)
    out["time"] = time.time() - t0
    with open(os.path.join(CACHECHI, RS._key(N, J2) + ".json"), "w") as fh:
        json.dump(out, fh, default=_f)
    return RS._key(N, J2), out["time"]


def run_chi_chain():
    todo = [(N, J2) for N in (48, 32, 16) for J2 in RS.J2S]
    todo = [t for t in todo if not os.path.exists(os.path.join(CACHECHI, RS._key(*t) + ".json"))]
    with ProcessPoolExecutor(JOBS) as ex:
        for k, dt in ex.map(compute_chi_chain, todo):
            print(f"[chi_chain] {k} {dt:.0f}s", flush=True)


# ----------------------------------------------------------------------------- inputs: chain, N <= 14 with statevectors
CACHESMALL = os.path.join(RES, "cache_small")


def _chain_small(args):
    RS.CACHE = CACHESMALL                 # same front end as run_study.inputs, separate cache
    return RS.compute_inputs(args)


def run_chain_small():
    """Statevectors (trial, ED ground state) for the chain at N <= 14; the .npz files are not in git."""
    os.makedirs(CACHESMALL, exist_ok=True)
    todo = [(N, J2) for N in RS.N_VAL for J2 in RS.J2S]
    todo = [t for t in todo if not os.path.exists(os.path.join(CACHESMALL, RS._key(*t) + ".npz"))]
    with ProcessPoolExecutor(JOBS) as ex:
        for k, dt in ex.map(_chain_small, todo):
            print(f"[chain_small] {k} {dt:.0f}s", flush=True)


def get_case(kind, size, J2, chi):
    """(est, ansatz, ground-or-None) for a chain of N = size or a cylinder of Lx = size."""
    if kind == "chain":
        base = CACHESMALL if os.path.exists(os.path.join(CACHESMALL, RS._key(size, J2) + ".npz")) else CACHE
        key = RS._key(size, J2)
        rec = json.load(open(os.path.join(base, key + ".json")))
        if str(chi) not in rec["ansatz"]:
            extra = json.load(open(os.path.join(CACHECHI, key + ".json")))
            rec["ansatz"][str(chi)] = extra["ansatz"][str(chi)]
    else:
        base, key = CACHE2D, _k2(size, J2)
        rec = json.load(open(os.path.join(base, key + ".json")))
    est = EnergyEstimates(**rec["est"])
    a = rec["ansatz"][str(chi)]
    sv = g = None
    p = os.path.join(base, key + ".npz")
    if os.path.exists(p):
        z = np.load(p)
        if f"trial_chi{chi}" in z:
            sv, g = z[f"trial_chi{chi}"], z["ground"]
    ans = Ansatz(gamma=a["gamma"], statevector=sv, cost=a["cost"], name=f"chi={chi}", meta=a["meta"])
    return rec, est, ans, g


# ----------------------------------------------------------------------------- QETU: validation
def _qv_case(args):
    kind, size, J2, eps, sched, chi = args
    rec, est, ans, g = get_case(kind, size, J2, chi)
    H = model(kind, size, J2)
    t0 = time.time()
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    out = dict(kind=kind, size=size, N=H.num_qubits, J2=J2, eps=eps, schedule=sched, order=2, chi=chi,
               plan=plan.summary(), plan_time=time.time() - t0, steps=plan.steps)
    if ans.statevector is not None and plan.resources["filter_rotations"] * 2 ** H.num_qubits <= 3e10:
        t1 = time.time()
        out["sim"] = validate(plan, H, ans, g)
        out["sim_time"] = time.time() - t1
    return out


def run_qetu_val():
    cases = [("chain", N, J2, eps, "qetu", 2) for N in RS.N_VAL for J2 in RS.J2S for eps in (1e-2, 1e-3)]
    cases += [("cyl", Lx, J2, eps, s, 2) for Lx in (3, 4) for J2 in J2S_2D for eps in (1e-2, 1e-3)
              for s in ("dyadic", "optimized", "qetu")]
    kinds = os.environ.get("KINDS", "chain,cyl").split(",")
    cases = [c for c in cases if c[0] in kinds]
    name = "ext_qetu_validate.json" if len(kinds) == 2 else f"ext_qetu_validate_{kinds[0]}.json"
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_qv_case, cases):
            s = r.get("sim", {})
            print(f"[qetu_val] {r['kind']} N={r['N']} J2={r['J2']} eps={r['eps']:g} {r['schedule']:9s} "
                  f"bound={r['plan']['infidelity_bound']:.2e} meas={s.get('infidelity_trotter', float('nan')):.2e} "
                  f"ok={s.get('bound_ok')}", flush=True)
            out.append(r)
    dump(name, out)


# ----------------------------------------------------------------------------- QETU: scaling (Table IV setting)
def _qs_case(args):
    N, J2, eps, order = args
    rec, est, ans, _ = RS.get_inputs(N, J2, 2)
    ans.statevector = None
    base = dict(N=N, J2=J2, eps=eps, schedule="qetu", order=order, chi=2, gamma=ans.gamma)
    try:
        plan = refine(j1j2_chain(N, 1.0, J2), est, ans, Config(eps=eps, schedule="qetu", order=order))
    except RuntimeError as e:
        return dict(base, plan=None, error=str(e))
    return dict(base, plan=plan.summary())


def run_qetu_scale():
    cases = [(N, J2, eps, o) for N in RS.N_SCALE for J2 in RS.J2S for eps in (1e-2, 1e-3) for o in (1, 2)]
    with ProcessPoolExecutor(JOBS) as ex:
        out = list(ex.map(_qs_case, cases))
    dump("ext_qetu_scaling.json", out)


# ----------------------------------------------------------------------------- regime: MPS alone vs filters
EPS_REGIME = (1e-2, 1e-3, 1e-4)


def _regime_case(args):
    kind, size, J2, eps, chi, sched = args
    rec, est, ans, _ = get_case(kind, size, J2, chi)
    ans.statevector = None
    H = model(kind, size, J2)
    base = dict(kind=kind, size=size, N=H.num_qubits, J2=J2, eps=eps, chi=chi, schedule=sched,
                gamma=ans.gamma, gamma_ed=rec.get("ed", {}).get("gamma", {}).get(str(chi)),
                trial_cost=ans.cost)
    t0 = time.time()
    try:
        plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    except (RuntimeError, ValueError) as e:
        return dict(base, plan=None, error=str(e), time=time.time() - t0)
    return dict(base, plan=plan.summary(), time=time.time() - t0)


def run_regime():
    cases = [("chain", N, J2, eps, chi, s) for N in (16, 32, 48) for J2 in RS.J2S
             for chi in (2, 4, 8, 16, 32) for eps in EPS_REGIME for s in ("auto", "qetu")]
    cases += [("cyl", Lx, J2, eps, chi, s) for Lx in LX for J2 in J2S_2D
              for chi in CHIS_2D for eps in EPS_REGIME for s in ("auto", "qetu")
              if os.path.exists(os.path.join(CACHE2D, _k2(Lx, J2) + ".json"))]
    kinds = os.environ.get("KINDS", "chain,cyl").split(",")
    cases = [c for c in cases if c[0] in kinds]
    name = "ext_regime.json" if len(kinds) == 2 else f"ext_regime_{kinds[0]}.json"
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_regime_case, cases):
            pl = r["plan"]
            print(f"[regime] {r['kind']} N={r['N']} J2={r['J2']} chi={r['chi']} eps={r['eps']:g} "
                  f"{r['schedule']:5s} gamma={r['gamma']:.4f} "
                  f"cnot={pl['res_expected_cnot'] if pl else float('nan'):.2e} ({r['time']:.0f}s)", flush=True)
            out.append(r)
    dump(name, out)


# ----------------------------------------------------------------------------- wider sensitivity study
def _sens_wide_case(args):
    kind, f, N, J2, sched = args
    from dataclasses import replace
    rec, est, ans, g = get_case("chain", N, J2, 2)
    H = j1j2_chain(N, 1.0, J2)
    e = replace(est, E1=est.E0 + f * est.gap, E1_err=0.0) if kind == "gap" else est
    a = replace(ans, gamma=min(1.0, f * ans.gamma)) if kind == "gamma" else ans
    try:
        plan = refine(H, e, a, Config(eps=1e-2, schedule=sched, order=2))
    except (RuntimeError, ValueError) as err:
        return dict(kind=kind, f=f, N=N, J2=J2, schedule=sched, plan=None, error=str(err))
    sim = validate(plan, H, ans, g)            # always the TRUE trial state and ground state
    return dict(kind=kind, f=f, N=N, J2=J2, schedule=sched, plan=plan.summary(), sim=sim)


def run_sens_wide():
    """Gap and overlap overestimation for every N <= 14, coupling and second-order schedule."""
    cases = [("gap", f, N, J2, s) for f in (1.0, 1.1, 1.25, 1.5, 2.0, 3.0)
             for N in RS.N_VAL for J2 in RS.J2S for s in ("dyadic", "optimized", "qetu")]
    cases += [("gamma", f, N, J2, s) for f in (1.02, 1.05)
              for N in RS.N_VAL for J2 in RS.J2S for s in ("dyadic", "optimized", "qetu")]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_sens_wide_case, cases):
            sm = r.get("sim") or {}
            print(f"[sens] {r['kind']} f={r['f']} N={r['N']} J2={r['J2']} {r['schedule']:9s} "
                  f"meas={sm.get('infidelity_trotter', float('nan')):.2e}", flush=True)
            out.append(r)
    dump("ext_sensitivity.json", out)


# ----------------------------------------------------------------------------- what drives the cost: gap and overlap
J2S_GAP = [0.0, 0.1, 0.2, 0.2411, 0.3, 0.4, 0.5]
N_DRIVE = [8, 10, 12, 14, 16, 24, 32, 48]
GAMMAS = [0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97]
SCHED_DRIVE = ("qetu", "dyadic")


def run_inputs_j2():
    """DMRG (+ ED for N <= 14) inputs for the extra couplings of the gap study, in run_study's cache."""
    todo = [(N, J2) for J2 in J2S_GAP for N in sorted(N_DRIVE, reverse=True)]
    todo = [t for t in todo if not os.path.exists(os.path.join(CACHE, RS._key(*t) + ".json"))]
    with ProcessPoolExecutor(JOBS) as ex:
        for k, dt in ex.map(RS.compute_inputs, todo):
            print(f"[inputs_j2] {k} {dt:.0f}s", flush=True)


def _chain_est(N, J2):
    rec = json.load(open(os.path.join(CACHE, RS._key(N, J2) + ".json")))
    return rec, EnergyEstimates(**rec["est"])


def _drive_plan(N, J2, est, gamma, eps, sched):
    """Filter-only plan (no trial-circuit cost) for an assumed overlap gamma."""
    from dataclasses import replace
    H = j1j2_chain(N, 1.0, J2)
    ans = Ansatz(gamma=gamma, cost={}, name=f"gamma={gamma:g}")
    try:
        plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    except (RuntimeError, ValueError) as err:
        return None, str(err)
    return plan, None


def _overlap_case(args):
    N, J2, eps, sched, gamma, gapmode = args
    from dataclasses import replace
    rec, est = _chain_est(N, J2)
    if gapmode == "common":                    # one gap for every N: the smallest true gap of this J2
        g = min(_chain_est(n, J2)[1].gap for n in N_DRIVE)
        est = replace(est, E1=est.E0 + g, E1_err=0.0)
    plan, err = _drive_plan(N, J2, est, gamma, eps, sched)
    return dict(study="overlap", N=N, J2=J2, eps=eps, schedule=sched, gamma=gamma, gapmode=gapmode,
                gap=est.E1 - est.E0 - est.E1_err, plan=plan.summary() if plan else None, error=err)


def _gap_case(args):
    N, J2, eps, sched, gamma, f = args
    from dataclasses import replace
    rec, est = _chain_est(N, J2)
    true_gap = est.gap
    est = replace(est, E1=est.E0 + f * true_gap, E1_err=0.0)
    plan, err = _drive_plan(N, J2, est, gamma, eps, sched)
    return dict(study="gap", N=N, J2=J2, eps=eps, schedule=sched, gamma=gamma, f=f, true_gap=true_gap,
                gap=f * true_gap, plan=plan.summary() if plan else None, error=err)


def run_drivers():
    """Study 1: overlap scan at fixed gap (true gap per N, or one common gap for all N).
    Study 2: gap scan at fixed overlap (J2 scan, and a fraction f of the true gap at fixed H)."""
    c1 = [(N, J2, eps, s, g, gm) for N in N_DRIVE for J2 in RS.J2S for eps in (1e-2, 1e-3)
          for s in SCHED_DRIVE for g in GAMMAS for gm in ("true", "common")]
    c2 = [(N, J2, eps, s, g, f) for N in N_DRIVE for J2 in J2S_GAP for eps in (1e-2, 1e-3)
          for s in SCHED_DRIVE for g in (0.1, 0.5) for f in (0.25, 0.5, 0.75, 1.0)]
    with ProcessPoolExecutor(JOBS) as ex:
        o1 = list(ex.map(_overlap_case, c1, chunksize=8))
        o2 = list(ex.map(_gap_case, c2, chunksize=8))
    dump("ext_drivers.json", o1 + o2)


def _synthetic_trial(N, J2, gamma):
    """sqrt(gamma) |E0> + sqrt(1-gamma) |xi>, |xi> the normalised excited part of the chi = 2 MPS."""
    rec, est, ans, g = get_case("chain", N, J2, 2)
    t = ans.statevector
    xi = t - np.vdot(g, t) * g
    xi /= np.linalg.norm(xi)
    ph = np.vdot(g, t) / abs(np.vdot(g, t))
    return est, np.sqrt(gamma) * g * ph + np.sqrt(1 - gamma) * xi, g


def _drive_val_case(args):
    study, N, J2, eps, sched, gamma, knob = args
    from dataclasses import replace
    est, psi, g = _synthetic_trial(N, J2, gamma)
    if study == "gap":
        est = replace(est, E1=est.E0 + knob * est.gap, E1_err=0.0)
    H = j1j2_chain(N, 1.0, J2)
    ans = Ansatz(gamma=gamma, statevector=psi, cost={})
    plan = refine(H, est, ans, Config(eps=eps, schedule=sched, order=2))
    sim = validate(plan, H, ans, g)
    return dict(study=study, N=N, J2=J2, eps=eps, schedule=sched, gamma=gamma, knob=knob,
                plan=plan.summary(), sim=sim)


def run_drive_val():
    """Simulate a subset of both studies with synthetic trial states of exactly the stated overlap."""
    c = [("overlap", N, J2, 1e-2, s, g, None) for N in RS.N_VAL for J2 in RS.J2S
         for s in SCHED_DRIVE for g in (0.03, 0.3, 0.9)]
    c += [("gap", N, J2, 1e-2, s, 0.5, f) for N in RS.N_VAL for J2 in RS.J2S
          for s in SCHED_DRIVE for f in (0.5, 1.0)]
    out = []
    with ProcessPoolExecutor(JOBS) as ex:
        for r in ex.map(_drive_val_case, c):
            print(f"[drive_val] {r['study']} N={r['N']} J2={r['J2']} {r['schedule']} gamma={r['gamma']} "
                  f"knob={r['knob']} bound={r['plan']['infidelity_bound']:.2e} "
                  f"meas={r['sim']['infidelity_trotter']:.2e} ok={r['sim']['bound_ok']}", flush=True)
            out.append(r)
    dump("ext_drivers_validate.json", out)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    steps = dict(inputs2d=run_inputs2d, chi_chain=run_chi_chain, chain_small=run_chain_small,
                 qetu_val=run_qetu_val, qetu_scale=run_qetu_scale, regime=run_regime,
                 sens_wide=run_sens_wide, inputs_j2=run_inputs_j2, drivers=run_drivers,
                 drive_val=run_drive_val)
    for k in (steps if what == "all" else what.split(",")):
        t = time.time()
        steps[k]()
        print(f"== {k} done in {time.time() - t:.0f}s", flush=True)
