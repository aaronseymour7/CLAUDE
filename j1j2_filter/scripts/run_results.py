"""
run_results.py -- collect every pipeline quantity needed for the results document.

For each (source, N, J2, eps) it runs the pipeline's own functions (build_ctx -> make_design -> find_guaranteed /
find_empirical) and records, in results/data/<source>_N<N>_J2<J2>_eps<eps>.json:

  inputs     E0, E1, Etop, W, scaled gap, gamma (+ source), alpha, <DMRG|ED>^2, MPO check, trial F
  design     source, pulses m, x = T*Delta/pi, T, eta target / certified, leakage budget
  per design ("guaranteed" = rigorous a-priori bound <= eps; "empirical" = measured 1-F <= eps, an oracle):
     n, pulses used, k_i, t_i, phases, dt, T (scaled and physical), eta, bound terms, P_succ lower bound,
     measured 1-F / F / P_succ / state error (when the point is simulated), full-transpile check (n <= 200)
     CX, RZ (all and non-Clifford), H, S, depth, approx. T-count -- trial, filter, total, all-to-all and line
  per-step    transpiled counts of ONE Trotter step of one pulse (and the 3-step average), all-to-all and line

All runs use the DMRG spectrum (DMRG gap, certified bandwidth W). N <= 14: ED cross-check and statevector simulation available;
N > 14: DMRG/MPS only (no ED, no statevector, no dense Hamiltonian), circuits are counted, not simulated.

usage:  python scripts/run_results.py [--worker i --of K] [--only ed|dmrg|large] [--list]
Resumable: finished points are skipped. Run from anywhere; set OMP_NUM_THREADS=1 when running several workers.
"""
import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

import numpy as np                                    # noqa: E402
import pipeline as P                                  # noqa: E402
from core.circuits import apply_filter_pulse          # noqa: E402
from core.resources import resource_costs, rotation_synthesis_t_count   # noqa: E402
from qiskit import QuantumCircuit                     # noqa: E402
from qiskit.transpiler import CouplingMap             # noqa: E402

import get_energies as GE                             # noqa: E402   (settings + RUN_LOG live in get_energies.py)

_ROLES = ("ground", "top", "excited")

OUT = ROOT / "results" / "data"
OUT.mkdir(parents=True, exist_ok=True)

SMALL_N = (4, 6, 8, 10, 12, 14)       # statevector-simulated (measured 1-F, P_succ); ED cross-check available
J2_ALL = (0.0, 0.2411, 0.4)
EPS_MAIN = 1e-2
EPS_SCAN = (1e-1, 3e-2, 3e-3, 1e-3)
SCAN_CASES = ((6, 0.0), (8, 0.0), (6, 0.4))
LARGE_N = (16, 18, 20, 24, 28, 32)                    # DMRG/MPS only: no ED, no statevector, circuits counted not simulated
SIM_MAX_FOR_N = {4: 200_000, 6: 200_000, 8: 200_000, 10: 200_000, 12: 200_000, 14: 150_000}   # n above this: not simulated


def jobs():
    """(source, N, J2, [(eps, [which...])]).  Everything uses the DMRG spectrum (certified bandwidth)."""
    out = []
    for N in SMALL_N:
        for J2 in J2_ALL:
            ep = [(EPS_MAIN, ["guaranteed", "empirical"])]
            if (N, J2) in SCAN_CASES:
                ep += [(e, ["guaranteed", "empirical"]) for e in EPS_SCAN]
            out.append(("dmrg", N, J2, ep))
    for N in LARGE_N:
        for J2 in J2_ALL:
            out.append(("dmrg", N, J2, [(1e-1, ["guaranteed"]), (EPS_MAIN, ["guaranteed"])]))
    return out


# ----------------------------------------------------------------------------- counting
COUNT_KEYS = ("cx", "h", "s", "rz", "rz_nonclifford", "total", "depth")


def pick(c):
    return {k: float(c[k]) for k in COUNT_KEYS} | {"t_count_est": float(c["t_count_est"])}


def step_counts(H_scaled, N):
    """Transpiled counts of one Trotter step of one pulse (phi = 0), all-to-all and on a line (ancilla = qubit N,
    the end of the line: pessimistic). steps=1 and the 3-step average (routing / cancellation across steps)."""
    out = {}
    for label, cmap in (("a2a", None), ("line", CouplingMap.from_line(N + 1))):
        per = {}
        for reps in (1, 3):
            qc = QuantumCircuit(N + 1)
            apply_filter_pulse(qc, list(range(N)), N, H_scaled, 1.0, 0.0, reps, P.ORDER)
            c = resource_costs(qc, verbose=False, coupling_map=cmap, synth_eps=P.SYNTH_EPS)
            per[reps] = {k: float(c[k]) / reps for k in COUNT_KEYS}
        out[label] = dict(one=per[1], avg3=per[3])
    return out


def totals(ctx, pt, step, topo):
    """Filter and total (with trial) counts for the success path = n * (per-step) + per-pulse Rz(2 phi) + 2 H."""
    s = step[topo]["one"] if topo == "a2a" else step[topo]["avg3"]
    n, nz = pt["n"], pt["n_nz"]
    tc = ctx["trial_cost"]
    f = {k: s[k] * n for k in COUNT_KEYS if k != "depth"}
    f["rz"] += nz                                    # one Rz(2 phi) per non-zero pulse
    f["rz_nonclifford"] += nz
    f["h"] += 2 * nz                                 # ancilla Hadamards
    f["depth"] = s["depth"] * n
    f["total"] = f["cx"] + f["h"] + f["s"] + f["rz"]
    t_total = {k: f[k] + float(tc[k]) for k in ("cx", "h", "s", "rz", "rz_nonclifford", "total") if k in tc}
    t_total["depth"] = f["depth"] + float(tc["depth"])
    t_total["t_count_est"] = rotation_synthesis_t_count(f["rz_nonclifford"], P.SYNTH_EPS) + float(tc["t_count_est"])
    f["t_count_est"] = rotation_synthesis_t_count(f["rz_nonclifford"], P.SYNTH_EPS)
    return dict(filter=f, total=t_total)


def describe(ctx, des, pt, step, which, W):
    d = dict(which=which, n=int(pt["n"]), n_pulses_used=int(pt["n_nz"]), n_pulses_designed=int(pt["n_designed"]),
             k=[int(x) for x in pt["k"]], tg=[float(x) for x in pt["tg"]], phases=[float(x) for x in pt["ph"]],
             dt=float(pt["dt"]), T_scaled=float(np.sum(pt["tg"])), T_phys=float(np.sum(pt["tg"]) / W),
             eta=float(pt["eta"]), grid_ok=bool(pt["grid_ok"]), eps_bound=float(pt["eps_bound"]),
             leak_bound=float(pt["leak"]), d_T=float(pt["d_T"]), eps_T=float(pt["eps_T"]),
             p_succ_lb=float(pt["p_g_lb"]), f0=float(pt["f0"]), simulated=bool(pt["simulated"]),
             kmax=int(np.max(pt["k"])))
    if pt["simulated"]:
        d.update(eps_actual=float(pt["eps_actual"]), F_ed=float(pt["F_ed"]), F_dmrg=float(pt["F_dmrg"]),
                 p_succ=float(pt["p_succ"]), p_succ_exact=float(pt["p_succ_exact"]), state_err=float(pt["state_err"]),
                 circ_dev=None if pt["circ_dev"] is None else float(pt["circ_dev"]), warn=pt["warn"])
    d["counts_a2a"] = totals(ctx, pt, step, "a2a")
    d["counts_line"] = totals(ctx, pt, step, "line")
    try:
        m = P.measured_cost(ctx, pt)                 # full transpile of the actual early-abort circuit (n <= 200)
    except Exception as e:
        m = dict(error=f"{type(e).__name__}: {e}")
    d["full_transpile_check"] = m
    return d


def validation(ctx, des, step, W, n=48):
    """Per-step x n estimate against a full transpile of the real circuit at a fixed small n."""
    try:
        pt = P.point(ctx, des, n, simulate=False)
        est = totals(ctx, pt, step, "a2a")["filter"]
        m = P.measured_cost(ctx, pt)
        return dict(n=n, est_cx=est["cx"], est_depth=est["depth"], full_cx=m["cx"], full_depth=m["depth"],
                    cx_err_pct=100 * (m["cx"] - est["cx"]) / est["cx"])
    except Exception as e:
        return dict(n=n, error=f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------------- one case
def run_point(source, N, J2, eps, whichs, cache):
    P.SPECTRUM_SOURCE = source
    P.SIM_MAX_N = SIM_MAX_FOR_N.get(N, 0)
    key = (source, N, J2)
    if key not in cache:
        t0 = time.time()
        GE.RUN_LOG.clear()
        ctx = P.run_quiet(P.build_ctx, N, P.J1, J2)
        ctx["dmrg_runs"] = [dict(role=_ROLES[i % 3], **x) for i, x in enumerate(GE.RUN_LOG)]
        step = step_counts(ctx["H_scaled"], N)
        cache.clear()
        cache[key] = dict(ctx=ctx, step=step, build_s=time.time() - t0, des={})
    c = cache[key]
    ctx, step = c["ctx"], c["step"]
    spec = ctx["spec"]
    W = float(spec["W"])
    P._STEP_CX["cx"] = ctx["step"]["cx"]
    rec = dict(source=source, N=N, J2=J2, eps=eps, order=P.ORDER, bandwidth=P.BANDWIDTH, sweep_L=P.SWEEP_L, seed=P.SEED,
               inputs=dict(E0=float(spec["E0"]), E1=float(spec["E1"]), Etop=float(spec["Etop"]), W=W,
                           W_cert=float(spec.get("W_cert", np.nan)), gap_scaled=float(spec["gap"]),
                           gap_raw=float(spec["gap"] * W), gamma=float(ctx["gamma"]), gamma_src=ctx["gamma_src"],
                           trial_F_ed=float(ctx["g_ed"]), trial_F_dmrg=float(ctx["g_dmrg"]), dmrg_ed_overlap=float(ctx["dmrg_ed"]),
                           alpha=float(ctx["alpha"]), alpha_tri=float(ctx["alpha_tri"]), e0_dev_indicative=float(spec["e0_dev_indicative"]),
                           W_exact_dmrg=float(spec["Etop"] - spec["E0"]), dense=bool(ctx["dense"]),
                           mpo_err=None if ctx.get("mpo_err") is None else float(ctx["mpo_err"]),
                           n_pauli_terms=int(len(ctx["H_scaled"]) - 1), trial_layers=P.SWEEP_L,
                           trial_cost=pick(ctx["trial_cost"]), build_s=c["build_s"], dmrg_runs=ctx.get("dmrg_runs", []),
                           e1_method=(spec.get("dmrg_res") or {}).get("e1_method"),
                           dmrg_vs_ed=({k: float(spec["ed"][k]) for k in ("E0_diff", "E1_diff", "E_top_diff") if k in spec["ed"]}
                                       if spec.get("ed") else None),
                           dmrg_converged=all(x["converged"] for x in ctx.get("dmrg_runs", []))),
               per_step=step, designs={})
    t0 = time.time()
    if eps not in c["des"]:
        c["des"][eps] = P.run_quiet(P.make_design, ctx, eps)
    des, cands, n_inf = c["des"][eps]
    rec["design"] = dict(source=des["source"], m=int(len(des["times"])), x=float(des["x"]), T_scaled=float(des["T"]),
                         T_phys=float(des["T"] / W), eta_cert=float(des["eta_cert"]),
                         eta_target=None if des.get("eta_target") is None else float(des["eta_target"]),
                         leak_budget=None if des.get("leak_budget") is None else float(des["leak_budget"]),
                         times=[float(x) for x in des["times"]], phases=[float(x) for x in des["phases"]],
                         candidates=[dict(eps_leak=float(cd["row"]["eps"]), x=float(cd["row"]["x"]), m=int(cd["row"]["m"]),
                                          T=float(cd["row"]["T"]), eta=float(cd["row"]["eta_cert"]), f0sq=float(cd["row"]["f0sq"]),
                                          n_req=int(cd["n_req"]), cost=float(cd["cost"])) for cd in cands],
                         design_s=time.time() - t0)
    if des["source"] == "none":
        rec["note"] = "gamma = 1: trial state already exact, no filter"
        return rec
    rec["validation"] = validation(ctx, des, step, W)
    for which in whichs:
        t1 = time.time()
        try:
            if which == "guaranteed":
                pt = P.find_guaranteed(ctx, des, eps)
            else:
                pt = P.find_empirical(ctx, des, eps, {})
        except Exception as e:
            rec["designs"][which] = dict(error=f"{type(e).__name__}: {e}", trace=traceback.format_exc()[-1500:])
            continue
        if pt is None:
            rec["designs"][which] = dict(error="not achievable (leakage alone exceeds the budget / not reached in sweep)")
            continue
        d = describe(ctx, des, pt, step, which, W)
        d["search_s"] = time.time() - t1
        rec["designs"][which] = d
    return rec


def fname(source, N, J2, eps):
    return OUT / f"{source}_N{N}_J2{J2:g}_eps{eps:g}.json"


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=1, default=lambda o: o.item() if hasattr(o, "item") else (o.tolist() if hasattr(o, "tolist") else str(o))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--only", default="all", choices=["all", "small", "large"])
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    J = jobs()
    if a.only == "small":
        J = [j for j in J if j[1] <= 14]
    elif a.only == "large":
        J = [j for j in J if j[1] > 14]
    J = sorted(J, key=lambda j: -(j[1] * 10 + len(j[3])))          # biggest first, spread over workers
    mine = [j for i, j in enumerate(J) if i % a.of == a.worker]
    if a.list:
        for j in mine:
            print(j[0], j[1], j[2], [e for e, _ in j[3]])
        return
    cache = {}
    for source, N, J2, eps_list in mine:
        for eps, whichs in eps_list:
            f = fname(source, N, J2, eps)
            if f.exists():
                continue
            t0 = time.time()
            print(f"[{time.strftime('%H:%M:%S')}] {source} N={N} J2={J2:g} eps={eps:g} {whichs}", flush=True)
            try:
                rec = run_point(source, N, J2, eps, whichs, cache)
            except Exception as e:
                traceback.print_exc()
                dump(f.with_suffix(".err"), dict(source=source, N=N, J2=J2, eps=eps, error=f"{type(e).__name__}: {e}",
                                                 trace=traceback.format_exc()))
                continue
            rec["wall_s"] = time.time() - t0
            dump(f, rec)
            g = rec["designs"].get("guaranteed", {})
            e = rec["designs"].get("empirical", {})
            print(f"   done {time.time() - t0:.0f}s  n_guar={g.get('n')} n_emp={e.get('n')}", flush=True)


if __name__ == "__main__":
    main()
