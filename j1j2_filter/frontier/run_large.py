"""Planned SBC / QETU filters on J1-J2 cylinders (4 x Lx) beyond a workstation, for Frontier (OLCF).

Two stages, so the CPU-heavy DMRG and the GPU-heavy statevector simulation can run as separate jobs:

  inputs    (CPU, one system per task)  DMRG energies E0 and E1 (lowest S^z = 1 state) on a bond-dimension
            ladder, the variational trial MPS for each chi, its overlap gamma with the DMRG ground state, and a
            reference ground state (the last rung of the ladder). Dense vectors are written as .npy so the GPU
            stage needs no DMRG.  -> <out>/inputs/cyl<Lx>_J2<J2>.{json,npz}

  circuits  (GPU, one case per GCD)     for each (J2, eps, schedule): pick the trial chi whose certified plan
            is cheapest (or --chi), plan the filter with sbpipe.refine from the DMRG inputs, and simulate the
            Trotterised planned circuit at the size-transfer step: dt_prev, dt_prev/2, dt_prev/4, dt_prev/8
            (stopping once a halved step passes), dt_prev = the coarsest passing step of the next-smaller
            cylinder (from --prev, a results file of this script or of sbpipe_study). Records infidelity to the
            reference ground state, success probability, post-selected energy and CNOT counts.
            -> <out>/circuits_cyl<Lx>.jsonl (one json row per case; same format as
               sbpipe_study/results/transfer20.jsonl)

  selftest  checks the simulator against sbpipe on 12 qubits.

Typical:
    python frontier/run_large.py inputs   --Lx 6 --J2 0 0.3
    python frontier/run_large.py circuits --Lx 6 --J2 0 0.3 --prev sbpipe_study/results/cross20.jsonl
The infidelity is measured against the DMRG reference, so it cannot resolve below that reference's own error
(stored as ref_err in the inputs json; keep it well below the smallest eps).
"""
import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import numpy as np

LY = 4
EPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)


def key(Lx, J2):
    return f"cyl{Lx}_J2{J2:g}"


def model(Lx, J2):
    from sbpipe.models import j1j2_square
    return j1j2_square(Lx, LY, 1.0, J2)


# ----------------------------------------------------------------------------- stage 1: DMRG inputs (CPU)
def inputs(Lx, J2, out, chis, ladder):
    from sbpipe.frontends import dmrg_quimb as DQ
    d = os.path.join(out, "inputs")
    os.makedirs(d, exist_ok=True)
    t0 = time.time()
    H = model(Lx, J2)
    N = H.num_qubits
    mpo = DQ.pauli_to_mpo(H)
    hist, psi, prev = [], None, None
    for chi in ladder:                                   # ground-state ladder, keeping the last two rungs
        E, st, _ = DQ._dmrg(mpo, [chi], 1e-12, 1e-11, 60, p0=psi)
        prev, psi = psi, st
        hist.append((chi, E, int(st.max_bond())))
        print(f"[{key(Lx, J2)}] ground chi={chi} E0={E:.10f} {time.time() - t0:.0f}s", flush=True)
    g = dict(E0=hist[-1][1], history=hist, E0_conv=abs(hist[-1][1] - hist[-2][1]) if len(hist) > 1 else float("nan"))
    e = DQ.sector_ground_state(mpo, DQ.total_sz(N), 1.0, 2.0, bond_dims=ladder)
    print(f"[{key(Lx, J2)}] DMRG E0={g['E0']:.10f} E1={e['E0']:.10f}  {time.time() - t0:.0f}s", flush=True)
    # reference error: infidelity between the last two rungs of the ladder
    ov = (abs(DQ.mps_overlap(prev, psi)) ** 2 / (abs(DQ.mps_overlap(prev, prev)) * abs(DQ.mps_overlap(psi, psi)))
          if prev is not None else 0.0)
    rec = dict(Lx=Lx, Ly=LY, N=N, J2=J2, ladder=list(ladder),
               est=dict(E0=g["E0"], E1=e["E0"], E0_err=max(g["E0_conv"], 1e-10), E1_err=max(e["E0_conv"], 1e-10),
                        E_max=None, source=f"DMRG chi<={max(ladder)}; E1 = lowest state with Q=1"),
               E0_history=g["history"], E1_history=e["history"], ref_err=float(1 - ov), gamma={}, E_trial={})
    np.save(os.path.join(d, key(Lx, J2) + "_ground.npy"), np.asarray(psi.to_dense(), complex).reshape(-1))
    for chi in chis:
        E_t, t = DQ.trial_mps(mpo, chi)
        gam = abs(DQ.mps_overlap(psi, t)) ** 2 / (abs(DQ.mps_overlap(psi, psi)) * abs(DQ.mps_overlap(t, t)))
        rec["gamma"][str(chi)] = float(gam)
        rec["E_trial"][str(chi)] = float(E_t)
        v = np.asarray(t.to_dense(), complex).reshape(-1)
        np.save(os.path.join(d, key(Lx, J2) + f"_trial{chi}.npy"), v / np.linalg.norm(v))
        print(f"[{key(Lx, J2)}] chi={chi} 1-gamma={1 - gam:.2e}", flush=True)
    rec["time"] = time.time() - t0
    json.dump(rec, open(os.path.join(d, key(Lx, J2) + ".json"), "w"), indent=1)


# ----------------------------------------------------------------------------- stage 2: circuits (GPU)
def coarsest_pass(row):
    c, best = sorted(row["curve"], key=lambda q: -q["dt"]), None
    for i in range(len(c) - 1, -1, -1):
        if c[i]["inf"] > row["eps"]:
            break
        best = c[i]["dt"]
    return best


def prev_steps(paths, Lx):
    """{(J2, eps, schedule): coarsest passing dt} of the cylinder 4 x (Lx - 1) from earlier results files."""
    P = {}
    for p in paths:
        for r in map(json.loads, open(p)):
            if r.get("kind", "cyl") == "cyl" and r["size"] == Lx - 1 and r["curve"]:
                dt = coarsest_pass(r)
                if dt is not None:
                    k = (r["J2"], r["eps"], r["schedule"])
                    P[k] = max(P.get(k, 0.0), dt)
    return P


def plan_for(H, rec, chi, eps, sched):
    from sbpipe import Ansatz, EnergyEstimates
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    from sbpipe.pipeline import Config, refine
    ans = Ansatz(gamma=rec["gamma"][str(chi)], statevector=None, cost=mps_circuit_cost(H.num_qubits, chi),
                 name=f"chi={chi}", meta=dict(chi=chi))
    return refine(H, EnergyEstimates(**rec["est"]), ans, Config(eps=eps, schedule=sched, order=2))


def choose_chi(H, rec, eps, sched):
    """Cheapest certified plan over the trial chis that do not meet eps alone (the study's rule)."""
    best = None
    for c in sorted(map(int, rec["gamma"])):
        if 1 - rec["gamma"][str(c)] <= eps:
            continue
        try:
            cost = plan_for(H, rec, c, eps, sched).resources["expected_cnot"]
        except RuntimeError:
            continue
        if best is None or cost < best[1]:
            best = (c, cost)
    return best[0] if best else None


def circuit_case(Lx, J2, eps, sched, chi, dt_prev, out):
    import gpusim as G
    from sbpipe.resources import filter_counts
    t0 = time.time()
    H = model(Lx, J2)
    d = os.path.join(out, "inputs")
    rec = json.load(open(os.path.join(d, key(Lx, J2) + ".json")))
    if chi is None:
        chi = choose_chi(H, rec, eps, sched)
    row = dict(kind="cyl", size=Lx, N=H.num_qubits, J2=J2, eps=eps, schedule=sched, chi=chi, dt_prev=dt_prev,
               ref_err=rec["ref_err"], curve=[])
    if chi is None:
        row["note"] = "every trial meets eps alone, or no plan"
        return row
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    cost = mps_circuit_cost(H.num_qubits, chi)
    plan = plan_for(H, rec, chi, eps, sched)
    des, Hs, W = plan.design, plan.H_scaled, plan.window["W"]
    row.update(design=des.schedule, m=des.m, W=W, gamma=rec["gamma"][str(chi)], trial_cnot=cost["cnot"],
               cert_steps=int(np.sum(plan.steps)), cert_cnot=plan.resources["expected_cnot"],
               cert_p=plan.resources["p_succ"])
    if des.m == 0:
        row["time"] = time.time() - t0
        return row
    S, Sh = G.BondStrang(Hs), G.BondStrang(H)
    psi = G.to_dev(np.load(os.path.join(d, key(Lx, J2) + f"_trial{chi}.npy")))
    gnd = np.load(os.path.join(d, key(Lx, J2) + "_ground.npy"), mmap_mode="r")
    anc = 2 if des.schedule == "qetu" else 1
    t = np.abs(des.times)
    for f in (1, 2, 4, 8):
        dtp = dt_prev / f
        k = np.maximum(1, np.ceil(t / (dtp * W) - 1e-12)).astype(int)
        if des.schedule == "qetu":
            dm = des.meta
            v = G.qetu_trotter(S, psi, dm["degree"], dm["tau"], dm["x_s"], k)
        else:
            v = G.filter_trotter(S, psi, des.times, des.phases, k)
        p = float(G.xp.vdot(v, v).real)
        E = Sh.energy(v)
        g = G.to_dev(gnd)
        inf = max(0.0, 1.0 - abs(complex(G.xp.vdot(g, v))) ** 2 / p / float(G.xp.vdot(g, g).real))
        del g, v
        cn = filter_counts(Hs, k, 2, anc_rot_per_pulse=anc)["cnot"]
        row["curve"].append(dict(dt=dtp, steps=int(k.sum()), filter_cnot=cn, inf=inf, p=p, E=E,
                                 cnot=(cn + cost["cnot"]) / p))
        print(f"[{key(Lx, J2)} eps={eps:g} {sched} chi={chi}] dt={dtp:g} steps={int(k.sum())} inf={inf:.2e} "
              f"p={p:.3f} {time.time() - t0:.0f}s", flush=True)
        if f >= 2 and inf <= eps:
            break
    row["time"] = time.time() - t0
    return row


def circuits(a):
    P = prev_steps(a.prev, a.Lx)
    cases = [(J2, e, s) for J2 in a.J2 for e in a.eps for s in a.schedule]
    rank, size = int(os.environ.get("SLURM_PROCID", 0)), int(os.environ.get("SLURM_NTASKS", 1))
    path = os.path.join(a.out, f"circuits_cyl{a.Lx}.jsonl" if size == 1 else f"circuits_cyl{a.Lx}.rank{rank}.jsonl")
    done = set()
    if os.path.exists(path):
        done = {(r["J2"], r["eps"], r["schedule"]) for r in map(json.loads, open(path))}
    mine = [c for i, c in enumerate(cases) if i % size == rank and c not in done]
    for J2, e, s in mine:
        dt = a.dt if a.dt else P.get((J2, e, s))
        if dt is None:
            print(f"[skip] no step for J2={J2} eps={e} {s} at 4x{a.Lx - 1}: pass --prev or --dt", flush=True)
            continue
        r = circuit_case(a.Lx, J2, e, s, a.chi, dt, a.out)
        with open(path, "a") as fh:
            fh.write(json.dumps(r) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=("inputs", "circuits", "selftest"))
    ap.add_argument("--Lx", type=int, default=6, help="cylinder length (4 x Lx sites); 6 -> 24 qubits, 7 -> 28")
    ap.add_argument("--J2", type=float, nargs="+", default=[0.0, 0.3])
    ap.add_argument("--eps", type=float, nargs="+", default=list(EPS))
    ap.add_argument("--schedule", nargs="+", default=["auto", "qetu"], choices=("auto", "qetu", "dyadic", "optimized"))
    ap.add_argument("--chis", type=int, nargs="+", default=[8, 16, 32, 64], help="trial bond dimensions (inputs)")
    ap.add_argument("--ladder", type=int, nargs="+", default=[16, 32, 64, 128, 256, 512, 1024],
                    help="DMRG bond-dimension ladder; the last rung is the reference ground state")
    ap.add_argument("--chi", type=int, default=None, help="fix the trial chi (default: cheapest certified plan)")
    ap.add_argument("--prev", nargs="*", default=[], help="results files holding the 4 x (Lx-1) curves")
    ap.add_argument("--dt", type=float, default=None, help="use this dt_prev for every case instead of --prev")
    ap.add_argument("--out", default=os.path.join(os.getcwd(), "frontier_out"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.stage == "selftest":
        import gpusim as G
        print(G.selftest())
        return
    if a.stage == "inputs":
        rank, size = int(os.environ.get("SLURM_PROCID", 0)), int(os.environ.get("SLURM_NTASKS", 1))
        for i, J2 in enumerate(a.J2):
            if i % size == rank:
                inputs(a.Lx, J2, a.out, a.chis, tuple(a.ladder))
        return
    circuits(a)


if __name__ == "__main__":
    main()
