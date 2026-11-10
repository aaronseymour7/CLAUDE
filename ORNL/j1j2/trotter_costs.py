"""Gate costs of the J1-J2 QPF pipeline with first- and second-order Trotter evolution.

    python j1j2/trotter_costs.py                  # N = 8..16 -> results/j1j2_trotter.json, results/j1j2_trotter.md
    python j1j2/trotter_costs.py --N 8 10

Uses the schedules of results/j1j2_pipeline.json (run j1j2/run_pipeline.py first): the J^2 ladder and, for the energy
stage, QPF(t_i, delta_i) (v3), QPF(t_i) (v4) and the constant schedule, each block iterated as in the pipeline. Every
QPF step is built as a Qiskit circuit (j1j2/trotter.py), with exp(-i t O (x) Z_a) Trotterised by LieTrotter (order 1)
or SuzukiTrotter(order=2), and transpiled at optimization level 3 to {cx, rz, sx, x}. The trial state is the
one-layer mps-to-circuit circuit, transpiled the same way.

Step sizes. Step i of a stage gets n_i = ceil(|t_i| / dt) Trotter steps, with dt_J = (pi/4) / k_J for the J^2 stage
and dt_E for the energy stage, on the grids K_J and DT_E. For each (order, schedule, eps) the pair (k_J, dt_E) and
the number of blocks R are chosen to minimise CNOTs per success = CNOTs per attempt / success probability, subject to
1 - F <= eps with the chosen steps AND with the next finer k_J and the next finer dt_E (so the choice is not a lucky
cancellation of Trotter errors). Fidelity F is to the exact ground state, simulated exactly (full 2^N space, no
symmetry restriction) with the same Trotter product as the circuits (checked against the transpiled circuits in
tests/test_trotter.py). CNOTs per success assume the whole circuit is rerun after a failed ancilla readout; aborting
at the first failed readout makes it cheaper, so these are upper bounds of the expected cost.
"""
import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import inputs, trotter  # noqa: E402
from qpf import design  # noqa: E402

K_J = [1, 2, 3, 4, 6, 8, 12, 16]
DT_E = [32 * 2 ** (-h / 2) for h in range(18)]           # 32 ... 0.088 (dt_E >= max t_i means one Trotter step)
EPS = (1e-2, 1e-3)
SCHEDULES = ("QPF (t_i, δ_i)", "QPF (t_i)", "Constant")


def schedules(r):
    T = r["energy"]["Constant"]["T_block"]
    out = {"QPF (t_i, δ_i)": (r["design_v3"]["times"], r["design_v3"]["phases"]),
           "QPF (t_i)": (r["design_v4"]["times"], r["design_v4"]["phases"])}
    out["Constant"] = design.constant(T, r["design_v4"]["m"])
    blocks = {}
    for k in SCHEDULES:
        rows = r["energy"][k]["rows"]
        need = next((x["block"] for x in rows if x["infid"] <= min(EPS)), rows[-1]["block"])
        blocks[k] = need + 1
    return {k: (np.asarray(t, float), np.asarray(d, float)) for k, (t, d) in out.items()}, blocks


def run_N(N, r, seed=0):
    t0 = time.time()
    ex = inputs.exact(N)
    est = inputs.dmrg_estimates(N, seed=seed)
    qc, psi = inputs.trial(est["psi0"])
    g0 = ex["g0"]
    idx = np.arange(1 << N)
    Hop = trotter.scaled_hamiltonian(N, inputs.bonds(N), est["E0"], est["Emax"])
    Jop = trotter.total_spin_sq(N)
    trial_cx, trial_nc = trotter.counts(qc)
    tJ = r["j2"]["times"]
    scheds, blocks = schedules(r)
    res = dict(N=N, trial=dict(cx=trial_cx, nonclifford_rz=trial_nc, gamma=r["trial"]["gamma"]),
               exact=dict(j2_success=r["j2"]["success"],
                          energy={k: dict(rows=[dict(block=x["block"], p=x["p"], infid=x["infid"])
                                                for x in r["energy"][k]["rows"]]) for k in SCHEDULES}),
               orders={})
    for order in (1, 2):
        # J^2 stage for every k_J, kept as the columns of one array
        VJ, pJ, cJ = [], [], []
        for k in K_J:
            v, ns = trotter.run(Jop, idx, psi, tJ, [0.0] * len(tJ), (np.pi / 4) / k, order)
            pJ.append(float(np.linalg.norm(v) ** 2))
            VJ.append(v / np.sqrt(pJ[-1]))
            c = [trotter.step_counts(Jop, ("J", N), t, 0.0, n, order) for t, n in zip(tJ, ns)]
            cJ.append(dict(ns=ns, cx=sum(x[0] for x in c), nonclifford_rz=sum(x[1] for x in c)))
        VJ = np.stack(VJ, axis=1)
        jstage = [dict(k=k, dt=(np.pi / 4) / k, p=pJ[i], infid=float(1 - abs(np.vdot(g0, VJ[:, i])) ** 2), **cJ[i])
                  for i, k in enumerate(K_J)]
        ores = dict(j2=jstage, energy={})
        for label, (times, phases) in scheds.items():
            grid = []                                         # grid[d][b] = per-k arrays after block b + 1
            prev = None
            seen = set()
            for dt in DT_E:
                ns_dt = tuple(trotter.steps_for(t, dt) for t in times)
                if ns_dt in seen:                             # same circuit as a coarser dt_E already on the grid
                    continue
                seen.add(ns_dt)
                d = len(grid)
                V, P, rows = VJ.copy(), np.ones(len(K_J)), []
                for b in range(blocks[label]):
                    V, ns = trotter.run(Hop, idx, V, times, phases, dt, order)
                    q = np.sum(abs(V) ** 2, axis=0)
                    V, P = V / np.sqrt(q), P * q
                    rows.append(dict(p=P.copy(), infid=1 - abs(g0.conj() @ V) ** 2))
                c = [trotter.step_counts(Hop, ("H", N), t, dl, n, order) for t, dl, n in zip(times, phases, ns)]
                grid.append(dict(dt=dt, ns=ns, cx=sum(x[0] for x in c), nonclifford_rz=sum(x[1] for x in c),
                                 rows=rows))
                cur = np.array([x["infid"] for x in rows])
                conv = prev is not None and np.all(abs(cur - prev) <= 0.01 * np.maximum(prev, 1e-6))
                if conv and d >= 2 and grid[-2].get("conv"):
                    break
                grid[-1]["conv"] = bool(conv)
                prev = cur
            ores["energy"][label] = dict(grid=grid, choice=choose(grid, jstage, trial_cx, trial_nc))
        res["orders"][order] = ores
        print(f"  N={N} order {order} done [{time.time() - t0:.0f} s]", flush=True)
    res["seconds"] = time.time() - t0
    return res


def choose(grid, jstage, trial_cx, trial_nc):
    out = {}
    for eps in EPS:
        best = None
        for d in range(len(grid) - 1):
            for i in range(len(jstage) - 1):
                for b in range(len(grid[d]["rows"])):
                    ok = all(grid[dd]["rows"][b]["infid"][ii] <= eps for dd, ii in ((d, i), (d + 1, i), (d, i + 1)))
                    if not ok:
                        continue
                    R = b + 1
                    P = jstage[i]["p"] * grid[d]["rows"][b]["p"][i]
                    cx = trial_cx + jstage[i]["cx"] + R * grid[d]["cx"]
                    nc = trial_nc + jstage[i]["nonclifford_rz"] + R * grid[d]["nonclifford_rz"]
                    if best is None or cx / P < best["cx_per_success"]:
                        best = dict(k_J=jstage[i]["k"], n_J=jstage[i]["ns"], dt_E=grid[d]["dt"], n_E=grid[d]["ns"],
                                    blocks=R, infid=float(grid[d]["rows"][b]["infid"][i]), success=float(P),
                                    cx_attempt=int(cx), nonclifford_rz_attempt=int(nc),
                                    cx_per_success=cx / P, nonclifford_rz_per_success=nc / P,
                                    cx_trial=trial_cx, cx_j2=jstage[i]["cx"], cx_block=grid[d]["cx"],
                                    nc_trial=trial_nc, nc_j2=jstage[i]["nonclifford_rz"],
                                    nc_block=grid[d]["nonclifford_rz"])
        out[str(eps)] = best
    return out


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def table(res):
    L = ["| N | order | ε | schedule | k_J, dt_E | blocks | success | 1−F | CNOT / attempt (trial + J² + blocks) "
         "| non-Clifford Rz / attempt | CNOT / success | non-Clifford Rz / success |",
         "|--|--|--|--|--|--|--|--|--|--|--|--|"]
    for r in res:
        for order, o in r["orders"].items():
            for eps in EPS:
                for label in SCHEDULES:
                    c = o["energy"][label]["choice"][str(eps)]
                    if c is None:
                        L.append(f"| {r['N']} | {order} | {eps:g} | {label} | not reached | | | | | | | |")
                        continue
                    L.append(f"| {r['N']} | {order} | {eps:g} | {label} | {c['k_J']}, {c['dt_E']:.3g} | {c['blocks']} "
                             f"| {c['success']:.3f} | {c['infid']:.1e} | {c['cx_attempt']} ({c['cx_trial']} + {c['cx_j2']} + "
                             f"{c['blocks']}×{c['cx_block']}) | {c['nonclifford_rz_attempt']} | {c['cx_per_success']:.0f} "
                             f"| {c['nonclifford_rz_per_success']:.0f} |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12, 14, 16])
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    pipe = {x["N"]: x for x in json.load(open(os.path.join(HERE, "results", "j1j2_pipeline.json")))}
    with ProcessPoolExecutor(a.workers) as ex:
        res = list(ex.map(run_N, a.N, [pipe[N] for N in a.N]))
    res = _clean(res)
    json.dump(res, open(os.path.join(HERE, "results", "j1j2_trotter.json"), "w"), indent=1)
    open(os.path.join(HERE, "results", "j1j2_trotter.md"), "w").write(table(res) + "\n")
    print(table(res))


if __name__ == "__main__":
    main()
