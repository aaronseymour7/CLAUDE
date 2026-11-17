"""Gate cost of different routes to the J1-J2 ground state, with Trotterised QPF circuits.

    python j1j2/route_costs.py                    # N = 8..16 -> results/j1j2_routes.{json,md}
    python j1j2/route_costs.py --N 8 10 --workers 2

Routes (all start from the mps-to-circuit trial with L layers, L in LAYERS):
  S2full -> H   J^2 ladder (pi/4, pi/8, ...) to J = 0, then the energy filter for the singlet gap   (the pipeline)
  S2(pi/4) -> H only the first J^2 pulse (t = pi/4: removes S = 1, 2, 5, 6), then the singlet-gap filter
  H -> S2(pi/4) the energy filter for the lowest (triplet) gap first, then the pi/4 J^2 pulse
  H only        the energy filter for the triplet gap, no J^2 stage
Energy filter: the authors' design (qpf.design.optimize), QPF(t, delta) (v3) and QPF(t) (v4), m chosen by design
cost, block iterated R times. Gaps from DMRG: singlet gap as in the pipeline; triplet gap from DMRG on
H + 2 lam |psi0><psi0| (projector on the compressed DMRG ground state), whose lowest level is the first triplet.

Circuits: second-order Trotter (SuzukiTrotter(order=2)); first order is worse at these targets (results/j1j2_trotter*).
One Trotter step size per stage, set by the longest pulse (others snapped to whole steps; j1j2/trotter.uniform_steps).
Counts from transpiling at optimization level 3 to {cx, rz, sx, x}. For each (route, schedule, layers, eps) the
J^2 step count k, the energy step count n_max and the number of blocks R minimise CNOTs per success subject to
1 - F <= eps at the chosen steps and at the next finer k and the next finer n_max. Exact fidelity to the ED ground
state, full 2^N space.

Reported per choice: CNOT and non-Clifford Rz per attempt and per success, total evolution time (energy stage in
scaled units of H_s = (H - E0)/(Emax - E0); J^2 stage in units of J^2, i.e. pi/4 per pulse), total Trotter steps.
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
from qpf import design, symmetry  # noqa: E402

ORDER = 2
EPS = (1e-2, 1e-3)
LAYERS = (1, 3)
K_J = [1, 2, 3, 4, 6, 8, 12]                          # Trotter steps of the pi/4 J^2 pulse
N_MAX = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]      # Trotter steps of the longest energy pulse
R_MAX = 5
M_RANGE = range(4, 13)
ROUTES = ("S2full -> H", "S2(pi/4) -> H", "H -> S2(pi/4)", "H only")
SCHEDS = ("QPF (t_i, δ_i)", "QPF (t_i)")


def dmrg_triplet(N, est, seed=0, bond_dims=(16, 32, 64, 128), projector_bond=16):
    lam = est["Emax"] - est["E0"]
    mpo = inputs._builder(N, 1.0, 0.2411).build_mpo().add_MPO(inputs._projector_mpo(est["psi0"], projector_bond) * (2 * lam))
    return inputs._dmrg(mpo, N, seed, bond_dims)[0]


def pick(gap, method):
    best = None
    for m in M_RANGE:
        q = design.optimize(gap, m, method=method)
        if method == "v3" and q["f0"] ** 2 < 0.85:              # infeasible SLSQP exit (all phases ~ pi/2)
            continue
        if best is None or q["cost"] < best["cost"] - 1e-9:
            best = dict(q, m=m)
    o = np.argsort(best["times"])
    return dict(times=np.asarray(best["times"])[o], phases=np.asarray(best["phases"])[o], m=best["m"],
                cost=best["cost"], f0=best["f0"])


class Stage:
    """A stage with a step-count grid: J^2 pulses (zero phase) or a designed energy block."""

    def __init__(self, kind, op, key, times, phases, grid, blocks=1):
        self.kind, self.op, self.key, self.times, self.phases = kind, op, key, np.asarray(times), np.asarray(phases)
        self.blocks, self.grid = blocks, []
        seen = set()
        for g in grid:                                         # drop grid points that give the same circuit
            ns, ts, dt = trotter.uniform_steps(self.times, g)
            if tuple(ns) not in seen:
                seen.add(tuple(ns))
                self.grid.append((g, ns, ts, dt))

    def apply(self, V, gi, idx):
        _, ns, ts, _ = self.grid[gi]
        return trotter.run_steps(self.op, idx, V, ts, self.phases, ns, ORDER)

    def counts(self, gi):
        _, ns, ts, _ = self.grid[gi]
        c = [trotter.step_counts(self.op, self.key, t, d, n, ORDER) for t, d, n in zip(ts, self.phases, ns)]
        return dict(cx=sum(x[0] for x in c), rz=sum(x[1] for x in c), steps=int(sum(ns)), ns=list(ns),
                    time=float(np.sum(np.abs(ts))))


def _norm(V):
    q = np.sum(abs(V) ** 2, axis=0)
    return V / np.sqrt(q), q


def evaluate(first, second, psi, g0, idx):
    """All (grid of first stage, block count if energy, grid of second stage) -> (infid, success).
    Returns dict[(gi, R, gj)] = (infid, P). For a J^2 stage R = 1; for an energy stage R = number of blocks."""
    out = {}
    for gi in range(len(first.grid)):
        V, P = psi[:, None], np.ones(1)
        outs = []
        for R in range(1, first.blocks + 1):
            V, q = _norm(first.apply(V, gi, idx))
            P = P * q
            outs.append((R, V[:, 0], float(P[0])))
        if second is None:
            for R, v, p in outs:
                out[(gi, R, None)] = (float(1 - abs(np.vdot(g0, v)) ** 2), p)
            continue
        W = np.stack([v for _, v, _ in outs], axis=1)
        P1 = np.array([p for *_, p in outs])
        conv_prev = None
        for gj in range(len(second.grid)):
            X, P2 = W, np.ones(W.shape[1])
            res = []
            for R2 in range(1, second.blocks + 1):
                X, q = _norm(second.apply(X, gj, idx))
                P2 = P2 * q
                res.append((R2, 1 - abs(g0.conj() @ X) ** 2, P1 * P2))
            for R2, inf, PP in res:
                for c, (R, _, _) in enumerate(outs):
                    out[(gi, R, gj, R2)] = (float(inf[c]), float(PP[c]))
    return out


def choose(first, second, table, trial_cx, trial_rz, eps):
    best = None
    nf, ns_ = len(first.grid), (len(second.grid) if second else 1)
    for key, (inf, P) in table.items():
        gi, R1 = key[0], key[1]
        gj, R2 = (key[2], key[3]) if second else (None, 1)
        if inf > eps:
            continue
        if second:
            nb = ([(gi + 1, R1, gj, R2)] if gi + 1 < nf else []) + ([(gi, R1, gj + 1, R2)] if gj + 1 < ns_ else [])
            if len(nb) < 2:
                continue
        else:
            nb = [(gi + 1, R1, None)] if gi + 1 < nf else []
            if not nb:
                continue
        ok = all(table.get(k, (1, 0))[0] <= eps for k in nb)
        if not ok:
            continue
        c1 = first.counts(gi)
        c2 = second.counts(gj) if second else dict(cx=0, rz=0, steps=0, time=0.0, ns=[])
        r1 = R1 if first.kind == "E" else 1
        r2 = R2 if second and second.kind == "E" else 1
        cx = trial_cx + r1 * c1["cx"] + r2 * c2["cx"]
        rz = trial_rz + r1 * c1["rz"] + r2 * c2["rz"]
        if best is None or cx / P < best["cx_per_success"]:
            tE = (r1 * c1["time"] if first.kind == "E" else 0) + (r2 * c2["time"] if second and second.kind == "E" else 0)
            tJ = (c1["time"] if first.kind == "J" else 0) + (c2["time"] if second and second.kind == "J" else 0)
            best = dict(infid=inf, success=P, cx_attempt=int(cx), rz_attempt=int(rz), cx_per_success=cx / P,
                        rz_per_success=rz / P, time_E=tE, time_J=tJ,
                        trotter_steps=int(r1 * c1["steps"] + r2 * c2["steps"]),
                        first=dict(kind=first.kind, grid=first.grid[gi][0], ns=c1["ns"], blocks=r1, cx=c1["cx"], rz=c1["rz"]),
                        second=None if not second else dict(kind=second.kind, grid=second.grid[gj][0], ns=c2["ns"],
                                                            blocks=r2, cx=c2["cx"], rz=c2["rz"]),
                        trial_cx=trial_cx, trial_rz=trial_rz)
    return best


def run_N(N, seed=0):
    t0 = time.time()
    ex = inputs.exact(N)
    est = inputs.dmrg_estimates(N, seed=seed)
    E1t = dmrg_triplet(N, est, seed)
    W = est["Emax"] - est["E0"]
    gaps = dict(singlet=(est["E1_singlet"] - est["E0"]) / W, triplet=(E1t - est["E0"]) / W)
    g0, idx = ex["g0"], np.arange(1 << N)
    Hop = trotter.scaled_hamiltonian(N, inputs.bonds(N), est["E0"], est["Emax"])
    Jop = trotter.total_spin_sq(N)
    designs = {(g, s): pick(gaps[g], m) for g in gaps for s, m in zip(SCHEDS, ("v3", "v4"))}
    full = symmetry.j2_ladder(0, N // 2)
    res = dict(N=N, gaps=gaps, gaps_exact=dict(singlet=(ex["E1_singlet"] - ex["E0"]) / (ex["Emax"] - ex["E0"]),
                                               triplet=(ex["E1"] - ex["E0"]) / (ex["Emax"] - ex["E0"])),
               designs={f"{g}|{s}": dict(m=d["m"], times=list(d["times"]), phases=list(d["phases"]), f0=d["f0"])
                        for (g, s), d in designs.items()},
               trials={}, routes={})
    for L in LAYERS:
        qc, psi = inputs.trial(est["psi0"], num_layers=L)
        tcx, trz = trotter.counts(qc)
        res["trials"][L] = dict(gamma=float(abs(np.vdot(g0, psi)) ** 2), cx=tcx, rz=trz)
        for route in ROUTES:
            for s in SCHEDS:
                gapk = "singlet" if route.startswith("S2") else "triplet"
                d = designs[(gapk, s)]
                E = Stage("E", Hop, ("H", N), d["times"], d["phases"], N_MAX, blocks=R_MAX)
                Jf = Stage("J", Jop, ("J", N), full, np.zeros(len(full)), K_J)
                J4 = Stage("J", Jop, ("J", N), [np.pi / 4], [0.0], K_J)
                first, second = {"S2full -> H": (Jf, E), "S2(pi/4) -> H": (J4, E),
                                 "H -> S2(pi/4)": (E, J4), "H only": (E, None)}[route]
                tab = evaluate(first, second, psi, g0, idx)
                for eps in EPS:
                    res["routes"][f"{L}|{route}|{s}|{eps}"] = choose(first, second, tab, tcx, trz, eps)
            print(f"  N={N} L={L} {route} [{time.time() - t0:.0f} s]", flush=True)
    res["seconds"] = time.time() - t0
    return res


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
    L = ["| N | ε | layers (γ) | route | schedule | steps: stage 1 → stage 2 (blocks) | success | CNOT / attempt "
         "| non-Clifford Rz / attempt | CNOT / success | Rz / success | time E (scaled) + J² | Trotter steps |",
         "|--|--|--|--|--|--|--|--|--|--|--|--|--|"]
    for r in res:
        for eps in EPS:
            for lay in LAYERS:
                for route in ROUTES:
                    for s in SCHEDS:
                        c = r["routes"].get(f"{lay}|{route}|{s}|{eps}")
                        g = r["trials"][str(lay)]["gamma"]
                        if c is None:
                            L.append(f"| {r['N']} | {eps:g} | {lay} ({g:.2f}) | {route} | {s} | not reached | | | | | | | |")
                            continue
                        st = lambda x: "" if x is None else f"{x['kind']}{tuple(x['ns'])}×{x['blocks']}"
                        L.append(f"| {r['N']} | {eps:g} | {lay} ({g:.2f}) | {route} | {s} | {st(c['first'])} → {st(c['second'])} "
                                 f"| {c['success']:.3f} | {c['cx_attempt']} | {c['rz_attempt']} | {c['cx_per_success']:.0f} "
                                 f"| {c['rz_per_success']:.0f} | {c['time_E']:.1f} + {c['time_J']:.2f} | {c['trotter_steps']} |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12, 14, 16])
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--out", default="j1j2_routes")
    a = ap.parse_args()
    with ProcessPoolExecutor(a.workers) as ex:
        res = _clean(list(ex.map(run_N, a.N)))
    json.dump(res, open(os.path.join(HERE, "results", a.out + ".json"), "w"), indent=1)
    open(os.path.join(HERE, "results", a.out + ".md"), "w").write(table(res) + "\n")
    print(table(res))


if __name__ == "__main__":
    main()
