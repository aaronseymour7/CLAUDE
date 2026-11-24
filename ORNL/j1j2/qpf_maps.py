"""Gap-vs-overlap cost maps with the published QPF (J^2 ladder, then the authors' designed energy pulses).

    python j1j2/qpf_maps.py --N 8 --out DIR            # one N, every (delta, knob, gamma, eps); writes DIR/qpf_N8.json
    python j1j2/qpf_maps.py --N 8 --gaps 0.5 1.0 --gammas 0.1 0.5

Setup (manuscript Section A, decisions D1-D5):
  chain    H = J1 sum_i (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, open, J1 = 1, J2 = 0.2411, S^z = 0 sector
  gap      D5: each stage filters its own operator, and its gap is the spacing of that operator's two lowest
           eigenvalues in the space the stage acts on. J^2 stage: 0 -> 2. Energy stage (after J^2): the singlet gap
           E1(S=0) - E0. delta is solved per N so that the singlet gap hits each target (branch right of the
           singlet-gap minimum). The global gap E1 - E0 (lowest triplet where it is below the singlet) is stored too.
  overlap  trial = knob applied to the exact ground state, gamma = |<E0|trial>|^2 (D4):
             staggered  exp(-i e sum_i (-1)^i Z_i)                (breaks S^2; the cost_study knob)
             singlet    prod_i exp(-+i e SWAP_{i,i+1}) brick layer (keeps S^2)
           e solved for each target gamma on the first descent.
  routes   "J2 -> H": the published QPF pipeline below. "H only": the energy stage alone, designed for the gap of H in
           the space the trial lives in (D5): the global gap E1 - E0 for the staggered knob, the singlet gap for the
           singlet knob (an S^2-preserving trial never sees a triplet).
  QPF      stage 1: J^2 halving ladder pi/4, pi/8, ... to J = 0 (qpf.symmetry.j2_ladder), zero phases.
           stage 2: energy pulses from the authors' design code (qpf.design.optimize) for the scaled singlet gap
           Delta_s / (Emax - E0); both published variants, (t_i, delta_i) = "v3" and (t_i) = "v4", m = 4..12 chosen by
           design cost (j1j2.route_costs.pick). The block may be repeated R times.
  circuit  second-order Trotter, one step size per stage set by the longest pulse (j1j2.trotter.uniform_steps);
           gate counts from transpiling each pulse at optimization level 3 to {cx, rz, sx, x}
           (j1j2.trotter.step_counts), plus the knob layer. The exact ground-state preparation is not counted.
           Also stored: the cost_study analytic model (6 CNOT + 3 Rz per controlled term, (2B - 2) n + 1 blocks for
           n merged Strang steps of B terms, + 1 Rz for a non-Clifford phase), about 0.7x the transpile (Qiskit
           synthesises the controlled bond with 8 CNOT). Every row carries gate_model; the (k_J, n_E, R) search is done
           separately for each model, and the other model's per-attempt count for the same circuit is kept in
           *_other_model. AA reflections are transpiled in both.
           Merges: the transpile is per pulse, and its CNOT count does not depend on the angles (checked), so it relies
           on no angle-specific simplification; the Rz count is the larger of the counts at the actual and at generic
           angles (j1j2.qpf_maps.robust_counts). The analytic model merges only the end terms of consecutive Strang
           half-steps, which are the same Pauli term (exp(-iaP) exp(-ibP) = exp(-i(a+b)P) for any a, b).
  success  1 - F <= eps at the chosen (k_J, n_E, R) and also at the next finer k_J and the next finer n_E (as
           route_costs). Measured costs only (D3).
  RUS      CNOT per success = C_attempt / p.
  AA       one ancilla per pulse (no mid-circuit measurement), good = all ancillas 0, k = round(pi / 4 theta - 1/2),
           sin^2 theta = p; cost = ((2k + 1) C_A + k C_R) / sin^2((2k + 1) theta), C_R = mcz(n_anc) + mcz(N + n_anc)
           from Qiskit's mcx. The (k_J, n_E, R) choice is made separately for RUS and AA.

Per (N, delta, knob, gamma, eps, schedule, protocol) the output stores the pulse list (times, phases, Trotter steps),
the step sizes dt, p, total evolution time (energy stage in units of 1/J1 and scaled units, J^2 stage separately),
Trotter steps, CNOT and non-Clifford Rz per attempt and per success, and the AA round count.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import eigsh

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import trotter  # noqa: E402
from j1j2.route_costs import pick  # noqa: E402
from qpf import symmetry  # noqa: E402

J1, J2 = 1.0, 0.2411
ORDER = 2
EPS = (1e-2, 1e-3)
K_J = [1, 2, 3, 4, 6, 8, 12]                            # Trotter steps of the pi/4 J^2 pulse
N_E = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]       # Trotter steps of the longest energy pulse
R_MAX = 4
SCHEDS = {"QPF (t_i, δ_i)": "v3", "QPF (t_i)": "v4"}
GAPS = (0.55, 0.8, 1.15, 1.65, 2.35, 3.2)               # singlet gap targets, units of J1 (coarse grid)
GAMMAS = (0.01, 0.03, 0.1, 0.3, 0.6, 0.9)
KNOBS = ("staggered", "singlet")


# ---------------------------------------------------------------- model in the S^z = 0 sector
def bonds(N, delta):
    return [(i, i + 1, J1 * (1 + delta * (-1) ** i)) for i in range(N - 1)] + [(i, i + 2, J2) for i in range(N - 2)]


class SectorSwapSum(trotter.SwapSum):
    """SwapSum whose permutations on a sorted basis subset are found with searchsorted (fast for N = 16)."""

    def perms(self, index):
        key = (len(index), int(index[0]), int(index[-1]))
        if key not in self._perm:
            out = []
            for i, j in self.pairs:
                flip = ((index >> i) & 1) != ((index >> j) & 1)
                t = np.where(flip, index ^ (1 << i) ^ (1 << j), index)
                out.append(np.searchsorted(index, t))
            self._perm[key] = out
        return self._perm[key]


def as_sector(op):
    return SectorSwapSum(op.N, op.pairs, op.weights, op.const)


def sector(N):
    s = np.arange(1 << N)
    pop = np.array([bin(x).count("1") for x in s])
    return s[pop == N // 2]


def matrix(op, idx):
    """Sparse matrix of a SwapSum on the basis idx."""
    d = len(idx)
    P = op.perms(idx)
    rows = np.concatenate(P)
    cols = np.tile(np.arange(d), len(P))
    vals = np.concatenate([np.full(d, w) for w in op.weights])
    return (sp.csr_matrix((vals, (rows, cols)), shape=(d, d)) + op.const * sp.identity(d, format="csr")).tocsr()


class Chain:
    def __init__(self, N, delta, idx=None, S2=None):
        self.N, self.delta = N, delta
        self.idx = sector(N) if idx is None else idx
        self.Hop = as_sector(trotter.scaled_hamiltonian(N, bonds(N, delta), 0.0, 1.0))      # H itself (W = 1)
        self.H = matrix(self.Hop, self.idx)
        self.S2 = matrix(as_sector(trotter.total_spin_sq(N)), self.idx) if S2 is None else S2

    def spectrum(self):
        w, V = eigsh(self.H, k=4, which="SA", tol=1e-12)
        o = np.argsort(w)
        w, V = w[o], V[:, o]
        self.E0, self.g0 = float(w[0]), V[:, 0].astype(complex)
        self.E1 = float(w[w > w[0] + 1e-8][0])
        self.Emax = float(eigsh(self.H, k=1, which="LA", tol=1e-10, return_eigenvectors=False)[0])
        self.W = self.Emax - self.E0
        self.gap_global = self.E1 - self.E0
        self.gap_singlet = self.singlet_gap()
        return self

    def singlet_gap(self, E0=None):
        E0 = self.E0 if E0 is None else E0
        W = 2 * (abs(E0) + 1)
        ws = np.sort(eigsh(self.H + W * self.S2, k=2, which="SA", tol=1e-12, return_eigenvectors=False))
        return float(ws[1] - ws[0])


def singlet_gap_of(N, delta, idx, S2):
    c = Chain(N, delta, idx, S2)
    w = eigsh(c.H, k=1, which="SA", tol=1e-12, return_eigenvectors=False)[0]
    return c.singlet_gap(w)


def solve_deltas(N, targets, idx, S2, lo=-0.3, hi=0.85, scan=24):
    """delta for each singlet-gap target on the branch right of the singlet-gap minimum."""
    ds = np.linspace(lo, hi, scan)
    g = np.array([singlet_gap_of(N, d, idx, S2) for d in ds])
    i0 = int(np.argmin(g))
    out = {}
    for t in targets:
        j = [k for k in range(i0, len(ds) - 1) if (g[k] - t) * (g[k + 1] - t) <= 0]
        if not j:
            out[t] = None
            continue
        k = j[0]
        out[t] = float(brentq(lambda d: singlet_gap_of(N, d, idx, S2) - t, ds[k], ds[k + 1], xtol=1e-7))
    return out, dict(delta=ds.tolist(), gap_singlet=g.tolist())


# ---------------------------------------------------------------- knobs
def staggered_diag(N, idx):
    return sum((-1) ** i * (1 - 2 * ((idx >> i) & 1)) for i in range(N)).astype(float)


def singlet_layer(N, idx, v, e):
    for i in range(N - 1):
        t = np.where(((idx >> i) & 1) != ((idx >> (i + 1)) & 1), idx ^ (1 << i) ^ (1 << (i + 1)), idx)
        p = np.searchsorted(idx, t)
        th = e if i % 2 == 0 else -e
        v = np.cos(th) * v - 1j * np.sin(th) * v[p]
    return v


def knob_fn(name, N, idx, g0):
    if name == "staggered":
        Z = staggered_diag(N, idx)
        return (lambda e: np.exp(-1j * e * Z) * g0), 1.5
    return (lambda e: singlet_layer(N, idx, g0, e)), np.pi / 2


def eps_for(f, g0, gamma, hi, n=600):
    es = np.linspace(0, hi, n)
    g = np.array([abs(np.vdot(g0, f(e))) ** 2 for e in es])
    hit = np.flatnonzero(g <= gamma)
    if not len(hit):
        return None
    i = hit[0]
    return float(brentq(lambda e: abs(np.vdot(g0, f(e))) ** 2 - gamma, es[i - 1], es[i], xtol=1e-13))


_KNOB_COUNTS = {}


def knob_counts(name, N, e):
    """(CNOT, non-Clifford Rz) of the knob layer (transpiled like the pulses)."""
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    key = (name, N, round(e, 12))
    if key not in _KNOB_COUNTS:
        qc = QuantumCircuit(N)
        if name == "staggered":
            for i in range(N):
                qc.rz(2 * e * (-1) ** i, i)
        else:
            for i in range(N - 1):        # exp(-i th SWAP) = e^{-i th/2} exp(-i th/2 (XX + YY + ZZ))
                th = e if i % 2 == 0 else -e
                A = SparsePauliOp.from_sparse_list([(p * 2, [i, i + 1], 0.5) for p in "XYZ"], num_qubits=N)
                qc.append(PauliEvolutionGate(A, time=th), range(N))
        _KNOB_COUNTS[key] = trotter.counts(qc)
    return _KNOB_COUNTS[key]


# ---------------------------------------------------------------- AA reflections
_MCZ = {}


def mcz_counts(nq):
    if nq not in _MCZ:
        from qiskit import QuantumCircuit
        qc = QuantumCircuit(nq)
        qc.x(range(nq))
        if nq == 1:
            qc.z(0)
        else:
            qc.h(nq - 1); qc.mcx(list(range(nq - 1)), nq - 1); qc.h(nq - 1)
        qc.x(range(nq))
        _MCZ[nq] = trotter.counts(qc)
    return _MCZ[nq]


def aa(p, cA, rzA, N, n_anc):
    """Cheapest number of Grover rounds k in 0..round(pi / 4 theta - 1/2) (k = 0 is plain RUS on the unitary
    circuit, so AA is used only where it beats RUS)."""
    th = np.arcsin(np.sqrt(min(p, 1.0)))
    k_opt = max(0, int(round(np.pi / (4 * th) - 0.5)))
    rc = np.add(mcz_counts(n_anc), mcz_counts(N + n_anc))
    best = None
    for k in range(k_opt + 1):
        P = np.sin((2 * k + 1) * th) ** 2
        cx = ((2 * k + 1) * cA + k * rc[0]) / P
        if best is None or cx < best["cx"]:
            best = dict(k=k, k_opt=k_opt, P_aa=float(P), cx=cx, rz=((2 * k + 1) * rzA + k * rc[1]) / P,
                        cx_reflection=int(rc[0]), rz_reflection=int(rc[1]))
    return best


# ---------------------------------------------------------------- stages
class Stage:
    def __init__(self, op, key, times, phases, grid):
        self.op, self.key = op, key
        self.times, self.phases = np.asarray(times, float), np.asarray(phases, float)
        self.grid, seen = [], set()
        for g in grid:
            ns, ts, dt = trotter.uniform_steps(self.times, g)
            if tuple(ns) not in seen:
                seen.add(tuple(ns))
                self.grid.append(dict(n_longest=g, ns=list(ns), times=ts, dt=dt))
        self._c = {}

    def apply(self, V, gi, idx):
        g = self.grid[gi]
        return trotter.run_steps(self.op, idx, V, g["times"], self.phases, g["ns"], ORDER)

    def counts(self, gi):
        if gi not in self._c:
            g = self.grid[gi]
            c = [robust_counts(self.op, self.key, t, d, n) for t, d, n in zip(g["times"], self.phases, g["ns"])]
            B = len(self.op.pairs)
            blocks = [(2 * B - 2) * n + 1 for n in g["ns"]]               # analytic model: merged Strang half-steps
            rz_ph = [int(abs(d / (np.pi / 2) - round(d / (np.pi / 2))) > 1e-9) for d in self.phases]
            self._c[gi] = dict(cx=sum(x[0] for x in c), rz=sum(x[1] for x in c), steps=int(sum(g["ns"])),
                               cx_an=6 * sum(blocks), rz_an=3 * sum(blocks) + sum(rz_ph),
                               time=float(np.sum(np.abs(g["times"]))), per_pulse=[list(x) for x in c])
        return self._c[gi]


def robust_counts(op, key, t, d, n):
    """Transpiled (CNOT, non-Clifford Rz) of one pulse that do not rely on a lucky angle: the larger of the counts at
    the actual angles and at generic nearby angles (t (1 + sqrt(2) 1e-3), delta + sqrt(3) 1e-3 unless delta = 0, which
    is no gate at all). The CNOT count is angle independent (checked); the Rz count moves by about 1%."""
    a = trotter.step_counts(op, key, t, d, n, ORDER)
    b = trotter.step_counts(op, key, t * (1 + np.sqrt(2) * 1e-3), d + (np.sqrt(3) * 1e-3 if d != 0 else 0.0), n, ORDER)
    return max(a[0], b[0]), max(a[1], b[1])


def _norm(V):
    q = np.sum(abs(V) ** 2, axis=0)
    return V / np.sqrt(q), q


def run_point(N, delta, gammas, knobs, idx, S2, log=print):
    t0 = time.time()
    ch = Chain(N, delta, idx, S2).spectrum()
    W = ch.W
    gs_scaled = ch.gap_singlet / W
    res = dict(N=N, delta=delta, E0=ch.E0, Emax=ch.Emax, W=W, gap_singlet=ch.gap_singlet, gap_global=ch.gap_global,
               gap_singlet_scaled=gs_scaled, gap_global_scaled=ch.gap_global / W, gap_J2=2.0, trials=[], designs={},
               jstage={}, results=[])
    # trials (columns)
    cols, meta = [], []
    for kn in knobs:
        f, hi = knob_fn(kn, N, idx, ch.g0)
        for gt in gammas:
            e = eps_for(f, ch.g0, gt, hi)
            if e is None:
                res["trials"].append(dict(knob=kn, gamma_target=gt, reached=False))
                continue
            v = f(e)
            g = float(abs(np.vdot(ch.g0, v)) ** 2)
            w = np.real(np.vdot(v, S2 @ v))
            kc = knob_counts(kn, N, e)
            meta.append(dict(knob=kn, gamma_target=gt, gamma=g, eps_knob=e, S2_expect=float(w), knob_cx=kc[0],
                             knob_rz=kc[1]))
            res["trials"].append(dict(meta[-1], reached=True))
            cols.append(v)
    if not cols:
        return res
    V0 = np.stack(cols, axis=1)
    Hs = SectorSwapSum(N, ch.Hop.pairs, ch.Hop.weights / W, (ch.Hop.const - ch.E0) / W)
    Jop = as_sector(trotter.total_spin_sq(N))
    jt = symmetry.j2_ladder(0, N // 2)
    J = Stage(Jop, ("J", N), jt, np.zeros(len(jt)), K_J)
    res["jstage"] = dict(times=list(map(float, jt)), grid=[dict(g, counts=J.counts(i)) for i, g in enumerate(J.grid)])
    VJ = [_norm(J.apply(V0, gi, idx)) for gi in range(len(J.grid))]
    log(f"  N={N} delta={delta:+.4f} J^2 stage [{time.time() - t0:.0f} s]")
    # routes: (name, uses J^2 stage, trial columns, energy-stage gap per D5)
    routes = [("J2 -> H", True, list(range(len(meta))), "singlet")]
    for kn in knobs:
        c_kn = [c for c, m in enumerate(meta) if m["knob"] == kn]
        if c_kn:
            routes.append(("H only", False, c_kn, "singlet" if kn == "singlet" else "global"))
    designs = {}
    for route, useJ, cs, gapk in routes:
        gap_s = (ch.gap_singlet if gapk == "singlet" else ch.gap_global) / W
        for sname, method in SCHEDS.items():
            key = (gapk, sname)
            if key not in designs:
                d = pick(gap_s, method)
                designs[key] = d
                res["designs"][f"{gapk}|{sname}"] = dict(gap=gapk, gap_scaled=gap_s, m=int(d["m"]),
                                                         times=d["times"].tolist(), phases=d["phases"].tolist(),
                                                         design_cost=float(d["cost"]), f0=float(d["f0"]))
            d = designs[key]
            E = Stage(Hs, ("H", N, round(delta, 10)), d["times"], d["phases"], N_E)
            if useJ:
                ins = [(x[:, cs], q[cs]) for x, q in VJ]
            else:
                ins = [(V0[:, cs], np.ones(len(cs)))]
            tab = evaluate_energy(E, ins, ch.g0, idx)
            nJ = len(ins)
            for c_local, c in enumerate(cs):
                m = meta[c]
                for eps in EPS:
                    MODELS = ("transpile", "analytic")
                    best = {(gm, pr): None for gm in MODELS for pr in ("rus", "aa")}
                    for (gi, gj, R), (inf, P) in tab.items():
                        if not valid(tab, gi, gj, R, c_local, eps, nJ, useJ):
                            continue
                        cJ = J.counts(gi) if useJ else dict(cx=0, rz=0, cx_an=0, rz_an=0, steps=0, time=0.0)
                        cE = E.counts(gj)
                        p = float(P[c_local])
                        n_anc = (len(jt) if useJ else 0) + R * len(d["times"])
                        att = {"transpile": (m["knob_cx"] + cJ["cx"] + R * cE["cx"], m["knob_rz"] + cJ["rz"] + R * cE["rz"]),
                               "analytic": (m["knob_cx"] + cJ["cx_an"] + R * cE["cx_an"],
                                            m["knob_rz"] + cJ["rz_an"] + R * cE["rz_an"])}
                        for gm in MODELS:
                            cx, rz = att[gm]
                            other = att["analytic" if gm == "transpile" else "transpile"]
                            cand = {"rus": dict(cx=cx / p, rz=rz / p), "aa": aa(p, cx, rz, N, n_anc)}
                            for proto in ("rus", "aa"):
                                bk = (gm, proto)
                                if best[bk] is not None and cand[proto]["cx"] >= best[bk]["cx_per_success"]:
                                    continue
                                g = E.grid[gj]
                                best[bk] = dict(
                                    k_J=J.grid[gi]["n_longest"] if useJ else 0, ns_J=J.grid[gi]["ns"] if useJ else [],
                                    dt_J=J.grid[gi]["dt"] if useJ else None,
                                    n_E=g["n_longest"], ns_E=g["ns"], dt_E_scaled=g["dt"], dt_E=g["dt"] / W, R=R,
                                    m_E=len(d["times"]), times_E_scaled=g["times"], phases_E=d["phases"].tolist(),
                                    infid=float(inf[c_local]), p=p, p_J=float(ins[gi][1][c_local]),
                                    cx_attempt=int(cx), rz_attempt=int(rz),
                                    cx_attempt_other_model=int(other[0]), rz_attempt_other_model=int(other[1]),
                                    cx_J=cJ["cx" if gm == "transpile" else "cx_an"],
                                    rz_J=cJ["rz" if gm == "transpile" else "rz_an"],
                                    cx_E_block=cE["cx" if gm == "transpile" else "cx_an"],
                                    rz_E_block=cE["rz" if gm == "transpile" else "rz_an"],
                                    time_E_scaled=R * cE["time"], time_E=R * cE["time"] / W, time_J=cJ["time"],
                                    trotter_steps_J=cJ["steps"], trotter_steps_E=R * cE["steps"],
                                    trotter_steps=cJ["steps"] + R * cE["steps"], n_anc_aa=n_anc,
                                    cx_per_success=cand[proto]["cx"], rz_per_success=cand[proto]["rz"],
                                    **({k: cand["aa"][k] for k in ("k", "k_opt", "P_aa", "cx_reflection", "rz_reflection")}
                                       if proto == "aa" else {}))
                    for (gm, proto), bb in best.items():
                        res["results"].append(dict(m, eps=eps, route=route, gap_used=gapk, gap_used_value=gap_s * W,
                                                   schedule=sname, gate_model=gm, protocol=proto, reached=bb is not None,
                                                   **(bb or {})))
            log(f"  N={N} delta={delta:+.4f} {route} [{gapk}] {sname} [{time.time() - t0:.0f} s]")
    res["seconds"] = time.time() - t0
    return res


def valid(tab, gi, gj, R, c, eps, nJ, useJ):
    """1 - F <= eps at (gi, gj, R) and at the next finer J^2 grid (if a J^2 stage is used) and next finer energy grid."""
    k = (gi, gj, R)
    if k not in tab or tab[k][0][c] > eps:
        return False
    nb = [(gi, gj + 1, R)] + ([(gi + 1, gj, R)] if useJ else [])
    return all(x in tab and tab[x][0][c] <= eps for x in nb)


def evaluate_energy(E, ins, g0, idx):
    """tab[(gi, gj, R)] = (1 - F, P) per column, for every first-stage output gi (stacked as columns), energy grid gj
    and block count R. Stops refining the energy grid once every (column, eps) has had a valid choice at two
    consecutive grid levels (finer steps only cost more), and stops adding blocks once every column is below min(EPS)."""
    nc = ins[0][0].shape[1]
    nJ = len(ins)
    X0 = np.concatenate([x for x, _ in ins], axis=1)
    P0 = np.concatenate([q for _, q in ins])
    tab, streak = {}, 0
    useJ = nJ > 1
    for gj in range(len(E.grid)):
        X, P = X0, P0.copy()
        for R in range(1, R_MAX + 1):
            X, q = _norm(E.apply(X, gj, idx))
            P = P * q
            inf = 1 - abs(g0.conj() @ X) ** 2
            for gi in range(nJ):
                sl = slice(gi * nc, (gi + 1) * nc)
                tab[(gi, gj, R)] = (inf[sl], P[sl])
            if np.all(inf <= min(EPS)):
                break
        if gj > 0:
            ok = all(any(valid(tab, gi, gj - 1, R, c, e, nJ, useJ) for gi in range(nJ) for R in range(1, R_MAX + 1))
                     for c in range(nc) for e in EPS)
            streak = streak + 1 if ok else 0
            if streak >= 2:
                break
    return tab


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    if isinstance(x, np.generic):
        return x.item()
    return x


def grid(N, gaps, out):
    """delta for each singlet-gap target (cached in out/grid_N{N}.json)."""
    path = os.path.join(out, f"grid_N{N}.json")
    g = json.load(open(path)) if os.path.exists(path) else dict(N=N, deltas={})
    todo = [t for t in gaps if f"{t:g}" not in g["deltas"]]
    if todo:
        idx = sector(N)
        S2 = matrix(as_sector(trotter.total_spin_sq(N)), idx)
        d, scan = solve_deltas(N, todo, idx, S2)
        g["deltas"].update({f"{t:g}": v for t, v in d.items()})
        g["gap_scan"] = scan
        json.dump(_clean(g), open(path, "w"), indent=1)
    return {t: g["deltas"][f"{t:g}"] for t in gaps}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--gaps", type=float, nargs="+", default=list(GAPS))
    ap.add_argument("--gammas", type=float, nargs="+", default=list(GAMMAS))
    ap.add_argument("--knobs", nargs="+", default=list(KNOBS))
    ap.add_argument("--out", default=os.path.join(HERE, "results", "qpf_maps"))
    ap.add_argument("--grid-only", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    deltas = grid(a.N, a.gaps, a.out)
    print(f"N={a.N}: deltas {deltas} [{time.time() - t0:.0f} s]", flush=True)
    if a.grid_only:
        return
    idx = sector(a.N)
    S2 = matrix(as_sector(trotter.total_spin_sq(a.N)), idx)
    for gt, dl in deltas.items():
        path = os.path.join(a.out, f"qpf_N{a.N}_gap{gt:g}.json")
        if dl is None or os.path.exists(path):
            continue
        r = run_point(a.N, dl, a.gammas, a.knobs, idx, S2, log=lambda s: print(s, flush=True))
        r.update(gap_target=gt, J2=J2, eps_list=list(EPS), K_J=K_J, N_E=N_E, R_MAX=R_MAX)
        json.dump(_clean(r), open(path + ".tmp", "w"))
        os.replace(path + ".tmp", path)
    print(f"N={a.N} done [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
