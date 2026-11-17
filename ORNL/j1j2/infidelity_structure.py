"""Gap-overlap study, one layer down: where the trial's infidelity sits, in S^2 and in H, and over how many states.

    python j1j2/infidelity_structure.py                 # N = 8..14 -> results/j1j2_infidelity_structure.{json,md}
    python j1j2/infidelity_structure.py --N 8 10

Chain: H = sum_i J1 (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2} (open, J2 = 0.2411); the gap is tuned by
the dimerisation delta, as in the gap-overlap study. Trials (all in the full 2^N space, Qiskit order, site i = bit i):

  mps1          the one-layer approximate mps-to-circuit circuit for the ground state (the MPS is the exact ground
                state written as an MPS; the circuit is what loses the overlap)
  staggered Rz  the study's overlap knob on top of it: exp(-i eps sum_i (-1)^i Z_i) |mps1>, eps set so that
                gamma = |<E0|trial>|^2 hits each target (first descent). Conserves S^z, breaks S^2.
  singlet bond  the SU(2)-preserving knob: one brick layer prod e^{-+i eps SWAP_{i,i+1}} on |mps1>.

H and S^2 commute and both conserve S^z, so in each S^z sector the common eigenbasis is found by diagonalising
H + a S^2 densely (a small and irrational to split degeneracies). A "level" is one (S, E) multiplet, i.e. its 2S+1
S^z members merged. For the infidelity 1 - gamma (all weight off the ground state) we report:
  by S        weight in each total spin S; number of S values holding 90% / 99% of 1 - gamma
  by level    number of H levels (multiplets) holding 90% / 99% of 1 - gamma, the participation ratio
              (sum w)^2 / sum w^2 of the level weights, the weight on the lowest triplet and the first excited
              singlet, and the same count restricted to singlets (what the energy filter still sees after J^2).
"""
import argparse
import json
import os
import sys
from itertools import combinations  # noqa: F401

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import inputs  # noqa: E402

J1, J2 = 1.0, 0.2411
DELTAS = (0.0, 0.2, 0.4, 0.6, 0.8)
GAMMAS = (0.8, 0.5, 0.2, 0.05)
A = 1e-3 * np.pi
TOL = 1e-7


def bonds(N, delta):
    return [(i, i + 1, J1 * (1 + delta * (-1) ** i)) for i in range(N - 1)] + [(i, i + 2, J2) for i in range(N - 2)]


def hamiltonian(N, delta):
    b = bonds(N, delta)
    H = inputs._swap_sum(N, [(i, j) for i, j, _ in b], [J / 2 for *_, J in b])
    return (H - 0.25 * sum(J for *_, J in b) * sp.identity(1 << N, format="csr")).tocsr()


class Spectrum:
    """Common eigenbasis of H and S^2, sector by sector."""

    def __init__(self, N, delta):
        self.N = N
        H, S2 = hamiltonian(N, delta), inputs.total_spin_sq(N)
        pop = np.array([bin(i).count("1") for i in range(1 << N)])
        self.sectors = []
        for nup in range(N + 1):
            idx = np.flatnonzero(pop == nup)
            h, s2 = H[idx][:, idx].toarray(), S2[idx][:, idx].toarray()
            _, V = np.linalg.eigh(h + A * s2)
            E = np.einsum("ij,ij->j", V.conj(), h @ V).real
            J = np.einsum("ij,ij->j", V.conj(), s2 @ V).real
            S = np.abs(np.round(-0.5 + np.sqrt(0.25 + np.maximum(J, 0)))) + 0.0
            self.sectors.append((idx, V, E, S))
        allE = np.concatenate([s[2] for s in self.sectors])
        allS = np.concatenate([s[3] for s in self.sectors])
        self.E0 = float(allE.min())
        k = np.argmin(allE)
        for idx, V, E, S in self.sectors:                           # ground state (S^z = 0 member)
            j = np.flatnonzero(abs(E - self.E0) < TOL)
            if len(j) and len(idx) == len(np.flatnonzero(pop == N // 2)):
                self.g0 = np.zeros(1 << N, complex)
                self.g0[idx] = V[:, j[0]]
        self.S0 = float(allS[k])
        ts = np.sort(allE[(allS == 1)])
        ss = np.unique(np.round(allE[(allS == 0)], 9))
        self.gap_triplet, self.gap_singlet = float(ts[0] - self.E0), float(ss[1] - self.E0)

    def levels(self, psi):
        """[(S, E, weight)] per multiplet."""
        rows = []
        for idx, V, E, S in self.sectors:
            w = abs(V.conj().T @ psi[idx]) ** 2
            rows.append(np.stack([S, E, w], 1))
        r = np.concatenate(rows)
        r = r[np.lexsort((r[:, 1], r[:, 0]))]
        out = []
        for s, e, w in r:
            if out and out[-1][0] == s and abs(out[-1][1] - e) < 1e-6:
                out[-1][2] += w
            else:
                out.append([s, e, w])
        return np.array(out)


def count_for(w, frac):
    w = np.sort(w)[::-1]
    c = np.cumsum(w)
    return int(np.searchsorted(c, frac * c[-1] - 1e-15) + 1) if c[-1] > 0 else 0


def analyse(sp_, psi):
    L = sp_.levels(psi)
    ground = (L[:, 0] == 0) & (abs(L[:, 1] - sp_.E0) < 1e-6)
    gamma = float(L[ground, 2].sum())
    ex = L[~ground]
    w = ex[:, 2]
    infid = float(w.sum())
    byS = {}
    for s in np.unique(ex[:, 0]):
        byS[float(s)] = float(w[ex[:, 0] == s].sum())
    Sw = np.array(list(byS.values()))
    sing = ex[ex[:, 0] == 0]
    trip = ex[ex[:, 0] == 1]
    out = dict(gamma=gamma, infid=infid,
               by_S={f"{k:g}": v for k, v in byS.items() if v > 1e-8},
               nS90=count_for(Sw, 0.9), nS99=count_for(Sw, 0.99),
               n90=count_for(w, 0.9), n99=count_for(w, 0.99), PR=float(infid ** 2 / np.sum(w ** 2)),
               n_levels_total=int(len(L)),
               w_first_singlet=float(sing[np.argmin(sing[:, 1]), 2]) if len(sing) else 0.0,
               w_lowest_triplet=float(trip[np.argmin(trip[:, 1]), 2]) if len(trip) else 0.0,
               singlet_infid=float(sing[:, 2].sum()),
               singlet_n90=count_for(sing[:, 2], 0.9), singlet_n99=count_for(sing[:, 2], 0.99),
               singlet_PR=float(sing[:, 2].sum() ** 2 / np.sum(sing[:, 2] ** 2)) if sing[:, 2].sum() > 0 else 0.0,
               top=[dict(S=float(s), dE=float(e - sp_.E0), w=float(x)) for s, e, x in ex[np.argsort(-w)][:5]])
    return out


def mps1(g0, N):
    import quimb.tensor as qtn
    from mps_to_circuit import mps_to_circuit
    from qiskit.quantum_info import Statevector
    v = g0.reshape([2] * N).transpose(list(range(N))[::-1]).ravel()          # quimb order: site 0 most significant
    p = qtn.MatrixProductState.from_dense(v, dims=[2] * N)
    p.permute_arrays("lrp")
    qc = mps_to_circuit(list(p.arrays), method="approximate", shape="lrp", num_layers=1)
    return np.asarray(Statevector(qc).data)


def staggered(N):
    s = np.arange(1 << N)
    return sum((-1) ** i * (1 - 2 * ((s >> i) & 1)) for i in range(N)).astype(float)    # sum_i (-1)^i Z_i


def singlet_layer(N, v, eps):
    for i in range(N - 1):
        p = inputs._swap(N, i, i + 1)
        th = eps if i % 2 == 0 else -eps
        v = np.cos(th) * v - 1j * np.sin(th) * v[p]
    return v


def eps_for(f, g0, gamma, hi):
    from scipy.optimize import brentq
    es = np.linspace(0, hi, 1501)
    g = np.array([abs(np.vdot(g0, f(e))) ** 2 for e in es])
    hit = np.flatnonzero(g <= gamma)
    if not len(hit) or not g[0] >= gamma:
        return None
    i = hit[0]
    return float(brentq(lambda e: abs(np.vdot(g0, f(e))) ** 2 - gamma, es[i - 1], es[i], xtol=1e-13))


def run(N):
    res = []
    M = staggered(N)
    for d in DELTAS:
        S = Spectrum(N, d)
        base = mps1(S.g0, N)
        row = dict(N=N, delta=d, gap_triplet=S.gap_triplet, gap_singlet=S.gap_singlet, trials=[])
        row["trials"].append(dict(trial="mps1", eps=0.0, **analyse(S, base)))
        for knob, f, hi in (("staggered Rz", lambda e: np.exp(-1j * e * M) * base, 1.5),
                            ("singlet bond", lambda e: singlet_layer(N, base, e), np.pi / 2)):
            for gt in GAMMAS:
                e = eps_for(f, S.g0, gt, hi)
                if e is None:
                    continue
                row["trials"].append(dict(trial=knob, eps=e, **analyse(S, f(e))))
        res.append(row)
        t0 = row["trials"][0]
        print(f"N={N} delta={d}: gaps T {S.gap_triplet:.3f} S {S.gap_singlet:.3f}; mps1 gamma {t0['gamma']:.3f} "
              f"n90 {t0['n90']} PR {t0['PR']:.1f} by_S {t0['by_S']}", flush=True)
    return res


def table(res):
    L = ["| N | δ | gap T, S | trial | γ | 1−γ by S (S=0 exc / 1 / 2 / ≥3) | # S for 90%, 99% | # levels for 90%, 99% "
         "| PR | 1st exc. singlet, lowest triplet | singlet part: # levels 90%, 99%, PR |",
         "|--|--|--|--|--|--|--|--|--|--|--|"]
    for row in res:
        for t in row["trials"]:
            b = {float(k): v for k, v in t["by_S"].items()}
            s3 = sum(v for k, v in b.items() if k >= 3)
            L.append(f"| {row['N']} | {row['delta']:g} | {row['gap_triplet']:.3f}, {row['gap_singlet']:.3f} | {t['trial']} "
                     f"| {t['gamma']:.3f} | {b.get(0.0, 0):.3f} / {b.get(1.0, 0):.3f} / {b.get(2.0, 0):.3f} / {s3:.3f} "
                     f"| {t['nS90']}, {t['nS99']} | {t['n90']}, {t['n99']} | {t['PR']:.1f} "
                     f"| {t['w_first_singlet']:.3f}, {t['w_lowest_triplet']:.3f} "
                     f"| {t['singlet_n90']}, {t['singlet_n99']}, {t['singlet_PR']:.1f} |")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12, 14])
    a = ap.parse_args()
    res = []
    for N in a.N:
        res += run(N)
        json.dump(res, open(os.path.join(HERE, "results", "j1j2_infidelity_structure.json"), "w"), indent=1)
        open(os.path.join(HERE, "results", "j1j2_infidelity_structure.md"), "w").write(table(res) + "\n")
    print(table(res))


if __name__ == "__main__":
    main()
