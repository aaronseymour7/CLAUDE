"""Where the infidelity of the one-layer mps-to-circuit trial sits, by total spin S and by energy level.

    python j1j2/trial_spectrum.py            # N = 8..16 -> results/j1j2_trial_spectrum.{json,md}

H and S^2 commute (every term is an isotropic S_i.S_j), and both conserve S^z. In each S^z sector the common
eigenbasis is found by diagonalising H + a S^2 (a = small irrational number, which splits the degeneracies) densely;
S comes from <S^2> = S(S+1) and E from <H>. The trial's weight on each (S, E) level is then exact. Levels are
counted within each S, so (S = 0, level 1) is the first excited singlet.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import inputs  # noqa: E402

A = 1e-3 * np.pi
TOL = 1e-6


def decompose(N, seed=0):
    est = inputs.dmrg_estimates(N, seed=seed)
    _, psi = inputs.trial(est["psi0"])
    H, S2 = inputs.hamiltonian(N), inputs.total_spin_sq(N)
    pop = np.array([bin(i).count("1") for i in range(1 << N)])
    levels = []                                                  # (Sz, S, E, weight)
    for nup in range(N + 1):
        idx = np.flatnonzero(pop == nup)
        w_sector = float(np.sum(abs(psi[idx]) ** 2))
        if w_sector < 1e-14:
            continue
        h, s2 = H[idx][:, idx].toarray(), S2[idx][:, idx].toarray()
        _, V = np.linalg.eigh(h + A * s2)
        E = np.einsum("ij,ij->j", V.conj(), h @ V).real
        J = np.einsum("ij,ij->j", V.conj(), s2 @ V).real
        S = np.abs(np.round((-1 + np.sqrt(1 + 4 * np.maximum(J, 0))) / 2 * 2) / 2) + 0.0
        w = abs(V.conj().T @ psi[idx]) ** 2
        for e, s, x in zip(E, S, w):
            levels.append((nup - N / 2, float(s), float(e), float(x)))
    lv = np.array(levels)
    E0 = lv[:, 2].min()
    out = dict(N=N, gamma=None, by_S={}, singlet_levels=[], by_Sz={})
    for s in sorted(set(lv[:, 1])):
        m = lv[:, 1] == s
        out["by_S"][s] = float(lv[m, 3].sum())
    for sz in sorted(set(lv[:, 0])):
        out["by_Sz"][sz] = float(lv[lv[:, 0] == sz, 3].sum())
    # distinct energy levels per S (merge the 2S+1 S^z members)
    for s in sorted(set(lv[:, 1]))[:3]:
        m = lv[:, 1] == s
        E, W = lv[m, 2], lv[m, 3]
        o = np.argsort(E)
        E, W = E[o], W[o]
        groups = []
        for e, w in zip(E, W):
            if groups and abs(e - groups[-1][0]) < TOL:
                groups[-1][1] += w
            else:
                groups.append([e, w])
        out[f"levels_S{s:g}"] = [dict(level=k, E=g[0], dE=g[0] - E0, weight=g[1]) for k, g in enumerate(groups[:6])]
        out[f"rest_S{s:g}"] = float(sum(g[1] for g in groups[6:]))
    out["gamma"] = out["levels_S0"][0]["weight"]
    out["E0"] = E0
    return out


def main():
    Ns = [int(x) for x in sys.argv[1:]] or [8, 10, 12, 14, 16]
    res = []
    for N in Ns:
        r = decompose(N)
        res.append(r)
        print(N, round(1 - r["gamma"], 4), {k: round(v, 4) for k, v in r["by_S"].items() if v > 1e-5},
              [(x["level"], round(x["dE"], 3), round(x["weight"], 4)) for x in r["levels_S0"][:4]],
              [(x["level"], round(x["dE"], 3), round(x["weight"], 4)) for x in r["levels_S1"][:3]], flush=True)
        json.dump(res, open(os.path.join(HERE, "results", "j1j2_trial_spectrum.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
