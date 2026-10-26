"""Model Hamiltonians as PauliSums, and rigorous spectral-range bounds for any PauliSum."""
from __future__ import annotations

from collections import defaultdict
from typing import Tuple

import numpy as np

from .pauli import PauliSum, local_matrix, support


def _two_site(N, i, j, p):
    lab = ["I"] * N
    lab[i] = lab[j] = p
    return "".join(lab)


def j1j2_chain(N: int, J1: float = 1.0, J2: float = 0.0, pbc: bool = False) -> PauliSum:
    """H = J1 sum S_i.S_{i+1} + J2 sum S_i.S_{i+2},  S = sigma/2.

    Term order (= Trotter order): nearest-neighbour bonds with even, then odd left site, then
    next-nearest-neighbour bonds likewise; XX, YY, ZZ within each bond. One group id per bond."""
    bonds = []
    for r, J in ((1, J1), (2, J2)):
        if J == 0.0:
            continue
        pairs = [(i, i + r) for i in range(N - r)]
        if pbc and N > 2 * r:
            pairs += [(i, (i + r) % N) for i in range(N - r, N)]
        for par in (0, 1):
            bonds += [(i, j, J) for (i, j) in pairs if i % 2 == par]
    labels, coeffs, groups = [], [], []
    for g, (i, j, J) in enumerate(bonds):
        for p in "XYZ":
            labels.append(_two_site(N, i, j, p))
            coeffs.append(0.25 * J)
            groups.append(g)
    return PauliSum(labels, np.array(coeffs), groups)


def heisenberg_ladder(L: int, J_leg: float = 1.0, J_rung: float = 1.0) -> PauliSum:
    """Open two-leg spin-1/2 Heisenberg ladder with L rungs, N = 2L qubits, S = sigma/2.

    Site (rung l, leg s) is qubit 2l + s, so leg bonds connect qubits at distance 2 and rung bonds
    neighbours. Term order (= Trotter order): rung bonds, leg bonds of even rungs, leg bonds of odd
    rungs; XX, YY, ZZ within each bond; one group id per bond."""
    bonds = [(2 * l, 2 * l + 1, J_rung) for l in range(L)]
    for par in (0, 1):
        for l in range(L - 1):
            if l % 2 == par:
                bonds += [(2 * l + s, 2 * l + 2 + s, J_leg) for s in (0, 1)]
    labels, coeffs, groups = [], [], []
    for g, (i, j, J) in enumerate(bonds):
        for p in "XYZ":
            labels.append(_two_site(2 * L, i, j, p))
            coeffs.append(0.25 * J)
            groups.append(g)
    return PauliSum(labels, np.array(coeffs), groups)


def weyl_bounds(H: PauliSum) -> Tuple[float, float]:
    """Rigorous [lambda_min_lb, lambda_max_ub] of H from Weyl's inequality over its groups:
    lambda_max(sum_g h_g) <= sum_g lambda_max(h_g) (and likewise for lambda_min).
    Each group is diagonalised exactly on its support. Cost linear in the number of groups."""
    by_group = defaultdict(dict)
    c_id = 0.0
    for lab, c, g in zip(H.labels, H.coeffs, H.groups):
        if set(lab) == {"I"}:
            c_id += c
        else:
            by_group[g][lab] = by_group[g].get(lab, 0.0) + c
    lo = hi = c_id
    for terms in by_group.values():
        qs = sorted({q for l in terms for q in support(l)})
        w = np.linalg.eigvalsh(local_matrix(terms, qs))
        lo += float(w[0])
        hi += float(w[-1])
    return lo, hi
