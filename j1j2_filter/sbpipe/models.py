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


def j1j2_square(Lx: int, Ly: int, J1: float = 1.0, J2: float = 0.0, pbc_y: bool = True) -> PauliSum:
    """J1-J2 Heisenberg model on an Lx x Ly square lattice, open in x and (by default) periodic in y
    (a cylinder). Site (x, y) is qubit x*Ly + y, so the MPS order runs up each column in turn.

    Term order (= Trotter order): horizontal, then vertical nearest-neighbour bonds, then the two
    diagonal next-nearest-neighbour families, each split by the parity of its first coordinate;
    XX, YY, ZZ within each bond. One group id per bond."""
    if pbc_y and Ly < 3:
        raise ValueError("a periodic direction needs Ly >= 3")
    N = Lx * Ly
    q = lambda x, y: x * Ly + (y % Ly)
    fam = []
    fam.append([(x, q(x, y), q(x + 1, y)) for x in range(Lx - 1) for y in range(Ly)])
    fam.append([(y, q(x, y), q(x, y + 1)) for x in range(Lx) for y in range(Ly if pbc_y else Ly - 1)])
    if J2 != 0.0:
        for dy in (1, -1):
            fam.append([(x, q(x, y), q(x + 1, y + dy)) for x in range(Lx - 1) for y in range(Ly)
                        if pbc_y or 0 <= y + dy < Ly])
    bonds = []
    for k, f in enumerate(fam):
        J = J1 if k < 2 else J2
        for par in (0, 1):
            bonds += [(min(i, j), max(i, j), J) for (c, i, j) in f if c % 2 == par]
    labels, coeffs, groups = [], [], []
    for g, (i, j, J) in enumerate(bonds):
        for p in "XYZ":
            labels.append(_two_site(N, i, j, p))
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
