"""The paper's test system: the Heisenberg model H = sum_<ij> sigma_i . sigma_j + h sum_i sigma_i^z on an
Lx x Ly periodic square lattice (4 x 4 in the paper), in the J_z = 0 sector, with the Neel initial state and the
total-spin operator J^2 (eigenvalues J(J+1)).

sigma_i . sigma_j = 2 SWAP_ij - 1, and J^2 = sum_{i<j} SWAP_ij + 3N/4 - N(N-1)/4 for spin-1/2.
States are bit strings, bit (N - 1 - i) = 1 meaning site i is up.
"""
from itertools import combinations

import numpy as np
import scipy.sparse as sp


def square_bonds(Lx, Ly, periodic=True):
    """Nearest-neighbour bonds (i, j) on an Lx x Ly lattice, site = x + Lx y; each bond once."""
    bonds = set()
    for y in range(Ly):
        for x in range(Lx):
            i = x + Lx * y
            for dx, dy in ((1, 0), (0, 1)):
                xx, yy = x + dx, y + dy
                if not periodic and (xx >= Lx or yy >= Ly):
                    continue
                j = (xx % Lx) + Lx * (yy % Ly)
                if i != j:
                    bonds.add((min(i, j), max(i, j)))
    return sorted(bonds)


class SpinSector:
    """Basis of N spin-1/2 with n_up up spins, and SWAP_ij as permutations of it."""

    def __init__(self, N, n_up=None):
        self.N = N
        self.n_up = N // 2 if n_up is None else n_up
        self.states = np.array(sorted(sum(1 << (N - 1 - q) for q in c) for c in combinations(range(N), self.n_up)))
        self.dim = len(self.states)
        self._index = {int(s): k for k, s in enumerate(self.states)}

    def index(self, s):
        return self._index[int(s)]

    def swap_perm(self, i, j):
        N, s = self.N, self.states
        mi, mj = 1 << (N - 1 - i), 1 << (N - 1 - j)
        flip = ((s & mi) > 0) != ((s & mj) > 0)
        t = np.where(flip, s ^ mi ^ mj, s)
        return np.array([self._index[int(x)] for x in t])

    def swap_sum(self, pairs, weights=None):
        """sum_k w_k SWAP_{pairs[k]} as a sparse matrix."""
        weights = np.ones(len(pairs)) if weights is None else weights
        rows, cols, vals = [], [], []
        ar = np.arange(self.dim)
        for (i, j), w in zip(pairs, weights):
            rows.append(self.swap_perm(i, j)); cols.append(ar); vals.append(np.full(self.dim, w))
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.dim, self.dim))

    def sz_total(self):
        bits = (self.states[:, None] >> (self.N - 1 - np.arange(self.N))[None, :]) & 1
        return 0.5 * (2 * bits - 1).sum(1)


def heisenberg(sector, bonds, h=0.0):
    """H = sum_<ij> sigma_i . sigma_j + h sum_i sigma_i^z (Pauli matrices, as in the paper's Eq. 3)."""
    H = 2 * sector.swap_sum(bonds) - len(bonds) * sp.identity(sector.dim, format="csr")
    if h:
        H = H + sp.diags(2 * h * sector.sz_total())
    return H.tocsr()


def total_spin_sq(sector):
    """J^2 for spin-1/2: sum_{i<j} SWAP_ij + 3N/4 - N(N-1)/4."""
    N = sector.N
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    c = 0.75 * N - 0.25 * N * (N - 1)
    return (sector.swap_sum(pairs) + c * sp.identity(sector.dim, format="csr")).tocsr()


def neel(sector, Lx, Ly):
    """Spins up on the even sublattice ((x + y) even), down on the odd one; J_z = 0."""
    s = 0
    for y in range(Ly):
        for x in range(Lx):
            if (x + y) % 2 == 0:
                s |= 1 << (sector.N - 1 - (x + Lx * y))
    v = np.zeros(sector.dim)
    v[sector.index(s)] = 1.0
    return v
