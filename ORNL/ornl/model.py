"""The open J1-J2 chain in its S^z = 0 sector, and the two knobs.

H = sum_i J1 (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, even N.

Gap knob      dimerisation delta. It keeps SU(2): the ground state is a singlet in S^z = 0 and the first
              excitation a triplet, also reachable in S^z = 0. delta < 0 closes the gap (weak end bonds, nearly
              free edge spins), delta -> 1 opens it to about 1.8 J1 (isolated dimers).
Overlap knob  |trial> = exp(-i eps sum_i (-1)^i Z_i) |E0>: one layer of staggered Rz on the exact ground state.
              It conserves S^z and moves the missing weight onto the low triplets, E1 included (the hard case).

Basis convention: a state is an integer whose bit (N - 1 - i) is site i (site 0 = most significant bit).
For two spin-1/2s, S_i.S_j = (2 SWAP_ij - 1)/4, and SWAP_ij is a permutation of the S^z = 0 basis (`perm`).
"""
from itertools import combinations

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigsh


def chain_bonds(N, delta, J1=1.0, J2=0.2411):
    """Bonds (i, j, J) in Trotter order: J1 even, J1 odd, J2 even, J2 odd."""
    bonds = []
    for r, J in ((1, J1), (2, J2)):
        for par in (0, 1):
            for i in range(par, N - r, 2):
                f = (1 + delta if i % 2 == 0 else 1 - delta) if r == 1 else 1.0
                bonds.append((i, i + r, J * f))
    return bonds


class Chain:
    """S^z = 0 sector of one chain: basis, bonds in Trotter order, spectrum, exact ground state, trial states.

    Attributes: N, delta, dim, states, bonds [(i, j, J)], perm (SWAP_ij as an index array per bond), M (staggered
    magnetisation, diagonal), E0, E1, Emax, gap, W = Emax - E0, gap_s = gap / W, g0 (ground state), c_id."""

    def __init__(self, N, delta, J1=1.0, J2=0.2411):
        assert N % 2 == 0, "even N only"
        self.N, self.delta, self.J1, self.J2 = N, delta, J1, J2
        self.states = np.array(sorted(sum(1 << (N - 1 - q) for q in c) for c in combinations(range(N), N // 2)))
        self.dim = len(self.states)
        bits = (self.states[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1
        self.M = ((1 - 2 * bits) * (-1) ** np.arange(N)).sum(1).astype(float)       # sum_i (-1)^i Z_i
        self.bonds = chain_bonds(N, delta, J1, J2)
        index = np.full(1 << N, -1)
        index[self.states] = np.arange(self.dim)
        self.index_of = lambda s: index[s]
        self.perm = []
        for i, j, _ in self.bonds:
            mi, mj = 1 << (N - 1 - i), 1 << (N - 1 - j)
            s = self.states
            flip = ((s & mi) > 0) != ((s & mj) > 0)
            self.perm.append(index[np.where(flip, s ^ mi ^ mj, s)])
        self.c_id = -0.25 * sum(J for _, _, J in self.bonds)                          # identity part of H
        self._spectrum()

    def hamiltonian(self):
        """Sparse H in the sector: sum_b J_b/4 (2 SWAP_b - 1)."""
        rows = np.concatenate(self.perm)
        cols = np.tile(np.arange(self.dim), len(self.bonds))
        vals = np.concatenate([np.full(self.dim, 0.5 * J) for _, _, J in self.bonds])
        H = csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim))
        return H + self.c_id * identity(self.dim, format="csr")

    def _spectrum(self):
        H = self.hamiltonian()
        if self.dim <= 2000:
            w, V = np.linalg.eigh(H.toarray())
            E0, E1, Emax, g0 = w[0], w[1], w[-1], V[:, 0]
        else:
            w, V = eigsh(H, k=2, which="SA", tol=1e-13)
            E0, E1, g0 = w[0], w[1], V[:, 0]
            Emax = eigsh(H, k=1, which="LA", tol=1e-10, return_eigenvectors=False)[0]
        self.E0, self.E1, self.Emax = E0, E1, Emax
        self.gap, self.W = E1 - E0, Emax - E0
        self.gap_s = self.gap / self.W
        self.g0 = g0 / np.linalg.norm(g0)

    # ---- overlap knob
    def trial(self, eps):
        return np.exp(-1j * eps * self.M) * self.g0

    def overlap(self, eps):
        return abs(np.vdot(self.g0, self.trial(eps))) ** 2

    def eps_for_overlap(self, gamma):
        """Smallest eps with overlap = gamma (the overlap falls monotonically up to its first minimum)."""
        es = np.linspace(0, 1.5, 1501)
        g = np.array([self.overlap(e) for e in es])
        rise = np.flatnonzero(np.diff(g) > 0)
        stop = rise[0] if len(rise) else len(es) - 1
        if not g[stop] <= gamma <= 1.0:
            return None
        return brentq(lambda e: self.overlap(e) - gamma, 0.0, es[stop], xtol=1e-14)

    # ---- full-register vectors (for Qiskit): sector vector -> 2^N vector in Qiskit's little-endian order
    def to_full(self, v):
        """Qiskit index = sum_i bit_i 2^i with bit_i the spin of site i, i.e. site i = qubit i."""
        N = self.N
        q = np.zeros(len(self.states), dtype=np.int64)
        for i in range(N):
            q |= ((self.states >> (N - 1 - i)) & 1) << i
        out = np.zeros(1 << N, dtype=complex)
        out[q] = v
        return out


def delta_for_gap(N, gap, J1=1.0, J2=0.2411, lo=-0.9, hi=0.95):
    """Dimerisation that gives the target physical gap (units of J1)."""
    return brentq(lambda d: Chain(N, d, J1, J2).gap - gap, lo, hi, xtol=1e-10)


def chain_for_gap(N, gap, cfg):
    """Chain at the target gap with the model settings of a StudyConfig."""
    d = delta_for_gap(N, gap, cfg.J1, cfg.J2, *cfg.delta_range)
    return Chain(N, d, cfg.J1, cfg.J2)
