"""The open J1-J2 chain in its S^z = 0 sector, and the two knobs.

H = sum_i J1 (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, even N.

Gap knob      dimerisation delta. It keeps SU(2): the ground state is a singlet in S^z = 0 and the first
              excitation a triplet, also reachable in S^z = 0. delta < 0 closes the gap (weak end bonds, nearly
              free edge spins), delta -> 1 opens it to about 1.8 J1 (isolated dimers).
Overlap knob  cfg.knob = "staggered_z" (default):
              |trial> = exp(-i eps sum_i (-1)^i Z_i) |base>: one layer of staggered Rz after the state preparation.
              |base> is the one-layer MPS-to-circuit approximation of |E0> (cfg.trial = "mps1", application; see ornl.prep),
              or |E0> itself (cfg.trial = "exact", the default for the study). gamma = |<E0|trial>|^2 includes the preparation loss, so at
              eps = 0 it is gamma_0 = |<E0|base>|^2 and targets above gamma_0 are out of reach.
              It conserves S^z but not S^2: 90-95% of the missing weight goes to S > 0, mostly the low triplets.
              cfg.knob = "singlet_bond": |trial> = prod_{odd i} e^{+i eps SWAP_{i,i+1}} prod_{even i} e^{-i eps SWAP_{i,i+1}} |base>,
              one brick layer of SU(2)-invariant two-qubit gates (3 CNOT each, counted in knob_cnot). The trial
              stays a singlet, so all missing weight sits on excited singlets; the overlap falls to ~0 at eps = pi/2.

Gap kinds     "triplet": E1 - E0 in the S^z = 0 sector (a triplet). "singlet": the first excited singlet, the gap
              that matters once the trial has no S > 0 weight (singlet knob, or after the S^2 stage). In the
              singlet case the dimerisation knob has a minimum at small negative delta (two dimer patterns
              compete), so delta_for_gap searches only to the right of it.

Basis convention: a state is an integer whose bit (N - 1 - i) is site i (site 0 = most significant bit).
For two spin-1/2s, S_i.S_j = (2 SWAP_ij - 1)/4, and SWAP_ij is a permutation of the S^z = 0 basis (`perm`).
"""
from itertools import combinations

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigsh


KNOBS = ("staggered_z", "singlet_bond")


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

    def __init__(self, N, delta, J1=1.0, J2=0.2411, knob="staggered_z", spectrum=True):
        assert N % 2 == 0, "even N only"
        if knob not in KNOBS:
            raise ValueError(f"knob must be one of {KNOBS}")
        self.N, self.delta, self.J1, self.J2, self.knob = N, delta, J1, J2, knob
        self.knob_cnot = 3 * (N - 1) if knob == "singlet_bond" else 0
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
        self.gap_kind = "triplet"
        if spectrum:
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
        self.gap_triplet = self.gap
        self.g0 = g0 / np.linalg.norm(g0)
        self.set_base(self.g0)

    def singlet_gap(self):
        """E1(S = 0) - E0 (the ground state is a singlet)."""
        if not hasattr(self, "_gap_singlet"):
            from .symmetry import lowest_singlets, s2_matrix
            lam = sum(abs(J) for _, _, J in self.bonds)                        # >= the spectral width
            w = lowest_singlets(self.hamiltonian(), s2_matrix(self), lam)
            self._gap_singlet = float(w[1] - w[0])
        return self._gap_singlet

    def use_gap(self, kind):
        """Make gap / gap_s the gap the filter has to resolve: "triplet" (E1 in the sector) or "singlet"."""
        self.gap_kind = kind
        self.gap = self.gap_triplet if kind == "triplet" else self.singlet_gap()
        self.gap_s = self.gap / self.W

    # ---- overlap knob
    def set_base(self, v, prep_cnot=0, label="exact"):
        """State the overlap knob acts on (default: the exact ground state), and the CNOTs that prepare it."""
        self.base, self.prep_cnot, self.trial_mode = v / np.linalg.norm(v), int(prep_cnot), label
        self.gamma0 = float(abs(np.vdot(self.g0, self.base)) ** 2)

    def trial(self, eps):
        if self.knob == "staggered_z":
            return np.exp(-1j * eps * self.M) * self.base
        v = np.array(self.base, dtype=complex)                  # singlet_bond: even bonds e^{-i eps SWAP}, odd e^{+i eps SWAP}
        for b, (i, j, _) in enumerate(self.bonds):
            if j == i + 1:
                th = eps if i % 2 == 0 else -eps
                v = np.cos(th) * v - 1j * np.sin(th) * v[self.perm[b]]
        return v

    def overlap(self, eps):
        return abs(np.vdot(self.g0, self.trial(eps))) ** 2

    def eps_for_overlap(self, gamma):
        """Smallest eps with overlap = gamma. Staggered knob: on the first descent only (it falls monotonically up to
        its first minimum); singlet knob: the first crossing on [0, pi/2]."""
        if self.knob == "singlet_bond":
            # the brick-layer overlap can have a shallow local minimum before it falls to ~0 at pi/2: take the
            # first grid crossing of gamma and refine it
            es = np.linspace(0, np.pi / 2, 1571)
            g = np.array([self.overlap(e) for e in es])
            if not gamma <= g[0]:
                return None
            hit = np.flatnonzero(g <= gamma)
            if not len(hit):
                return None
            i = hit[0]
            if g[i] == gamma or i == 0:
                return float(es[i])
            return brentq(lambda e: self.overlap(e) - gamma, es[i - 1], es[i], xtol=1e-14)
        es = np.linspace(0, 1.5, 1501)
        g = np.array([self.overlap(e) for e in es])
        rise = np.flatnonzero(np.diff(g) > 0)
        stop = rise[0] if len(rise) else len(es) - 1
        if not g[stop] <= gamma <= g[0]:                    # g[0] = gamma_0, the overlap of the prepared state
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


def singlet_gap(N, delta, J1=1.0, J2=0.2411):
    return Chain(N, delta, J1, J2, spectrum=False).singlet_gap()


def singlet_dip(N, J1=1.0, J2=0.2411, lo=-0.5, hi=0.2):
    """delta where the singlet gap is smallest (two dimer patterns cross), and that gap."""
    r = minimize_scalar(lambda d: singlet_gap(N, d, J1, J2), bounds=(lo, hi), method="bounded",
                        options=dict(xatol=1e-6))
    return float(r.x), float(r.fun)


def delta_for_gap(N, gap, J1=1.0, J2=0.2411, lo=-0.9, hi=0.95, kind="triplet"):
    """Dimerisation that gives the target physical gap (units of J1). kind = "singlet" searches the branch to the
    right of the singlet dip; returns None if the target is below the dip or above the value at hi."""
    if kind == "triplet":
        return brentq(lambda d: Chain(N, d, J1, J2).gap - gap, lo, hi, xtol=1e-10)
    d0, g0 = singlet_dip(N, J1, J2)
    f = lambda d: singlet_gap(N, d, J1, J2) - gap
    if gap < g0 or f(hi) < 0:
        return None
    return brentq(f, d0, hi, xtol=1e-10)


def filter_gap_kind(cfg):
    """The gap the energy filter must resolve: singlet when the state entering it has no S > 0 weight."""
    return "singlet" if (cfg.knob == "singlet_bond" or cfg.s2_stage) else "triplet"


def chain_for_gap(N, gap, cfg):
    """Chain at the target gap (of kind cfg.gap_kind) with the model settings of a StudyConfig, or None if a
    singlet target is out of reach. cfg.trial = "mps1" prepares the trial base with the one-layer MPS-to-circuit
    (ornl.prep); "exact" uses the exact ground state (no prep cost). ch.gap / ch.gap_s are then set to the gap the
    filter resolves (filter_gap_kind)."""
    d = delta_for_gap(N, gap, cfg.J1, cfg.J2, *cfg.delta_range, kind=cfg.gap_kind)
    if d is None:
        return None
    ch = Chain(N, d, cfg.J1, cfg.J2, knob=cfg.knob)
    ch.use_gap(filter_gap_kind(cfg))
    if cfg.trial == "mps1":
        from .prep import prepare
        v, gates, c, _ = prepare(ch)
        ch.set_base(v, c, "mps1")
        ch.prep_gates = gates
    elif cfg.trial != "exact":
        raise ValueError("cfg.trial must be 'mps1' or 'exact'")
    return ch
