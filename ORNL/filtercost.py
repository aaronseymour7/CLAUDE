"""Measured cost of a Trotterised cosine ground-state filter on the 1D J1-J2 chain.

Everything here is measured by exact state-vector simulation; nothing is a bound.

Model     H = sum_i J1 (1 + (-1)^i delta) S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, open chain, even N.
Gap knob  the dimerisation delta. It keeps SU(2): the ground state is a singlet in S^z = 0 and the first
          excitation a triplet, so E1 is reachable from S^z = 0.
Overlap   one layer of staggered Rz on the exact ground state, |trial> = exp(-i eps sum_i (-1)^i Z_i) |E0>.
knob      It conserves S^z and moves the missing weight onto the low triplets (E1 included): the hard case.
Filter    one-ancilla cosine filter. Pulse t: H_anc, controlled e^{-i t H_s}, H_anc; post-selecting the ancilla
          on 0 applies cos(t H_s) with H_s = (H - E0) / (Emax - E0), spectrum in [0, 1].
          Dyadic schedule: t = tau 2^k for k < K, the ladder repeated r times, tau = pi / (1 + Delta_s).
Trotter   second-order Strang, n steps per tau, so pulse k has n 2^k steps. The Trotterised pulse is
          (S^m + S^-m) / 2 with S one Strang step, which is what the code applies.
Measured  for every trial, the cheapest (K, r, n) whose simulated 1 - F <= eps, where F is the fidelity of the
cost      post-selected state with the exact ground state. A choice only counts if the same (K, r) still meets
          eps at the next finer n, so a lucky Trotter cancellation is not counted.
Protocols repeat until success (RUS): cost = CNOT per attempt / p.
          amplitude amplification (AA): (2k+1) uses of the attempt with success sin^2((2k+1) theta),
          sin theta = sqrt(p); cost = CNOT per attempt * (2k+1) / sin^2((2k+1) theta).
          RUS and AA are optimised separately: each picks its own cheapest (K, r, n).
CNOTs     each Heisenberg bond controlled by the ancilla is one 6-CNOT block. Strang merges the outer half
          steps, so a pulse of k steps costs 6 (2 k (B - 1) + 1) CNOTs with B = 2N - 3 bonds.

The simulation runs in the S^z = 0 sector (dimension C(N, N/2), 12870 at N = 16), where a Heisenberg bond is
J/4 (2 SWAP - 1) and every Trotter factor is a + b SWAP, i.e. one gather. Only numpy and scipy are used.
"""
from itertools import combinations

import numpy as np
from scipy.optimize import brentq
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigsh

J1, J2 = 1.0, 0.2411            # J2 near the critical point of the uniform chain


# ---------------------------------------------------------------- model --------------------------------------------
class Chain:
    """S^z = 0 sector of an open J1-J2 chain with dimerisation delta: basis, bonds, spectrum, ground state."""

    def __init__(self, N, delta, J1=J1, J2=J2):
        assert N % 2 == 0
        self.N, self.delta = N, delta
        self.states = np.array(sorted(sum(1 << (N - 1 - q) for q in c) for c in combinations(range(N), N // 2)))
        self.dim = len(self.states)
        bits = (self.states[:, None] >> (N - 1 - np.arange(N))[None, :]) & 1
        self.M = ((1 - 2 * bits) * (-1) ** np.arange(N)).sum(1).astype(float)       # sum_i (-1)^i Z_i
        # bonds in Trotter order: J1 even, J1 odd, J2 even, J2 odd
        self.bonds = []
        for r, J in ((1, J1), (2, J2)):
            for par in (0, 1):
                for i in range(par, N - r, 2):
                    f = (1 + delta if i % 2 == 0 else 1 - delta) if r == 1 else 1.0
                    self.bonds.append((i, i + r, J * f))
        index = np.full(1 << N, -1)
        index[self.states] = np.arange(self.dim)
        self.index_of = lambda s: index[s]
        self.perm = []                                                                 # SWAP_ij as a gather
        for i, j, _ in self.bonds:
            mi, mj = 1 << (N - 1 - i), 1 << (N - 1 - j)
            s = self.states
            flip = ((s & mi) > 0) != ((s & mj) > 0)
            self.perm.append(index[np.where(flip, s ^ mi ^ mj, s)])
        self._spectrum()

    def hamiltonian(self):
        rows = np.concatenate(self.perm)
        cols = np.tile(np.arange(self.dim), len(self.bonds))
        vals = np.concatenate([np.full(self.dim, 0.5 * J) for _, _, J in self.bonds])
        H = csr_matrix((vals, (rows, cols)), shape=(self.dim, self.dim))
        self.c_id = -0.25 * sum(J for _, _, J in self.bonds)                            # identity part of H
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
        self.gap_s = self.gap / self.W                                                 # scaled gap Delta_s
        self.g0 = g0 / np.linalg.norm(g0)

    # overlap knob
    def trial(self, eps):
        return np.exp(-1j * eps * self.M) * self.g0

    def overlap(self, eps):
        return abs(np.vdot(self.g0, self.trial(eps))) ** 2

    def eps_for_overlap(self, gamma):
        """Smallest eps with overlap = gamma (overlap falls monotonically up to its first minimum)."""
        es = np.linspace(0, 1.5, 1501)
        g = np.array([self.overlap(e) for e in es])
        rise = np.flatnonzero(np.diff(g) > 0)
        stop = rise[0] if len(rise) else len(es) - 1
        if not g[stop] <= gamma <= 1.0:
            return None
        return brentq(lambda e: self.overlap(e) - gamma, 0.0, es[stop], xtol=1e-14)

    def min_overlap(self):
        es = np.linspace(0, 1.5, 1501)
        g = np.array([self.overlap(e) for e in es])
        rise = np.flatnonzero(np.diff(g) > 0)
        return g[rise[0] if len(rise) else -1]

    # Trotterised evolution
    def strang(self, V, dt, steps):
        """steps Strang steps of e^{-i dt H_s} on V (dim x m). Bond factor exp(-i a (J/4W)(2 SWAP - 1)): the
        constant part is collected into one phase, the rest is cos(c) - i sin(c) SWAP with c = a J / (2 W)."""
        half = [(np.cos(0.5 * J * dt / (2 * self.W)), -1j * np.sin(0.5 * J * dt / (2 * self.W)), p)
                for (_, _, J), p in zip(self.bonds, self.perm)]
        seq = half + half[::-1]
        V = np.array(V, dtype=complex)
        for _ in range(steps):
            for a, b, p in seq:
                V = a * V + b * V[p]
        c0 = (self.c_id - self.E0) / self.W                                            # identity part of H_s
        return np.exp(-1j * c0 * dt * steps) * V


def delta_for_gap(N, gap, lo=-0.9, hi=0.95, J1=J1, J2=J2):
    """Dimerisation that gives the target physical gap (units of J1)."""
    return brentq(lambda d: Chain(N, d, J1, J2).gap - gap, lo, hi, xtol=1e-10)


# ---------------------------------------------------------------- costs --------------------------------------------
def cnot_per_attempt(N, K, r, n):
    k = np.tile(n * 2 ** np.arange(K), r)
    B = 2 * N - 3
    return int(np.sum(6 * (2 * k * (B - 1) + 1)))


def aa_cost(c_attempt, p):
    """Expected cost per success with amplitude amplification; k = 0 is plain repetition."""
    p = np.asarray(p, float)
    th = np.arcsin(np.sqrt(np.clip(p, 1e-300, 1)))
    k = np.maximum(0, np.round(np.pi / (4 * th) - 0.5))
    return np.minimum(c_attempt * (2 * k + 1) / np.sin((2 * k + 1) * th) ** 2, c_attempt / p)   # AA only if it helps


def cost(c_attempt, p, protocol):
    return c_attempt / p if protocol == "rus" else float(aa_cost(c_attempt, p))


# ---------------------------------------------------------------- Trotter spectrum ---------------------------------
def step_spectrum(ch, dt):
    """Exact eigen-decomposition of one Strang step S(dt), so that every Trotterised filter can be evaluated at once:
    the pulse with m steps is (S^m + S^-m)/2 = O cos(m phi) O^T.

    How: S = A_1 ... A_B A_B ... A_1 with every bond factor A_b complex symmetric, so S = P P^T is complex symmetric
    and unitary, hence S = O e^{-i phi} O^T with O real orthogonal. O is found with a real symmetric eigensolver
    applied to Re S + alpha Im S. Spin flip commutes with every SWAP, so S is block diagonal in the spin-flip
    sectors, which halves the matrix size (at N = 16: 6435 instead of 12870).
    Returns (phi, project): eigenphases per step, and a function that writes vectors (dim x m) in the eigenbasis."""
    d, N = ch.dim, ch.N
    flip = ch.index_of(((1 << N) - 1) ^ ch.states)
    a = np.flatnonzero(np.arange(d) < flip)                       # one state of each spin-flip pair
    b = flip[a]
    half = [(np.cos(0.5 * J * dt / (2 * ch.W)), -1j * np.sin(0.5 * J * dt / (2 * ch.W)), p)
            for (_, _, J), p in zip(ch.bonds, ch.perm)]
    c0 = (ch.c_id - ch.E0) / ch.W
    phis, Os = [], []
    for sign in (1.0, -1.0):
        G = np.zeros((d, len(a)), dtype=complex)                   # columns: spin-flip (anti)symmetric states
        G[a, np.arange(len(a))] = 1 / np.sqrt(2)
        G[b, np.arange(len(a))] = sign / np.sqrt(2)
        for x, y, p in half:                                       # G = P^T Q
            Gp = G[p]
            G *= x
            Gp *= y
            G += Gp
            del Gp
        S = (G.T @ G) * np.exp(-1j * c0 * dt)                      # S restricted to the sector
        del G
        _, O = np.linalg.eigh(S.real + 0.1234 * S.imag)
        lam = np.einsum("ij,ij->j", O, S.real @ O) + 1j * np.einsum("ij,ij->j", O, S.imag @ O)
        phis.append(-np.angle(lam))
        Os.append((sign, O))

    def project(v):
        return np.concatenate([O.T @ ((v[a] + sign * v[b]) / np.sqrt(2)) for sign, O in Os], axis=0)
    return np.concatenate(phis), project


# ---------------------------------------------------------------- measured search ---------------------------------
N_GRID = (1, 2, 3, 4, 6, 8, 12, 16)          # Strang steps per tau
PROTOCOLS = ("rus", "aa")
K_MAX, R_MAX = 20, 64


def _scan(ch, trials, n):
    """Every dyadic filter (K <= K_MAX, r <= R_MAX) at n Strang steps per tau, for every trial.
    Returns {K: (p[r, j], infid[r, j])}, exact for the Trotterised circuit."""
    tau = np.pi / (1 + ch.gap_s)
    phi, project = step_spectrum(ch, tau / n)
    theta = n * phi                                                # phase per tau; pulse k: cos(2^k theta)
    c = project(trials)                                            # trial in the step eigenbasis
    u = project(ch.g0)
    w = np.abs(c) ** 2
    out, B = {}, np.ones_like(theta)
    for K in range(1, K_MAX + 1):
        B = B * np.cos(2.0 ** (K - 1) * theta)
        ps, infs, Br = [], [], np.ones_like(theta)
        for r in range(1, R_MAX + 1):
            Br = Br * B
            p = (Br ** 2) @ w
            amp = np.abs((u * Br) @ c) ** 2
            ps.append(p); infs.append(np.maximum(0.0, 1 - amp / np.maximum(p, 1e-300)))
        out[K] = (np.array(ps), np.array(infs))
    return out


def measure_point(ch, gammas, eps_list, verbose=False):
    """Measured cheapest filter for every target overlap and eps, under RUS and AA, at one Hamiltonian.
    Returns (rows, eps_knob): rows[j][eps][protocol] = dict(K, r, n, cnot, p, infid, cost, T_scaled) or None.
    The Strang step count n is raised along N_GRID until a level brings no cheaper choice."""
    eps_knob = [ch.eps_for_overlap(g) for g in gammas]
    ok = [j for j, e in enumerate(eps_knob) if e is not None]
    trials = np.stack([ch.trial(eps_knob[j]) for j in ok], axis=1)
    nkeys = len(ok) * len(eps_list) * len(PROTOCOLS)
    levels, best = [], {}
    for i, n in enumerate(N_GRID):
        levels.append((n, _scan(ch, trials, n)))
        new_best = _select(ch, levels, ok, eps_list)
        improved = any(k not in best or new_best[k]["cost"] < best[k]["cost"] for k in new_best)
        level_best = _select(ch, levels, ok, eps_list, unverified_last=True)
        could_improve = any(k not in new_best or level_best[k]["cost"] < new_best[k]["cost"] for k in level_best)
        best = new_best
        if verbose:
            print(f"    n={n}: {len(best)}/{nkeys} verified", flush=True)
        if i > 0 and len(best) == nkeys and not improved and not could_improve:
            break
    rows = [{e: {pr: None for pr in PROTOCOLS} for e in eps_list} for _ in gammas]
    for (jj, e, pr), b in best.items():
        rows[ok[jj]][e][pr] = b
    return rows, eps_knob


def _select(ch, levels, ok, eps_list, unverified_last=False):
    """Cheapest choice per (trial, eps, protocol). A choice at level i only counts if the same (K, r) also meets
    eps at level i + 1 (finer Trotter step)."""
    best = {}
    last = len(levels) - 1
    tau = np.pi / (1 + ch.gap_s)
    for i, (n, data) in enumerate(levels):
        check = levels[i + 1][1] if i < last else None
        if check is None and not unverified_last:
            continue
        for K, (p, inf) in data.items():
            for jj in range(len(ok)):
                for e in eps_list:
                    ok_r = inf[:, jj] <= e
                    if check is not None:
                        ok_r &= check[K][1][:, jj] <= e
                    hit = np.flatnonzero(ok_r)
                    if not len(hit):
                        continue
                    r = hit[0] + 1                                 # fewest repetitions that work (and verify)
                    c = cnot_per_attempt(ch.N, K, r, n)
                    for pr in PROTOCOLS:
                        cc = cost(c, p[r - 1, jj], pr)
                        key = (jj, e, pr)
                        if key not in best or cc < best[key]["cost"]:
                            best[key] = dict(K=K, r=int(r), n=n, cnot=c, p=float(p[r - 1, jj]),
                                             infid=float(inf[r - 1, jj]), cost=float(cc),
                                             trotter_steps=int(n * r * (2 ** K - 1)),
                                             T_scaled=float(tau * r * (2 ** K - 1)))
    return best
