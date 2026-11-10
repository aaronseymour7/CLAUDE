"""Total spin: the S^2 operator in the S^z = 0 sector, the singlet gap, and the SBC S^2 projection stage.

S^2 = c + sum_{i<j} SWAP_ij,  c = 3N/4 - N(N-1)/4,  eigenvalues J(J+1).

S^2 stage (Stetcu-Baroni-Carlson symmetry projection, all phases zero): pulses cos(t_k S^2), t_k = pi / 2^(k+1),
k = 1..L. J(J+1) = 2m, and cos(pi m / 2^k) = 0 exactly when m = 2^(k-1) x odd, so the ladder removes every J whose
v2(J(J+1)/2) < L; L = 1 + max_J v2(J(J+1)/2) over J <= N/2 removes all J > 0 (L = 3 up to N = 16).

Each pulse is H . controlled e^{-i t S^2 Z_a} . H on its own ancilla, like the energy pulses. e^{-i t S^2} is
Trotterised with second-order Strang steps over all N(N-1)/2 pairs; every factor e^{-i b SWAP_ij} is SU(2)
invariant, so the Trotterised pulse still conserves S^2 exactly. It is just not an exact function of S^2, which
leaves a little J > 0 weight and slightly rotates the singlet sector; the measured fidelity includes both.
Gate cost per controlled pair is the same 6-CNOT block as a Heisenberg bond (all-to-all connectivity, as for the
J2 bonds), with the same merging of consecutive Strang steps.
"""
from functools import lru_cache

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def pairs(N):
    return [(i, j) for i in range(N) for j in range(i + 1, N)]


def pair_perms(ch):
    """SWAP_ij of every pair as an index array on the sector basis."""
    N, s = ch.N, ch.states
    out = []
    for i, j in pairs(N):
        mi, mj = 1 << (N - 1 - i), 1 << (N - 1 - j)
        flip = ((s & mi) > 0) != ((s & mj) > 0)
        out.append(ch.index_of(np.where(flip, s ^ mi ^ mj, s)))
    return out


def s2_const(N):
    return 0.75 * N - 0.25 * N * (N - 1)


def s2_matrix(ch):
    """Sparse S^2 in the sector."""
    P = pair_perms(ch)
    d = ch.dim
    S = sp.csr_matrix((np.ones(d * len(P)), (np.concatenate(P), np.tile(np.arange(d), len(P)))), shape=(d, d))
    return S + s2_const(ch.N) * sp.identity(d, format="csr")


def lowest_singlets(H, S2, W, k=2):
    """k lowest S = 0 levels: the lowest of H + W S^2 (every J >= 1 level is pushed up by at least 2W)."""
    A = H + W * S2
    if A.shape[0] <= 2000:
        w = np.linalg.eigvalsh(A.toarray())[:k]
    else:
        w = np.sort(eigsh(A, k=k, which="SA", tol=1e-12, return_eigenvectors=False))
    return w


@lru_cache(maxsize=None)
def ladder(N):
    """S^2 pulse times pi / 2^(k+1), k = 1..L, removing every J = 1..N/2."""
    def v2(m):
        return (m & -m).bit_length() - 1
    L = 1 + max(v2(J * (J + 1) // 2) for J in range(1, N // 2 + 1))
    return tuple(np.pi / 2 ** (k + 1) for k in range(1, L + 1))


def steps(n, k):
    """Strang steps for pulse k (0-based): n for the longest, halved with the time (at least 1)."""
    return max(1, int(np.ceil(n / 2 ** k)))


def _strang(perms, V, dt, m, c):
    """m Strang steps of e^{-i dt S^2} on V (dim or dim x cols)."""
    b = dt / 2
    a, s = np.cos(b), -1j * np.sin(b)
    seq = list(perms) + list(perms)[::-1]
    for _ in range(m):
        for p in seq:
            V = a * V + s * V[p]
    return np.exp(-1j * c * dt * m) * V


def project_singlet(ch, V, S2=None):
    """Exact projector onto J = 0: prod_{J=1..N/2} (S^2 - J(J+1)) / (-J(J+1)). This is also the ideal S^2 stage
    (every pulse is cos(t * 0) = 1 on J = 0 and the ladder zeroes every J > 0)."""
    S2 = s2_matrix(ch) if S2 is None else S2
    V = np.array(V, dtype=complex)
    for J in range(1, ch.N // 2 + 1):
        x = J * (J + 1)
        V = (S2 @ V - x * V) / (-x)
    return V


def stage(ch, V, n, perms=None):
    """Apply the Trotterised S^2 stage (n Strang steps for the longest pulse) to V (dim or dim x cols).
    Returns the unnormalised post-selected output."""
    V = np.array(V, dtype=complex)
    c = s2_const(ch.N)
    perms = pair_perms(ch) if perms is None else perms
    for k, t in enumerate(ladder(ch.N)):
        m = steps(n, k)
        V = 0.5 * (_strang(perms, V, t / m, m, c) + _strang(perms, V, -t / m, m, c))
    return V


def stage_blocks(N, n):
    """Controlled pair blocks of the whole stage (merged Strang steps, like the energy pulses)."""
    P = len(pairs(N))
    return sum(steps(n, k) * (2 * P - 2) + 1 for k in range(len(ladder(N))))


def stage_cnots(N, n, bond_cnots=6):
    return bond_cnots * stage_blocks(N, n)


def singlet_weight(ch, V):
    """Weight of V (dim or dim x cols) in the J = 0 sector."""
    return np.sum(np.abs(project_singlet(ch, V)) ** 2, axis=0)


def _infid(a, b):
    a, b = a.reshape(len(a), -1), b.reshape(len(b), -1)
    na, nb = np.sum(abs(a) ** 2, 0), np.sum(abs(b) ** 2, 0)
    return float(np.max(1 - np.abs(np.sum(a.conj() * b, 0)) ** 2 / (na * nb)))
