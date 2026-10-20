"""s2_ops.py -- the total-spin operator S^2 as a filter operator (Stetcu-Baroni-Carlson symmetry projection) and its Trotter certificate.

S^2 = 3N/4 + sum_{i<j} sigma_i.sigma_j / 2 = const + sum_{i<j} P_ij  (P_ij = swap),  eigenvalues J(J+1).
Jucys-Murphy blocks X_k = sum_{i<k} P_ik (k = 2..N) commute with each other, so  exp(-i t sum_k X_k) = prod_k exp(-i t X_k)  EXACTLY;
a Lie-Trotter circuit that orders the Pauli terms block by block therefore has error  sum_k alpha_k t^2/(2 m)  with
alpha_k = sum_{g in block k} || [H_g, sum_{g' > g, g' in block k} H_g'] ||   (terms of different blocks never enter).
Inside a block the only anticommuting pairs are (i,k,a),(j,k,b) with i != j, a != b (they share qubit k with different Paulis).
"""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from qiskit.quantum_info import SparsePauliOp

PAULIS = "XYZ"


def s2_terms(N):
    """List of (block k, qubits (i,k), pauli, coeff) in block-by-block order (k = 1..N-1 zero-based)."""
    return [(k, (i, k), P, 0.5) for k in range(1, N) for i in range(k) for P in PAULIS]


def s2_pauli_op(N):
    """sum_{i<j} sigma_i.sigma_j / 2  (the 3N/4 constant is dropped: it only shifts all J(J+1) by a known number), block order."""
    return SparsePauliOp.from_sparse_list([(P + P, list(q), c) for _, q, P, c in s2_terms(N)], num_qubits=N)


def s2_sparse(N):
    """S^2 - 3N/4 + N(N-1)/4 ... implemented as sum_{i<j} P_ij ; full 2^N space.  Eigenvalues: J(J+1) - 3N/4 + N(N-1)/4 ... we return S^2 itself."""
    s = np.arange(1 << N, dtype=np.int64); R, C, V = [], [], []; diag = np.zeros(1 << N)
    for i in range(N):
        for j in range(i + 1, N):
            bi, bj = (s >> i) & 1, (s >> j) & 1
            same = bi == bj
            diag += np.where(same, 1.0, 0.0)                       # P_ij |same> = |same>
            src = np.where(~same)[0]
            R.append(src ^ ((1 << i) | (1 << j))); C.append(src); V.append(np.ones(len(src)))     # P_ij swaps unequal pairs
    # S^2 = sum_{i<j} P_ij + 3N/4 - N(N-1)/2 * 1/2 * ... : use S^2 = 3N/4 + sum_{i<j}(2 P_ij - 1)/2 = 3N/4 - N(N-1)/4 + sum_{i<j} P_ij
    const = 3 * N / 4 - N * (N - 1) / 4
    return (sp.diags(diag + const) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(1 << N, 1 << N))).tocsr()


def j_weights(vec, N, S2=None):
    """Weight of a state in each total-spin sector J = 0 .. N/2 (even N), via the exact projectors prod_{J' != J}(S^2 - J'(J'+1)) / (J(J+1) - J'(J'+1))."""
    S2 = s2_sparse(N) if S2 is None else S2
    v = np.asarray(vec, complex); v = v / np.linalg.norm(v); out = {}
    Js = list(range(0, N // 2 + 1))
    for J in Js:
        w = v.copy()
        for Jp in Js:
            if Jp != J:
                w = (S2 @ w - Jp * (Jp + 1) * w) / (J * (J + 1) - Jp * (Jp + 1))
        out[J] = float(np.real(np.vdot(v, w)))
    return out


def zero_pulses(N, tol=1e-9):
    """Smallest set of times t_i = pi/2^(i+1), i = 1..m (zero phases) with prod_i cos(t_i J(J+1)) = 0 for every J = 1..N/2 (even N)."""
    Js = np.arange(1, N // 2 + 1)
    for m in range(1, 12):
        t = np.pi / 2.0 ** (np.arange(1, m + 1) + 1)
        f = np.prod(np.cos(np.outer(Js * (Js + 1), t)), axis=1)
        if np.all(np.abs(f) < tol):
            return t
    raise RuntimeError("no zero-pulse set found")


def alpha_blocks(N, exact=True, nmax_exact=12):
    """Block-aware Trotter constant alpha_S2 = sum_k alpha_k (exact norms of the inner sums for N <= nmax_exact, else the triangle bound).
    Returns (alpha_triangle, alpha_exact_or_None)."""
    terms = s2_terms(N); tri = 0.0; ex = 0.0 if (exact and N <= nmax_exact) else None
    for k in range(1, N):
        blk = [(i, P) for i in range(k) for P in PAULIS]              # block order
        for gi, (i, a) in enumerate(blk):
            later = [(j, b) for (j, b) in blk[gi + 1:] if j != i and b != a]    # anticommuting later terms
            if not later:
                continue
            tri += 2 * 0.5 * 0.5 * len(later)                         # 2|c_g| * sum|c'|
            if ex is not None:
                # B = sum_later (1/2) sigma_j^b sigma_k^b on qubits {j<k} U {k}; build exactly on k+1 qubits
                B = SparsePauliOp.from_sparse_list([(b + b, [j, k], 0.5) for (j, b) in later], num_qubits=k + 1)
                M = B.to_matrix(sparse=True)
                nrm = abs(sla.eigsh(M, k=1, which="LM", return_eigenvectors=False)[0]) if M.shape[0] > 4 else np.max(np.abs(np.linalg.eigvalsh(M.toarray())))
                ex += 2 * 0.5 * nrm
    return tri, ex
