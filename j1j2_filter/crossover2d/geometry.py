"""geometry.py -- open W x L square-lattice J1-J2 Heisenberg clusters (H = sum J S.S), qubit q = x*W + y (columns of height W, snake order: an MPS cut across
the middle crosses W sites).  Bonds: nearest neighbours (J1) and both diagonals (J2), open in both directions."""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from qiskit.quantum_info import SparsePauliOp


def bonds(W, L):
    idx = lambda x, y: x * W + y
    b1, b2 = [], []
    for x in range(L):
        for y in range(W):
            if y + 1 < W: b1.append((idx(x, y), idx(x, y + 1)))
            if x + 1 < L: b1.append((idx(x, y), idx(x + 1, y)))
            if x + 1 < L and y + 1 < W: b2.append((idx(x, y), idx(x + 1, y + 1)))
            if x + 1 < L and y - 1 >= 0: b2.append((idx(x, y), idx(x + 1, y - 1)))
    return sorted(set(map(lambda t: tuple(sorted(t)), b1))), sorted(set(map(lambda t: tuple(sorted(t)), b2)))


def pauli_hamiltonian(W, L, j1=1.0, j2=0.0):
    b1, b2 = bonds(W, L); terms = []
    for J, bl in ((j1, b1), (j2, b2)):
        if J == 0:
            continue
        for i, j in bl:
            for P in "XYZ":
                terms.append((P + P, [i, j], 0.25 * J))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=W * L)


def sector_matrix(W, L, nup, j1=1.0, j2=0.0):
    N = W * L; b1, b2 = bonds(W, L)
    s = np.arange(1 << N, dtype=np.int64); pc = np.zeros(1 << N, dtype=np.int8)
    for q in range(N):
        pc += ((s >> q) & 1).astype(np.int8)
    st = s[pc == nup]; del s, pc
    d = len(st); diag = np.zeros(d); R, C, V = [], [], []
    for J, bl in ((j1, b1), (j2, b2)):
        if J == 0:
            continue
        for i, j in bl:
            bi, bj = (st >> i) & 1, (st >> j) & 1
            diag += J * np.where(bi == bj, 0.25, -0.25)
            src = np.where(bi != bj)[0]; tgt = np.searchsorted(st, st[src] ^ ((1 << i) | (1 << j)))
            R.append(tgt.astype(np.int32)); C.append(src.astype(np.int32)); V.append(np.full(len(src), 0.5 * J))
    return (sp.diags(diag) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(d, d))).tocsr(), st
