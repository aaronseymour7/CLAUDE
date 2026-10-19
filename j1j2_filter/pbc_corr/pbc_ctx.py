"""pbc_ctx.py -- periodic J1-J2 chain inputs for the cosine-filter pipeline (exact classical inputs by ED, N <= ~20).

Everything downstream (floor design, certified step count n, exact post-selected filter, gate counts) is the pipeline's own code; only the
Hamiltonian (periodic bonds), the classical inputs (E0, E1, ground state by sector ED instead of DMRG) and the trial circuit (L = 1 brickwork from
an SVD-compressed MPS of the exact ground state) are new.  Convention as elsewhere: H = sum J S.S, qubit q = chain site q (qiskit little-endian).
"""
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from qiskit.quantum_info import SparsePauliOp, Statevector      # noqa: E402
from hamiltonians import scale_hamiltonian, reorder_axes        # noqa: E402


def pbc_hamiltonian(N, j1=1.0, j2=0.0):
    """Periodic J1-J2 chain as a Pauli sum (S.S = (XX+YY+ZZ)/4); r=1 bonds first, then r=2 bonds, X,Y,Z consecutive per bond (as the open chain)."""
    terms = []
    for r, J in ((1, j1), (2, j2)):
        if J == 0:
            continue
        for i in range(N):
            for P in "XYZ":
                terms.append((P + P, [i, (i + r) % N], 0.25 * J))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=N)


def _states(N, nup):
    return np.sort(np.array([sum(1 << b for b in c) for c in combinations(range(N), nup)], dtype=np.int64))


def sector_matrix(N, nup, j1, j2):
    """J1-J2 periodic Hamiltonian in the fixed-S^z sector (csr)."""
    st = _states(N, nup); d = len(st); diag = np.zeros(d); R, Cc, V = [], [], []
    for r, J in ((1, j1), (2, j2)):
        if J == 0:
            continue
        for i in range(N):
            j = (i + r) % N; bi, bj = (st >> i) & 1, (st >> j) & 1
            diag += J * np.where(bi == bj, 0.25, -0.25)
            src = np.where(bi != bj)[0]; tgt = np.searchsorted(st, st[src] ^ ((1 << i) | (1 << j)))
            R.append(tgt); Cc.append(src); V.append(np.full(len(src), 0.5 * J))
    return (sp.diags(diag) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(Cc))), shape=(d, d))).tocsr(), st


def full_sparse(N, j1, j2, shift=0.0, scale=1.0):
    """(H - shift)/scale on the full 2^N space (csr), built directly (no qiskit)."""
    dim = 1 << N; s = np.arange(dim, dtype=np.int64); diag = np.zeros(dim); R, Cc, V = [], [], []
    for r, J in ((1, j1), (2, j2)):
        if J == 0:
            continue
        for i in range(N):
            j = (i + r) % N; bi, bj = (s >> i) & 1, (s >> j) & 1
            diag += J * np.where(bi == bj, 0.25, -0.25)
            src = np.where(bi != bj)[0]
            R.append(src ^ ((1 << i) | (1 << j))); Cc.append(src); V.append(np.full(len(src), 0.5 * J / scale))
    diag = (diag - shift) / scale
    return (sp.diags(diag) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(Cc))), shape=(dim, dim))).tocsr()


def ed_inputs(N, j1, j2):
    """E0, E1 (lowest level of any symmetry other than the ground state), exact ground state in qiskit ordering (full 2^N vector)."""
    A0, st0 = sector_matrix(N, N // 2, j1, j2); A1, _ = sector_matrix(N, N // 2 + 1, j1, j2)
    w0, v0 = sla.eigsh(A0, k=3, which="SA", tol=1e-13); o = np.argsort(w0); w0, v0 = w0[o], v0[:, o]
    w1 = sla.eigsh(A1, k=1, which="SA", return_eigenvectors=False, tol=1e-13)
    if w0[1] - w0[0] < 1e-8:
        raise ValueError("degenerate ground state")
    psi = np.zeros(1 << N, dtype=complex); psi[st0] = v0[:, 0]
    E0 = float(w0[0]); E1 = float(min(w0[1], w1[0]))
    return dict(E0=E0, E1=E1, psi=psi, E1_singlet=float(w0[1]), E1_triplet=float(w1[0]))


def correlations(vec, N):
    """Translation-averaged <Sz_i Sz_{i+r}>, r = 1..N/2, of a qiskit-ordered state."""
    p = np.abs(vec) ** 2; p = p / p.sum(); s = np.arange(1 << N, dtype=np.int64)
    z = [((s >> i) & 1) - 0.5 for i in range(N)]
    return [float(p @ (sum(z[i] * z[(i + r) % N] for i in range(N)) / N)) for r in range(1, N // 2 + 1)]


def build_ctx(N, j1, j2, chi=96, L=1, verbose=False):
    import quimb.tensor as qtn
    from mps_to_circuit import mps_to_circuit
    from core.trotter import trotter_alpha
    ed = ed_inputs(N, j1, j2)
    W = N * (j1 + j2) / 4.0 - ed["E0"]                    # Weyl: every bond J S.S has lambda_max = J/4; N bonds of each range
    gap = (ed["E1"] - ed["E0"]) / W
    H = pbc_hamiltonian(N, j1, j2); Hsc = scale_hamiltonian(H, ed["E0"], W)
    mps = qtn.MatrixProductState.from_dense(reorder_axes(ed["psi"]), dims=[2] * N)
    mps.compress(cutoff=1e-12, max_bond=chi)
    trial_qc = mps_to_circuit(mps.arrays, method="approximate", shape="lpr", num_layers=L)
    trial_vec = Statevector(trial_qc).data
    gamma = float(abs(np.vdot(ed["psi"], trial_vec)) ** 2)
    return dict(N=N, j1=j1, j2=j2, ed=ed, W=W, gap=gap, H_qk=H, H_scaled=Hsc, trial_qc=trial_qc, trial_vec=trial_vec, gamma=gamma,
                alpha=trotter_alpha(Hsc), energies=np.array([0.0, gap, 1.0]), chi=int(max(mps.bond_sizes())),
                Hs=full_sparse(N, j1, j2, shift=ed["E0"], scale=W))
