"""check_x1_symmetry.py -- the first-order Trotter correction X1 = -(i/2) sum_{a>b}[H_a, H_b] (term order of the circuit) commutes with S^2 and S^z_tot.
Prints ||X1||_F, ||[X1,S^2]||_F and ||[X1,Sz]||_F for the open chain at N = 6 (J2 = 0.4) and N = 8 (J2 = 0.2411).   usage: python trotter2/check_x1_symmetry.py"""
import sys
from pathlib import Path
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent), str(ROOT / "s2proj")]
from hamiltonians import j1j2_hamiltonian, scale_hamiltonian
from s2_ops import s2_sparse
from qiskit.quantum_info import SparsePauliOp
for N, J2 in ((6, 0.4), (8, 0.2411)):
    Hs = scale_hamiltonian(j1j2_hamiltonian(N, j1=1.0, j2=J2), -1.0, 5.0); dim = 1 << N
    mats = [float(np.real(c)) * sp.csr_matrix(SparsePauliOp(l).to_matrix(sparse=True)) for l, c in zip(Hs.paulis.to_labels(), Hs.coeffs) if set(l) != {"I"}]
    L = sp.csr_matrix((dim, dim), dtype=complex); S = sp.csr_matrix((dim, dim), dtype=complex)
    for Ha in mats:
        S = S + (Ha @ L - L @ Ha); L = L + Ha
    X1 = (-0.5j) * S; S2 = s2_sparse(N)
    Sz = SparsePauliOp.from_sparse_list([("Z", [i], 0.5) for i in range(N)], num_qubits=N).to_matrix(sparse=True)
    nrm = sp.linalg.norm
    print(f"N={N} J2={J2}: ||X1||_F = {nrm(X1):.3f}   ||[X1, S^2]||_F = {nrm(X1 @ S2 - S2 @ X1):.1e}   ||[X1, Sz]||_F = {nrm(X1 @ Sz - Sz @ X1):.1e}")
