"""validate_alpha2.py -- brute-force check of the second-order bound  ||S2(t/k)^k - e^{-itH}|| <= t^3 B2 / k^2  on small chains (dense matrices)."""
import sys
from pathlib import Path
import numpy as np, scipy.linalg as sl
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); sys.path.insert(0, str(ROOT / "pbc_corr"))
from alpha2 import alpha2, alpha2_constants
from pbc_ctx import pbc_hamiltonian
from hamiltonians import scale_hamiltonian
from core.trotter import trotter_alpha
from qiskit.quantum_info import SparsePauliOp

worst = 0.0
for (N, J2, pbc) in ((4, 0.0, False), (6, 0.0, False), (6, 0.5, True), (8, 0.2411, False), (8, 0.3, True)):
    H = pbc_hamiltonian(N, 1.0, J2, pbc); Hs = scale_hamiltonian(H, -1.0, 6.0)         # arbitrary shift / width: bound is for the Pauli sum anyway
    terms = [(p, c) for p, c in zip(Hs.paulis, np.real(Hs.coeffs)) if set(p.to_label()) != {"I"}]
    mats = [SparsePauliOp(p, [c]).to_matrix() for p, c in terms]; Hd = sum(mats)
    B2 = alpha2(Hs); a1 = trotter_alpha(Hs)
    for t in (3.0, 10.0):
        Uex = sl.expm(-1j * t * Hd)
        for k in (1, 2, 4, 8, 16):
            half = np.eye(1 << N, dtype=complex)
            for M in mats:
                half = sl.expm(-1j * (t / k) * M / 2) @ half
            half_r = np.eye(1 << N, dtype=complex)
            for M in mats[::-1]:
                half_r = sl.expm(-1j * (t / k) * M / 2) @ half_r
            step = half_r @ half                                                  # S2 = (reversed half) (forward half): palindromic
            err = np.linalg.norm(np.linalg.matrix_power(step, k) - Uex, 2); bnd = t ** 3 * B2 / k ** 2
            e1 = np.linalg.norm(np.linalg.matrix_power(sl.expm(-1j * (t / k) * Hd) * 0 + (lambda s: s)(np.eye(1)), 1), 2) if False else None
            worst = max(worst, err / bnd)
            print(f"N={N} J2={J2} pbc={pbc} t={t:g} k={k:2d}: err {err:.3e}  bound {bnd:.3e}  ratio {err / bnd:.3f}  {'OK' if err <= bnd else 'VIOLATED'}")
print("max err/bound over all cases:", round(worst, 3))
