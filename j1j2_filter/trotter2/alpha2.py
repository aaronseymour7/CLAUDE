"""alpha2.py -- second-order (Strang) product-formula error constant with commutator scaling, evaluated exactly for local Pauli Hamiltonians.

Childs, Su, Tran, Wiebe, Zhu, PRX 11, 011020 (2021): for H = sum_g H_g and the symmetric formula S2(t) = prod_{g=Gamma..1} e^{-i H_g t/2} prod_{g=1..Gamma} e^{-i H_g t/2},
   || S2(t) - e^{-itH} ||  <=  t^3 [ (1/12) sum_a || [B_a, [B_a, H_a]] ||  +  (1/24) sum_a || [H_a, [H_a, B_a]] || ],   B_a = sum_{b > a} H_b.
(Their bound has two independent inner sums over b > a; they coincide here, so [B_a,[B_a,H_a]].)  Cost: k steps of t/k give error <= t^3 B2 / k^2 with B2 = (1/12) c1 + (1/24) c2.
For Pauli terms H_a = c_a P_a:  [H_a,[H_a,B_a]] = 4 c_a^2 B_a^{anti}  (B^{anti} = later terms anticommuting with P_a), so c2 = sum_a 4 c_a^2 ||B_a^{anti}||.
The orientation of the term order in the circuit is not fixed by the formula, so the constant is taken as the max over the forward and the reversed order.
I verified the formula (as recalled) numerically in validate_alpha2.py before use."""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from qiskit.quantum_info import SparsePauliOp, Pauli
from core.trotter import _strip_identity, _pauli_label_on


def _supp(P):
    return {q for q in range(len(P)) if P.x[q] or P.z[q]}


def _norm_local(op, hermitian=True):
    """Spectral norm of a SparsePauliOp acting on few qubits (restricted to its support)."""
    op = op.simplify(atol=1e-14)
    if len(op) == 0:
        return 0.0
    qs = sorted({q for P in op.paulis for q in _supp(P)})
    if not qs:
        return float(np.max(np.abs(op.coeffs)))
    M = SparsePauliOp([_pauli_label_on(P, qs) for P in op.paulis], op.coeffs).to_matrix(sparse=True)
    if M.shape[0] <= 512:
        return float(np.max(np.abs(np.linalg.eigvalsh(M.toarray()) if hermitian else np.linalg.svd(M.toarray(), compute_uv=False))))
    return float(abs(sla.eigsh(M, k=1, which="LM", return_eigenvectors=False)[0])) if hermitian else float(sla.svds(M, k=1, return_singular_vectors=False)[0])


def alpha2_constants(H, reverse=False):
    """(c1, c2): c1 = sum_a ||[B_a,[B_a,H_a]]||,  c2 = sum_a ||[H_a,[H_a,B_a]]||  (exact, local)."""
    H = _strip_identity(H); P, c = list(H.paulis), np.real(H.coeffs)
    if reverse:
        P, c = P[::-1], c[::-1]
    n = len(P); supp = [_supp(p) for p in P]; c1 = c2 = 0.0
    for a in range(n - 1):
        anti = [b for b in range(a + 1, n) if not P[a].commutes(P[b])]
        if not anti:
            continue
        Banti = SparsePauliOp([P[b] for b in anti], [c[b] for b in anti])
        c2 += 4.0 * c[a] ** 2 * _norm_local(Banti)
        # D = [B_a, H_a] = 2 c_a sum_{b anti} c_b P_b P_a  (a Pauli sum), then E = [B_a, D]
        D = SparsePauliOp([P[b].dot(P[a]) for b in anti], [2.0 * c[a] * c[b] for b in anti]).simplify(atol=1e-14)
        sD = set().union(*[_supp(p) for p in D.paulis]) if len(D) else set()
        later = [b for b in range(a + 1, n) if supp[b] & sD]
        if not later:
            continue
        Bloc = SparsePauliOp([P[b] for b in later], [c[b] for b in later])
        E = (Bloc.dot(D) - D.dot(Bloc)).simplify(atol=1e-14)
        c1 += _norm_local(E)
    return c1, c2


def alpha2(H):
    """B2 = (1/12) c1 + (1/24) c2, max over the two term orders:  error of k steps of the Strang formula over total time t  <=  t^3 B2 / k^2."""
    best = 0.0
    for rev in (False, True):
        c1, c2 = alpha2_constants(H, reverse=rev); best = max(best, c1 / 12.0 + c2 / 24.0)
    return best
