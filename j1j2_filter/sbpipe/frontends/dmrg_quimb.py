"""DMRG front end (quimb): PauliSum -> MPO -> E0, E1, reference state, low-bond-dimension trial ansatz.

This is one possible producer of the pipeline inputs (EnergyEstimates, Ansatz); nothing in the
pipeline depends on it.

* E0 and a reference MPS psi_ref from two-site DMRG with increasing bond dimension.
* E1 as the lowest state in a symmetry sector chosen by a penalty mu (Q - q)^2, Q a conserved
  charge given as a PauliSum (for spin chains: total S^z, q = 1, which selects the lowest
  triplet). This is valid as the gap only if the first excitation lies in that sector; the
  study checks it with exact diagonalisation where possible.
* Ansatz: the variationally optimal MPS at a small bond dimension chi (DMRG with chi fixed).
  gamma = |<psi_ref|psi_chi>|^2 by MPS contraction. An MPS of bond dimension 2^b is exactly a
  staircase of (b+1)-qubit unitaries, which is the circuit cost model used here.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
import quimb.tensor as qtn

from ..inputs import Ansatz, EnergyEstimates
from ..pauli import PauliSum

_P = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def _product_mpo(label: str, coeff: float) -> qtn.MatrixProductOperator:
    N = len(label)
    arrs = []
    for q, c in enumerate(label):
        m = _P[c] * (coeff if q == 0 else 1.0)
        if q == 0 or q == N - 1:
            arrs.append(m.reshape(1, 2, 2))
        else:
            arrs.append(m.reshape(1, 1, 2, 2))
    return qtn.MatrixProductOperator(arrs, shape="lrud")


def pauli_to_mpo(H: PauliSum, chunk: int = 16, cutoff: float = 1e-13) -> qtn.MatrixProductOperator:
    """Generic MPO of a PauliSum by summing product MPOs and compressing every `chunk` terms."""
    mpo = None
    for n, (lab, c) in enumerate(zip(H.labels, H.coeffs)):
        t = _product_mpo(lab, float(c))
        mpo = t if mpo is None else mpo + t
        if (n + 1) % chunk == 0:
            mpo.compress(cutoff=cutoff)
    mpo.compress(cutoff=cutoff)
    return mpo


def total_sz(N: int) -> PauliSum:
    labels = ["I" * q + "Z" + "I" * (N - q - 1) for q in range(N)]
    return PauliSum(labels, np.full(N, 0.5))


def _dmrg(mpo, bond_dims, cutoff, tol, max_sweeps, p0=None):
    d = qtn.DMRG2(mpo, bond_dims=bond_dims, cutoffs=cutoff, p0=p0)
    d.solve(tol=tol, max_sweeps=max_sweeps, verbosity=0)
    return float(np.real(d.energy)), d.state.copy(), d


def ground_state(H_mpo, bond_dims=(16, 32, 64, 128, 256), cutoff=1e-12, tol=1e-11, max_sweeps=60):
    """E0, psi_ref and a convergence-based error bar |E(chi_max) - E(chi_prev)| (heuristic)."""
    E_prev, psi = None, None
    hist = []
    for chi in bond_dims:
        E, psi, d = _dmrg(H_mpo, [chi], cutoff, tol, max_sweeps, p0=psi)
        hist.append((chi, E, int(psi.max_bond())))
        E_prev = E
    err = abs(hist[-1][1] - hist[-2][1]) if len(hist) > 1 else float("nan")
    return dict(E0=hist[-1][1], psi=psi, E0_conv=err, history=hist)


def sector_ground_state(H_mpo, Q: PauliSum, q: float, mu: float, **kw):
    """Lowest state of H + mu (Q - q)^2, i.e. the lowest state with charge q if mu exceeds the
    energy gained by leaving the sector."""
    N = Q.num_qubits
    Qm = pauli_to_mpo(Q)
    Q2 = Qm.apply(Qm)
    pen = Q2 + (-2.0 * q) * Qm + (q * q) * qtn.MPO_identity(N, phys_dim=2)
    Hp = H_mpo + mu * pen
    Hp.compress(cutoff=1e-13)
    return ground_state(Hp, **kw)


def trial_mps(H_mpo, chi: int, tol=1e-10, max_sweeps=60, seed: int = 0):
    np.random.seed(seed)
    E, psi, _ = _dmrg(H_mpo, [chi], 1e-14, tol, max_sweeps)
    return E, psi


def mps_circuit_cost(N: int, chi: int) -> Dict[str, float]:
    """Exact sequential (staircase) circuit of an MPS with bond dimension chi = 2^b:
    N - b unitaries on b+1 qubits. Generic n-qubit unitary: CNOT <= (23/48) 4^n - (3/2) 2^n + 4/3
    (Shende-Bullock-Markov quantum Shannon decomposition; 3, 20, 100 for n = 2, 3, 4) and 4^n - 1
    rotation parameters. chi = 1 is a product state: N single-qubit gates, no CNOT."""
    b = int(np.ceil(np.log2(max(chi, 1))))
    if b == 0:
        return dict(cnot=0, rotations=3 * N, unitaries=N, width=1)
    n = b + 1
    cnot_u = round((23 / 48) * 4 ** n - 1.5 * 2 ** n + 4 / 3)
    return dict(cnot=(N - b) * cnot_u, rotations=(N - b) * (4 ** n - 1), unitaries=N - b, width=n)


def mps_overlap(a, b) -> complex:
    """<a|b> by a left-to-right transfer-matrix sweep (no large einsum; any length)."""
    a, b = a.copy(), b.copy()
    a.permute_arrays("lrp")
    b.permute_arrays("lrp")
    A = [np.asarray(x) for x in a.arrays]
    B = [np.asarray(x) for x in b.arrays]
    env = np.einsum("rp,sp->rs", A[0].conj(), B[0])
    for x, y in zip(A[1:-1], B[1:-1]):
        env = np.einsum("rs,rRp,sSp->RS", env, x.conj(), y)
    return complex(np.einsum("rs,rp,sp->", env, A[-1].conj(), B[-1]))


def ansatz_from_dmrg(H_mpo, psi_ref, chi: int, N: int, dense_max_n: int = 16) -> Ansatz:
    """Optimal bond-dimension-chi MPS as an Ansatz; gamma is its overlap with psi_ref."""
    E_t, psi_t = trial_mps(H_mpo, chi)
    gamma = float(abs(mps_overlap(psi_ref, psi_t)) ** 2
                  / (abs(mps_overlap(psi_ref, psi_ref)) * abs(mps_overlap(psi_t, psi_t))))
    sv = np.asarray(psi_t.to_dense()).reshape(-1) if N <= dense_max_n else None
    return Ansatz(gamma=gamma, statevector=sv, cost=mps_circuit_cost(N, chi),
                  name=f"DMRG MPS chi={chi}",
                  meta=dict(chi=chi, E_trial=E_t, max_bond=int(psi_t.max_bond())))


def dmrg_inputs(H: PauliSum, chi_trial=2, Q: Optional[PauliSum] = None, q: float = 1.0,
                mu: float = 2.0, bond_dims=(16, 32, 64, 128, 256), E0_err_floor: float = 1e-10):
    """Run the front end and return (EnergyEstimates, Ansatz or list of Ansatz, extras).

    E0_err / E1_err are the change between the last two bond dimensions (floored at
    E0_err_floor): a convergence heuristic, not a rigorous bound."""
    N = H.num_qubits
    H_mpo = pauli_to_mpo(H)
    g = ground_state(H_mpo, bond_dims=bond_dims)
    Q = Q if Q is not None else total_sz(N)
    e = sector_ground_state(H_mpo, Q, q, mu, bond_dims=bond_dims)
    est = EnergyEstimates(E0=g["E0"], E1=e["E0"], E0_err=max(g["E0_conv"], E0_err_floor),
                          E1_err=max(e["E0_conv"], E0_err_floor),
                          source=f"DMRG chi<={max(bond_dims)}; E1 = lowest state with Q={q:g}")
    chis = [chi_trial] if np.isscalar(chi_trial) else list(chi_trial)
    ans = [ansatz_from_dmrg(H_mpo, g["psi"], c, N) for c in chis]
    extras = dict(psi_ref=g["psi"], ref_history=g["history"], E1_history=e["history"], H_mpo=H_mpo)
    return est, (ans[0] if np.isscalar(chi_trial) else ans), extras
