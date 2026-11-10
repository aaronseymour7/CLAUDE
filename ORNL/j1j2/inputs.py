"""Inputs for running the QPF pipeline on the open J1-J2 Heisenberg chain.

H = J1 sum_i S_i.S_{i+1} + J2 sum_i S_i.S_{i+2}, open boundaries (no dimerisation).

Everything here is in the full 2^N space in Qiskit's ordering (bit i of a basis index = site i = qubit i), so the
trial statevector from the mps-to-circuit circuit can be used as is.

  operators     H and J^2 as sparse matrices (S_i.S_j = (2 SWAP_ij - 1)/4, J^2 = sum_{i<j} SWAP_ij + 3N/4 - N(N-1)/4)
  dmrg_estimates DMRG (quimb) spectral estimates: E0 from H, Emax from -H, and E1 = the lowest excited singlet,
                from DMRG on H + lam J^2 + 2 lam |psi0><psi0| (lam = Emax - E0): the J^2 penalty pushes every
                J > 0 level above the singlets and the projector on the (compressed) DMRG ground state pushes it up,
                so the lowest level left is the first excited singlet -- the gap the energy filter faces after the
                J^2 projection.
  trial         the default mps-to-circuit circuit (method "approximate", one layer) for the DMRG ground state,
                and its statevector.
  exact         ED reference (only for measuring fidelity and the error of the estimates).
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def bonds(N, J1=1.0, J2=0.2411):
    return [(i, i + 1, J1) for i in range(N - 1)] + [(i, i + 2, J2) for i in range(N - 2)]


def _swap(N, i, j):
    s = np.arange(1 << N)
    flip = ((s >> i) & 1) != ((s >> j) & 1)
    return np.where(flip, s ^ (1 << i) ^ (1 << j), s)


def _swap_sum(N, pairs, weights):
    d = 1 << N
    rows = np.concatenate([_swap(N, i, j) for i, j in pairs])
    cols = np.tile(np.arange(d), len(pairs))
    vals = np.concatenate([np.full(d, w, float) for w in weights])
    return sp.csr_matrix((vals, (rows, cols)), shape=(d, d))


def hamiltonian(N, J1=1.0, J2=0.2411):
    b = bonds(N, J1, J2)
    H = _swap_sum(N, [(i, j) for i, j, _ in b], [J / 2 for *_, J in b])
    return (H - 0.25 * sum(J for *_, J in b) * sp.identity(1 << N, format="csr")).tocsr()


def total_spin_sq(N):
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    S = _swap_sum(N, pairs, np.ones(len(pairs)))
    return (S + (0.75 * N - 0.25 * N * (N - 1)) * sp.identity(1 << N, format="csr")).tocsr()


def exact(N, J1=1.0, J2=0.2411):
    """ED reference: E0, ground state, Emax, lowest triplet-or-higher E1 and lowest excited singlet."""
    H = hamiltonian(N, J1, J2)
    w, V = eigsh(H, k=2, which="SA", tol=1e-12)
    o = np.argsort(w)
    Emax = float(eigsh(H, k=1, which="LA", tol=1e-10, return_eigenvectors=False)[0])
    W = Emax - w[o[0]]
    ws = np.sort(eigsh(H + W * total_spin_sq(N), k=2, which="SA", tol=1e-12, return_eigenvectors=False))
    return dict(E0=float(w[o[0]]), g0=V[:, o[0]], E1=float(w[o[1]]), E1_singlet=float(ws[1]), Emax=Emax, H=H)


# ---------------------------------------------------------------- DMRG
def _builder(N, J1, J2, sign=1.0, s2=0.0):
    import quimb.operator as qop
    E = {}
    for i, j, J in bonds(N, J1, J2):
        E[(i, j)] = E.get((i, j), 0.0) + sign * J
    if s2:
        for i in range(N):
            for j in range(i + 1, N):
                E[(i, j)] = E.get((i, j), 0.0) + 2 * s2          # s2 (J^2 - 3N/4) = s2 sum_{i<j} 2 S_i.S_j
    return qop.heisenberg_from_edges(list(E), j=lambda i, j: E[(i, j)] if (i, j) in E else E[(j, i)])


def _seeded_mps(N, seed, D=8):
    import quimb.tensor as qtn
    rng = np.random.default_rng(seed)
    arrs = [rng.normal(size=(D, 2))] + [rng.normal(size=(D, D, 2)) for _ in range(N - 2)] + [rng.normal(size=(D, 2))]
    p0 = qtn.MatrixProductState(arrs, shape="lrp")
    return p0 / p0.norm()


def _dmrg(mpo, N, seed, bond_dims):
    import quimb.tensor as qtn
    d = qtn.DMRG2(mpo, bond_dims=list(bond_dims), cutoffs=1e-12, p0=_seeded_mps(N, seed))
    np.random.seed(seed)
    d.solve(tol=1e-10, verbosity=0, max_sweeps=40)
    return float(np.real(d.energy)), d.state


def _projector_mpo(psi, max_bond):
    import quimb.tensor as qtn
    p = psi.copy()
    p.compress(max_bond=max_bond)
    p /= p.norm()
    p.permute_arrays("lrp")
    arrs = []
    n = p.L
    for k, A in enumerate(p.arrays):
        if k in (0, n - 1):
            W = np.einsum("ru,sd->rsud", A, A.conj()).reshape(A.shape[0] ** 2, 2, 2)
        else:
            W = np.einsum("lru,msd->lmrsud", A, A.conj()).reshape(A.shape[0] ** 2, A.shape[1] ** 2, 2, 2)
        arrs.append(W)
    return qtn.MatrixProductOperator(arrs, shape="lrud")


def dmrg_estimates(N, J1=1.0, J2=0.2411, seed=0, bond_dims=(16, 32, 64, 128), projector_bond=16):
    """DMRG E0, Emax, E1 (lowest excited singlet) and the DMRG ground-state MPS."""
    E0, psi0 = _dmrg(_builder(N, J1, J2).build_mpo(), N, seed, bond_dims)
    Emax = -_dmrg(_builder(N, J1, J2, sign=-1.0).build_mpo(), N, seed, bond_dims)[0]
    lam = Emax - E0
    Hpen = _builder(N, J1, J2, s2=lam).build_mpo().add_MPO(_projector_mpo(psi0, projector_bond) * (2 * lam))
    E1p, _ = _dmrg(Hpen, N, seed, bond_dims)
    E1 = E1p + lam * 0.75 * N                                     # remove the constant -lam 3N/4 of the J^2 penalty
    return dict(E0=E0, Emax=Emax, E1_singlet=E1, psi0=psi0)


# ---------------------------------------------------------------- trial state
def mps_statevector(psi):
    """Dense vector of a quimb MPS in Qiskit ordering (site i = bit i)."""
    N = psi.L
    v = np.asarray(psi.to_dense()).ravel()                         # quimb: site 0 = most significant
    return v.reshape([2] * N).transpose(list(range(N))[::-1]).ravel()


def trial(psi0, num_layers=1):
    """The default mps-to-circuit circuit for psi0 and its statevector (Qiskit ordering)."""
    from mps_to_circuit import mps_to_circuit
    from qiskit.quantum_info import Statevector
    mps = psi0.copy()
    mps.permute_arrays("lrp")
    qc = mps_to_circuit(list(mps.arrays), method="approximate", shape="lrp", num_layers=num_layers)
    return qc, np.asarray(Statevector(qc).data)
