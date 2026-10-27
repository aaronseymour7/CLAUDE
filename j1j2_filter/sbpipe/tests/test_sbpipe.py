"""Checks of every component against an independent reference.  Run:  python -m pytest sbpipe/tests -q"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import pytest
from scipy.linalg import expm

import sbpipe
from sbpipe import design as D
from sbpipe.pauli import PauliSum, apply_pauli, pauli_sparse
from sbpipe.resources import filter_counts
from sbpipe.trotter import alpha1, beta2, product_formula

_P = {"I": np.eye(2), "X": np.array([[0, 1], [1, 0]]), "Y": np.array([[0, -1j], [1j, 0]]),
      "Z": np.diag([1, -1])}


def kron_label(lab):
    m = np.array([[1.0 + 0j]])
    for c in lab:
        m = np.kron(m, _P[c])
    return m


def random_paulisum(N, n, rng):
    labs = ["".join(rng.choice(list("IXYZ"), N)) for _ in range(n)]
    labs = [l if set(l) != {"I"} else "X" + l[1:] for l in labs]
    return PauliSum(labs, rng.normal(size=n))


def test_pauli_matrices_and_action():
    rng = np.random.default_rng(0)
    for lab in ["XYZ", "IZY", "YYI", "ZIX"]:
        M = kron_label(lab)
        assert np.allclose(pauli_sparse(lab).toarray(), M)
        v = rng.normal(size=8) + 1j * rng.normal(size=8)
        assert np.allclose(apply_pauli(lab, v), M @ v)


def test_weyl_bound_j1j2():
    for N, J2 in [(6, 0.0), (7, 0.4), (8, 0.2411)]:
        H = sbpipe.j1j2_chain(N, 1.0, J2)
        w = np.linalg.eigvalsh(H.to_sparse().toarray())
        lo, hi = sbpipe.weyl_bounds(H)
        assert lo <= w[0] + 1e-12 and hi >= w[-1] - 1e-12
        assert abs(hi - w[-1]) < 1e-9          # attained by the ferromagnet for J1, J2 >= 0


@pytest.mark.parametrize("order", [1, 2])
def test_trotter_bounds_hold(order):
    rng = np.random.default_rng(order)
    cases = [sbpipe.j1j2_chain(5, 1.0, 0.4), random_paulisum(4, 7, rng)]
    worst = 0.0
    for H in cases:
        M = H.to_sparse().toarray()
        c = alpha1(H) if order == 1 else beta2(H)
        for t in [0.05, 0.3, 1.0]:
            for k in [1, 2, 5]:
                U = expm(-1j * t * M)
                S = np.column_stack([product_formula(H, e, t, k, order) for e in np.eye(len(M))])
                err = np.linalg.norm(S - U, 2)
                bnd = c * t ** 2 / (2 * k) if order == 1 else c * t ** 3 / k ** 2
                assert err <= bnd * (1 + 1e-9) + 1e-13
                worst = max(worst, err / bnd)
    assert worst > 0.05                         # bound is not vacuous on these small cases


def test_dyadic_block_is_dirichlet_kernel():
    tau, K = 0.7, 5
    E = np.linspace(0.01, 1, 400)
    F = D.filter_values(tau * 2.0 ** np.arange(K), np.zeros(K), E)
    assert np.allclose(F, np.sin(2 ** K * E * tau) / (2 ** K * np.sin(E * tau)))


@pytest.mark.parametrize("sched", ["uniform", "dyadic"])
def test_design_certificate_vs_dense_grid(sched):
    d = D.design_filter(sched, delta=0.05, w0=1e-4, eps_leak=1e-3, gamma=0.9)
    E = np.linspace(d.delta, 1.0, 2_000_001)
    sup = np.max(np.abs(d.values(E)))
    fpass = np.min(np.abs(d.values(np.linspace(0, d.w0, 101))))
    assert d.f_pass <= fpass + 1e-12
    assert sup / fpass <= d.eta + 1e-12
    assert d.eta <= D.eta_target(1e-3, 0.9)


def test_strang_counts_against_explicit_sequence():
    H = sbpipe.j1j2_chain(4, 1.0, 0.4)
    G = len(H.without_identity())
    for k in [1, 3]:
        c = filter_counts(H, [k], 2)
        # explicit: k steps of forward+backward half sweeps, merge adjacent identical terms
        seq = []
        for _ in range(k):
            seq += list(range(G)) + list(range(G))[::-1]
        merged = [seq[0]]
        for s in seq[1:]:
            if s != merged[-1]:
                merged.append(s)
        assert c["rotations"] == len(merged) + 1      # + merged ancilla Rz
        assert c["cnot"] == 4 * len(merged)           # weight-2 terms (x) Z_anc: 4 CNOT each


def test_cnot_count_against_qiskit():
    qiskit = pytest.importorskip("qiskit")
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    from qiskit.synthesis import LieTrotter

    H = sbpipe.j1j2_chain(4, 1.0, 0.4).without_identity()
    # qiskit labels are little-endian: reverse, and append Z on the ancilla (qubit N = leftmost)
    op = SparsePauliOp(["Z" + l[::-1] for l in H.labels], H.coeffs)
    qc = QuantumCircuit(5)
    qc.append(PauliEvolutionGate(op, time=0.3, synthesis=LieTrotter(reps=1)), range(5))
    ops = transpile(qc, basis_gates=["cx", "rz", "h", "sx", "x"], optimization_level=0).count_ops()
    c = filter_counts(H, [1], 1)
    # the ladder model is an upper bound: qiskit shares CNOTs between consecutive terms of a bond
    # (4 + 2 + 2 per bond instead of 3 x 4)
    assert ops["cx"] <= c["cnot"]
    assert ops["cx"] == 8 * (len(H) // 3)


def test_end_to_end_certificate_holds():
    N = 6
    H = sbpipe.j1j2_chain(N, 1.0, 0.2411)
    w, v = np.linalg.eigh(H.to_sparse().toarray())
    g = v[:, 0]
    rng = np.random.default_rng(3)
    psi = g + 0.25 * rng.normal(size=2 ** N)
    psi /= np.linalg.norm(psi)
    gam = abs(g @ psi) ** 2
    est = sbpipe.EnergyEstimates(E0=w[0], E1=w[1], E0_err=1e-6, E1_err=1e-6)
    for order in (1, 2):
        plan = sbpipe.refine(H, est, sbpipe.Ansatz(gam, psi),
                             sbpipe.Config(eps=1e-2, schedule="dyadic", order=order))
        val = sbpipe.validate(plan, H, sbpipe.Ansatz(gam, psi), g)
        assert val["bound_ok"]
        assert val["infidelity_trotter"] <= 1e-2


def _qsp_phases(p, d, seed=0):
    """Phases phi_0..phi_d with <0| e^{i phi_0 X} prod_j [e^{i theta Z} e^{i phi_j X}] |0> = p(cos theta)."""
    from scipy.optimize import least_squares
    xg = np.cos(np.linspace(0.02, np.pi / 2 - 0.02, 4 * d + 8))

    def blk(ph, th):
        c, s = np.cos(ph), 1j * np.sin(ph)
        U = np.array([[c[0], s[0]], [s[0], c[0]]])
        for j in range(1, len(ph)):
            U = U @ np.diag([np.exp(1j * th), np.exp(-1j * th)]) @ np.array([[c[j], s[j]], [s[j], c[j]]])
        return U[0, 0]

    def res(ph):
        r = np.array([blk(ph, np.arccos(x)) - p(x) for x in xg])
        return np.r_[r.real, r.imag]
    best = None
    for s in range(30):
        ph0 = np.random.default_rng(seed + s).normal(size=d + 1) * 0.5
        o = least_squares(res, (ph0 + ph0[::-1]) / 2)
        if best is None or o.cost < best.cost:
            best = o
        if best.cost < 1e-20:
            break
    return best.x, best.cost


def test_qetu_certificate_and_circuit():
    """The QETU design's closed-form eta and F_pass hold on a dense grid, and an explicit
    ancilla circuit (QSP phases, Trotterised e^{-i tau H (x) Z}) reproduces p(X~) psi."""
    from sbpipe.simulate import qetu_trotter
    for delta, w0, eta_t in [(0.05, 1e-6, 1e-3), (0.2, 1e-3, 1e-2), (0.01, 0.0, 1e-5)]:
        des = D.qetu(delta, w0, eta_t)
        assert des.eta <= eta_t and des.m % 2 == 0
        E = np.linspace(delta, 1.0, 200001)
        assert np.max(np.abs(des.values(E))) <= des.eta * des.f_pass * (1 + 1e-9)
        assert np.min(np.abs(des.values(np.linspace(0, w0, 101)))) >= des.f_pass * (1 - 1e-12)
    # explicit circuit for a small degree on a random 3-qubit Hamiltonian scaled into [0, 1]
    rng = np.random.default_rng(3)
    H = random_paulisum(3, 6, rng)
    w = np.linalg.eigvalsh(H.to_sparse().toarray())
    Hs = PauliSum(H.labels + ["III"], np.r_[H.coeffs / (w[-1] - w[0]), -w[0] / (w[-1] - w[0])])
    d, tau, x_s, k = 6, np.pi / 2, 0.8, 2
    p = lambda x: D.qetu_values(d, 1.0, x_s, np.arccos(x))      # p as a function of x
    ph, cost = _qsp_phases(p, d)
    assert cost < 1e-14
    psi = rng.normal(size=8) + 1j * rng.normal(size=8)
    psi /= np.linalg.norm(psi)
    # ancilla is the first tensor factor; query = S(tau) on |0>, S(-tau) on |1>
    S = lambda t: np.column_stack([product_formula(Hs, e, t, k, 2) for e in np.eye(8)])
    Q = np.kron(np.diag([1, 0]), S(tau)) + np.kron(np.diag([0, 1]), S(-tau))
    rx = lambda a: np.kron(np.array([[np.cos(a), 1j * np.sin(a)], [1j * np.sin(a), np.cos(a)]]), np.eye(8))
    U = rx(ph[0])
    for j in range(1, d + 1):
        U = U @ Q @ rx(ph[j])
    out = (U @ np.kron([1, 0], psi))[:8]
    ref = qetu_trotter(Hs, psi, d, tau, x_s, [k] * d)
    assert np.allclose(out, ref, atol=1e-6)
