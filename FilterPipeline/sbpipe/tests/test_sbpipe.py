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


def test_operator_norm_never_underestimates_up_to_10_qubits():
    """Dense bound (with its 1e-11 margin) for supports of 9 and 10 qubits against sparse Lanczos on the
    same matrix, and the triangle fallback beyond the limit."""
    from scipy.sparse.linalg import eigsh

    from sbpipe.pauli import local_matrix, operator_norm
    rng = np.random.default_rng(3)
    for nq in (9, 10):
        terms = {}
        while len({q for l in terms for q, c in enumerate(l) if c != "I"}) < nq:
            lab = "".join(rng.choice(list("IXYZ"), nq))
            if set(lab) != {"I"}:
                terms[lab] = float(rng.normal())
        qs = list(range(nq))
        lam = abs(eigsh(local_matrix(terms, qs, sparse=True), k=1, which="LM",
                        return_eigenvectors=False, tol=1e-13)[0])
        b, exact = operator_norm(terms)
        assert exact and lam <= b <= lam * (1 + 1e-9)
    big = {"X" * 11: 1.0, "Z" * 11: 0.5}
    b, exact = operator_norm(big)
    assert not exact and b == pytest.approx(1.5)


def test_eckart_gamma_is_a_lower_bound():
    H = sbpipe.j1j2_chain(8, 1.0, 0.2411)
    w, v = np.linalg.eigh(H.to_sparse().toarray())
    est = sbpipe.EnergyEstimates(E0=w[0], E1=w[np.argmax(w > w[0] + 1e-8)])
    rng = np.random.default_rng(5)
    for mix in (0.02, 0.1, 0.3):
        psi = v[:, 0] + mix * (v[:, 1:9] @ rng.normal(size=8))
        psi /= np.linalg.norm(psi)
        g_true = abs(np.vdot(v[:, 0], psi)) ** 2
        E_t = float(np.vdot(psi, H.to_sparse() @ psi).real)
        g_eck = sbpipe.eckart_gamma(E_t, est)
        assert g_eck <= g_true + 1e-12
        assert g_eck > 0.0
    # a margin on the gap lowers the bound further
    est2 = sbpipe.EnergyEstimates(E0=est.E0, E1=est.E1, E0_err=1e-3, E1_err=1e-3)
    assert sbpipe.eckart_gamma(E_t, est2) < g_eck


def test_heisenberg_ladder_spectrum_and_weyl():
    H = sbpipe.heisenberg_ladder(3, 1.0, 0.7)
    w = np.linalg.eigvalsh(H.to_sparse().toarray())
    lo, hi = sbpipe.weyl_bounds(H)
    assert lo <= w[0] + 1e-12 and abs(hi - w[-1]) < 1e-9
    # a 2-rung ladder is a 4-site ring, E0 = -(2 J_leg + ...) check against dense diagonalisation of S.S terms
    H2 = sbpipe.heisenberg_ladder(2, 1.0, 1.0)
    w2 = np.linalg.eigvalsh(H2.to_sparse().toarray())
    assert w2[0] == pytest.approx(-2.0)          # 4-site Heisenberg ring: E0 = -2 J


def test_qsvt_degree_scales_with_inverse_gap():
    from sbpipe.qsvt import qsvt_filter_cost
    H = sbpipe.j1j2_chain(8, 1.0, 0.2411)
    a = qsvt_filter_cost(H, 0.2, 0.5, 1e-3, 1e-2, {})
    b = qsvt_filter_cost(H, 0.1, 0.5, 1e-3, 1e-2, {})
    assert b["degree"] == pytest.approx(2 * a["degree"], rel=0.02)
    assert a["t_success"] > a["t_attempt"] > 0


# ----------------------------------------------------------------------------- J^2 projection stage
def _sa2_dense(N, axis):
    S = sum(0.5 * kron_label("I" * q + axis + "I" * (N - q - 1)) for q in range(N))
    return S @ S


def test_spin_projection_ladder_is_exact():
    from sbpipe.symmetry import j2_max, projection_times, projection_values
    for N in (2, 4, 6, 8, 12, 24, 48):
        t = projection_times(N)
        R = len(t)
        assert 2 ** (R + 1) > j2_max(N) and (R == 1 or 2 ** R <= j2_max(N))
        v = projection_values(N)
        assert v[0] == pytest.approx(1.0)
        assert max(abs(x) for j, x in v.items() if j > 0) < 1e-9
    assert len(projection_times(48)) == 9


def test_spin_matrices_casimir_and_commutator():
    from sbpipe.symmetry import spin_matrices
    for j in (0, 1, 2, 5):
        Sx, Sy, Sz = spin_matrices(j)
        assert np.allclose(Sx @ Sx + Sy @ Sy + Sz @ Sz, j * (j + 1) * np.eye(2 * j + 1))
        assert np.allclose(Sx @ Sy - Sy @ Sx, 1j * Sz)


def test_evolve_sa2_matches_dense_expm():
    from sbpipe.symmetry import evolve_sa2
    rng = np.random.default_rng(11)
    N = 4
    psi = rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N)
    for ax in "xyz":
        M = _sa2_dense(N, ax.upper())
        assert np.allclose(evolve_sa2(psi, N, ax, 0.37), expm(-0.37j * M) @ psi)


@pytest.mark.parametrize("order", [1, 2])
def test_spin_trotter_constants_bound_the_error(order):
    from sbpipe.symmetry import spin_product_formula, spin_trotter_constants
    N = 4
    a1, b2 = spin_trotter_constants(N)
    J2 = sum(_sa2_dense(N, a) for a in "XYZ")
    for t in (0.1, 0.4, 0.8):
        for k in (1, 3):
            cols = [spin_product_formula(np.eye(2 ** N)[:, i].astype(complex), N, t, k, order)
                    for i in range(2 ** N)]
            err = np.linalg.norm(np.array(cols).T - expm(-1j * t * J2), 2)
            bound = a1 * t ** 2 / (2 * k) if order == 1 else b2 * t ** 3 / k ** 2
            assert err <= bound + 1e-12


def test_spin_stage_end_to_end_certificate():
    from sbpipe.symmetry import refine_spin, validate_spin
    N = 6
    H = sbpipe.j1j2_chain(N, 1.0, 0.2411)
    Hd = H.to_sparse().toarray()
    from sbpipe.symmetry import j2_sparse
    J2 = j2_sparse(N).toarray()
    w, v = np.linalg.eigh(Hd + 2.0 * J2)            # singlet excitation = second level of H + 2 J^2
    E0, E1s = np.linalg.eigvalsh(Hd)[0], w[1]
    g = v[:, 0]
    rng = np.random.default_rng(2)
    psi = g + 0.4 * rng.normal(size=2 ** N)         # contains triplets and higher spin components
    psi /= np.linalg.norm(psi)
    gamma = float(abs(np.vdot(g, psi)) ** 2)
    est = sbpipe.EnergyEstimates(E0=E0, E1=E1s)
    ans = sbpipe.Ansatz(gamma=gamma, statevector=psi, cost={"cnot": 10, "rotations": 20})
    plan = refine_spin(H, est, ans, sbpipe.Config(eps=1e-2, schedule="dyadic", order=2))
    out = validate_spin(plan, H, ans, g)
    assert out["bound_ok"]
    assert out["j2_trial"] > 0.5 and out["j2_after_projection"] < 1e-4
    plain = sbpipe.refine(H, sbpipe.EnergyEstimates(E0=E0, E1=np.linalg.eigvalsh(Hd)[1]), ans,
                          sbpipe.Config(eps=1e-2, schedule="dyadic", order=2))
    assert plan.design.T < plain.design.T            # larger gap, shorter energy filter
