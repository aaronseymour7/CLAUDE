"""Checks of the QPF engine against the paper's definitions and its 4 x 4 Heisenberg numbers.
Run from ORNL/: python -m pytest tests -q"""
import os
import sys

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qpf import core, design, models, symmetry  # noqa: E402

Y = np.array([[0, -1j], [1j, 0]])


def random_hermitian(d, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    return (A + A.conj().T) / 2


def test_ancilla_step_is_paper_eq1():
    """exp(-i (tO + delta) (x) Y_a) |psi>|0> = cos(tO+delta)|psi>|0> + sin(tO+delta)|psi>|1>."""
    d, t, delta = 6, 0.7, 0.3
    O = random_hermitian(d)
    psi = np.random.default_rng(1).normal(size=d) + 0j
    U = expm(-1j * np.kron(t * O + delta * np.eye(d), Y))          # system (x) ancilla
    out = U @ np.kron(psi, [1, 0])
    a0, a1 = core.ancilla_step(O, psi, t, delta)
    assert np.allclose(out[0::2], a0) and np.allclose(out[1::2], a1)
    w, V = np.linalg.eigh(O)
    assert np.allclose(a0, V @ (np.cos(t * w + delta) * (V.conj().T @ psi)))


def test_run_is_product_of_cosines_and_probabilities():
    O = random_hermitian(8, 2)
    psi = np.random.default_rng(3).normal(size=8) + 0j
    psi /= np.linalg.norm(psi)
    times, phases = [0.4, 1.1, 2.3], [0.1, -0.2, 0.05]
    r = core.run(O, psi, times, phases)
    w, V = np.linalg.eigh(O)
    c = V.conj().T @ psi
    out = V @ (core.filter_function(w, times, phases) * c)
    assert abs(r.p_total - np.vdot(out, out).real) < 1e-12
    assert abs(abs(np.vdot(out / np.linalg.norm(out), r.state)) - 1) < 1e-12
    d = core.run(core.Diagonal(w), c, times, phases)                # eigenbasis form used for design
    assert abs(d.p_total - r.p_total) < 1e-12


def test_converged_success_probability():
    """P = f(0)^2 |<Psi0|psi>|^2 once every excited component is gone (O with a large gap, many steps)."""
    w = np.array([0.0, 0.5, 0.7, 1.0])
    c = np.array([0.6, 0.5, 0.4, np.sqrt(1 - 0.36 - 0.25 - 0.16)])
    times, phases = [np.pi / 1.4] * 40, [0.05] * 40
    r = core.run(core.Diagonal(w), c, times, phases)
    assert abs(r.state[0]) > 1 - 1e-9
    assert abs(r.p_total - core.success_probability_limit(0.36, phases)) < 1e-9


def test_halving_ladder():
    assert symmetry.halving_times(1.0, 3) == [1.0, 0.5, 0.25]
    assert symmetry.removes([2, 6, 30, 42], [np.pi / 4]).all()       # J = 1, 2, 5, 6 with t1 = pi/4
    assert np.allclose(symmetry.j2_ladder(0, 8), [np.pi / 4, np.pi / 8, np.pi / 16])   # paper: three steps
    ts = symmetry.j2_ladder(3, 6)                                     # J^2 - 12 (paper's 10B)
    assert symmetry.removes(symmetry.j2_values(3, 6), ts).all()


def test_j2_projection_small_lattice():
    sec = models.SpinSector(8)
    J2 = models.total_spin_sq(sec)
    w, V = np.linalg.eigh(J2.toarray())
    assert np.allclose(w, np.round(w)) and set(np.round(w)) <= {J * (J + 1) for J in range(5)}
    psi = models.neel(sec, 4, 2)
    r = core.run(J2, psi, symmetry.j2_ladder(0, 4))
    P0 = V[:, np.round(w) == 0]
    assert abs(r.p_total - np.sum(abs(P0.T @ psi) ** 2)) < 1e-10      # success = J = 0 weight
    assert np.linalg.norm(J2 @ r.state) < 1e-9


@pytest.fixture(scope="module")
def heis4x4():
    sec = models.SpinSector(16)
    H = models.heisenberg(sec, models.square_bonds(4, 4))
    J2 = models.total_spin_sq(sec)
    w, V = eigsh(H, k=2, which="SA", tol=1e-12)
    o = np.argsort(w)
    Emax = eigsh(H, k=1, which="LA", return_eigenvectors=False)[0]
    return sec, H, J2, w[o], V[:, o], Emax


def test_paper_heisenberg_numbers(heis4x4):
    sec, H, J2, w, V, Emax = heis4x4
    W = Emax - w[0]
    assert len(models.square_bonds(4, 4)) == 32 and sec.dim == 12870      # paper: 12,870 states with J_z = 0
    Hs = core.scale(H, w[0], Emax)
    psi = models.neel(sec, 4, 4)
    assert abs(core.Operator(Hs).expect(psi) - 0.17) < 0.005              # paper: scaled energy ~0.17
    assert abs((w[1] - w[0]) / W - 0.03) < 0.002                           # paper: gap ~0.03
    ws = np.sort(eigsh(H + W * J2, k=2, which="SA", return_eigenvectors=False))
    assert abs((ws[1] - ws[0]) / W - 0.15) < 0.01                         # paper: ~0.15 after J^2 projection
    r = core.run(J2, psi, symmetry.j2_ladder(0, 8))
    assert abs(r.p_total - 0.11) < 0.005                                   # paper: ~11%
    assert abs(core.Operator(Hs).expect(r.state) - 0.06) < 0.005          # paper: E ~0.06


def test_design_reproduces_authors_notebook():
    """reference/optimize_times_phases.ipynb, its example cell (method v3, 8-10 pulses): fun = 1.319e-01,
    fdiff = 2.498529e-03 at overlap 0.34, total time pi/gap = 8.183140."""
    wmin, wmax = -3.90983149, 7.92112829
    gap = (0.63219691 - wmin) / abs(wmax - wmin)
    for m in (8, 9, 10):
        r = design.optimize(gap, m, method="v3")
        assert abs(r["T"] - 8.183140) < 1e-6 and abs(r["times"].sum() - r["T"]) < 1e-8
        assert abs(r["cost"] - 1.319e-01) < 5e-4
        fd, _ = design.evaluate(r["times"], r["phases"], r["E"], design.notebook_trial(0.34, 1000))
        assert abs(fd - 2.498529e-03) < 1e-8


def test_design_paper_properties():
    q = design.optimize(0.15, 6, method="v3")
    z = design.optimize(0.15, 6, method="v4")
    assert q["cost"] < z["cost"]                                           # phases help
    assert q["f0"] ** 2 >= 0.9 - 1e-9                                      # success kept within 10%
    assert (q["phases"] > 0).all()                                         # paper: phases all positive
    o = np.argsort(q["times"])
    assert (np.diff(q["phases"][o]) > 0).all()                             # ... and larger for longer times
    assert q["times"].max() <= q["T"] / 3 + 1e-9
    E = design.baroni_spectrum(0.15)
    assert design.baroni_cost(*design.constant(np.pi / 0.15, 6), E) > z["cost"]   # optimised beats constant
    g, _ = design.gaussian(np.pi / 0.15, 7, seed=3)
    assert abs(np.abs(g).sum() - np.pi / 0.15) < 1e-9
