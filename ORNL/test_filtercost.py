"""Checks of the engine. Run with `python test_filtercost.py` (or pytest). Takes a few seconds."""
import numpy as np
from scipy.linalg import expm

import filtercost as fc


def test_strang_is_second_order():
    ch = fc.Chain(8, 0.1)
    Hs = (ch.hamiltonian().toarray() - ch.E0 * np.eye(ch.dim)) / ch.W
    v = ch.trial(0.3)[:, None]
    exact = expm(-0.7j * Hs) @ v
    err = [np.linalg.norm(ch.strang(v, 0.7 / k, k) - exact) for k in (10, 100)]
    assert 90 < err[0] / err[1] < 110                                  # error ~ dt^2


def test_eigenphase_filter_equals_circuit():
    ch = fc.Chain(8, -0.1)
    dt = 0.9
    phi, project = fc.step_spectrum(ch, dt)
    v = ch.trial(0.25)
    for m in (1, 6, 37):
        circuit = 0.5 * (ch.strang(v[:, None], dt, m) + ch.strang(v[:, None], -dt, m))[:, 0]
        eig = np.cos(m * phi) * project(v)
        assert abs(np.vdot(ch.g0, circuit) - project(ch.g0).conj() @ eig) < 1e-12
        assert abs(np.linalg.norm(circuit) ** 2 - np.sum(abs(eig) ** 2)) < 1e-12


def test_knobs_hit_targets():
    d = fc.delta_for_gap(10, 0.1)
    ch = fc.Chain(10, d)
    assert abs(ch.gap - 0.1) < 1e-8
    assert abs(ch.overlap(ch.eps_for_overlap(0.3)) - 0.3) < 1e-10


def test_cnot_count():
    # 6-CNOT controlled bond blocks, outer Strang half-steps merged (same rule as sbpipe.resources.filter_counts)
    assert fc.cnot_per_attempt(10, 2, 1, 1) == 588
    assert fc.cnot_per_attempt(10, 8, 2, 3) == 293856


def test_aa_cost_matches_grover():
    p = 0.02
    th = np.arcsin(np.sqrt(p))
    k = int(round(np.pi / (4 * th) - 0.5))
    R = np.array([[np.cos(2 * th), -np.sin(2 * th)], [np.sin(2 * th), np.cos(2 * th)]])   # one Grover iterate
    s = np.linalg.matrix_power(R, k) @ np.array([np.cos(th), np.sin(th)])
    assert abs(fc.aa_cost(1.0, p) - (2 * k + 1) / s[1] ** 2) < 1e-12
    assert fc.aa_cost(1.0, p) < 1 / p


def test_measured_point_meets_target():
    ch = fc.Chain(8, 0.0)
    rows, _ = fc.measure_point(ch, [0.1, 0.9], [1e-2, 1e-3])
    for r in rows:
        for e, v in r.items():
            for pr in fc.PROTOCOLS:
                assert v[pr]["infid"] <= e
            assert v["aa"]["cost"] <= v["rus"]["cost"] * (1 + 1e-12)


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f(); print("ok", name)
