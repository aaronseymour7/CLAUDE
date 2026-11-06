"""Engine checks (numpy/scipy only). Run from ORNL/: python -m pytest tests -q"""
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ornl import Chain, StudyConfig, chain_for_gap, costs, measure_point, trotter  # noqa: E402


def test_strang_is_second_order():
    ch = Chain(8, 0.1)
    Hs = (ch.hamiltonian().toarray() - ch.E0 * np.eye(ch.dim)) / ch.W
    v = ch.trial(0.3)
    exact = expm(-0.7j * Hs) @ v
    err = [np.linalg.norm(trotter.strang(ch, v, 0.7 / k, k) - exact) for k in (10, 100)]
    assert 90 < err[0] / err[1] < 110                                  # error ~ dt^2


def test_eigenphase_filter_equals_gate_by_gate():
    ch = Chain(8, -0.1)
    dt = 0.9
    phi, project = trotter.step_spectrum(ch, dt)
    v = ch.trial(0.25)
    for m in (1, 6, 37):
        gates = trotter.pulse(ch, v, dt, m)
        eig = np.cos(m * phi) * project(v)
        assert abs(np.vdot(ch.g0, gates) - project(ch.g0).conj() @ eig) < 1e-12
        assert abs(np.linalg.norm(gates) ** 2 - np.sum(abs(eig) ** 2)) < 1e-12


def test_knobs_hit_targets():
    ch = chain_for_gap(10, 0.1, StudyConfig())
    assert abs(ch.gap - 0.1) < 1e-8
    assert abs(ch.overlap(ch.eps_for_overlap(0.3)) - 0.3) < 1e-10


def test_cnot_count():
    assert costs.cnot_per_attempt(10, 2, 1, 1) == 588                 # same rule as sbpipe.resources
    assert costs.cnot_per_attempt(10, 8, 2, 3) == 293856


def test_aa_cost_matches_grover():
    p = 0.02
    th = np.arcsin(np.sqrt(p))
    k = int(round(np.pi / (4 * th) - 0.5))
    R = np.array([[np.cos(2 * th), -np.sin(2 * th)], [np.sin(2 * th), np.cos(2 * th)]])
    s = np.linalg.matrix_power(R, k) @ np.array([np.cos(th), np.sin(th)])
    assert abs(costs.aa_cost(1.0, p) - (2 * k + 1) / s[1] ** 2) < 1e-12
    assert costs.aa_cost(1.0, p, c_reflect=0.1) > costs.aa_cost(1.0, p)
    assert costs.aa_cost(1.0, p) < 1 / p


def test_measured_point_meets_target_and_matches_gates():
    cfg = StudyConfig()
    ch = Chain(8, 0.0)
    rows, eps_knob = measure_point(ch, [0.1, 0.9], [1e-2, 1e-3], cfg)
    for r, e_k in zip(rows, eps_knob):
        for e, v in r.items():
            for pr in cfg.protocols:
                assert v[pr]["infid"] <= e
            assert v["aa"]["cost"] <= v["rus"]["cost"] * (1 + 1e-12)
            b = v["rus"]
            _, p, inf = trotter.dyadic_filter(ch, ch.trial(e_k), b["K"], b["r"], b["n"])
            assert abs(p - b["p"]) < 1e-10 and abs(inf - b["infid"]) < 1e-10


def test_config_overrides():
    cfg = StudyConfig.from_cli(None, ["Ns=[8,10]", "J2=0.3"])
    assert cfg.Ns == (8, 10) and cfg.J2 == 0.3


def test_aa_cost_no_overflow_at_tiny_p():
    from ornl import costs
    c = costs.cost(603979080, 1.5e-30, "aa", 1000)
    assert 0 < c < 603979080 / 1.5e-30                              # positive, and below RUS (~ pi / (4 sqrt p) rounds)
