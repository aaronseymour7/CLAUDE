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
    rows, eps_knob, _ = measure_point(ch, [0.1, 0.9], [1e-2, 1e-3], cfg)
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
    for p in (1.5e-30, 1e-130, 1e-285):
        c = costs.cost(603979080, p, "aa", 1000)
        assert 0 < c < 603979080 / p                                # positive, and below RUS (~ pi / (4 sqrt p) rounds)


# ---- total spin: S^2 stage and the singlet knob
from ornl import symmetry  # noqa: E402


def test_s2_ladder_and_projector():
    ch = Chain(8, 0.1)
    S2 = symmetry.s2_matrix(ch).toarray()
    w = np.linalg.eigvalsh(S2)
    J = (-1 + np.sqrt(1 + 4 * w)) / 2
    assert np.allclose(J, np.round(J), atol=1e-9) and abs(np.vdot(ch.g0, S2 @ ch.g0)) < 1e-10   # singlet ground state
    assert symmetry.ladder(16) == (np.pi / 4, np.pi / 8, np.pi / 16)
    for N in (8, 10, 12, 14, 16):                                      # the ladder zeroes every J > 0
        t = np.array(symmetry.ladder(N))
        for Jv in range(1, N // 2 + 1):
            assert abs(np.prod(np.cos(t * Jv * (Jv + 1)))) < 1e-12
    v = ch.trial(0.4)
    P = symmetry.project_singlet(ch, v)
    assert np.linalg.norm(S2 @ P) < 1e-9 and abs(np.vdot(ch.g0, P) - np.vdot(ch.g0, v)) < 1e-12


def test_s2_stage_trotter_converges_and_counts():
    ch = Chain(8, 0.0)
    v = ch.trial(0.5)
    ideal = symmetry.project_singlet(ch, v)
    inf = [symmetry._infid(symmetry.stage(ch, v, n), ideal) for n in (4, 8, 16)]
    assert inf[0] > 8 * inf[1] > 64 * inf[2] and inf[2] < 1e-3          # infidelity ~ dt^4 (second order)
    # gate count: (2P - 2) m + 1 blocks per pulse, P = 28 pairs at N = 8; J <= 4 needs two pulses (8 and 4 steps)
    assert len(symmetry.ladder(8)) == 2
    assert symmetry.stage_cnots(8, 8) == 6 * ((54 * 8 + 1) + (54 * 4 + 1))


def test_singlet_knob_keeps_singlet_and_hits_target():
    ch = Chain(10, 0.0, knob="singlet_bond")
    e = ch.eps_for_overlap(0.01)
    v = ch.trial(e)
    assert abs(ch.overlap(e) - 0.01) < 1e-10
    assert abs(symmetry.singlet_weight(ch, v) - 1) < 1e-10           # all weight stays in J = 0
    assert ch.knob_cnot == 27


def test_singlet_gap_grid():
    cfg = StudyConfig(knob="singlet_bond", gap_kind="singlet")
    ch = chain_for_gap(10, 0.6, cfg)
    assert abs(ch.singlet_gap() - 0.6) < 1e-8 and ch.gap_kind == "singlet" and abs(ch.gap - 0.6) < 1e-8
    assert ch.delta > chain_for_gap.__globals__["singlet_dip"](10)[0]
    assert chain_for_gap(10, 0.05, cfg) is None                       # below the singlet dip


def test_measured_point_with_s2_stage_matches_gates():
    cfg = StudyConfig(s2_stage=True)
    ch = Chain(8, 0.0)
    ch.use_gap("singlet")
    rows, eps_knob, s2 = measure_point(ch, [0.1], [1e-2], cfg)
    b = rows[0][1e-2]["rus"]
    v = symmetry.stage(ch, ch.trial(eps_knob[0]), b["s2_n"])
    assert abs(np.sum(abs(v) ** 2) - b["s2_p"]) < 1e-12
    _, p, inf = trotter.dyadic_filter(ch, v, b["K"], b["r"], b["n"])
    assert abs(p - b["p"]) < 1e-10 and abs(inf - b["infid"]) < 1e-10 and inf <= 1e-2
    assert b["cnot"] == costs.cnot_per_attempt(8, b["K"], b["r"], b["n"]) + symmetry.stage_cnots(8, b["s2_n"])
    assert s2["stage_infid"] == sorted(s2["stage_infid"], reverse=True)
