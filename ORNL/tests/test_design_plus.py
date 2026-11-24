import numpy as np

from qpf import design, design_plus as dp


def test_v4_cost_matches_authors():
    q = design.optimize(0.3, 4, method="v4")
    assert abs(dp.v4_cost(0.3, q["times"]) - q["cost"]) < 1e-9


def test_minimax_beats_v4_worst_case():
    g = 0.3
    v4 = design.optimize(g, 4, method="v4")
    mm = dp.design(g, "minimax", c=1.0, m=4, n_random=4)
    assert abs(mm["times"].sum() - np.pi / g) < 1e-8
    assert mm["times"].max() <= np.pi / g / 3 + 1e-9
    assert mm["max_leak"] <= dp.max_leak(g, v4["times"]) * 1.001
    assert dp.design(g, "minimax", c=1.5, n_random=4)["max_leak"] < 1e-3


def test_phases_respect_success_floor():
    r = dp.design_phases(0.3, c=1.0, p0=0.9, m=4, n_random=2)
    assert r["f0"] ** 2 >= 0.9 * (1 - 1e-6)
    assert r["max_leak"] < dp.design(0.3, "minimax", c=1.0, m=4, n_random=4)["max_leak"]
