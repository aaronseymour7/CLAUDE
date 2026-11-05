"""MPS mode vs the exact engine (skipped if quimb is not installed). Run from ORNL/: python -m pytest tests -q"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("quimb")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ornl import Chain, trotter  # noqa: E402
from ornl.mps import MPSChain  # noqa: E402


@pytest.mark.parametrize("reference", ["ed", "dmrg"])
def test_mps_matches_exact(reference):
    ch = Chain(8, 0.15)
    m = MPSChain(8, 0.15, reference=reference)
    assert abs(m.E0 - ch.E0) < 1e-9 and abs(m.E1 - ch.E1) < 1e-8 and abs(m.Emax - ch.Emax) < 1e-8
    e = ch.eps_for_overlap(0.2)
    assert abs(m.overlap(e) - 0.2) < 1e-8
    for K, r, n in [(3, 1, 1), (2, 2, 2)]:
        _, p, inf = trotter.dyadic_filter(ch, ch.trial(e), K, r, n)
        _, pm, infm = m.filter(m.trial(e), K, r, n)
        assert abs(pm - p) < 1e-7 and abs(infm - inf) < 1e-7


def test_mps_pulse_is_the_post_selected_branch():
    ch = Chain(6, -0.2)
    m = MPSChain(6, -0.2, reference="ed")
    e = ch.eps_for_overlap(0.4)
    w = m.pulse(m.trial(e), 0.7, 3).to_dense().ravel()[ch.states]
    we = trotter.pulse(ch, ch.trial(e), 0.7, 3)
    assert abs(abs(np.vdot(w, we)) - np.vdot(we, we).real) < 1e-9
