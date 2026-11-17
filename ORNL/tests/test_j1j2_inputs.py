"""J1-J2 pipeline inputs: operators, DMRG estimates vs ED, trial ordering. Needs quimb and mps-to-circuit."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
pytest.importorskip("quimb")
from j1j2 import inputs  # noqa: E402


def test_operators_and_dmrg_match_ed():
    N = 8
    ex = inputs.exact(N)
    J2 = inputs.total_spin_sq(N)
    assert abs(ex["g0"].conj() @ (J2 @ ex["g0"])) < 1e-10                  # singlet ground state
    est = inputs.dmrg_estimates(N)
    for k in ("E0", "Emax", "E1_singlet"):
        assert abs(est[k] - ex[k]) < 1e-6, k
    v = inputs.mps_statevector(est["psi0"])
    assert abs(np.vdot(ex["g0"], v)) ** 2 > 1 - 1e-8                       # Qiskit ordering is consistent
    pytest.importorskip("mps_to_circuit")
    _, t = inputs.trial(est["psi0"])
    assert abs(np.linalg.norm(t) - 1) < 1e-10 and abs(np.vdot(v, t)) ** 2 > 0.5
