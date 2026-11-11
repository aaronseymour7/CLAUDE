"""The Trotter simulator applies exactly the product the transpiled Qiskit circuits apply."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from j1j2 import inputs, trotter  # noqa: E402
from qpf import core  # noqa: E402


@pytest.mark.parametrize("order", [1, 2])
def test_circuit_matches_simulator(order):
    from qiskit import transpile
    from qiskit.quantum_info import Statevector
    N = 5
    ex = inputs.exact(N)
    rng = np.random.default_rng(1)
    v = rng.normal(size=1 << N) + 1j * rng.normal(size=1 << N)
    v /= np.linalg.norm(v)
    idx = np.arange(1 << N)
    for op, t, d in ((trotter.scaled_hamiltonian(N, inputs.bonds(N), ex["E0"], ex["Emax"]), 4.3, 0.3),
                     (trotter.total_spin_sq(N), np.pi / 4, 0.0)):
        for n in (1, 3):
            tq = transpile(trotter.step_circuit(op, t, d, n, order), basis_gates=["cx", "rz", "sx", "x"],
                           optimization_level=3, seed_transpiler=0)
            out = Statevector(np.kron([1, 0], v)).evolve(tq).data[: 1 << N]      # ancilla (qubit N) reads 0
            assert np.allclose(out, trotter.pulse(op, idx, v, t, d, n, order), atol=1e-10)


def test_trotter_converges_to_exact_step():
    N = 6
    ex = inputs.exact(N)
    op = trotter.scaled_hamiltonian(N, inputs.bonds(N), ex["E0"], ex["Emax"])
    v = np.random.default_rng(2).normal(size=1 << N)
    v /= np.linalg.norm(v)
    ideal = core.Operator(core.scale(ex["H"], ex["E0"], ex["Emax"])).cos(5.0, 0.1, v)
    err = [np.linalg.norm(trotter.pulse(op, np.arange(1 << N), v, 5.0, 0.1, n, 2) - ideal) for n in (10, 20)]
    assert err[1] < err[0] / 3.5                                                   # second order: ~ 1/4
