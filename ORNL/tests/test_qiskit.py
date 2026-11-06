"""Engine vs real Qiskit circuits (skipped if qiskit is not installed). Run from ORNL/: python -m pytest tests -q"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("qiskit")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qiskit import QuantumCircuit  # noqa: E402
from qiskit.quantum_info import Operator, SparsePauliOp  # noqa: E402
from scipy.linalg import expm  # noqa: E402

from ornl import StudyConfig, chain_for_gap, measure_point  # noqa: E402
from ornl import qiskit_check as qk  # noqa: E402


def test_bond_block_is_exact_and_6_cnots():
    phi = 0.37
    qc = QuantumCircuit(3)
    qk.bond_block(qc, 0, 1, 2, phi)
    target = expm(-1j * phi * SparsePauliOp(["ZXX", "ZYY", "ZZZ"]).to_matrix())
    U = Operator(qc).data
    ph = U[0, 0] / target[0, 0]
    assert np.allclose(U, ph * target, atol=1e-10) and abs(abs(ph) - 1) < 1e-10
    assert qc.count_ops()["cx"] == 6


@pytest.mark.parametrize("compile", ["bond", "pauli"])
def test_circuit_matches_engine(compile):
    cfg = StudyConfig()
    ch = chain_for_gap(6, 0.8, cfg)
    rows, eps_knob = measure_point(ch, [0.1, 0.5], [1e-2], cfg)
    for r, e in zip(rows, eps_knob):
        b = r[1e-2]["rus"]
        v = qk.verify(ch, ch.trial(e), b["K"], b["r"], b["n"], compile=compile, aa=(compile == "bond"))
        assert abs(v["p"] - b["p"]) < 1e-9 and abs(v["infid"] - b["infid"]) < 1e-9
        if compile == "bond":
            assert v["cnot_circuit"] == v["cnot_formula"]
            assert abs(v["aa_success"] - v["aa_success_formula"]) < 1e-9


def test_mps1_trial_circuit_matches_engine():
    cfg = StudyConfig(trial="mps1")
    ch = chain_for_gap(6, 0.8, cfg)
    assert 0.9 < ch.gamma0 < 1
    rows, eps_knob = measure_point(ch, [0.3], [1e-2], cfg)
    b, e = rows[0][1e-2]["rus"], eps_knob[0]
    v = qk.verify(ch, ch.trial(e), b["K"], b["r"], b["n"], aa=False, eps=e)
    assert abs(v["gamma_circuit"] - 0.3) < 1e-9
    assert abs(v["p"] - b["p"]) < 1e-9 and abs(v["infid"] - b["infid"]) < 1e-9
    assert v["prep_cnot_circuit"] == ch.prep_cnot
