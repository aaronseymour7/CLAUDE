"""Simple, transparent resource model for the Trotterised Stetcu-Baroni filter.

Circuit (one ancilla, qubit N): for each pulse i,
    H_anc . [controlled evolution e^{-i t_i H_s (x) Z_anc}, k_i product-formula steps] . Rz_anc(2 phi_i) . H_anc,
then measure the ancilla (post-select 0; abort the attempt on 1).

Counting rules (all-to-all connectivity, Clifford+Rz):
* a term c P (x) Z_anc with |P| = w is a weight-(w+1) Pauli rotation: 2w CNOTs, one Rz and
  single-qubit Clifford basis changes (not counted);
* Lie step: one rotation per non-identity term; Strang step: 2*Gamma - 1 rotations (the two middle
  half-steps merge), and adjacent steps share their outer half-step, so k Strang steps cost
  k (2 Gamma - 2) + 1 rotations;
* the identity term of H_s becomes c_0 t Z_anc, which commutes with every other factor, so it merges
  with Rz(2 phi_i) into one ancilla rotation per pulse;
* arbitrary-angle Rz cost T_rot = a log2(1/delta) + b T gates at synthesis error delta
  (default a = 3, b = 0: the leading term of Ross-Selinger gridsynth);
* expected cost per success = (trial + filter) / p_succ for the repeat-until-success protocol.
  An early-abort protocol (measure after each pulse) is cheaper; it is not modelled.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .pauli import PauliSum, weight


@dataclass
class GateModel:
    t_per_rot_a: float = 3.0
    t_per_rot_b: float = 0.0

    def t_per_rotation(self, delta: float) -> float:
        return self.t_per_rot_a * np.log2(1.0 / delta) + self.t_per_rot_b


def per_step_counts(H: PauliSum):
    h = H.without_identity()
    w = np.array([weight(l) for l in h.labels])
    return dict(terms=len(w), cnot_lie=int(np.sum(2 * w)), rot_lie=len(w),
                cnot_per_term=2 * w)


def filter_counts(H: PauliSum, steps, order: int, phases=None):
    """CNOTs and arbitrary rotations of the filter (no trial circuit)."""
    s = per_step_counts(H)
    G, cw = s["terms"], s["cnot_per_term"]
    k = np.asarray(steps, int)
    if order == 1:
        rot = int(np.sum(k) * G)
        cnot = int(np.sum(k) * cw.sum())
    else:
        # per pulse: k(2G-2)+1 rotations; the first and last term are applied k+1 / k times, the rest 2k
        rot = int(np.sum(k * (2 * G - 2) + 1))
        inner = cw[1:-1].sum() if G > 2 else 0
        cnot = int(np.sum(2 * k * inner + (k + 1) * cw[0] + k * cw[-1])) if G > 1 else int(np.sum(k) * cw.sum())
    rot += len(k)                                    # one merged ancilla Rz per pulse
    return dict(cnot=cnot, rotations=rot, steps=int(np.sum(k)),
                qubits=H.num_qubits + (1 if len(k) else 0))


def total_costs(H, steps, order, trial_cost, p_succ, synth_delta, model: GateModel = GateModel()):
    f = filter_counts(H, steps, order)
    t_rot = model.t_per_rotation(synth_delta) if synth_delta > 0 else 0.0
    tc = dict(cnot=trial_cost.get("cnot", 0), rotations=trial_cost.get("rotations", 0))
    attempt_cnot = f["cnot"] + tc["cnot"]
    attempt_rot = f["rotations"] + tc["rotations"]
    return dict(
        qubits=f["qubits"], steps=f["steps"],
        filter_cnot=f["cnot"], filter_rotations=f["rotations"],
        trial_cnot=tc["cnot"], trial_rotations=tc["rotations"],
        attempt_cnot=attempt_cnot, attempt_rotations=attempt_rot,
        synth_delta=synth_delta, t_per_rotation=t_rot,
        attempt_t=attempt_rot * t_rot,
        p_succ=p_succ,
        expected_cnot=attempt_cnot / p_succ, expected_t=attempt_rot * t_rot / p_succ,
    )
