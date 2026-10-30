"""Simple, transparent resource model for the Trotterised Stetcu-Baroni filter.

Circuit (one ancilla, qubit N): for each pulse i,
    H_anc . [controlled evolution e^{-i t_i H_s (x) Z_anc}, k_i product-formula steps] . Rz_anc(2 phi_i) . H_anc,
then measure the ancilla (post-select 0; abort the attempt on 1).

Counting rules (all-to-all connectivity, Clifford+Rz):
* compile = "bond" (default): the commuting terms of one Heisenberg bond, c (XX + YY + ZZ) (x) Z_anc,
  form one exact three-qubit block. A Bell-basis change on the pair (one CNOT each side) makes it
  diagonal, exp(-i c t (Z_i Z_a + Z_j Z_a - Z_i Z_j Z_a)), and a four-CNOT parity network on the
  ancilla applies the three phases: 6 CNOTs and 3 Rz per block instead of 12 CNOTs (the uncontrolled
  bond is the 3-CNOT KAK gate, so this is the controlled-gate analogue). Any other group falls back to
  the ladder rule below. The block is the same unitary as its three term rotations, so the Trotter
  circuit, its error and every simulated infidelity are unchanged;
* compile = "ladder": a term c P (x) Z_anc with |P| = w is a weight-(w+1) Pauli rotation: 2w CNOTs,
  one Rz and single-qubit Clifford basis changes (not counted);
* Lie step: one rotation per non-identity term; Strang step: 2*Gamma - 1 rotations (the two middle
  half-steps merge), and adjacent steps share their outer half-step, so k Strang steps cost
  k (2 Gamma - 2) + 1 rotations;
* the identity term of H_s becomes c_0 t Z_anc, which commutes with every other factor, so it merges
  with Rz(2 phi_i) into one ancilla rotation per pulse;
* arbitrary-angle Rz cost T_rot = a log2(1/delta) + b T gates at synthesis error delta
  (default a = 3, b = 0: the leading term of Ross-Selinger gridsynth);
* expected cost per success = (trial + filter) / p_succ for the repeat-until-success protocol, with
  p_succ the certified lower bound (sqrt(p_g) - eps_U)^2 of Theorem 1 (pipeline._plan_for).
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


COMPILE = "bond"
BOND_BLOCK_CNOT = 6


def _is_heisenberg_bond(labels):
    sup = {tuple(i for i, c in enumerate(l) if c != "I") for l in labels}
    if len(sup) != 1 or len(next(iter(sup))) != 2 or len(labels) != 3:
        return False
    i, j = next(iter(sup))
    return sorted(l[i] + l[j] for l in labels) == ["XX", "YY", "ZZ"]


def per_block_counts(H: PauliSum, compile: str = None):
    """Trotter factors as compiled: per block (consecutive terms sharing a group id) its CNOTs and Rz."""
    compile = compile or COMPILE
    h = H.without_identity()
    if compile == "ladder":
        w = np.array([weight(l) for l in h.labels])
        return 2 * w, np.ones(len(w), int)
    blocks, cur = [], []
    for l, g in zip(h.labels, h.groups):
        if cur and g != cur[0][1]:
            blocks.append([x for x, _ in cur]); cur = []
        cur.append((l, g))
    if cur:
        blocks.append([x for x, _ in cur])
    cw, rz = [], []
    for b in blocks:
        if _is_heisenberg_bond(b):
            cw.append(BOND_BLOCK_CNOT); rz.append(3)
        else:
            cw.append(sum(2 * weight(l) for l in b)); rz.append(len(b))
    return np.array(cw), np.array(rz)


def filter_counts(H: PauliSum, steps, order: int, phases=None, anc_rot_per_pulse: int = 1,
                  compile: str = None):
    """CNOTs and arbitrary rotations of the filter (no trial circuit). A Stetcu-Baroni pulse has one
    merged ancilla Rz; a QETU query has two (the identity term's Z rotation and the X-rotation phase).
    Strang merging acts on the compiled blocks (bonds by default, see the module docstring)."""
    cw, rw = per_block_counts(H, compile)
    G = len(cw)
    k = np.asarray(steps, int)
    if order == 1:
        rot = int(np.sum(k) * rw.sum())
        cnot = int(np.sum(k) * cw.sum())
    else:
        # per pulse: the first and last block are applied k+1 / k times, the rest 2k
        if G > 1:
            inner, rinner = (cw[1:-1].sum(), rw[1:-1].sum()) if G > 2 else (0, 0)
            cnot = int(np.sum(2 * k * inner + (k + 1) * cw[0] + k * cw[-1]))
            rot = int(np.sum(2 * k * rinner + (k + 1) * rw[0] + k * rw[-1]))
        else:
            cnot, rot = int(np.sum(k) * cw.sum()), int(np.sum(k) * rw.sum())
    rot += anc_rot_per_pulse * len(k)                # ancilla rotations per pulse / query
    return dict(cnot=cnot, rotations=rot, steps=int(np.sum(k)),
                qubits=H.num_qubits + (1 if len(k) else 0))


def total_costs(H, steps, order, trial_cost, p_succ, synth_delta, model: GateModel = GateModel(),
                anc_rot_per_pulse: int = 1, compile: str = None):
    f = filter_counts(H, steps, order, anc_rot_per_pulse=anc_rot_per_pulse, compile=compile)
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
