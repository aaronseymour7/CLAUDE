"""Leading-order T-count model of a qubitization / QSVT ground-state filter, for comparison with the
Trotterised Stetcu-Baroni filter. It is a model with stated assumptions, not a compiled circuit.

Assumptions (every one favours the QSVT side):
* H = sum_a c_a P_a is block-encoded by a Pauli LCU with 1-norm lambda = sum |c_a| and Gamma terms;
  the walk operator (PREP, SELECT, PREP^dagger, reflection) is the query.
* The filter is a polynomial of degree d = c_deg (lambda / g) ln(2 / eta) with c_deg = 1, where g is the
  gap and eta the suppression ratio that gives the same leakage as the SBC plan. Known constructions
  (Lin-Tong, Dong-Lin-Tong) scale like this with constants of order one or larger; parity and
  normalisation overheads are ignored.
* One attempt succeeds with probability gamma (the trial overlap); amplitude amplification is not
  used for either method.
* Costs per query: SELECT by unary iteration 4 (Gamma - 1) T gates; PREP and PREP^dagger together
  2 (Gamma - 1) arbitrary rotations; the reflection 4 ceil(log2 Gamma) T gates; each of the d QSVT
  phase rotations is arbitrary. An arbitrary rotation at accuracy delta costs 3 log2(1/delta) T
  gates, as in resources.py, and delta spends 10% of the same error budget.
* Only T gates are counted (Clifford CNOTs are not), and the qubit count is N + ceil(log2 Gamma) + 3.
"""
from __future__ import annotations

import math
from typing import Dict

from .pauli import PauliSum


def qsvt_filter_cost(H: PauliSum, gap: float, gamma: float, eta: float, eps: float,
                     trial_cost: Dict[str, float], c_deg: float = 1.0) -> Dict[str, float]:
    h = H.without_identity()
    Gamma = len(h)
    lam = float(sum(abs(h.coeffs)))
    d = max(1, math.ceil(c_deg * lam / gap * math.log(2.0 / eta)))
    rot_q = 2 * (Gamma - 1) + 1                       # PREP, PREP^dagger and one QSVT phase per query
    R = d * rot_q + trial_cost.get("rotations", 0)
    budget = 0.1 * math.sqrt(eps) * math.sqrt(gamma)  # share of the error budget for synthesis
    delta = budget / R
    t_rot = 3.0 * math.log2(1.0 / delta)
    t_query = 4 * (Gamma - 1) + 4 * math.ceil(math.log2(max(Gamma, 2)))
    t_attempt = d * (t_query + rot_q * t_rot) + trial_cost.get("rotations", 0) * t_rot
    return dict(lam=lam, Gamma=Gamma, degree=d, rotations=R, delta=delta, t_per_rotation=t_rot,
                t_attempt=t_attempt, t_success=t_attempt / gamma,
                qubits=H.num_qubits + math.ceil(math.log2(max(Gamma, 2))) + 3)
