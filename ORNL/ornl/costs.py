"""Gate counts and expected cost per success.

CNOTs of one attempt (filter only, trial preparation not counted):
  each ancilla-controlled Heisenberg bond is one block of `bond_cnots` (6: Bell change, 1 CNOT each side, and a
  4-CNOT parity network for exp(-i b (Z_i Z_a + Z_j Z_a - Z_i Z_j Z_a))). A Strang step has 2B - 1 blocks (the two
  middle half steps merge) and consecutive steps share their outer half step, so a pulse of k steps has
  k (2B - 2) + 1 blocks, B = 2N - 3 bonds. `ornl.qiskit_check` builds this circuit and counts it.

Protocols:
  rus  repeat until success: C_attempt / p
  aa   amplitude amplification: k rounds, A used 2k + 1 times, the two reflections k times; success
       sin^2((2k+1) theta), sin theta = sqrt(p). Reflections: R_good flips "all K r ancillas = 0" (one ancilla per
       pulse: nothing may be measured), R_0 flips |0...0> on all N + K r qubits; an n-qubit controlled Z is taken
       as mcz_cnots_per_qubit * n CNOTs. AA is only used where it beats rus.
"""
import numpy as np


def cnot_per_attempt(N, K, r, n, bond_cnots=6):
    k = np.tile(n * 2 ** np.arange(K), r)                 # Strang steps per pulse
    B = 2 * N - 3
    return int(np.sum(bond_cnots * (k * (2 * B - 2) + 1)))


def cnot_mcz(n, per_qubit=12):
    return 0 if n <= 1 else per_qubit * n


def reflection_cnots(N, K, r, per_qubit=12):
    m = K * r
    return cnot_mcz(m, per_qubit) + cnot_mcz(N + m, per_qubit)


def aa_rounds(p):
    th = np.arcsin(np.sqrt(np.clip(np.asarray(p, float), 1e-300, 1)))
    return np.maximum(0, np.round(np.pi / (4 * th) - 0.5)).astype(int), th


def aa_cost(c_attempt, p, c_reflect=0.0):
    k, th = aa_rounds(p)
    k = k.astype(float)                                   # int64 overflows at tiny p (k ~ 1e15) and the cost goes negative
    c_attempt = float(c_attempt)
    aa = (c_attempt * (2 * k + 1) + c_reflect * k) / np.sin((2 * k + 1) * th) ** 2
    return np.minimum(aa, c_attempt / np.asarray(p, float))


def cost(c_attempt, p, protocol, c_reflect=0.0):
    if protocol == "rus":
        return c_attempt / p
    if protocol == "aa":
        return float(aa_cost(c_attempt, p, c_reflect))
    raise ValueError(protocol)
