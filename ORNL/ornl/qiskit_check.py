"""Real circuits in Qiskit, to check the fast engine at small N (noiseless state-vector simulation).

Register: system qubits 0..N-1 (site i = qubit i), then one ancilla per pulse (K r of them), so that the whole
attempt is one unitary A with no mid-circuit measurement; post-selecting "all ancillas 0" at the end is the same
as measuring after each pulse. This is also the form amplitude amplification needs.

One pulse on ancilla a, time t = m dt:
    H(a) . [m merged Strang steps of exp(-i dt H_s (x) Z_a)] . Rz_a(-2 t E0 / W) . H(a)
H_s = (H - E0)/W. Its Pauli form is sum_b (J_b / 4W)(XX + YY + ZZ)_b - (E0/W) I; the identity part, multiplied
by Z_a, is the ancilla Rz.

Controlled bond exp(-i phi (XX+YY+ZZ) (x) Z_a), two compilations:
  "bond"   Bell-basis change on (i, j), synthesised by Qiskit's KAK TwoQubitBasisDecomposer (1 CNOT each side),
           which turns XX+YY+ZZ into Z_i + Z_j - Z_i Z_j; then a 4-CNOT parity network onto the ancilla.
           6 CNOTs per block, the count `ornl.costs` uses.
  "pauli"  one PauliEvolutionGate per bond (XXZ, YYZ, ZZZ commute, so it is exact), left to Qiskit's
           transpiler at optimization_level=3 (block consolidation + KAK resynthesis of 2-qubit blocks).
           Shows what a generic compile gives.
"""
import warnings

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import CXGate, PauliEvolutionGate, StatePreparation
from qiskit.quantum_info import Operator, SparsePauliOp, Statevector
from qiskit.synthesis import TwoQubitBasisDecomposer

from .costs import aa_rounds, cnot_per_attempt, cnot_mcz
from .trotter import strang_sequence

warnings.filterwarnings("ignore", category=Warning, module="scipy.sparse")

KAK = TwoQubitBasisDecomposer(CXGate())

# Bell change: after it, XX + YY + ZZ is diag(1, 1, 1, -3) = Z_0 + Z_1 - Z_0 Z_1
_bell = QuantumCircuit(2)
_bell.cx(0, 1)
_bell.h(0)
BELL = KAK(Operator(_bell).data)                 # KAK-synthesised (1 CNOT)
BELL_DG = BELL.inverse()
BOND_PAULI = SparsePauliOp(["ZXX", "ZYY", "ZZZ"])  # qubit order (i, j, a) -> labels read a j i


def bond_block(qc, i, j, a, phi, compile="bond"):
    """Append exp(-i phi (X_iX_j + Y_iY_j + Z_iZ_j) Z_a)."""
    if compile == "pauli":
        qc.append(PauliEvolutionGate(BOND_PAULI, time=phi), [i, j, a])
        return
    qc.compose(BELL, [i, j], inplace=True)
    qc.cx(i, a); qc.rz(2 * phi, a)               # exp(-i phi Z_i Z_a)
    qc.cx(j, a); qc.rz(-2 * phi, a)              # exp(+i phi Z_i Z_j Z_a)
    qc.cx(i, a); qc.rz(2 * phi, a)               # exp(-i phi Z_j Z_a)
    qc.cx(j, a)
    qc.compose(BELL_DG, [i, j], inplace=True)


def append_pulse(qc, ch, a, dt, m, compile="bond"):
    """One filter pulse of m Strang steps of length dt, on ancilla a."""
    qc.h(a)
    for b, w in strang_sequence(len(ch.bonds), m):
        i, j, J = ch.bonds[b]
        bond_block(qc, i, j, a, w * (dt / 2) * J / (4 * ch.W), compile)
    qc.rz(-2 * m * dt * ch.E0 / ch.W, a)
    qc.h(a)


def filter_circuit(ch, K, r, n, compile="bond"):
    """The filter part of one attempt: K r pulses, each on its own ancilla (qubits N .. N + K r - 1)."""
    qc = QuantumCircuit(ch.N + K * r)
    tau = np.pi / (1 + ch.gap_s)
    a = ch.N
    for _ in range(r):
        for k in range(K):
            append_pulse(qc, ch, a, tau / n, n * 2 ** k, compile)
            a += 1
    return qc


def attempt_circuit(ch, psi, K, r, n, compile="bond"):
    """A: prepare the trial (exact StatePreparation, not costed) and run the filter."""
    qc = QuantumCircuit(ch.N + K * r)
    qc.append(StatePreparation(ch.to_full(psi)), range(ch.N))
    qc.compose(filter_circuit(ch, K, r, n, compile), inplace=True)
    return qc


def flip_zero(qc, qubits):
    """Sign flip of |0...0> on `qubits`: X layer, multi-controlled Z, X layer."""
    qubits = list(qubits)
    qc.x(qubits)
    if len(qubits) == 1:
        qc.z(qubits[0])
    else:
        qc.h(qubits[-1]); qc.mcx(qubits[:-1], qubits[-1]); qc.h(qubits[-1])
    qc.x(qubits)


def aa_circuit(A, n_sys, k):
    """A followed by k Grover rounds Q = -A R_0 A^dagger R_good (R_good: ancillas all 0; R_0: everything 0)."""
    nq = A.num_qubits
    qc = A.copy()
    Adg = A.inverse()
    for _ in range(k):
        flip_zero(qc, range(n_sys, nq))
        qc.compose(Adg, inplace=True)
        flip_zero(qc, range(nq))
        qc.compose(A, inplace=True)
    return qc


def postselect(state, ch, n_anc):
    """System state with all ancillas 0 (unnormalised), p, and 1 - F against the exact ground state."""
    v = np.asarray(state).reshape(2 ** n_anc, 2 ** ch.N)[0]      # ancillas are the high qubits
    p = float(np.vdot(v, v).real)
    g = ch.to_full(ch.g0)
    return v, p, 1 - abs(np.vdot(g, v)) ** 2 / p


def cnot_count(qc, transpile_first=True, seed=1):
    if transpile_first:
        qc = transpile(qc, basis_gates=["cx", "rz", "sx", "x"], optimization_level=3, seed_transpiler=seed)
    return int(qc.count_ops().get("cx", 0))


def verify(ch, psi, K, r, n, compile="bond", aa=True, count_transpiled=False):
    """Simulate the real circuits and return a dict to compare with the engine:
    p, 1-F, CNOTs (built vs formula); with aa: rounds k, success after k rounds vs sin^2((2k+1) theta),
    and the CNOTs of the two reflections vs the cost model."""
    m = K * r
    A = attempt_circuit(ch, psi, K, r, n, compile)
    _, p, infid = postselect(Statevector(A).data, ch, m)
    filt = filter_circuit(ch, K, r, n, compile)
    out = dict(N=ch.N, K=K, r=r, n=n, qubits=A.num_qubits, compile=compile, p=p, infid=infid,
               cnot_formula=cnot_per_attempt(ch.N, K, r, n),
               cnot_circuit=cnot_count(filt, transpile_first=(compile == "pauli" or count_transpiled)))
    if aa:
        k, th = aa_rounds(p)
        k = int(k)
        _, p_aa, inf_aa = postselect(Statevector(aa_circuit(A, ch.N, k)).data, ch, m)
        rg, r0 = QuantumCircuit(A.num_qubits), QuantumCircuit(A.num_qubits)
        flip_zero(rg, range(ch.N, A.num_qubits)); flip_zero(r0, range(A.num_qubits))
        out.update(aa_k=k, aa_success=p_aa, aa_success_formula=float(np.sin((2 * k + 1) * th) ** 2),
                   aa_infid=inf_aa, reflect_cnot_circuit=cnot_count(rg) + cnot_count(r0),
                   reflect_cnot_model=cnot_mcz(m) + cnot_mcz(ch.N + m))
    return out
