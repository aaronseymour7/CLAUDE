"""Trotterised QPF steps for operators that are sums of SWAPs, as Qiskit circuits and as a fast exact simulator.

The QPF step exp(-i (t O + delta) (x) Y_a) is built as V . exp(-i (t O + delta) (x) Z_a) . V^dagger with
V = S H on the ancilla (S H Z H S^dagger = Y). With O = c + sum_k w_k SWAP_k and SWAP = (I + XX + YY + ZZ)/2,

    O (x) Z_a = (c + sum_k w_k / 2) Z_a + sum_k (w_k / 2) (X X Z + Y Y Z + Z Z Z)_(i_k, j_k, a)

and exp(-i t O (x) Z_a) is Trotterised with Qiskit's PauliEvolutionGate: LieTrotter (first order) or
SuzukiTrotter(order=2) with `reps` steps; the delta part is one ancilla Rz. Within one SWAP the three Pauli terms
commute, so the Trotter error comes only from the order of the SWAP terms (H: J1 bonds, then J2 bonds; J^2: all
pairs i < j), and the simulator below applies exactly the same product: each controlled SWAP factor
exp(-i s w dt SWAP) = cos(w dt) - i s sin(w dt) SWAP on the branch with ancilla Z = s = +-1. Reading the ancilla 0
keeps (e^{-i delta} U(t) + e^{i delta} U(-t)) / 2, with U(+-t) the Trotter product of the two branches.

Gate counts: the circuits are transpiled with optimization_level=3 to {cx, rz, sx, x}; level 3 consolidates
two-qubit blocks and resynthesises them with Qiskit's two-qubit (KAK) basis decomposer. An Rz is non-Clifford when
its angle is not a multiple of pi/2. No coupling map (all-to-all connectivity). Qiskit's SuzukiTrotter does not merge
the half-steps of neighbouring repetitions, so a second-order step costs about twice a first-order one.
"""
import numpy as np

CLIFFORD_TOL = 1e-9


# ---------------------------------------------------------------- operators as SWAP sums
class SwapSum:
    """O = const + sum_k w_k SWAP_(i_k, j_k) on N qubits (Qiskit ordering, bit i = site i)."""

    def __init__(self, N, pairs, weights, const):
        self.N, self.pairs, self.weights, self.const = N, list(pairs), np.asarray(weights, float), float(const)
        self._perm = {}

    def perms(self, index):
        """SWAP permutations on the basis `index` (array of basis integers, sorted)."""
        key = (len(index), int(index[0]), int(index[-1]))
        if key not in self._perm:
            pos = {int(s): k for k, s in enumerate(index)} if len(index) < (1 << self.N) else None
            out = []
            for i, j in self.pairs:
                flip = ((index >> i) & 1) != ((index >> j) & 1)
                t = np.where(flip, index ^ (1 << i) ^ (1 << j), index)
                out.append(t if pos is None else np.array([pos[int(x)] for x in t]))
            self._perm[key] = out
        return self._perm[key]


def scaled_hamiltonian(N, bonds, E0, Emax):
    """H_s = (H - E0) / (Emax - E0) with H = sum_b J_b (2 SWAP_b - 1) / 4."""
    W = Emax - E0
    c_id = -0.25 * sum(J for *_, J in bonds)
    return SwapSum(N, [(i, j) for i, j, _ in bonds], [J / (2 * W) for *_, J in bonds], (c_id - E0) / W)


def total_spin_sq(N):
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    return SwapSum(N, pairs, np.ones(len(pairs)), 0.75 * N - 0.25 * N * (N - 1))


# ---------------------------------------------------------------- simulator
def _branch(op, perms, V, t, n, order, s):
    """Trotter product for exp(-i s t O) on V (s = +-1)."""
    dt = t / n
    if order == 1:
        seq = [(p, w * dt) for p, w in zip(perms, op.weights)]
    elif order == 2:
        half = [(p, w * dt / 2) for p, w in zip(perms, op.weights)]
        seq = half + half[::-1]
    else:
        raise ValueError("order must be 1 or 2")
    for _ in range(n):
        for p, a in seq:
            V = np.cos(a) * V - 1j * s * np.sin(a) * V[p]
    return np.exp(-1j * s * op.const * t) * V


def pulse(op, index, V, t, delta, n, order):
    """Unnormalised post-selected output of one Trotterised QPF step on V (vectors on the basis `index`)."""
    perms = op.perms(index)
    return 0.5 * (np.exp(-1j * delta) * _branch(op, perms, V, t, n, order, 1)
                  + np.exp(1j * delta) * _branch(op, perms, V, t, n, order, -1))


def steps_for(t, dt):
    return max(1, int(np.ceil(abs(t) / dt - 1e-12)))


def run(op, index, V, times, phases, dt, order):
    """All steps, n_i = ceil(|t_i| / dt) Trotter steps each. Returns (unnormalised state, list of n_i)."""
    ns = []
    for t, d in zip(times, phases):
        n = steps_for(t, dt)
        ns.append(n)
        V = pulse(op, index, V, t, d, n, order)
    return V, ns


# ---------------------------------------------------------------- Qiskit circuits and counts
def step_circuit(op, t, delta, n, order):
    """The QPF step exp(-i (t O + delta) (x) Y_a) with a Trotterised evolution; ancilla = qubit N."""
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import PauliEvolutionGate
    from qiskit.quantum_info import SparsePauliOp
    from qiskit.synthesis import LieTrotter, SuzukiTrotter
    N = op.N
    terms = [("Z", [N], op.const + 0.5 * op.weights.sum())]
    for (i, j), w in zip(op.pairs, op.weights):
        terms += [(p * 2 + "Z", [i, j, N], w / 2) for p in "XYZ"]
    A = SparsePauliOp.from_sparse_list(terms, num_qubits=N + 1)
    synth = LieTrotter(reps=n) if order == 1 else SuzukiTrotter(order=2, reps=n)
    qc = QuantumCircuit(N + 1)
    qc.sdg(N); qc.h(N)
    qc.append(PauliEvolutionGate(A, time=t, synthesis=synth), range(N + 1))
    qc.rz(2 * delta, N)
    qc.h(N); qc.s(N)
    return qc


def counts(qc):
    """(CNOT, non-Clifford Rz) after transpiling to {cx, rz, sx, x} at optimization level 3."""
    from qiskit import transpile
    tq = transpile(qc, basis_gates=["cx", "rz", "sx", "x"], optimization_level=3, seed_transpiler=0)
    rz = [float(g.operation.params[0]) for g in tq.data if g.operation.name == "rz"]
    nc = sum(abs(a / (np.pi / 2) - round(a / (np.pi / 2))) > CLIFFORD_TOL for a in rz)
    return int(tq.count_ops().get("cx", 0)), int(nc)


_COUNTS = {}
N_DIRECT = 32


def step_counts(op, key, t, delta, n, order):
    """(CNOT, non-Clifford Rz) of one transpiled QPF step. Circuits with n <= N_DIRECT Trotter steps are transpiled
    as they are; above that the counts are extrapolated linearly from the transpiled n = N_DIRECT / 2 and N_DIRECT
    circuits with the same step dt = t / n (both counts grow linearly in n to within a few per cent; the
    extrapolation is flagged in the results)."""
    k = (key, round(t, 12), round(delta, 12), n, order)
    if k not in _COUNTS:
        if n <= N_DIRECT:
            _COUNTS[k] = counts(step_circuit(op, t, delta, n, order))
        else:
            a = step_counts(op, key, t * (N_DIRECT // 2) / n, delta, N_DIRECT // 2, order)
            b = step_counts(op, key, t * N_DIRECT / n, delta, N_DIRECT, order)
            s = (n - N_DIRECT) / (N_DIRECT // 2)
            _COUNTS[k] = tuple(int(round(b[i] + s * (b[i] - a[i]))) for i in (0, 1))
    return _COUNTS[k]
