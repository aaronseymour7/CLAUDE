"""
resources.py -- CX counts for the Trotterised filter (needs qiskit; the rest of the study does not).

One Trotter step of one pulse = prod over bonds of exp(-i dt * c_g * S_i.S_j (x) Z_anc), written as three
Pauli rotations (XX, YY, ZZ on the bond) tensored with Z on the ancilla, transpiled with qiskit (level 3)
to {cx, rz, h, s}.  Total CX = n * (CX per step) -- ancilla Hadamards/Rz and the phase rotations add no CX.
Connectivity: all-to-all, and a line 0-1-...-(N-1)-anc (pessimistic: the ancilla touches every system qubit).
Reference: generic (exact) state preparation of an N-qubit state with qiskit's StatePreparation.
"""
import sys, json, glob, pathlib
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate, StatePreparation
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import workflow as wf

BASIS = ["cx", "rz", "h", "s", "sdg"]


def step_circuit(N, J2, W, steps=1, dt=0.1):
    qc = QuantumCircuit(N + 1)
    for _ in range(steps):
        for (i, j, J) in wf.bonds(N, 1.0, J2):
            c = J / W
            for P in "XYZ":
                lab = ["I"] * (N + 1)
                lab[i] = lab[j] = P
                lab[N] = "Z"
                op = SparsePauliOp("".join(lab[::-1]), c / 4)
                qc.append(PauliEvolutionGate(op, time=dt), range(N + 1))
    return qc


def cx_per_step(N, J2, W, cmap=None):
    out = {}
    for k in (1, 3):
        rc = transpile(step_circuit(N, J2, W, k), basis_gates=BASIS, coupling_map=cmap,
                       optimization_level=3, seed_transpiler=0)
        out[k] = rc.count_ops().get("cx", 0)
    return out[1], out[3] / 3          # single step, and amortised over 3 steps (cross-step cancellation)


def generic_prep_cx(N):
    rng = np.random.default_rng(0)
    v = rng.standard_normal(2 ** N) + 1j * rng.standard_normal(2 ** N)
    qc = QuantumCircuit(N)
    qc.append(StatePreparation(v / np.linalg.norm(v)), range(N))
    rc = transpile(qc, basis_gates=BASIS, optimization_level=1, seed_transpiler=0)
    return rc.count_ops().get("cx", 0)


if __name__ == "__main__":
    cases = [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "results" / "case_*chi2_eps0.01.json")))]
    cases = [c for c in cases if "at_n_bound" in c]
    rows = []
    gp = {}
    for c in cases:
        N, J2 = c["N"], c["J2"]
        ch = wf.Chain(N, J2)
        a2a1, a2a3 = cx_per_step(N, J2, ch.W)
        line1, line3 = cx_per_step(N, J2, ch.W, CouplingMap.from_line(N + 1))
        if N not in gp and N <= 12:
            gp[N] = generic_prep_cx(N)
        a = c["at_n_bound"]
        rows.append(dict(N=N, J2=J2, n_bound=c["n_bound"], n_emp=c["n_emp"], psucc=a["psucc"],
                         cx_step_a2a=a2a1, cx_step_a2a_amortised=a2a3, cx_step_line=line1, cx_step_line_amortised=line3,
                         cx_bound=c["n_bound"] * a2a1, cx_emp=c["n_emp"] * a2a1,
                         cx_bound_line=c["n_bound"] * line1, cx_emp_line=c["n_emp"] * line1,
                         generic_prep_cx=gp.get(N)))
        print(rows[-1], flush=True)
    (HERE / "results" / "resources.json").write_text(json.dumps(rows, indent=1, default=float))
