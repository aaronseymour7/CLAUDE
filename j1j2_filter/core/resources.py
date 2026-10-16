"""Gate-count / T-count resource estimates."""


import numpy as np
from qiskit import QuantumCircuit, transpile

from core.circuits import apply_filter_pulse, flatten_success_path


BASIS = ("h", "s", "cx", "rz")
RZZ_BASIS = ("rzz", "rz", "sx", "x")     # native two-qubit gate = RZZ(theta) = exp(-i theta/2 Z(x)Z); CX is not in the basis


def count_nonclifford_rz(circ, atol=1e-9):
    """Number of rz gates whose angle is not a multiple of pi/2."""
    cnt = 0
    for inst in circ.data:
        if inst.operation.name == "rz":
            x = float(inst.operation.params[0]) / (np.pi / 2)
            if abs(x - round(x)) > atol:
                cnt += 1
    return cnt


def count_nonclifford_rzz(circ, atol=1e-9):
    """Number of rzz gates whose angle is not a multiple of pi/2 (RZZ(pi/2) is Clifford)."""
    cnt = 0
    for inst in circ.data:
        if inst.operation.name == "rzz":
            x = float(inst.operation.params[0]) / (np.pi / 2)
            if abs(x - round(x)) > atol:
                cnt += 1
    return cnt


def rotation_synthesis_t_count(n_rot, eps=1e-3):
    """APPROXIMATE T-count to synthesize n_rot arbitrary Rz rotations to
    per-rotation precision eps, using ~ 1.15 log2(1/eps) + 9.2 T per rotation
    (my recollection of the Kliuchnikov/Ross-Selinger-type scaling -- verify
    the constants before quoting). Only for fault-tolerant cost estimates."""
    return float(n_rot * (1.15 * np.log2(1.0 / eps) + 9.2))


def resource_costs(circuit, label="", flatten=False, verbose=True, basis=BASIS,
                   coupling_map=None, synth_eps=1e-3, seed_transpiler=0, optimization_level=3):
    """Transpile to `basis` (opt level 3) and count gates. flatten=True for the
    early-abort circuit. coupling_map: optional qiskit CouplingMap (routing
    overhead shows up as extra cx). Also returns non-Clifford rz count and an
    approximate T-count for rotation synthesis. Pass the trial circuit
    separately (or inside the circuit via trial_prep) to include it (A11)."""
    circ = flatten_success_path(circuit) if flatten else circuit
    rc = transpile(circ, basis_gates=list(basis) + ["measure", "reset"],
                   coupling_map=coupling_map, optimization_level=optimization_level,
                   seed_transpiler=seed_transpiler)
    counts = rc.count_ops()
    tidy = {g: int(counts.get(g, 0)) for g in basis}
    mid = {g: int(counts[g]) for g in ("measure", "reset") if g in counts}
    other = {g: int(c) for g, c in counts.items()
             if g not in basis and g not in ("measure", "reset", "barrier")}
    n_nc = count_nonclifford_rz(rc)
    n_nzz = count_nonclifford_rzz(rc) if "rzz" in basis else 0
    out = {**tidy, "total": sum(tidy.values()), "depth": rc.depth(),
           "measure": mid.get("measure", 0), "reset": mid.get("reset", 0),
           "rz_nonclifford": n_nc, "rzz_nonclifford": n_nzz, "nc_total": n_nc + n_nzz,
           "t_count_est": rotation_synthesis_t_count(n_nc + n_nzz, synth_eps),
           "other": other}
    if verbose:
        print(f"\n--- resource costs: {label} ---")
        print("  " + "  ".join(f"{g.upper()}={tidy[g]}" for g in basis)
              + f"  | total={out['total']}  depth={out['depth']}"
              + (f"  | mid-circuit: {mid}" if mid else ""))
        print(f"  non-Clifford rz={n_nc}  approx T-count (eps={synth_eps:g}/rot)"
              f"={out['t_count_est']:.0f}")
        if other:
            print(f"  [warning] non-basis ops left after transpile: {other}")
    return out


def cx_depth_per_step(H, order=1, basis=BASIS, optimization_level=3):
    """(CX count, depth) of ONE Trotter step of one pulse, transpiled once.
    Ignores cancellation across step boundaries -- validate against a full
    resource_costs() transpile at more than one n."""
    N = H.num_qubits
    qc = QuantumCircuit(N + 1)
    apply_filter_pulse(qc, list(range(N)), N, H, 1.0, 0.0,
                       trotter_steps=1, order=order)
    rc = transpile(qc, basis_gates=list(basis),
                   optimization_level=optimization_level, seed_transpiler=0)
    return rc.count_ops().get("cx", 0), rc.depth()


def direct_step_counts(H, N, order=1, basis=BASIS):
    """Counts of ONE Trotter step of one pulse compiled DIRECTLY (qiskit optimization level 1: every Pauli rotation exp(-i theta P (x) Z)
    becomes a CX ladder around one Rz, no block re-synthesis), all-to-all.  Level 3 (the default elsewhere) re-synthesises two-qubit
    blocks: ~1.7x fewer CX but ~3x more non-Clifford Rz.  Direct compilation gives, for the Heisenberg chain with ancilla,
    CX = 4 x (#Pauli terms) and non-Clifford Rz = #Pauli terms + 1 per step (verified against a full transpile)."""
    qc = QuantumCircuit(N + 1)
    apply_filter_pulse(qc, list(range(N)), N, H, 1.0, 0.0, trotter_steps=1, order=order)
    c = resource_costs(qc, verbose=False, optimization_level=1, basis=basis)
    keys = ("rzz", "rzz_nonclifford", "rz", "rz_nonclifford", "nc_total", "sx", "x", "depth") if "rzz" in basis \
        else ("cx", "h", "s", "rz", "rz_nonclifford", "depth")
    return {k: float(c[k]) for k in keys}
