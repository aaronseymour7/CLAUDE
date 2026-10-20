"""add_direct_counts.py -- add the DIRECT-compilation (qiskit level 1) gate counts to results/data/*.json without rerunning DMRG/design.

Per-step counts come from core.resources.direct_step_counts on the scaled Pauli Hamiltonian rebuilt from the stored (N, J2, E0, W); the
filter counts of every design are n x (per-step) + one Rz(2 phi) per pulse.  For N <= 12 the per-step x n estimate is also checked against
a full transpile of a 48-step circuit.  Idempotent.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from hamiltonians import j1j2_hamiltonian, scale_hamiltonian
from core.resources import direct_step_counts, resource_costs
from core.circuits import build_early_abort_circuit


def filter_counts(step, n, n_pulses):
    return dict(cx=step["cx"] * n, rz_nonclifford=step["rz_nonclifford"] * n + n_pulses, rz=step["rz"] * n + n_pulses, depth=step["depth"] * n)


def validate(H, step, N):
    ks = [12, 12, 12, 12]
    rng = np.random.default_rng(3)
    qc, _ = build_early_abort_circuit(H, np.array(ks) * 0.37, rng.uniform(0.02, 0.2, 4), trial_prep=None, trotter_steps=ks, order=1)
    c = resource_costs(qc, flatten=True, verbose=False, optimization_level=1)
    n = sum(ks)
    return dict(n=n, est_cx=step["cx"] * n, full_cx=int(c["cx"]), est_nc=step["rz_nonclifford"] * n + 4, full_nc=int(c["rz_nonclifford"]))


if __name__ == "__main__":
    for f in sorted(glob.glob(str(ROOT / "results" / "data" / "dmrg_N*_eps*.json"))):
        r = json.load(open(f)); i = r["inputs"]
        H = scale_hamiltonian(j1j2_hamiltonian(r["N"], 1.0, r["J2"]), i["E0"], i["W"])
        step = r["per_step"].get("direct") or direct_step_counts(H, r["N"])
        r["per_step"]["direct"] = step
        for d in r["designs"].values():
            if d and "n" in d:
                d["counts_direct"] = dict(filter=filter_counts(step, d["n"], d["n_pulses_used"]))
        if r["N"] <= 12 and "validation_direct" not in r:
            r["validation_direct"] = validate(H, step, r["N"])
        json.dump(r, open(f, "w"), indent=1)
        print(Path(f).name, "cx/step", step["cx"], "ncRz/step", step["rz_nonclifford"])
