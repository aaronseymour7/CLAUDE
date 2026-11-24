"""Trial states from the DMRG MPS compiled by mps_to_circuit (as in j1j2_filter/pipeline.py), saved as statevectors.
Runs under Python <= 3.12 (mps-to-circuit does not install on 3.13).

    python3.12 j1j2/mps_trial_gen.py --N 8 10 12 --J2 0.3 --layers 1 --out DIR

DMRG: quimb DMRG2 on j1j2_filter.hamiltonians.MPO_ham_j1j2 (open chain, J1 = 1), bond dims 8..64, cutoff 1e-12.
Circuit: mps_to_circuit(psi.arrays, method="approximate", shape="lpr", num_layers=L) (the pipeline's call, L = 1).
Saved: the circuit's statevector (qiskit order), CNOT count after transpiling to {cx, u}, and the DMRG energy."""
import argparse
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "j1j2_filter"))
import quimb.tensor as qtn  # noqa: E402
from hamiltonians import MPO_ham_j1j2  # noqa: E402
from mps_to_circuit import mps_to_circuit  # noqa: E402
from qiskit import transpile  # noqa: E402
from qiskit.quantum_info import Statevector  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", required=True)
    ap.add_argument("--J2", type=float, nargs="+", required=True)
    ap.add_argument("--layers", type=int, default=1)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for N in a.N:
        for J2 in a.J2:
            H = MPO_ham_j1j2(N, j1=1.0, j2=J2)
            d = qtn.DMRG2(H, bond_dims=[8, 16, 32, 64], cutoffs=1e-12)
            d.solve(tol=1e-10, max_sweeps=50, verbosity=0)
            psi = d.state
            psi.normalize()
            qc = mps_to_circuit(psi.arrays, method="approximate", shape="lpr", num_layers=a.layers)
            v = Statevector(qc).data
            cx = transpile(qc, basis_gates=["cx", "u"], optimization_level=1).count_ops().get("cx", 0)
            np.savez(os.path.join(a.out, f"mps_trial_N{N}_J2{J2:g}_L{a.layers}.npz"), psi=v, cx=cx,
                     E_dmrg=float(np.real(d.energy)), psi_mps=np.asarray(psi.to_dense()).reshape(-1))
            print(f"N={N} J2={J2:g}: E_dmrg={float(np.real(d.energy)):.8f} trial CX={cx}", flush=True)


if __name__ == "__main__":
    main()
