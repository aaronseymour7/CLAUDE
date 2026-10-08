"""
dmrg_inputs.py -- the classical front end of the workflow (DMRG), run in the python-3.12 venv (quimb + qiskit + mps-to-circuit).

  dmrg_reference(N, J2)      -> E0, E1 (penalty method), E_top, and the high-bond-dimension ground-state MPS
  trial_dmrg(N, J2, chi)     -> DMRG state at maximum bond dimension chi and its exact circuit (mps-to-circuit "exact")
  trial_approx(N, J2, L)     -> L-layer "approximate" brickwork circuit compiled from the reference MPS
Statevectors are returned in the site order used by workflow.py (site 0 = most significant bit).
"""
import sys, pathlib, io, contextlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "j1j2_filter"))
sys.path.insert(0, str(HERE))
import quimb.tensor as qtn                                   # noqa: E402
from qiskit import transpile                                 # noqa: E402
from qiskit.quantum_info import Statevector                  # noqa: E402
from mps_to_circuit import mps_to_circuit                    # noqa: E402
from hamiltonians import MPO_ham_j1j2, reorder_axes          # noqa: E402
import get_energies as GE                                    # noqa: E402

BASIS = ["cx", "rz", "sx", "x"]


def _quiet(f, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return f(*a, **k)


def dense(mps):
    return np.asarray(mps.to_dense()).reshape(-1)


def dmrg_reference(N, J2, j1=1.0):
    """DMRG E0, E1 (penalty method), E_top and the reference ground-state MPS (bond dim up to 64)."""
    H = MPO_ham_j1j2(N, j1, J2)
    r = _quiet(GE.get_energies, H, N=N, j1=j1, j2=J2, run_ed_check=False)
    return dict(E0=float(r["E0"]), E1=float(r["E1"]), Etop=float(r["E_top"]), psi0=r["psi0"], H=H,
                chi_ref=max(max(t.shape) for t in r["psi0"].arrays))


def circuit_state(qc):
    return reorder_axes(Statevector(qc).data)


def cx_count(qc):
    return int(transpile(qc, basis_gates=BASIS, optimization_level=3, seed_transpiler=0).count_ops().get("cx", 0))


def trial_dmrg(N, J2, chi, j1=1.0):
    """Genuine DMRG at maximum bond dimension chi, compiled exactly (sequential isometries)."""
    H = MPO_ham_j1j2(N, j1, J2)
    E, psi = _quiet(GE.run_dmrg, H, bond_dims=[chi], cutoff=0.0)
    psi.normalize()
    qc = mps_to_circuit(psi.arrays, method="exact", shape="lpr")
    v = circuit_state(qc)
    ov = abs(np.vdot(dense(psi), v)) ** 2
    return dict(E=float(E), psi=v / np.linalg.norm(v), cx=cx_count(qc), mps_overlap=float(ov),
                chi=max(max(t.shape) for t in psi.arrays))


def trial_approx(psi_ref, L):
    qc = mps_to_circuit(psi_ref.arrays, method="approximate", shape="lpr", num_layers=L)
    v = circuit_state(qc)
    return dict(psi=v / np.linalg.norm(v), cx=cx_count(qc), L=L)
