"""sbpipe -- a generic, modular Stetcu-Baroni filter pipeline.

    energy estimates (EnergyEstimates) + trial state (Ansatz) + Hamiltonian (PauliSum)
        -> refine(): scaled window, certified cosine filter, Trotter steps, certificate, resources
        -> validate(): statevector check of the planned circuit (small N)

Front ends that produce the inputs live in sbpipe.frontends (DMRG via quimb is provided).
"""
import os as _os

# The many tiny dense eigenproblems in trotter.py are ~50x slower with multithreaded BLAS.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")
try:                                    # numpy may already be loaded: limit its BLAS pool directly
    from threadpoolctl import threadpool_limits as _tpl
    _tpl(int(_os.environ.get("SBPIPE_BLAS_THREADS", "1")))
except ImportError:
    pass

from .inputs import Ansatz, EnergyEstimates
from .models import j1j2_chain, weyl_bounds
from .pauli import PauliSum
from .pipeline import Config, Plan, certificate, refine, validate
from .resources import GateModel

__all__ = ["Ansatz", "EnergyEstimates", "PauliSum", "j1j2_chain", "weyl_bounds", "Config",
           "Plan", "refine", "validate", "certificate", "GateModel"]
