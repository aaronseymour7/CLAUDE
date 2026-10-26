"""Pipeline inputs: classical energy estimates and a trial state (ansatz).

Both are plain containers, so any classical method (DMRG, VMC/RBM, coupled cluster, ED,
a previous quantum run) can feed the filter. The pipeline only reads the fields below.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import numpy as np


@dataclass
class EnergyEstimates:
    """What the filter needs to know about the spectrum, in the Hamiltonian's own units.

    E0, E0_err : ground-state energy and an error bar; the pipeline assumes E0 lies in
                 [E0 - E0_err, E0 + E0_err]. A variational estimate (DMRG, VMC) is an upper
                 bound, so E0_err only has to cover the downward side; it is used symmetrically
                 to keep the interface simple.
    E1, E1_err : first excited energy (any symmetry the trial state can overlap with) and its
                 error bar; the pipeline assumes E1 >= E1 - E1_err.
    E_max      : an upper bound on the top of the spectrum. If None, the pipeline uses the
                 rigorous Weyl bound computed from the Hamiltonian's term groups.
    source     : free-text provenance ("DMRG chi<=256", "ED", ...), copied to the report.
    """

    E0: float
    E1: float
    E0_err: float = 0.0
    E1_err: float = 0.0
    E_max: Optional[float] = None
    source: str = ""

    @property
    def gap(self) -> float:
        return self.E1 - self.E0

    def window(self, E_max_bound: float):
        """Scaled window: shift = lower end of the E0 interval, W = E_max - shift.
        Returns dict(shift, W, w0, delta): pass band [0, w0], stop band [delta, 1]."""
        lo = self.E0 - self.E0_err
        E_top = self.E_max if self.E_max is not None else E_max_bound
        W = E_top - lo
        w0 = 2.0 * self.E0_err / W
        delta = (self.E1 - self.E1_err - lo) / W
        if not (0.0 <= w0 < delta < 1.0):
            raise ValueError(f"inconsistent energy window: pass band [0,{w0:.3g}], "
                             f"stop band [{delta:.3g}, 1]")
        return dict(shift=lo, W=W, w0=w0, delta=delta, E_top=E_top)


@dataclass
class Ansatz:
    """A trial state |psi> for the filter.

    gamma       : estimated weight of |psi> on the ground space, |<E0|psi>|^2. The certificate
                  is conditional on this value; pass a conservative (lower) estimate.
    statevector : optional dense vector (big-endian, see pauli.py), only needed for simulation.
    cost        : gate counts of the circuit that prepares |psi>, e.g. {"cnot": 45, "rotations": 90}.
                  Added once per filter attempt in the resource estimate.
    name        : label for reports.
    meta        : anything else (bond dimension, optimiser settings, ...).
    """

    gamma: float
    statevector: Optional[np.ndarray] = None
    cost: Dict[str, float] = field(default_factory=dict)
    name: str = "ansatz"
    meta: Dict = field(default_factory=dict)

    def __post_init__(self):
        if not 0.0 < self.gamma <= 1.0 + 1e-12:
            raise ValueError(f"gamma must be in (0, 1], got {self.gamma}")
        self.gamma = min(float(self.gamma), 1.0)
        if self.statevector is not None:
            v = np.asarray(self.statevector, dtype=complex).reshape(-1)
            self.statevector = v / np.linalg.norm(v)
