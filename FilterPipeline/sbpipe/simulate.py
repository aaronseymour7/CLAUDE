"""Statevector simulation of the post-selected filter, exact and Trotterised (small N only)."""
from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as sla

from .pauli import PauliSum
from .trotter import product_formula


def filter_exact(Hs_sparse, psi, times, phases):
    """Unnormalised cos(H_s t_m + phi_m) ... cos(H_s t_1 + phi_1) psi."""
    v = np.array(psi, dtype=complex, copy=True)
    for t, ph in zip(times, phases):
        a = sla.expm_multiply(-1j * t * Hs_sparse, v)
        b = sla.expm_multiply(1j * t * Hs_sparse, v)
        v = 0.5 * (np.exp(-1j * ph) * a + np.exp(1j * ph) * b)
    return v


def filter_trotter(Hs: PauliSum, psi, times, phases, steps, order):
    """Ancilla-0 block of the Trotterised circuit: (e^{-i phi} S(t) + e^{i phi} S(-t)) / 2 per pulse,
    S the product formula. This is exactly what the controlled circuit implements."""
    v = np.array(psi, dtype=complex, copy=True)
    for t, ph, k in zip(times, phases, steps):
        a = product_formula(Hs, v, t, int(k), order)
        b = product_formula(Hs, v, -t, int(k), order)
        v = 0.5 * (np.exp(-1j * ph) * a + np.exp(1j * ph) * b)
    return v


def state_metrics(v, ground, H_sparse=None):
    """Success probability, infidelity to the ground state, energy of the normalised state."""
    p = float(np.vdot(v, v).real)
    u = v / np.sqrt(p)
    F = float(abs(np.vdot(ground, u)) ** 2)
    out = dict(p_succ=p, infidelity=max(0.0, 1.0 - F))
    if H_sparse is not None:
        out["energy"] = float(np.vdot(u, H_sparse @ u).real)
    return out


def trace_distance_pure(u, v):
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    return float(np.sqrt(max(0.0, 1.0 - abs(np.vdot(u, v)) ** 2)))
