"""Quantum-number projection (paper, "Projecting on quantum numbers").

All phases are zero. The first time t1 sets the resolution (it zeroes the lowest unwanted value), and each next
time is half the previous one, t_{i+1} = t_i / 2, which zeroes larger and larger values, including ones missed
before. The total time approaches 2 t1. A different target is selected by shifting the operator by a constant
(e.g. J^2 - 12 for J = 3). With exact zeros the step removes the unwanted quantum numbers entirely, keeps the
relative amplitudes inside the target space, and succeeds with probability equal to the target-space weight.

For J^2 (eigenvalues J(J+1)) and target J = 0 the paper takes t1 = pi/4: cos(J(J+1) pi/4) = 0 for
J = 1, 2, 5, 6, 9, 10, ...
"""
import numpy as np


def halving_times(t1, n):
    """t1, t1/2, ..., t1/2^(n-1)."""
    return [t1 / 2 ** k for k in range(n)]


def removes(values, times, tol=1e-12):
    """Which of `values` (eigenvalues of the shifted operator) the zero-phase ladder sends exactly to zero."""
    v = np.asarray(values, float)
    return np.abs(np.prod([np.cos(t * v) for t in times], axis=0)) < tol


def ladder_for(values, t1, n_max=64):
    """Shortest halving ladder from t1 that zeroes every value in `values` (the unwanted eigenvalues of the shifted
    operator; the target has value 0, where every factor is 1). Raises if n_max steps are not enough."""
    for n in range(1, n_max + 1):
        ts = halving_times(t1, n)
        if removes(values, ts).all():
            return ts
    raise ValueError("no exact-zero ladder within n_max steps; choose another t1")


def j2_values(J_target, J_max, half_integer=False):
    """Eigenvalues of J^2 - J_t(J_t+1) for every J != J_target (J = 0..J_max, or 1/2..J_max)."""
    start = 0.5 if half_integer else 0.0
    Js = np.arange(start, J_max + 0.5, 1.0)
    return np.array([J * (J + 1) - J_target * (J_target + 1) for J in Js if abs(J - J_target) > 1e-9])


def j2_ladder(J_target, J_max, t1=np.pi / 4, half_integer=False):
    """Paper's J^2 projection: halving ladder from t1 = pi/4 on J^2 - J_t(J_t+1)."""
    return ladder_for(j2_values(J_target, J_max, half_integer), t1)
