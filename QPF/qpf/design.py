"""Energy-projection schedules (paper, energy filter, Figs. 1-2, and its Supplemental Material).

The paper designs the filter from three pieces of information: the shifted ground state is near 0, an estimate of
the gap Delta_e, and a rough estimate alpha_e of the overlap with the target.

Pseudo-Hamiltonian and initial state (Supplement): 1000 levels with equal amplitudes; a fraction alpha_e of them
(50 of 1000 for alpha_e = 0.05) at zero energy, "the remainder distributed with energies from the gap to one on a
grid but randomly displaced within a grid spacing".

Objective (Supplement): a target function that is one at zero energy and zero from the gap to one; times and
phases are optimised "to produce the smallest mean squared difference between the state produced by the algorithm
acting on the initial state and the target state all at energy zero". With equal amplitudes c this is
    L = mean_k ( f(E_k) - target(E_k) )^2 ,   f(E) = prod_i cos(t_i E + delta_i),
so the zero-energy levels pull f(0) = prod_i cos(delta_i) towards 1 (which is why the phases come out small) and
the excited levels push f towards 0. Main text: the total time is limited to pi/Delta_e (Fig. 1, Heisenberg),
six (t_i, delta_i) for the Heisenberg case, steps applied shortest time first.

Comparison schedules from the paper:
  constant  m equal times summing to T ("similar to the Ge et al. algorithm"); m chosen to give the lowest energy
            at fixed T.
  gaussian  7 times drawn from a normal distribution, scaled so that the sum of |t_i| is T, averaged over draws
            ("Gaussian-sampled times ... Rodeo algorithm"; negative times are allowed, as in the Supplement).

Where the paper is silent we choose (SPEC.md): the displacement is uniform inside the grid cell; the total time
is used in full (sum t_i = T); the optimiser is L-BFGS-B with random restarts, times positive.
"""
import numpy as np
from scipy.optimize import minimize

from .core import filter_function


def pseudo_spectrum(gap, alpha=0.05, n_levels=1000, seed=0):
    """Energies of the Supplement's test spectrum: round(alpha n) levels at 0, the rest on a grid over [gap, 1],
    each displaced at random within its grid cell. Returns (E, target) with target 1 at E = 0, else 0."""
    rng = np.random.default_rng(seed)
    n0 = int(round(alpha * n_levels))
    n1 = n_levels - n0
    spacing = (1.0 - gap) / n1
    E1 = gap + (np.arange(n1) + rng.uniform(0, 1, n1)) * spacing
    E = np.concatenate([np.zeros(n0), E1])
    return E, (E == 0).astype(float)


def mse(E, target, times, phases):
    """The Supplement's objective: mean squared difference between f(E) (equal amplitudes) and the target."""
    return float(np.mean((filter_function(E, times, phases) - target) ** 2))


def ground_fraction(E, target, times, phases):
    """Weight of the zero-energy levels in the normalised filtered state (equal initial amplitudes)."""
    f2 = filter_function(E, times, phases) ** 2
    return float(np.sum(f2 * target) / np.sum(f2))


def _unpack(x, m, T, phases):
    u = x[:m]
    t = np.exp(u - u.max())
    t = T * t / t.sum()                                   # t_i >= 0, sum t_i = T
    d = 0.5 * np.pi * np.tanh(x[m:]) if phases else np.zeros(m)
    return t, d


def optimize(gap, alpha=0.05, m=6, T=None, phases=True, n_levels=1000, seed=0, restarts=30):
    """The paper's optimised schedule (QPF (t_i) with phases=False, QPF (t_i and delta_i) with phases=True).
    Returns dict(times, phases, T, loss, fraction, f0, E, target), times sorted shortest first."""
    T = np.pi / gap if T is None else T
    E, target = pseudo_spectrum(gap, alpha, n_levels, seed)
    rng = np.random.default_rng(seed + 1)

    def loss(x):
        t, d = _unpack(x, m, T, phases)
        return mse(E, target, t, d)

    best = None
    for _ in range(restarts):
        x0 = np.concatenate([rng.normal(0, 1, m), rng.normal(0, 0.05, m) if phases else []])
        r = minimize(loss, x0, method="L-BFGS-B")
        if best is None or r.fun < best.fun:
            best = r
    t, d = _unpack(best.x, m, T, phases)
    order = np.argsort(t)                                 # shortest first, as in the paper's convergence plots
    t, d = t[order], d[order]
    return dict(times=t, phases=d, T=T, loss=float(best.fun), fraction=ground_fraction(E, target, t, d),
                f0=float(np.prod(np.cos(d))), E=E, target=target)


def constant(T, m):
    """m equal times summing to T, zero phases."""
    return np.full(m, T / m), np.zeros(m)


def gaussian(T, m=7, seed=0):
    """m times from N(0, 1), scaled so that sum |t_i| = T, zero phases."""
    t = np.random.default_rng(seed).normal(0, 1, m)
    return T * t / np.abs(t).sum(), np.zeros(m)
