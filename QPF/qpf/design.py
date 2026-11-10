"""Energy-projection schedules.

Reference implementation: the authors' design notebook (A. Baroni, `reference/optimize_times_phases.ipynb`),
reproduced here exactly as `optimize` (method "v3" times + phases, "v3b", "v4" times only):

  pseudo-spectrum  E = [0, linspace(gap, 1, N - 1)], N = 1000, one ground level, uniform trial state
  objective        C = |1 - prod_i cos(delta_i) / sqrt( sum_n prod_i cos^2(t_i E_n + delta_i) )|
                   (one minus the ground amplitude of the normalised filtered uniform state)
  constraints      sum_i t_i = T (= pi / gap); 0 <= t_i <= T/3; -pi/2 <= delta_i <= pi/2;
                   v3: prod_i cos^2(delta_i) >= 0.9 (success kept within 10% of the zero-phase value)
                   v3b: prod_i cos^2(E_0 t_i + delta_i) >= 0.9 (same when E_0 = 0); v4: no phases
  start            halving times t_{i+1} = t_i / 2 scaled to T, zero phases; SLSQP, maxiter 5000, ftol 1e-12

The paper's own wording (main text: maximise the squared overlap with the ground state; Supplement: 1000 levels
with a fraction alpha_e at zero energy, randomly displaced within a grid spacing, mean squared difference to a
target) is kept as `optimize_supplement` for comparison. The notebook does not use alpha_e or random
displacements in the design: the trial overlap only enters the evaluation.

Comparison schedules from the paper:
  constant  m equal times summing to T ("similar to the Ge et al. algorithm"); m chosen for the lowest energy.
  gaussian  7 times from a normal distribution, scaled so that sum |t_i| = T, averaged over draws.
"""
import numpy as np
from scipy.optimize import minimize

from .core import filter_function


def baroni_spectrum(gap, n_levels=1000):
    """The notebook's design spectrum: 0, then n_levels - 1 evenly spaced levels on [gap, 1]."""
    E = np.zeros(n_levels)
    E[1:] = np.linspace(gap, 1, n_levels - 1)
    return E


def baroni_cost(times, phases, E):
    """The notebook's objective (new_func_v3; new_func_v4 is the same with zero phases)."""
    num = np.prod(np.cos(phases))
    den = np.prod(np.cos(E[:, None] * times + phases) ** 2, axis=1).sum()
    return abs(1 - num / np.sqrt(den))


def halving_guess(T, m):
    t = 1.0 / 2.0 ** np.arange(m)
    return T * t / t.sum()


def optimize(gap, m=6, T=None, method="v3", n_levels=1000, p_min=0.9, maxiter=5000, ftol=1e-12, E=None):
    """The authors' design (see module docstring). Returns dict(times, phases, T, cost, f0, success, E)."""
    T = np.pi / gap if T is None else float(T)
    E = baroni_spectrum(gap, n_levels) if E is None else np.asarray(E, float)
    split = lambda x: (x[:m], x[m:])
    if method == "v4":
        fun = lambda x: baroni_cost(split(x)[0], np.zeros(m), E)
    else:
        fun = lambda x: baroni_cost(*split(x), E)
    cons = [{"type": "eq", "fun": lambda x: abs(np.sum(split(x)[0]) - T)}]
    if method == "v3":
        cons.append({"type": "ineq", "fun": lambda x: np.prod(np.cos(split(x)[1]) ** 2) - p_min})
    elif method == "v3b":
        cons.append({"type": "ineq", "fun": lambda x: np.prod(np.cos(E[0] * split(x)[0] + split(x)[1]) ** 2) - p_min})
    elif method != "v4":
        raise ValueError("method must be 'v3', 'v3b' or 'v4'")
    x0 = np.concatenate([halving_guess(T, m), np.zeros(m)])
    bnds = [(0.0, T / 3.0)] * m + [(-np.pi / 2, np.pi / 2)] * m
    r = minimize(fun, x0, method="SLSQP", bounds=bnds, constraints=cons, options=dict(maxiter=maxiter, ftol=ftol),
                 tol=1e-13)
    t, d = split(r.x)
    d = np.zeros(m) if method == "v4" else d
    return dict(times=t.copy(), phases=d.copy(), T=T, cost=float(r.fun), f0=float(np.prod(np.cos(d))),
                success=bool(r.success), E=E)


def evaluate(times, phases, E, trial):
    """The notebook's evaluation (apply_filter + fdiff): filter the trial amplitudes on the design spectrum,
    normalise, and return (sum (gs - filtered)^2 with gs = level 0, filtered amplitudes)."""
    f = np.asarray(trial, float) * filter_function(E, times, phases)
    f = f / np.sqrt(np.sum(f ** 2))
    gs = np.zeros(len(E)); gs[0] = 1.0
    return float(np.sum((gs - f) ** 2)), f


def notebook_trial(overlap, n_levels):
    """The notebook's evaluation trial: amplitude `overlap` on level 0, the rest spread uniformly."""
    v = np.full(n_levels, np.sqrt((1 - overlap ** 2) / (n_levels - 1)))
    v[0] = overlap
    return v


# ---- the paper's text, read literally (kept for comparison)
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


def optimize_supplement(gap, alpha=0.05, m=6, T=None, phases=True, n_levels=1000, seed=0, restarts=30):
    """The Supplement's wording read literally (MSE to a step, fraction alpha at zero, random displacements).
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
