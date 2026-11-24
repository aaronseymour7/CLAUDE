"""Improved energy-pulse design for QPF: same filter, prod_i cos(t_i H), better objective and optimizer.

This is an OPTION next to the authors' design (`qpf.design.optimize`, method "v4"), not a replacement. The circuit
is unchanged: pulses cos(t_i H) with zero phases, one ancilla and one measurement per pulse. Only the times differ.

What the authors' v4 does (qpf.design, reproduced exactly):
  objective   |1 - 1 / sqrt(sum_n prod_i cos^2(t_i E_n))| on E = {0} U linspace(gap, 1, 999)
              = a monotone function of the MEAN of f(E)^2 = prod_i cos^2(t_i E) over [gap, 1] (uniform trial)
  constraints sum t_i = T = pi / gap, 0 <= t_i <= T / 3
  optimizer   SLSQP from one start (halving times), m scanned 4..12
  The phases are passed to SLSQP but v4 ignores them.

What goes wrong, on the scaled design spectrum:
  1. The optimizer is not the bottleneck. Multistart SLSQP on the same objective and constraints lowers the v4 cost
     by at most ~1% (scaled gaps 0.06-0.42). Extra pulses go to t = 0 and the cost stops improving past m = 4-7,
     because the objective and the cap give them nothing to do. Zero-length pulses are still kept, and in the
     circuit each one costs a Trotter step, an ancilla and a measurement (e.g. N = 12, singlet gap 3.2).
     Without the T / 3 cap, the same objective does reach 1.2-1.5x lower v4 costs.
  2. The objective is an average. The mean of f^2 is ~1e-3, but the worst level in the band sits at the gap edge
     with f^2 = 0.03-0.08. On a real spectrum the lowest excited states sit right at that edge and carry most of the
     leakage, so one block suppresses them by only 12-30x, and the pipeline repeats the block (R > 1) to make up.
  3. At T = pi / gap the cap t_i <= T / 3 binds: putting a zero of f at the gap edge needs one pulse of length
     pi / (2 gap) = T / 2. With the cap, even the minimax optimum only reaches max f^2 ~ 0.03.
  4. T is fixed at pi / gap. The suppression a run needs depends on gamma and eps, so the total time is better
     treated as a knob, like the block count R is in the pipeline.

The option here:
  objective   soft-max of log f(E)^2 over [gap, 1] (log-sum-exp with a temperature schedule beta = 2 ... 2048),
              i.e. minimise the WORST-case leakage max_{E in [gap, 1]} f(E)^2 (minimax). It needs no trial and no
              overlap: for any trial with ground-state weight gamma, one block gives
              1 - F <= (1 - gamma) / gamma * max f^2  (exact evolution).
              beta = 1 on the authors' 999-level grid is exactly the v4 objective (monotone transform), which gives
              the "v4, multistart" option: the authors' objective and constraints with a better optimizer.
  constraints sum t_i = T = c pi / gap (c = 1 is the authors' T), 0 <= t_i <= T / 3 (the authors' cap)
  optimizer   SLSQP on the smooth soft-max with analytic gradients, 24 starts per m (halving, equal, and random
              Dirichlet splits), beta raised in stages; the best of m = 4..12 by exact max f^2 on a 20x finer grid,
              ties within 2% going to the smaller m. Pulses shorter than 1e-3 T are dropped.
"""
import numpy as np
from scipy.optimize import minimize

BETAS = (2.0, 8.0, 32.0, 128.0, 512.0, 2048.0)
M_RANGE = range(4, 13)


def log_f2(E, t):
    """log prod_i cos^2(t_i E) per energy, and the per-pulse tan(t_i E) for the gradient."""
    x = np.outer(E, t)
    c2 = np.cos(x) ** 2
    return np.log(np.maximum(c2, 1e-300)).sum(1), np.tan(x)


def softmax_obj(t, E, beta):
    """(1 / beta) log mean_k exp(beta log f(E_k)^2) and its gradient in t. beta = 1: log mean f^2 (v4);
    beta -> inf: log max f^2 (minimax)."""
    ell, tan = log_f2(E, t)
    a = beta * ell
    amax = a.max()
    w = np.exp(a - amax)
    Z = w.sum()
    w /= Z
    val = (amax + np.log(Z) - np.log(len(E))) / beta
    grad = (w[:, None] * (-2.0 * E[:, None] * tan)).sum(0)
    return val, grad


def logmean_obj(t, E, floor):
    """mean_k log(f(E_k)^2 + floor) and its gradient: the log-average (geometric-mean) leakage. Every level counts
    in decibels, so it rewards pushing many levels deep; the floor (e.g. 1e-12) stops exact zeros from going to -inf."""
    ell, tan = log_f2(E, t)
    f2 = np.exp(ell)
    val = np.mean(np.log(f2 + floor))
    w = f2 / (f2 + floor) / len(E)
    grad = (w[:, None] * (-2.0 * E[:, None] * tan)).sum(0)
    return val, grad


def _starts(T, m, n_random, rng):
    out = [1.0 / 2.0 ** np.arange(m), np.ones(m)]
    out += [rng.dirichlet(np.ones(m)) for _ in range(n_random)]
    return [T * s / s.sum() for s in out]


def _solve(E, T, m, cap, betas, starts, fun=softmax_obj):
    bnds = [(0.0, cap)] * m
    cons = [{"type": "eq", "fun": lambda t: np.sum(t) - T, "jac": lambda t: np.ones_like(t)}]
    best = None
    for t in starts:
        t = np.clip(t, 0, cap)
        t = t * T / t.sum()
        for b in betas:
            r = minimize(fun, t, args=(E, b), jac=True, method="SLSQP", bounds=bnds, constraints=cons,
                         options=dict(maxiter=500, ftol=1e-12))
            t = np.clip(r.x, 0, cap)
        val = fun(t, E, betas[-1])[0]
        if best is None or val < best[0]:
            best = (val, t)
    return best[1]


def max_leak(gap, t, n=20001):
    """max_{E in [gap, 1]} prod_i cos^2(t_i E) on a fine grid."""
    return float(np.exp(log_f2(np.linspace(gap, 1, n), np.asarray(t, float))[0].max()))


def mean_leak(gap, t, n_levels=1000):
    """Mean of f^2 over the authors' 999 levels on [gap, 1] (the v4 objective up to a monotone map)."""
    E = np.linspace(gap, 1, n_levels - 1)
    return float(np.mean(np.exp(log_f2(E, np.asarray(t, float))[0])))


def v4_cost(gap, t, n_levels=1000):
    """The authors' v4 cost |1 - 1/sqrt(sum_n f(E_n)^2)| on {0} U linspace(gap, 1, n_levels - 1)."""
    return abs(1 - 1 / np.sqrt(1 + (n_levels - 1) * mean_leak(gap, t, n_levels)))


def design(gap, objective="minimax", c=1.0, m=None, cap_frac=1.0 / 3.0, n_random=22, seed=0, drop=1e-3,
           floor=1e-12):
    """Times for prod_i cos(t_i H) on the scaled spectrum [0] U [gap, 1].

    objective  "minimax" (worst-case leakage over [gap, 1]), "v4" (the authors' objective, multistart; optimised as
               log mean f^2, which has the same minimiser), or "logmean" (mean of log(f^2 + floor) over [gap, 1]).
    c          total time T = c pi / gap (the authors use c = 1).
    m          number of pulses; None scans 4..12 (best by the objective's own figure of merit).
    Returns dict(times (sorted, shortest first), phases (zeros), m, T, max_leak, mean_leak, v4_cost, objective, c).
    """
    T = c * np.pi / gap
    cap = cap_frac * T
    rng = np.random.default_rng(seed)
    if objective == "minimax":
        E = np.linspace(gap, 1, int(max(400, 10 * T * (1 - gap))))
        betas = BETAS
        score = lambda t: max_leak(gap, t)
    elif objective == "v4":
        E = np.linspace(gap, 1, 999)
        betas = (1.0,)
        score = lambda t: mean_leak(gap, t)
    elif objective == "logmean":
        E = np.linspace(gap, 1, int(max(999, 10 * T * (1 - gap))))
        betas = (floor,)
        Ef = np.linspace(gap, 1, 20001)
        score = lambda t: np.exp(logmean_obj(t, Ef, floor)[0])
    else:
        raise ValueError("objective must be 'minimax', 'v4' or 'logmean'")
    fun = logmean_obj if objective == "logmean" else softmax_obj
    best = None
    for mm in ([m] if m else M_RANGE):
        if mm * cap < T - 1e-12:                       # cap makes sum t = T infeasible
            continue
        t = _solve(E, T, mm, cap, betas, _starts(T, mm, n_random, rng), fun)
        t = np.sort(t[t > drop * T])
        t = t * T / t.sum()
        s = score(t)
        if best is None or s < best[0] * (1 - 0.02):
            best = (s, t)
    t = best[1]
    return dict(times=t, phases=np.zeros(len(t)), m=len(t), T=T, c=c, objective=objective,
                max_leak=max_leak(gap, t), mean_leak=mean_leak(gap, t), v4_cost=v4_cost(gap, t),
                logmean_leak=float(np.exp(logmean_obj(t, np.linspace(gap, 1, 20001), floor)[0])))


# ---------------------------------------------------------------- with phases
def _ratio_obj(x, E, beta, m):
    """Soft-max over E of log f(E)^2 - log f(0)^2 for f = prod cos(t_i E + delta_i), and its gradient."""
    t, d = x[:m], x[m:]
    X = np.outer(E, t) + d
    ell = np.log(np.maximum(np.cos(X) ** 2, 1e-300)).sum(1) - np.log(np.cos(d) ** 2).sum()
    a = beta * ell
    amax = a.max()
    w = np.exp(a - amax)
    Z = w.sum()
    w /= Z
    tan = np.tan(X)
    gt = (w[:, None] * (-2 * E[:, None] * tan)).sum(0)
    gd = (w[:, None] * (-2 * tan)).sum(0) + 2 * np.tan(d)
    return (amax + np.log(Z) - np.log(len(E))) / beta, np.concatenate([gt, gd])


def max_leak_ratio(gap, t, d, n=20001):
    """max_{E in [gap, 1]} f(E)^2 / f(0)^2 for f = prod cos(t_i E + delta_i)."""
    E = np.linspace(gap, 1, n)
    f2 = np.prod(np.cos(np.outer(E, t) + d) ** 2, axis=1)
    return float(f2.max() / np.prod(np.cos(d) ** 2))


def design_phases(gap, c=1.0, p0=0.9, m=None, cap_frac=1.0 / 3.0, n_random=12, seed=0, drop=1e-3,
                  m_range=range(4, 11)):
    """Minimax with phases: minimise max_{E in [gap, 1]} f(E)^2 / f(0)^2 subject to f(0)^2 = prod cos^2(delta_i) >= p0
    (the success cost of the phases, 1 / p0 at most), sum t_i = T = c pi / gap, 0 <= t_i <= cap_frac T,
    |delta_i| <= 1.4. Without the p0 constraint the optimum drives f(0) -> 0 (the authors' v3 has the same constraint
    on prod cos^2 delta_i). Returns dict like `design`, plus f0 = prod cos(delta_i) and max_leak = the ratio."""
    T = c * np.pi / gap
    cap = cap_frac * T
    rng = np.random.default_rng(seed)
    E = np.linspace(gap, 1, int(max(400, 10 * T * (1 - gap))))
    best = None
    for mm in ([m] if m else m_range):
        if mm * cap < T - 1e-12:
            continue
        cons = [{"type": "eq", "fun": lambda x, mm=mm: x[:mm].sum() - T,
                 "jac": lambda x, mm=mm: np.r_[np.ones(mm), np.zeros(mm)]},
                {"type": "ineq", "fun": lambda x, mm=mm: np.log(np.cos(x[mm:]) ** 2).sum() - np.log(p0),
                 "jac": lambda x, mm=mm: np.r_[np.zeros(mm), -2 * np.tan(x[mm:])]}]
        bnds = [(0.0, cap)] * mm + [(-1.4, 1.4)] * mm
        starts = [design(gap, "minimax", c, m=mm, cap_frac=cap_frac, n_random=4, seed=seed)["times"]]
        starts += [rng.dirichlet(np.ones(mm)) * T for _ in range(n_random)]
        for t in starts:
            t = np.clip(np.r_[t, np.zeros(mm - len(t))], 0, cap)
            t = t * T / t.sum()
            x = np.r_[t, rng.normal(0, 0.05, mm)]
            for b in BETAS:
                r = minimize(_ratio_obj, x, args=(E, b, mm), jac=True, method="SLSQP", bounds=bnds,
                             constraints=cons, options=dict(maxiter=500, ftol=1e-12))
                x = r.x
            t, d = np.clip(x[:mm], 0, cap), x[mm:]
            if np.prod(np.cos(d) ** 2) < p0 * (1 - 1e-6):
                continue
            keep = t > drop * T
            t, d = t[keep], d[keep]
            s = max_leak_ratio(gap, t, d)
            if best is None or s < best[0] * (1 - 0.02):
                best = (s, t, d)
    s, t, d = best
    o = np.argsort(t)
    t, d = t[o] * T / t.sum(), d[o]
    return dict(times=t, phases=d, m=len(t), T=T, c=c, p0=p0, objective="minimax+phases",
                max_leak=max_leak_ratio(gap, t, d), f0=float(np.prod(np.cos(d))))
