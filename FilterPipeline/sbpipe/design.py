"""Stetcu-Baroni cosine filters F(E) = prod_i cos(E t_i + phi_i) on the scaled spectrum.

The scaled Hamiltonian has its (uncertain) ground energy in the pass band [0, w0] and every
other level in the stop band [delta, 1]. A design is certified by
    eta   = sup_{[delta,1]} |F| / F_pass,      F_pass = F(0) - L w0 <= min_{[0,w0]} |F|,
with L = sum |t_i| (|F'| <= L). The sup is a rigorous interval branch-and-bound
(``floor.certify``), so eta is an upper bound, not a grid estimate.

Three schedules are provided:

* ``uniform``  m identical pulses t = pi/(1+delta), phi = 0: F = cos^m(E t). The naive
  repeated-cosine filter; T ~ 2 ln(1/eta) / (pi delta^2).
* ``dyadic``   r repetitions of the halving ladder t = tau, 2 tau, ..., 2^{K-1} tau
  (tau = pi/(1+delta), phi = 0), the time schedule Stetcu-Baroni-Carlson use for quantum-number
  projection, here applied to energy. Each block is a Dirichlet kernel
  sin(2^K E tau) / (2^K sin(E tau)), so T ~ e ln(1/eta) / delta.
* ``optimized`` times and phases optimised for maximal F(0)^2 at certified eta (the
  Stetcu-Baroni-Carlson recipe for energy filters), via the repository's ``floor.py``.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from . import floor as _floor  # certified minimax designer

eta_target = _floor.eta_target
filter_values = _floor.filter_values


@dataclass
class FilterDesign:
    schedule: str
    times: np.ndarray
    phases: np.ndarray
    delta: float
    w0: float
    eta: float            # certified sup_stop |F| / F_pass
    f_pass: float         # certified lower bound on |F| over the pass band
    meta: Dict = field(default_factory=dict)

    @property
    def T(self) -> float:
        return float(np.sum(np.abs(self.times)))

    @property
    def m(self) -> int:
        return len(self.times)

    def values(self, E):
        return filter_values(self.times, self.phases, E)


def certify(times, phases, delta, w0=0.0):
    c = _floor.certify(np.asarray(times, float), np.asarray(phases, float), delta,
                       e0_slack=w0, hi=1.0, rtol=0.005)
    return c["eta"], c["f0_lb"], c


def _make(schedule, times, phases, delta, w0, **meta):
    eta, fp, c = certify(times, phases, delta, w0)
    return FilterDesign(schedule, np.asarray(times, float), np.asarray(phases, float),
                        delta, w0, eta, fp, dict(meta, converged=c["converged"]))


def uniform(delta, w0, eta_t, m_max=10_000_000):
    """cos^m(E t), t = pi/(1+delta). Certified in closed form: for E in [delta, 1], Et lies in
    [theta, pi - theta] with theta = pi delta/(1+delta), so |cos Et| <= cos theta; on the pass band
    |F| >= cos(w0 t)^m (w0 t < pi/2)."""
    t = math.pi / (1.0 + delta)
    c1 = math.cos(math.pi * delta / (1.0 + delta))
    cp = math.cos(w0 * t)
    if cp <= c1:
        return None
    m = max(1, math.ceil(math.log(eta_t) / math.log(c1 / cp)))
    if m > m_max:
        return None
    eta = (c1 / cp) ** m
    return FilterDesign("uniform", np.full(m, t), np.zeros(m), delta, w0, eta, cp ** m,
                        dict(certificate="closed form"))


def dyadic(delta, w0, eta_t, K_max=40):
    """r repetitions of the ladder tau 2^k, k < K, tau = pi/(1+delta); the (K, r) with the smallest T.
    Certified in closed form: one block is sin(2^K E tau)/(2^K sin(E tau)); on [delta, 1],
    E tau lies in [theta, pi - theta], so |block| <= 1/(2^K sin theta). On the pass band
    |F| >= 1 - L w0 with L = sum t_i (|F'| <= L, F(0) = 1)."""
    tau = math.pi / (1.0 + delta)
    sin_t = math.sin(math.pi * delta / (1.0 + delta))
    best = None
    for K in range(1, K_max + 1):
        e_blk = 1.0 / (2 ** K * sin_t)
        if e_blk >= 1.0:
            continue
        Tb = tau * (2 ** K - 1)
        r = max(1, math.ceil(math.log(eta_t) / math.log(e_blk)))      # first guess, ignores w0
        while True:
            L = r * Tb
            fp = 1.0 - L * w0
            if fp <= 0:
                r = None
                break
            if e_blk ** r / fp <= eta_t:
                break
            r += 1
        if r is None:
            continue
        T = r * Tb
        if best is None or T < best[0]:
            best = (T, K, r, e_blk ** r / (1.0 - T * w0), 1.0 - T * w0)
        elif r == 1:
            break                                                        # larger K only gets longer
    if best is None:
        return None
    T, K, r, eta, fp = best
    block = tau * 2.0 ** np.arange(K)
    return FilterDesign("dyadic", np.tile(block, r), np.zeros(K * r), delta, w0, eta, fp,
                        dict(K=K, r=r, certificate="closed form"))


def optimized(delta, w0, eps_leak, gamma, m_list=(4, 6, 8), x_list=None, seed=0):
    """Certified minimax design (floor.py): maximise F(0)^2 subject to certified
    eta <= eta_target(eps_leak, gamma); picks (x = T delta / pi, m) by T^2 / F(0)^2."""
    x_list = x_list or _floor.DEFAULT_X
    rows, _ = _floor.floor_vs_precision([eps_leak], gamma, delta, e0_slack=w0,
                                        x_list=x_list, m_list=m_list, seed=seed,
                                        verbose=False, t_max_frac=0.5)
    r = rows[0]
    if (not r.get("feasible") or r.get("m", 0) == 0) and x_list is _floor.DEFAULT_X:
        # long filters are needed for tiny overlaps or tight targets: widen the search once
        rows, _ = _floor.floor_vs_precision([eps_leak], gamma, delta, e0_slack=w0,
                                            x_list=(3.0, 4.0, 5.0, 6.0, 8.0),
                                            m_list=(8, 10, 12), seed=seed,
                                            verbose=False, t_max_frac=0.5)
        r = rows[0]
    if not r.get("feasible") or r.get("m", 0) == 0:
        return None
    return _make("optimized", r["times"], r["phases"], delta, w0, x=r["x"])


SCHEDULES = ("uniform", "dyadic", "optimized")


def no_filter(delta, w0):
    """Empty filter (F = 1): eta = 1, so the certificate's leakage equals 1 - gamma."""
    return FilterDesign("none", np.zeros(0), np.zeros(0), delta, w0, 1.0, 1.0,
                        dict(note="trial state used as is"))


def design_filter(schedule: str, delta: float, w0: float, eps_leak: float, gamma: float,
                  **kw) -> Optional[FilterDesign]:
    """Design a filter whose exact (un-Trotterised) leakage is <= eps_leak for a trial state
    with ground-space weight gamma. Returns None if no design was found."""
    eta_t = eta_target(eps_leak, gamma)
    if not np.isfinite(eta_t):
        return no_filter(delta, w0)
    if schedule == "uniform":
        return uniform(delta, w0, eta_t, **kw)
    if schedule == "dyadic":
        return dyadic(delta, w0, eta_t, **kw)
    if schedule == "optimized":
        return optimized(delta, w0, eps_leak, gamma, **kw)
    raise ValueError(f"unknown schedule {schedule!r}; choose from {SCHEDULES}")
