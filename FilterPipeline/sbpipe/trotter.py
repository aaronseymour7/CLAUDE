"""Commutator-scaling Trotter error constants (Childs, Su, Tran, Wiebe, Zhu, PRX 11, 011020 (2021)).

For H = sum_a H_a (Pauli terms, stored order) and B_a = sum_{b>a} H_b:

  first order  (Prop. 15)  ||S1(t) - e^{-itH}|| <= alpha1 t^2 / 2,
                           alpha1 = sum_a ||[B_a, H_a]||
  second order (Prop. 16)  ||S2(t) - e^{-itH}|| <= beta2 t^3,
                           beta2 = (1/12) sum_a ||[B_a,[B_a,H_a]]|| + (1/24) sum_a ||[H_a,[H_a,B_a]]||

so k steps over time t give alpha1 t^2/(2k) and beta2 t^3/k^2. Controlling every term on an
ancilla Z (H_a -> H_a (x) Z) leaves both constants unchanged because Z^2 = 1. Identity terms
commute with everything and are dropped.

Every norm is evaluated on the few qubits the commutator acts on: exactly (dense) when that
support has <= 10 qubits, otherwise by the triangle inequality, so the constants are rigorous
upper bounds in every case and the cost is linear in the number of terms for local H.
Both term orders are evaluated and the larger value is returned (the bound then holds for
either order a compiler might realise).
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

from .pauli import PauliSum, apply_rotation, commutes, operator_norm, pauli_mul, support


def _terms(H: PauliSum, reverse: bool):
    h = H.without_identity()
    L, c = list(h.labels), list(h.coeffs)
    if reverse:
        L, c = L[::-1], c[::-1]
    return L, c


def _neighbours(L: List[str]) -> List[List[int]]:
    """For each term, the indices of terms that share a qubit with it (candidate non-commuters)."""
    on_qubit = defaultdict(list)
    for j, lab in enumerate(L):
        for q in support(lab):
            on_qubit[q].append(j)
    out = []
    for j, lab in enumerate(L):
        s = set()
        for q in support(lab):
            s.update(on_qubit[q])
        out.append(sorted(s))
    return out


def _add(acc: Dict[str, complex], lab: str, c: complex):
    acc[lab] = acc.get(lab, 0.0) + c


def alpha1(H: PauliSum, both_orders: bool = True) -> float:
    vals = [_alpha1(H, False)] + ([_alpha1(H, True)] if both_orders else [])
    return max(vals)


def _alpha1(H, reverse):
    L, c = _terms(H, reverse)
    nb = _neighbours(L)
    tot = 0.0
    for a in range(len(L)):
        acc = {}
        for b in nb[a]:
            if b > a and not commutes(L[a], L[b]):
                _add(acc, L[b], c[b])
        if acc:
            tot += 2.0 * abs(c[a]) * operator_norm(acc)[0]   # ||[P_a, Q]|| = 2||Q|| for anticommuting Q
    return tot


def beta2(H: PauliSum, both_orders: bool = True) -> float:
    vals = [_beta2(H, False)] + ([_beta2(H, True)] if both_orders else [])
    return max(vals)


def _beta2(H, reverse):
    L, c = _terms(H, reverse)
    nb = _neighbours(L)
    on_qubit = defaultdict(list)
    for j, lab in enumerate(L):
        for q in support(lab):
            on_qubit[q].append(j)
    s1 = s2 = 0.0
    for a in range(len(L)):
        # [H_a, [H_a, B_a]] = sum_{b>a, anticomm} 4 c_a^2 c_b P_b
        acc2 = {}
        inner = {}                               # [B_a, H_a] = sum_b c_b c_a 2 P_b P_a
        for b in nb[a]:
            if b > a and not commutes(L[a], L[b]):
                _add(acc2, L[b], 4.0 * c[a] ** 2 * c[b])
                ph, Q = pauli_mul(L[b], L[a])
                _add(inner, Q, 2.0 * c[b] * c[a] * ph)
        if acc2:
            s2 += operator_norm(acc2)[0]
        if not inner:
            continue
        # [B_a, inner] = sum_{c>a} sum_Q c_c [P_c, Q]
        acc1 = {}
        for Q, cq in inner.items():
            cand = set()
            for q in support(Q):
                cand.update(on_qubit[q])
            for j in cand:
                if j > a and not commutes(L[j], Q):
                    ph, R = pauli_mul(L[j], Q)
                    _add(acc1, R, 2.0 * c[j] * cq * ph)
        s1 += operator_norm(acc1)[0]
    return s1 / 12.0 + s2 / 24.0


def trotter_error(order: int, const: float, times, steps) -> float:
    """Operator-norm error of the product formula summed over pulses."""
    t = np.abs(np.asarray(times, float))
    k = np.asarray(steps, float)
    if order == 1:
        return float(np.sum(const * t ** 2 / (2.0 * k)))
    if order == 2:
        return float(np.sum(const * t ** 3 / k ** 2))
    raise ValueError("order must be 1 or 2")


def steps_for_budget(order: int, const: float, times, budget: float) -> np.ndarray:
    """Per-pulse step counts k_i = ceil(t_i / dt) with a common dt chosen so that the summed
    error bound is <= budget. A common dt (k_i proportional to t_i) minimises sum k_i for both
    orders (Lagrange / Cauchy-Schwarz). Times are left unchanged, so no re-certification is needed."""
    t = np.abs(np.asarray(times, float))
    T = t.sum()
    if T == 0 or const == 0:
        return np.ones(len(t), dtype=int)
    if order == 1:
        dt = 2.0 * budget / (const * T)          # sum alpha t_i dt / 2 = alpha T dt / 2
    else:
        dt = np.sqrt(budget / (const * T))       # sum beta t_i dt^2 = beta T dt^2
    k = np.maximum(1, np.ceil(t / dt - 1e-12)).astype(int)
    assert trotter_error(order, const, t, k) <= budget * (1 + 1e-9)
    return k


# ----------------------------------------------------------------------------- statevector product formulas
def product_formula(H: PauliSum, psi: np.ndarray, t: float, k: int, order: int) -> np.ndarray:
    """Apply k steps of the first- (Lie) or second-order (Strang) formula for e^{-itH} to psi.
    Lie step: e^{-i dt H_last} ... e^{-i dt H_1} (H_1 first). Identity terms are applied exactly."""
    dt = t / k
    terms = list(zip(H.labels, H.coeffs))
    out = np.array(psi, dtype=complex, copy=True)
    for _ in range(k):
        if order == 1:
            for lab, cf in terms:
                out = apply_rotation(lab, cf * dt, out)
        else:
            for lab, cf in terms:
                out = apply_rotation(lab, cf * dt / 2, out)
            for lab, cf in reversed(terms):
                out = apply_rotation(lab, cf * dt / 2, out)
    return out
