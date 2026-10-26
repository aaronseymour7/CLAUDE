"""Total-spin (J^2) projection stage in front of the energy filter (Stetcu-Baroni-Carlson, main text).

An SU(2)-symmetric H commutes with J^2 = (sum_i S_i)^2, whose eigenvalues are k = j(j+1) for
j = 0, 1, ..., N/2 (N even). Removing every j > 0 first leaves the singlet sector, in which the gap to
the first excitation is larger than the triplet gap that limits the plain filter (1.4x to 2.7x for the
J1-J2 chains studied here), so the energy filter that follows is shorter.

Exact projection with a halving ladder. Every k = j(j+1) is an even integer, and
    prod_{r=0}^{R-1} cos(pi k / 2^{r+2}) = sin(pi k / 2) / (2^R sin(pi k / 2^{R+1}))      (k even)
vanishes unless 2^{R+1} divides k. With 2^{R+1} > W_O = (N/2)(N/2 + 1) only k = 0 survives, where
the product is 1. So the smallest R with 2^{R+1} > W_O (9 pulses for N = 48) with times t_r = pi / 2^{r+2} (operator J^2,
phases 0) remove every j > 0 exactly and leave the j = 0 sector untouched: eta_O = 0, F_O(0) = 1.

Combined certificate. The two stages act on commuting operators, so the product F_O(J^2) F_H(H)
multiplies a component (E, j) by F_O(j(j+1)) F_H(E). Components with j > 0 are annihilated, the ground
state keeps |F_H(E0)| >= F_pass, and the remaining components have E in the stop band of the singlet
sector. The certificate of ``pipeline.certificate`` therefore applies unchanged, with eta, F_pass of
the energy stage and with E1 the lowest level of the J = 0 sector. The implementation error of the
two stages adds (telescoping product of ancilla-0 blocks).

Implementation of e^{-i t J^2 (x) Z_a}. J^2 = Sx^2 + Sy^2 + Sz^2 with S_a = (1/2) sum_i sigma^a_i. Each
S_a^2 = N/4 + (1/2) sum_{i<j} sigma^a_i sigma^a_j is a sum of commuting two-body terms, so e^{-i tau S_a^2}
is implemented exactly with N(N-1)/2 Pauli rotations (weight 3 once controlled by the ancilla: 4 CNOTs
and one rotation in the ladder model of resources.py). Only the three non-commuting groups S_x^2, S_y^2,
S_z^2 are Trotterised. Their commutator constants are evaluated exactly: S_a act block-diagonally on the
spin-j irreps (dimension 2j+1 <= N+1), so the operator norm of any polynomial in them equals
the maximum of its norm over j = 0..N/2, which are small dense matrices.

Costs are an upper-bound ladder model; cheaper implementations of e^{-i tau S_z^2} exist (a Hamming-weight
register and one phase) and are not modelled.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache
from typing import Dict, List, Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

from . import design as D
from .inputs import Ansatz, EnergyEstimates
from .models import weyl_bounds
from .pauli import PauliSum, pauli_sparse
from .pipeline import Config, certificate, trotter_constant
from .resources import filter_counts
from .trotter import steps_for_budget, trotter_error


# ----------------------------------------------------------------------------- the projection ladder
def j2_max(N: int) -> int:
    """Largest eigenvalue j(j+1) of J^2 for N spin-1/2 (N even)."""
    if N % 2:
        raise ValueError("the J = 0 projection needs an even number of spins")
    return (N // 2) * (N // 2 + 1)


def projection_times(N: int) -> np.ndarray:
    """Times t_r = pi / 2^{r+2}, r < R, of the exact J = 0 projection for N spins (operator J^2)."""
    R = 1
    while 2 ** (R + 1) <= j2_max(N):          # smallest R with 2^{R+1} > W_O
        R += 1
    return math.pi / 2.0 ** (np.arange(R) + 2)


def projection_values(N: int, times: Optional[Sequence[float]] = None) -> Dict[int, float]:
    """prod_r cos(t_r j(j+1)) for j = 0..N/2 (1 at j = 0, 0 for j > 0 up to rounding)."""
    t = projection_times(N) if times is None else np.asarray(times, float)
    return {j: float(np.prod(np.cos(t * j * (j + 1)))) for j in range(N // 2 + 1)}


# ----------------------------------------------------------------------------- spin-j matrices, constants
def spin_matrices(j: float):
    d = int(round(2 * j + 1))
    m = j - np.arange(d)
    Sz = np.diag(m).astype(complex)
    Sp = np.zeros((d, d), complex)
    for k in range(1, d):
        Sp[k - 1, k] = math.sqrt(j * (j + 1) - m[k] * (m[k] + 1))
    Sx = 0.5 * (Sp + Sp.conj().T)
    Sy = -0.5j * (Sp - Sp.conj().T)
    return Sx, Sy, Sz


def _comm(A, B):
    return A @ B - B @ A


def _alpha_beta_over_j(N: int, order_groups):
    """(alpha1, beta2) for the group order given as a list of axes, e.g. ['x','y','z'], with every
    operator norm taken as the maximum over the spin-j irreps, j = 0..N/2."""
    jmax = N // 2
    a_terms = [0.0] * (len(order_groups) - 1)
    b1_terms = [0.0] * (len(order_groups) - 1)
    b2_terms = [0.0] * (len(order_groups) - 1)
    for j in range(jmax + 1):
        S = dict(zip("xyz", spin_matrices(j)))
        Hs = [S[a] @ S[a] for a in order_groups]
        for a in range(len(Hs) - 1):
            Ba = sum(Hs[a + 1:])
            Ha = Hs[a]
            c = _comm(Ba, Ha)
            a_terms[a] = max(a_terms[a], np.linalg.norm(c, 2))
            b1_terms[a] = max(b1_terms[a], np.linalg.norm(_comm(Ba, c), 2))
            b2_terms[a] = max(b2_terms[a], np.linalg.norm(_comm(Ha, _comm(Ha, Ba)), 2))
    return float(sum(a_terms)), float(sum(b1_terms) / 12.0 + sum(b2_terms) / 24.0)


@lru_cache(maxsize=None)
def spin_trotter_constants(N: int):
    """(alpha1, beta2) of the splitting J^2 = Sx^2 + Sy^2 + Sz^2 for N spins, the larger over the two
    group orders (forward and reversed), exact over all spin sectors."""
    f = _alpha_beta_over_j(N, ["x", "y", "z"])
    r = _alpha_beta_over_j(N, ["z", "y", "x"])
    return max(f[0], r[0]), max(f[1], r[1])


# ----------------------------------------------------------------------------- resources of the O stage
def spin_layer_counts(N: int, model: str) -> Dict[str, int]:
    """Cost of one layer e^{-i tau S_a^2 (x) Z_a} (one axis, exact), ancilla-controlled.

    ``pairs``  : N(N-1)/2 commuting weight-3 Pauli rotations, ladder model: 4 CNOTs + 1 rotation each.
    ``hamming``: compute the Hamming weight w of the N qubits into a b-bit register, b = ceil(log2(N+1)),
                 with 2 (N - b) Toffolis (compute and uncompute; 6 CNOTs and 4 T gates each), apply the
                 phase exp(-i tau (N - 2w)^2 Z_a / 4), which is quadratic in the bits of w and costs b single
                 and b(b-1)/2 pair rotations controlled by Z_a (2 and 4 CNOTs), and b extra ancillas. The
                 Hadamard/phase basis changes for the x and y axes are Cliffords and are not counted."""
    if model == "pairs":
        pairs = N * (N - 1) // 2
        return dict(cnot=4 * pairs, rotations=pairs, toffoli=0, ancillas=0)
    if model == "hamming":
        b = math.ceil(math.log2(N + 1))
        tof = 2 * (N - b)
        pr = b * (b - 1) // 2
        return dict(cnot=6 * tof + 2 * b + 4 * pr, rotations=b + pr, toffoli=tof, ancillas=b)
    raise ValueError(f"unknown model {model!r}")


def spin_stage_counts(N: int, steps, order: int, model: str = "pairs") -> Dict[str, int]:
    """CNOTs, arbitrary rotations and Toffolis of the J^2 pulses. A Lie step has 3 layers; k Strang steps
    have 4k + 1 (the x half-steps of adjacent steps merge); one merged ancilla rotation per pulse."""
    k = np.asarray(steps, int)
    layers = int(np.sum(3 * k)) if order == 1 else int(np.sum(4 * k + 1))
    L = spin_layer_counts(N, model)
    return dict(cnot=L["cnot"] * layers, rotations=L["rotations"] * layers + len(k),
                toffoli=L["toffoli"] * layers, ancillas=L["ancillas"], steps=int(k.sum()))


# ----------------------------------------------------------------------------- planning
@dataclass
class SpinPlan:
    config: Config
    window: Dict
    design: D.FilterDesign
    steps: np.ndarray            # energy stage
    spin_times: np.ndarray
    spin_steps: np.ndarray
    trotter_const: float
    spin_const: float
    budget_share: float          # share of the Trotter budget given to the energy stage
    cert: Dict
    resources: Dict
    H_scaled: PauliSum

    def summary(self) -> Dict:
        d = self.design
        return dict(schedule=d.schedule, order=self.config.order, eps=self.config.eps,
                    delta=self.window["delta"], w0=self.window["w0"], W=self.window["W"], m=d.m,
                    T_scaled=d.T, eta=d.eta, f_pass=d.f_pass, spin_pulses=len(self.spin_times),
                    spin_steps=int(self.spin_steps.sum()), energy_steps=int(self.steps.sum()),
                    budget_share=self.budget_share, trotter_const=self.trotter_const,
                    spin_const=self.spin_const, **self.cert,
                    **{f"res_{k}": v for k, v in self.resources.items()})


def refine_spin(H: PauliSum, est: EnergyEstimates, ansatz: Ansatz, cfg: Config = Config(),
                shares=(0.05, 0.1, 0.25, 0.5, 0.7, 0.85, 0.95, 0.99), model: str = "pairs") -> SpinPlan:
    """Plan a J = 0 projection followed by a certified energy filter. ``est.E1`` must be the lowest
    level of the J = 0 sector (above E0) and the certificate is conditional on it, like the plain pipeline.
    The (1 - synth_frac) share of the error budget that goes to the product formulas is split between the
    two stages on a grid of ``shares``; the split with the fewest expected CNOTs is kept. ``model`` selects the
    cost model of the all-to-all J^2 layers (see ``spin_layer_counts``)."""
    N = H.num_qubits
    lo, hi = weyl_bounds(H)
    win = est.window(hi)
    Hs = H.shifted_scaled(win["shift"], win["W"])
    constH = trotter_constant(Hs, cfg.order)
    tO = projection_times(N)
    a1, b2 = spin_trotter_constants(N)
    constO = a1 if cfg.order == 1 else b2
    gamma = ansatz.gamma
    R_trial = ansatz.cost.get("rotations", 0)
    best = None
    scheds = ("dyadic", "optimized") if cfg.schedule == "auto" else (cfg.schedule,)
    for sched, f in [(s_, f_) for s_ in scheds for f_ in cfg.leak_fracs]:
        des = D.design_filter(sched, win["delta"], win["w0"], f * cfg.eps, gamma, **cfg.design_kw)
        if des is None or des.m == 0:
            continue
        c0 = certificate(gamma, des.eta, des.f_pass, 0.0)
        budget = (math.sqrt(cfg.eps) - math.sqrt(c0["leak"])) * math.sqrt(c0["p_g"])
        if budget <= 0:
            continue
        BT = (1 - cfg.synth_frac) * budget
        for share in shares:
            sH = steps_for_budget(cfg.order, constH, des.times, share * BT)
            sO = steps_for_budget(cfg.order, constO, tO, (1 - share) * BT)
            epsT = trotter_error(cfg.order, constH, des.times, sH) + trotter_error(cfg.order, constO, tO, sO)
            fh = filter_counts(Hs, sH, cfg.order)
            fo = spin_stage_counts(N, sO, cfg.order, model)
            rot = fh["rotations"] + fo["rotations"] + R_trial
            delta_syn = (budget - epsT) / rot
            eps_U = epsT + rot * delta_syn
            cert = dict(c0)
            Dist = math.sqrt(cert["leak"]) + eps_U / math.sqrt(cert["p_g"])
            cert.update(eps_U=eps_U, trace_dist_bound=Dist, infidelity_bound=min(1.0, Dist ** 2),
                        eps_trotter=epsT, synth_delta=delta_syn)
            assert cert["infidelity_bound"] <= cfg.eps * (1 + 1e-9)
            p_succ = (math.sqrt(cert["p_g"]) - eps_U) ** 2
            cert["p_succ_bound"] = p_succ
            t_rot = cfg.gate_model.t_per_rotation(delta_syn)
            attempt_cnot = fh["cnot"] + fo["cnot"] + ansatz.cost.get("cnot", 0)
            res = dict(qubits=fh["qubits"], steps=fh["steps"] + fo["steps"],
                       filter_cnot=fh["cnot"], spin_cnot=fo["cnot"], trial_cnot=ansatz.cost.get("cnot", 0),
                       filter_rotations=fh["rotations"], spin_rotations=fo["rotations"],
                       spin_toffoli=fo["toffoli"], qubits_total=fh["qubits"] + fo["ancillas"],
                       attempt_cnot=attempt_cnot, attempt_rotations=rot, synth_delta=delta_syn,
                       t_per_rotation=t_rot, attempt_t=rot * t_rot + 4 * fo["toffoli"], p_succ=p_succ,
                       expected_cnot=attempt_cnot / p_succ, expected_t=(rot * t_rot + 4 * fo["toffoli"]) / p_succ)
            if best is None or res["expected_cnot"] < best[5]["expected_cnot"]:
                best = (des, sH, sO, share, cert, res)
    if best is None:
        raise RuntimeError("no feasible filter design for these inputs")
    des, sH, sO, share, cert, res = best
    return SpinPlan(cfg, win, des, sH, tO, sO, constH, constO, share, cert, res, Hs)


# ----------------------------------------------------------------------------- simulation (small N)
def _single_qubit_all(psi: np.ndarray, U: np.ndarray, N: int) -> np.ndarray:
    v = psi.reshape((2,) * N)
    for ax in range(N):
        v = np.moveaxis(np.tensordot(U, v, axes=([1], [ax])), 0, ax)
    return v.reshape(-1)


_HAD = np.array([[1, 1], [1, -1]], complex) / math.sqrt(2)
# U sigma^y U^dagger = sigma^z for U = H S^dagger (checked in the tests)
_SDG_H = _HAD @ np.diag([1, -1j])


def _sz2_phases(N: int, tau: float) -> np.ndarray:
    idx = np.arange(2 ** N)
    pc = np.zeros_like(idx)
    for b in range(N):
        pc += (idx >> b) & 1
    sz = 0.5 * (N - 2 * pc)
    return np.exp(-1j * tau * sz ** 2)


def evolve_sa2(psi: np.ndarray, N: int, axis: str, tau: float) -> np.ndarray:
    """e^{-i tau S_a^2} applied to a statevector (exact: diagonal in the rotated basis)."""
    if axis == "z":
        return _sz2_phases(N, tau) * psi
    U = _HAD if axis == "x" else _SDG_H
    return _single_qubit_all(_sz2_phases(N, tau) * _single_qubit_all(psi, U, N), U.conj().T, N)


def spin_product_formula(psi: np.ndarray, N: int, t: float, k: int, order: int) -> np.ndarray:
    """k steps of the Lie (order 1, x first) or Strang (order 2) formula for e^{-i t J^2}."""
    dt = t / k
    v = np.array(psi, complex, copy=True)
    for _ in range(k):
        if order == 1:
            for a in "xyz":
                v = evolve_sa2(v, N, a, dt)
        else:
            for a in "xyz":
                v = evolve_sa2(v, N, a, dt / 2)
            for a in "zyx":
                v = evolve_sa2(v, N, a, dt / 2)
    return v


def spin_filter_trotter(psi, N, times, steps, order):
    """Ancilla-0 block of the Trotterised J^2 pulses: (S(t) + S(-t)) / 2 per pulse (phases 0)."""
    v = np.array(psi, complex, copy=True)
    for t, k in zip(times, steps):
        v = 0.5 * (spin_product_formula(v, N, t, int(k), order) + spin_product_formula(v, N, -t, int(k), order))
    return v


def j2_sparse(N: int) -> sp.csr_matrix:
    out = sp.csr_matrix((2 ** N, 2 ** N), dtype=complex)
    for a in "XYZ":
        S = sum(0.5 * pauli_sparse("I" * q + a + "I" * (N - q - 1)) for q in range(N))
        out = out + S @ S
    return out.tocsr()


def validate_spin(plan: SpinPlan, H: PauliSum, ansatz: Ansatz, ground: np.ndarray) -> Dict:
    """Simulate the planned two-stage circuit on the ansatz statevector (N <= ~14) and compare with the
    certificate. Also reports the weight left outside J = 0 after the (Trotterised) projection."""
    from .simulate import filter_trotter, state_metrics

    N = H.num_qubits
    psi = ansatz.statevector
    J2 = j2_sparse(N)
    Hsp = H.to_sparse()
    d = plan.design
    v1 = spin_filter_trotter(psi, N, plan.spin_times, plan.spin_steps, plan.config.order)
    v2 = filter_trotter(plan.H_scaled, v1, d.times, d.phases, plan.steps, plan.config.order)
    u1 = v1 / np.linalg.norm(v1)
    m = state_metrics(v2, ground, Hsp)
    # exact stages for comparison
    ve1 = psi.astype(complex)
    for t in plan.spin_times:
        ve1 = 0.5 * (sla.expm_multiply(-1j * t * J2, ve1) + sla.expm_multiply(1j * t * J2, ve1))
    from .simulate import filter_exact
    ve2 = filter_exact(plan.H_scaled.to_sparse(), ve1, d.times, d.phases)
    me = state_metrics(ve2, ground, Hsp)
    return dict(
        gamma_true=float(abs(np.vdot(ground, psi)) ** 2),
        j2_trial=float(np.vdot(psi, J2 @ psi).real), j2_after_projection=float(np.vdot(u1, J2 @ u1).real),
        j2_exact_stage=float(np.vdot(ve1 / np.linalg.norm(ve1), J2 @ (ve1 / np.linalg.norm(ve1))).real),
        infidelity_trotter=m["infidelity"], infidelity_exact=me["infidelity"],
        p_succ_trotter=m["p_succ"], p_succ_exact=me["p_succ"],
        bound_ok=bool(m["infidelity"] <= plan.cert["infidelity_bound"] + 1e-12
                      and me["infidelity"] <= plan.cert["leak"] + 1e-12))
