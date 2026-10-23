"""refine(): energy estimates + ansatz -> certified Stetcu-Baroni filter plan -> resources (-> simulation).

Certificate (proof in the paper, Appendix A; same structure as j1j2_filter/CERTIFICATION.tex):
  ell  = (1-gamma) eta^2 / (gamma + (1-gamma) eta^2)      leakage of the exact filter
  p_g  = gamma F_pass^2                                    lower bound on the success probability
  eps_U = eps_Trotter + R delta                            operator-norm error of the implemented block
  D    <= sqrt(ell) + eps_U / sqrt(p_g)                    trace distance to the ground space
  1 - F <= D^2,   <H> - E0 <= (E_top - E0) D^2.
It holds whenever E0 is in [E0 - E0_err, E0 + E0_err], every other level is >= E1 - E1_err,
E_top <= E_max, and the trial state's ground-space weight is >= gamma.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, Optional, Sequence

import numpy as np

from . import design as D
from .inputs import Ansatz, EnergyEstimates
from .models import weyl_bounds
from .pauli import PauliSum
from .resources import GateModel, filter_counts, total_costs
from .trotter import alpha1, beta2, steps_for_budget, trotter_error


@dataclass
class Config:
    eps: float = 1e-2                       # target infidelity 1 - F
    schedule: str = "optimized"             # uniform | dyadic | optimized | auto (cheaper of the last two)
    order: int = 2                          # Trotter order (1 or 2)
    leak_fracs: Sequence[float] = (0.5, 0.25, 0.1, 0.05)   # leakage share of eps to scan
    synth_frac: float = 0.1                 # share of the unitary-error budget for Rz synthesis
    gate_model: GateModel = field(default_factory=GateModel)
    design_kw: Dict = field(default_factory=dict)


@dataclass
class Plan:
    config: Config
    window: Dict
    design: D.FilterDesign
    steps: np.ndarray
    trotter_const: float
    cert: Dict
    resources: Dict
    H_scaled: PauliSum
    candidates: list

    def summary(self) -> Dict:
        d = self.design
        return dict(
            schedule=d.schedule, order=self.config.order, eps=self.config.eps,
            delta=self.window["delta"], w0=self.window["w0"], W=self.window["W"],
            m=d.m, T_scaled=d.T, T_phys=d.T / self.window["W"], eta=d.eta, f_pass=d.f_pass,
            trotter_const=self.trotter_const, **self.cert,
            **{f"res_{k}": v for k, v in self.resources.items()})


def certificate(gamma, eta, f_pass, eps_U):
    leak = (1 - gamma) * eta ** 2 / (gamma + (1 - gamma) * eta ** 2) if gamma < 1 else 0.0
    p_g = gamma * f_pass ** 2
    Dist = math.sqrt(leak) + eps_U / math.sqrt(p_g)
    return dict(leak=leak, p_g=p_g, eps_U=eps_U, trace_dist_bound=Dist,
                infidelity_bound=min(1.0, Dist ** 2))


def _plan_for(H, Hs, win, des, gamma, cfg, const, trial_cost):
    """Steps, certificate and costs for one design. The synthesis budget covers the trial circuit's
    rotations too: an error r on the trial state moves F(H_s)|psi> by at most r (||F|| <= 1)."""
    eta, fp = des.eta, des.f_pass
    c0 = certificate(gamma, eta, fp, 0.0)
    p_eff = c0["p_g"] if des.m else 1.0          # no filter: nothing is post-selected
    budget = (math.sqrt(cfg.eps) - math.sqrt(c0["leak"])) * math.sqrt(p_eff)
    if budget <= 0:
        return None
    R_trial = trial_cost.get("rotations", 0)
    if des.m == 0:
        steps = np.zeros(0, int)
        eps_T = 0.0
    else:
        steps = steps_for_budget(cfg.order, const, des.times, (1 - cfg.synth_frac) * budget)
        eps_T = trotter_error(cfg.order, const, des.times, steps)
    R = filter_counts(Hs, steps, cfg.order)["rotations"] + R_trial
    delta_syn = (budget - eps_T) / R if R else 0.0   # all remaining budget to synthesis
    eps_U = eps_T + R * delta_syn
    cert = certificate(gamma, eta, fp, 0.0)
    cert["p_g"] = p_eff
    D = math.sqrt(cert["leak"]) + eps_U / math.sqrt(p_eff)
    cert.update(eps_U=eps_U, trace_dist_bound=D, infidelity_bound=min(1.0, D ** 2),
                eps_trotter=eps_T, synth_delta=delta_syn)
    assert cert["infidelity_bound"] <= cfg.eps * (1 + 1e-9)
    res = total_costs(Hs, steps, cfg.order, trial_cost, p_eff, delta_syn, cfg.gate_model)
    return steps, cert, res


_CONST_CACHE: Dict = {}


def trotter_constant(Hs: PauliSum, order: int) -> float:
    """alpha1 (order 1) or beta2 (order 2) of the scaled Hamiltonian, cached."""
    key = (order, tuple(Hs.labels), Hs.without_identity().coeffs.tobytes())
    if key not in _CONST_CACHE:
        _CONST_CACHE[key] = alpha1(Hs) if order == 1 else beta2(Hs)
    return _CONST_CACHE[key]


def refine(H: PauliSum, est: EnergyEstimates, ansatz: Ansatz, cfg: Config = Config()) -> Plan:
    """Design, certify and cost a filter that takes `ansatz` to infidelity <= cfg.eps."""
    lo, hi = weyl_bounds(H)
    win = est.window(hi)
    if est.E_max is None:
        win["E_max_source"] = "Weyl bound"
    Hs = H.shifted_scaled(win["shift"], win["W"])
    const = trotter_constant(Hs, cfg.order)
    best, cands = None, []
    if 1.0 - ansatz.gamma <= cfg.eps:       # the trial state alone meets the target
        des = D.no_filter(win["delta"], win["w0"])
        steps, cert, res = _plan_for(H, Hs, win, des, ansatz.gamma, cfg, const, ansatz.cost)
        cands.append(dict(leak_frac=None, feasible=True, T=0.0, m=0,
                          expected_cnot=res["expected_cnot"], steps=0))
        best = (des, steps, cert, res)
    fracs = [] if ansatz.gamma >= 1 - 1e-12 else cfg.leak_fracs
    scheds = ("dyadic", "optimized") if cfg.schedule == "auto" else (cfg.schedule,)
    for sched, f in [(s, f) for s in scheds for f in fracs]:
        eps_leak = f * cfg.eps
        des = D.design_filter(sched, win["delta"], win["w0"], eps_leak, ansatz.gamma,
                              **cfg.design_kw)
        if des is None:
            cands.append(dict(schedule=sched, leak_frac=f, feasible=False))
            continue
        out = _plan_for(H, Hs, win, des, ansatz.gamma, cfg, const, ansatz.cost)
        if out is None:
            cands.append(dict(schedule=sched, leak_frac=f, feasible=False))
            continue
        steps, cert, res = out
        cands.append(dict(schedule=sched, leak_frac=f, feasible=True, T=des.T, m=des.m,
                          expected_cnot=res["expected_cnot"], steps=res["steps"]))
        if best is None or res["expected_cnot"] < best[3]["expected_cnot"]:
            best = (des, steps, cert, res)
    if best is None:
        raise RuntimeError("no feasible filter design for these inputs")
    des, steps, cert, res = best
    return Plan(cfg, win, des, steps, const, cert, res, Hs, cands)


def validate(plan: Plan, H: PauliSum, ansatz: Ansatz, ground: Optional[np.ndarray] = None) -> Dict:
    """Simulate the planned circuit (exact filter and Trotterised filter) on the ansatz statevector
    and compare with the certificate. Dense; intended for N <= ~14."""
    from scipy.sparse.linalg import eigsh

    from .simulate import filter_exact, filter_trotter, state_metrics, trace_distance_pure

    if ansatz.statevector is None:
        raise ValueError("ansatz has no statevector")
    Hsp = H.to_sparse()
    if ground is None:
        w, v = eigsh(Hsp, k=1, which="SA")
        ground = v[:, 0]
    Hs_sp = plan.H_scaled.to_sparse()
    d = plan.design
    psi = ansatz.statevector
    ve = filter_exact(Hs_sp, psi, d.times, d.phases)
    vt = filter_trotter(plan.H_scaled, psi, d.times, d.phases, plan.steps, plan.config.order)
    me, mt = state_metrics(ve, ground, Hsp), state_metrics(vt, ground, Hsp)
    trot_dist = trace_distance_pure(ve, vt)
    out = dict(
        gamma_true=float(abs(np.vdot(ground, psi)) ** 2),
        infidelity_exact=me["infidelity"], infidelity_trotter=mt["infidelity"],
        p_succ_exact=me["p_succ"], p_succ_trotter=mt["p_succ"],
        energy_trotter=mt["energy"], energy_trial=float(np.vdot(psi, Hsp @ psi).real),
        trotter_state_dist=trot_dist,
        trotter_dist_bound=plan.cert["eps_trotter"] / math.sqrt(plan.cert["p_g"]),   # synthesis not simulated
    )
    out["bound_ok"] = (out["infidelity_trotter"] <= plan.cert["infidelity_bound"] + 1e-12
                       and trot_dist <= out["trotter_dist_bound"] + 1e-12
                       and me["infidelity"] <= plan.cert["leak"] + 1e-12)
    return out
