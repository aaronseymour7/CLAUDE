"""Minimal end-to-end example: DMRG inputs -> certified Stetcu-Baroni filter -> resources -> check.

    python sbpipe_study/example_j1j2.py [N] [J2] [chi] [eps]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sbpipe
from sbpipe.frontends.dmrg_quimb import dmrg_inputs

N = int(sys.argv[1]) if len(sys.argv) > 1 else 12
J2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2411
chi = int(sys.argv[3]) if len(sys.argv) > 3 else 2
eps = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-2

H = sbpipe.j1j2_chain(N, J1=1.0, J2=J2)
est, ansatz, _ = dmrg_inputs(H, chi_trial=chi, bond_dims=(16, 32, 64, 128))
print(f"DMRG: E0 = {est.E0:.10f}  gap = {est.gap:.6f}  trial gamma (chi={chi}) = {ansatz.gamma:.4f}")

plan = sbpipe.refine(H, est, ansatz, sbpipe.Config(eps=eps, schedule="auto", order=2))
s = plan.summary()
print(f"filter: {s['schedule']}, m = {s['m']} pulses, T = {s['T_phys']:.2f}/J1, certified eta = {s['eta']:.3g}")
print(f"certificate: 1-F <= {s['infidelity_bound']:.3g}  (leakage {s['leak']:.2e}, "
      f"Trotter operator error {s['eps_trotter']:.2e}, success prob >= {s['p_g']:.3f})")
print(f"resources: {s['res_qubits']} qubits, {s['res_steps']} Strang steps, "
      f"{s['res_attempt_cnot']:,} CNOT/attempt, {s['res_expected_cnot']:,.0f} CNOT and "
      f"{s['res_expected_t']:.3g} T per success")

if N <= 14:
    v = sbpipe.validate(plan, H, ansatz)
    print(f"simulated: 1-F = {v['infidelity_trotter']:.3g} (trial {1 - v['gamma_true']:.3g}), "
          f"success prob = {v['p_succ_trotter']:.3f}, certificate holds: {v['bound_ok']}")
