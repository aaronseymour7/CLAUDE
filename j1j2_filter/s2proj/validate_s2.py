"""validate_s2.py -- checks behind study.py:  (V1) block-aware Trotter bound holds for the S^2 evolution; (V2) RZZ per S^2 Trotter step by transpile;
(V3) end-to-end exact simulation: S^2 pulses then energy pulses on the L=1 trial state (no Trotter): J>0 weight, infidelity, success probability."""
import sys, os
from pathlib import Path
import numpy as np, scipy.linalg as sl
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
from s2_ops import *                                                     # noqa: E402
from qiskit import QuantumCircuit                                        # noqa: E402
from core.circuits import apply_filter_pulse                              # noqa: E402
from core.resources import resource_costs, RZZ_BASIS                      # noqa: E402

print("V1: Trotter error of exp(-i t S^2) (block order) vs certified alpha t^2/(2k)")
for N in (6, 8):
    tri, ex = alpha_blocks(N); H = s2_pauli_op(N); terms = [(P, q, c) for _, q, P, c in s2_terms(N)]
    from qiskit.quantum_info import SparsePauliOp
    mats = [SparsePauliOp.from_sparse_list([(P + P, list(q), c)], num_qubits=N).to_matrix() for P, q, c in terms]
    Hd = sum(mats)
    for t in (np.pi / 4, np.pi / 8):
        Uex = sl.expm(-1j * t * Hd)
        for k in (1, 2, 4, 8):
            step = np.eye(1 << N, dtype=complex)
            for M in mats:
                step = sl.expm(-1j * (t / k) * M) @ step
            U = np.linalg.matrix_power(step, k); err = np.linalg.norm(U - Uex, 2); bnd = ex * t ** 2 / (2 * k)
            print(f"  N={N} t={t / np.pi:.3f}pi k={k}: error {err:.3e} <= bound {bnd:.3e}  {'OK' if err <= bnd else 'VIOLATED'}  (ratio {err / bnd:.2f})")
print("V2: RZZ per controlled S^2 Trotter step (level 3)")
for N in (6, 8):
    H = s2_pauli_op(N); qc = QuantumCircuit(N + 1); apply_filter_pulse(qc, list(range(N)), N, H, 1.0, 0.0, 1, 1)
    c = resource_costs(qc, verbose=False, basis=RZZ_BASIS); print(f"  N={N}: terms {len(H)}, RZZ per step {c['rzz']:.0f} (model 2 x terms = {2 * len(H)}), non-Clifford rotations {c['nc_total']:.0f}")
print("V3: end-to-end exact (no Trotter) on the L=1 trial state, open chain N=8, J2=0.2411")
import pipeline as P; from pbc_ctx import build_ctx; from core.simulate import apply_pulse_exact; from s2_study import candidates
ctx = build_ctx(8, 1.0, 0.2411, pbc=False); e = ctx["ed"]; W = ctx["W"]; N = 8
gap_s = (e["E_singlet_true"] - e["E0"]) / W; S2 = s2_sparse(N); wJ = j_weights(ctx["trial_vec"], N); pS = wJ[0]; gamma = ctx["gamma"]


def run(pulses, v):
    state = np.concatenate([v, np.zeros(1 << N, complex)]); p = 1.0
    for Hm, t, ph in pulses:
        state = apply_pulse_exact(state, Hm, t, ph); a0 = state[:1 << N]; q = float(np.vdot(a0, a0).real); p *= q; state = np.concatenate([a0 / np.sqrt(q), np.zeros(1 << N, complex)])
    return p, state[:1 << N]


rows, _ = __import__("floor").floor_vs_precision(eps_list=[1e-3], gamma=min(gamma / pS, 1 - 1e-9), delta=gap_s, hi=1.0, m_list=[4, 6, 8], n_jobs=1, verbose=False, **P.FLOOR_KW)
r = [x for x in rows if x.get("feasible")][0]
Epulses = [(ctx["Hs"], t, ph) for t, ph in zip(r["times"], r["phases"])]; Spulses = [(S2, t, 0.0) for t in zero_pulses(N)]
for name, pulses in (("S^2 stage only", Spulses), ("S^2 + energy stage", Spulses + Epulses)):
    p, v = run(pulses, ctx["trial_vec"]); w = j_weights(v, N)
    F = abs(np.vdot(e["psi"], v)) ** 2 / np.vdot(v, v).real
    print(f"  {name}: success prob {p:.4f}, J=0 weight {w[0]:.12f}, weights J>0 max {max(w[j] for j in w if j > 0):.1e}, 1-F = {1 - F:.2e}")
print(f"  prediction: p_S = {pS:.4f}; p_g = p_S * gamma_B * F(0)^2 = {pS * (gamma / pS) * abs(np.prod(np.cos(r['phases']))) ** 2:.4f}; leakage bound for the energy stage: eta = {r['eta_cert']:.3f} -> leak <= {(1 - gamma / pS) * r['eta_cert'] ** 2 / (gamma / pS + (1 - gamma / pS) * r['eta_cert'] ** 2):.2e}")
