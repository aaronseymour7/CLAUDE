"""mixing_theory.py -- first-order perturbation theory for the Trotter error of the filter, tested against the measured error.

Lie-Trotter step (term a applied after b for a > b):  V0(tau) = prod_a e^{-i tau H_a} = e^{-i tau H_eff},   H_eff = H + tau X1 + O(tau^2),
    X1 = -(i/2) sum_{a>b} [H_a, H_b]      (BCH).  All H_a are real, so X1 is imaginary Hermitian: it has ZERO diagonal elements in a real eigenbasis,
    i.e. no O(tau) energy shift (confirmed: quasi-energy shift ~ tau^2, dephasing exponent 2.00) -- its O(tau) effect is level MIXING.
Ancilla branches: V1 = conj(V0) = e^{+it(H - tau X1)}.  The block of the pulse is  A~ = (e^{-i phi} V0(t) + e^{i phi} V1(t))/2  vs exact A = cos(Ht + phi).
First-order perturbation theory in the eigenbasis of H (levels E_n, X_mn = <m|X1|n>), for k steps of tau = t/k:
    dV0_mn = -tau X_mn (e^{-iE_n t} - e^{-iE_m t}) / (E_m - E_n),    dV1_mn = +tau X_mn (e^{+iE_n t} - e^{+iE_m t}) / (E_m - E_n),    dA_mn = (e^{-i phi} dV0_mn + e^{i phi} dV1_mn)/2.
|dV_mn| <= 2 tau |X_mn| / |E_m - E_n|: it does NOT grow with t (non-resonant levels), whereas the commutator bound alpha t^2/(2k) = alpha t tau / 2 does.
The final unnormalised error vector is  d psi = sum_i F_{>i}(E_m) [dA_i a_{<i}]_m,  a_{<i} = c_n prod_{j<i} F_j(E_n).

usage: python trotter2/mixing_theory.py N J2      ->  trotter2/mixing_N{N}_J2{J2}.json   (compares with the measured first-order state error of dephasing.py)"""
import sys, os, json
from pathlib import Path
import numpy as np, scipy.sparse as sp
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT), str(ROOT.parent), str(ROOT / "s2proj")]
os.chdir(ROOT)
import pipeline as P
from s2_ops import s2_sparse
from qiskit.quantum_info import SparsePauliOp
N, J2 = int(sys.argv[1]), float(sys.argv[2])
dd = json.load(open(HERE / f"dephasing_N{N}_J2{J2:g}.json")); pul = dd["pulses"]
tg, ph, k1 = np.array(pul["t"]), np.array(pul["phi"]), np.array(pul["k1"], int)
ctx = P.build_ctx(N, 1.0, J2); Hsp = ctx["H_scaled"]; psi = np.asarray(ctx["trial_vec"], complex); psi /= np.linalg.norm(psi); dim = 1 << N
Hm = ctx["Hs"].toarray()
E_, V = np.linalg.eigh(Hm + 1e-6 * s2_sparse(N).toarray()); E = np.einsum("in,ij,jn->n", V.conj(), Hm, V).real
c = V.conj().T @ psi
labs = [(l, float(np.real(x))) for l, x in zip(Hsp.paulis.to_labels(), Hsp.coeffs) if set(l) != {"I"}]
mats = [cc * sp.csr_matrix(SparsePauliOp(l).to_matrix(sparse=True)) for l, cc in labs]
L = sp.csr_matrix((dim, dim), dtype=complex); S = sp.csr_matrix((dim, dim), dtype=complex)
for Ha in mats:                                              # sum_{a>b}[H_a,H_b] = sum_a [H_a, L_{a-1}],  L = running sum of earlier terms
    S = S + (Ha @ L - L @ Ha); L = L + Ha
X1 = (-0.5j) * S; X = V.conj().T @ (X1.toarray() @ V)        # eigenbasis
print(f"max|diag X1| = {np.max(np.abs(np.diag(X))):.2e}  (O(tau) energy shift; expect ~0 up to degenerate-subspace rotation)   max|X1 - X1^dagger| = {np.max(np.abs(X - X.conj().T)):.1e}", flush=True)
dE = E[:, None] - E[None, :]                                  # E_m - E_n   (rows m, cols n)
safe = np.where(np.abs(dE) < 1e-9, 1.0, dE)
def dA_matrix(t, phi, tau):
    em, en = np.exp(-1j * E * t)[:, None], np.exp(-1j * E * t)[None, :]
    d0 = -tau * X * (en - em) / safe;  d1 = tau * X * (np.conj(en) - np.conj(em)) / safe
    zero = np.abs(dE) < 1e-9                                  # degenerate pairs: limit (t-linear); X vanishes between symmetry-equivalent states, keep for completeness
    d0 = np.where(zero, -1j * tau * X * t * em, d0); d1 = np.where(zero, 1j * tau * X * t * np.conj(em), d1)
    return 0.5 * (np.exp(-1j * phi) * d0 + np.exp(1j * phi) * d1)
Fcum = [np.ones(dim)]
for t, p_ in zip(tg, ph):
    Fcum.append(Fcum[-1] * np.cos(E * t + p_))
Ffull = Fcum[-1]
a_exact = c * Ffull; norm_exact = np.linalg.norm(a_exact)
dpsi = np.zeros(dim, complex); per_source = np.zeros(dim)
contrib = []
for i, (t, p_, k) in enumerate(zip(tg, ph, k1)):
    tau = t / k; dA = dA_matrix(t, p_, tau); a_before = c * Fcum[i]
    after = np.ones(dim)
    for j in range(i + 1, len(tg)):
        after = after * np.cos(E * tg[j] + ph[j])
    M = (dA * a_before[None, :]) * after[:, None]            # [m, n]: contribution of source level n to final component m
    dpsi += M.sum(axis=1); contrib.append(float(np.linalg.norm(M.sum(axis=1)) ** 2))
    per_source += np.sum(np.abs(M) ** 2, axis=0)
# remove the component along the exact state (overall scalar), normalise: trace distance^2
lam = np.vdot(a_exact, dpsi) / np.vdot(a_exact, a_exact); perp = dpsi - lam * a_exact
pred_td = float(np.linalg.norm(perp) / norm_exact)
meas_td = dd["trotter"]["order1"]["state_err_trace"]
order = np.argsort(per_source)[::-1]; cs = np.cumsum(per_source[order]) / per_source.sum()
# T-dependence at fixed tau: scale all pulse lengths by lam_T keeping tau = t/k fixed (k scales too) and see how the predicted error moves (bound would grow as lam_T^2 on the same n/T)
out = dict(N=N, J2=J2, predicted_trace_distance=pred_td, measured_trace_distance=meas_td, ratio=pred_td / meas_td, max_abs_diag_X1=float(np.max(np.abs(np.diag(X)))),
           n_source_levels_90pct=int(np.searchsorted(cs, 0.9) + 1), n_source_states_for_99pct=int(np.searchsorted(cs, 0.99) + 1), per_pulse_share=[x / sum(contrib) for x in contrib],
           alpha=float(ctx["alpha"]), bound_trace_distance=float(dd["trotter"]["order1"].get("bound", np.nan)) if "bound" in dd["trotter"]["order1"] else None)
json.dump(out, open(HERE / f"mixing_N{N}_J2{J2:g}.json", "w"), indent=1)
print(f"N={N} J2={J2:g}: first-order perturbation theory predicts trace distance {pred_td:.3e}, measured {meas_td:.3e} (ratio {pred_td / meas_td:.3f}); "
      f"{out['n_source_states_for_99pct']} source eigenstates carry 99% of the predicted error; per-pulse share {np.round(out['per_pulse_share'], 3).tolist()}", flush=True)
