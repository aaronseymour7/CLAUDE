"""dephasing.py -- spectral decomposition of the trial state and of the Trotter error (ED, N <= 12).

usage: python trotter2/dephasing.py N J2        (stored guaranteed design, eps = 1e-2; second-order k from trotter2/circuit_N*_J2*.json)

Writes trotter2/dephasing_N{N}_J2{J2}.json with
  levels      distinct (E, S) levels of H_s: energy, total spin S, multiplicity, trial weight p, weight after the full filter, per-pulse cumulative weights
  counts      number of levels / eigenstates carrying more than a threshold of the trial weight (before the filter, after it), participation ratios
  trotter     per eigenlevel error of the post-selected filter state:  ||P_L (psi_trotter - psi_exact)||^2  (unnormalised, same normalisation), for first and second order
              split into  (a) the ground level (amplitude/phase error of the target),  (b) levels that carry filter weight,  (c) levels the exact filter suppresses (leakage)
  floquet     quasi-energy shifts of one Trotter step e^{-i tau H_eff} (ancilla-0 branch) vs the exact levels, for several tau: ground state and low levels,
              with fitted exponents (first order: diagonal shift ~ tau^2, no O(tau) term because H is real; second order: ~ tau^4)
"""
import sys, os, json, time
from pathlib import Path
import numpy as np, scipy.sparse as sp
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT), str(ROOT.parent), str(ROOT / "s2proj")]
os.chdir(ROOT)
import pipeline as P
from s2_ops import s2_sparse
from core.trotter import _trotter_run
from qiskit.quantum_info import SparsePauliOp

N, J2 = int(sys.argv[1]), float(sys.argv[2])
d = json.load(open(ROOT / "results/data" / f"dmrg_N{N}_J2{J2:g}_eps0.01.json")); g = d["designs"]["guaranteed"]
tg, ph, k1 = np.array(g["tg"]), np.array(g["phases"]), np.array(g["k"], int)
k2 = np.array(json.load(open(HERE / f"circuit_N{N}_J2{J2:g}.json"))["k2"], int)
ctx = P.build_ctx(N, 1.0, J2); Hsp = ctx["H_scaled"]; psi = np.asarray(ctx["trial_vec"], complex); psi /= np.linalg.norm(psi)
dim = 1 << N
Hm = ctx["Hs"].toarray()                                        # scaled H, spectrum in [0, 1]
S2 = s2_sparse(N).toarray()
dl = 1e-6                                                       # splits S sectors inside degenerate energy levels
E_, V = np.linalg.eigh(Hm + dl * S2)
S2v = np.einsum("in,ij,jn->n", V.conj(), S2, V).real
E = np.einsum("in,ij,jn->n", V.conj(), Hm, V).real             # unperturbed energies of those eigenvectors
Sval = np.round((-1 + np.sqrt(1 + 4 * np.maximum(S2v, 0))) / 2).astype(int)
c = V.conj().T @ psi; p = np.abs(c) ** 2

# distinct levels: same energy (1e-8) and same S
key = [(round(E[i], 8), int(Sval[i])) for i in range(dim)]
lev = {}
for i, kk in enumerate(key):
    lev.setdefault(kk, []).append(i)
keys = sorted(lev); idx = [np.array(lev[kk]) for kk in keys]
Lv = np.array([kk[0] for kk in keys]); LS = np.array([kk[1] for kk in keys]); mult = np.array([len(x) for x in idx])
pL = np.array([p[x].sum() for x in idx])

def filt(Es, upto=None):
    f = np.ones_like(Es)
    for t, ph_ in list(zip(tg, ph))[:upto]:
        f = f * np.cos(Es * t + ph_)
    return f
Fexact = filt(E)
pull = [np.array([(np.abs(c[x]) ** 2 * filt(E[x], m) ** 2).sum() for x in idx]) for m in range(len(tg) + 1)]   # unnormalised weights after m pulses
pfin = pull[-1]; psucc = float(pfin.sum())
norm = lambda w: w / w.sum()
def counts(w, name):
    w = norm(w); o = {}
    for th in (1e-2, 1e-3, 1e-4, 1e-6, 1e-8):
        o[f"levels>{th:g}"] = int((w > th).sum()); o[f"states>{th:g}"] = int(sum(mult[w > th]))
    o["participation_levels"] = float(1 / (w ** 2).sum())
    o["levels_for_99pct"] = int(np.searchsorted(np.cumsum(np.sort(w)[::-1]), 0.99) + 1)
    o["levels_for_99.99pct"] = int(np.searchsorted(np.cumsum(np.sort(w)[::-1]), 0.9999) + 1)
    return o
out = dict(N=N, J2=J2, dim=dim, n_levels=len(keys), gamma=float(pL[0]), gap_scaled=float(Lv[1] - Lv[0]), psucc_exact=psucc,
           pulses=dict(t=tg.tolist(), phi=ph.tolist(), k1=k1.tolist(), k2=k2.tolist()),
           levels=dict(E=Lv.tolist(), S=LS.tolist(), mult=mult.tolist(), p_trial=pL.tolist(), p_after_filter=norm(pfin).tolist(),
                       p_unnorm_after_pulse=[w.tolist() for w in pull]),
           weight_by_S={int(s): float(pL[LS == s].sum()) for s in sorted(set(LS))},
           counts_trial=counts(pL, "trial"), counts_after_filter=counts(pfin, "filtered"))
# which levels have weight at the first pulses (the part that can dephase)
out["weight_above_gap"] = float(pL[Lv > Lv[0] + 0.5 * out["gap_scaled"]].sum())
print(f"N={N} J2={J2}: {len(keys)} distinct levels, gamma={pL[0]:.4f}, levels>1e-4 trial: {out['counts_trial']['levels>0.0001']} -> after filter {out['counts_after_filter']['levels>0.0001']}", flush=True)

# ---- Trotter error per level (post-selected, unnormalised)
terms_m = [(sp.csr_matrix(SparsePauliOp("Z" + lab).to_matrix(sparse=True)), float(np.real(cc))) for lab, cc in zip(Hsp.paulis.to_labels(), Hsp.coeffs)]
hd = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
def run_order(order, k):
    state = psi.copy(); amp = 1.0
    for t_i, phi_i, k_i in zip(tg, ph, k):
        s = np.stack([state, np.zeros(dim, complex)]); rz = np.diag([np.exp(-1j * phi_i), np.exp(1j * phi_i)])
        flat = ((rz @ hd) @ s).reshape(-1)
        for _ in range(int(k_i)):
            for seq in ((terms_m,) if order == 1 else (terms_m, terms_m[::-1])):
                for Pm, cc in seq:
                    th = cc * t_i / k_i / (1 if order == 1 else 2)
                    flat = np.cos(th) * flat - 1j * np.sin(th) * (Pm @ flat)
        a0 = (hd @ flat.reshape(2, dim))[0]; pr = float(np.vdot(a0, a0).real); amp *= np.sqrt(pr); state = a0 / np.sqrt(pr)
    return state * amp                                                  # unnormalised post-selected vector
exact_vec = V @ (c * Fexact)
out["trotter"] = {}
for order, k in ((1, k1), (2, k2)):
    t0 = time.time(); vt = run_order(order, k)
    lam = np.vdot(exact_vec, vt) / np.vdot(exact_vec, exact_vec)        # overall scalar (post-selection amplitude) does not change the normalised state
    dv = V.conj().T @ (vt - lam * exact_vec)
    errL = np.array([(np.abs(dv[x]) ** 2).sum() for x in idx]) / float(np.vdot(vt, vt).real); tot = float(errL.sum())   # = (state error)^2 resolved by level, normalised state
    wL = pfin
    grp = dict(ground=float(errL[0]), carrying_weight=float(errL[(wL / wL.sum() > 1e-6) & (np.arange(len(keys)) > 0)].sum()),
               suppressed=float(errL[(wL / wL.sum() <= 1e-6) & (np.arange(len(keys)) > 0)].sum()))
    # relative phase of the ground-level amplitude, and of the excited levels that carry weight, after the full Trotter run
    ph_gs = float(np.angle(np.vdot(V[:, idx[0][0]], vt) / np.vdot(V[:, idx[0][0]], exact_vec)))
    amp_ratio = float(np.abs(np.vdot(V[:, idx[0][0]], vt)) / np.abs(np.vdot(V[:, idx[0][0]], exact_vec)))
    out["trotter"][f"order{order}"] = dict(n=int(k.sum()), err_norm2_total=tot, err_by_level=errL.tolist(), err_groups=grp,
                                            ground_phase_error_rad=ph_gs, ground_amp_ratio=amp_ratio, state_err_trace=float(np.sqrt(max(0, 1 - abs(np.vdot(vt, exact_vec)) ** 2 / (np.vdot(vt, vt).real * np.vdot(exact_vec, exact_vec).real)))), wall_s=time.time() - t0)
    print(f"order {order} (n={int(k.sum())}): state-error^2={tot:.3e} (trace dist {np.sqrt(tot):.2e})  ground {grp['ground']:.2e}  carrying-weight {grp['carrying_weight']:.2e}  suppressed-levels {grp['suppressed']:.2e}", flush=True)

# ---- Floquet / quasi-energy shifts of one Trotter step (ancilla-0 branch, system only), dense
def step_matrix(tau, order):
    M = np.eye(dim, dtype=complex)
    seqs = ((terms_m,) if order == 1 else (terms_m, terms_m[::-1]))
    Hterms = [(sp.csr_matrix(SparsePauliOp(lab).to_matrix(sparse=True)), float(np.real(cc))) for lab, cc in zip(Hsp.paulis.to_labels(), Hsp.coeffs)]
    for seq in ((Hterms,) if order == 1 else (Hterms, Hterms[::-1])):
        for Pm, cc in seq:
            th = cc * tau / (1 if order == 1 else 2)
            M = np.cos(th) * M - 1j * np.sin(th) * (Pm @ M)
    return M
if N <= 10:
    taus = [0.4, 0.2, 0.1, 0.05]; Esorted = np.sort(E); fl = {}
    for order in (1, 2):
        rows = []
        for tau in taus:
            lam = np.linalg.eigvals(step_matrix(tau, order)); eps = np.sort((-np.angle(lam) / tau))
            # E ranges within [0,1] and tau*E << pi: the principal branch is unambiguous; the identity term only shifts E by a constant, sort-match
            dE = eps - Esorted
            rows.append(dict(tau=tau, gs=float(dE[0]), median_low50=float(np.median(np.abs(dE[:50]))), max_abs=float(np.max(np.abs(dE))), mean=float(dE.mean())))
        gs = np.array([abs(r["gs"]) for r in rows]); tt = np.array(taus)
        fl[f"order{order}"] = dict(rows=rows, exponent_gs=float(np.polyfit(np.log(tt[-3:]), np.log(np.maximum(gs[-3:], 1e-300)), 1)[0]),
                                   exponent_max=float(np.polyfit(np.log(tt[-3:]), np.log([r["max_abs"] for r in rows[-3:]]), 1)[0]))
        print(f"floquet order {order}: |dE_gs| = {gs}, exponent {fl[f'order{order}']['exponent_gs']:.2f}; max|dE| exponent {fl[f'order{order}']['exponent_max']:.2f}", flush=True)
    out["floquet"] = fl
json.dump(out, open(HERE / f"dephasing_N{N}_J2{J2:g}.json", "w"), indent=1, default=float)
print("DONE", flush=True)
