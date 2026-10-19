"""study.py -- is the S^2 (total-spin) projection stage useful?  Certified cost comparison, 1D J1-J2 chains (open and periodic), exact classical inputs.

A (no projection) : energy filter must resolve the lowest excitation of ANY spin:  Delta_A = min(triplet gap, singlet gap), trial overlap gamma.
B (S^2 stage first): zero-phase S^2 pulses t = pi/4, pi/8, (pi/16)  remove every J > 0 exactly (see s2_ops.zero_pulses).  The energy filter then
                     only sees the S = 0 sector:  Delta_B = first excited SINGLET gap, gamma_B = gamma / p_S  (p_S = trial weight in J = 0).
Both stages are Trotterised with the certificate  1 - F <= (sqrt(leak) + eps_T / sqrt(p_g))^2,  eps_T = sum over ALL pulses alpha t^2 /(2 k),  p_g = gamma F(0)^2
(the S^2 pulses have F = 1 at J = 0, so they cost no success probability beyond p_S).  Optimal split of the eps_T allowance between the stages and
RZZ per Trotter step = 2 x (number of Pauli terms) (qiskit level 3, validated by transpile).
usage: python s2proj/s2_study.py --bc obc pbc --N 8 10 12 --J2 0 0.2411 0.4 --eps 1e-2 [--worker i --of K]"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
import pipeline as P                                                      # noqa: E402
from pbc_ctx import build_ctx                                             # noqa: E402
from s2_ops import s2_sparse, j_weights, zero_pulses, alpha_blocks        # noqa: E402
from floor import floor_vs_precision                                      # noqa: E402
from core.trotter import total_error_bound                                # noqa: E402


def candidates(gamma_design, gap, alpha_E, p_true_factor, eps):
    """Floor-design candidates for leakage budgets eps*FLOOR_EPS_FRACS: (leak, p_g_true, allowed eps_T, T, row).
    gamma_design: overlap used for the leakage/eta target (gamma_B = gamma/p_S for B); p_true_factor: multiplies gamma_design*f0^2 to get the TRUE success probability (p_S for B)."""
    rows, _ = floor_vs_precision(eps_list=[f * eps for f in P.FLOOR_EPS_FRACS], gamma=gamma_design, delta=gap, hi=1.0, m_list=P.FLOOR_M,
                                 n_jobs=P.FLOOR_JOBS, verbose=False, **P.FLOOR_KW)
    out = []
    for r in rows:
        if not r.get("feasible") or r.get("m", 0) == 0:
            continue
        eb = total_error_bound(gamma_design, gap, alpha_E, r["times"], r["phases"], np.ones(len(r["times"])), hi=1.0, eta=r["eta_cert"])
        dT = np.sqrt(eps) - np.sqrt(eb["leak"])
        p_true = p_true_factor * eb["p_g_lb"]
        if dT <= 0 or p_true <= 0:
            continue
        out.append(dict(leak=float(eb["leak"]), p_g=float(p_true), allowed=float(np.sqrt(p_true) * dT), T=float(r["T"]), m=int(r["m"]), eta=float(r["eta_cert"])))
    return out


_STEP_CACHE = HERE / "results" / "s2_step_rzz.json"


def s2_rzz_per_step(N):
    """RZZ per controlled S^2 Lie-Trotter step (block order, qiskit level 3, all-to-all), by transpile; cached."""
    cache = json.load(open(_STEP_CACHE)) if _STEP_CACHE.exists() else {}
    if str(N) not in cache:
        from qiskit import QuantumCircuit
        from core.circuits import apply_filter_pulse
        from core.resources import resource_costs, RZZ_BASIS
        from s2_ops import s2_pauli_op
        H = s2_pauli_op(N); qc = QuantumCircuit(N + 1); apply_filter_pulse(qc, list(range(N)), N, H, 1.0, 0.0, 1, 1)
        c = resource_costs(qc, verbose=False, basis=RZZ_BASIS); cache = json.load(open(_STEP_CACHE)) if _STEP_CACHE.exists() else {}
        cache[str(N)] = dict(rzz=float(c["rzz"]), nc=float(c["nc_total"]), terms=len(H)); _STEP_CACHE.write_text(json.dumps(cache))
    return cache[str(N)]["rzz"]


def study_point(N, J2, pbc, eps):
    t0 = time.time(); ctx = build_ctx(N, 1.0, J2, pbc=pbc); e = ctx["ed"]; W = ctx["W"]
    gap_t = (e["E1_triplet"] - e["E0"]) / W; gap_s = (e["E_singlet_true"] - e["E0"]) / W
    gamma = ctx["gamma"]; alpha_E = ctx["alpha"]; n_terms = len(ctx["H_scaled"]) - 1
    wJ = j_weights(ctx["trial_vec"], N); pS = wJ[0]
    cE = 2.0 * n_terms                                                        # RZZ per energy Trotter step (level 3)
    tri, ex = alpha_blocks(N); aS = ex if ex is not None else tri
    tS = zero_pulses(N); cS = s2_rzz_per_step(N); aS_term = aS * float(tS.sum()) ** 2 / 2        # measured RZZ per S^2 step, a_S = alpha (sum t)^2 / 2
    res = dict(N=N, J2=J2, pbc=pbc, eps=eps, gamma=gamma, p_S=pS, J_weights={str(k): v for k, v in wJ.items()}, gap_triplet=gap_t, gap_singlet=gap_s, W=W, alpha_E=alpha_E,
               alpha_S2=aS, alpha_S2_exact=ex, alpha_S2_triangle=tri, S2_times_over_pi=(tS / np.pi).tolist(), terms_E=n_terms, terms_S2=3 * N * (N - 1) // 2, rzz_per_step_E=cE, rzz_per_step_S2=cS)
    # A: no projection
    candA = candidates(gamma, min(gap_t, gap_s), alpha_E, 1.0, eps)
    bestA = None
    if candA:
        for c in candA:
            a_E = alpha_E * c["T"] ** 2 / 2; c["n_E"] = a_E / c["allowed"]; c["rzz"] = cE * c["n_E"]; c["rzz_per_success"] = c["rzz"] / c["p_g"]
        bestA = min(candA, key=lambda c: c["rzz_per_success"])
    # B: S^2 stage + energy filter on the S = 0 sector
    gB = min(gamma / pS, 1.0 - 1e-9)
    candB = candidates(gB, gap_s, alpha_E, pS, eps)
    bestB = None
    if candB:
        for c in candB:
            a_E = alpha_E * c["T"] ** 2 / 2; sE, sS = np.sqrt(cE * a_E), np.sqrt(cS * aS_term)
            c["n_E"] = a_E / (c["allowed"] * sE / (sE + sS)); c["n_S"] = aS_term / (c["allowed"] * sS / (sE + sS))
            c["rzz_E"] = cE * c["n_E"]; c["rzz_S"] = cS * c["n_S"]; c["rzz"] = c["rzz_E"] + c["rzz_S"]; c["rzz_per_success"] = c["rzz"] / c["p_g"]
        bestB = min(candB, key=lambda c: c["rzz_per_success"])
    res.update(A=bestA, B=bestB, seconds=time.time() - t0)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--bc", nargs="+", default=["obc"]); ap.add_argument("--N", type=int, nargs="+", default=[8]); ap.add_argument("--J2", type=float, nargs="+", default=[0.0])
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-2]); ap.add_argument("--worker", type=int, default=0); ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--outdir", default="results"); a = ap.parse_args(); out = HERE / a.outdir; out.mkdir(exist_ok=True)
    jobs = sorted([(bc, N, J2, eps) for bc in a.bc for N in a.N for J2 in a.J2 for eps in a.eps], key=lambda j: (-j[1], j[0], j[2]))
    for i, (bc, N, J2, eps) in enumerate(jobs):
        if i % a.of != a.worker:
            continue
        f = out / f"s2_{bc}_N{N}_J2{J2:g}_eps{eps:g}.json"
        if f.exists():
            continue
        r = study_point(N, J2, bc == "pbc", eps); f.write_text(json.dumps(r, indent=1))
        A, B = r["A"], r["B"]
        print(f"{bc} N={N} J2={J2:g} eps={eps:g}: gamma={r['gamma']:.3f} p_S={r['p_S']:.3f} Delta_t={r['gap_triplet']:.4f} Delta_s={r['gap_singlet']:.4f} | "
              + (f"A: {A['rzz']:.2e} RZZ (n={A['n_E']:.0f})" if A else "A: infeasible") + " | " + (f"B: {B['rzz']:.2e} RZZ (E: n={B['n_E']:.0f}, S2: n={B['n_S']:.0f})" if B else "B: infeasible")
              + f" ({r['seconds']:.0f}s)", flush=True)
