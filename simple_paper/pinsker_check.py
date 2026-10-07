"""
pinsker_check.py -- numerical illustration of the (optional) Pinsker remark.

Energy-resolved distributions over eigenstates:  p_k = |<E_k|psi>|^2   (trial),
q_k = p_k F(E_k)^2 / P_succ  (ideally filtered).  Then
   KL(q||p) = sum q_k ln(F_k^2/P) <= ln(1/P)               (|F| <= 1)
   TV(p,q) <= sqrt(KL(q||p)/2)                              (Pinsker)
   TV(p,q) >= q_0 - p_0 = (1 - l) - gamma
=> P_succ <= exp(-2 (1 - l - gamma)^2): a (weak) cap on how likely post-selection can succeed.
"""
import sys, json, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import workflow as wf

RES = pathlib.Path(__file__).resolve().parent / "results"
out = []
for N, J2 in [(4, 0.0), (6, 0.0), (8, 0.0), (6, 0.4), (8, 0.4)]:
    f = RES / f"case_N{N}_J2{J2}_chi2_eps0.01.json"
    if not f.exists():
        continue
    c = json.loads(f.read_text())
    if "times" not in c:
        continue
    ch = wf.Chain(N, J2)
    w, v = np.linalg.eigh(ch.Hs.toarray())
    psi = wf.mps_truncate(ch.g, N, c["chi"])
    p = np.abs(v.conj().T @ psi) ** 2
    F = wf.FL.filter_values(c["times"], c["phases"], w)
    P = float(np.sum(p * F ** 2))
    q = p * F ** 2 / P
    tv = 0.5 * np.abs(p - q).sum()
    m = (q > 0) & (p > 0)
    kl = float(np.sum(q[m] * np.log(q[m] / p[m])))
    l = 1 - q[0]
    cap = float(np.exp(-2 * (q[0] - p[0]) ** 2))
    out.append(dict(N=N, J2=J2, gamma=float(p[0]), ell=float(l), P=P, TV=float(tv), KL=kl,
                    pinsker_rhs=float(np.sqrt(kl / 2)), KL_le_ln_invP=bool(kl <= np.log(1 / P) + 1e-12),
                    pinsker_ok=bool(tv <= np.sqrt(kl / 2) + 1e-12), P_cap=cap, P_le_cap=bool(P <= cap + 1e-12)))
    print(out[-1])
(RES / "pinsker.json").write_text(json.dumps(out, indent=1))
