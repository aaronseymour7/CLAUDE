"""
pinsker_check.py -- numerical illustration of the Pinsker remark (python-3.12 venv; uses the designs in results/main_*.json).
Energy-resolved distributions over eigenstates: p_k = |<E_k|psi>|^2 (trial), q_k = p_k F(E_k)^2 / P (ideally filtered):
   KL(q||p) <= ln(1/P),  TV(p,q) <= sqrt(KL/2),  TV(p,q) >= q_0 - p_0 = 1 - l - gamma  =>  P <= exp(-2 (1 - l - gamma)^2)
"""
import sys, json, pathlib
import numpy as np
HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import workflow as wf
import dmrg_inputs as DI

out = []
for N in (4, 6, 8):
    for J2 in (0.0, 0.4):
        d = json.load(open(HERE / "results" / f"main_N{N}_J2{J2}.json"))
        c = [x for x in d["cases"] if abs(x["eps"] - 0.01) < 1e-12][0]
        if "times" not in c:
            continue
        ch = wf.Chain(N, J2, E0_used=d["inputs"]["E0_dmrg"], W_used=d["inputs"]["W_dmrg"])
        w, v = np.linalg.eigh(ch.Hs.toarray())
        psi = DI.trial_dmrg(N, J2, 2)["psi"]
        p = np.abs(v.conj().T @ psi) ** 2
        F = wf.FL.filter_values(c["times"], c["phases"], w)
        P = float(np.sum(p * F ** 2)); q = p * F ** 2 / P
        tv = 0.5 * np.abs(p - q).sum(); m = (q > 0) & (p > 0)
        kl = float(np.sum(q[m] * np.log(q[m] / p[m])))
        cap = float(np.exp(-2 * (q[0] - p[0]) ** 2))
        out.append(dict(N=N, J2=J2, gamma=float(p[0]), ell=float(1 - q[0]), P=P, TV=float(tv), KL=kl,
                        KL_le_ln_invP=bool(kl <= np.log(1 / P) + 1e-12), pinsker_ok=bool(tv <= np.sqrt(kl / 2) + 1e-12),
                        P_cap=cap, P_le_cap=bool(P <= cap + 1e-12)))
        print(out[-1], flush=True)
(HERE / "results" / "pinsker.json").write_text(json.dumps(out, indent=1))
