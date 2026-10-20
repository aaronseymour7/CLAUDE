"""audit_results.py -- independent re-derivation of every stored result from its stored inputs.

For every results/data/*.json (guaranteed and empirical designs) it re-computes, with a clean re-implementation of the formulas in
CERTIFICATION.tex (no call into the pipeline's bound code):

  W = Weyl bound - E0;  Delta = g/W;  eta (dense-grid maximum <= stored certified eta <= 1.01 x dense);  F(0);
  leakage l = (1-g)eta^2/(g+(1-g)eta^2);  p_g = gamma F(0)^2;  eps_T = sum alpha t_i^2/(2 k_i);  d_T = min(1, eps_T/sqrt(p_g));
  bound = min(1, (sqrt(l) + d_T)^2)  (must equal the stored bound and be <= eps for guaranteed designs);
  sum k_i = n, t_i = k_i dt, T = sum t_i, T_phys = T/W;
  gate counts (RZZ basis; legacy CX basis too) = n x (per-step transpile) + trial circuit (+ per-pulse terms);
  alpha re-evaluated from the Pauli Hamiltonian (exact local commutator sum, both orders);
  every simulated design: measured 1-F <= bound, p_succ >= p_g, trace distance to the exact filter <= eps_T/sqrt(p_exact).

usage: python scripts/audit_results.py [--alpha]     (--alpha recomputes alpha for every point, ~2 min)
Exit status 1 if any check fails.
"""
import glob, json, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
RES = ROOT / "results"
fails, checks = [], 0


def check(name, ok, detail=""):
    global checks
    checks += 1
    if not ok:
        fails.append(f"{name}: {detail}")


def F(t, ph, E):
    return np.prod(np.cos(np.outer(E, t) + ph), axis=1)


do_alpha = "--alpha" in sys.argv
if do_alpha:
    from hamiltonians import j1j2_hamiltonian, scale_hamiltonian
    from core.trotter import alpha_commutator

files = sorted(glob.glob(str(RES / "data" / "dmrg_N*_eps*.json")))
nd = 0
for f in files:
    r = json.load(open(f)); N, J2, eps = r["N"], r["J2"], r["eps"]; i = r["inputs"]; tag = f"N{N} J2={J2:g} eps={eps:g}"
    # ---- inputs
    weyl = (N - 1) / 4 + (N - 2) * J2 / 4
    check(f"{tag} W=Weyl-E0", abs(i["W"] - (weyl - i["E0"])) < 1e-9, f"{i['W']} vs {weyl - i['E0']}")
    check(f"{tag} Etop<=Weyl", i["Etop"] <= weyl + 1e-9, f"{i['Etop']} > {weyl}")
    check(f"{tag} W>=Etop-E0", i["W"] >= i["Etop"] - i["E0"] - 1e-9)
    check(f"{tag} Delta=g/W", abs(i["gap_scaled"] - i["gap_raw"] / i["W"]) < 1e-12)
    check(f"{tag} gap=E1-E0", abs(i["gap_raw"] - (i["E1"] - i["E0"])) < 1e-12)
    check(f"{tag} 0<gamma<=1", 0 < i["gamma"] <= 1 + 1e-9)
    check(f"{tag} dmrg converged", all(x["converged"] for x in i["dmrg_runs"]))
    if do_alpha:
        Hq = j1j2_hamiltonian(N, 1.0, J2)
        Hs = scale_hamiltonian(Hq, i["E0"], i["W"])
        a = max(alpha_commutator(Hs, False), alpha_commutator(Hs, True))
        check(f"{tag} alpha", abs(a - i["alpha"]) < 1e-9 * max(1, a), f"{a} vs {i['alpha']}")
    st = r["per_step"]["a2a"]["one"]
    for which, d in r["designs"].items():
        if not d or "n" not in d:
            continue
        nd += 1
        t = f"{tag} {which}"
        k, tg, ph = np.array(d["k"]), np.array(d["tg"]), np.array(d["phases"])
        check(f"{t} sum k = n", int(k.sum()) == d["n"], f"{k.sum()} vs {d['n']}")
        check(f"{t} t=k dt", np.allclose(tg, k * d["dt"], rtol=1e-9, atol=0))
        check(f"{t} T", abs(d["T_scaled"] - tg.sum()) < 1e-9 * max(1, tg.sum()) and abs(d["T_phys"] - tg.sum() / i["W"]) < 1e-9 * max(1, d["T_phys"]))
        check(f"{t} n_pulses", d["n_pulses_used"] == len(k) and (k > 0).all())
        delta = i["gap_scaled"]
        E = np.linspace(delta, 1.0, 400_001)
        dense = float(np.max(np.abs(F(tg, ph, E))))
        f0 = abs(float(np.prod(np.cos(ph))))
        check(f"{t} F0", abs(f0 - d["f0"]) < 1e-12)
        eta_dense = dense / f0
        check(f"{t} eta sound", d["eta"] >= eta_dense * (1 - 1e-9), f"stored {d['eta']} < dense {eta_dense}")
        check(f"{t} eta tight", d["eta"] <= 1.02 * eta_dense + 1e-9, f"stored {d['eta']} vs dense {eta_dense}")
        g_, eta = i["gamma"], d["eta"]
        leak = (1 - g_) * eta ** 2 / (g_ + (1 - g_) * eta ** 2) if g_ < 1 else 0.0
        pg = g_ * f0 ** 2
        epsT = float(np.sum(i["alpha"] * tg ** 2 / (2.0 * k)))
        dT = min(1.0, epsT / np.sqrt(pg))
        bound = min(1.0, (np.sqrt(leak) + dT) ** 2)
        check(f"{t} leak", abs(leak - d["leak_bound"]) < 1e-10, f"{leak} vs {d['leak_bound']}")
        check(f"{t} p_g", abs(pg - d["p_succ_lb"]) < 1e-10)
        check(f"{t} eps_T", abs(epsT - d["eps_T"]) < 1e-9 * max(1, epsT))
        check(f"{t} d_T", abs(dT - d["d_T"]) < 1e-9)
        check(f"{t} bound", abs(bound - d["eps_bound"]) < 1e-9, f"{bound} vs {d['eps_bound']}")
        if which == "guaranteed":
            check(f"{t} bound<=eps", bound <= eps + 1e-12, f"{bound} > {eps}")
        # ---- counts (opt-level-3 per-step transpile x n)
        c = d["counts_a2a"]["filter"]
        check(f"{t} CX = n x step", abs(c["cx"] - d["n"] * st["cx"]) < 1e-6, f"{c['cx']} vs {d['n'] * st['cx']}")
        check(f"{t} nonCliffRz = n x step + pulses", abs(c["rz_nonclifford"] - (d["n"] * st["rz_nonclifford"] + len(k))) < 1e-6)
        # ---- direct (level-1) compilation: analytic per-step counts and totals
        sd = r["per_step"]["direct"]; nt = i["n_pauli_terms"]
        check(f"{t} direct CX/step = 4 x terms", sd["cx"] == 4 * nt, f"{sd['cx']} vs {4 * nt}")
        check(f"{t} direct ncRz/step = terms+1", sd["rz_nonclifford"] == nt + 1, f"{sd['rz_nonclifford']} vs {nt + 1}")
        cd = d["counts_direct"]["filter"]
        check(f"{t} direct CX = n x step", abs(cd["cx"] - d["n"] * sd["cx"]) < 1e-6)
        check(f"{t} direct ncRz = n x step + pulses", abs(cd["rz_nonclifford"] - (d["n"] * sd["rz_nonclifford"] + len(k))) < 1e-6)
        check(f"{t} direct CX <= level-3 CX x 2 and ncRz < level-3 ncRz", cd["cx"] > c["cx"] and cd["rz_nonclifford"] < c["rz_nonclifford"])
        # ---- RZZ basis (the reported gate basis): level-3 and direct
        sz = r["per_step"]["rzz"]; z3s, zds = sz["a2a"]["one"], sz["direct"]
        z3 = d["counts_rzz_a2a"]["filter"]; zd = d["counts_rzz_direct"]["filter"]
        check(f"{t} RZZ = n x step (level 3)", abs(z3["rzz"] - d["n"] * z3s["rzz"]) < 1e-6)
        check(f"{t} nc rot = n x step + pulses (level 3)", abs(z3["nc_total"] - (d["n"] * z3s["nc_total"] + len(k))) < 1e-6)
        check(f"{t} RZZ direct/step = 4 x terms", zds["rzz"] == 4 * nt, f"{zds['rzz']} vs {4 * nt}")
        check(f"{t} nc rot direct/step = terms+1", zds["nc_total"] == nt + 1 and zds["rzz_nonclifford"] == 0, f"{zds['nc_total']} vs {nt + 1}")
        check(f"{t} RZZ direct = n x step", abs(zd["rzz"] - d["n"] * zds["rzz"]) < 1e-6 and abs(zd["nc_total"] - (d["n"] * zds["nc_total"] + len(k))) < 1e-6)
        check(f"{t} level-3 RZZ/step = 2 x terms, non-Clifford RZZ = terms", z3s["rzz"] == 2 * nt and z3s["rzz_nonclifford"] == nt, f"{z3s['rzz']}, {z3s['rzz_nonclifford']} vs {nt}")
        check(f"{t} direct RZZ = CX (one-to-one)", abs(zd["rzz"] - cd["cx"]) < 1e-6)
        tcz = i["trial_cost_rzz"]; tz = d["counts_rzz_a2a"]["total"]
        check(f"{t} RZZ total = filter + trial", abs(tz["rzz"] - (z3["rzz"] + tcz["rzz"])) < 1e-6 and abs(tz["nc_total"] - (z3["nc_total"] + tcz["nc_total"])) < 1e-6)
        check(f"{t} nc_total = rz_nc + rzz_nc", abs(z3s["nc_total"] - z3s["rz_nonclifford"] - z3s["rzz_nonclifford"]) < 1e-9)
        # ---- simulated designs
        if d.get("simulated"):
            check(f"{t} measured<=bound(guaranteed)", which != "guaranteed" or d["eps_actual"] <= d["eps_bound"] + 1e-12, f"{d['eps_actual']} > {d['eps_bound']}")
            check(f"{t} p_succ>=p_g", d["p_succ_exact"] >= d["p_succ_lb"] - 1e-12, f"{d['p_succ_exact']} < {d['p_succ_lb']}")
            check(f"{t} state_err<=epsT/sqrt(p)", d["state_err"] <= d["eps_T"] / np.sqrt(d["p_succ_exact"]) + 1e-9, f"{d['state_err']} > {d['eps_T'] / np.sqrt(d['p_succ_exact'])}")
            check(f"{t} F_ed+eps=1", abs(d["F_ed"] + d["eps_actual"] - 1) < 1e-9 if np.isfinite(d["F_ed"]) else True)
            if d.get("circ_dev") is not None:
                check(f"{t} circuit==numpy sim", d["circ_dev"] < 1e-9, f"circ_dev {d['circ_dev']}")
        vr = r.get("validation_rzz")
        if vr and "l3" in vr:
            check(f"{tag} RZZ per-step x n vs full transpile (level 3)", vr["l3"]["est_rzz"] == vr["l3"]["full_rzz"], str(vr["l3"]))
            check(f"{tag} RZZ per-step x n vs full transpile (direct)", vr["direct"]["est_rzz"] == vr["direct"]["full_rzz"] and vr["direct"]["est_nc"] == vr["direct"]["full_nc"], str(vr["direct"]))
            check(f"{tag} level-3 nc estimate within 25% of full", abs(vr["l3"]["est_nc"] - vr["l3"]["full_nc"]) / vr["l3"]["full_nc"] < 0.25, str(vr["l3"]))
        vdr = r.get("validation_direct")
        if vdr:
            check(f"{tag} direct per-step x n vs full transpile", vdr["est_cx"] == vdr["full_cx"] and vdr["est_nc"] == vdr["full_nc"], str(vdr))
        v = r.get("validation")
        if v:
            check(f"{tag} per-step x n vs full transpile CX", abs(v["cx_err_pct"]) < 0.5, f"{v['cx_err_pct']}%")

print(f"{len(files)} data files, {nd} designs, {checks} checks")
if fails:
    print(f"FAILED {len(fails)}:")
    for x in fails[:40]:
        print("  ", x)
    sys.exit(1)
print("ALL CHECKS PASSED")
