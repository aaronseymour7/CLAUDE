"""sym_cost.py -- cost of the energy filter started from the SU(2)-symmetric trial state (p_S = 1, gap = singlet gap, NO S^2 stage) vs the baseline (A) of s2_study.
Reads results/s2_*_eps*.json (baseline A, gaps, alpha_E, W), optimises the symmetric trial (sym_trial.py, 2 layers by default), designs the filter with the
pipeline's floor design for (gamma_sym, Delta_singlet) and reports RZZ per success.   usage: python s2proj/sym_cost.py [--layers 2] [--worker i --of K]"""
import argparse, glob, json, os, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
from pbc_ctx import ed_inputs                                    # noqa: E402
from sym_trial import optimise                                   # noqa: E402
from s2_study import candidates                                  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--layers", type=int, default=2); ap.add_argument("--worker", type=int, default=0); ap.add_argument("--of", type=int, default=1)
a = ap.parse_args(); out = HERE / "results"
files = sorted(glob.glob(str(out / "s2_*_N*_eps*.json")), key=lambda f: -json.load(open(f))["N"])
for i, f in enumerate(files):
    if i % a.of != a.worker:
        continue
    r = json.load(open(f)); N, J2, pbc, eps = r["N"], r["J2"], r["pbc"], r["eps"]
    g = out / f"sym_{'pbc' if pbc else 'obc'}_N{N}_J2{J2:g}_eps{eps:g}_L{a.layers}.json"
    if g.exists():
        continue
    psi = ed_inputs(N, 1.0, J2, pbc)["psi"]; gam, th = optimise(psi, N, a.layers, pbc)
    gap = r["gap_singlet"]; cE = r["rzz_per_step_E"]; best = None
    cands = candidates(min(gam, 1 - 1e-9), gap, r["alpha_E"], 1.0, eps)
    for c in cands:
        c["n_E"] = r["alpha_E"] * c["T"] ** 2 / 2 / c["allowed"]; c["rzz"] = cE * c["n_E"]; c["rzz_per_success"] = c["rzz"] / c["p_g"]
    if cands:
        best = min(cands, key=lambda c: c["rzz_per_success"])
    res = dict(N=N, J2=J2, pbc=pbc, eps=eps, layers=a.layers, gamma_sym=gam, gamma_baseline=r["gamma"], gap_singlet=gap, gap_triplet=r["gap_triplet"], A=r["A"], B_s2_stage=r["B"], SYM=best)
    g.write_text(json.dumps(res, indent=1))
    A = r["A"]
    print(f"{'pbc' if pbc else 'obc'} N={N} J2={J2:g}: gamma {r['gamma']:.3f} -> sym {gam:.3f}; Delta_s/Delta_t={gap / r['gap_triplet']:.2f}; A {A['rzz'] if A else float('nan'):.2e}  "
          f"S2-stage {r['B']['rzz'] if r['B'] else float('nan'):.2e}  SYM {best['rzz'] if best else float('nan'):.2e}", flush=True)
