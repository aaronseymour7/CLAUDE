"""sym_attrib.py -- attribution of the symmetric-trial gain:  total = (better trial overlap) x (larger gap).
For each case (energy-filter RZZ, certified, eps from the s2_study files):
  A        generic L=1 trial (the pipeline's), gap = min(triplet, singlet)                       [baseline]
  A2       generic L=2 trial (same number of layers as the symmetric ansatz), same gap            [like-for-like baseline]
  SYMg     symmetric 2-layer trial but the FULL-spectrum gap (overlap benefit only; not a valid design for the symmetric state, used for attribution)
  SYM      symmetric 2-layer trial, singlet gap (what the symmetry buys: only valid because the circuit preserves J = 0 exactly)
usage: python s2proj/sym_attrib.py --N 8 12 16 --worker i --of K"""
import argparse, glob, json, os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
from pbc_ctx import build_ctx                                       # noqa: E402
from s2_study import candidates                                     # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--N", type=int, nargs="+", default=[8, 12, 16]); ap.add_argument("--worker", type=int, default=0); ap.add_argument("--of", type=int, default=1)
a = ap.parse_args(); out = HERE / "results"; i = 0
for f in sorted(glob.glob(str(out / "sym_*_L2.json")), key=lambda f: (-json.load(open(f))["N"], f)):
    s = json.load(open(f))
    if s["N"] not in a.N:
        continue
    i += 1
    if (i - 1) % a.of != a.worker:
        continue
    g = out / f"attrib_{'pbc' if s['pbc'] else 'obc'}_N{s['N']}_J2{s['J2']:g}_eps{s['eps']:g}.json"
    if g.exists():
        continue
    r = json.load(open(str(out / f"s2_{'pbc' if s['pbc'] else 'obc'}_N{s['N']}_J2{s['J2']:g}_eps{s['eps']:g}.json")))
    ctx = build_ctx(s["N"], 1.0, s["J2"], pbc=s["pbc"], L=2); gamma2 = ctx["gamma"]; cE = r["rzz_per_step_E"]; aE = r["alpha_E"]
    gapA = min(r["gap_triplet"], r["gap_singlet"]); eps = s["eps"]

    def best(gam, gap):
        c = candidates(min(gam, 1 - 1e-9), gap, aE, 1.0, eps)
        for x in c:
            x["rzz"] = cE * aE * x["T"] ** 2 / 2 / x["allowed"]; x["rps"] = x["rzz"] / x["p_g"]
        return min(c, key=lambda x: x["rps"]) if c else None
    A2, SYMg = best(gamma2, gapA), best(s["gamma_sym"], gapA)
    res = dict(N=s["N"], J2=s["J2"], pbc=s["pbc"], eps=eps, gamma_L1=r["gamma"], gamma_L2_generic=gamma2, gamma_sym=s["gamma_sym"], A=r["A"], A2=A2, SYMg=SYMg, SYM=s["SYM"])
    g.write_text(json.dumps(res, indent=1))
    f_ = lambda x: f"{x['rzz']:.2e}" if x else "infeas"
    print(f"{'pbc' if s['pbc'] else 'obc'} N={s['N']} J2={s['J2']:g}: gamma L1/L2gen/sym = {r['gamma']:.3f}/{gamma2:.3f}/{s['gamma_sym']:.3f} | A {f_(r['A'])}  A2 {f_(A2)}  SYMg {f_(SYMg)}  SYM {f_(s['SYM'])}", flush=True)
