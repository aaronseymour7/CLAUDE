"""model2d.py -- crossover estimate: certified cosine-filter cost vs sequential-MPS state-preparation cost for the open L x L square-lattice Heisenberg model (J2 = 0).

MEASURED (ed_measure.py, exact diagonalisation of open W x L clusters, N <= 24):  E0, triplet/singlet gaps, Schmidt spectra -> MPS bond dimension chi(W).
EXACT (any size):  the Trotter constant alpha of the scaled Hamiltonian (alpha_commutator), the number of Pauli terms, the Weyl width.
MODELLED (explicit scenarios, flagged):
  * gaps for large L:  triplet (Anderson tower) Delta_t = a_t / N,   singlet Delta_s = c_s / L   (constants from the 4 x L clusters),
  * trial overlap gamma = exp(-kappa N)  (kappa = 0.02 / 0.05 per site; 1D symmetric-trial values were 0.018-0.026, generic L = 1 trial 0.03-0.05),
  * chi(W) ~ exp(a W + b) extrapolated from W = 1..4.
Filter: the pipeline's floor design for (gamma, Delta), certificate eps = 1e-2, RZZ per Trotter step = 2 x terms, cost per success = RZZ / p_g.
  A = no symmetry (resolve the triplet gap);  S = symmetric (S = 0) trial: the energy filter only resolves the singlet gap (valid because the circuit preserves J).
MPS: one isometry per site, (23/48) 4^n RZZ with n = ceil(log2 chi) + 1, eps = 1e-2.
usage: python crossover2d/model2d.py designs --L 4 6 8 10 12 16 20 --kappa 0.02 0.05 [--worker i --of K] ;  python crossover2d/model2d.py report"""
import argparse, glob, json, os, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "s2proj")); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
from mps_cost import chi_per_cut, mps_rzz                                      # noqa: E402

ED = [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "results" / "ed_*_J20.json")))]


def fits(a_t=None, c_s=None):
    E = [r for r in ED if r["W"] >= 2]
    A = np.array([[r["N"], r["W"] + r["L"], 1.0] for r in E]); y = np.array([r["E0"] for r in E])
    (a, b, c), *_ = np.linalg.lstsq(A, y, rcond=None)
    sq = [r for r in ED if r["W"] == 4]
    at = a_t or np.mean([r["N"] * r["gap_triplet"] for r in sq if r["N"] >= 20] or [8.0])
    cs = c_s or float(min(r["gap_singlet"] * r["L"] for r in sq))                    # conservative: the 4x4 value
    return dict(E0=(a, b, c), a_t=float(at), c_s=float(cs))


def chi_growth(eps=1e-2):
    best = {}
    for r in ED:
        c = chi_per_cut(r["spectra"], eps); key = r["W"]
        if key not in best or r["N"] > best[key][0]:
            best[key] = (r["N"], max(c), mps_rzz(c) / (r["N"] * 23 / 48 * 4.0 ** (np.ceil(np.log2(max(c))) + 1)))
    ws = sorted(w for w in best if w >= 2); slope, icpt = np.polyfit(ws, [np.log(best[w][1]) for w in ws], 1)
    return best, float(slope), float(icpt), float(np.mean([best[w][2] for w in ws]))


def geometry(L):
    B1 = 2 * L * (L - 1); return dict(N=L * L, B1=B1, terms=3 * B1)


def inputs(L, fit):
    g = geometry(L); a, b, c = fit["E0"]; E0 = a * g["N"] + b * 2 * L + c
    Wd = g["B1"] / 4.0 - E0
    from geometry import pauli_hamiltonian
    from core.trotter import alpha_commutator
    f = HERE / "results" / f"alpha_{L}.json"
    if not f.exists():
        H = pauli_hamiltonian(L, L, 1.0, 0.0); f.write_text(json.dumps(float(max(alpha_commutator(H, False), alpha_commutator(H, True)))))
    alpha_u = json.load(open(f))
    return dict(g, E0=E0, Wd=Wd, alpha_unscaled=alpha_u, alpha_scaled=alpha_u / Wd ** 2,
                gap_t=fit["a_t"] / g["N"] / Wd, gap_s=fit["c_s"] / L / Wd)


def get_design(gamma, delta, eps=1e-2):
    from s2_study import candidates
    import pipeline as P
    P.FLOOR_M = [6, 8, 10, 12]
    f = HERE / "results" / f"design_g{gamma:.5g}_d{delta:.4e}_e{eps:g}.json"
    if not f.exists():
        f.write_text(json.dumps(candidates(min(gamma, 1 - 1e-9), delta, 1.0, 1.0, eps)))      # alpha enters later (candidates only uses it for an unused estimate)
    return json.load(open(f))


def filter_cost(inp, gamma, mode, eps=1e-2):
    delta = inp["gap_s"] if mode == "S" else inp["gap_t"]
    cands = get_design(gamma, delta, eps)
    if not cands:
        return None
    best = None
    for c in cands:
        n = inp["alpha_scaled"] * c["T"] ** 2 / 2 / c["allowed"]; rzz = 2.0 * inp["terms"] * n
        rec = dict(rzz=rzz, rzz_per_success=rzz / c["p_g"], n=n, T=c["T"], m=c["m"], p_g=c["p_g"], delta=delta)
        if best is None or rec["rzz_per_success"] < best["rzz_per_success"]:
            best = rec
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("--L", type=int, nargs="+", default=[4, 6, 8, 10, 12, 16, 20]); ap.add_argument("--kappa", type=float, nargs="+", default=[0.02, 0.05])
    ap.add_argument("--worker", type=int, default=0); ap.add_argument("--of", type=int, default=1); a = ap.parse_args(); fit = fits()
    if a.cmd == "designs":
        jobs = [(L, k, m) for L in a.L for k in a.kappa for m in ("S", "A")]
        for i, (L, k, m) in enumerate(jobs):
            if i % a.of:
                continue
            inp = inputs(L, fit); r = filter_cost(inp, np.exp(-k * inp["N"]), m)
            print(f"L={L} kappa={k} mode={m}: N={inp['N']} Delta={inp['gap_s'] if m == 'S' else inp['gap_t']:.3e} gamma={np.exp(-k * inp['N']):.3f} -> " + (f"RZZ per success {r['rzz_per_success']:.3e} (n={r['n']:.3e}, T={r['T']:.1f}, m={r['m']})" if r else "INFEASIBLE with <= 12 pulses"), flush=True)
