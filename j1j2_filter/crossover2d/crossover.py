"""crossover.py -- filter vs sequential-MPS state preparation, open L x L Heisenberg (J2 = 0).  All inputs documented in model2d.py / README; calibration of the filter
cost model against the pipeline's certified designs:  n = alpha_scaled T^2 / (2 allowed),  T = x pi / Delta (x ~ 1.0-1.25),  allowed = 0.056 sqrt(p_g),  p_g = 0.78 gamma
(reproduces the actual certified step counts within 5% for gamma 0.3-0.98, Delta 0.006-0.3, eps = 1e-2).
usage: python crossover2d/crossover.py"""
import sys, json, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.argv = ["x"]
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import model2d as M
from mps_cost import mps_rzz

FIT = M.fits(); A_T, C_S = FIT["a_t"], FIT["c_s"]; X = 1.1; F0SQ = 0.78; ALLOW = 0.056


def alpha_unscaled(L, terms):
    f = HERE / "results" / f"alpha_{L}.json"
    return json.load(open(f)) if f.exists() else terms * (0.553 - 0.48 / L)


def geom(L):
    g = M.geometry(L); a, b, c = FIT["E0"]; E0 = a * g["N"] + b * 2 * L + c; Wd = g["B1"] / 4 - E0
    return dict(g, Wd=Wd, alpha_scaled=alpha_unscaled(L, g["terms"]) / Wd ** 2)


def filter_cost(L, mode, gamma, a_t=A_T, c_s=C_S, x=X):
    g = geom(L); gap = (c_s / L if mode == "S" else a_t / g["N"]) / g["Wd"]
    T = x * np.pi / gap; pg = F0SQ * gamma; allowed = ALLOW * np.sqrt(pg)
    n = g["alpha_scaled"] * T ** 2 / 2 / allowed; rzz = 2.0 * g["terms"] * n
    return dict(rzz=rzz, per_success=rzz / pg, n=n, T=T, gap=gap)


def fit_chi(target_key):
    ED = [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "results" / "ed_*_J20.json")))]

    def gamma_of(sp, chi):
        g = 1.0
        for s in sp:
            s = np.sort(np.array(s))[::-1]; g *= max(0.0, 1 - s[chi:].sum())
        return g
    best = {}
    for r in ED:
        for chi in range(1, 4097):
            if gamma_of(r["spectra"], chi) >= target_key:
                break
        if r["W"] not in best or r["N"] > best[r["W"]][0]:
            best[r["W"]] = (r["N"], chi)
    ws = [w for w in best if w >= 2]; return np.polyfit(ws, [np.log(best[w][1]) for w in ws], 1)


S99, I99 = fit_chi(0.99); S30, I30 = fit_chi(0.30)
mps_cell = lambda W, sl, ic: (23 / 48) * 4.0 ** (np.ceil(np.log2(max(2.0, np.exp(sl * W + ic)))) + 1) * 0.59


def mps_pure(L):
    return L * L * mps_cell(L, S99, I99)


def mps_trial(L):
    return L * L * mps_cell(L, S30, I30)


if __name__ == "__main__":
    Ls = [4, 5, 6, 8, 10, 12, 14, 16, 20, 24]
    print("L   N    chi99  MPS(0.99)    | hybrid: trial(g>=0.3) + filter (per success)         | filter only, generic ansatz gamma=exp(-kappa N)")
    print("                              | A: triplet gap      S: singlet gap (sym trial)       | A kappa=.02     A kappa=.05")
    rows = []
    for L in Ls:
        N = L * L; chi = np.exp(S99 * L + I99); fa = filter_cost(L, "A", 0.3); fs = filter_cost(L, "S", 0.3)
        k2 = filter_cost(L, "A", np.exp(-0.02 * N))["per_success"]; k5 = filter_cost(L, "A", np.exp(-0.05 * N))["per_success"]
        rows.append((L, N, chi, mps_pure(L), mps_trial(L) + fa["per_success"], mps_trial(L) + fs["per_success"], k2, k5))
        print(f"{L:<3d} {N:<4d} {chi:7.0f} {mps_pure(L):10.2e}   | {rows[-1][4]:10.2e} ({rows[-1][4] / rows[-1][3]:6.2f}x)  {rows[-1][5]:10.2e} ({rows[-1][5] / rows[-1][3]:6.2f}x) | {k2:10.2e}  {k5:10.2e}")
    cross = lambda col: next((r[0] for r in rows if r[col] < r[3]), None)
    print("\nfirst L where the filter is cheaper than the pure MPS circuit:  hybrid A:", cross(4), " hybrid S:", cross(5), " generic kappa=0.02:", cross(6), " kappa=0.05:", cross(7))
    # sensitivity of the hybrid-S crossover to the modelled gap constants and chi growth
    print("\nsensitivity of the first-crossover L (hybrid A / hybrid S):")
    for lab, kw in (("central", {}), ("a_t=6 (smaller triplet gap)", dict(a_t=6.0)), ("a_t=15", dict(a_t=15.0)), ("c_s=4", dict(c_s=4.0)), ("c_s=8", dict(c_s=8.0)), ("x=1.5", dict(x=1.5))):
        res = []
        for mode in ("A", "S"):
            c = next((L for L in range(4, 41) if mps_trial(L) + filter_cost(L, mode, 0.3, **kw)["per_success"] < mps_pure(L)), None); res.append(c)
        print(f"   {lab:32s}: A -> L={res[0]},  S -> L={res[1]}")
    fig, ax = plt.subplots(figsize=(7.2, 4.8)); Lr = np.array([r[0] for r in rows]); Nr = Lr ** 2
    ax.semilogy(Nr, [r[3] for r in rows], "k-o", label="pure MPS circuit (fidelity 0.99)")
    ax.semilogy(Nr, [r[4] for r in rows], "C0--s", label="hybrid: MPS trial (gamma 0.3) + filter, triplet gap")
    ax.semilogy(Nr, [r[5] for r in rows], "C0-s", label="hybrid + symmetric trial (singlet gap)")
    ax.semilogy(Nr, [r[6] for r in rows], ":", color="C3", label="filter, fixed ansatz, gamma = exp(-0.02 N)"); ax.semilogy(Nr, [r[7] for r in rows], ":", color="C1", label="filter, fixed ansatz, gamma = exp(-0.05 N)")
    ax.set_xlabel("number of spins N = L x L"); ax.set_ylabel("RZZ gates (per success for the filter)"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
    fig.tight_layout(); (HERE.parent / "results" / "figures").mkdir(exist_ok=True, parents=True); fig.savefig(HERE.parent / "results" / "figures" / "crossover2d.png", dpi=150)
