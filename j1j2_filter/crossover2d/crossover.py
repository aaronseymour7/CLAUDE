"""crossover.py -- filter vs sequential-MPS state preparation, open L x L J1-J2 Heisenberg, parametrised by J2.  Documentation of every input in README.md.
Filter cost model, calibrated against the pipeline's certified designs (within 5%, gamma 0.3-0.98, Delta 0.006-0.3, eps = 1e-2):
    n = alpha_scaled T^2 / (2 allowed),  T = x pi / Delta (x ~ 1.1),  allowed = 0.056 sqrt(p_g),  p_g = 0.78 gamma,  RZZ = 2 terms n,  cost per success = RZZ / p_g.
Modes: A = no symmetry (resolve min(triplet, singlet) gap);  S = exactly-S=0 trial (resolve the singlet gap only).
Trial: hybrid = truncated MPS giving gamma = 0.3 (bond dimension fitted from W = 2..4 clusters).   MPS competitor: sequential isometries, fidelity 0.99.
usage: python crossover2d/crossover.py --j2 0.0|0.5 [--singlet inv_L|const]"""
import sys, json, glob, argparse
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent)); sys.path.insert(0, str(HERE.parent.parent))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from mps_cost import mps_rzz
from geometry import pauli_hamiltonian

X = 1.1; F0SQ = 0.78; ALLOW = 0.056


def load_ed(j2):
    return [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "results" / f"ed_*_J2{j2:g}.json")))]


def geometry(L, j2):
    B1 = 2 * L * (L - 1); B2 = 2 * (L - 1) ** 2
    return dict(N=L * L, B1=B1, B2=B2, terms=3 * (B1 + (B2 if j2 else 0)), weyl=(B1 + j2 * B2) / 4.0)


def fits_j2(j2):
    ED = load_ed(j2); E = [r for r in ED if r["W"] >= 2]
    A = np.array([[r["N"], r["W"] + r["L"], 1.0] for r in E]); y = np.array([r["E0"] for r in E]); (a, b, c), *_ = np.linalg.lstsq(A, y, rcond=None)
    sq = [r for r in ED if r["W"] == 4]
    return dict(E0=(a, b, c), a_t=float(np.mean([r["N"] * r["gap_triplet"] for r in sq])), c_s=float(min(r["gap_singlet"] * r["L"] for r in sq)),
                s_const=float(np.mean([r["gap_singlet"] for r in sq])), sq=[(r["W"], r["L"], r["gap_triplet"], r["gap_singlet"]) for r in sq])


def alpha_per_term(L, j2):
    """exact for L <= 12 (cached), else p(L) = p_inf - b/L fitted to L = 8..12."""
    def ex(l):
        f = HERE / "results" / f"alpha_{l}_J2{j2:g}.json"
        if not f.exists():
            from core.trotter import alpha_commutator
            H = pauli_hamiltonian(l, l, 1.0, j2); f.write_text(json.dumps(float(max(alpha_commutator(H, False), alpha_commutator(H, True))) / len(H)))
        return json.load(open(f))
    if L <= 12:
        return ex(L)
    p8, p12 = ex(8), ex(12); b = (p12 - p8) / (1 / 8 - 1 / 12) * -1; pinf = p12 + b / 12
    return pinf - b / L


def gaps(L, fit, singlet):
    N = L * L; gt = fit["a_t"] / N; gs = fit["c_s"] / L if singlet == "inv_L" else fit["s_const"]
    return gt, gs


def filter_cost(L, mode, gamma, j2, fit, singlet="inv_L", x=X, a_t=None):
    g = geometry(L, j2); a, b, c = fit["E0"]; E0 = a * g["N"] + b * 2 * L + c; Wd = g["weyl"] - E0
    gt, gs = gaps(L, fit, singlet)
    if a_t is not None:
        gt = a_t / g["N"]
    gap_J = gs if mode == "S" else min(gt, gs); gap = gap_J / Wd
    alpha_scaled = alpha_per_term(L, j2) * g["terms"] / Wd ** 2
    T = x * np.pi / gap; pg = F0SQ * gamma; allowed = ALLOW * np.sqrt(pg)
    n = alpha_scaled * T ** 2 / 2 / allowed; rzz = 2.0 * g["terms"] * n
    return dict(rzz=rzz, per_success=rzz / pg, n=n, T=T, gap_J=gap_J, Wd=Wd)


def fit_chi(j2, target):
    best = {}
    for r in load_ed(j2):
        sp = r["spectra"]
        for chi in range(1, 4097):
            g = 1.0
            for s in sp:
                s = np.sort(np.array(s))[::-1]; g *= max(0.0, 1 - s[chi:].sum())
            if g >= target:
                break
        if r["W"] not in best or r["N"] > best[r["W"]][0]:
            best[r["W"]] = (r["N"], chi)
    ws = [w for w in best if w >= 2]; return np.polyfit(ws, [np.log(best[w][1]) for w in ws], 1), best


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--j2", type=float, default=0.0); ap.add_argument("--singlet", default="inv_L"); a = ap.parse_args(); j2 = a.j2
    fit = fits_j2(j2); (S99, I99), b99 = fit_chi(j2, 0.99); (S30, I30), b30 = fit_chi(j2, 0.30)
    cell = lambda W, sl, ic: (23 / 48) * 4.0 ** (np.ceil(np.log2(max(2.0, np.exp(sl * W + ic)))) + 1) * 0.59
    mps_pure = lambda L: L * L * cell(L, S99, I99); mps_trial = lambda L: L * L * cell(L, S30, I30)
    print(f"J2={j2:g}: W=4 clusters (W,L,triplet gap,singlet gap): {[(w, l, round(t, 3), round(s, 3)) for w, l, t, s in fit['sq']]}")
    print(f"   a_t = N*Delta_t = {fit['a_t']:.2f};  singlet model '{a.singlet}': c_s = {fit['c_s']:.2f} (Delta_s = c_s/L) or constant {fit['s_const']:.3f};   chi(0.99) x{np.exp(S99):.2f}/width, chi(0.3) x{np.exp(S30):.2f}/width  (W=4: chi {b99[4][1]} / {b30[4][1]})")
    Ls = [4, 5, 6, 8, 10, 12, 14, 16, 20, 24]; rows = []
    print(f"\n{'L':<3s} {'N':<4s} {'chi99':>8s} {'MPS(0.99)':>10s} | {'hybrid A':>10s} {'ratio':>7s}  {'hybrid S':>10s} {'ratio':>7s} | Delta_t, Delta_s (J units)")
    for L in Ls:
        N = L * L; fa = filter_cost(L, "A", 0.3, j2, fit, a.singlet); fs = filter_cost(L, "S", 0.3, j2, fit, a.singlet); gt, gs = gaps(L, fit, a.singlet)
        rows.append((L, N, mps_pure(L), mps_trial(L) + fa["per_success"], mps_trial(L) + fs["per_success"]))
        print(f"{L:<3d} {N:<4d} {np.exp(S99 * L + I99):8.0f} {mps_pure(L):10.2e} | {rows[-1][3]:10.2e} {rows[-1][3] / rows[-1][2]:7.2f}  {rows[-1][4]:10.2e} {rows[-1][4] / rows[-1][2]:7.2f} | {gt:.3f}, {gs:.3f}")
    cross = lambda c: next((r[0] for r in rows if r[c] < r[2]), None)
    print(f"\nfirst L (on this grid) where the hybrid is cheaper than the pure MPS circuit:  A: {cross(3)}   S: {cross(4)}")
    for lab, kw in (("central", {}), ("a_t x0.6", dict(a_t=0.6 * fit["a_t"])), ("a_t x1.7", dict(a_t=1.7 * fit["a_t"])), ("x=1.5", dict(x=1.5))):
        r = [next((L for L in range(4, 61) if mps_trial(L) + filter_cost(L, m, 0.3, j2, fit, a.singlet, **kw)["per_success"] < mps_pure(L)), None) for m in "AS"]
        print(f"   sensitivity {lab:10s}: A -> L={r[0]},  S -> L={r[1]}")
    fig, ax = plt.subplots(figsize=(7.2, 4.6)); Nr = [r[1] for r in rows]
    ax.semilogy(Nr, [r[2] for r in rows], "k-o", label="pure MPS circuit (fidelity 0.99)"); ax.semilogy(Nr, [r[3] for r in rows], "C0--s", label="MPS trial (gamma 0.3) + filter, no symmetry")
    ax.semilogy(Nr, [r[4] for r in rows], "C0-s", label="MPS trial + filter, symmetric trial (singlet gap)")
    ax.set_title(f"J2 = {j2:g}, singlet gap model: {a.singlet}", fontsize=9); ax.set_xlabel("N = L x L"); ax.set_ylabel("RZZ gates (per success for the filter)"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(HERE.parent / "results" / "figures" / f"crossover2d_J2{j2:g}_{a.singlet}.png", dpi=150)
