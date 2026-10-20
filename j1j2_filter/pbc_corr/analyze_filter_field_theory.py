"""analyze_filter_field_theory.py -- field-theory data from FILTER-PREPARED periodic ground states, and extrapolation to sizes the filter cannot purify.

Inputs: results/pbc_N*_J2*_eps*.json (filter runs; 'error' records = not purifiable) and results/field_theory_ed.json (exact reference).
(1) For each (J2, eps): c from the half-chain entropy slope, c*v from the Casimir energy of the filtered state, eta_eff from N|C(N/2)|, compared with exact.
(2) Extrapolation test: observables from purifiable sizes only, fitted per family (N = 4n and 4n+2) with the field-theory forms
        S_half = a ln N + b ;  [N C(N/2)]^2 = a ln N + b ;  C(1), E0/N = a + b/N^2
    and predicted for the larger, non-purifiable sizes; compared with exact values there.
usage: python pbc_corr/analyze_filter_field_theory.py"""
import json, glob
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
ED = {float(k): {r["N"]: r for r in v} for k, v in json.load(open(HERE / "results" / "field_theory_ed.json")).items()}
F = {}
for f in glob.glob(str(HERE / "results" / "pbc_*.json")):
    r = json.load(open(f)); F[(r["J2"], r["eps"], r["N"])] = r
NS = (8, 10, 12, 14, 16, 18, 20)


def obs_exact(J2, N):
    e = ED[J2][N]; C = e["C"]
    return dict(S=e["S_half"], X=(-1) ** (N // 2) * N * C[-1], C1=C[0], eN=e["E0"] / N)


def obs_filter(r):
    N = r["N"]; C = r["C_filter"]
    return dict(S=r["S_filter"][-1], X=(-1) ** (N // 2) * N * C[-1], C1=C[0], eN=r["E_filter"] / N)


def fit_predict(Ns, ys, Nt, form):
    Ns, ys = np.array(Ns, float), np.array(ys, float)
    if form == "log":                       # y = a ln N + b
        A = lambda n: np.vstack([np.log(n), np.ones_like(n)]).T
        t = lambda y: y; inv = lambda y: y
    elif form == "log2":                    # y^2 = a ln N + b   (multiplicative sqrt(ln N))
        A = lambda n: np.vstack([np.log(n), np.ones_like(n)]).T
        ys = ys ** 2; inv = lambda y: np.sqrt(max(y, 0))
    else:                                   # y = a + b/N^2
        A = lambda n: np.vstack([np.ones_like(n), 1 / n ** 2]).T
        inv = lambda y: y
    cof = np.linalg.lstsq(A(Ns), ys, rcond=None)[0]
    return [inv(float((A(np.array([n], float)) @ cof)[0])) for n in Nt]


FORM = dict(S="log", X="log2", C1="inv2", eN="inv2")
report = {}
for J2 in (0.0, 0.2411, 0.4):
    for eps in (1e-3, 1e-2):
        ok = {N: F[(J2, eps, N)] for N in NS if (J2, eps, N) in F and "error" not in F[(J2, eps, N)]}
        if not ok:
            continue
        pur = sorted(ok); bad = [N for N in NS if N not in ok]
        rec = dict(purifiable=pur, not_purifiable=bad, direct={}, extrap={})
        # (1) direct comparison filter vs exact on purifiable sizes
        for key in ("S", "X", "C1", "eN"):
            rec["direct"][key] = {N: (obs_filter(ok[N])[key], obs_exact(J2, N)[key]) for N in pur}
        # (2) extrapolation test: per family, using purifiable sizes only, predict the non-purifiable ones
        for key in ("S", "X", "C1", "eN"):
            res = {}
            for fam in (0, 2):
                fn = [N for N in pur if N % 4 == fam]; tg = [N for N in bad if N % 4 == fam]
                if len(fn) >= 2 and tg:
                    pred = fit_predict(fn, [obs_filter(ok[N])[key] for N in fn], tg, FORM[key] if key != "eN" else "inv2")
                    for N, p in zip(tg, pred):
                        ex = obs_exact(J2, N)[key]; res[N] = (p, ex)
            rec["extrap"][key] = res
        report[(J2, eps)] = rec
for (J2, eps), rec in report.items():
    print(f"\n=== J2 = {J2:g}, eps = {eps:g}: purifiable N = {rec['purifiable']}, not purifiable N = {rec['not_purifiable']}")
    for key, lab in (("S", "S_half"), ("X", "(-1)^(N/2) N C(N/2)"), ("C1", "C(1)"), ("eN", "E0/N")):
        d = rec["direct"][key]
        print(f"  {lab:22s} filter vs exact, rel. err.: " + ", ".join(f"N={N}: {abs(a - b) / abs(b):.1e}" for N, (a, b) in d.items()))
        if rec["extrap"][key]:
            print(f"  {'':22s} extrapolated vs exact      : " + ", ".join(f"N={N}: {p:.5f} vs {e:.5f} ({abs(p - e) / abs(e):.1e})" for N, (p, e) in rec["extrap"][key].items()))
# entropy-based c and Casimir c*v from the filter states vs exact, on all purifiable sizes
print("\n--- central charge c from entropy slope (N, N+4) and c*v from Casimir fit (N>=10) ---")
for J2 in (0.0, 0.2411):
    for eps in (1e-2, 1e-3):
        ok = {N: F[(J2, eps, N)] for N in NS if (J2, eps, N) in F and "error" not in F[(J2, eps, N)]}
        ce = []
        for a in ok:
            if a + 4 in ok:
                Sf = (ok[a + 4]["S_filter"][-1] - ok[a]["S_filter"][-1]); Se = (ED[J2][a + 4]["S_half"] - ED[J2][a]["S_half"])
                ce.append((a, a + 4, 3 * Sf / np.log((a + 4) / a), 3 * Se / np.log((a + 4) / a)))
        line = f"J2={J2:g} eps={eps:g}: c(entropy) filter/exact: " + ", ".join(f"({a},{b}) {x:.3f}/{y:.3f}" for a, b, x, y in ce)
        Ns = [N for N in ok if N >= 10]
        if len(Ns) >= 4:
            A = np.vstack([Ns, 1 / np.array(Ns, float), 1 / np.array(Ns, float) ** 3]).T
            cf = np.linalg.lstsq(A, [ok[N]["E_filter"] for N in Ns], rcond=None)[0]; ee = np.linalg.lstsq(A, [ED[J2][N]["E0"] for N in Ns], rcond=None)[0]
            line += f" | c*v: filter {-6 * cf[1] / np.pi:.4f} exact {-6 * ee[1] / np.pi:.4f}"
        print(line)
json.dump({f"{k[0]:g}_{k[1]:g}": {"purifiable": v["purifiable"], "not_purifiable": v["not_purifiable"]} for k, v in report.items()}, open(HERE / "results" / "filter_reach.json", "w"))
# figure: extrapolation test for the three cases at the tightest eps with a non-purifiable range
fig, ax = plt.subplots(2, 3, figsize=(14, 7.5))
cases = [(0.0, 1e-3), (0.2411, 1e-3), (0.4, 1e-2)]
for c_, (J2, eps) in enumerate(cases):
    rec = report[(J2, eps)]
    for r_, (key, lab) in enumerate((("S", r"$S_{N/2}$"), ("X", r"$(-1)^{N/2}N\,C^{zz}(N/2)$"))):
        a = ax[r_, c_]
        a.plot(NS, [obs_exact(J2, N)[key] for N in NS], "k-", lw=0.8, label="exact (ED)")
        pn = rec["purifiable"]; a.plot(pn, [obs_filter(F[(J2, eps, N)])[key] for N in pn], "o", color="C0", label=f"filter-prepared (eps={eps:g})")
        for N, (p, e) in rec["extrap"][key].items():
            a.plot([N], [p], "s", color="C3", label="extrapolated" if N == min(rec["extrap"][key]) else None)
        a.set_title(f"J2={J2:g}, purifiable up to N={max(pn)}", fontsize=9); a.set_xlabel("N"); a.set_ylabel(lab); a.grid(alpha=.3); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(HERE.parent / "results" / "figures" / "filter_field_theory.png", dpi=150)
