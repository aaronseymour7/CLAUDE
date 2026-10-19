"""analyze_field_theory.py -- field-theory quantities from results/field_theory_ed.json (exact periodic chains N=8..20) -> table + figure."""
import json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
D = json.load(open(HERE / "results" / "field_theory_ed.json"))
fig, ax = plt.subplots(2, 3, figsize=(14, 7.5)); ax = ax.ravel()
rep = {}
for J2s, rows in D.items():
    J2 = float(J2s); N = np.array([r["N"] for r in rows]); E0 = np.array([r["E0"] for r in rows]); gt = np.array([r["gap_t"] for r in rows]); gs = np.array([r["gap_s"] for r in rows])
    S = np.array([r["S_half"] for r in rows]); Cm = np.array([abs(r["C"][-1]) for r in rows])
    v = N * gt / np.pi                                                      # x_triplet = 1/2
    # Casimir: E0 = N e_inf - (pi c v / 6)/N + b/N^3 using N >= 10;  c = (c v)/v with v from the largest-N gap
    A = np.vstack([N, 1 / N, 1 / N ** 3]).T[2:]; cof = np.linalg.lstsq(A, E0[2:], rcond=None)[0]; cv = -6 * cof[1] / np.pi
    c_casimir = cv / v[-1]
    # entanglement: slope between sizes of the same family (N, N+4)
    c_ent = []
    for a in range(len(N)):
        for b in range(a + 1, len(N)):
            if (N[b] - N[a]) == 4:
                c_ent.append((int(N[a]), int(N[b]), 3 * (S[b] - S[a]) / np.log(N[b] / N[a])))
    eta = []
    for a in range(len(N)):
        for b in range(a + 1, len(N)):
            if (N[b] - N[a]) == 4:
                eta.append((int(N[a]), int(N[b]), 1 - np.log((N[b] * Cm[b]) / (N[a] * Cm[a])) / np.log(N[b] / N[a])))
    rep[J2] = dict(v=v.tolist(), c_casimir=float(c_casimir), c_ent=c_ent, eta=eta, ratio=(gs / gt).tolist(), e_inf=float(cof[0]), cv=float(cv))
    ax[0].plot(N, v, "o-", label=f"J2={J2:g}"); ax[1].plot(N, gs / gt, "o-", label=f"J2={J2:g}")
    ax[2].plot(N, S, "o-", label=f"J2={J2:g}"); ax[3].plot(N, N * Cm, "o-", label=f"J2={J2:g}")
    big = rows[-1]; ls = np.arange(1, len(big["S_l"]) + 1); Sl = np.array([big["S_l"][str(l)] for l in ls]); Nn = big["N"]
    ax[4].plot(np.log(Nn / np.pi * np.sin(np.pi * ls / Nn)), Sl, "o-", label=f"J2={J2:g} (N={Nn})")
    ax[5].plot(N, [-(E0[i] / N[i]) for i in range(len(N))], "o-", label=f"J2={J2:g}")
ax[0].set_ylabel(r"$v_N=N\,\mathrm{gap}_{S=1}/\pi$"); ax[0].axhline(np.pi / 2, color="gray", ls=":", label=r"$\pi/2$ (J2=0 exact)")
ax[1].set_ylabel("singlet gap / triplet gap"); ax[1].axhline(1, color="gray", ls=":")
ax[2].set_ylabel(r"$S(N/2)$ (half-chain entropy)"); ax[3].set_ylabel(r"$N\,|C^{zz}(N/2)|$")
ax[4].set_xlabel(r"$\ln[(N/\pi)\sin(\pi\ell/N)]$"); ax[4].set_ylabel(r"$S(\ell)$, N=20 ($c/3$ = slope)"); ax[5].set_ylabel("$-E_0/N$")
for i, a in enumerate(ax):
    if i != 4: a.set_xlabel("N")
    a.grid(alpha=.3); a.legend(fontsize=7)
fig.tight_layout(); fig.savefig(HERE.parent / "results" / "figures" / "field_theory_ed.png", dpi=150)
for J2, r in rep.items():
    print(f"\nJ2 = {J2:g}")
    print("  v_N from triplet gap (N=8..20):", np.round(r["v"], 4), "  [pi/2 = 1.5708]")
    print("  singlet/triplet gap ratio     :", np.round(r["ratio"], 3))
    print("  c from Casimir energy (fit N>=12, v from N=20 gap): %.3f   (c*v = %.4f, e_inf = %.5f)" % (r["c_casimir"], r["cv"], r["e_inf"]))
    print("  c from entropy slope (N,N+4)  :", [(a, b, round(c, 3)) for a, b, c in r["c_ent"]])
    print("  eta_eff from N|C(N/2)| slope  :", [(a, b, round(e, 3)) for a, b, e in r["eta"]])
