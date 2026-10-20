"""analyze_dephasing.py -- table + figure from trotter2/dephasing_N*_J2*.json  ->  trotter2/dephasing_table.md, trotter2/dephasing.png"""
import json, glob
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
R = sorted([json.load(open(f)) for f in glob.glob(str(HERE / "dephasing_N*_J2*.json"))], key=lambda r: (r["N"], r["J2"]))
L = ["| N | J2 | states 2^N | distinct levels | γ | levels with >1e-4 of trial weight | …after filter | levels for 99.99% (trial → filtered) | weight above gap/2 | state error² 1st | 2nd | share in suppressed levels 1st / 2nd | floquet exponent (gs) 1st / 2nd |", "|" + "--|" * 13]
for r in R:
    t1, t2 = r["trotter"]["order1"], r["trotter"]["order2"]; ct, cf = r["counts_trial"], r["counts_after_filter"]
    fl = r.get("floquet"); fe = f"{fl['order1']['exponent_gs']:.2f} / {fl['order2']['exponent_gs']:.2f}" if fl else "–"
    L.append(f"| {r['N']} | {r['J2']:g} | {r['dim']} | {r['n_levels']} | {r['gamma']:.3f} | {ct['levels>0.0001']} | {cf['levels>0.0001']} | {ct['levels_for_99.99pct']} → {cf['levels_for_99.99pct']} | {r['weight_above_gap']:.3f} | "
             f"{t1['err_norm2_total']:.1e} | {t2['err_norm2_total']:.1e} | {t1['err_groups']['suppressed'] / t1['err_norm2_total']:.2f} / {t2['err_groups']['suppressed'] / t2['err_norm2_total']:.2f} | {fe} |")
open(HERE / "dephasing_table.md", "w").write("\n".join(L) + "\n"); print("\n".join(L))

sel = [r for r in R if r["N"] == max(x["N"] for x in R if x["J2"] == 0.4) and r["J2"] == 0.4] or R[-1:]
r = sel[0]; E = np.array(r["levels"]["E"]); S = np.array(r["levels"]["S"]); p0 = np.array(r["levels"]["p_trial"]); p1 = np.array(r["levels"]["p_after_filter"])
fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
cm = {0: "#1f4e79", 1: "#c0392b", 2: "#2e8b57"}
for s in sorted(set(S)):
    m = S == s
    ax[0].semilogy(E[m], np.maximum(p0[m], 1e-16), "o", ms=4, color=cm.get(s, "#888"), label=f"S={s}")
    ax[1].semilogy(E[m], np.maximum(p1[m], 1e-16), "o", ms=4, color=cm.get(s, "#888"))
ax[0].set_title(f"trial state weight per level (N={r['N']}, J2={r['J2']:g})", fontsize=9); ax[1].set_title("after the full filter (normalised)", fontsize=9)
for a in ax[:2]:
    a.set_xlabel("scaled level energy"); a.set_ylim(1e-12, 2)
ax[0].legend(fontsize=7, frameon=False)
e1 = np.array(r["trotter"]["order1"]["err_by_level"]); e2 = np.array(r["trotter"]["order2"]["err_by_level"])
ax[2].semilogy(E, np.maximum(e1, 1e-20), "o", ms=4, color="#1f4e79", label=f"first order (n={r['trotter']['order1']['n']})")
ax[2].semilogy(E, np.maximum(e2, 1e-20), "s", ms=4, color="#c0392b", label=f"second order (n={r['trotter']['order2']['n']})")
ax[2].set_xlabel("scaled level energy"); ax[2].set_title("state error² resolved by level", fontsize=9); ax[2].legend(fontsize=7, frameon=False); ax[2].set_ylim(1e-18, 1e-3)
fig.tight_layout(); fig.savefig(HERE / "dephasing.png", dpi=170); print("wrote figure")
