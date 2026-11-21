"""Median CNOT-per-success ratio (improved option / authors' v4) vs overlap, per route, protocol and eps."""
import sys, os
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

d, out = sys.argv[1], sys.argv[2]
w = pd.read_csv(os.path.join(d, "design_compare_wide.csv"))
w["r"] = w["minimax any"] / w["QPF (t_i) v4"]
C = {"J2 -> H": "#2a78d6", "H only": "#eb6834"}
fig, axs = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
for ax, proto in zip(axs, ("rus", "aa")):
    for route, col in C.items():
        for eps, ls in ((1e-2, "-"), (1e-3, "--")):
            s = w[(w.route == route) & (w.protocol == proto) & (w.eps == eps)].groupby("gamma_target")["r"]
            m, lo, hi = s.median(), s.quantile(0.1), s.quantile(0.9)
            ax.plot(m.index, m.values, ls, color=col, lw=2, marker="o", ms=5, label=f"{route}, ε = {eps:g}")
            if eps == 1e-2:
                ax.fill_between(m.index, lo.values, hi.values, color=col, alpha=0.12, lw=0)
    ax.axhline(1, color="#52514e", lw=1)
    ax.set_xscale("log"); ax.set_xlabel("overlap γ"); ax.set_title(proto.upper(), fontsize=11); ax.grid(alpha=0.25)
axs[0].set_ylabel("CNOT per success: minimax option / authors' v4")
fig.legend(*axs[0].get_legend_handles_labels(), fontsize=8, frameon=False, loc="lower center", ncol=4)
fig.suptitle("Median over N = 8–14 and the 6 singlet gaps (band: 10–90%, ε = 10⁻²)", fontsize=10)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(out, dpi=150)
