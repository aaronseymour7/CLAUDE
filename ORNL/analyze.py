"""Fits and figures from results/map_N*.json.

    python analyze.py                 # reads results/, writes figs/ and results/fits.csv

Cost kinds:  attempt  CNOTs of one attempt (the filter chosen for repeat-until-success)
             rus      expected CNOTs per success, repeat until success
             aa       expected CNOTs per success, amplitude amplification
Fits:        ln C = c + a ln(1/gap) + b ln(1/gamma)              per N
             ln C = c + alpha ln N + a ln(1/gap) + b ln(1/gamma)  all N together
over gamma <= GAMMA_FIT_MAX and 1 - gamma > 10 eps (near gamma = 1 the filter is barely needed).
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

KINDS = ("attempt", "rus", "aa")
LABEL = dict(attempt="CNOT per attempt", rus="CNOT per success (RUS)", aa="CNOT per success (AA)")
GAMMA_FIT_MAX = 0.95


def load(res_dir):
    maps = {}
    for f in sorted(glob.glob(os.path.join(res_dir, "map_N*.json"))):
        d = json.load(open(f))
        maps[d["N"]] = d
    return dict(sorted(maps.items()))


def grid(d, eps, kind):
    """Arrays (gap, gamma, cost) of shape (n_gap, n_gamma); NaN where missing."""
    rows = [d["rows"][k] for k in sorted(d["rows"], key=float)]
    G = np.array([[r["gap"]] * len(d["gammas"]) for r in rows])
    Y = np.array([[np.nan if g is None else g for g in r["gamma"]] for r in rows])
    C = np.full(G.shape, np.nan)
    for i, r in enumerate(rows):
        for j, b in enumerate(r["best"]):
            x = b[f"{eps:g}"]["rus" if kind == "attempt" else kind]
            if x is not None:
                C[i, j] = x["cnot"] if kind == "attempt" else x["cost"]
    return G, Y, C


def fit(X_cols, y):
    X = np.column_stack([np.ones(len(y))] + X_cols)
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ coef
    se = np.sqrt(np.diag(res @ res / (len(y) - X.shape[1]) * np.linalg.inv(X.T @ X)))
    return coef, se, 1 - res @ res / np.sum((y - y.mean()) ** 2)


def fit_mask(Y, C, eps):
    return np.isfinite(C) & np.isfinite(Y) & (Y <= GAMMA_FIT_MAX + 1e-9) & (1 - Y > 10 * eps)


def main(res_dir="results", fig_dir="figs"):
    os.makedirs(fig_dir, exist_ok=True)
    maps = load(res_dir)
    Ns = list(maps)
    eps_list = maps[Ns[0]]["eps"]
    lines = ["eps,kind,N,a_gap,se_a,b_gamma,se_b,alpha_N,se_alpha,R2,points"]
    fits = {}
    for eps in eps_list:
        for kind in KINDS:
            allX, ally = [], []
            for N in Ns:
                G, Y, C = grid(maps[N], eps, kind)
                m = fit_mask(Y, C, eps)
                coef, se, r2 = fit([-np.log(G[m]), -np.log(Y[m])], np.log(C[m]))
                fits[eps, kind, N] = (coef[1], se[1], coef[2], se[2])
                lines.append(f"{eps:g},{kind},{N},{coef[1]:.3f},{se[1]:.3f},{coef[2]:.3f},{se[2]:.3f},,,{r2:.3f},{m.sum()}")
                allX.append(np.column_stack([np.log(np.full(m.sum(), N)), -np.log(G[m]), -np.log(Y[m])]))
                ally.append(np.log(C[m]))
            if len(Ns) > 1:
                X = np.vstack(allX)
                coef, se, r2 = fit([X[:, 0], X[:, 1], X[:, 2]], np.concatenate(ally))
                lines.append(f"{eps:g},{kind},all,{coef[2]:.3f},{se[2]:.3f},{coef[3]:.3f},{se[3]:.3f},"
                             f"{coef[1]:.3f},{se[1]:.3f},{r2:.3f},{len(X)}")
    open(os.path.join(res_dir, "fits.csv"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    colors = plt.cm.viridis(np.linspace(0, 0.9, len(Ns)))

    # 1. the two knobs: gap vs dimerisation, overlap vs eps
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for N, c in zip(Ns, colors):
        rows = [maps[N]["rows"][k] for k in sorted(maps[N]["rows"], key=float)]
        ax[0].semilogy([r["delta"] for r in rows], [r["gap"] for r in rows], "o-", color=c, ms=4, label=f"N={N}")
        mid = rows[len(rows) // 2]
        e = [x for x in mid["eps_knob"] if x is not None]
        g = [y for y in mid["gamma"] if y is not None]
        ax[1].semilogy(e, g, "o-", color=c, ms=4, label=f"N={N}, gap {mid['gap']:.2f}")
    ax[0].set(xlabel=r"dimerisation $\delta$", ylabel=r"gap $E_1-E_0$ ($J_1$)", title="gap knob (grid points)")
    ax[1].set(xlabel=r"staggered $R_z$ angle $\epsilon$", ylabel=r"overlap $\gamma$", title="overlap knob (grid points)")
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "fig_knobs.pdf")); plt.close(fig)

    for eps in eps_list:
        # 2. cost maps
        fig, ax = plt.subplots(len(KINDS), len(Ns), figsize=(3.2 * len(Ns), 2.8 * len(KINDS)), squeeze=False)
        allC = np.concatenate([grid(maps[N], eps, k)[2].ravel() for N in Ns for k in KINDS])
        vmin, vmax = np.log10(np.nanmin(allC)), np.log10(np.nanmax(allC))
        for i, kind in enumerate(KINDS):
            for j, N in enumerate(Ns):
                G, Y, C = grid(maps[N], eps, kind)
                m = ax[i, j].pcolormesh(G[:, 0], np.array(maps[N]["gammas"]), np.log10(C).T, shading="nearest",
                                        vmin=vmin, vmax=vmax, cmap="magma")
                ax[i, j].set_xscale("log"); ax[i, j].set_yscale("log")
                ax[i, j].set_title(f"N={N}: {LABEL[kind]}", fontsize=8)
                ax[i, j].set_xlabel("gap ($J_1$)"); ax[i, j].set_ylabel(r"$\gamma$")
        fig.colorbar(m, ax=ax, label="log10 CNOT")
        fig.suptitle(rf"measured cost, $1-F\leq{eps:g}$")
        fig.savefig(os.path.join(fig_dir, f"fig_maps_eps{eps:g}.pdf")); plt.close(fig)

        # 3. cuts: cost vs gap at fixed gamma, cost vs gamma at fixed gap
        gam_cut, gap_cut = (0.01, 0.1, 0.5, 0.9), (0.01, 0.1, 1.0)
        fig, ax = plt.subplots(2, len(gam_cut), figsize=(3.6 * len(gam_cut), 6.4), squeeze=False)
        for N, c in zip(Ns, colors):
            for kind, ls in (("rus", "-"), ("aa", "--"), ("attempt", ":")):
                G, Y, C = grid(maps[N], eps, kind)
                gams = np.array(maps[N]["gammas"])
                for j, gc in enumerate(gam_cut):
                    jj = np.argmin(abs(gams - gc))
                    ax[0, j].loglog(G[:, jj], C[:, jj], ls, color=c, label=f"N={N} {kind}" if j == 0 else None)
                    ax[0, j].set_title(rf"$\gamma$ = {gams[jj]:g}", fontsize=9); ax[0, j].set_xlabel("gap ($J_1$)")
                for j, gc in enumerate(gap_cut):
                    ii = np.argmin(abs(G[:, 0] - gc))
                    ax[1, j].loglog(Y[ii], C[ii], ls, color=c)
                    ax[1, j].set_title(f"gap = {G[ii, 0]:.3g}", fontsize=9); ax[1, j].set_xlabel(r"$\gamma$")
        ax[1, -1].axis("off")
        ax[0, 0].legend(fontsize=6, ncol=2)
        for a in ax.ravel():
            a.grid(alpha=0.3, which="both")
        ax[0, 0].set_ylabel("CNOT"); ax[1, 0].set_ylabel("CNOT")
        fig.suptitle(rf"$1-F\leq{eps:g}$: RUS solid, AA dashed, per attempt dotted")
        fig.tight_layout(); fig.savefig(os.path.join(fig_dir, f"fig_cuts_eps{eps:g}.pdf")); plt.close(fig)

    # 4. exponents vs N
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    for eps, mk in zip(eps_list, ("o", "s", "^")):
        for kind, c in zip(KINDS, ("C0", "C1", "C2")):
            a = np.array([fits[eps, kind, N] for N in Ns])
            ax[0].errorbar(Ns, a[:, 0], a[:, 1], marker=mk, color=c, capsize=3, label=f"{kind}, eps={eps:g}")
            ax[1].errorbar(Ns, a[:, 2], a[:, 3], marker=mk, color=c, capsize=3, label=f"{kind}, eps={eps:g}")
    ax[0].set(xlabel="N", ylabel="a  (cost ~ gap^-a)", title="gap exponent")
    ax[1].set(xlabel="N", ylabel=r"b  (cost ~ $\gamma^{-b}$)", title="overlap exponent")
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "fig_exponents.pdf")); plt.close(fig)

    # 5. local exponents: the same two-variable fit on a sliding window along one axis (all points of the other)
    fig, ax = plt.subplots(len(eps_list), 2, figsize=(10, 3.6 * len(eps_list)), squeeze=False)
    for row, eps in enumerate(eps_list):
        for N, c in zip(Ns, colors):
            for kind, ls in (("rus", "-"), ("aa", "--"), ("attempt", ":")):
                G, Y, C = grid(maps[N], eps, kind)
                m = fit_mask(Y, C, eps)
                xs, a_loc = [], []
                for i in range(G.shape[0] - 3):                    # 4 gap values per window
                    w = m.copy(); w[:i] = False; w[i + 4:] = False
                    coef = fit([-np.log(G[w]), -np.log(Y[w])], np.log(C[w]))[0]
                    xs.append(np.exp(np.log(G[i:i + 4, 0]).mean())); a_loc.append(coef[1])
                ax[row, 0].semilogx(xs, a_loc, ls, color=c, label=f"N={N} {kind}")
                gams = np.array(maps[N]["gammas"])
                xs, b_loc = [], []
                for j in range(len(gams) - 2):                     # 3 overlap values per window
                    w = m.copy(); w[:, :j] = False; w[:, j + 3:] = False
                    if len(np.unique(gams[np.any(w, axis=0)])) < 2:
                        continue
                    coef = fit([-np.log(G[w]), -np.log(Y[w])], np.log(C[w]))[0]
                    xs.append(np.exp(np.log(gams[j:j + 3]).mean())); b_loc.append(coef[2])
                ax[row, 1].semilogx(xs, b_loc, ls, color=c)
        ax[row, 0].set(xlabel="gap ($J_1$)", ylabel="local gap exponent a", title=rf"$1-F\leq{eps:g}$")
        ax[row, 1].set(xlabel=r"$\gamma$", ylabel="local overlap exponent b", title=rf"$1-F\leq{eps:g}$")
        for a in ax[row]:
            a.grid(alpha=0.3, which="both")
    ax[0, 0].legend(fontsize=6, ncol=3)
    fig.suptitle("RUS solid, AA dashed, per attempt dotted")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "fig_local_exponents.pdf")); plt.close(fig)
    print(f"figures in {fig_dir}/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--figs", default="figs")
    a = ap.parse_args()
    main(a.results, a.figs)
