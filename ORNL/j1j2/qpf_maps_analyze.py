"""Tables, fits and figures from the QPF gap-vs-overlap runs (j1j2/qpf_maps.py).

    python j1j2/qpf_maps_analyze.py RAW_DIR OUT_DIR

Writes OUT_DIR/qpf_points.csv (one row per N, delta, knob, gamma, eps, route, schedule, protocol; schedule "best" =
the cheaper of the two published designs at that point), qpf_fits.csv and the figures:
Figures use the analytic gate model (GM); the CSVs hold both (column gate_model).
  fig_qpf_maps_{rus,aa}.png      CNOT per success over (singlet gap / W, gamma), rows eps, columns N (J2 -> H, staggered)
  fig_qpf_cuts.png               cuts at fixed gamma and fixed gap, both knobs and both routes
  fig_qpf_shot_vs_success.png    per-attempt CNOT, evolution time and CNOT per success vs gamma and vs gap (A3)
  fig_qpf_runtime.png            T vs L / Delta, Trotter dt vs eps, Trotter steps vs L / Delta (A4)
  fig_qpf_exponents.png          fitted gap and overlap exponents vs N
Fit: ln C = c + a ln(W / Delta) + b ln(1 / gamma) (+ alpha ln N for "all"), Delta = the gap the energy stage
resolves (D5), W = Emax - E0.
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

LIST_COLS = ("ns_J", "ns_E", "times_E_scaled", "phases_E")


def load(raw):
    rows = []
    for f in sorted(glob.glob(os.path.join(raw, "qpf_N*_gap*.json"))):
        p = json.load(open(f))
        base = {k: p[k] for k in ("N", "delta", "gap_target", "E0", "Emax", "W", "gap_singlet", "gap_global")}
        for r in p["results"]:
            row = dict(base, **{k: v for k, v in r.items()})
            for k in LIST_COLS:
                if k in row and isinstance(row[k], list):
                    row[k] = ";".join(f"{x:.10g}" if isinstance(x, float) else str(x) for x in row[k])
            rows.append(row)
    d = pd.DataFrame(rows)
    d["gap_over_W"] = d.gap_used_value / d.W
    d["L"] = np.log((1 - d.gamma) / (d.gamma * d.eps))
    d["T_total"] = d.time_E + d.time_J.fillna(0)            # energy stage in 1/J1; J^2 stage in units of 1/(J^2 unit)
    d["cx_E_attempt"] = d.R * d.cx_E_block                  # energy stage only, per attempt
    return d


def add_best(d):
    ok = d[d.reached == True].copy()  # noqa: E712
    key = ["N", "delta", "knob", "gamma_target", "eps", "route", "gate_model", "protocol"]
    i = ok.groupby(key).cx_per_success.idxmin()
    b = ok.loc[i].copy()
    b["schedule_chosen"] = b.schedule
    b["schedule"] = "best"
    return pd.concat([d, b], ignore_index=True)


def fit(s, col, withN=False):
    s = s[np.isfinite(s[col]) & (s[col] > 0)]
    X = [np.ones(len(s)), np.log(1 / s.gap_over_W), np.log(1 / s.gamma)]
    if withN:
        X.append(np.log(s.N))
    X = np.column_stack(X)
    y = np.log(s[col].values)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ c
    dof = max(1, len(y) - X.shape[1])
    cov = np.linalg.pinv(X.T @ X) * (r @ r) / dof
    se = np.sqrt(np.diag(cov))
    R2 = 1 - (r @ r) / np.sum((y - y.mean()) ** 2)
    out = dict(a_gap=c[1], se_a=se[1], b_gamma=c[2], se_b=se[2], R2=R2, points=len(y), C0=np.exp(c[0]))
    if withN:
        out.update(alpha_N=c[3], se_alpha=se[3])
    return out


def fits(d):
    rows = []
    b = d[(d.schedule == "best")]
    for (gm, route, knob, eps), s in b.groupby(["gate_model", "route", "knob", "eps"]):
        for kind, col, proto in (("attempt", "cx_attempt", "rus"), ("rus", "cx_per_success", "rus"),
                                 ("aa", "cx_per_success", "aa"), ("energy_stage_attempt", "cx_E_attempt", "rus"),
                                 ("time_E", "time_E", "rus")):
            t = s[s.protocol == proto]
            for N, u in list(t.groupby("N")) + [("all", t)]:
                if len(u) < 5:
                    continue
                rows.append(dict(gate_model=gm, route=route, knob=knob, eps=eps, kind=kind, N=N, **fit(u, col, withN=(N == "all"))))
    return pd.DataFrame(rows)


GM = "analytic"


def log_edges(x):
    """Cell edges for log-spaced centres (geometric midpoints, half a step beyond the ends)."""
    lx = np.log(np.asarray(x, float))
    if len(lx) == 1:
        return np.exp([lx[0] - 0.1, lx[0] + 0.1])
    m = (lx[1:] + lx[:-1]) / 2
    return np.exp(np.concatenate([[2 * lx[0] - m[0]], m, [2 * lx[-1] - m[-1]]]))


def maps(d, path, proto, route="J2 -> H", knob="staggered"):
    d = d[d.gate_model == GM]
    s = d[(d.schedule == "best") & (d.protocol == proto) & (d.route == route) & (d.knob == knob)]
    Ns, epss = sorted(s.N.unique()), sorted(s.eps.unique(), reverse=True)
    from matplotlib.colors import LogNorm
    norm = LogNorm(s.cx_per_success.min(), s.cx_per_success.max())
    fig, ax = plt.subplots(len(epss), len(Ns), figsize=(3.4 * len(Ns) + 1, 3.1 * len(epss)), squeeze=False,
                           sharex=True, sharey=True, constrained_layout=True)
    pc = None
    for i, e in enumerate(epss):
        for j, N in enumerate(Ns):
            a = ax[i, j]
            t = s[(s.N == N) & (s.eps == e)].pivot_table(index="gamma_target", columns="gap_over_W", values="cx_per_success")
            if t.empty:
                continue
            pc = a.pcolormesh(log_edges(t.columns.values), log_edges(t.index.values), t.values, norm=norm,
                              cmap="magma_r", shading="flat")
            lv = [x for x in (3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6) if np.nanmin(t.values) < x < np.nanmax(t.values)]
            if lv and t.shape[0] > 1 and t.shape[1] > 1:
                cs = a.contour(t.columns, t.index, t.values, levels=lv, colors="k", linewidths=0.6)
                a.clabel(cs, fmt=lambda v: f"{v:.0e}", fontsize=6)
            a.set(xscale="log", yscale="log", title=f"N = {N}, ε = {e:g}")
            if i == len(epss) - 1:
                a.set_xlabel("singlet gap Δ / W")
            if j == 0:
                a.set_ylabel("overlap γ")
    if pc is not None:
        fig.colorbar(pc, ax=ax, label=f"CNOT per success ({proto.upper()})", shrink=0.7)
    fig.suptitle(f"QPF ({route}, {knob} knob, best of the two published designs, {GM} gate model): {proto.upper()}", fontsize=10)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def cuts(d, path, gamma_fix=0.1, eps=0.01):
    d = d[d.gate_model == GM]
    s = d[(d.schedule == "best") & (d.eps == eps)]
    Ns = sorted(s.N.unique())
    gaps = sorted(s.gap_target.unique())
    gap_fix = gaps[len(gaps) // 2]
    col = dict(zip(Ns, plt.cm.viridis(np.linspace(0, 0.85, len(Ns)))))
    combos = [("J2 -> H", "staggered", "o-"), ("J2 -> H", "singlet", "s--"), ("H only", "staggered", "^:"),
              ("H only", "singlet", "d-.")]
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    for r, proto in enumerate(("rus", "aa")):
        for route, knob, st in combos:
            for N in Ns:
                u = s[(s.protocol == proto) & (s.route == route) & (s.knob == knob) & (s.N == N)]
                v = u[np.isclose(u.gamma_target, gamma_fix)].sort_values("gap_singlet")
                ax[r, 0].loglog(v.gap_singlet, v.cx_per_success, st, color=col[N], ms=4)
                v = u[np.isclose(u.gap_target, gap_fix)].sort_values("gamma")
                ax[r, 1].loglog(v.gamma, v.cx_per_success, st, color=col[N], ms=4)
        ax[r, 0].set(xlabel="singlet gap Δ_S (J1)", ylabel=f"CNOT per success ({proto.upper()})",
                     title=f"γ = {gamma_fix}, ε = {eps:g}")
        ax[r, 1].set(xlabel="overlap γ", title=f"singlet gap = {gap_fix} J1, ε = {eps:g}")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=col[N], lw=2, label=f"N = {N}") for N in Ns] + \
        [Line2D([], [], color="0.3", ls=st[1:], marker=st[0], label=f"{route}, {knob}") for route, knob, st in combos]
    ax[0, 1].legend(handles=h, fontsize=7, ncol=2)
    for a in ax.flat:
        a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def shot_vs_success(d, path, route="J2 -> H", knob="staggered", eps=0.01):
    d = d[d.gate_model == GM]
    s = d[(d.schedule == "best") & (d.eps == eps) & (d.protocol == "rus") & (d.route == route) & (d.knob == knob)]
    Ns = sorted(s.N.unique())
    gaps = sorted(s.gap_target.unique())
    gap_fix = gaps[len(gaps) // 2]
    col = dict(zip(Ns, plt.cm.viridis(np.linspace(0, 0.85, len(Ns)))))
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    for N in Ns:
        u = s[s.N == N]
        a = u[np.isclose(u.gap_target, gap_fix)].sort_values("gamma")
        b = u[np.isclose(u.gamma_target, 0.1)].sort_values("gap_singlet")
        for r, (v, x) in enumerate(((a, "gamma"), (b, "gap_singlet"))):
            ax[r, 0].loglog(v[x], v.cx_attempt, "o-", color=col[N], ms=4, label=f"N = {N}")
            ax[r, 1].loglog(v[x], v.time_E, "o-", color=col[N], ms=4)
            ax[r, 2].loglog(v[x], v.cx_per_success, "o-", color=col[N], ms=4)
    for j, yl in enumerate(("CNOT per attempt", "energy-stage evolution time (1/J1)", "CNOT per success (RUS)")):
        ax[0, j].set(xlabel="overlap γ", ylabel=yl, title=f"singlet gap = {gap_fix} J1")
        ax[1, j].set(xlabel="singlet gap Δ_S (J1)", ylabel=yl, title="γ = 0.1")
    ax[0, 0].legend(fontsize=7)
    for a in ax.flat:
        a.grid(alpha=0.3, which="both")
    fig.suptitle(f"Per shot vs per success: QPF {route}, {knob} knob, ε = {eps:g}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def runtime(d, path, route="J2 -> H", knob="staggered"):
    d = d[d.gate_model == GM]
    s = d[(d.schedule == "best") & (d.protocol == "rus") & (d.route == route) & (d.knob == knob)]
    Ns = sorted(s.N.unique())
    col = dict(zip(Ns, plt.cm.viridis(np.linspace(0, 0.85, len(Ns)))))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for N in Ns:
        for e, m in ((0.01, "o"), (0.001, "s")):
            u = s[(s.N == N) & (s.eps == e)]
            ax[0].loglog(u.L / u.gap_used_value, u.time_E, m, color=col[N], ms=4, alpha=0.8)
            ax[2].loglog(u.L / u.gap_used_value, u.trotter_steps, m, color=col[N], ms=4, alpha=0.8)
    for N in Ns:
        u = s[s.N == N]
        ax[1].scatter(u.eps * (1 + 0.06 * (N - 12)), u.dt_E, color=col[N], s=10, alpha=0.6, label=f"N = {N}")
    ax[0].set(xlabel="L / Δ, L = ln[(1−γ)/(γε)]", ylabel="energy-stage evolution time T (1/J1)",
              title="T vs L/Δ (o ε = 1e-2, □ 1e-3)")
    ax[1].set(xscale="log", yscale="log", xlabel="ε", ylabel="Trotter step dt of the energy stage (1/J1)", title="dt vs ε")
    ax[2].set(xlabel="L / Δ", ylabel="Trotter steps (J² + energy)", title="steps vs L/Δ")
    ax[1].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def exponents(f, path):
    if f.empty:
        return
    f = f[f.gate_model == GM]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    sty = {"attempt": "o-", "rus": "s-", "aa": "d--"}
    for (route, knob, eps, kind), u in f[f.N != "all"].groupby(["route", "knob", "eps", "kind"]):
        if eps != 0.01 or route != "J2 -> H" or kind not in sty:
            continue
        u = u.sort_values("N")
        lab = f"{knob}, {kind}"
        ax[0].errorbar(u.N.astype(int), u.a_gap, u.se_a, fmt=sty[kind], ms=4, capsize=2, label=lab)
        ax[1].errorbar(u.N.astype(int), u.b_gamma, u.se_b, fmt=sty[kind], ms=4, capsize=2, label=lab)
    ax[0].set(xlabel="N", ylabel="gap exponent a", title="ln C = c + a ln(W/Δ) + b ln(1/γ), QPF J2 -> H, ε = 1e-2")
    ax[1].set(xlabel="N", ylabel="overlap exponent b")
    ax[1].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    raw, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    d = add_best(load(raw))
    d.to_csv(os.path.join(out, "qpf_points.csv"), index=False)
    f = fits(d)
    f.to_csv(os.path.join(out, "qpf_fits.csv"), index=False)
    for proto in ("rus", "aa"):
        maps(d, os.path.join(out, f"fig_qpf_maps_{proto}.png"), proto)
        maps(d, os.path.join(out, f"fig_qpf_maps_{proto}_singletknob.png"), proto, knob="singlet")
    cuts(d, os.path.join(out, "fig_qpf_cuts.png"))
    cuts(d, os.path.join(out, "fig_qpf_cuts_eps0.001.png"), eps=0.001)
    shot_vs_success(d, os.path.join(out, "fig_qpf_shot_vs_success.png"))
    runtime(d, os.path.join(out, "fig_qpf_runtime.png"))
    exponents(f, os.path.join(out, "fig_qpf_exponents.png"))
    print(f"{len(d)} rows, {d.reached.sum()} reached; N = {sorted(d.N.unique())}")


if __name__ == "__main__":
    main()
