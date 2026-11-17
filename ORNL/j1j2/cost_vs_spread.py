"""Resource cost vs spread of the infidelity, for the trials of the gap-overlap study (j1j2/infidelity_structure.py).

    python j1j2/cost_vs_spread.py            # N = 8..14 -> results/j1j2_cost_vs_spread.{json,md}, figs/j1j2_cost_vs_spread.png

For every trial (mps1, staggered-Rz knob, singlet-bond knob, at several overlaps gamma) and every chain (N, delta) the
QPF pipeline is run exactly in the common eigenbasis of H and S^2 (every step is diagonal there):

  route A  S^2 first: the J^2 ladder (pi/4, pi/8, ...) to J = 0, then the energy filter designed for the singlet gap
  route B  H only:    the energy filter designed for the lowest gap (a triplet), no J^2 stage

The energy filter is one block of the authors' QPF(t) design (v4, m = M_H steps, total time T = pi / gap_scaled on
H_s = (H - E0) / (Emax - E0)), iterated until 1 - F <= eps. Cost = total evolution time (scaled units, J^2 stage
included) / success probability, i.e. expected evolution time per success. With a fixed Trotter step this is
proportional to the CNOT count (see j1j2/trotter_costs.py for the conversion).

Spread measures (of the infidelity entering the energy filter): number of (S, E) levels holding 90% of it, its
participation ratio, and the fraction of it within 2 gaps of E0.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import infidelity_structure as IS  # noqa: E402
from qpf import design, symmetry  # noqa: E402

DELTAS = (0.0, 0.2, 0.4, 0.6)
GAMMAS = (0.9, 0.8, 0.65, 0.5, 0.35, 0.2, 0.1, 0.05)
M_H = 8
EPS = (1e-2, 1e-3)
R_MAX = 60


def spread(E, w, E0, gap):
    tot = w.sum()
    if tot <= 0:
        return dict(weight=0.0, n90=0, PR=0.0, near_gap=0.0, mean_dE_over_gap=0.0)
    return dict(weight=float(tot), n90=IS.count_for(w, 0.9), PR=float(tot ** 2 / np.sum(w ** 2)),
                near_gap=float(w[E - E0 <= 2 * gap + 1e-9].sum() / tot),
                mean_dE_over_gap=float(np.sum(w * (E - E0)) / tot / gap))


def pipeline(L, Spec, route, filt, T_block, tJ):
    S, E, w = L[:, 0], L[:, 1], L[:, 2].copy()
    ground = (S == 0) & (abs(E - Spec.E0) < 1e-6)
    P, T = 1.0, 0.0
    if route == "A":
        fac = np.prod(np.cos(np.outer(S * (S + 1), tJ)) ** 2, axis=1)
        w = w * fac
        P, T = w.sum(), float(np.sum(tJ))
        w = w / P
    f2 = filt(E) ** 2
    out = {}
    infid = 1 - w[ground].sum()
    for r in range(R_MAX + 1):
        for eps in EPS:
            if str(eps) not in out and infid <= eps:
                out[str(eps)] = dict(blocks=r, T=T, success=float(P), cost=T / P)
        if len(out) == len(EPS):
            break
        q = (w * f2).sum()
        w, P, T = w * f2 / q, P * q, T + T_block
        infid = 1 - w[ground].sum()
    return out


def run(N):
    rows = []
    M = IS.staggered(N)
    tJ = np.asarray(symmetry.j2_ladder(0, N // 2))
    for d in DELTAS:
        Spec = IS.Spectrum(N, d)
        W = Spec.Emax - Spec.E0
        base = IS.mps1(Spec.g0, N)
        filt, T = {}, {}
        for route, gap in (("A", Spec.gap_singlet), ("B", Spec.gap_triplet)):
            q = design.optimize(gap / W, M_H, method="v4")
            filt[route] = IS.energy_filter(Spec, q)
            T[route] = float(np.sum(np.abs(q["times"])))
        trials = [("mps1", 0.0, base)]
        for knob, fk, hi in (("staggered Rz", lambda e: np.exp(-1j * e * M) * base, 1.5),
                             ("singlet bond", lambda e: IS.singlet_layer(N, base, e), np.pi / 2)):
            for gt in GAMMAS:
                e = IS.eps_for(fk, Spec.g0, gt, hi)
                if e is not None:
                    trials.append((knob, e, fk(e)))
        for name, e, psi in trials:
            L = Spec.levels(psi)
            S, E, w = L[:, 0], L[:, 1], L[:, 2]
            ground = (S == 0) & (abs(E - Spec.E0) < 1e-6)
            ex = ~ground
            sing = ex & (S == 0)
            rows.append(dict(N=N, delta=d, gap_singlet=Spec.gap_singlet, gap_triplet=Spec.gap_triplet, trial=name,
                             eps_knob=e, gamma=float(w[ground].sum()),
                             spread_all=spread(E[ex], w[ex], Spec.E0, Spec.gap_triplet),
                             spread_singlet=spread(E[sing], w[sing], Spec.E0, Spec.gap_singlet),
                             A=pipeline(L, Spec, "A", filt["A"], T["A"], tJ),
                             B=pipeline(L, Spec, "B", filt["B"], T["B"], tJ)))
        print(f"N={N} delta={d} done ({len(trials)} trials)", flush=True)
    return rows


def figure(rows, path, N=None, deltas=(0.2, 0.6), eps="0.001"):
    """Cost vs overlap for both knobs (equal gamma, different spread), each point labelled by its spread."""
    N = N or max(r["N"] for r in rows)
    fig, axs = plt.subplots(len(deltas), 2, figsize=(11, 4.2 * len(deltas)), squeeze=False)
    sty = {"staggered Rz": ("o-", "#c0392b"), "singlet bond": ("s-", "#1f618d")}
    for i, d in enumerate(deltas):
        rs = [r for r in rows if r["N"] == N and abs(r["delta"] - d) < 1e-9]
        for j, (route, sk, title) in enumerate((("A", "spread_singlet", "S² first, then H"),
                                                ("B", "spread_all", "H only"))):
            a = axs[i, j]
            for knob, (m, c) in sty.items():
                pts = sorted((r["gamma"], r[route][eps]["cost"], r[sk]["n90"], r[route][eps]["blocks"])
                             for r in rs if r["trial"] == knob and eps in r[route])
                if not pts:
                    continue
                g, cost, n, b = map(np.array, zip(*pts))
                a.plot(g, cost, m, color=c, ms=5, label=knob)
                for x, y, k, bb in zip(g, cost, n, b):
                    a.annotate(f"{k}|{bb}", (x, y), textcoords="offset points", xytext=(4, 4 if c == "#c0392b" else -10),
                               fontsize=6.5, color=c)
            for r in rs:
                if r["trial"] == "mps1" and eps in r[route]:
                    a.plot(r["gamma"], r[route][eps]["cost"], "k*", ms=10, label="mps1 (no knob)")
            gs = np.linspace(0.05, 1, 50)
            a.set(xscale="log", yscale="log", xlabel="overlap γ", ylabel=f"evolution time per success, 1−F ≤ {float(eps):g}",
                  title=f"N = {N}, δ = {d:g} (singlet gap {rs[0]['gap_singlet']:.2f}, triplet {rs[0]['gap_triplet']:.2f}): {title}")
            a.title.set_fontsize(9)
            a.grid(alpha=0.3, which="both")
            a.legend(fontsize=7, title="labels: # levels for 90% of infidelity | blocks", title_fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def table(rows):
    L = ["| N | δ | trial | γ | spread into H filter, route A: 1−γ after S², # levels 90%, PR, near-gap | route A cost (blocks, P) "
         "@1e-2 / 1e-3 | spread, route B: 1−γ, # levels 90%, PR, near-gap | route B cost (blocks, P) @1e-2 / 1e-3 |",
         "|--|--|--|--|--|--|--|--|"]
    for r in rows:
        a, b = r["spread_singlet"], r["spread_all"]
        c = lambda R: " / ".join(f"{R[e]['cost']:.0f} ({R[e]['blocks']}, {R[e]['success']:.2f})" if e in R else "—"
                                 for e in map(str, EPS))
        L.append(f"| {r['N']} | {r['delta']:g} | {r['trial']} | {r['gamma']:.3f} | {a['weight']:.3f}, {a['n90']}, "
                 f"{a['PR']:.1f}, {a['near_gap']:.2f} | {c(r['A'])} | {b['weight']:.3f}, {b['n90']}, {b['PR']:.1f}, "
                 f"{b['near_gap']:.2f} | {c(r['B'])} |")
    return "\n".join(L)


def main():
    if sys.argv[1:] == ["--plot"]:
        figure(json.load(open(os.path.join(HERE, "results", "j1j2_cost_vs_spread.json"))),
               os.path.join(HERE, "figs", "j1j2_cost_vs_spread.png"))
        return
    Ns = [int(x) for x in sys.argv[1:]] or [8, 10, 12, 14]
    rows = []
    for N in Ns:
        rows += run(N)
        json.dump(rows, open(os.path.join(HERE, "results", "j1j2_cost_vs_spread.json"), "w"), indent=1)
        open(os.path.join(HERE, "results", "j1j2_cost_vs_spread.md"), "w").write(table(rows) + "\n")
    figure(rows, os.path.join(HERE, "figs", "j1j2_cost_vs_spread.png"))


if __name__ == "__main__":
    main()
