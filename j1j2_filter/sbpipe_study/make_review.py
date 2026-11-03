"""Tables, figures and macros for the review additions to paper/paper_measured.tex.

    python sbpipe_study/make_review.py            # everything whose inputs exist

Writes paper/generated/macros_review.tex, paper/generated/tab_fitci.tex and paper/figs/fig_*_rev.pdf.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
import math
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import make_figures as M
from derive_gapov_filter import load_cases

RES, GEN, FIG = M.RES, M.GEN, M.FIG
MAC = {}
COL = {"dyadic": M.SCHED_COL["dyadic"], "optimized": M.SCHED_COL["optimized"], "qetu": M.QETU_COL}


def macro(k, v):
    MAC[k] = v


def write_macros():
    p = os.path.join(GEN, "macros_review.tex")
    old = {}
    if os.path.exists(p):
        for l in open(p):
            if l.startswith("\\newcommand{\\"):
                k = l[len("\\newcommand{\\"):l.index("}")]
                old[k] = l
    for k, v in MAC.items():
        old[k] = f"\\newcommand{{\\{k}}}{{{v}}}\n"
    with open(p, "w") as fh:
        fh.write("".join(old[k] for k in sorted(old)))


def jsonl(name):
    p = os.path.join(RES, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


# ----------------------------------------------------------------------------- 1. fits with confidence intervals
KEYS = (("meas_dyadic", "dyadic"), ("meas_optimized", "optimized"), ("meas_qetu", "QETU"))


def _cost(r, m, what):
    p = m["p"]
    fil = m["cnot"] * p - r["trial_cnot"]
    return dict(succ=m["cnot"], succ_filter=fil / p, shot_filter=fil)[what]


def _design(G, key, what):
    X, y, grp = [], [], []
    for r in G:
        m = r[key]
        if not m:
            continue
        X.append([1, np.log(r["N"]), -np.log(r["gap"]), np.log(np.log(1 / r["eps"])), -np.log(r["gamma"])])
        y.append(np.log(_cost(r, m, what)))
        grp.append((r["N"], r["gap_target"]))
    return np.array(X), np.array(y), grp


def _fit_ci(X, y, grp, nboot=2000, seed=0):
    """Least squares in logs; 95% intervals from a cluster bootstrap over the (N, gap) systems."""
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    rng = np.random.default_rng(seed)
    keys = sorted(set(grp))
    idx = {k: [i for i, g in enumerate(grp) if g == k] for k in keys}
    B = []
    for _ in range(nboot):
        pick = [i for k in rng.choice(len(keys), len(keys)) for i in idx[keys[k]]]
        b, *_ = np.linalg.lstsq(X[pick], y[pick], rcond=None)
        B.append(b)
    lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
    return c, lo, hi, float(np.exp(np.sqrt(np.mean(res ** 2)))), float(np.exp(np.abs(res).max()))


def fits():
    G = load_cases()
    rows = []
    for key, name in KEYS:
        for what, wname in (("shot_filter", "filter, per shot"), ("succ_filter", "filter, per success"),
                            ("succ", "filter + trial, per success")):
            X, y, grp = _design(G, key, what)
            c, lo, hi, rms, mx = _fit_ci(X, y, grp)
            cell = lambda i: f"${c[i]:.2f}_{{[{lo[i]:.2f},\\,{hi[i]:.2f}]}}$"
            rows.append(f"{name} & {wname} & ${M.sci(np.exp(c[0]), 1)}$ & {cell(1)} & {cell(2)} & {cell(3)} & "
                        f"{cell(4)} & {rms:.2f} \\\\")
            t = {"meas_dyadic": "D", "meas_optimized": "O", "meas_qetu": "Q"}[key]
            w = {"shot_filter": "Shot", "succ_filter": "Succ", "succ": "Tot"}[what]
            for i, nm in ((1, "N"), (2, "gap"), (3, "eps"), (4, "gam")):
                macro(f"ci{t}{w}{nm}", f"{c[i]:.2f}")
                macro(f"ci{t}{w}{nm}lo", f"{lo[i]:.2f}")
                macro(f"ci{t}{w}{nm}hi", f"{hi[i]:.2f}")
            macro(f"ci{t}{w}rms", f"{rms:.2f}")
            macro(f"ci{t}{w}a", f"{np.exp(c[0]):.3g}")
        rows.append("\\addlinespace")
    with open(os.path.join(GEN, "tab_fitci.tex"), "w") as fh:
        fh.write("\n".join(rows[:-1]) + "\n")
    macro("ciNcases", str(len(G)))
    macro("ciNsys", str(len({(r["N"], r["gap_target"]) for r in G})))
    # the per-shot ranges quoted in the text, with the exponent they imply
    A = [r for r in G if r["study"] == "A" and r["eps"] == 1e-3 and r["meas_dyadic"]]
    B = [r for r in G if r["study"] == "B" and r["eps"] == 1e-3 and r["meas_dyadic"]]
    return G


# ----------------------------------------------------------------------------- 2. measured step vs bounds
def _psd_shift_energy(N, J2, E0):
    """Energy above the zero of the positive-semidefinite bond terms h_b + 3J/4 (Mizuta-Kuwahara)."""
    return E0 + 0.75 * ((N - 1) * 1.0 + (N - 2) * J2)


def mk_dt(T, budget, N, J2, E0, p=2, cp=2):
    """Coarsest common step allowed by the explicit low-energy bound of Mizuta and Kuwahara (Corollary S4)
    for a state inside the ground space (the most favourable case), Strang (p = 2, c_p = 2 exponentials per
    group), Gamma commuting groups (2 for J2 = 0, 4 otherwise, the smallest edge colouring), k = 2,
    g = 2 J1 + 2 J2. The bound holds only for dt <= 1/(2 e c_p Gamma k g)."""
    k, g = 2, 2.0 + 2.0 * J2
    Gam = 2 if J2 == 0 else 4
    A = 2 * cp * Gam * k * g
    C = 2 * cp * Gam / ((1 - math.exp(-1)) * (p + 1))
    D = _psd_shift_energy(N, J2, E0)
    dt = 1.0 / (math.e * A)
    for _ in range(50):
        n = max(1.0, T / dt)
        aux = budget / (2 * n)
        Dp = D + 4 * k * g * math.log(2 ** (-p + 1) / k * N / aux)
        new = math.sqrt(budget / (2 * T * C * A ** p * Dp))
        new = min(new, 1.0 / (math.e * A))
        if abs(new - dt) < 1e-12 * dt:
            break
        dt = new
    return dt, 1.0 / (math.e * A)


def bounds():
    from sbpipe import j1j2_chain
    from sbpipe.trotter import beta2
    from scipy.sparse.linalg import eigsh
    G = load_cases()
    sysinfo = {}
    for r in G:
        key = (r["N"], round(r["J2"], 7))
        if key in sysinfo:
            continue
        H = j1j2_chain(r["N"], 1.0, r["J2"])
        E0 = float(eigsh(H.to_sparse(), k=1, which="SA", return_eigenvectors=False, tol=1e-10)[0])
        Hs = H.shifted_scaled(E0, r["W"])
        sysinfo[key] = dict(E0=E0, beta=beta2(Hs))
    pts = []
    for r in G:
        s = sysinfo[(r["N"], round(r["J2"], 7))]
        for key, fam in (("meas_dyadic", "dyadic"), ("meas_optimized", "optimized"), ("meas_qetu", "qetu")):
            m = r[key]
            if not m:
                continue
            if fam == "qetu":
                T, dt = m["d"] * math.pi / 2, math.pi / 2 / m["k"]
            else:
                T, dt = m["T"], m["dt"]
            budget = 0.5 * math.sqrt(r["eps"]) * math.sqrt(m["p"])
            dt_c = math.sqrt(budget / (s["beta"] * T)) / r["W"]          # physical units
            dt_mk, dt_valid = mk_dt(T / r["W"], budget, r["N"], r["J2"], s["E0"])
            pts.append(dict(fam=fam, eps=r["eps"], N=r["N"], gap=r["gap"], dt=dt / r["W"], dt_comm=dt_c,
                            dt_mk=dt_mk, dt_valid=dt_valid))
    json.dump(pts, open(os.path.join(RES, "review_bounds.json"), "w"))
    D = [p for p in pts if p["fam"] == "dyadic"]
    rc = np.array([p["dt"] / p["dt_comm"] for p in D])
    rm = np.array([p["dt"] / p["dt_mk"] for p in D])
    macro("bdCommMed", f"{np.median(rc):.0f}"); macro("bdCommLo", f"{rc.min():.0f}"); macro("bdCommHi", f"{rc.max():.0f}")
    macro("bdMKMed", f"{np.median(rm):.0f}"); macro("bdMKLo", f"{rm.min():.0f}"); macro("bdMKHi", f"{rm.max():.0f}")
    macro("bdMKvalid", f"{min(p['dt_valid'] for p in D):.4f}")
    macro("bdMKvalidHi", f"{max(p['dt_valid'] for p in D):.4f}")
    for fam in ("optimized", "qetu"):
        F = [p for p in pts if p["fam"] == fam]
        macro(f"bdComm{fam.capitalize()}Med", f"{np.median([p['dt'] / p['dt_comm'] for p in F]):.0f}")
    # figure: step against eps, measured and both bounds (dyadic and optimized SBC, QETU)
    fig, ax = plt.subplots(1, 1, figsize=(3.4, 3.2))
    E = sorted({p["eps"] for p in pts})
    off = {"dyadic": 0.88, "optimized": 1.0, "qetu": 1.13}
    mk = {"dyadic": "o", "optimized": "s", "qetu": "^"}
    for fam in ("dyadic", "optimized", "qetu"):
        F = [p for p in pts if p["fam"] == fam]
        for k, lab, mfc in (("dt", "measured", None), ("dt_comm", "commutator bound", "white")):
            for e in E:
                v = [p[k] for p in F if p["eps"] == e]
                x = np.full(len(v), e * off[fam]) * np.exp(np.random.default_rng(1).normal(0, 0.03, len(v)))
                ax.plot(x, v, mk[fam], ms=2.4, color=COL[fam], alpha=0.5 if k == "dt" else 0.35,
                        mfc=COL[fam] if mfc is None else "white", mew=0.6)
    for e in E:
        v = [p["dt_mk"] for p in D if p["eps"] == e]
        ax.plot([e * 0.8, e * 1.25], [np.median(v)] * 2, color=M.INK2, lw=1.2, ls="--")
    ax.axhline(min(p["dt_valid"] for p in D), color=M.INK2, lw=0.6, ls=":")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(2.5e-5, 4e-2)
    ax.set_xlabel(r"target infidelity $\epsilon$"); ax.set_ylabel(r"Strang step $\delta t$ ($1/J_1$)")
    h = [plt.Line2D([], [], color=COL[f], marker=mk[f], ls="", ms=3, label=l)
         for f, l in (("dyadic", "dyadic SBC"), ("optimized", "optimized SBC"), ("qetu", "QETU"))]
    h += [plt.Line2D([], [], color=M.INK2, marker="o", ls="", ms=3, label="measured (filled)"),
          plt.Line2D([], [], color=M.INK2, marker="o", mfc="white", ls="", ms=3, label="commutator bound (open)"),
          plt.Line2D([], [], color=M.INK2, ls="--", lw=1.2, label="low-energy bound (median)")]
    ax.set_ylim(3e-7, 4)
    ax.annotate("", xy=(5e-5, 6e-3), xytext=(5e-5, 2.5e-1), arrowprops=dict(arrowstyle="<->", lw=0.6, color=M.INK2))
    ax.text(4.3e-5, 4e-2, rf"$\times${MAC['bdCommMed']}", fontsize=6, ha="right", va="center", color=M.INK2)
    ax.text(3.3e-2, 1.4e-3, r"low-energy bound needs $\delta t$ below", fontsize=5.5, ha="right", va="top", color=M.INK2)
    fig.legend(handles=h, fontsize=5.6, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.55, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    fig.savefig(os.path.join(FIG, "fig_bounds_rev.pdf"))
    plt.close(fig)
    print("bounds: measured/commutator median", np.median(rc), "range", rc.min(), rc.max(),
          "; measured/MK median", np.median(rm), "range", rm.min(), rm.max())


# ----------------------------------------------------------------------------- 3. harder MPS baselines
def prep_cnot(n):
    return round((23 / 24) * 2 ** n - 1.5 * 2 ** (n / 2) + 4 / 3)


def baselines():
    """Re-score the benchmark against the MPS alone with (a) every staircase unitary costed at half its
    quantum-Shannon count (an optimistic stand-in for isometry decompositions) for both routes, and
    (b) generic state preparation of the exact ground state as a third route (Plesch-Brukner)."""
    A1 = jsonl("emp_actual1d.jsonl")
    A2 = [r for r in json.load(open(os.path.join(RES, "emp_actual.json")))]
    out = {}
    for tag, A in (("chain", A1 + [r for r in A2 if r["kind"] == "chain"]), ("cyl", [r for r in A2 if r["kind"] == "cyl"])):
        by = defaultdict(list)
        for r in A:
            by[(r["kind"], r["size"], r["J2"], r["eps"])].append(r)
        res = {}
        for scale in (1.0, 0.5):
            wins = tot = gen = 0
            for key, R in by.items():
                alone = [scale * r["trial_cnot"] for r in R if r["mps_alone"]]
                filt = []
                for r in R:
                    if r["mps_alone"]:
                        continue
                    for k in ("dyadic", "optimized"):
                        m = r.get(k)
                        if m:
                            fc = m["cnot"] * m["p"] - r["trial_cnot"]
                            filt.append((scale * r["trial_cnot"] + fc) / m["p"])
                if not filt:
                    continue
                tot += 1
                best_alone = min(alone) if alone else np.inf
                if min(filt) < best_alone:
                    wins += 1
                if prep_cnot(R[0]["N"]) < min(min(filt), best_alone):
                    gen += 1
            res[scale] = (wins, tot, gen)
        out[tag] = res
        macro(f"bl{tag.capitalize()}Wins", str(res[1.0][0])); macro(f"bl{tag.capitalize()}Tot", str(res[1.0][1]))
        macro(f"bl{tag.capitalize()}WinsHalf", str(res[0.5][0]))
        macro(f"bl{tag.capitalize()}Gen", str(res[1.0][2])); macro(f"bl{tag.capitalize()}GenHalf", str(res[0.5][2]))
    for n in (12, 14, 16):
        macro(f"genPrep{['', 'Twelve', '', 'Fourteen', '', 'Sixteen'][n - 11]}", f"{prep_cnot(n):,}".replace(",", "\\,"))
    print("baselines", out)



# ----------------------------------------------------------------------------- 4. truncation trials
def _fit_rows(R, fam):
    X, y, grp = {"shot": [], "succ": [], "tot": []}, {"shot": [], "succ": [], "tot": []}, []
    rows = []
    for r in R:
        m = r.get(fam)
        if not m:
            continue
        rows.append(r)
    out = {}
    for what in ("shot", "succ", "tot"):
        XX, yy, gg = [], [], []
        for r in rows:
            m = r[fam]
            c = dict(shot=m["cnot_filter"], succ=m["cnot_filter"] / m["p"], tot=m["cnot"])[what]
            XX.append([1, np.log(r["N"]), -np.log(r["gap"]), np.log(np.log(1 / r["eps"])), -np.log(r["gamma"])])
            yy.append(np.log(c))
            gg.append((r["N"], r["gap_target"]))
        out[what] = _fit_ci(np.array(XX), np.array(yy), gg)
    return out, len(rows)


def trunc():
    R = [r for r in jsonl("review_trunc.jsonl") if 1 - r["gamma"] > r["eps"]]
    if not R:
        return
    fams = (("dyadic", "D"), ("optimized", "O"), ("qetu_pi2", "Q"))
    res = {}
    for fam, t in fams:
        f, n = _fit_rows(R, fam)
        res[t] = f
        for what, w in (("shot", "Shot"), ("succ", "Succ"), ("tot", "Tot")):
            c, lo, hi, rms, mx = f[what]
            for i, nm in ((1, "N"), (2, "gap"), (3, "eps"), (4, "gam")):
                macro(f"tr{t}{w}{nm}", f"{c[i]:.2f}")
                macro(f"tr{t}{w}{nm}lo", f"{lo[i]:.2f}")
                macro(f"tr{t}{w}{nm}hi", f"{hi[i]:.2f}")
            macro(f"tr{t}{w}rms", f"{rms:.2f}")
    macro("trNcases", str(len(R)))
    macro("trGamLo", f"{min(r['gamma'] for r in R):.2f}"); macro("trGamHi", f"{max(r['gamma'] for r in R):.3f}")
    rE_t = [r["resid_E"] for r in R]
    rot = json.load(open(os.path.join(RES, "review_rot_resid.json"))) if os.path.exists(
        os.path.join(RES, "review_rot_resid.json")) else []
    rE_r = [r["resid_E"] for r in rot]
    macro("trResLo", f"{min(rE_t):.1f}"); macro("trResHi", f"{max(rE_t):.1f}")
    if rE_r:
        macro("rotResLo", f"{min(rE_r):.0f}"); macro("rotResHi", f"{max(rE_r):.0f}")
    # direct comparison at matched (N, gap, eps): log-cost residual of each truncation case against the
    # rotation-trial fit (controlled study), i.e. would the rotation law have predicted it?
    G = load_cases()
    X, y, grp = _design(G, "meas_dyadic", "shot_filter")
    c, *_ = _fit_ci(X, y, grp, nboot=1)
    ratio = []
    for r in R:
        m = r.get("dyadic")
        if not m:
            continue
        x = np.array([1, np.log(r["N"]), -np.log(r["gap"]), np.log(np.log(1 / r["eps"])), -np.log(r["gamma"])])
        ratio.append(m["cnot_filter"] / np.exp(x @ c))
    ratio = np.array(ratio)
    macro("trPredMed", f"{np.median(ratio):.2f}"); macro("trPredLo", f"{ratio.min():.2f}")
    macro("trPredHi", f"{ratio.max():.2f}")
    print("trunc:", len(R), "cases; dyadic shot", np.round(res["D"]["shot"][0], 2), "succ", np.round(res["D"]["succ"][0], 2),
          "; rotation-law prediction ratio median", np.median(ratio), ratio.min(), ratio.max())
    txt = r"""The rotation lowers the overlap by spin flips, which put the missing weight mostly at high energy: the residual
of a rotation trial has a mean excitation energy of \rotResLo--\rotResHi{} gaps. To check that the split between gap
and overlap is not specific to that kind of residual, we repeat study A with a different mechanism. The trial is the
exact ground state truncated to Schmidt rank one on a chosen set of bonds (at most five), an MPS of bond dimension at
most 16 costed like the $\chi=16$ trial, whose error is an entanglement truncation and whose residual sits low in the
spectrum, \trResLo--\trResHi{} gaps above $E_0$. For each $N$, $\Delta\in\{0.40,0.50\}$ and target $\gamma$ the bond
set whose overlap lies closest to the target is used, giving $\gamma=\trGamLo$--\trGamHi{} (\trNcases{} cases with
$1-\gamma>\epsilon$; \texttt{run\_review.py trunc}).

The same law holds. Fitting Eq.~\eqref{eq:fit} to the truncation trials gives, for the dyadic filter, a per-shot
overlap exponent $\nu=\trDShotgam$ $[\trDShotgamlo,\trDShotgamhi]$ and a per-success exponent \trDSuccgam{}
$[\trDSuccgamlo,\trDSuccgamhi]$, against \ciDShotgam{} and \ciDSuccgam{} for the rotation trials, with
$\lambda=\trDShoteps$ for the target. (The gap exponent is poorly constrained here, since only two gaps are used.)
Predicting each truncation case's per-shot CNOTs from the rotation-trial fit of Table~\ref{tab:gofit} is off by a
median factor \trPredMed{} (range \trPredLo--\trPredHi). The overlap therefore enters through its size, not through
where in the spectrum the missing weight sits: the filter has to suppress every level above the gap, wherever the
weight is."""
    open(os.path.join(GEN, "sec_trunc.tex"), "w").write(txt + "\n")



# ----------------------------------------------------------------------------- 5. monotonicity in the step
def mono():
    R = jsonl("review_mono.jsonl")
    if not R:
        return
    n = mono_ok = 0
    ratio, worst, qn, qok = [], [], 0, 0
    for r in R:
        for fam in ("dyadic", "optimized"):
            d = r.get(fam)
            if not d:
                continue
            sc = sorted(d["scan"], key=lambda x: -x[0])          # coarse to fine
            ok = [x[2] <= r["eps"] for x in sc]
            n += 1
            first = ok.index(True) if True in ok else None
            if first is not None and all(ok[first:]):
                mono_ok += 1
            if first is not None:
                ratio.append(sc[first][0] / d["dt_bisect"])
                worst.append(max(x[2] for x in sc[first:]) / r["eps"])
        q = r.get("qetu")
        if q:
            qn += 1
            ok = [x[1] <= r["eps"] for x in q["scan"]]
            first = ok.index(True) if True in ok else None
            qok += int(first is not None and all(ok[first:]))
    ratio, worst = np.array(ratio), np.array(worst)
    macro("monoN", str(n)); macro("monoOK", str(mono_ok)); macro("monoCases", str(len(R)))
    macro("monoQN", str(qn)); macro("monoQOK", str(qok))
    macro("monoRatioMed", f"{np.median(ratio):.2f}"); macro("monoRatioHi", f"{ratio.max():.2f}")
    macro("monoWorst", f"{worst.max():.2f}")
    macro("monoCoarser", str(int(np.sum(ratio > 1.05))))
    print("mono:", mono_ok, "/", n, "SBC scans monotone;", qok, "/", qn, "QETU; dense/bisect", np.median(ratio), ratio.max(),
          "; worst finer infidelity/eps", worst.max())
    txt = r"""We checked this on \monoCases{} cases ($N=8$ and 10; three overlaps at $\Delta=0.40$ and three gaps at
$\gamma=0.1$; all three targets) by scanning the common step on a dense logarithmic grid from one step per pulse down
to an eighth of the bisected step. For the dyadic filter and for QETU (one to six steps per query) passing was
monotone in every scan: once a step passes, every finer step passes. The exceptions are all in the optimized SBC
filter, \monoOK{} of \monoN{} SBC scans being monotone. Its pulse times are unequal, so the number of steps in each
pulse is rounded separately and the realized error does not fall smoothly with the step. In four of the six exceptions
the failing step lies above the bisected one, which only means the bisection could have stopped at a coarser step
(the dense scan found a passing step up to \monoRatioHi{} times coarser; median ratio \monoRatioMed). In two
($N=8$ and 10 at $\Delta=0.30$) a step finer than the bisected one fails, at 0.72 and 0.91 times it, with infidelity
up to \monoWorst{}$\,\epsilon$. The step we report is always one that was evaluated and passed, so the tables are
unaffected, but for the optimized filter a finer step is not automatically safe and the chosen step must be checked
directly."""
    open(os.path.join(GEN, "sec_mono.tex"), "w").write(txt + "\n")



# ----------------------------------------------------------------------------- 6. wrong inputs: E0, Emax, gamma margin
def sens2():
    p = os.path.join(RES, "review_sens2.json")
    if not os.path.exists(p):
        return
    C = json.load(open(p)) + [dict(r, kind="gamma") for r in json.load(open(os.path.join(RES, "ext_sensitivity.json")))
                              if r["kind"] == "gamma"]
    scheds = ("dyadic", "optimized", "qetu")
    lines = []
    spec = (("E0", (-0.25, -0.1, 0.1, 0.25), r"$E_0$ shift / $\Delta$"),
            ("Emax", (0.9, 0.95, 1.0, 1.25), r"$E_{\max}$ / exact"),
            ("gamma", (1.02, 1.05), r"$\gamma$, no margin"),
            ("gammaM", (1.0, 1.02, 1.05, 1.1), r"$\gamma$, margin 1.1"))
    for kind, fs, name in spec:
        for f in fs:
            cells = []
            for sc in scheds:
                rs = [r for r in C if r["kind"] == kind and r["f"] == f and r["schedule"] == sc]
                ok = [r for r in rs if r.get("sim")]
                miss = [r for r in ok if r["sim"]["infidelity_trotter"] > 1e-2]
                ref = len(rs) - len(ok)
                worst = max((r["sim"]["infidelity_trotter"] for r in ok), default=float("nan"))
                cells.append(f"{len(miss)}/{len(rs)}" + (f"$^{{{ref}}}$" if ref else "") + f" & ${M.sci(worst)}$")
            lines.append(f"{name} & {f:+g} & " if kind == "E0" else f"{name} & {f:g} & ")
            lines[-1] += " & ".join(cells) + " \\\\"
    with open(os.path.join(GEN, "tab_sens2.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    miss = lambda kind, f, sc=None: sum(1 for r in C if r["kind"] == kind and r["f"] == f and r.get("sim")
                                        and (sc is None or r["schedule"] == sc)
                                        and r["sim"]["infidelity_trotter"] > 1e-2)
    tot = lambda kind, f: sum(1 for r in C if r["kind"] == kind and r["f"] == f)
    for kind, fs in (("E0", (-0.25, -0.1, 0.1, 0.25)), ("Emax", (0.9, 0.95, 1.0, 1.25)), ("gammaM", (1.0, 1.02, 1.05, 1.1))):
        for f in fs:
            tag = {-0.25: "mQ", -0.1: "mT", 0.1: "pT", 0.25: "pQ", 0.9: "Nine", 0.95: "NineFive", 1.0: "One",
                   1.25: "OneQ", 1.02: "OneOTwo", 1.05: "OneOFive", 1.1: "OneOne"}[f]
            macro(f"sII{kind.replace('0', 'Zero')}{tag}Miss", str(miss(kind, f)))
            macro(f"sII{kind.replace('0', 'Zero')}{tag}N", str(tot(kind, f)))
    macro("sIIrefused", str(sum(1 for r in json.load(open(p)) if not r.get("sim"))))
    print("sens2:", {(k, f): (miss(k, f), tot(k, f)) for k, fs, _ in spec for f in fs})


# ----------------------------------------------------------------------------- 7. certified step rule
def eckart():
    R = jsonl("review_eckart.jsonl")
    if not R:
        return
    ok = cert = 0
    ratio, shots, redesign, ratio_q, ratio_s = [], [], 0, [], []
    for r in R:
        if not r["attempts"]:
            continue
        a = r["attempts"][-1]
        c = a["curve"][-1]
        if c["bound"] > r["eps"]:
            continue
        cert += 1
        ok += int(c["inf"] <= r["eps"])
        redesign += int(a["frac"] < 1)
        cost = (a["trial_cnot"] + c["filter_cnot"]) / c["p"]
        a1 = r["attempts"][0]
        cv = a1["curve"]
        orc = None
        for i in range(len(cv)):
            if all(q["inf"] <= r["eps"] for q in cv[i:]):
                orc = cv[i]
                break
        if orc is not None:
            rr = cost / ((a1["trial_cnot"] + orc["filter_cnot"]) / orc["p"])
            ratio.append(rr)
            (ratio_q if r["schedule"] == "qetu" else ratio_s).append(rr)
        margin = r["eps"] * r["gap_dmrg"] - (c["E"] - r["E0_dmrg"])
        shots.append(3 * 4 * sum(c["var"]) / margin ** 2 if margin > 0 else np.inf)
    shots = np.array(shots)
    macro("ekN", str(len(R))); macro("ekCert", str(cert)); macro("ekOK", str(ok)); macro("ekRedesign", str(redesign))
    macro("ekRatioS", f"{np.median(ratio_s):.1f}"); macro("ekRatioQ", f"{np.median(ratio_q):.1f}")
    macro("ekRatioHi", f"{max(ratio):.1f}")
    e = np.array([r["eps"] for r in R if r["attempts"] and r["attempts"][-1]["curve"][-1]["bound"] <= r["eps"]])
    macro("ekShotsMed", M.sci(np.median(shots[np.isfinite(shots)])))
    lo = shots[(e == 1e-2) & np.isfinite(shots)]
    if len(lo):
        macro("ekShotsLoose", M.sci(np.median(lo)))
    print("eckart:", cert, "/", len(R), "certified;", ok, "within eps;", redesign, "redesigned; cost ratio SBC",
          np.median(ratio_s), "QETU", np.median(ratio_q), "max", max(ratio), "; shots median", np.median(shots))



# ----------------------------------------------------------------------------- 8. SBC/QETU crossover target versus N
def qcross():
    from derive_gapov_filter import load_cases
    Q = [json.loads(l) for l in open(os.path.join(RES, "qetu_tau.jsonl"))]
    rows = {(r["N"], round(r["J2"], 7), r["gamma_target"], r["eps"]): r for r in load_cases()}
    best = lambda q: min([x["cnot"] for x in (q["tau_pi2"], q["3pi/4"], q["pi/(1+delta)"]) if x] or [np.inf])
    D = {}
    for q in Q:
        if q["src"] != "gapov":
            continue
        r = rows[(q["N"], round(q["J2"], 7), q["gamma_target"], q["eps"])]
        s = min(r["meas_dyadic"]["cnot"], (r["meas_optimized"] or {}).get("cnot", np.inf))
        D.setdefault((q["N"], q["eps"]), []).append(best(q) / s)
    lines = []
    names = {8: "Eight", 10: "Ten", 12: "Twelve", 14: "Fourteen", 16: "Sixteen"}
    for N in sorted({k[0] for k in D}):
        med = [np.median(D[(N, e)]) for e in (1e-2, 1e-3, 1e-4)]
        a, b = np.polyfit([-2, -3, -4], np.log(med), 1)
        xc = -b / a
        macro("qx" + names[N], M.sci(10 ** xc, 1))
        lines.append(f"{N} & " + " & ".join(f"{m:.2f}" for m in med) + f" & ${M.sci(10 ** xc, 1)}$ \\\\")
        print("qcross", N, np.round(med, 2), 10 ** xc)
    open(os.path.join(GEN, "tab_qcross.tex"), "w").write("\n".join(lines) + "\n")

if __name__ == "__main__":
    what = sys.argv[1:] or ["fits", "bounds"]
    for w in what:
        globals()[w]()
    write_macros()
