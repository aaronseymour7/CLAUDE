"""The QPF pipeline, as the paper prescribes, on the open J1-J2 chain.

    python j1j2/run_pipeline.py                      # N = 8..16 -> results/j1j2_pipeline.json, figs/j1j2_*.png
    python j1j2/run_pipeline.py --N 8 10 --seed 1

Per N:
  1. classical inputs: DMRG E0, Emax and the lowest excited singlet E1 (j1j2/inputs.py); the trial state is the
     default one-layer mps-to-circuit circuit for the DMRG ground state.
  2. J^2 projection onto J = 0 (qpf.symmetry.j2_ladder, zero phases, t = pi/4, pi/8, ...). This also removes every
     S^z != 0 component the trial circuit leaks into.
  3. energy projection on H_s = (H - E0) / (Emax - E0) with the DMRG E0, Emax: the authors' design
     (qpf.design.optimize) for Delta_e = (E1 - E0) / (Emax - E0) from DMRG, total time pi / Delta_e, with
     method v3 (times and phases, f(0)^2 >= 0.9) and v4 (times only); m is the number of steps with the lowest
     design cost in M_RANGE. As in the paper, the block of steps is iterated ("iterating with the same times")
     until 1 - F <= EPS_STOP or R_MAX blocks. Constant and Gaussian schedules of the same total time per block are
     run for comparison.
Exact evolution throughout (no Trotter). Fidelity, energy error and the error of the estimates are measured
against ED, which the pipeline itself never uses.
"""
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from j1j2 import inputs  # noqa: E402
from qpf import core, design, symmetry  # noqa: E402

J1, J2 = 1.0, 0.2411
M_RANGE = range(4, 13)
R_MAX, EPS_STOP = 12, 1e-4
N_GAUSS = 10


def pick_m(gap, method):
    best = None
    for m in M_RANGE:
        q = design.optimize(gap, m, method=method)
        if best is None or q["cost"] < best["cost"] - 1e-9:
            best = dict(q, m=m)
    o = np.argsort(best["times"])                                   # paper: shortest time first
    best["times"], best["phases"] = best["times"][o], best["phases"][o]
    return best


def iterate(Hs, H, g0, E0, W, psi, p0, t0, times, phases):
    """Apply the block repeatedly; record after every block."""
    Op = core.Operator(Hs)
    v, p, t = psi, p0, t0
    rows = [dict(block=0, t=t, p=p, infid=1 - abs(np.vdot(g0, v)) ** 2, dE=core.Operator(H).expect(v) - E0)]
    for r in range(1, R_MAX + 1):
        run = core.run(Op, v, times, phases)
        v, p, t = run.state, p * run.p_total, t + float(np.sum(np.abs(times)))
        rows.append(dict(block=r, t=t, p=p, infid=1 - abs(np.vdot(g0, v)) ** 2,
                         dE=core.Operator(H).expect(v) - E0))
        if rows[-1]["infid"] <= EPS_STOP:
            break
    for x in rows:
        x["dE_scaled"] = x["dE"] / W
    return rows


def first(rows, eps, key="t"):
    for x in rows:
        if x["infid"] <= eps:
            return x[key]
    return None


def run_N(N, seed):
    t_start = time.time()
    ex = inputs.exact(N, J1, J2)
    est = inputs.dmrg_estimates(N, J1, J2, seed=seed)
    qc, psi = inputs.trial(est["psi0"])
    H, g0 = ex["H"], ex["g0"]
    W_est = est["Emax"] - est["E0"]
    gap_est = (est["E1_singlet"] - est["E0"]) / W_est
    Hs = core.scale(H, est["E0"], est["Emax"])
    out = dict(N=N, seed=seed,
               exact=dict(E0=ex["E0"], Emax=ex["Emax"], E1=ex["E1"], E1_singlet=ex["E1_singlet"],
                          gap_singlet_scaled=(ex["E1_singlet"] - ex["E0"]) / (ex["Emax"] - ex["E0"]),
                          gap_any_scaled=(ex["E1"] - ex["E0"]) / (ex["Emax"] - ex["E0"])),
               dmrg=dict(E0=est["E0"], Emax=est["Emax"], E1_singlet=est["E1_singlet"], gap_scaled=gap_est),
               trial=dict(gamma=abs(np.vdot(g0, psi)) ** 2, gamma_dmrg=abs(np.vdot(inputs.mps_statevector(est["psi0"]), psi)) ** 2,
                          dE=core.Operator(H).expect(psi) - ex["E0"], two_qubit_gates=int(qc.num_nonlocal_gates()),
                          leak_out_of_Sz0=float(np.sum(abs(psi[np.array([bin(i).count("1") for i in range(1 << N)]) != N // 2]) ** 2))))

    # 2. J^2 projection
    J2op = inputs.total_spin_sq(N)
    ts = symmetry.j2_ladder(0, N // 2)
    rj = core.run(J2op, psi, ts)
    out["j2"] = dict(times=ts, T=float(sum(ts)), success=rj.p_total, gamma_after=abs(np.vdot(g0, rj.state)) ** 2,
                     dE=core.Operator(H).expect(rj.state) - ex["E0"])

    # 3. energy projection. After the J^2 projection the state lies in J = 0, inside S^z = 0, which H conserves:
    #    the rest of the run is restricted to S^z = 0 (exact; only to save time).
    sz0 = np.flatnonzero(np.array([bin(i).count("1") for i in range(1 << N)]) == N // 2)
    out["j2"]["weight_outside_Sz0"] = float(np.sum(abs(np.delete(rj.state, sz0)) ** 2))
    H, Hs, g0 = H[sz0][:, sz0], Hs[sz0][:, sz0], g0[sz0]
    psiJ = rj.state[sz0] / np.linalg.norm(rj.state[sz0])
    T = np.pi / gap_est
    scheds = {}
    for label, method in (("QPF (t_i, δ_i)", "v3"), ("QPF (t_i)", "v4")):
        q = pick_m(gap_est, method)
        scheds[label] = [(q["times"], q["phases"])]
        out[f"design_{method}"] = dict(m=q["m"], cost=q["cost"], f0=q["f0"], times=list(q["times"]),
                                       phases=list(q["phases"]))
    m_c = out["design_v4"]["m"]
    scheds["Constant"] = [design.constant(T, m_c)]
    scheds["Gaussian"] = [design.gaussian(T, 7, seed=k) for k in range(N_GAUSS)]
    out["energy"] = {}
    for label, runs in scheds.items():
        allrows = [iterate(Hs, H, g0, ex["E0"], ex["Emax"] - ex["E0"], psiJ, rj.p_total, out["j2"]["T"], t, d)
                   for t, d in runs]
        rows = allrows[0] if len(allrows) == 1 else _mean_rows(allrows)
        out["energy"][label] = dict(T_block=T, rows=rows,
                                    t_at_1e2=first(rows, 1e-2), t_at_1e3=first(rows, 1e-3),
                                    p_at_1e2=first(rows, 1e-2, "p"), p_at_1e3=first(rows, 1e-3, "p"))
    out["seconds"] = time.time() - t_start
    return out


def _mean_rows(allrows):
    n = min(len(r) for r in allrows)
    return [{k: float(np.mean([r[i][k] for r in allrows])) for k in allrows[0][0]} for i in range(n)]


def figure(res, path):
    Ns = [r["N"] for r in res]
    fig, axs = plt.subplots(2, len(Ns), figsize=(3.4 * len(Ns), 5.6), sharey="row", squeeze=False)
    style = {"Constant": ("o", "#c0392b"), "Gaussian": ("^", "#8e44ad"), "QPF (t_i)": ("s", "#1e8449"),
             "QPF (t_i, δ_i)": ("D", "#000000")}
    for c, r in enumerate(res):
        for label, e in r["energy"].items():
            m, col = style[label]
            t = [x["t"] for x in e["rows"]]
            axs[0, c].plot(t, [x["p"] for x in e["rows"]], m + "-", color=col, ms=3.5, label=label)
            axs[1, c].semilogy(t, [max(x["infid"], 1e-16) for x in e["rows"]], m + "-", color=col, ms=3.5)
        axs[0, c].set_title(f"N = {r['N']}  (γ = {r['trial']['gamma']:.2f}, Δe = {r['dmrg']['gap_scaled']:.3f})", fontsize=9)
        axs[1, c].set_xlabel("total time (scaled units), incl. J² stage")
        for a in axs[:, c]:
            a.grid(alpha=0.3, which="both")
    axs[0, 0].set_ylabel("success probability")
    axs[1, 0].set_ylabel("1 − F")
    axs[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[8, 10, 12, 14, 16])
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    res = []
    for N in a.N:
        r = run_N(N, a.seed)
        res.append(r)
        e = r["energy"]
        print(f"N={N:2d}  gamma={r['trial']['gamma']:.3f}  J2: p={r['j2']['success']:.3f} gamma->{r['j2']['gamma_after']:.3f}  "
              f"gap_e {r['dmrg']['gap_scaled']:.4f} (exact {r['exact']['gap_singlet_scaled']:.4f})  "
              + "  ".join(f"{k}: T(1e-3)={v['t_at_1e3'] if v['t_at_1e3'] is None else round(v['t_at_1e3'], 1)} p={v['p_at_1e3'] if v['p_at_1e3'] is None else round(v['p_at_1e3'], 3)}"
                          for k, v in e.items())
              + f"  [{r['seconds']:.0f} s]", flush=True)
        os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
        json.dump(res, open(os.path.join(HERE, "results", "j1j2_pipeline.json"), "w"), indent=1, default=float)
    os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
    figure(res, os.path.join(HERE, "figs", "j1j2_pipeline.png"))


if __name__ == "__main__":
    main()
