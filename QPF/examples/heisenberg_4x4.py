"""The paper's Heisenberg example: 4 x 4 periodic lattice, h = 0, J_z = 0 sector, Neel initial state.

    python examples/heisenberg_4x4.py          # ~2 min; writes figs/fig1_filters.png, figs/fig2_heisenberg.png
                                               # and results/heisenberg_4x4.json

Steps, as in the paper:
  1. J^2 projection onto J = 0: halving ladder from t1 = pi/4, zero phases (three steps for J <= 8).
  2. Energy projection on the scaled H (spectrum 0..1, target at 0), total time pi/Delta with Delta = 0.15:
     constant times (m chosen for the lowest final energy), Gaussian-sampled times (7 steps, averaged over draws),
     QPF optimised times, QPF optimised times and phases (six steps each, designed for Delta_e = 0.15 and a 5%
     ground-state probability).
Paper values to compare: Neel scaled energy ~0.17; J^2 projection needs three steps, succeeds with ~11% and
lowers the energy to ~0.06; the relevant gap goes from ~0.03 to ~0.15.
"""
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import eigsh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from qpf import core, design, models, symmetry  # noqa: E402

DELTA_E, ALPHA_E, M_QPF, M_GAUSS, N_GAUSS = 0.15, 0.05, 6, 7, 20
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def system():
    sec = models.SpinSector(16)
    H = models.heisenberg(sec, models.square_bonds(4, 4))
    J2 = models.total_spin_sq(sec)
    w, V = eigsh(H, k=2, which="SA", tol=1e-12)
    o = np.argsort(w)
    E0, g0 = w[o[0]], V[:, o[0]]
    Emax = eigsh(H, k=1, which="LA", return_eigenvectors=False)[0]
    W = Emax - E0
    ws = np.sort(eigsh(H + W * J2, k=2, which="SA", return_eigenvectors=False))
    return dict(sec=sec, H=H, J2=J2, E0=E0, Emax=Emax, W=W, g0=g0, gap_all=(w[o[1]] - E0) / W,
                gap_J0=(ws[1] - ws[0]) / W, Hs=core.scale(H, E0, Emax))


def trajectory(Hs, g0, psi, times, phases, p0):
    """Cumulative time t*Delta/pi, scaled energy, success probability and ground fidelity after each step."""
    r = core.run(Hs, psi, times, phases, keep=True)
    Op = core.Operator(Hs)
    T = np.cumsum(np.abs(times)) * DELTA_E / np.pi            # running sum of |t_i| (Gaussian times can be negative)
    E = [Op.expect(v) for v in r.trajectory]
    P = p0 * np.cumprod(r.p_steps)
    F = [abs(np.vdot(g0, v)) ** 2 for v in r.trajectory]
    return dict(t=[0.0] + list(T), E=[Op.expect(psi)] + E, P=[p0] + list(P), F=[abs(np.vdot(g0, psi)) ** 2] + F)


def main():
    t0 = time.time()
    s = system()
    psi = models.neel(s["sec"], 4, 4)
    Hop = core.Operator(s["Hs"])
    out = dict(neel_energy=Hop.expect(psi), neel_overlap_sq=abs(np.vdot(s["g0"], psi)) ** 2,
               gap_all=s["gap_all"], gap_J0=s["gap_J0"])

    # 1. J^2 projection
    ts = symmetry.j2_ladder(0, 8)
    r = core.run(s["J2"], psi, ts)
    out.update(j2_times=ts, j2_success=r.p_total, j2_energy=Hop.expect(r.state),
               j2_J2=core.Operator(s["J2"]).expect(r.state), j2_overlap_sq=abs(np.vdot(s["g0"], r.state)) ** 2)
    psiJ, pJ = r.state, r.p_total
    print(f"Neel: E = {out['neel_energy']:.3f}; gap {out['gap_all']:.3f} -> {out['gap_J0']:.3f} in J = 0")
    print(f"J^2 projection, {len(ts)} steps: success {pJ:.4f}, E = {out['j2_energy']:.4f}, <J^2> = {out['j2_J2']:.1e}")

    # 2. energy projection, total time pi / Delta_e
    T = np.pi / DELTA_E
    sched = {}
    best_m = min(range(1, 31), key=lambda m: Hop.expect(core.run(s["Hs"], psiJ, *design.constant(T, m)).state))
    sched["Constant"] = [design.constant(T, best_m)]
    sched["Gaussian"] = [design.gaussian(T, M_GAUSS, seed=k) for k in range(N_GAUSS)]
    q_t = design.optimize(DELTA_E, ALPHA_E, M_QPF, T, phases=False)
    q_td = design.optimize(DELTA_E, ALPHA_E, M_QPF, T, phases=True)
    sched["QPF (t_i)"] = [(q_t["times"], q_t["phases"])]
    sched["QPF (t_i, δ_i)"] = [(q_td["times"], q_td["phases"])]
    out["constant_m"] = best_m
    out["qpf_t"] = dict(times=list(q_t["times"]), design_loss=q_t["loss"], design_fraction=q_t["fraction"])
    out["qpf_td"] = dict(times=list(q_td["times"]), phases=list(q_td["phases"]), design_loss=q_td["loss"], design_fraction=q_td["fraction"],
                         f0=float(np.prod(np.cos(q_td["phases"]))))

    traj = {}
    for name, runs in sched.items():
        tr = [trajectory(s["Hs"], s["g0"], psiJ, t, d, pJ) for t, d in runs]
        traj[name] = {k: list(np.mean([x[k] for x in tr], axis=0)) for k in tr[0]}
        e = traj[name]
        print(f"{name:16s} steps {len(runs[0][0]):2d}: E = {e['E'][-1]:.2e}, success {e['P'][-1]:.4f}, "
              f"fidelity {e['F'][-1]:.6f}")
    out["energy_projection"] = {k: dict(final_energy=v["E"][-1], final_success=v["P"][-1], final_fidelity=v["F"][-1])
                                for k, v in traj.items()}

    os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    style = {"Constant": ("o", "#c0392b"), "Gaussian": ("^", "#8e44ad"), "QPF (t_i)": ("s", "#1e8449"),
             "QPF (t_i, δ_i)": ("D", "#000000")}

    # Fig. 1: filter functions on [0, 1]
    E = np.linspace(0, 1, 2001)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for name, runs in sched.items():
        t, d = runs[0]
        ax.plot(E, core.filter_function(E, t, d), color=style[name][1], label=name, lw=1.2)
    ax.axvline(DELTA_E, color="gray", ls=":", lw=0.8)
    ax.set_xlabel("Energy (scaled)"); ax.set_ylabel("f"); ax.legend(fontsize=8)
    ax.set_title(f"Filter functions, Δe = {DELTA_E}, αe = {ALPHA_E}, total time π/Δe", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "fig1_filters.png"), dpi=150); plt.close(fig)

    # Fig. 2: success probability and energy vs t Delta / pi
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(5.5, 6), sharex=True, gridspec_kw=dict(height_ratios=(1, 2)))
    for name, e in traj.items():
        m, c = style[name]
        a1.plot(e["t"], e["P"], m + "-", color=c, ms=4, label=name)
        a2.semilogy(e["t"], np.maximum(e["E"], 1e-12), m + "-", color=c, ms=4)
    a1.set_ylabel("Prob"); a1.legend(fontsize=8); a2.set_ylabel("Energy (scaled)"); a2.set_xlabel("tΔ/π")
    a1.set_title("Heisenberg 4×4 after J² projection", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figs", "fig2_heisenberg.png"), dpi=150); plt.close(fig)

    json.dump(dict(out, trajectories=traj), open(os.path.join(HERE, "results", "heisenberg_4x4.json"), "w"),
              indent=1, default=float)
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
