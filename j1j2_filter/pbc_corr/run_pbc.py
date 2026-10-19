"""run_pbc.py -- filter-prepared ground states of the PERIODIC J1-J2 chain and their SzSz correlations (benchmark of the pipeline).

For each (N, J2): exact classical inputs (sector ED) -> trial circuit (L = 1) -> pipeline floor design + certified n -> exact cos-filter at the
certified design (tg, ph) applied to the trial state -> translation-averaged <Sz_i Sz_{i+r}>, r = 1..N/2, compared with the exact ground state.
The circuit-level (Trotter) state is within the certified trace distance of the exact-filter state; for N <= --trotter-max it is also simulated.
Gate counts (RZZ basis) = n x per-step + trial circuit.   Output: results/pbc_N{N}_J2{J2}_eps{eps}.json  (resumable)

usage: python pbc_corr/run_pbc.py --N 8 12 --J2 0.2 0.24 0.28 --eps 1e-3 [--trotter-max 10] [--outdir results]
"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
import pipeline as P                                                   # noqa: E402
from pbc_ctx import build_ctx, correlations                            # noqa: E402
from core.simulate import postselected_run_exact                        # noqa: E402
from core.circuits import apply_filter_pulse                            # noqa: E402
from core.resources import resource_costs, RZZ_BASIS                    # noqa: E402
from core.trotter import _trotter_run, trotter_bounds                   # noqa: E402
from qiskit import QuantumCircuit                                       # noqa: E402

P.SIM_MAX_N = 0                                                         # the pipeline's own Trotter simulation is not used for the design search


def step_counts(ctx):
    N = ctx["N"]; out = {}
    for tag, lvl in (("l3", 3), ("direct", 1)):
        qc = QuantumCircuit(N + 1); apply_filter_pulse(qc, list(range(N)), N, ctx["H_scaled"], 1.0, 0.0, 1, P.ORDER)
        c = resource_costs(qc, verbose=False, basis=RZZ_BASIS, optimization_level=lvl)
        out[tag] = {k: float(c[k]) for k in ("rzz", "rzz_nonclifford", "rz_nonclifford", "nc_total", "depth")}
    t = resource_costs(ctx["trial_qc"], verbose=False, basis=RZZ_BASIS)
    out["trial"] = {k: float(t[k]) for k in ("rzz", "rzz_nonclifford", "rz_nonclifford", "nc_total", "depth")}
    return out


def run_point(N, J2, eps, trotter_max, outdir):
    f = outdir / f"pbc_N{N}_J2{J2:g}_eps{eps:g}.json"
    if f.exists():
        return
    t0 = time.time(); ctx = build_ctx(N, 1.0, J2); t_ctx = time.time() - t0
    des, cands, n_inf = P.make_design(ctx, eps)
    if des["source"] == "none":
        raise RuntimeError("trial state already exact")
    pt = P.find_guaranteed(ctx, des, eps)
    if pt is None:
        json.dump(dict(N=N, J2=J2, eps=eps, error="no certified design", gamma=ctx["gamma"], gap_scaled=ctx["gap"]), open(f, "w")); return
    t1 = time.time(); p_ex, sv = postselected_run_exact(ctx["trial_qc"], pt["tg"], pt["ph"], ctx["Hs"]); t_filt = time.time() - t1
    psi = ctx["ed"]["psi"]; F = float(abs(np.vdot(psi, sv)) ** 2 / np.vdot(sv, sv).real)
    rec = dict(N=N, J2=J2, eps=eps, E0=ctx["ed"]["E0"], E1=ctx["ed"]["E1"], E1_singlet=ctx["ed"]["E1_singlet"], E1_triplet=ctx["ed"]["E1_triplet"],
               W=ctx["W"], gap_scaled=ctx["gap"], gamma=ctx["gamma"], alpha=ctx["alpha"], chi=ctx["chi"],
               n=int(pt["n"]), n_pulses=int(pt["n_nz"]), k=[int(x) for x in pt["k"]], tg=[float(x) for x in pt["tg"]], phases=[float(x) for x in pt["ph"]],
               T_scaled=float(np.sum(pt["tg"])), T_phys=float(np.sum(pt["tg"]) / ctx["W"]), eta=float(pt["eta"]), eps_bound=float(pt["eps_bound"]),
               p_g_lb=float(pt["p_g_lb"]), p_succ_exact=float(p_ex), F_exact_filter=F, infidelity=1 - F,
               C_exact=correlations(psi, N), C_filter=correlations(sv, N), C_trial=correlations(ctx["trial_vec"], N),
               E_filter=float(np.real(np.vdot(sv, ctx["H_qk"].to_matrix(sparse=True) @ sv)) / np.vdot(sv, sv).real) if N <= 14 else None,
               seconds=dict(ctx=t_ctx, filter=t_filt))
    sc = step_counts(ctx); n, nz = pt["n"], pt["n_nz"]
    rec["counts"] = dict(per_step=sc, rzz_l3=n * sc["l3"]["rzz"] + sc["trial"]["rzz"], nc_l3=n * sc["l3"]["nc_total"] + nz + sc["trial"]["nc_total"],
                         rzz_direct=n * sc["direct"]["rzz"], nc_direct=n * sc["direct"]["nc_total"] + nz)
    if N <= trotter_max:
        t2 = time.time(); svt, pt_ = _trotter_run(ctx["H_scaled"], ctx["trial_vec"], pt["tg"], pt["ph"], pt["k"])
        ov = abs(np.vdot(sv, svt)) / (np.linalg.norm(sv) * np.linalg.norm(svt))
        rec["trotter"] = dict(C=correlations(svt, N), state_err=float(np.sqrt(max(0.0, 1 - ov ** 2))), bound=float(trotter_bounds(ctx["alpha"], pt["tg"], pt["k"], p_ex)["bound_state"]),
                              infidelity=float(1 - abs(np.vdot(psi, svt)) ** 2 / np.vdot(svt, svt).real), seconds=time.time() - t2)
    f.write_text(json.dumps(rec, indent=1))
    ce, cf = np.array(rec["C_exact"]), np.array(rec["C_filter"])
    print(f"N={N} J2={J2:g} eps={eps:g}: gamma={ctx['gamma']:.3f} Delta={ctx['gap']:.4f} n={rec['n']} m={rec['n_pulses']} 1-F={rec['infidelity']:.1e} "
          f"max|dC|/|C|={np.max(abs(cf - ce) / abs(ce)):.1e} ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--N", type=int, nargs="+", default=[8]); ap.add_argument("--J2", type=float, nargs="+", default=[0.24])
    ap.add_argument("--eps", type=float, nargs="+", default=[1e-3]); ap.add_argument("--trotter-max", type=int, default=0); ap.add_argument("--outdir", default="results")
    a = ap.parse_args(); out = HERE / a.outdir; out.mkdir(exist_ok=True)
    for N in a.N:
        for J2 in a.J2:
            for eps in a.eps:
                run_point(N, J2, eps, a.trotter_max, out)
