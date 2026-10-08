"""
run_extras.py -- supporting experiments (python-3.12 venv).   usage: python run_extras.py trials|sens|design|baseline

  trials   : which trial state?  DMRG at chi = 1,2,4,8 (exact circuits) and L-layer approximate circuits,
             N=8: gamma, trial CX, certified filter, total CX  (eps = 1e-2)
  sens     : what if the classical inputs are wrong?  gap over/under-estimated, gamma over-estimated
  design   : sensitivity to the design choices (leakage share eps_l/eps, number of pulses m)
  baseline : direct MPS-circuit preparation at equal fidelity (mps-to-circuit exact + approximate, generic
             state preparation), measured against the exact ground state
"""
import sys, json, pathlib, time
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import workflow as wf
import dmrg_inputs as DI
import run_experiments as R
import resources as RS

OUT = R.OUT
dump = lambda name, obj: (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def setup(N, J2):
    ref = DI.dmrg_reference(N, J2)
    Wd = ref["Etop"] - ref["E0"]
    ch = wf.Chain(N, J2, E0_used=ref["E0"], W_used=Wd)
    rd = DI.dense(ref["psi0"]); rd = rd / np.linalg.norm(rd)
    return ref, ch, rd, Wd, (ref["E1"] - ref["E0"]) / Wd


def trials():
    out = []
    eps = 0.01
    for N, J2 in [(8, 0.0), (8, 0.4)]:
        ref, ch, rd, Wd, dl = setup(N, J2)
        s1, _ = RS.cx_per_step(N, J2, Wd)
        steps = dict(a2a=int(s1))
        fams = [("dmrg", chi, DI.trial_dmrg(N, J2, chi)) for chi in (1, 2, 4, 8)] + \
               [("approx", L, DI.trial_approx(ref["psi0"], L)) for L in (1, 2, 3)]
        for kind, p, tr in fams:
            psi = tr["psi"]
            gam = float(abs(np.vdot(rd, psi)) ** 2)
            cx = R.Ctx(ch, psi, gam, dl)
            c = R.run_eps(cx, tr, eps, steps, sweep=False)
            out.append(dict(N=N, J2=J2, kind=kind, param=p, gamma_est=gam, gamma_true=cx.gamma_true, trial_cx=tr["cx"],
                            D_trial=c["D_trial"], note=c.get("note"), n_bound=c.get("n_bound"), n_hyb=c.get("n_hyb"),
                            n_meas=c.get("n_meas"), cx_total_bound=c.get("cx_total_bound"), cx_total_hyb=c.get("cx_total_hyb"),
                            cx_total_meas=c.get("cx_total_meas"), T=(c.get("bound_point") or {}).get("T")))
            print(out[-1], flush=True)
    dump("extra_trials", out)


def sens():
    out = []
    eps = 0.01
    for N, J2 in [(8, 0.0), (12, 0.0)]:
        ref, ch, rd, Wd, dl = setup(N, J2)
        tr = DI.trial_dmrg(N, J2, 2)
        psi = tr["psi"]
        gt = float(abs(np.vdot(ch.g, psi)) ** 2)
        variants = [("gap", d, gt, dl * (1 + d)) for d in (-0.2, 0.0, 0.1, 0.25, 0.5, 1.0)] + \
                   [("gamma", f, gt + f * (1 - gt), dl) for f in (0.25, 0.5, 0.75)]
        for kind, p, gam, dlt in variants:
            cx = R.Ctx(ch, psi, gam, dlt)
            base = R.design(gam, dlt, eps)
            nb = R.bound_n(cx, base, eps)
            o, pk = R.point(cx, base, nb)
            m = R.measure(cx, dict(o), pk) if cx.cost(nb) <= R.SIM_BUDGET else dict(o)
            row = dict(N=N, J2=J2, kind=kind, param=p, gamma_used=gam, gamma_true=gt, delta_used=dlt, delta_true=float(ch.gap),
                       n_bound=nb, bound_claimed=o["bound"], bound_with_truth=o["bound_truth"], eta_used=o["eta"],
                       eta_true=o["eta_true"], meas=m.get("meas"), meas_leak=m.get("meas_leak"), sqrt_eps=float(np.sqrt(eps)))
            row["claim_violated"] = bool(row["meas"] > row["bound_claimed"] + 1e-12) if row["meas"] is not None else None
            row["target_missed"] = bool(row["meas"] > row["sqrt_eps"]) if row["meas"] is not None else None
            out.append(row)
            print(row, flush=True)
    dump("extra_sens", out)


def design_sens():
    out = []
    eps = 0.01
    N, J2 = 8, 0.0
    ref, ch, rd, Wd, dl = setup(N, J2)
    tr = DI.trial_dmrg(N, J2, 2)
    gam = float(abs(np.vdot(rd, tr["psi"])) ** 2)
    cx = R.Ctx(ch, tr["psi"], gam, dl)
    s1, _ = RS.cx_per_step(N, J2, Wd)
    for frac in (0.05, 0.1, 0.25, 0.5, 0.75):
        for ml in ((4,), (6,), (8,)):
            try:
                rows_, _ = R.FL.floor_vs_precision([eps * frac], gam, dl, 0.0, x_list=R.XL, m_list=ml, hi=1.0, verbose=False)
                if not rows_[0]["feasible"]:
                    out.append(dict(frac=frac, m=ml[0], error="infeasible")); print(out[-1], flush=True); continue
                base = rows_[0]
                nb = R.bound_n(cx, base, eps)
                o, _ = R.point(cx, base, nb)
                out.append(dict(frac=frac, m=ml[0], m_used=int(base["m"]), x=float(base["x"]), T=float(base["T"]), f0sq=float(base["f0sq"]),
                                eta=float(base["eta_cert"]), n_bound=nb, cx_filter=nb * int(s1), bound_leak=o["bound_leak"]))
            except Exception as e:
                out.append(dict(frac=frac, m=ml[0], error=repr(e)))
            print(out[-1], flush=True)
    dump("extra_design", out)


def baseline():
    from qiskit import QuantumCircuit
    from qiskit.circuit.library import StatePreparation
    from mps_to_circuit import mps_to_circuit
    out = []
    for N in (4, 6, 8, 10, 12, 14, 16):
        for J2 in (0.0, 0.4):
            ref, ch, rd, Wd, dl = setup(N, J2)
            row = dict(N=N, J2=J2, exact=[], approx=[])
            for chi in ((1, 2, 4, 8, 16) if N <= 12 else (1, 2, 4, 8)):
                psi = ref["psi0"].copy(); psi.compress(max_bond=chi, cutoff=0.0); psi.normalize()
                try:
                    qc = mps_to_circuit(psi.arrays, method="exact", shape="lpr")
                    v = DI.circuit_state(qc)
                    row["exact"].append(dict(chi=chi, cx=DI.cx_count(qc), infid=float(1 - abs(np.vdot(ch.g, v)) ** 2)))
                except Exception as e:
                    row["exact"].append(dict(chi=chi, error=repr(e)[:80]))
            for L in range(1, 9):
                try:
                    t = DI.trial_approx(ref["psi0"], L)
                    row["approx"].append(dict(L=L, cx=t["cx"], infid=float(1 - abs(np.vdot(ch.g, t["psi"])) ** 2)))
                except Exception as e:
                    break
            if N <= 12:
                qc = QuantumCircuit(N)
                from hamiltonians import reorder_axes
                qc.append(StatePreparation(reorder_axes(ch.g).astype(complex)), range(N))
                row["generic"] = dict(cx=DI.cx_count(qc))
            else:
                row["generic"] = dict(cx=None)
            out.append(row)
            print(N, J2, [(e.get("chi"), e.get("cx"), round(e.get("infid", -1), 6)) for e in row["exact"]], flush=True)
    dump("extra_baseline", out)


if __name__ == "__main__":
    {"trials": trials, "sens": sens, "design": design_sens, "baseline": baseline}[sys.argv[1]]()
