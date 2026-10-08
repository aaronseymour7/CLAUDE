"""
run_experiments.py -- main experiments (run in the python-3.12 venv: needs quimb, qiskit, mps-to-circuit).

Per (N, J2):
  1. classical front end: DMRG E0, E1 (penalty), E_top, reference MPS (bond dim up to 64); DMRG trial state at
     bond dimension chi=2, compiled to an exact circuit (mps-to-circuit)
  2. the workflow sees ONLY DMRG quantities: W, Delta from the DMRG energies, gamma from the overlap of the trial
     with the DMRG reference. Exact diagonalisation is used afterwards to validate (ground state, true gap, true gamma).
  3. for each target eps: design the filter, find n_bound (smallest n with certified D_bar <= sqrt(eps)),
     n_hyb (non-oracle hybrid rule: certified leakage + n-vs-2n Trotter estimate), n_meas (oracle),
     simulate, count CX (filter step transpiled once; total = trial + n * step).
usage: python run_experiments.py [N ...]      (default: all sizes, 4 worker processes)
"""
import sys, json, time, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import workflow as wf
import dmrg_inputs as DI
import resources as RS
FL = wf.FL

OUT = pathlib.Path(__file__).resolve().parent / "results"
XL, ML = (0.6, 1.0, 1.5, 2.0), (4, 6)
SIM_BUDGET = 1.5e10        # n * bonds * 2 * 2^N : cap on a single simulated filter run


def design(gamma, delta, eps, frac=0.25, xl=XL, ml=ML):
    """Filter design with a widening fallback (more pulses / longer evolution) if no design is feasible."""
    for xs, ms in [(xl, ml), (xl + (2.5, 3.0), tuple(sorted(set(ml) | {8}))), (xl + (2.5, 3.0, 4.0), tuple(sorted(set(ml) | {8, 10})))]:
        rows, _ = FL.floor_vs_precision([eps * frac], gamma, delta, 0.0, x_list=xs, m_list=ms, hi=1.0, verbose=False)
        if rows[0]["feasible"]:
            return rows[0]
    raise RuntimeError("no feasible filter")


class Ctx:
    """Everything the workflow knows (est.) and everything only ED knows (true)."""

    def __init__(self, ch, psi, gamma_est, delta_est):
        self.ch, self.psi = ch, psi
        self.gamma_est, self.delta_est = gamma_est, delta_est
        self.gamma_true = float(abs(np.vdot(ch.g, psi)) ** 2)
        self.cost = lambda n, bonds=len(ch.bonds): n * bonds * 2 * 2 ** ch.N


def point(cx, base, n):
    """Bound for the base design snapped to n steps (no simulation)."""
    ch = cx.ch
    g = FL.grid_design(base["times"], base["phases"], base["T"], n, cx.delta_est, base["eta_target"], 0.0, 1.0)
    c = g["cert"]
    ks, dt, ph = g["k"], g["dt"], g["phases"]
    tt = ks * dt
    Db, Dbl, Dbt = wf.trace_distance_bound(cx.gamma_est, c["eta"], c["f0sq"], ch.alpha, tt.sum(), n)
    eta_true = FL.certify(tt, ph, ch.gap, 0.0, 1.0, rtol=0.005)["eta"] if abs(ch.gap - cx.delta_est) > 1e-12 else c["eta"]
    Dt_, _, _ = wf.trace_distance_bound(cx.gamma_true, eta_true, c["f0sq"], ch.alpha, tt.sum(), n)
    out = dict(n=int(n), k=[int(k) for k in ks], T=float(tt.sum()), eta=float(c["eta"]), f0sq=float(c["f0sq"]),
               feasible=bool(g["feasible"]), bound=float(Db), bound_leak=float(Dbl), bound_trot=float(Dbt),
               bound_truth=float(Dt_), eta_true=float(eta_true))
    return out, (ks, dt, ph, tt)


def measure(cx, out, pack, with_rich=False):
    ch, psi = cx.ch, cx.psi
    ks, dt, ph, tt = pack
    a = wf.run_filter_exact(ch, psi, tt, ph)
    b = wf.run_filter_trotter(ch, psi, ks, dt, ph)
    out.update(meas=wf.dist(b, ch.g), meas_leak=wf.dist(a, ch.g), meas_trot=wf.dist(a, b),
               psucc=float(np.linalg.norm(b) ** 2), psucc_exact=float(np.linalg.norm(a) ** 2),
               psucc_lb=float(cx.gamma_true * out["f0sq"]))
    bn = b / np.linalg.norm(b)
    out["energy_err"] = float(np.vdot(bn, ch.H @ bn).real - ch.E0_exact)
    out["energy_err_bound"] = float(ch.W * out["bound_truth"] ** 2)
    out["triangle_ok"] = bool(out["meas"] <= out["meas_leak"] + out["meas_trot"] + 1e-12)
    out["bound_ok"] = bool(out["meas"] <= out["bound_truth"] + 1e-12)
    out["leak_ok"] = bool(out["meas_leak"] <= np.sqrt(wf.leakage_bound(cx.gamma_true, out["eta_true"])) + 1e-12)
    pg = cx.gamma_true * out["f0sq"]
    out["trot_ok"] = bool(out["meas_trot"] <= min(1.0, ch.alpha * out["T"] ** 2 / (2 * out["n"]) / np.sqrt(pg)) + 1e-12)
    if with_rich:
        b2 = wf.run_filter_trotter(ch, psi, 2 * ks, dt / 2, ph)
        out["rich_trot"] = 2 * wf.dist(b, b2)
    return out


def bound_n(cx, base, eps, nmax=2_000_000):
    tgt, n, lo = np.sqrt(eps), 8, None
    while n <= nmax:
        o, _ = point(cx, base, n)
        if o["feasible"] and o["bound"] <= tgt:
            break
        lo, n = n, int(n * 1.3) + 1
    else:
        return None
    hi, lo = n, lo or 4
    while hi - lo > max(1, int(0.02 * hi)):
        mid = (lo + hi) // 2
        o, _ = point(cx, base, mid)
        if o["feasible"] and o["bound"] <= tgt:
            hi = mid
        else:
            lo = mid
    return hi


def measured_n(cx, base, eps, nmax):
    tgt, n = np.sqrt(eps), 2
    while n <= nmax:
        o, pk = point(cx, base, n)
        if o["feasible"] and measure(cx, o, pk)["meas"] <= tgt:
            return n
        n = int(n * 1.25) + 1
    return None


def hybrid_n(cx, base, eps, nmax):
    """Non-oracle rule: smallest n on a geometric grid with  sqrt(l_bar) + 2 D(b_n, b_2n) <= sqrt(eps),
    confirmed at the next grid point. Never touches the exact ground state."""
    tgt, n, hit = np.sqrt(eps), 4, None
    hist = []
    while n <= nmax:
        o, pk = point(cx, base, n)
        if o["feasible"]:
            ks, dt, ph, _ = pk
            b = wf.run_filter_trotter(cx.ch, cx.psi, ks, dt, ph)
            b2 = wf.run_filter_trotter(cx.ch, cx.psi, 2 * ks, dt / 2, ph)
            crit = o["bound_leak"] + 2 * wf.dist(b, b2)
            hist.append((n, crit))
            if crit <= tgt:
                if hit is not None:
                    return hit
                hit = n
            else:
                hit = None
        n = int(n * 1.25) + 1
    return hit


def hybrid_record(cx, base, nh):
    o, pk = point(cx, base, nh)
    return measure(cx, o, pk)


def run_eps(cx, trial, eps, steps, sweep=True, nsweep=9):
    ch, gamma = cx.ch, cx.gamma_est
    case = dict(eps=eps, D_trial=float(np.sqrt(1 - cx.gamma_true)))
    if 1 - cx.gamma_est <= eps:
        case["note"] = "trial state already within eps: no filter needed"
        case["cx_total_bound"] = trial["cx"]
        return case
    t0 = time.time()
    base = design(gamma, cx.delta_est, eps)
    case.update(m=int(base["m"]), x=float(base["x"]), T_design=float(base["T"]), eta_target=float(base["eta_target"]),
                eta_cert=float(base["eta_cert"]), f0sq_design=float(base["f0sq"]),
                times=[float(v) for v in base["times"]], phases=[float(v) for v in base["phases"]])
    nb = bound_n(cx, base, eps)
    case["n_bound"] = nb
    if nb is None:
        return case
    o, pk = point(cx, base, nb)
    case["bound_point"] = o
    case["cx_filter_bound"] = nb * steps["a2a"]
    case["cx_total_bound"] = nb * steps["a2a"] + trial["cx"]
    cheap = cx.cost(nb) <= SIM_BUDGET
    if cheap:
        case["at_n_bound"] = measure(cx, dict(o), pk, with_rich=nb <= 6000)
    nm_max = nb if cheap else min(nb, 2000)
    nmeas = measured_n(cx, base, eps, nm_max)
    case["n_meas"] = nmeas
    nh = hybrid_n(cx, base, eps, nm_max)
    case["n_hyb"] = nh
    if nh is not None:
        h = hybrid_record(cx, base, nh)
        case["hybrid"] = {k: h[k] for k in ("n", "bound_leak", "meas", "meas_leak", "meas_trot", "psucc", "eta") if k in h}
        case["cx_filter_hyb"] = nh * steps["a2a"]
        case["cx_total_hyb"] = nh * steps["a2a"] + trial["cx"]
    if nmeas is not None:
        case["cx_filter_meas"] = nmeas * steps["a2a"]
        case["cx_total_meas"] = nmeas * steps["a2a"] + trial["cx"]
    if sweep and nmeas is not None:
        ns = sorted(set(int(round(v)) for v in np.geomspace(max(8, nmeas // 2), 2 * nb, nsweep)))
        ns = [v for v in ns if v <= 6000] or ns[:3]
        sw = []
        for n in ns:
            o2, pk2 = point(cx, base, n)
            if o2["feasible"]:
                sw.append(measure(cx, o2, pk2, with_rich=True))
        case["sweep"] = sw
    case["seconds"] = time.time() - t0
    return case


def run_size(N, J2, chi=2, epss=(0.1, 0.01, 0.001)):
    t0 = time.time()
    ref = DI.dmrg_reference(N, J2)
    E0d, E1d, Etd = ref["E0"], ref["E1"], ref["Etop"]
    Wd = Etd - E0d
    tr = DI.trial_dmrg(N, J2, chi)
    ch = wf.Chain(N, J2, E0_used=E0d, W_used=Wd)
    psi = tr["psi"]
    ref_dense = DI.dense(ref["psi0"]); ref_dense = ref_dense / np.linalg.norm(ref_dense)
    gamma_est = float(abs(np.vdot(ref_dense, psi)) ** 2)
    delta_est = (E1d - E0d) / Wd
    cx = Ctx(ch, psi, gamma_est, delta_est)
    ch_alpha = ch.alpha
    s1, s3 = RS.cx_per_step(N, J2, Wd)
    from qiskit.transpiler import CouplingMap
    l1, l3 = RS.cx_per_step(N, J2, Wd, CouplingMap.from_line(N + 1))
    steps = dict(a2a=int(s1), line=float(l3))
    inputs = dict(N=N, J2=J2, chi=chi, E0_dmrg=E0d, E1_dmrg=E1d, Etop_dmrg=Etd, E0_ed=float(ch.E0_exact), E1_ed=float(ch.E1),
                  Etop_ed=float(ch.Etop), W_dmrg=Wd, W_ed=float(ch.W_exact), gap_dmrg=delta_est, gap_true=float(ch.gap),
                  gamma_est=gamma_est, gamma_true=cx.gamma_true, alpha=float(ch_alpha), chi_ref=ref["chi_ref"],
                  dmrg_E_trial=tr["E"], trial_cx=tr["cx"], trial_mps_overlap=tr["mps_overlap"], bonds=len(ch.bonds),
                  cx_step_a2a=steps["a2a"], cx_step_line=steps["line"])
    cases = []
    for eps in epss:
        c = run_eps(cx, tr, eps, steps, sweep=(eps == 0.01 and N <= 12))
        cases.append(c)
        print(f"   N={N} J2={J2} eps={eps}: n_bound={c.get('n_bound')} n_hyb={c.get('n_hyb')} n_meas={c.get('n_meas')}", flush=True)
    inputs["seconds"] = time.time() - t0
    return dict(inputs=inputs, cases=cases)


def job(args):
    N, J2 = args
    f = OUT / f"main_N{N}_J2{J2}.json"
    if f.exists():
        return str(f)
    print(f"[{time.strftime('%H:%M:%S')}] start N={N} J2={J2}", flush=True)
    epss = (0.1, 0.01, 0.001) if N <= 8 else (0.1, 0.01)
    try:
        r = run_size(N, J2, epss=epss)
    except Exception:
        import traceback; traceback.print_exc()
        return None
    f.write_text(json.dumps(r, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(f"[{time.strftime('%H:%M:%S')}] done N={N} J2={J2} ({r['inputs']['seconds']:.0f}s)", flush=True)
    return str(f)


if __name__ == "__main__":
    import multiprocessing as mp
    sizes = [int(a) for a in sys.argv[1:]] or [4, 6, 8, 10, 12, 14, 16]
    jobs = [(N, J2) for N in sorted(sizes) for J2 in (0.0, 0.4)]
    with mp.get_context("spawn").Pool(4) as p:
        list(p.imap_unordered(job, jobs))
