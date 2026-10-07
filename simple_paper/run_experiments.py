"""
run_experiments.py -- the empirical side of the paper.

For every case (N, J2, chi, eps):
  1. exact spectrum + ground state, chi-truncated trial state, gamma
  2. filter design with leakage budget eps_l = eps/4 (so sqrt(l) <= sqrt(eps)/2)  [floor.py]
  3. smallest n whose a-priori bound  sqrt(l) + eps_T/sqrt(p_g)  is <= sqrt(eps)   ("bound n")
  4. simulate the Trotterised filter, measure the three trace distances, and sweep n
Usage: python run_experiments.py [main|chi|all]
"""
import sys, json, time, pathlib
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import workflow as wf
FL = wf.FL

OUT = pathlib.Path(__file__).resolve().parent / "results"
XL, ML = (0.6, 1.0, 1.5, 2.0), (4, 6)


def design(ch, gamma, eps):
    rows, _ = FL.floor_vs_precision([eps / 4], gamma, ch.gap, 0.0, x_list=XL, m_list=ML, hi=1.0, verbose=False)
    r = rows[0]
    assert r["feasible"], "no feasible filter"
    return r


def point(ch, psi, gamma, base, n, with_rich=False):
    """Bound and measurement for the base design snapped to n Trotter steps."""
    g = FL.grid_design(base["times"], base["phases"], base["T"], n, ch.gap, base["eta_target"], 0.0, 1.0)
    c = g["cert"]
    ks, dt, ph = g["k"], g["dt"], g["phases"]
    tt = ks * dt
    Db, Dbl, Dbt = wf.trace_distance_bound(gamma, c["eta"], c["f0sq"], ch.alpha, tt.sum(), n)
    out = dict(n=int(n), k=[int(k) for k in ks], T=float(tt.sum()), eta=float(c["eta"]), f0sq=float(c["f0sq"]),
               feasible=bool(g["feasible"]), bound=Db, bound_leak=Dbl, bound_trot=Dbt)
    return out, (ks, dt, ph, tt)


def measure(ch, psi, gamma, out, pack, with_rich=False):
    ks, dt, ph, tt = pack
    a = wf.run_filter_exact(ch, psi, tt, ph)
    b = wf.run_filter_trotter(ch, psi, ks, dt, ph)
    out.update(meas=wf.dist(b, ch.g), meas_leak=wf.dist(a, ch.g), meas_trot=wf.dist(a, b),
               psucc=float(np.linalg.norm(b) ** 2), psucc_exact=float(np.linalg.norm(a) ** 2),
               psucc_lb=float(gamma * out["f0sq"]))
    bn = b / np.linalg.norm(b)
    out["energy_err"] = float(ch.W * np.vdot(bn, ch.Hs @ bn).real)           # <H> - E0
    out["energy_err_bound"] = float(ch.W * out["bound"] ** 2)
    out["triangle_ok"] = bool(out["meas"] <= out["meas_leak"] + out["meas_trot"] + 1e-12)
    out["bound_ok"] = bool(out["meas"] <= out["bound"] + 1e-12)
    out["leak_ok"] = bool(out["meas_leak"] <= out["bound_leak"] + 1e-12)
    out["trot_ok"] = bool(out["meas_trot"] <= out["bound_trot"] + 1e-12)
    if with_rich:                                     # Richardson-type estimate, no access to |E0>
        b2 = wf.run_filter_trotter(ch, psi, 2 * ks, dt / 2, ph)
        out["rich_trot"] = 2 * wf.dist(b, b2)
    return out


def bound_n(ch, gamma, base, eps, nmax=400_000):
    """Smallest n (coarse geometric scan + bisection) with certified bound <= sqrt(eps)."""
    tgt = np.sqrt(eps)
    n, lo = 8, None
    while n <= nmax:
        o, _ = point(ch, None, gamma, base, n)
        if o["feasible"] and o["bound"] <= tgt:
            break
        lo, n = n, int(n * 1.3) + 1
    else:
        return None
    hi = n
    lo = lo or 4
    while hi - lo > max(1, int(0.02 * hi)):
        mid = (lo + hi) // 2
        o, _ = point(ch, None, gamma, base, mid)
        if o["feasible"] and o["bound"] <= tgt:
            hi = mid
        else:
            lo = mid
    return hi


def empirical_n(ch, psi, gamma, base, eps, nstart=2, nmax=None):
    """Oracle: smallest n on the same snapped designs with measured D <= sqrt(eps) (stays below from there on, checked at 3 larger n)."""
    tgt = np.sqrt(eps)
    n = nstart
    while n <= nmax:
        o, pk = point(ch, psi, gamma, base, n)
        if o["feasible"]:
            m = measure(ch, psi, gamma, o, pk)["meas"]
            if m <= tgt:
                return n
        n = int(n * 1.25) + 1
    return None


def run_case(N, J2, chi, eps, sweep=True, nsweep=9, rich=True):
    t0 = time.time()
    ch = wf.Chain(N, J2)
    psi = wf.mps_truncate(ch.g, N, chi)
    gamma = float(abs(np.vdot(ch.g, psi)) ** 2)
    case = dict(N=N, J2=J2, chi=chi, eps=eps, gamma=gamma, E0=ch.E0, E1=ch.E1, Etop=ch.Etop, W=ch.W,
                gap_scaled=ch.gap, alpha=ch.alpha, D_trial=float(np.sqrt(1 - gamma)))
    if 1 - gamma <= eps:
        case.update(note="trial state already within eps: no filter needed")
        return case
    base = design(ch, gamma, eps)
    case.update(m=int(base["m"]), x=float(base["x"]), T_design=float(base["T"]),
                times=[float(v) for v in base["times"]], phases=[float(v) for v in base["phases"]],
                eta_target=float(base["eta_target"]), eta_cert=float(base["eta_cert"]), f0sq_design=float(base["f0sq"]))
    nb = bound_n(ch, gamma, base, eps)
    case["n_bound"] = nb
    if nb is not None:
        o, pk = point(ch, psi, gamma, base, nb)
        case["at_n_bound"] = measure(ch, psi, gamma, o, pk, with_rich=rich and nb <= 6000)
        ne = empirical_n(ch, psi, gamma, base, eps, nmax=nb)
        case["n_emp"] = ne
        if sweep:
            ns = sorted(set(int(round(v)) for v in np.geomspace(max(8, (ne or 8) // 2), 2 * nb, nsweep)))
            ns = [v for v in ns if v <= 6000] or ns[:3]
            sw = []
            for n in ns:
                o, pk = point(ch, psi, gamma, base, n)
                if not o["feasible"]:
                    continue
                sw.append(measure(ch, psi, gamma, o, pk, with_rich=rich))
            case["sweep"] = sw
    case["seconds"] = time.time() - t0
    return case


MAIN = [(N, J2, 2, eps) for N in (4, 6, 8, 10, 12) for J2 in (0.0, 0.4) for eps in (0.1, 0.01)]
CHI = [(8, 0.0, chi, 0.01) for chi in (1, 2, 3, 4)]

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    jobs = (MAIN if which in ("main", "all") else []) + (CHI if which in ("chi", "all") else [])
    for (N, J2, chi, eps) in jobs:
        f = OUT / f"case_N{N}_J2{J2}_chi{chi}_eps{eps}.json"
        if f.exists():
            continue
        print(f"[{time.strftime('%H:%M:%S')}] N={N} J2={J2} chi={chi} eps={eps}", flush=True)
        try:
            r = run_case(N, J2, chi, eps, sweep=(eps == 0.01 and chi == 2))
        except Exception as e:                       # keep going, report at the end
            import traceback; traceback.print_exc()
            r = dict(N=N, J2=J2, chi=chi, eps=eps, error=repr(e))
        f.write_text(json.dumps(r, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
        print("   done", {k: r.get(k) for k in ("gamma", "n_bound", "n_emp", "seconds")}, flush=True)
