"""Measured cheapest dyadic filter at one Hamiltonian.

For every target overlap, eps and protocol: the cheapest (K, r, n) whose simulated 1 - F <= eps, where F is the
fidelity of the post-selected state with the exact ground state. A choice only counts if the same (K, r) still
meets eps at the next finer n (no lucky Trotter cancellations). n is raised along cfg.n_grid until a level brings
no cheaper choice.
"""
import numpy as np

from .costs import cnot_per_attempt, cost, reflection_cnots
from .trotter import step_spectrum


def scan(ch, trials, n, cfg):
    """Every dyadic filter (K <= cfg.K_max, r <= cfg.r_max) at n Strang steps per tau, for every trial column.
    Returns {K: (p[r, j], infid[r, j])}, exact for the Trotterised circuit."""
    tau = np.pi / (1 + ch.gap_s)
    phi, project = step_spectrum(ch, tau / n)
    theta = n * phi                                                # phase per tau; pulse k: cos(2^k theta)
    c = project(trials)
    u = project(ch.g0)
    w = np.abs(c) ** 2
    out, B = {}, np.ones_like(theta)
    for K in range(1, cfg.K_max + 1):
        B = B * np.cos(2.0 ** (K - 1) * theta)
        ps, infs, Br = [], [], np.ones_like(theta)
        for _ in range(cfg.r_max):
            Br = Br * B
            p = (Br ** 2) @ w
            amp = np.abs((u * Br) @ c) ** 2
            ps.append(p); infs.append(np.maximum(0.0, 1 - amp / np.maximum(p, 1e-300)))
        out[K] = (np.array(ps), np.array(infs))
    return out


def _select(ch, levels, ok, eps_list, cfg, unverified_last=False):
    best = {}
    last = len(levels) - 1
    tau = np.pi / (1 + ch.gap_s)
    for i, (n, data) in enumerate(levels):
        check = levels[i + 1][1] if i < last else None
        if check is None and not unverified_last:
            continue
        for K, (p, inf) in data.items():
            for jj in range(len(ok)):
                for e in eps_list:
                    good = inf[:, jj] <= e
                    if check is not None:
                        good &= check[K][1][:, jj] <= e
                    hit = np.flatnonzero(good)
                    if not len(hit):
                        continue
                    r = hit[0] + 1                                 # fewest repetitions that work (and verify)
                    c = cnot_per_attempt(ch.N, K, r, n, cfg.bond_cnots)
                    refl = reflection_cnots(ch.N, K, r, cfg.mcz_cnots_per_qubit)
                    for pr in cfg.protocols:
                        cc = cost(c, p[r - 1, jj], pr, refl)
                        key = (jj, e, pr)
                        if key not in best or cc < best[key]["cost"]:
                            best[key] = dict(K=K, r=int(r), n=n, cnot=c, p=float(p[r - 1, jj]),
                                             infid=float(inf[r - 1, jj]), cost=float(cc),
                                             trotter_steps=int(n * r * (2 ** K - 1)),
                                             T_scaled=float(tau * r * (2 ** K - 1)))
    return best


def measure_point(ch, gammas, eps_list, cfg, verbose=False):
    """Returns (rows, eps_knob): rows[j][eps][protocol] = dict(K, r, n, cnot, p, infid, cost, trotter_steps,
    T_scaled) or None if the overlap is out of reach; eps_knob[j] is the staggered-Rz angle for gammas[j]."""
    eps_knob = [ch.eps_for_overlap(g) for g in gammas]
    ok = [j for j, e in enumerate(eps_knob) if e is not None]
    trials = np.stack([ch.trial(eps_knob[j]) for j in ok], axis=1)
    nkeys = len(ok) * len(eps_list) * len(cfg.protocols)
    levels, best = [], {}
    for i, n in enumerate(cfg.n_grid):
        levels.append((n, scan(ch, trials, n, cfg)))
        new_best = _select(ch, levels, ok, eps_list, cfg)
        improved = any(k not in best or new_best[k]["cost"] < best[k]["cost"] for k in new_best)
        level_best = _select(ch, levels, ok, eps_list, cfg, unverified_last=True)
        could_improve = any(k not in new_best or level_best[k]["cost"] < new_best[k]["cost"] for k in level_best)
        best = new_best
        if verbose:
            print(f"    n={n}: {len(best)}/{nkeys} verified", flush=True)
        if i > 0 and len(best) == nkeys and not improved and not could_improve:
            break
    rows = [{e: {pr: None for pr in cfg.protocols} for e in eps_list} for _ in gammas]
    for (jj, e, pr), b in best.items():
        rows[ok[jj]][e][pr] = b
    return rows, eps_knob
