# ======================================================================
# v5: minimax cost (worst-case leakage), pasted after the FilterBuilder cell
# and the setup cell (needs: energies, totaltime, overlap, FilterBuilder).
#
#   v4:  C = |1 - 1/sqrt(sum_n prod_i cos^2(t_i E_n))|  -> mean of f^2 over the band
#   v5:  C = max_{E in [gap, 1]} prod_i cos^2(t_i E)      -> worst leaked level
#
# Same filter (times only, zero phases), same constraints as v4:
# sum t_i = T, 0 <= t_i <= T/3 (set cap_frac=1.0 to drop the cap, c > 1 for T = c*totaltime).
# The max is minimised through its smooth soft-max
#   L_beta = (1/beta) log mean_k exp(beta log f(E_k)^2),
# beta = 1 is log(mean f^2) (= v4), beta -> inf is log(max f^2) (= v5).
# ======================================================================
import numpy as np
from scipy import optimize as opt
import matplotlib.pyplot as plt


def softmax_logf2(times, E, beta):
    """(1/beta) log mean exp(beta log f(E)^2) and its gradient in the times."""
    x = np.outer(E, times)
    l = np.log(np.maximum(np.cos(x) ** 2, 1e-300)).sum(1)
    a = beta * l
    amax = a.max()
    w = np.exp(a - amax)
    Z = w.sum()
    w /= Z
    val = (amax + np.log(Z) - np.log(len(E))) / beta
    grad = (w[:, None] * (-2.0 * E[:, None] * np.tan(x))).sum(0)
    return val, grad


def new_func_v5(timesphases, energies, overlap):
    """v5 cost in the notebook's signature: worst leaked level max_{E>0} prod cos^2(t_i E) (phases ignored)."""
    n = len(timesphases) // 2
    E = energies[energies > 0]
    return float(np.prod(np.cos(np.outer(E, timesphases[:n])) ** 2, axis=1).max())


def build_v5(builder, c=1.0, cap_frac=1.0 / 3.0, n_starts=24, seed=0,
             betas=(2, 8, 32, 128, 512, 2048), verbose=True, label="v5"):
    """Minimax design for each ntimes in builder.a..builder.b. Returns the same list of dicts as builder.build().
    betas=(1,) instead gives the v4 objective (log mean f^2 on the notebook's levels) with this multistart
    optimizer, e.g. v4 without the cap: build_v5(builder, cap_frac=1.0, betas=(1,), label="v4")."""
    T = c * builder.total_time
    cap = cap_frac * T
    gap = builder.energies[builder.energies > 0].min()
    if tuple(betas) == (1,):
        E = builder.energies[builder.energies > 0]          # v4: exactly the notebook's levels
    else:
        E = np.linspace(gap, builder.energies.max(), max(len(builder.energies), int(10 * T)))
    rng = np.random.default_rng(seed)
    results = []
    for m in range(builder.a, builder.b + 1):
        if m * cap < T:
            continue
        starts = [1 / 2.0 ** np.arange(m), np.ones(m)] + [rng.dirichlet(np.ones(m)) for _ in range(n_starts)]
        cons = [{"type": "eq", "fun": lambda t: t.sum() - T, "jac": lambda t: np.ones_like(t)}]
        best = None
        for s in starts:
            t = np.clip(T * s / s.sum(), 0, cap)
            t *= T / t.sum()
            for b in betas:
                r = opt.minimize(softmax_logf2, t, args=(E, b), jac=True, method="SLSQP",
                                 bounds=[(0, cap)] * m, constraints=cons, options=dict(maxiter=500, ftol=1e-12))
                t = np.clip(r.x, 0, cap)
            v = softmax_logf2(t, E, betas[-1])[0]
            if best is None or v < best[0]:
                best = (v, t)
        t = np.sort(best[1])
        fun = new_func_v5(np.r_[t, np.zeros(m)], builder.energies, builder.overlap)
        results.append({"ntimes": m, "times": t.copy(), "phases": np.zeros(m), "fun": fun, "success": True,
                        "message": label, "result": None})
        if verbose:
            print(f"{label} ntimes={m:2d}  time={t.sum():.6f}  max f^2 above gap={fun:.3e}")
    return results


def leak_stats(times, E):
    f2 = np.prod(np.cos(np.outer(E, times)) ** 2, axis=1)
    return f2.max(), f2.mean()


# ---------------------------------------------------------------- run v4 and v5 on the same setup
builder = FilterBuilder(total_time=totaltime, energies=energies, overlap=overlap, a=4, b=10)
res_v4 = builder.build("v4")
designs = {
    "v4 (authors: mean)": res_v4,
    "v4 mean, same T, no cap": build_v5(builder, cap_frac=1.0, betas=(1,), label="v4 nocap"),
    "v5 minimax, same T and cap": build_v5(builder),                       # authors' T and constraints
    "v5 minimax, same T, no cap": build_v5(builder, cap_frac=1.0),         # T = pi/gap, t_i <= T
    "v5 minimax, T = 1.5 pi/gap": build_v5(builder, c=1.5),                # cap kept, longer total time
}

gap = energies[energies > 0].min()
Eband = np.linspace(gap, 1, 20001)
nonzero = lambda t: np.sort(t[t > 1e-3 * t.sum()])                         # zero-length pulses do nothing
best = {}
for name, res in designs.items():
    crit = (lambda t: leak_stats(t, energies[1:])[1]) if name.startswith("v4") else (lambda t: leak_stats(t, Eband)[0])
    best[name] = nonzero(min(res, key=lambda r: crit(r["times"]))["times"])

# two trials with overlap `overlap`: the notebook's uniform one, and one whose excited weight sits on the
# 20 levels just above the gap (closer to a real low-energy trial state)
N = len(energies)
gs_state = np.zeros(N); gs_state[0] = 1.0
trials = {"uniform": np.full(N, np.sqrt((1 - overlap ** 2) / (N - 1))), "near gap": np.zeros(N)}
trials["near gap"][1:21] = np.sqrt((1 - overlap ** 2) / 20)
for v in trials.values():
    v[0] = overlap

print("\n                                m   T       max f^2>gap   mean f^2    fdiff uniform   fdiff near-gap")
for name, t in best.items():
    mx, mn = leak_stats(t, Eband)
    fd = [np.sum((gs_state - FilterBuilder.apply_filter(t, np.zeros(len(t)), energies, v.copy())[0]) ** 2)
          for v in trials.values()]
    print(f"{name:36s} {len(t):2d}  {t.sum():6.2f}   {mx:10.2e}   {mn:9.2e}    {fd[0]:10.2e}      {fd[1]:10.2e}")

# ---------------------------------------------------------------- plot: filter vs step, linear and log
COL = ["C0", "C2", "C1", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]   # one per design; zip() drops any design without a colour
Eplot = np.linspace(0, 1, 6001)
fig, ax = plt.subplots(1, 2, figsize=(14, 4.8))
for (name, t), col in zip(best.items(), COL):
    f2 = np.prod(np.cos(np.outer(Eplot, t)) ** 2, axis=1)
    ax[0].plot(Eplot, f2, color=col, lw=2, label=f"{name} (m={len(t)}, max>gap {leak_stats(t, Eband)[0]:.1e})")
    ax[1].semilogy(Eplot, np.maximum(f2, 1e-12), color=col, lw=1.4, label=name)
ax[0].plot([0, gap, gap, 1], [1, 1, 0, 0], "k:", lw=1.2, label="ideal step")
for a in ax:
    a.axvline(gap, color="gray", ls="--", lw=0.8)
    a.set_xlabel("Energy (scaled)")
    a.grid(alpha=0.3)
ax[0].set_ylabel(r"$f(E)^2 = \prod_i \cos^2(t_i E)$")
ax[1].set_ylim(1e-8, 2)
ax[0].legend(fontsize=8)
ax[0].set_title(f"gap = {gap:.3f}, authors' T = {totaltime:.3f}")
ax[1].set_title("same, log scale")
plt.tight_layout()
plt.show()

# ---------------------------------------------------------------- what survives of the near-gap trial
fig, ax = plt.subplots(figsize=(10, 4.2))
for (name, t), col in zip(best.items(), COL):
    f0, _ = FilterBuilder.apply_filter(t, np.zeros(len(t)), energies, trials["near gap"].copy())
    ax.semilogy(energies[1:21], f0[1:21] ** 2, "o", ms=5, color=col,
                label=f"{name}: ground weight {f0[0] ** 2:.3f}")
ax.set_xlabel("Energy of the excited level (just above the gap)")
ax.set_ylabel("weight left after filtering")
ax.set_title(f"near-gap trial (overlap {overlap}): leftover excited weight")
ax.legend(fontsize=8); ax.grid(alpha=0.3); plt.tight_layout(); plt.show()
