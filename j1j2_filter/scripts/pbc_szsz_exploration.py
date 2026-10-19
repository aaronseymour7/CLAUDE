"""pbc_szsz_exploration.py -- exact (ED) periodic-chain SzSz correlations near J2c: which finite-size estimator of J2c do they give?
(a) N*(-1)^(N/2)*<Sz_0 Sz_{N/2}> versus J2 for N = 8..22 (two families N=4n and 4n+2): crossings between same-family sizes.
(b) single-N log-slope estimator (Eggert-type): q(r) = (-1)^r C(r) * chord(r), chord = (N/pi) sin(pi r/N); q^2 is linear in ln(chord) with a slope that
    vanishes where the marginal coupling vanishes (J2c).  J2*(N) = zero of the slope, using even r and odd r separately (r >= 4).
Data are cached in results/pbc_szsz_ed.json.   usage: python scripts/pbc_szsz_exploration.py [--recompute]"""
import sys, json
from itertools import combinations
from pathlib import Path
import numpy as np, scipy.sparse.linalg as sla
from scipy.interpolate import CubicSpline
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; CACHE = ROOT / "results" / "pbc_szsz_ed.json"
RECOMPUTE = "--recompute" in sys.argv
J2S = np.round(np.arange(0.16, 0.3201, 0.01), 4); NS = (8, 10, 12, 14, 16, 18, 20, 22)
if RECOMPUTE or not CACHE.exists():
    sys.argv = ["x"]; exec(open(ROOT / "scripts" / "plot_low_levels.py").read().split("coarse =")[0])   # sector_ops
    data = {}
    for N in NS:
        st = np.sort(np.array([sum(1 << b for b in c) for c in combinations(range(N), N // 2)], dtype=np.int64)); A, B = sector_ops(N, N // 2, True)
        bits = [((st >> i) & 1) - 0.5 for i in range(N)]
        dd = [sum(bits[i] * bits[(i + r) % N] for i in range(N)) / N for r in range(1, N // 2 + 1)]
        C = []
        for J2 in J2S:
            w, v = sla.eigsh(A + J2 * B, k=1, which="SA", tol=1e-12); p = np.abs(v[:, 0]) ** 2
            C.append([float(p @ d) for d in dd])
        data[N] = C; print("N", N, flush=True)
    json.dump({"J2": J2S.tolist(), "C": {str(k): v for k, v in data.items()}}, open(CACHE, "w"))
D = json.load(open(CACHE)); J2S = np.array(D["J2"]); C = {int(k): np.array(v) for k, v in D["C"].items()}
J = np.linspace(J2S[0], J2S[-1], 3201)


def zero_slope(N, parity, rmin=4):
    rs = [r for r in range(rmin, N // 2 + 1) if r % 2 == parity]
    if len(rs) < 3:
        return None
    chord = np.array([N / np.pi * np.sin(np.pi * r / N) for r in rs]); sl = []
    for k in range(len(J2S)):
        q = np.array([(-1) ** r * C[N][k][r - 1] for r in rs]) * chord
        sl.append(np.polyfit(np.log(chord), q ** 2, 1)[0])
    s = CubicSpline(J2S, sl)(J); z = np.where(np.diff(np.sign(s)) != 0)[0]
    return float(J[z[0]]) if len(z) else None


fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for N in NS:
    y = N * (-1) ** (N // 2) * C[N][:, N // 2 - 1]
    ax[0].plot(J2S, y, "-" if N % 4 == 0 else "--", label=f"N={N}")
ax[0].axvline(0.2411, color="gray", ls=":"); ax[0].set_xlabel("$J_2/J_1$"); ax[0].set_ylabel(r"$N(-1)^{N/2}\langle S^z_0S^z_{N/2}\rangle$")
ax[0].set_title("(a) solid: N=4n, dashed: N=4n+2; same-family crossings at 0.277-0.282"); ax[0].legend(fontsize=7, ncol=2); ax[0].title.set_fontsize(8)
for par, lab, mk in ((0, "even r", "o"), (1, "odd r", "s")):
    pts = [(N, zero_slope(N, par)) for N in NS]; pts = [(n, z) for n, z in pts if z is not None]
    ax[1].plot([p[0] for p in pts], [p[1] for p in pts], mk + "-", label=lab)
    print(lab, [(n, round(z, 4)) for n, z in pts])
ax[1].axhline(0.2411, color="gray", ls=":", label="J2c = 0.2411"); ax[1].set_xlabel("N"); ax[1].set_ylabel("$J_2^*(N)$ (zero log-slope)")
ax[1].set_title("(b) Eggert-type log-slope estimator, even vs odd distances", fontsize=8); ax[1].legend(fontsize=8)
for a in ax: a.grid(alpha=.3)
fig.tight_layout(); (ROOT / "results" / "figures").mkdir(exist_ok=True); fig.savefig(ROOT / "results" / "figures" / "pbc_szsz_exploration.png", dpi=150)
