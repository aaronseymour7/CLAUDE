"""plot_low_levels.py -- ED of the open J1-J2 chain: lowest S=0 level, first excited S=0 level, lowest S=1 level vs J2.
S=0 levels = S^z=0 levels with no S^z=1 partner; S=1 min = lowest S^z=1 level.   usage: python scripts/plot_low_levels.py [N ...]"""
import sys, json
from pathlib import Path
import numpy as np, scipy.sparse.linalg as sla
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
import scipy.sparse as sp


def sector_ops(N, nup, pbc=False):
    """H1 (nearest-neighbour) and H2 (next-nearest) of sum S.S, built directly in the S^z=nup-N/2 sector (open chain, or periodic with pbc=True)."""
    states = np.array([s for s in range(2 ** N) if bin(s).count("1") == nup], dtype=np.int64) if N <= 16 else None
    if states is None:   # enumerate combinations of up-spin positions (fast for N=20)
        from itertools import combinations
        states = np.sort(np.array([sum(1 << p for p in c) for c in combinations(range(N), nup)], dtype=np.int64))
    d = len(states); ops = []
    for r in (1, 2):
        diag = np.zeros(d); rows, cols, vals = [], [], []
        for i in range(N if pbc else N - r):
            j = (i + r) % N; bi, bj = (states >> i) & 1, (states >> j) & 1
            diag += np.where(bi == bj, 0.25, -0.25)
            flip = bi != bj; src = np.where(flip)[0]
            tgt = np.searchsorted(states, states[src] ^ ((1 << i) | (1 << j)))
            rows.append(tgt); cols.append(src); vals.append(np.full(len(src), 0.5))
        ops.append((sp.diags(diag) + sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(d, d))).tocsr())
    return ops
PBC = "--pbc" in sys.argv
Ns = [int(a) for a in sys.argv[1:] if not a.startswith("--")] or [8, 12, 16]
SUF = "_pbc" if PBC else ""
coarse = np.arange(0, 1.0001, 0.025); fine = np.concatenate([np.arange(0.55, 0.7501, 0.005), np.arange(0.18, 0.3001, 0.005)] if PBC else [np.arange(0.55, 0.7501, 0.005)])
J2S = sorted(set(np.round(np.concatenate([coarse, fine, [0.2411]]), 4)))
data = {}
for N in Ns:
    A0, B0 = sector_ops(N, N // 2, PBC); A1, B1 = sector_ops(N, N // 2 + 1, PBC)
    rows = []
    for J2 in J2S:
        ev = lambda A, B: np.sort(sla.eigsh(A + J2 * B, k=12, which="SA", return_eigenvectors=False, tol=1e-11))
        l0, l1 = ev(A0, B0), ev(A1, B1)
        s_ = [e for e in l0 if np.min(abs(l1 - e)) > 1e-7]
        rows.append((float(J2), float(s_[0]), float(s_[1]), float(l1[0])))
    data[N] = rows; print("done N", N, flush=True)
out = ROOT / "results"; (out / "figures").mkdir(exist_ok=True)
json.dump(data, open(out / f"low_levels_ed{SUF}.json", "w"))
fig, ax = plt.subplots(3, len(Ns), figsize=(4.6 * len(Ns), 10), squeeze=False)
for c, N in enumerate(Ns):
    J, e0, s1, t = np.array(data[N]).T
    ax[0, c].plot(J, e0, "k-", label="S=0 min ($E_0$)"); ax[0, c].plot(J, s1, "C3--", label="S=0 first excited"); ax[0, c].plot(J, t, "C0-", label="S=1 min")
    ax[0, c].set_title(f"N={N}, " + ("periodic" if PBC else "open") + " chain: energies"); ax[0, c].set_ylabel("E"); ax[0, c].legend(fontsize=8)
    ax[1, c].plot(J, s1 - e0, "C3--", label="S=0 first excited $-E_0$"); ax[1, c].plot(J, t - e0, "C0-", label="S=1 min $-E_0$")
    ax[1, c].set_ylabel("excitation energy"); ax[1, c].set_title("gaps above $E_0$"); ax[1, c].legend(fontsize=8)
    ax[2, c].plot(J, s1 - t, "C2-", label="S=0 first excited $-$ S=1 min"); ax[2, c].axhline(0, color="k", lw=0.8)
    ax[2, c].set_ylabel("level difference"); ax[2, c].set_title("crossing = sign change"); ax[2, c].legend(fontsize=8)
    for a in ax[:, c]:
        a.axvline(0.2411, color="gray", ls=":", lw=1); a.set_xlabel("$J_2/J_1$"); a.grid(alpha=.3)
fig.tight_layout(); fig.savefig(out / "figures" / f"low_levels{SUF}.png", dpi=150); print("saved", out / "figures" / f"low_levels{SUF}.png")
