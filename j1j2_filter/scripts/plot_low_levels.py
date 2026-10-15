"""plot_low_levels.py -- ED of the open J1-J2 chain: lowest S=0 level, first excited S=0 level, lowest S=1 level vs J2.
S=0 levels = S^z=0 levels with no S^z=1 partner; S=1 min = lowest S^z=1 level.   usage: python scripts/plot_low_levels.py [N ...]"""
import sys, json
from pathlib import Path
import numpy as np, scipy.sparse.linalg as sla
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from hamiltonians import j1j2_hamiltonian
Ns = [int(a) for a in sys.argv[1:]] or [8, 12, 16]
J2S = sorted(set(np.round(np.arange(0, 1.0001, 0.025), 4)) | {0.2411})
data = {}
for N in Ns:
    pop = np.array([bin(s).count("1") for s in range(2 ** N)])
    i0, i1 = np.where(pop == N // 2)[0], np.where(pop == N // 2 + 1)[0]
    H1 = j1j2_hamiltonian(N, 1., 0.).to_matrix(sparse=True).tocsr(); H2 = (j1j2_hamiltonian(N, 1., 1.).to_matrix(sparse=True) - H1).tocsr()
    rows = []
    for J2 in J2S:
        H = (H1 + J2 * H2).tocsr()
        ev = lambda idx: np.sort(sla.eigsh(H[idx][:, idx], k=10, which="SA", return_eigenvectors=False, tol=1e-12))
        l0, l1 = ev(i0), ev(i1)
        s = [e for e in l0 if np.min(abs(l1 - e)) > 1e-7]
        rows.append((float(J2), float(s[0]), float(s[1]), float(l1[0])))
    data[N] = rows
out = ROOT / "results"; (out / "figures").mkdir(exist_ok=True)
json.dump(data, open(out / "low_levels_ed.json", "w"))
fig, ax = plt.subplots(2, len(Ns), figsize=(4.6 * len(Ns), 7), squeeze=False)
for c, N in enumerate(Ns):
    J, e0, s1, t = np.array(data[N]).T
    ax[0, c].plot(J, e0, "k-", label="S=0 min ($E_0$)"); ax[0, c].plot(J, s1, "C3--", label="S=0 first excited"); ax[0, c].plot(J, t, "C0-", label="S=1 min")
    ax[0, c].set_title(f"N={N}, open chain: energies"); ax[0, c].set_ylabel("E"); ax[0, c].legend(fontsize=8)
    ax[1, c].plot(J, s1 - e0, "C3--", label="S=0 first excited $-E_0$"); ax[1, c].plot(J, t - e0, "C0-", label="S=1 min $-E_0$")
    ax[1, c].set_ylabel("excitation energy"); ax[1, c].set_title("gaps above $E_0$"); ax[1, c].legend(fontsize=8)
    for a in ax[:, c]:
        a.axvline(0.2411, color="gray", ls=":", lw=1); a.set_xlabel("$J_2/J_1$"); a.grid(alpha=.3)
fig.tight_layout(); fig.savefig(out / "figures" / "low_levels.png", dpi=150); print("saved", out / "figures" / "low_levels.png")
