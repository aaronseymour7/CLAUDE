"""tower_ed.py -- 2D J1-J2 Heisenberg (square lattice, torus) by sector ED: spin-sector structure relevant to the symmetry-projected filter.

For each cluster Lx x Ly and J2 it reports, in units of the Weyl spectral width W = (#bonds weighted by J)/4 - E0 (as in the 1D study):
  * gap to the lowest triplet (what a filter without spin projection must resolve),
  * gap to the lowest excited SINGLET (what remains after projecting onto S = 0, the paper's J^2 stage),
  * the lowest energy of each total-spin sector S = 0..4 (the 'tower of states'): E(S) - E(0) vs S(S+1).
Sector matrices are built directly (no qiskit).   usage: python lattice2d/tower_ed.py LxLy [LxLy ...] [--j2 0 0.5]"""
import sys, json, time
from pathlib import Path
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
HERE = Path(__file__).resolve().parent


def bonds(Lx, Ly):
    idx = lambda x, y: (x % Lx) * Ly + (y % Ly)
    b1, b2 = set(), set()
    for x in range(Lx):
        for y in range(Ly):
            for dx, dy, S in ((1, 0, b1), (0, 1, b1), (1, 1, b2), (1, -1, b2)):
                i, j = idx(x, y), idx(x + dx, y + dy)
                if i != j:
                    S.add((min(i, j), max(i, j)))
    return sorted(b1), sorted(b2 - b1)


def sector(N, nup, b1, b2, j2):
    s = np.arange(1 << N, dtype=np.int64)
    pc = np.zeros(1 << N, dtype=np.int8)
    for q in range(N):
        pc += ((s >> q) & 1).astype(np.int8)
    st = s[pc == nup]; d = len(st); diag = np.zeros(d); R, C, V = [], [], []
    for J, bl in ((1.0, b1), (j2, b2)):
        if J == 0:
            continue
        for i, j in bl:
            bi, bj = (st >> i) & 1, (st >> j) & 1
            diag += J * np.where(bi == bj, 0.25, -0.25)
            src = np.where(bi != bj)[0]; tgt = np.searchsorted(st, st[src] ^ ((1 << i) | (1 << j)))
            R.append(tgt); C.append(src); V.append(np.full(len(src), 0.5 * J))
    return (sp.diags(diag) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(d, d))).tocsr()


def run(Lx, Ly, j2):
    N = Lx * Ly; b1, b2 = bonds(Lx, Ly); t0 = time.time()
    low = {}
    for m in range(0, 5):
        if N // 2 - m < 0:
            break
        A = sector(N, N // 2 + m, b1, b2, j2)
        k = 10 if m <= 1 else 3
        low[m] = np.sort(sla.eigsh(A, k=min(k, A.shape[0] - 2), which="SA", return_eigenvectors=False, tol=1e-9))
    singlets = [e for e in low[0] if np.min(abs(low[1] - e)) > 1e-6]
    E0 = float(singlets[0]); weyl = (len(b1) + j2 * len(b2)) / 4.0; W = weyl - E0
    tower = {m: float(low[m][0] - E0) for m in low}                      # lowest level with S >= m
    out = dict(N=N, Lx=Lx, Ly=Ly, J2=j2, E0=E0, W=float(W), gap_triplet=float(low[1][0] - E0), gap_singlet=float(singlets[1] - E0),
               tower=tower, seconds=time.time() - t0)
    print(f"{Lx}x{Ly} J2={j2:g}: W={W:.3f}  triplet gap/W={out['gap_triplet'] / W:.4f}  singlet gap/W={out['gap_singlet'] / W:.4f}  ratio={out['gap_singlet'] / out['gap_triplet']:.2f}  "
          f"tower E(S)-E0: " + ", ".join(f"S={m}: {v:.4f}" for m, v in tower.items()) + f"  ({out['seconds']:.0f}s)", flush=True)
    return out


if __name__ == "__main__":
    args = sys.argv[1:]; j2s = [0.0]
    if "--j2" in args:
        k = args.index("--j2"); j2s = [float(a) for a in args[k + 1:]]; args = args[:k]
    res = []
    for a in args:
        Lx, Ly = (int(v) for v in a.split("x"))
        for j2 in j2s:
            res.append(run(Lx, Ly, j2))
    json.dump(res, open(HERE / "results" / ("tower_" + "_".join(args) + ".json"), "w"))
