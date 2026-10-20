"""true_singlet_gap.py -- periodic chain: lowest excited SINGLET gap vs lowest TRIPLET gap (levels of the S^z=0 sector with no S^z=1 partner are singlets).
Writes results/true_singlet_gap.json.   usage: python pbc_corr/true_singlet_gap.py [J2 ...]"""
import sys, json
from pathlib import Path
import numpy as np, scipy.sparse.linalg as sla
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from pbc_ctx import sector_matrix
J2S = [float(a) for a in sys.argv[1:]] or [0.0, 0.2411, 0.4]
out = {}
for J2 in J2S:
    rows = []
    for N in (8, 10, 12, 14, 16, 18, 20):
        A0, _ = sector_matrix(N, N // 2, 1.0, J2); A1, _ = sector_matrix(N, N // 2 + 1, 1.0, J2)
        k0 = min(10, A0.shape[0] - 2)
        w0 = np.sort(sla.eigsh(A0, k=k0, which="SA", return_eigenvectors=False, tol=1e-10)); w1 = np.sort(sla.eigsh(A1, k=6, which="SA", return_eigenvectors=False, tol=1e-10))
        sing = [e for e in w0 if np.min(abs(w1 - e)) > 1e-7]
        rows.append(dict(N=N, E0=float(sing[0]), gap_singlet=float(sing[1] - sing[0]), gap_triplet=float(w1[0] - sing[0])))
        print("J2", J2, "N", N, "singlet/triplet = %.4f" % (rows[-1]["gap_singlet"] / rows[-1]["gap_triplet"]), flush=True)
    out[str(J2)] = rows
json.dump(out, open(HERE / "results" / "true_singlet_gap.json", "w"))
