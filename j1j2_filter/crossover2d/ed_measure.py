"""ed_measure.py -- exact data for open W x L J1-J2 clusters: E0, triplet and singlet gaps, exact Trotter constant of the scaled Hamiltonian, and the Schmidt spectra
at every cut of the snake-ordered chain (-> MPS bond dimension chi(eps) needed for infidelity eps; cost of the sequential MPS state-preparation circuit).
usage: python crossover2d/ed_measure.py 1x16 2x6 3x4 4x4 ...  [--j2 0.0]      (WxL)"""
import sys, json, time
from pathlib import Path
import numpy as np, scipy.sparse.linalg as sla
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from geometry import bonds, pauli_hamiltonian, sector_matrix
from hamiltonians import scale_hamiltonian
from core.trotter import alpha_commutator

args = sys.argv[1:]; j2 = 0.0
if "--j2" in args:
    k = args.index("--j2"); j2 = float(args[k + 1]); args = args[:k] + args[k + 2:]
for spec in args:
    W, L = (int(v) for v in spec.split("x")); N = W * L; t0 = time.time()
    f = HERE / "results" / f"ed_{W}x{L}_J2{j2:g}.json"
    if f.exists():
        continue
    A0, st0 = sector_matrix(W, L, N // 2, 1.0, j2); A1, _ = sector_matrix(W, L, N // 2 + 1, 1.0, j2)
    w0, v0 = sla.eigsh(A0, k=min(6, A0.shape[0] - 2), which="SA", tol=1e-10); o = np.argsort(w0); w0, v0 = w0[o], v0[:, o]
    w1 = np.sort(sla.eigsh(A1, k=min(3, A1.shape[0] - 2), which="SA", return_eigenvectors=False, tol=1e-10))
    sing = [e for e in w0 if np.min(abs(w1 - e)) > 1e-6]
    E0 = float(w0[0]); b1, b2 = bonds(W, L); weyl = (len(b1) + j2 * len(b2)) / 4.0; Wd = weyl - E0
    psi = np.zeros(1 << N); psi[st0] = v0[:, 0]
    spectra = []
    for cut in range(1, N):
        s = np.linalg.svd(psi.reshape(1 << cut, 1 << (N - cut)) if cut <= N - cut else psi.reshape(1 << cut, 1 << (N - cut)).T, compute_uv=False) ** 2
        spectra.append(s[: min(len(s), 4096)].tolist())
    H = scale_hamiltonian(pauli_hamiltonian(W, L, 1.0, j2), E0, Wd)
    alpha = max(alpha_commutator(H, False), alpha_commutator(H, True))
    out = dict(W=W, L=L, N=N, J2=j2, E0=E0, Weyl_width=float(Wd), gap_triplet=float(w1[0] - E0), gap_singlet=float(sing[1] - E0) if len(sing) > 1 else None,
               alpha_scaled=float(alpha), terms=len(H) - 1, bonds=[len(b1), len(b2)], spectra=spectra, seconds=time.time() - t0)
    f.write_text(json.dumps(out))
    print(f"{W}x{L} J2={j2:g}: E0/N={E0 / N:.5f} W={Wd:.2f} triplet gap {out['gap_triplet']:.4f} singlet gap {out['gap_singlet']} alpha_scaled {alpha:.3f} terms {out['terms']} ({time.time() - t0:.0f}s)", flush=True)
