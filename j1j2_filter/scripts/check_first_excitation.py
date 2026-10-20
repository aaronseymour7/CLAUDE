"""check_first_excitation.py -- exact diagonalisation: what is the total spin S of the first excited level of the open J1-J2 chain?

SU(2) symmetry: every eigenvalue of H appears in the S^z = 0 sector once per multiplet with S >= 0, in S^z = 1 once per multiplet with S >= 1,
in S^z = 2 once per multiplet with S >= 2.  Hence, with E0 the lowest S^z=0 level:
  * the lowest excitation is a SINGLET   if the second S^z=0 level lies strictly below the lowest S^z=1 level;
  * it is a TRIPLET (S = 1)              if the lowest S^z=1 level is <= every other S^z=0 level above E0 and below the lowest S^z=2 level;
  * it has S >= 2                        if the lowest S^z=2 level is the lowest excitation.
Also reports whether E0 is non-degenerate and a singlet (E0 absent from the S^z=1 sector).
usage: python scripts/check_first_excitation.py [N ...]      (default 8 10 12 14 16)
"""
import sys
from pathlib import Path
import numpy as np
import scipy.sparse.linalg as sla

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent))
from hamiltonians import j1j2_hamiltonian

TOL = 1e-7
J2S = [round(0.05 * i, 4) for i in range(21)] + [0.2411]
J2S = sorted(set(J2S))


def low_levels(H, idx, k=6):
    Hs = H[idx][:, idx].tocsr()
    d = Hs.shape[0]
    if d <= 1500:
        return np.sort(np.linalg.eigvalsh(Hs.toarray()))[:k]
    return np.sort(sla.eigsh(Hs, k=k, which="SA", return_eigenvectors=False, tol=1e-12))


def scan(N):
    pop = np.array([bin(s).count("1") for s in range(2 ** N)])
    sec = {m: np.where(pop == N // 2 + m)[0] for m in (0, 1, 2)}
    H1 = j1j2_hamiltonian(N, 1.0, 0.0).to_matrix(sparse=True).tocsr()
    H2 = (j1j2_hamiltonian(N, 1.0, 1.0).to_matrix(sparse=True) - H1).tocsr()
    rows = []
    for J2 in J2S:
        H = (H1 + J2 * H2).tocsr()
        l0, l1, l2 = (low_levels(H, sec[m]) for m in (0, 1, 2))
        e0 = l0[0]
        s_next = l0[1]                 # second S^z=0 level
        t1, q1 = l1[0], l2[0]
        e0_singlet = t1 > e0 + TOL     # ground level not repeated in S^z=1
        e0_nondeg = (l0[1] - e0) > TOL
        # first excitation: smallest of {second Sz=0 level, lowest Sz=1 level, lowest Sz=2 level}
        if s_next < t1 - TOL:
            kind = "SINGLET"
        elif t1 <= q1 + TOL and abs(s_next - t1) <= TOL:
            kind = "triplet" if t1 < q1 - TOL else "triplet/quintet-degenerate"
        elif q1 < t1 - TOL:
            kind = "S>=2"
        else:
            kind = "triplet (partner level not found among the lowest 6 Sz=0 levels)"
        rows.append((J2, e0, t1 - e0, s_next - e0, q1 - e0, kind, e0_singlet and e0_nondeg))
    return rows


if __name__ == "__main__":
    Ns = [int(a) for a in sys.argv[1:]] or [8, 10, 12, 14, 16]
    allok = True
    for N in Ns:
        rows = scan(N)
        print(f"\nN = {N}   (gap_t = lowest S^z=1 level - E0;  gap_s2 = second S^z=0 level - E0;  gap_q = lowest S^z=2 level - E0)")
        print(f"{'J2':>7} {'E0':>12} {'gap_t':>10} {'gap_s2':>10} {'gap_q':>10}   first excitation   E0 singlet+nondeg")
        for J2, e0, gt, gs, gq, kind, ok0 in rows:
            print(f"{J2:7.4f} {e0:12.6f} {gt:10.6f} {gs:10.6f} {gq:10.6f}   {kind:<16} {ok0}")
            allok &= kind == "triplet" and ok0
    print("\nALL CASES: first excitation is a triplet and the ground state is a non-degenerate singlet" if allok else "\nNOT ALL CASES TRIPLET -- see rows above")
