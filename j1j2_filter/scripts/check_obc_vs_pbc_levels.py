"""check_obc_vs_pbc_levels.py -- why the lowest singlet and triplet cross in periodic but not in open J1-J2 chains.
(1) ratio (first excited S=0 gap)/(S=1 gap) at J2 = 0 and 0.2411 for N = 8..20 (from results/low_levels_ed*.json): CFT predicts 2 for open chains and 1 for periodic chains at J2c.
(2) N=12 quantum numbers (reflection R, momentum k) of the ground state, lowest triplet and first excited singlet, open and periodic.
usage: python scripts/check_obc_vs_pbc_levels.py"""
import sys, json
from itertools import combinations
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / "scripts"))
d, p = (json.load(open(ROOT / "results" / f)) for f in ("low_levels_ed.json", "low_levels_ed_pbc.json"))
print("(1) (S=0 first-excited gap)/(S=1 gap):")
for N in d:
    r = {round(x[0], 4): x for x in d[N]}; q = {round(x[0], 4): x for x in p.get(N, [])}
    print(f"  N={N}: " + " | ".join(f"{tag} J2={J}: {(x[J][2] - x[J][1]) / (x[J][3] - x[J][1]):.3f}" for tag, x in (("OBC", r), ("PBC", q)) for J in (0.0, 0.2411) if J in x))
sys.argv = ["x"]; exec(open(ROOT / "scripts" / "plot_low_levels.py").read().split("coarse =")[0])   # sector_ops


def sector(N, nup):
    return np.sort(np.array([sum(1 << b for b in c) for c in combinations(range(N), nup)], dtype=np.int64))


def perm_op(states, perm):
    idx = {s: i for i, s in enumerate(states)}; P = np.zeros((len(states),) * 2)
    for i, s in enumerate(states):
        P[idx[sum(1 << perm[b] for b in range(len(perm)) if (s >> b) & 1)], i] = 1
    return P


N = 12; s0, s1 = sector(N, 6), sector(N, 7); refl = [N - 1 - i for i in range(N)]; shift = [(i + 1) % N for i in range(N)]
for pbc in (True, False):
    print(f"\n(2) {'PBC' if pbc else 'OBC'} N=12 (R = reflection, k = momentum / pi)")
    A0, B0 = sector_ops(N, 6, pbc); A1, B1 = sector_ops(N, 7, pbc)
    R0, R1, T0, T1 = perm_op(s0, refl), perm_op(s1, refl), perm_op(s0, shift), perm_op(s1, shift)
    for J2 in (0.1, 0.2411, 0.4, 0.6):
        w0, v0 = np.linalg.eigh((A0 + J2 * B0).toarray()); w1, v1 = np.linalg.eigh((A1 + J2 * B1).toarray())
        sg = [i for i in range(len(w0)) if np.min(abs(w1 - w0[i])) > 1e-8][:2]
        qn = lambda v, R, T: f"R={v @ R @ v:+.0f}" + (f", k={np.angle(v.conj() @ T @ v) / np.pi:.2f}" if pbc else "")
        print(f"  J2={J2}: GS {qn(v0[:, sg[0]], R0, T0)} | triplet gap {w1[0] - w0[sg[0]]:.4f} {qn(v1[:, 0], R1, T1)} | 1st exc. singlet gap {w0[sg[1]] - w0[sg[0]]:.4f} {qn(v0[:, sg[1]], R0, T0)}")
