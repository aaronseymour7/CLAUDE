"""summarize_low_levels.py -- numbers quoted in WRITEUP.md section 6, read from results/low_levels_ed.json (open) and results/low_levels_ed_pbc.json (periodic).

Columns of both files: J2, E(S=0 minimum), E(S=0 first excited), E(S=1 minimum).  Exact diagonalisation, J1 = 1, J2 grid step as in the files.
Prints, for the open chain, the smallest value over the grid of  E(S=0 first excited) - E(S=1 min)  (positive = the triplet is the lowest excitation) and its location;
for the periodic chain, the J2 at which that difference changes sign (linear interpolation between grid points).
usage: python scripts/summarize_low_levels.py"""
import json
from pathlib import Path
import numpy as np
RES = Path(__file__).resolve().parents[1] / "results"
d = json.load(open(RES / "low_levels_ed.json"))
print("open chain: lowest excitation = triplet?  (difference = E(S=0 first excited) - E(S=1 min))")
for N, rows in d.items():
    a = np.array(rows); J = a[:, 0]; diff = a[:, 2] - a[:, 3]
    near = int(np.argmin(np.abs(J - 0.2411)))
    print(f"  N={N}: {len(J)} grid points in [{J[0]:.2f}, {J[-1]:.2f}]; difference > 0 at every grid point: {bool((diff > 0).all())}; "
          f"smallest difference {diff.min():.4f} at J2 = {J[diff.argmin()]:.3f}; at J2 = {J[near]:.3f}: {diff[near]:.3f}")
d = json.load(open(RES / "low_levels_ed_pbc.json"))
print("periodic chain: sign change of E(S=0 first excited) - E(S=1 min)")
for N, rows in d.items():
    a = np.array(rows); J = a[:, 0]; diff = a[:, 2] - a[:, 3]
    idx = np.where(np.diff(np.sign(diff)) != 0)[0]
    cs = [float(J[i] - diff[i] * (J[i + 1] - J[i]) / (diff[i + 1] - diff[i])) for i in idx]
    print(f"  N={N}: crossing at J2 = {[round(c, 4) for c in cs]}   (grid step {J[1] - J[0]:.3f}, linear interpolation)")
