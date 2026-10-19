"""analyze.py -- summary of results/s2_*.json: is the S^2 stage useful?  Columns: gap gain, p_S, RZZ without (A) and with (B) the S^2 stage, B/A,
share of B spent on S^2, and B0 = cost of the energy stage alone with the B design and the FULL error allowance (an S^2 stage that costs nothing): B0/A is the
best any implementation of the projection could achieve with this design."""
import json, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
R = [json.load(open(f)) for f in sorted(glob.glob(str(HERE / "results" / "s2_*_N*_eps*.json")))]
print(f"{'bc':3s} {'N':>3s} {'J2':>6s} {'gamma':>6s} {'p_S':>5s} {'D_s/D_t':>7s} | {'A: RZZ':>9s} {'B: RZZ':>9s} {'B/A':>6s} {'S2 share':>8s} | {'B0/A (free S2)':>14s}")
rows = []
for r in sorted(R, key=lambda r: (not r["pbc"], r["J2"], r["N"])):
    A, B = r["A"], r["B"]; g = r["gap_singlet"] / r["gap_triplet"]
    if not A or not B:
        print(f"{'pbc' if r['pbc'] else 'obc':3s} {r['N']:>3d} {r['J2']:>6g} {r['gamma']:6.3f} {r['p_S']:5.3f} {g:7.2f} | {'infeasible' if not A else '%.2e' % A['rzz']:>9s} {'infeasible' if not B else '%.2e' % B['rzz']:>9s}"); continue
    tS = np.array(r["S2_times_over_pi"]) * np.pi
    aE = r["alpha_E"] * B["T"] ** 2 / 2; sE = np.sqrt(r["rzz_per_step_E"] * aE)
    B0 = sE ** 2 / B["allowed"]
    rows.append((r, B["rzz"] / A["rzz"], B0 / A["rzz"]))
    print(f"{'pbc' if r['pbc'] else 'obc':3s} {r['N']:>3d} {r['J2']:>6g} {r['gamma']:6.3f} {r['p_S']:5.3f} {g:7.2f} | {A['rzz']:9.2e} {B['rzz']:9.2e} {B['rzz'] / A['rzz']:6.1f} {B['rzz_S'] / B['rzz']:8.2f} | {B0 / A['rzz']:14.2f}")
if rows:
    print(f"\nB/A over {len(rows)} feasible cases: min {min(x[1] for x in rows):.2f}, median {np.median([x[1] for x in rows]):.2f}, max {max(x[1] for x in rows):.2f};  "
          f"B0/A (free S^2 stage): min {min(x[2] for x in rows):.2f}, median {np.median([x[2] for x in rows]):.2f}, max {max(x[2] for x in rows):.2f}")
