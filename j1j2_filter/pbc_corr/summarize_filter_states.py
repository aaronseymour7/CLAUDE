"""summarize_filter_states.py -- accuracy of the IDEAL (exact, untrotterised) filter states on periodic chains, from pbc_corr/results/pbc_N*_J2*_eps0.01.json.

Each file holds, for the certified filter design at eps = 1e-2: the exact-filter infidelity, the energy of the filter state, and spin correlators C and block entropies S of the
exact ground state, the trial state and the filter state (C_exact, C_filter, C_trial, S_exact, S_filter, S_trial).  The Trotterised circuit is NOT simulated here.
Prints per case: relative energy error, maximum relative correlator error and maximum relative entropy error of the filter state against the exact ground state.
Cases without the 'C_filter' key are those for which no feasible design was found (the file records an 'error').
usage: python pbc_corr/summarize_filter_states.py"""
import glob, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
rows, missing = [], []
for f in sorted(glob.glob(str(HERE / "results" / "pbc_N*_eps0.01.json")), key=lambda f: (json.load(open(f))["N"], json.load(open(f))["J2"])):
    d = json.load(open(f))
    if "C_filter" not in d:
        missing.append((d["N"], d["J2"], d.get("error")))
        continue
    Ce, Cf = np.array(d["C_exact"]), np.array(d["C_filter"]); Se, Sf = np.array(d["S_exact"]), np.array(d["S_filter"])
    r = (d["N"], d["J2"], d["gamma"], d["infidelity"], abs(d["E_filter"] - d["E0"]) / abs(d["E0"]),
         float(np.max(np.abs(Cf - Ce) / np.maximum(np.abs(Ce), 1e-12))), float(np.max(np.abs(Sf - Se) / np.maximum(Se, 1e-12))))
    rows.append(r)
    print("N=%2d J2=%-6g trial overlap^2 %.3f  1-F(exact filter) %.1e  rel. energy err %.1e  max rel. correlator err %.1e  max rel. entropy err %.1e" % r)
print(f"{len(rows)} cases with a filter state; max relative errors: correlator {max(r[5] for r in rows):.3f}, entropy {max(r[6] for r in rows):.3f}, energy {max(r[4] for r in rows):.1e}")
print("cases without a feasible design:", [(n, j) for n, j, _ in missing])
