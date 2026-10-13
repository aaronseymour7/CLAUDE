"""check_dmrg_status.py -- DMRG convergence (ground / top / excited) for every (N, J2) of the study; writes results/dmrg_status.json.
Needed for the points whose data files were produced before convergence logging existed (N <= 14).  usage: python scripts/check_dmrg_status.py N [N ...]"""
import json, sys, os
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts")); os.chdir(ROOT)
import run_results as RR
import pipeline as P
from hamiltonians import MPO_ham_j1j2
out_f = ROOT / "results" / "dmrg_status.json"
status = json.load(open(out_f)) if out_f.exists() else {}
for N in [int(a) for a in sys.argv[1:]]:
    for J2 in RR.J2_ALL:
        key = f"N{N}_J2{J2:g}"
        if key in status:
            continue
        RR.DMRG_LOG.clear()
        H = MPO_ham_j1j2(N, P.J1, J2)
        P.run_quiet(RR.GE.get_energies, H, N=N, j1=P.J1, j2=J2, run_ed_check=False)
        status[key] = list(RR.DMRG_LOG)
        cur = json.load(open(out_f)) if out_f.exists() else {}
        cur.update({key: status[key]})
        out_f.write_text(json.dumps(cur, indent=1))
        print(key, [(x["role"], x["converged"]) for x in status[key]], flush=True)
