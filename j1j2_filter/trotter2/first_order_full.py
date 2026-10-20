"""first_order_full.py -- FULL transpile (RZZ basis, level 3 and direct) of the first-order certified filter circuit of a stored design,
so that first- and second-order counts in the results are measured the same way (not per-step x n).

usage: python trotter2/first_order_full.py N J2     ->  trotter2/first_{N}_{J2}.json
Same circuit builder, basis and conventions as second_order_circuit.py (filter only, early-abort form, no trial circuit)."""
import sys, os, json, time
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT), str(ROOT.parent)]
os.chdir(ROOT)
from hamiltonians import j1j2_hamiltonian, scale_hamiltonian
from core.circuits import build_early_abort_circuit
from core.resources import resource_costs, RZZ_BASIS

N, J2 = int(sys.argv[1]), float(sys.argv[2])
d = json.load(open(ROOT / "results/data" / f"dmrg_N{N}_J2{J2:g}_eps0.01.json")); g = d["designs"]["guaranteed"]
inp = d["inputs"]
Hs = scale_hamiltonian(j1j2_hamiltonian(N, j1=1.0, j2=J2), inp["E0"], inp["W"])      # same scaling as the pipeline (shift E0, certified W)
out = dict(N=N, J2=J2, n=g["n"], k=g["k"])
for lev, key in ((3, "l3"), (1, "direct")):
    t = time.time()
    qc, _ = build_early_abort_circuit(Hs, g["tg"], g["phases"], trial_prep=None, trotter_steps=[int(x) for x in g["k"]], order=1)
    c = resource_costs(qc, flatten=True, verbose=False, basis=RZZ_BASIS, optimization_level=lev)
    out[key] = {x: c[x] for x in ("rzz", "rzz_nonclifford", "rz_nonclifford", "nc_total", "depth")}
    print(N, J2, key, out[key], f"({time.time() - t:.0f}s)", flush=True)
json.dump(out, open(HERE / f"first_{N}_{J2:g}.json", "w"), indent=1, default=float)
