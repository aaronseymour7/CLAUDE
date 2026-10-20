"""compare.py -- certified Trotter step count: first order n1 = alpha T^2/(2 sqrt(p_lb) d_T) (stored 'guaranteed' design) vs second order n2 = sqrt(B2 T^3/(sqrt(p_lb) d_T)).
Run per size:  NONLY=10 python trotter2/compare.py"""
import sys, os, json, glob, numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
for p in (HERE, ROOT, ROOT.parent, ROOT / "pbc_corr"): sys.path.insert(0, str(p))
from alpha2 import alpha2
from pbc_ctx import pbc_hamiltonian
from core.trotter import trotter_alpha
from qiskit.quantum_info import SparsePauliOp
cache, out = {}, []
for f in sorted(glob.glob(str(ROOT / "results/data/*.json"))):
    d = json.load(open(f)); N, J2 = d["N"], d["J2"]
    if N != int(os.environ["NONLY"]): continue
    g = d["designs"]["guaranteed"]; W = d["inputs"]["W_cert"]
    if (N, J2) not in cache:
        H = pbc_hamiltonian(N, 1.0, J2, False); Hs = SparsePauliOp(H.paulis, np.real(H.coeffs) / W)
        cache[(N, J2)] = (trotter_alpha(Hs), alpha2(Hs))
    a1, B2 = cache[(N, J2)]; T = g["T_scaled"]; allowed = g["d_T"] * np.sqrt(g["p_succ_lb"])
    n1 = int(np.ceil(a1 * T ** 2 / 2 / allowed)); n2 = int(np.ceil(np.sqrt(B2 * T ** 3 / allowed)))
    r = dict(N=N, J2=J2, eps=d["eps"], T=T, n_stored=g["n"], n1=n1, n2=n2, a1=a1, B2=B2, rzz_ratio=2 * n2 / n1)
    out.append(r); print(json.dumps(r), flush=True)
json.dump(out, open(HERE / f"compare_N{os.environ['NONLY']}.json", "w"), indent=1)
