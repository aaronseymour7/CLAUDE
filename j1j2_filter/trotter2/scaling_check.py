"""scaling_check.py -- how the normalised post-selected state error scales with the step count, first vs second order (N=6, stored design).

usage: python trotter2/scaling_check.py N J2   ->  prints trace distance of the Trotterised filter state to the exact filter state for k-scalings s = 1/4 .. 4 of the design's step counts,
and the fitted exponent  (expected: first order ~ 1/n, second order ~ 1/n^2)."""
import sys, os, json
from pathlib import Path
import numpy as np, scipy.sparse as sp
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT), str(ROOT.parent)]
os.chdir(ROOT)
import pipeline as P
from qiskit.quantum_info import SparsePauliOp
N, J2 = int(sys.argv[1]), float(sys.argv[2])
d = json.load(open(ROOT / "results/data" / f"dmrg_N{N}_J2{J2:g}_eps0.01.json")); g = d["designs"]["guaranteed"]
tg, ph, k1 = np.array(g["tg"]), np.array(g["phases"]), np.array(g["k"], int)
k2 = np.array(json.load(open(HERE / f"circuit_N{N}_J2{J2:g}.json"))["k2"], int)
ctx = P.build_ctx(N, 1.0, J2); Hsp = ctx["H_scaled"]; psi = np.asarray(ctx["trial_vec"], complex); psi /= np.linalg.norm(psi); dim = 1 << N
w, V = np.linalg.eigh(ctx["Hs"].toarray()); c = V.conj().T @ psi
F = np.ones_like(w)
for t, p_ in zip(tg, ph):
    F = F * np.cos(w * t + p_)
exact = V @ (c * F)
terms = [(sp.csr_matrix(SparsePauliOp("Z" + lab).to_matrix(sparse=True)), float(np.real(cc))) for lab, cc in zip(Hsp.paulis.to_labels(), Hsp.coeffs)]
hd = np.array([[1, 1], [1, -1]], complex) / np.sqrt(2)
def run(order, k):
    state = psi.copy(); amp = 1.0
    for t_i, phi_i, k_i in zip(tg, ph, k):
        s = np.stack([state, np.zeros(dim, complex)]); flat = ((np.diag([np.exp(-1j * phi_i), np.exp(1j * phi_i)]) @ hd) @ s).reshape(-1)
        for _ in range(int(k_i)):
            for seq in ((terms,) if order == 1 else (terms, terms[::-1])):
                for Pm, cc in seq:
                    th = cc * t_i / k_i / (1 if order == 1 else 2); flat = np.cos(th) * flat - 1j * np.sin(th) * (Pm @ flat)
        a0 = (hd @ flat.reshape(2, dim))[0]; pr = float(np.vdot(a0, a0).real); amp *= np.sqrt(pr); state = a0 / np.sqrt(pr)
    return state * amp
res = {}
for order, k0 in ((1, k1), (2, k2)):
    rows = []
    for s in (0.25, 0.5, 1, 2, 4):
        k = np.maximum(1, np.round(k0 * s)).astype(int); vt = run(order, k)
        td = float(np.sqrt(max(0.0, 1 - abs(np.vdot(vt, exact)) ** 2 / (np.vdot(vt, vt).real * np.vdot(exact, exact).real))))
        rows.append((int(k.sum()), td))
    n = np.log([r[0] for r in rows]); y = np.log([r[1] for r in rows])
    res[order] = dict(rows=rows, exponent=float(np.polyfit(n, y, 1)[0]))
    print(f"order {order}: " + "  ".join(f"n={a}: {b:.2e}" for a, b in rows) + f"   fitted exponent d log(err)/d log(n) = {res[order]['exponent']:.2f}", flush=True)
json.dump(res, open(HERE / f"scaling_check_N{N}_J2{J2:g}.json", "w"), indent=1)
