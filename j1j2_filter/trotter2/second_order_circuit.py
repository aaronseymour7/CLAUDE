"""second_order_circuit.py -- build, count and simulate the SECOND-ORDER (Strang) filter circuit for a stored design.

usage: python trotter2/second_order_circuit.py N J2     (eps = 1e-2, 'guaranteed' design of results/data/dmrg_N{N}_J2{J2}_eps0.01.json)

For the stored design (pulse times tg_i, phases) it
  1. keeps the leakage / p_g certificate unchanged and replaces the first-order Trotter error  sum_i alpha t_i^2/(2 k_i)
     by the second-order  sum_i B2 t_i^3 / k_i^2  (trotter2/alpha2.py), choosing integer k_i = ceil(s t_i) with the smallest s whose
     total bound  (sqrt(leak) + eps_T / sqrt(p_g))^2  <= eps;
  2. builds the actual circuit (qiskit SuzukiTrotter order 2, preserve_order) with those k_i, transpiles it to the RZZ basis
     (level 3 and level 1 = 'direct') and counts RZZ / non-Clifford rotations; compares with the first-order stored design;
  3. simulates the post-selected second-order filter against ED: measured 1-F, p_succ, state error vs. the bound, and the
     smallest s that empirically meets eps (second-order oracle); cross-checks the simulation against the qiskit circuit (N <= 8)."""
import sys, os, json, time
from pathlib import Path
import numpy as np, scipy.sparse as sp
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT), str(ROOT.parent)]
os.chdir(ROOT)
import pipeline as P
from alpha2 import alpha2
from core.trotter import total_error_bound, _trotter_run
from core.circuits import build_early_abort_circuit, apply_filter_pulse
from core.resources import resource_costs, RZZ_BASIS
from core.simulate import postselected_run, postselected_run_exact, state_metrics
from qiskit.quantum_info import SparsePauliOp

N, J2 = int(sys.argv[1]), float(sys.argv[2]); EPS = 1e-2
name = ROOT / "results/data" / f"dmrg_N{N}_J2{J2:g}_eps0.01.json"
d = json.load(open(name)); g = d["designs"]["guaranteed"]; e = d["designs"]["empirical"]
tg, ph, k1 = np.array(g["tg"]), np.array(g["phases"]), np.array(g["k"])
t0 = time.time(); ctx = P.build_ctx(N, 1.0, J2); Hs = ctx["H_scaled"]
B2 = alpha2(Hs); print(f"B2 = {B2:.5f} ({time.time() - t0:.0f}s)", flush=True)

eb = total_error_bound(ctx["gamma"], ctx["gap"], ctx["alpha"], tg, ph, k1, hi=1.0, eta=g["eta"])
leak, p_g = eb["leak"], eb["p_g_lb"]
bound2 = lambda k: (np.sqrt(leak) + np.sum(B2 * tg ** 3 / np.asarray(k, float) ** 2) / np.sqrt(p_g)) ** 2
grid = lambda s: np.maximum(1, np.ceil(s * tg)).astype(int)
lo, hi = 1e-3, 1e4
for _ in range(80):
    mid = np.sqrt(lo * hi)
    lo, hi = (mid, hi) if bound2(grid(mid)) > EPS else (lo, mid)
k2 = grid(hi); n2 = int(k2.sum()); print(f"first order: n={g['n']} k={list(k1)}\nsecond order: n={n2} k={list(k2)}  bound={bound2(k2):.5f}", flush=True)

# ---- circuits and RZZ counts
def counts(k, order, level):
    qc, _ = build_early_abort_circuit(Hs, tg, ph, trial_prep=None, trotter_steps=[int(x) for x in k], order=order)
    c = resource_costs(qc, flatten=True, verbose=False, basis=RZZ_BASIS, optimization_level=level)
    return {x: c[x] for x in ("rzz", "rzz_nonclifford", "rz_nonclifford", "nc_total", "depth")}
terms = d["per_step"]["rzz"]["a2a"]["one"]["rzz"] / 2
rec = dict(N=N, J2=J2, eps=EPS, B2=B2, k1=[int(x) for x in k1], n1=int(g["n"]), k2=[int(x) for x in k2], n2=n2, bound2=float(bound2(k2)),
           first_order_rzz_level3=2 * terms * g["n"], first_order_rzz_direct=4 * terms * g["n"])
for order, lev, key in ((2, 3, "second_l3"), (2, 1, "second_direct")):
    t = time.time(); rec[key] = counts(k2, order, lev); print(key, rec[key], f"({time.time() - t:.0f}s)", flush=True)
rec["rzz_ratio_l3"] = rec["second_l3"]["rzz"] / rec["first_order_rzz_level3"]
rec["rzz_ratio_direct"] = rec["second_direct"]["rzz"] / rec["first_order_rzz_direct"]
rec["rzz_per_step_second_l3_vs_2x"] = rec["second_l3"]["rzz"] / (2 * 2 * terms * n2)      # < 1 means adjacent half-steps were merged

# ---- second-order simulation (same Pauli order as the circuit: forward half-step then reversed half-step)
terms_m = [(sp.csr_matrix(SparsePauliOp("Z" + lab).to_matrix(sparse=True)), float(np.real(c))) for lab, c in zip(Hs.paulis.to_labels(), Hs.coeffs)]
def run2(psi, k):
    dim = 2 ** N; h = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2); state = np.asarray(psi, complex); p_succ = 1.0
    for t_i, phi_i, k_i in zip(tg, ph, k):
        s = np.stack([state, np.zeros(dim, complex)]); rz = np.diag([np.exp(-1j * phi_i), np.exp(1j * phi_i)])
        flat = ((rz @ h) @ s).reshape(-1)
        for _ in range(int(k_i)):
            for seq in (terms_m, terms_m[::-1]):
                for Pm, c in seq:
                    th = c * t_i / k_i / 2
                    flat = np.cos(th) * flat - 1j * np.sin(th) * (Pm @ flat)
        a0 = (h @ flat.reshape(2, dim))[0]; p = float(np.vdot(a0, a0).real); p_succ *= p; state = a0 / np.sqrt(p)
    return state, p_succ
def measure(k):
    sv, p = run2(ctx["trial_vec"], k); p_ex, sv_ex = postselected_run_exact(ctx["trial_qc"], tg, ph, ctx["Hs"])
    _, Fd, Fe = state_metrics(sv, ctx["H_qk"], ctx["psi0_dmrg"], ctx["psi0_ed"])
    ov = abs(np.vdot(sv_ex, sv)) / (np.linalg.norm(sv_ex) * np.linalg.norm(sv))
    return dict(eps_actual=1.0 - (Fe if np.isfinite(Fe) else Fd), p_succ=p, p_succ_exact=p_ex, state_err=float(np.sqrt(max(0.0, 1 - ov ** 2))),
                state_err_bound=float(np.sum(B2 * tg ** 3 / k ** 2) / np.sqrt(p_ex)), sv=sv)
m = measure(k2); sv2 = m.pop("sv"); rec["sim_second"] = m
rec["bound_holds"] = bool(m["eps_actual"] <= rec["bound2"] + 1e-12 and m["state_err"] <= m["state_err_bound"] + 1e-9)
if N <= 8:
    _, svc = postselected_run(ctx["trial_qc"], tg, ph, lambda qc, s, a, t, phi, i: apply_filter_pulse(qc, s, a, Hs, t, phi, int(k2[i]), order=2))
    rec["circ_vs_sim_dev"] = float(1.0 - abs(np.vdot(sv2, svc)) ** 2)
print("sim", {x: (round(v, 6) if isinstance(v, float) else v) for x, v in rec["sim_second"].items()}, "holds", rec["bound_holds"], rec.get("circ_vs_sim_dev"), flush=True)
# oracle: smallest s with measured 1-F <= eps (scan then bisect)
ok = lambda s: measure(grid(s))["eps_actual"] <= EPS
lo, hi = 1e-2, hi
for _ in range(14):
    mid = np.sqrt(lo * hi); lo, hi = (lo, mid) if ok(mid) else (mid, hi)
ke = grid(hi); rec["empirical2"] = dict(k=[int(x) for x in ke], n=int(ke.sum()), **{x: v for x, v in measure(ke).items() if x != "sv"})
rec["empirical2"]["rzz_l3"] = counts(ke, 2, 3)["rzz"]
rec["empirical1"] = dict(n=int(e["n"]), k=e["k"], eps_actual=e["eps_actual"], rzz_l3=2 * terms * e["n"])
rec["wall_s"] = time.time() - t0
json.dump(rec, open(HERE / f"circuit_N{N}_J2{J2:g}.json", "w"), indent=1, default=float)
print("DONE", json.dumps({x: rec[x] for x in ("n1", "n2", "rzz_ratio_l3", "rzz_ratio_direct", "bound_holds")}), "emp1 n", rec["empirical1"]["n"], "emp2 n", rec["empirical2"]["n"], flush=True)
