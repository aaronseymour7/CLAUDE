"""Regenerate pipeline_walkthrough.ipynb (cells live here so the notebook diff stays reviewable).

    python j1j2_filter/notebooks/make_notebook.py && jupyter nbconvert --to notebook --execute --inplace j1j2_filter/notebooks/pipeline_walkthrough.ipynb
"""
import os

import nbformat as nbf

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip()))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip()))

md(r"""
# J1–J2 filter pipeline, step by step

Four stages, each a **build cell** followed by a **test cell** of `assert`s:

1. **Hamiltonian** – `sbpipe.j1j2_chain` as a Pauli sum
2. **DMRG spectrum** – E0, gap (E1), E_max from `sbpipe.frontends.dmrg_quimb`
3. **Trial circuit** – the bond-dimension-χ DMRG trial state, compiled to a staircase of 2-qubit unitaries (qiskit)
4. **Filter** – certified Stetcu–Baroni cosine filter (`sbpipe.refine`) and its simulation (`sbpipe.validate`)

Every test compares against an independent reference, mostly exact diagonalisation (ED), so keep `N` ≤ 12.
Each stage only depends on the cells above it, so you can rerun a stage and its test on its own after changing the parameters.
""")

code(r"""
import os, sys, time
sys.path.insert(0, os.path.abspath(".."))          # j1j2_filter/, so `import sbpipe` finds j1j2_filter/sbpipe

import numpy as np
import scipy.sparse.linalg as sla

import sbpipe
from sbpipe.frontends.dmrg_quimb import dmrg_inputs, ground_state, pauli_to_mpo, total_sz, mps_circuit_cost

# ---- parameters ---------------------------------------------------------------------------
N         = 10          # sites (= qubits); ED tests need N <= ~12
J1, J2    = 1.0, 0.2411 # J2 = 0.2411: near the critical point of the open chain
CHI_TRIAL = 2           # bond dimension of the trial MPS (2 -> one bond qubit, 2-qubit gates)
BOND_DIMS = (16, 32, 64, 128)
EPS       = 1e-2        # target infidelity 1 - F
SCHEDULE  = "dyadic"    # "uniform" | "dyadic" | "optimized" | "auto"  (optimized takes ~1 min)
ORDER     = 2           # Trotter order (1 or 2)
TOL       = 1e-8        # numerical tolerance used by the tests

def passed(msg):
    print("✓", msg)
""")

# ------------------------------------------------------------------------------------------
md(r"""
## 1. Hamiltonian

$H = J_1 \sum_i \mathbf S_i\cdot\mathbf S_{i+1} + J_2 \sum_i \mathbf S_i\cdot\mathbf S_{i+2}$, $\mathbf S=\boldsymbol\sigma/2$, open chain,
stored as a `PauliSum` (one term group per bond, which the Weyl bound and the Trotter ordering use).
""")
code(r"""
H = sbpipe.j1j2_chain(N, J1=J1, J2=J2)
H_sp = H.to_sparse()                     # 2^N x 2^N, big-endian (site 0 = most significant bit)
weyl_lo, weyl_hi = sbpipe.weyl_bounds(H)

print(f"N = {N}, J1 = {J1}, J2 = {J2}: {len(H)} Pauli terms, {len(set(H.groups))} bonds")
print("first terms:", list(zip(H.labels[:3], H.coeffs[:3])))
print(f"Weyl spectral bounds: [{weyl_lo:.4f}, {weyl_hi:.4f}]")
""")
code(r"""
# --- test 1: Hamiltonian ---
n_bonds = (N - 1) + (N - 2 if J2 != 0 else 0)
assert len(H) == 3 * n_bonds and len(set(H.groups)) == n_bonds, "wrong number of terms / bonds"
assert H.num_qubits == N
passed("term count: XX, YY, ZZ on every J1 and J2 bond")

Hd = H_sp.toarray()
assert np.allclose(Hd, Hd.conj().T, atol=1e-12), "H is not Hermitian"
passed("Hermitian")

# independent construction: sum_bonds J S_i.S_j from Kronecker products
sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2; sz = np.diag([0.5, -0.5])
def site_op(op, i):
    return np.kron(np.kron(np.eye(2 ** i), op), np.eye(2 ** (N - i - 1)))
S = [[site_op(s, i) for s in (sx, sy, sz)] for i in range(N)]
H_ref = sum(J1 * sum(a @ b for a, b in zip(S[i], S[i + 1])) for i in range(N - 1))
H_ref = H_ref + sum(J2 * sum(a @ b for a, b in zip(S[i], S[i + 2])) for i in range(N - 2))
assert np.allclose(Hd, H_ref, atol=1e-12), "PauliSum disagrees with the Kronecker-product Hamiltonian"
passed("matches an independent Kronecker-product construction")

Sz_tot = total_sz(N).to_sparse().toarray()
Sx_tot = sum(S[i][0] for i in range(N))
assert np.allclose(Hd @ Sz_tot, Sz_tot @ Hd, atol=1e-12) and np.allclose(Hd @ Sx_tot, Sx_tot @ Hd, atol=1e-12)
passed("SU(2) symmetric: [H, S^z_tot] = [H, S^x_tot] = 0")

evals_ed = np.linalg.eigvalsh(Hd)        # full ED spectrum, reused by the later tests
assert weyl_lo <= evals_ed[0] + TOL and evals_ed[-1] <= weyl_hi + TOL, "Weyl bounds do not bracket the spectrum"
passed(f"Weyl bounds bracket the ED spectrum [{evals_ed[0]:.4f}, {evals_ed[-1]:.4f}]")
""")

# ------------------------------------------------------------------------------------------
md(r"""
## 2. DMRG energy spectrum: E0, gap, E_max

`dmrg_inputs` runs two-site DMRG (quimb) at increasing bond dimension for E0, and again with a penalty
$\mu (S^z_{tot}-1)^2$ for E1, the lowest triplet. Their error bars are the change between the last two bond
dimensions. It also returns the reference MPS and the trial ansatz used in section 3.
E_max is the top of the spectrum: DMRG on $-H$, compared with the rigorous Weyl bound (which the pipeline uses
when `E_max` is not given).
""")
code(r"""
t0 = time.time()
est, ansatz_dmrg, extras = dmrg_inputs(H, chi_trial=CHI_TRIAL, bond_dims=BOND_DIMS)
E_max_dmrg = -ground_state(pauli_to_mpo(H.shifted_scaled(0.0, -1.0)), bond_dims=BOND_DIMS)["E0"]

print(f"DMRG ({time.time() - t0:.1f} s)   source: {est.source}")
print(f"  E0    = {est.E0:.10f}  +- {est.E0_err:.1e}")
print(f"  E1    = {est.E1:.10f}  +- {est.E1_err:.1e}")
print(f"  gap   = {est.gap:.6f}")
print(f"  E_max = {E_max_dmrg:.6f} (DMRG on -H),  Weyl bound {weyl_hi:.6f}")
print("  E0 convergence (chi, E, max bond):", [(c, round(e, 10), b) for c, e, b in extras["ref_history"]])
""")
code(r"""
# --- test 2: DMRG spectrum vs exact diagonalisation ---
E0_ed, E_top_ed = evals_ed[0], evals_ed[-1]
assert abs(est.E0 - E0_ed) < 1e-7, f"DMRG E0 off by {est.E0 - E0_ed:.2e}"
assert est.E0 >= E0_ed - 1e-9, "DMRG E0 below the exact ground energy (not variational?)"
passed(f"E0 matches ED ({est.E0 - E0_ed:+.1e})")

# E1: lowest level with S^z_tot = 1, from ED restricted to that sector
mags = Sz_tot.diagonal().real
idx1 = np.flatnonzero(np.isclose(mags, 1.0))
E1_sector_ed = np.linalg.eigvalsh(Hd[np.ix_(idx1, idx1)])[0]
assert abs(est.E1 - E1_sector_ed) < 1e-7, f"DMRG E1 off by {est.E1 - E1_sector_ed:.2e}"
passed(f"E1 matches the lowest S^z=1 level from ED ({est.E1 - E1_sector_ed:+.1e})")

# the certificate assumes every other level is >= E1 - E1_err: check the true first excitation
E1_any_ed = evals_ed[evals_ed > E0_ed + 1e-9][0]
assert E1_any_ed >= est.E1 - est.E1_err - 1e-9, "a level below the DMRG E1 exists: the gap is overestimated"
assert np.sum(np.isclose(evals_ed, E0_ed, atol=1e-9)) == 1, "ground state is degenerate"
passed(f"gap {est.gap:.6f} is the true first gap (ED {E1_any_ed - E0_ed:.6f}); non-degenerate ground state")

assert abs(E_max_dmrg - E_top_ed) < 1e-7 and weyl_hi >= E_top_ed - TOL
passed(f"E_max: DMRG on -H matches ED top {E_top_ed:.6f}; Weyl bound {weyl_hi:.6f} is valid")

win = est.window(weyl_hi)
assert 0 <= win["w0"] < win["delta"] < 1
passed(f"scaled window valid: pass band [0, {win['w0']:.1e}], stop band [{win['delta']:.4f}, 1]")
""")

# ------------------------------------------------------------------------------------------
md(r"""
## 3. Trial circuit

The trial state is the variationally optimal MPS at bond dimension χ = `CHI_TRIAL` (from `dmrg_inputs`;
γ = its ground-state weight). An MPS with χ = 2 is exactly a staircase of 2-qubit unitaries: put it in
right-canonical form $B^{[i]}$, then $U_i$ on qubits $(i, i+1)$ maps $|l\rangle|0\rangle \mapsto \sum_{s,r} B^{[i]}_{l s r}|s\rangle|r\rangle$
(the bond index rides on the next qubit). This is the circuit `mps_circuit_cost` counts. The filter in section 4 runs on the
statevector this circuit prepares.
""")
code(r"""
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector


def right_canonical_mps(psi, n, chi_max):
    # sequential SVDs from the right; returns tensors B[i] of shape (Dl, 2, Dr) and the discarded weight
    Bs, rest, disc = [], psi.reshape(-1, 1), 0.0
    for i in range(n - 1, 0, -1):
        M = rest.reshape(2 ** i, -1)                      # (left sites) x (site i, right bond)
        U, s, Vh = np.linalg.svd(M, full_matrices=False)
        k = min(chi_max, int(np.sum(s > 1e-12 * s[0])))
        disc += float(np.sum(s[k:] ** 2))
        Bs.append(Vh[:k].reshape(k, 2, -1))
        rest = U[:, :k] * s[:k]
    Bs.append(rest.reshape(1, 2, -1))
    return Bs[::-1], disc


def complete_unitary(cols, dim):
    # unitary whose columns at positions cols.keys() are the given orthonormal vectors
    U = np.zeros((dim, dim), complex)
    for j, v in cols.items():
        U[:, j] = v
    free = [j for j in range(dim) if j not in cols]
    if free:
        A = np.column_stack([U[:, list(cols)], np.eye(dim)])
        Q, _ = np.linalg.qr(A)
        U[:, free] = Q[:, len(cols):len(cols) + len(free)]
    return U


def mps_staircase_circuit(Bs):
    # qubit i = site i; qiskit is little-endian, so a 2-qubit gate on (i, i+1) gets qargs [i+1, i]
    n = len(Bs)
    qc = QuantumCircuit(n, name="trial")
    for i, B in enumerate(Bs[:-1]):
        Dl, _, Dr = B.shape
        assert Dr <= 2, "bond dimension > 2 needs wider gates"
        Bp = np.zeros((Dl, 2, 2), complex); Bp[:, :, :Dr] = B
        cols = {2 * l: Bp[l].reshape(4) for l in range(Dl)}     # |l>|0>  ->  sum B[l,s,r] |s>|r>
        qc.unitary(complete_unitary(cols, 4), [i + 1, i], label=f"U{i}")
    B = Bs[-1]
    qc.unitary(complete_unitary({l: B[l, :, 0] for l in range(B.shape[0])}, 2), [n - 1], label=f"U{n-1}")
    return qc


psi_trial = ansatz_dmrg.statevector                       # dense MPS, big-endian
Bs, discarded = right_canonical_mps(psi_trial, N, CHI_TRIAL)
trial_qc = mps_staircase_circuit(Bs)
trial_tqc = transpile(trial_qc, basis_gates=["cx", "u"], optimization_level=1)
psi_circ = Statevector(trial_qc).reverse_qargs().data      # back to big-endian

# the ansatz used by the filter is the state the circuit prepares
ansatz = sbpipe.Ansatz(gamma=ansatz_dmrg.gamma, statevector=psi_circ, cost=ansatz_dmrg.cost,
                       name=f"staircase circuit, chi={CHI_TRIAL}", meta=ansatz_dmrg.meta)

print(f"trial: {ansatz_dmrg.name}, gamma (vs DMRG reference) = {ansatz.gamma:.6f}")
print(f"circuit: {trial_qc.size()} unitaries, after transpile {dict(trial_tqc.count_ops())}, depth {trial_tqc.depth()}")
print(f"cost model (mps_circuit_cost): {ansatz.cost}")
print(trial_qc.draw(fold=120) if N <= 8 else "(set N <= 8 to draw the circuit)")
""")
code(r"""
# --- test 3: trial state and circuit ---
assert abs(np.linalg.norm(psi_trial) - 1) < 1e-12
assert discarded < 1e-20, f"trial state is not a bond-dimension-{CHI_TRIAL} MPS (discarded weight {discarded:.1e})"
passed(f"trial state is normalised and exactly a chi={CHI_TRIAL} MPS")

for i, B in enumerate(Bs[1:], 1):     # right-canonical: B B^dagger = 1 on the left bond
    Bm = B.reshape(B.shape[0], -1)
    assert np.allclose(Bm @ Bm.conj().T, np.eye(B.shape[0]), atol=1e-10), f"site {i} not right-canonical"
for inst in trial_qc.data:
    U = inst.operation.to_matrix()
    assert np.allclose(U @ U.conj().T, np.eye(len(U)), atol=1e-10)
passed("MPS tensors right-canonical; every circuit gate unitary")

fid_circ = abs(np.vdot(psi_trial, psi_circ)) ** 2
assert fid_circ > 1 - 1e-10, f"circuit prepares the wrong state (fidelity {fid_circ:.12f})"
passed(f"circuit statevector = DMRG trial MPS (fidelity {fid_circ:.12f})")

n_cx = trial_tqc.count_ops().get("cx", 0)
assert n_cx <= mps_circuit_cost(N, CHI_TRIAL)["cnot"], "circuit uses more CNOTs than the cost model"
passed(f"{n_cx} CNOTs <= cost model {ansatz.cost['cnot']}")

w_ed, v_ed = sla.eigsh(H_sp, k=1, which="SA"); ground_ed = v_ed[:, 0]
gamma_ed = abs(np.vdot(ground_ed, psi_circ)) ** 2
assert abs(gamma_ed - ansatz.gamma) < 1e-6, f"gamma {ansatz.gamma:.8f} vs ED overlap {gamma_ed:.8f}"
assert 0 < gamma_ed < 1
passed(f"gamma = |<E0|trial>|^2 = {gamma_ed:.6f} agrees with the DMRG estimate")

E_trial = float(np.vdot(psi_circ, H_sp @ psi_circ).real)
assert E_trial >= E0_ed - TOL and abs(E_trial - ansatz.meta["E_trial"]) < 1e-7
passed(f"trial energy {E_trial:.6f} is variational (>= E0 {E0_ed:.6f}) and matches DMRG")
""")

# ------------------------------------------------------------------------------------------
md(r"""
## 4. Filter

`sbpipe.refine` scales the spectrum to $H_s = (H - \text{shift})/W$ (E0 at 0, gap edge at δ, top at 1), designs
$F(E)=\prod_i \cos(E t_i+\phi_i)$ with the chosen schedule (`dyadic` = the Stetcu–Baroni–Carlson halving ladder),
picks Trotter steps per pulse, and certifies $1-F \le$ `EPS`. `sbpipe.validate` simulates the exact and the
Trotterised post-selected filter on the trial circuit's statevector.
""")
code(r"""
t0 = time.time()
plan = sbpipe.refine(H, est, ansatz, sbpipe.Config(eps=EPS, schedule=SCHEDULE, order=ORDER))
s = plan.summary()
res = sbpipe.validate(plan, H, ansatz, ground=ground_ed)

print(f"design ({time.time() - t0:.1f} s): {s['schedule']}, m = {s['m']} pulses, T = {s['T_phys']:.2f}/J1, "
      f"eta = {s['eta']:.3g}, F_pass >= {s['f_pass']:.4f}")
print("  times  (scaled):", np.round(plan.design.times, 3))
print("  phases        :", np.round(plan.design.phases, 3))
print(f"certificate: 1-F <= {s['infidelity_bound']:.3g}  (leakage {s['leak']:.2e}, Trotter {s['eps_trotter']:.2e}), "
      f"p_succ >= {s['p_succ_bound']:.3f}")
print(f"resources: {s['res_qubits']} qubits, {s['res_steps']} Trotter steps, {s['res_attempt_cnot']:,} CNOT/attempt")
print(f"simulated: 1-F = {res['infidelity_trotter']:.2e} (trial {1 - res['gamma_true']:.3f}), "
      f"p_succ = {res['p_succ_trotter']:.4f}, E = {res['energy_trotter']:.8f} (E0 = {E0_ed:.8f})")
""")
code(r"""
# --- test 4: filter ---
d = plan.design
Hs_d = plan.H_scaled.to_sparse().toarray()
es, Vs = np.linalg.eigh(Hs_d)
assert abs(es[0]) <= win["w0"] + TOL and es[1] >= win["delta"] - TOL and es[-1] <= 1 + TOL
passed(f"scaled spectrum: E0 -> {es[0]:.1e}, first excitation -> {es[1]:.4f} >= delta, top -> {es[-1]:.4f} <= 1")

E_stop = np.linspace(win["delta"], 1.0, 200_001)
F_stop = np.abs(d.values(E_stop)).max()
assert F_stop <= d.eta * d.f_pass * (1 + 1e-9), "filter exceeds its certified stop-band bound"
assert abs(d.values(np.array([0.0]))[0]) >= d.f_pass - 1e-12
passed(f"certified eta holds on a 2e5-point grid: max|F| on stop band {F_stop:.2e} <= eta F_pass {d.eta * d.f_pass:.2e}")

# exact filter = cos-product applied in the eigenbasis of H_s
from sbpipe.simulate import filter_exact
c = Vs.conj().T @ ansatz.statevector
v_eig = Vs @ (d.values(es) * c)
v_sim = filter_exact(plan.H_scaled.to_sparse(), ansatz.statevector, d.times, d.phases)
assert np.allclose(v_sim, v_eig, atol=1e-8), "filter_exact disagrees with the spectral decomposition"
p_spec = float(np.sum(np.abs(d.values(es) * c) ** 2))
assert abs(p_spec - res["p_succ_exact"]) < 1e-8
passed(f"exact filter matches eigen-decomposition; success probability {p_spec:.6f}")

assert res["bound_ok"], "simulation violates the certificate"
assert res["infidelity_trotter"] <= s["infidelity_bound"] <= EPS * (1 + 1e-9)
assert res["infidelity_trotter"] < 1 - res["gamma_true"], "filter did not improve the fidelity"
passed(f"Trotterised infidelity {res['infidelity_trotter']:.2e} <= certificate {s['infidelity_bound']:.2e} <= eps {EPS}")

assert res["p_succ_trotter"] >= s["p_succ_bound"] - 1e-12
passed(f"success probability {res['p_succ_trotter']:.4f} >= certified bound {s['p_succ_bound']:.4f}")

assert res["trotter_state_dist"] <= res["trotter_dist_bound"]
passed(f"Trotter state error {res['trotter_state_dist']:.1e} <= bound {res['trotter_dist_bound']:.1e}")

assert E0_ed - TOL <= res["energy_trotter"] < E_trial
passed(f"energy after filter {res['energy_trotter']:.8f} in [E0, E_trial): error {res['energy_trotter'] - E0_ed:.1e}")
""")

md(r"""
### Run everything
`Kernel → Restart & Run All`, or from the repo root:
`jupyter nbconvert --to notebook --execute --inplace j1j2_filter/notebooks/pipeline_walkthrough.ipynb`.
Any failed `assert` stops at the stage that broke.
""")

nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = dict(name="python3", display_name="Python 3", language="python")
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pipeline_walkthrough.ipynb")
nbf.write(nb, out)
print("wrote", out)
