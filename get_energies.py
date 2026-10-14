"""
DMRG energies/states (E0, E1, E_top) with an ED cross-check.

Changes from the audited version
--------------------------------
* ED Hamiltonian is built from the PAULI form (hamiltonians.j1j2_hamiltonian),
  so it is independent of the MPO used by DMRG (audit A7); the old
  `build_sparse_j1j2` (same MPO, dense) is removed.
* ED_CHECK_MAX_N = 12 for DENSE ED (2^24 dense was impossible); sparse Lanczos
  ED is used up to ED_SPARSE_MAX_N (audit A12).
* ED eigenvectors are returned in MPS site order (axis-reversed), as before.
* Callers should call get_energies ONCE and pass the result around (A9).
"""
import numpy as np
import quimb.tensor as qtn
from scipy.sparse.linalg import eigsh
from qiskit import transpile
from qiskit.quantum_info import Statevector

from hamiltonians import MPO_ham_j1j2, j1j2_hamiltonian, reorder_axes
from mps_to_circuit import mps_to_circuit

# DMRG accuracy settings.  The old values (bond dims <= 64, cutoff 1e-10, energy tol 1e-10, quimb's default
# local-solver tolerance 1e-3, penalty-MPO compression 1e-12) left errors of ~1e-9 in E0/E1/gap and ~1e-10 in
# the ground-state overlap, even for N <= 14.  These settings reach ~1e-13 (energies) / ~1e-14 (overlap) for N <= 14.
BOND_DIMS = [8, 16, 32, 64, 128, 256]
CUTOFF = 1e-14                 # SVD truncation per bond (1D chains: weight discarded ~ energy error)
MAX_SWEEPS = 100
TOL = 1e-12                    # absolute energy change between sweeps
LOCAL_EIG_TOL = 1e-10          # inner Lanczos tolerance (quimb default 1e-3). Much below 1e-10 makes ARPACK fail.
LOCAL_EIG_NCV = 8
PENALTY_COMPRESS_CUTOFF = 1e-15   # compression of H + lambda|psi0><psi0| (was 1e-12, which limited E1 to ~1e-11)
E1_METHOD = "local"            # "local": projector penalty applied inside each local DMRG solve (any N, any symmetry);
                               # "mpo": projector penalty as a chi^2-bond MPO (original; runs out of memory for N >= 16)
RUN_LOG = []                   # one dict per run_dmrg call: converged, energy, max bond, sweeps, local tol used
ED_CHECK_MAX_N = 12      # dense eigh
ED_SPARSE_MAX_N = 22     # Lanczos (low end + high end)


def build_projector_mpo(mps):
    n = mps.L
    arrays = mps.arrays
    new_arrays = []
    for i in range(n):
        A = arrays[i]
        Ac = np.conj(A)
        if i == 0:
            p, r = A.shape
            T = np.einsum('pr,qR->rRpq', A, Ac).reshape(r * r, p, p)
        elif i == n - 1:
            l, p = A.shape
            T = np.einsum('lp,Lq->lLpq', A, Ac).reshape(l * l, p, p)
        else:
            l, p, r = A.shape
            T = np.einsum('lpr,LqR->lLrRpq', A, Ac).reshape(l * l, r * r, p, p)
        new_arrays.append(T)
    return qtn.MatrixProductOperator(new_arrays, shape='lrud')


class OrthoDMRG2(qtn.DMRG2):
    """Two-site DMRG for the lowest state of H + lam * sum_k |phi_k><phi_k|, i.e. the lowest state (of ANY symmetry)
    essentially orthogonal to the given states phi_k.  The projector enters each local two-site problem as the rank-1
    term lam |w><w|, w = phi contracted with the current left/right environments.  This is mathematically the same
    penalty as building H + lam|phi><phi| as an MPO, but never forms the chi^2-bond MPO (cost O(chi^3) per solve)."""

    def __init__(self, ham, ortho_states, penalty, **kw):
        super().__init__(ham, **kw)
        self._pen = float(penalty)
        self._ortho = []
        for psi in ortho_states:
            p = psi.copy()
            p.retag_({p.site_tag(i): f"_O{i}" for i in range(p.L)})
            p.reindex_({p.site_ind(j): self._b.site_ind(j) for j in range(p.L)})  # physical inds coincide with the bra's
            p.add_tag("_ORTHO")
            self._ortho.append(p)

    def sweep(self, direction, canonize=True, verbosity=0, **update_opts):
        # plain DMRG2 sweep; the overlap environments are rebuilt per local solve in form_local_ops
        return super().sweep(direction, canonize=canonize, verbosity=verbosity, **update_opts)

    def _block_overlaps(self, i, lix):
        """w = <bra with sites i,i+1 removed | phi> as a tensor on the open indices lix of the 2-site block."""
        tags = {self._b.site_tag(i), self._b.site_tag(i + 1)}
        rest = [t for t in self._b if not any(g in t.tags for g in tags)]
        ws = []
        for p in self._ortho:
            t = qtn.TensorNetwork(rest + list(p.tensors)) ^ all
            assert set(t.inds) == set(lix), (t.inds, lix)
            ws.append(np.asarray(t.transpose(*lix).data).reshape(-1))
        return ws

    def form_local_ops(self, i, dims, lix, uix):
        Heff, Neff = super().form_local_ops(i, dims, lix, uix)
        ws = self._block_overlaps(i, lix)
        if isinstance(Heff, np.ndarray):
            Heff = Heff + self._pen * sum(np.outer(w, w.conj()) for w in ws)
        else:
            from scipy.sparse.linalg import LinearOperator
            pen, H0 = self._pen, Heff
            Heff = LinearOperator(H0.shape, dtype=complex,
                                  matvec=lambda x: H0 @ x + pen * sum(w * (w.conj() @ x) for w in ws))
        return Heff, Neff


def sz_sector_mpo(H, N, sz=1.0, factor=None):
    """H + lam (S^z_tot - sz)^2 as an MPO (bond dimension of H + 3).  Its ground state is the lowest state with
    S^z_tot = sz.  For the singlet-ground-state J1-J2 chain (even N) the lowest excitation is a triplet, so for sz = 1
    that energy is E1; this avoids the chi^2-bond projector MPO of the penalty method.  lam is set in get_energies."""
    spin = qtn.SpinHam1D(S=1 / 2)
    spin += 1.0, "Z"                                   # S^z on every site
    Sz = spin.build_mpo(N)
    Sz2 = Sz.apply(Sz)
    P = Sz2 + (-2.0 * sz) * Sz + (sz ** 2) * qtn.MPO_identity(N)
    return H + factor * P


def run_dmrg(mpo, bond_dims=BOND_DIMS, cutoff=CUTOFF, tol=TOL,
             max_sweeps=MAX_SWEEPS, verbosity=0, local_eig_tol=LOCAL_EIG_TOL,
             local_eig_ncv=LOCAL_EIG_NCV, ortho_states=None, ortho_penalty=0.0):
    """Two-site DMRG.  If the inner ARPACK solve fails to converge at local_eig_tol, retry with a looser one
    (x100, x1e4) and record which tolerance was used in RUN_LOG."""
    from scipy.sparse.linalg import ArpackNoConvergence
    last = None
    for factor in (1.0, 1e2, 1e4):
        try:
            if ortho_states:
                dmrg = OrthoDMRG2(mpo, ortho_states, ortho_penalty, bond_dims=bond_dims, cutoffs=cutoff)
            else:
                dmrg = qtn.DMRG2(mpo, bond_dims=bond_dims, cutoffs=cutoff)
            dmrg.opts["local_eig_tol"] = local_eig_tol * factor
            dmrg.opts["local_eig_ncv"] = local_eig_ncv
            converged = dmrg.solve(tol=tol, max_sweeps=max_sweeps, verbosity=verbosity)
            break
        except ArpackNoConvergence as e:
            last = e
    else:
        raise last
    if not converged:
        print(f"  [warning] DMRG did not report convergence within {max_sweeps} sweeps")
    RUN_LOG.append(dict(converged=bool(converged), energy=float(np.real(dmrg.energy)),
                        max_bond=int(max(max(t.shape) for t in dmrg.state.arrays)),
                        sweeps=len(dmrg.energies), local_eig_tol=float(local_eig_tol * factor)))
    return dmrg.energy, dmrg.state


def ed_ground_and_excited(N, j1=1.0, j2=0.0, want_states=True):
    """Exact E0, E1, E_top from the PAULI Hamiltonian (independent of the MPO).
    States are returned in MPS site order."""
    if N > ED_SPARSE_MAX_N:
        raise ValueError(f"N={N} > ED_SPARSE_MAX_N={ED_SPARSE_MAX_N}")
    Hs = j1j2_hamiltonian(N, j1=j1, j2=j2).to_matrix(sparse=True)
    if N <= ED_CHECK_MAX_N:
        w, v = np.linalg.eigh(Hs.toarray())
        lo_w, lo_v = w[:3], v[:, :3]
        top_w, top_v = w[-1], v[:, -1]
    else:
        lw, lv = eigsh(Hs, k=3, which="SA")
        o = np.argsort(lw)
        lo_w, lo_v = lw[o], lv[:, o]
        tw, tv = eigsh(Hs, k=1, which="LA")
        top_w, top_v = tw[0], tv[:, 0]
    out = dict(E0_ed=float(lo_w[0]), E1_ed=float(lo_w[1]), E_top_ed=float(top_w))
    if want_states:
        out["psi0_ed"] = reorder_axes(lo_v[:, 0])
        out["psi1_ed"] = reorder_axes(lo_v[:, 1])
        out["psitop_ed"] = reorder_axes(top_v)
    return out


def ed_cross_check(N, j1=1.0, j2=0.0, E0_dmrg=None, E1_dmrg=None, Etop_dmrg=None,
                   rtol=1e-6, atol=1e-8, want_states=True):
    if N > ED_SPARSE_MAX_N:
        raise ValueError(f"N={N} exceeds ED_SPARSE_MAX_N={ED_SPARSE_MAX_N}; ED skipped.")
    print(f"Running ED cross-check (N={N}, dim=2^{N}={2**N})...")
    ed = ed_ground_and_excited(N, j1, j2, want_states=want_states)

    def _cmp(name, dmrg_val, ed_val):
        if dmrg_val is None:
            return None
        diff = abs(dmrg_val - ed_val)
        ok = diff <= atol + rtol * abs(ed_val)
        print(f"  [{'OK' if ok else 'MISMATCH'}] {name}: DMRG={dmrg_val:.10f}  "
              f"ED={ed_val:.10f}  |diff|={diff:.2e}")
        ed[f"{name}_dmrg"], ed[f"{name}_diff"], ed[f"{name}_ok"] = dmrg_val, diff, ok
        return ok

    checks = [_cmp("E0", E0_dmrg, ed["E0_ed"]),
              _cmp("E1", E1_dmrg, ed["E1_ed"]),
              _cmp("E_top", Etop_dmrg, ed["E_top_ed"])]
    ed["all_ok"] = all(c for c in checks if c is not None)
    if not ed["all_ok"]:
        print("  [warning] DMRG vs ED mismatch -- check bond dims / convergence / penalty.")
    return ed


def get_energies(H, N=None, j1=1.0, j2=0.0, safety_factor=1.1, run_ed_check=True):
    """Run once per configuration and pass the dict around (audit A9)."""
    if H is None:
        H = qtn.MPO_ham_heis(8, j=1, cyclic=False)
    if N is None:
        N = H.L

    print("Running DMRG for ground state...")
    E0, psi0 = run_dmrg(H)
    psi0.normalize()

    print("Running DMRG for highest excited state (via -H)...")
    E_top_neg, psi_top = run_dmrg(-1.0 * H)
    E_top = -E_top_neg

    bandwidth = E_top - E0
    if bandwidth <= 0:
        raise RuntimeError(f"Non-positive bandwidth (E_top={E_top:.6g}, E0={E0:.6g})")
    penalty = safety_factor * bandwidth
    print(f"  spectral bandwidth = {bandwidth:.6g}, PENALTY = {penalty:.6g}")


    e1_method = E1_METHOD
    if e1_method == "local":
        print("Running DMRG for first excited state (projector penalty, local rank-1 form)...")
        _, psi1 = run_dmrg(H, ortho_states=[psi0], ortho_penalty=penalty)
    else:
        print("Running DMRG for first excited state (projector penalty, MPO form)...")
        P0 = build_projector_mpo(psi0)
        H_exc = H + penalty * P0
        H_exc.compress(cutoff=PENALTY_COMPRESS_CUTOFF)
        _, psi1 = run_dmrg(H_exc)
    psi1.normalize()
    E1 = np.real(psi1.H @ H.apply(psi1))

    ed_result = None
    if run_ed_check and N <= ED_SPARSE_MAX_N:
        ed_result = ed_cross_check(N, j1, j2, E0_dmrg=E0, E1_dmrg=E1, Etop_dmrg=E_top)
    elif run_ed_check:
        print(f"  [skip] N={N} > ED_SPARSE_MAX_N; no ED cross-check run.")

    return dict(E0=E0, E1=E1, E_top=E_top, psi0=psi0, psi1=psi1, ed=ed_result, e1_method=e1_method)


# ----------------------------------------------------------------------
# Trial-circuit scoring helpers (unchanged logic)
# ----------------------------------------------------------------------
def circuit_fidelity_and_energy(qc, psi0, H, reverse_qubit_order=True,
                                contract_optimize="greedy", psi0_ed=None):
    N = psi0.L
    if qc.num_qubits != N:
        raise ValueError(f"circuit has {qc.num_qubits} qubits, psi0 has {N} sites")
    vec = Statevector(qc).data
    tensor = vec.reshape([2] * N)
    if reverse_qubit_order:
        tensor = np.transpose(tensor, axes=range(N - 1, -1, -1))
    tensor = np.ascontiguousarray(tensor)

    psi0_dense = psi0.to_dense().reshape([2] * N)
    fidelity = float(np.abs(np.vdot(psi0_dense, tensor)) ** 2)

    ket = qtn.Tensor(tensor, inds=H.lower_inds, tags="KET")
    bra = qtn.Tensor(tensor.conj(), inds=H.upper_inds, tags="BRA")
    energy = float(np.real((H & ket & bra).contract(optimize=contract_optimize)))
    out = dict(fidelity=fidelity, energy=energy)
    if psi0_ed is not None:
        out["fidelity_ed"] = float(np.abs(np.vdot(np.asarray(psi0_ed).reshape([2] * N),
                                                  tensor)) ** 2)
    return out


BASIS = ["cx", "rz", "sx", "x"]


def count_gates_optimized(qc, basis_gates=BASIS, optimization_level=3):
    qct = transpile(qc, basis_gates=basis_gates, optimization_level=optimization_level)
    ops = dict(qct.count_ops())
    return dict(cnot=ops.get("cx", 0), rz=ops.get("rz", 0),
                clifford_total=sum(v for k, v in ops.items() if k != "rz"),
                depth=qct.depth(), all_counts=ops)


def resource_and_accuracy_sweep(mps_arrays, psi0, H, layer_range, shape="lpr",
                                psi0_ed=None, E0=None, E0_ed=None, **kwargs):
    if psi0_ed is not None:
        ceil = float(np.abs(np.vdot(np.asarray(psi0_ed).reshape(-1),
                                    psi0.to_dense().reshape(-1))) ** 2)
        print(f"  [baseline] DMRG-vs-ED fidelity ceiling = {ceil:.8f}")
    rows = []
    for L in layer_range:
        qc = mps_to_circuit(mps_arrays, method="approximate", shape=shape,
                            num_layers=L, **kwargs)
        g = count_gates_optimized(qc)
        a = circuit_fidelity_and_energy(qc, psi0, H, psi0_ed=psi0_ed)
        row = dict(num_layers=L, cnot=g["cnot"], rz=g["rz"],
                   fidelity=a["fidelity"], energy=a["energy"])
        line = (f"L={L:2d}  CNOT={g['cnot']:4d}  RZ={g['rz']:4d}  "
                f"fidelity={a['fidelity']:.6f}  E={a['energy']:.6f}")
        if psi0_ed is not None:
            row["fidelity_ed"] = a["fidelity_ed"]
            if E0_ed is not None:
                row["energy_err_ed"] = a["energy"] - E0_ed
            line += f"  fidelity_ed={a['fidelity_ed']:.6f}"
        print(line)
        rows.append(row)
    return rows
