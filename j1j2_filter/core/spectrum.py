"""Spectral inputs: certified window, gap, scaling."""


import numpy as np
from get_energies import get_energies
from hamiltonians import (
    MPO_ham_j1j2,
    energy_variance,
    j1j2_hamiltonian,
    pauli_l1_norm,
    reorder_axes,
    symmetry_resolved_gaps,
)


try:
    from get_energies import ED_CHECK_MAX_N
except ImportError:
    ED_CHECK_MAX_N = 14


ED_DENSE_MAX_N = min(ED_CHECK_MAX_N, 14)   # dense eigh / dense Pauli matrices
SECTOR_MAX_N = 12                          # dense symmetry-resolved ED


def get_spectrum(N, j1=1.0, j2=0.0, source="dmrg", res=None, H_mpo=None,
                 bandwidth="certified", gap_mode="any", use_dense=True):
    """Spectrum for FilterBuilder, rescaled to (E - shift)/W with shift = E0.

    res: output of get_energies (pass it in -> DMRG runs once, A9). If None and
         source='dmrg', get_energies is called here.
    bandwidth='certified': W = (c_I + sum|c_beta|) - E0, an upper bound on
         E_top - E0 independent of DMRG's E_top (A3). The price is a smaller
         scaled gap (longer T). 'dmrg': W = E_top - E0 from DMRG (uncertified).
    gap_mode='any': DMRG first excited state (conservative for trial circuits
         that are not exactly symmetric). 'sector': lowest level in the ground
         state's (S, parity) sector (needs N <= SECTOR_MAX_N).

    The shift is the DMRG E0 (variational, >= E0_exact by a tiny amount at the accuracies used); no slack is
    applied. sqrt(var)/W is returned as e0_dev_indicative for reporting only."""
    if source not in ("dmrg", "ed"):
        raise ValueError("source must be 'dmrg' or 'ed'")
    if bandwidth not in ("certified", "dmrg"):
        raise ValueError("bandwidth must be 'certified' or 'dmrg'")
    H_p = j1j2_hamiltonian(N, j1, j2)
    l1 = pauli_l1_norm(H_p)
    cI = float(sum(np.real(c) for lab, c in zip(H_p.paulis.to_labels(), H_p.coeffs)
                   if set(lab) == {"I"}))

    if source == "ed":
        if N > ED_DENSE_MAX_N:
            raise ValueError(f"source='ed' requires N <= {ED_DENSE_MAX_N}")
        evals = np.linalg.eigvalsh(H_p.to_matrix())
        E0, E1, Etop = float(evals[0]), float(evals[1]), float(evals[-1])
        W = Etop - E0
        return dict(energies=(evals - E0) / W, gap=(E1 - E0) / W,
                    totaltime=np.pi * W / (E1 - E0), shift=E0, W=W, E0=E0,
                    E1=E1, Etop=Etop, l1=l1, source="ed", ed=None,
                    e0_dev_indicative=0.0, var=0.0, gap_mode="any", dmrg_res=None)

    if res is None:
        if H_mpo is None:
            H_mpo = MPO_ham_j1j2(N, j1=j1, j2=j2)
        res = get_energies(H_mpo, N=N, j1=j1, j2=j2)
    E0, E1, Etop = (float(np.real(res[k])) for k in ("E0", "E1", "E_top"))

    W_cert = (cI + l1) - E0
    W = max(W_cert, Etop - E0) if bandwidth == "certified" else Etop - E0
    gap_any = E1 - E0
    gap_sector = float("nan")
    if N <= SECTOR_MAX_N:
        sym = symmetry_resolved_gaps(H_p)
        gap_sector = sym["gap_sector"]
        if abs(sym["gap_any"] - gap_any) > 1e-6:
            print(f"  [warning] DMRG gap {gap_any:.8f} != exact lowest gap "
                  f"{sym['gap_any']:.8f} (excited-state DMRG not certified)")
    if gap_mode == "sector":
        if not np.isfinite(gap_sector):
            raise ValueError("gap_mode='sector' needs N <= SECTOR_MAX_N and a "
                             "level in the ground-state sector")
        gap_raw = gap_sector
    else:
        gap_raw = gap_any

    var = 0.0
    if res.get("psi0") is not None and use_dense:
        v = reorder_axes(np.asarray(res["psi0"].to_dense()).reshape(-1))
        var = max(energy_variance(H_p.to_matrix(sparse=True), v), 0.0)
    elif res.get("psi0") is not None:           # MPS/MPO only: Var = <H^2> - E^2 = ||H psi||^2 - E^2
        if H_mpo is None:
            H_mpo = MPO_ham_j1j2(N, j1=j1, j2=j2)
        Hpsi = H_mpo.apply(res["psi0"])
        var = max(float(np.real(Hpsi.H @ Hpsi)) - E0 ** 2, 0.0)
    gap = gap_raw / W
    print(f"[spectrum] E0={E0:.8f} E1={E1:.8f} Etop={Etop:.8f}  "
          f"W({bandwidth})={W:.6f} (W_cert={W_cert:.6f})  gap_any={gap_any:.6f} "
          f"gap_sector={gap_sector:.6f}  scaled gap={gap:.6f}  "
          f"sqrt(var)={np.sqrt(var):.2e}")
    return dict(energies=np.array([0.0, gap, 1.0]), gap=gap,
                totaltime=np.pi / gap, shift=E0, W=W, W_cert=W_cert, E0=E0,
                E1=E1, Etop=Etop, l1=l1, gap_any_raw=gap_any,
                gap_sector_raw=gap_sector, gap_mode=gap_mode, source="dmrg",
                ed=res.get("ed"), e0_dev_indicative=float(np.sqrt(var) / W), var=var,
                dmrg_res=res)
