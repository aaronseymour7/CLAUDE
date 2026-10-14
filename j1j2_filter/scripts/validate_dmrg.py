"""
validate_dmrg.py -- accuracy of the DMRG inputs (get_energies.py) beyond what run_results.py records.

  ed       DMRG vs exact diagonalisation (sparse Lanczos; testing only, never used by the N > 14 pipeline):
           |dE0|, |dE1|, |dEtop|, gap relative error, 1 - |<psi0_DMRG|psi0_ED>|^2.
  consist  no exact reference (N too large): two independent runs from different random starts; run-to-run
           differences of E0, E1, E_top, 1 - |<psi0_a|psi0_b>|^2, and the energy variances <H^2> - <H>^2 of the
           ground and first-excited MPS (a state with variance s^2 contains an eigenvalue within s of its energy).

  gamma    the trial-circuit ground-state weight used by the N > 14 path (MPS contraction against the DMRG state, no ED, no
           statevector) against the dense overlap with the exact ground state (testing only).

usage: python scripts/validate_dmrg.py ed 16 18 20   |   gamma 16 18   |   consist 24 28 32
Writes/updates results/dmrg_validation.json.  Set OMP_NUM_THREADS=1.
"""
import contextlib, io, json, os, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); os.chdir(ROOT)
import get_energies as GE                              # noqa: E402
from hamiltonians import MPO_ham_j1j2                  # noqa: E402

OUT = ROOT / "results" / "dmrg_validation.json"
J2_ALL = (0.0, 0.2411, 0.4)


def run(N, J2):
    GE.RUN_LOG.clear()
    t = time.time()
    H = MPO_ham_j1j2(N, 1.0, J2)
    with contextlib.redirect_stdout(io.StringIO()):
        r = GE.get_energies(H, N=N, j1=1.0, j2=J2, run_ed_check=False)
    return H, r, time.time() - t, all(x["converged"] for x in GE.RUN_LOG), [x["max_bond"] for x in GE.RUN_LOG]


def variance(H, psi, E):
    Hp = H.apply(psi)
    return max(float(np.real(Hp.H @ Hp)) - float(E) ** 2, 0.0)


def save(key, rec):
    cur = json.loads(OUT.read_text()) if OUT.exists() else {}
    cur[key] = rec
    OUT.write_text(json.dumps(cur, indent=1))


def main():
    mode, Ns = sys.argv[1], [int(a) for a in sys.argv[2:]]
    for N in Ns:
        for J2 in J2_ALL:
            key = f"{mode}_N{N}_J2{J2:g}"
            if OUT.exists() and key in json.loads(OUT.read_text()):
                continue
            H, r, secs, conv, chi = run(N, J2)
            if mode == "gamma":
                import pipeline as P
                from mps_to_circuit import mps_to_circuit
                from qiskit.quantum_info import Statevector
                from hamiltonians import reorder_axes
                qc = mps_to_circuit(r["psi0"].arrays, method="approximate", shape="lpr", num_layers=P.SWEEP_L)
                ed = GE.ed_ground_and_excited(N, 1.0, J2, want_states=True)
                g_ed = float(abs(np.vdot(ed["psi0_ed"], reorder_axes(Statevector(qc).data))) ** 2)
                g_mps = P.mps_gamma(qc, r["psi0"])
                rec = dict(N=N, J2=J2, gamma_ed_dense=g_ed, gamma_mps_contraction=g_mps, diff=g_mps - g_ed, secs=secs)
            elif mode == "ed":
                ed = GE.ed_ground_and_excited(N, 1.0, J2, want_states=True)
                d = np.asarray(r["psi0"].to_dense()).reshape(-1)
                d = d / np.linalg.norm(d)
                gap = ed["E1_ed"] - ed["E0_ed"]
                rec = dict(N=N, J2=J2, e1_method=r["e1_method"], dE0=abs(r["E0"].real - ed["E0_ed"]), dE1=abs(r["E1"] - ed["E1_ed"]),
                           dEtop=abs(r["E_top"].real - ed["E_top_ed"]),
                           gap_rel=abs((r["E1"] - r["E0"].real) - gap) / gap,
                           one_minus_overlap=float(1 - abs(np.vdot(ed["psi0_ed"], d)) ** 2), converged=conv, max_bond=chi, secs=secs)
            else:
                H2, r2, secs2, conv2, chi2 = run(N, J2)
                rec = dict(N=N, J2=J2, e1_method=r["e1_method"], E0=float(r["E0"].real), E1=float(r["E1"]), Etop=float(r["E_top"].real),
                           dE0=abs(r["E0"].real - r2["E0"].real), dE1=abs(r["E1"] - r2["E1"]), dEtop=abs(r["E_top"].real - r2["E_top"].real),
                           one_minus_overlap=float(1 - abs(r["psi0"].H @ r2["psi0"]) ** 2),
                           var0=variance(H, r["psi0"], r["E0"].real), var0_b=variance(H, r2["psi0"], r2["E0"].real),
                           var1=variance(H, r["psi1"], r["E1"]), var1_b=variance(H, r2["psi1"], r2["E1"]),
                           converged=conv and conv2, max_bond=chi, secs=secs + secs2)
            save(key, rec)
            print(key, {k: (f"{v:.1e}" if isinstance(v, float) and abs(v) < 1e-2 else v) for k, v in rec.items() if k not in ("N", "J2")}, flush=True)


if __name__ == "__main__":
    main()
