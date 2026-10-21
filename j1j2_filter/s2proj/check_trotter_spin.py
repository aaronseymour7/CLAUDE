"""check_trotter_spin.py -- does the TROTTERISED filter pulse preserve total spin?  (numerical check of the claim in WRITEUP.md section 5)

Starts from the exact ground state (a J = 0 state) of the open J1-J2 chain, applies two filter pulses built with apply_filter_pulse (first or second order, 3 steps per pulse)
on system + ancilla, and prints the weight outside J = 0 in each ancilla branch.  Rounding-level values (~1e-16) mean the circuit, including its Trotter error, preserves J = 0.
This is a check on the ground state as input, at N = 6 and 8, not a general proof; the supporting argument is that every bond's XX, YY, ZZ factors commute with each other.
usage: python s2proj/check_trotter_spin.py"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path[:0] = [str(ROOT), str(ROOT.parent), str(HERE)]
from hamiltonians import j1j2_hamiltonian, scale_hamiltonian
from core.circuits import apply_filter_pulse
from s2_ops import j_weights
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector
for N, J2 in ((6, 0.4), (8, 0.2411)):
    H = j1j2_hamiltonian(N, j1=1.0, j2=J2); w, v = np.linalg.eigh(H.to_matrix()); psi = v[:, 0]
    Hs = scale_hamiltonian(H, w[0], w[-1] - w[0])
    for order in (1, 2):
        qc = QuantumCircuit(N + 1); qc.initialize(np.concatenate([psi, np.zeros(2 ** N)]), range(N + 1))
        for t, ph in ((1.7, 0.1), (0.9, 0.3)):
            apply_filter_pulse(qc, list(range(N)), N, Hs, t, ph, trotter_steps=3, order=order)
        sv = Statevector(qc).data.reshape(2, 2 ** N)
        res = []
        for a in (0, 1):
            nrm = np.linalg.norm(sv[a])
            if nrm > 1e-12:
                Wt = j_weights(sv[a], N); res.append(f"ancilla {a}: branch probability {nrm ** 2:.3f}, weight with J>0 = {sum(x for k, x in Wt.items() if k > 0):.1e}")
        print(f"N={N} J2={J2} order {order}: " + "; ".join(res))
