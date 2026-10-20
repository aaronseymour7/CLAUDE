# trotter2 -- second-order (Strang) certificate, circuits and counts

Replaces the first-order Trotter error alpha*T^2/(2n) in the filter certificate by the commutator-scaling
second-order bound  err <= T^3 B2 / n^2,  B2 = (1/12) sum_a ||[B_a,[B_a,H_a]]|| + (1/24) sum_a ||[H_a,[H_a,B_a]]||
= Childs-Su-Tran-Wiebe-Zhu, arXiv:1912.08854v3, Proposition 16, Eq. (152) (first order: Proposition 15, Eq. (145));
evaluated exactly for the Pauli sum (`alpha2.py`), max over the two term orders. Everything else in the certificate
(leakage, p_g, trace-distance form) is unchanged.

* `validate_two_term.py`, `validate_alpha2.py`: numerical checks of the bound / my implementation (two-term: max err/bound 0.998 over 1800 random cases; chains N<=8: <= 0.28).
* `compare.py` -> `compare_N*.json`, `table.md`: certified step counts first vs second order from the stored designs.
* `second_order_circuit.py` (`run_all_circuits.sh`) -> `circuit_N*_J2*.json`, `circuit_table.md`: builds the second-order circuit (qiskit SuzukiTrotter order 2), full transpile to the RZZ basis, simulation against ED, second-order oracle.
* `first_order_full.py` (`run_first_all.sh`) -> `first_N_J2.json`: full transpile of the first-order certified circuit, same builder and basis, so both orders are measured the same way.
* Results: `results/RESULTS.md` Section 10 (built by `scripts/build_results.py`) and `WRITEUP.md` Section 4.

Caveats: alpha2 is O(N^3+) in python, so only N<=12 is computed; pulse designs are the stored first-order designs (not re-optimised); counts are filter only (no trial circuit).
