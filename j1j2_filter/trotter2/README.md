# trotter2 -- second-order (Strang) certificate

Replaces the first-order Trotter error alpha*T^2/(2n) in the filter certificate by the commutator-scaling
second-order bound  err <= T^3 B2 / n^2,  B2 = (1/12) sum_a ||[B_a,[B_a,H_a]]|| + (1/24) sum_a ||[H_a,[H_a,B_a]]||,
evaluated exactly for the Pauli sum (`alpha2.py`), max over the two term orders. Everything else in the certificate
(leakage, p_g, trace-distance form) is unchanged: n2 = sqrt(B2 T^3 / (d_T sqrt(p_lb))).

* `validate_alpha2.py`: dense-matrix check, N<=8, t in {3,10}, k in 1..16: measured error <= 0.28 x bound everywhere.
* `compare.py` -> `compare_N*.json`, `table.md`: certified step counts first vs second order on the stored designs
  (first-order recomputation reproduces the stored n to within +1).

Caveats: the 1/12, 1/24 constants and the form of the bound are from memory of Childs-Su-Tran-Wiebe-Zhu, PRX 11, 011020
(arXiv:1912.08854); arxiv.org was unreachable here, so they are NOT checked against the paper. Numerically the bound held
(N<=8) but that is not a proof. alpha2 is O(N^3+) in python, so only N<=12 is computed. RZZ ratio assumes 2 sweeps per
step (no merging of adjacent half-steps; real cost is somewhat lower). The new-circuit gate counts were not built or run.
