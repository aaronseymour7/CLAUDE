# nqs2d -- can a neural-network quantum state serve the filter pipeline in 2D?

A small feasibility test for the question "where does a neural network fit in the filter paper?". Open W x L J1-J2 clusters (same
geometry and qubit order as `../crossover2d/geometry.py`), checked against exact diagonalisation:

* `train`: complex RBM (Marshall sign, stochastic reconfiguration with exact local energies) in the S^z = 0 sector and, for the gap, in
  S^z = 1; best bond-dimension-chi MPS by truncating the exact ground state (snake order).
* `refcheck`: the RBM as the classical *reference* state: overlap of each MPS trial measured against the RBM vs against the exact state.
* `cost`: every overlap fed into `sbpipe.refine` with exact energies (eps = 1e-2, second order, cheaper of dyadic/optimized). The RBM
  has no known compact circuit, so its row is the filter cost alone plus the loading budget at which it would tie the best MPS-trial plan.

    python nqs2d/nqs_feasibility.py train 4x4 --j2 0.5 --alpha 1 2 4 --iters 800
    python nqs2d/nqs_feasibility.py train 4x5 --j2 0.5 --alpha 1 2 --iters 500
    python nqs2d/nqs_feasibility.py refcheck
    python nqs2d/nqs_feasibility.py cost [--only 4x5_J20]

Run from `j1j2_filter/`. Needs numpy and scipy only. Run times on one core each: 4x4 training about 1.5 h for alpha = 1, 2, 4; 4x5 about
1.5 h for alpha = 1, 2; `cost` about 30-40 min per lattice (the optimized filter designs dominate). Statevectors (`results/*.npy`) are not
committed.

## Results (eps = 1e-2)

| lattice, J2 | best RBM (real params) | energy rel. err | gamma | gap RBM / exact | gap, alpha = 1 |
|--|--|--|--|--|--|
| 4x4, 0.5 | alpha 4 (2208) | 7.4e-4 | 0.9984 | 0.5510 / 0.5497 | +1.9% |
| 4x4, 0   | alpha 4 (2208) | 6.8e-5 | 0.9999 | 0.5025 / 0.5023 | +0.3% |
| 4x5, 0.5 | alpha 2 (1720) | 3.1e-3 | 0.9885 | 0.4371 / 0.4327 | +4.1% |
| 4x5, 0   | alpha 2 (1720) | 5.3e-4 | 0.9988 | 0.4281 / 0.4275 | -0.1% |

Expected CNOTs per success, J2 = 0.5: MPS chi = 16 + filter 2.5e5 (4x4) and 6.8e5 (4x5); chi = 32 alone meets eps with 2.0e4 and 2.8e4.
The best RBM needs no filter on 4x4 but would have to be loaded in < 2.0e4 CNOTs to tie; on 4x5 it still needs a filter and loses to
chi = 32 even if loading were free (3.4e5). Free trial on 4x4, J2 = 0.5: 5.5e6 / 1.4e6 / 2.8e5 CNOTs per success at gamma = 0.2 / 0.5 / 0.95.

* The variational gap is overestimated in 3 of 4 cases for alpha = 1 (the dangerous direction for the certificate).
* Overlaps of MPS trials measured against the RBM match the exact ones to within 0.01 for gamma <= 0.95 (worst 0.793 vs 0.783, 4x5 J2 = 0.5).
* These clusters are easy (Marshall sign carries 99.8% of the weight at J2 = 0.5; chi = 32 is near exact): the test checks the plumbing
  and error directions, not a network advantage, which could only appear beyond exact diagonalisation.
