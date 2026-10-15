# Results: resource costs of the cosine-filter pipeline for the open J1–J2 chain

*Produced by `scripts/run_results.py` + `scripts/build_results.py` from the pipeline in this directory (`pipeline.py`, `floor.py`, `core/*`). Every number below is read from `results/data/*.json`; raw tables are in `results/tables/*.csv`.*

## 0. How to read this document

* **Pipeline.** DMRG (`get_energies`) → spectrum (`get_spectrum`) → L-layer approximate trial circuit (`mps-to-circuit`, L = 1) → filter design (`floor.py`, certified η) → snapped Trotter grid → counts. Settings are those in `pipeline.py`: first-order Lie–Trotter, ε targets as given, `SWEEP_L = 1`, `BANDWIDTH = "certified"` for the DMRG path, `SYNTH_EPS = 1e-3`.
* **Two designs per case.** *Guaranteed* = smallest n whose a-priori rigorous bound on 1−F is ≤ ε. *Empirical* = smallest n whose **measured** 1−F (against the exact ground state) is ≤ ε. The empirical design uses the exact ground state, so it is an oracle lower envelope, not something a practitioner can know in advance; it needs simulation and is therefore only available for N ≤ 14.
* **Everything uses the DMRG spectrum.** E0, E1 (orthogonality-penalty DMRG: lowest state of any symmetry orthogonal to ψ0) and E_top from DMRG; the bandwidth is the **certified** W = Σ_bonds λmax(bond) − E0 (Weyl's inequality: λmax of a sum ≤ sum of λmax; a bond J S·S has λmax = J/4), which equals E_top − E0 for J1, J2 ≥ 0 because the fully polarized state attains it (Section 6). The DMRG gap, shift and γ are not certified: every guarantee below is a *conditional certificate* (exact only if the DMRG reference is exact; see `CERTIFICATION.tex`, Section 1 and Section 12).
* **N ≤ 14 versus N > 14.** For N ≤ 14 an exact-diagonalisation cross-check and statevector simulation of the filter are available: γ is the overlap with the ED ground state, and measured 1−F, P_succ and the state error are reported. For **N > 14 there is no ED, no statevector and no dense or sparse 2^N operator**: γ is the overlap of the DMRG MPS with the trial circuit simulated as an MPS (Aer matrix-product-state), the energy variance comes from MPO contractions, α is the same exact commutator sum as for small N (Table 1b), and circuits are **counted, not simulated** (no measured 1−F or P_succ; only the rigorous bound).
* **Time.** `T (scaled)` is the total evolution time Σtᵢ in units of 1/W; `T (phys)` = T/W in units of 1/J1 (H = Σ J S·S, J1 = 1). Pulses with kᵢ = 0 on the grid are dropped.
* **Counting and the two compilations.** One Trotter step of one pulse (every Pauli term of H, controlled by the ancilla) is transpiled once and multiplied by n; per pulse one Rz(2φ) is added. Two compilations of the same circuit are reported because they trade CX against non-Clifford rotations: **direct** (qiskit `optimization_level=1`, every Pauli rotation becomes a CX ladder around one Rz): CX = 4 × (#Pauli terms) per step and non-Clifford Rz = #Pauli terms + 1 per step; **level 3** (block re-synthesis): about 1.7× fewer CX but about 3× more non-Clifford Rz. `RZ (non-Cliff.)` counts the Rz gates whose angle is not a multiple of π/2 (the cost metric for fault-tolerant hardware; no T-count is reported because the synthesis constants are unverified). Per-step × n was checked against a full transpile of a 48-step circuit: CX exact for both compilations, non-Clifford Rz exact for direct and up to 4% too high for level 3 (Section 9). **Level-3 totals include the trial circuit**; direct counts are for the filter only (the trial circuit is 45–93 CX). Connectivity is all-to-all unless marked *line* (ancilla at the end of a line: pessimistic).
* **Success path.** All counts are for the all-zeros (success) path. The expected cost per accepted run is total ÷ P_succ (early abort makes it lower).

Data present: 66 (N, J2, ε) points, N from 4 to 32; 0 failed point(s).

## Summary

* **N ≤ 14, guaranteed designs, ε = 10⁻²:** 76–22,902 Trotter steps, 1,794–4,007,889 total CX, T (phys) 1.9–13.3 (1/J1).
* **Guaranteed versus empirical (oracle):** the guaranteed design uses 15–559× the steps and 15–556× the CX of the smallest n that empirically meets ε (median 191×).
* **N = 16–32 (DMRG/MPS only, counted not simulated), ε = 10⁻²:** 17,333–224,514 steps, 2.79e+06–9.59e+07 total CX (level 3; direct compilation 4.78e+06–1.64e+08), 1.22e+06–4.13e+07 non-Clifford RZ (direct; level 3: 3.64e+06–1.24e+08), T (phys) 7.9–20.5 (converged-DMRG points only; 0 further point(s) have an unconverged excited-state DMRG and are flagged in Section 7; 0 point(s) designed by the fallback FilterBuilder, marked † there, are excluded from these ranges and from the fits).
* **Scaling, J2 = 0:** guaranteed total CX ~ N^5.1 (R² 0.999), steps ~ N^4.0, CX per step ~ N^1.1 (fits over N = 6–24).
* **Scaling, J2 = 0.2411:** guaranteed total CX ~ N^5.2 (R² 0.998), steps ~ N^4.1, CX per step ~ N^1.2 (fits over N = 6–24).
* **Scaling, J2 = 0.4:** guaranteed total CX ~ N^4.9 (R² 0.999), steps ~ N^3.8, CX per step ~ N^1.1 (fits over N = 6–32).
* **Certified bandwidth:** the Weyl bound W = Σ_bonds λmax − E0 equals E_top − E0 (ratio 1.000–1.000); W cancels in n and in T (phys), see Section 6.
* **Bound violations among simulated points:** 0.

## 1. Case inputs (DMRG spectrum, certified bandwidth)

| N | J2 | E0 (DMRG) | E1 | E_top | \|E0−E0(ED)\| | \|E1−E1(ED)\| | E_top−E0 | W (certified) | W / (E_top−E0) | Δ (scaled gap) | γ (trial) | γ source | trial-circuit infidelity 1−γ | α used | Pauli terms | trial CX | trial RZ | trial RZ (non-Cliff.) | trial depth | E1 method | DMRG converged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | -1.6160254038 | -0.9571067812 | 0.7500000000 | 2.7e-15 | 1.1e-16 | 2.3660 | 2.3660 | 1.00 | 0.27849 | 0.95534 | ED | 4.47e-02 | 0.2679 | 9 | 9 | 61 | 59 | 71 | local | yes |
| 6 | 0 | -2.4935771339 | -2.0019953569 | 1.2500000000 | 4.4e-16 | 4.4e-16 | 3.7436 | 3.7436 | 1.00 | 0.13131 | 0.89729 | ED | 1.03e-01 | 0.2141 | 15 | 15 | 98 | 96 | 111 | local | yes |
| 8 | 0 | -3.3749325987 | -2.9822404878 | 1.7500000000 | 7.6e-15 | 1.3e-15 | 5.1249 | 5.1249 | 1.00 | 0.07662 | 0.83650 | ED | 1.64e-01 | 0.1713 | 21 | 21 | 140 | 135 | 159 | local | yes |
| 10 | 0 | -4.2580352073 | -3.9306735895 | 2.2500000000 | 1.2e-14 | 7.1e-15 | 6.5080 | 6.5080 | 1.00 | 0.05030 | 0.77610 | ED | 2.24e-01 | 0.1417 | 27 | 27 | 179 | 173 | 203 | local | yes |
| 12 | 0 | -5.1420906328 | -4.8611479370 | 2.7500000000 | 2.0e-14 | 5.0e-14 | 7.8921 | 7.8921 | 1.00 | 0.03560 | 0.71708 | ED | 2.83e-01 | 0.1204 | 33 | 33 | 218 | 212 | 247 | local | yes |
| 14 | 0 | -6.0267246619 | -5.7804926045 | 3.2500000000 | 1.0e-13 | 1.0e-13 | 9.2767 | 9.2767 | 1.00 | 0.02654 | 0.66091 | ED | 3.39e-01 | 0.1046 | 39 | 39 | 252 | 248 | 291 | local | yes |
| 16 | 0 | -6.9117371456 | -6.6924604290 | 3.7500000000 | n/a (no ED) | n/a (no ED) | 10.6617 | 10.6617 | 1.00 | 0.02057 | 0.60469 | DMRG-MPS overlap (uncertified) | 3.95e-01 | 0.0924 | 45 | 45 | 299 | 288 | 335 | local | yes |
| 18 | 0 | -7.7970110685 | -7.5992848841 | 4.2500000000 | n/a (no ED) | n/a (no ED) | 12.0470 | 12.0470 | 1.00 | 0.01641 | 0.55466 | DMRG-MPS overlap (uncertified) | 4.45e-01 | 0.0827 | 51 | 51 | 335 | 326 | 379 | local | yes |
| 20 | 0 | -8.6824733344 | -8.5023786980 | 4.7500000000 | n/a (no ED) | n/a (no ED) | 13.4325 | 13.4325 | 1.00 | 0.01341 | 0.51386 | DMRG-MPS overlap (uncertified) | 4.86e-01 | 0.0748 | 57 | 57 | 372 | 362 | 423 | local | yes |
| 24 | 0 | -10.4537857604 | -10.3008340576 | 5.7500000000 | n/a (no ED) | n/a (no ED) | 16.2038 | 16.2038 | 1.00 | 0.00944 | 0.43229 | DMRG-MPS overlap (uncertified) | 5.68e-01 | 0.0628 | 69 | 69 | 451 | 436 | 511 | local | yes |
| 28 | 0 | -12.2254405486 | -12.0924250046 | 6.7500000000 | n/a (no ED) | n/a (no ED) | 18.9754 | 18.9754 | 1.00 | 0.00701 | 0.36523 | DMRG-MPS overlap (uncertified) | 6.35e-01 | 0.0542 | 81 | 81 | 520 | 509 | 599 | local | yes |
| 32 | 0 | -13.9973156182 | -13.8795758237 | 7.7500000000 | n/a (no ED) | n/a (no ED) | 21.7473 | 21.7473 | 1.00 | 0.00541 | 0.30813 | DMRG-MPS overlap (uncertified) | 6.92e-01 | 0.0476 | 93 | 93 | 592 | 578 | 687 | local | yes |
| 4 | 0.2411 | -1.5387409981 | -0.8776800957 | 0.8705500000 | 1.6e-15 | 7.8e-16 | 2.4093 | 2.4093 | 1.00 | 0.27438 | 0.98067 | ED | 1.93e-02 | 0.3977 | 15 | 9 | 60 | 58 | 67 | local | yes |
| 6 | 0.2411 | -2.3353474058 | -1.8475542730 | 1.4911000000 | 1.1e-16 | 6.7e-16 | 3.8264 | 3.8264 | 1.00 | 0.12748 | 0.94908 | ED | 5.09e-02 | 0.3295 | 27 | 15 | 100 | 97 | 115 | local | yes |
| 8 | 0.2411 | -3.1350011019 | -2.7491233495 | 2.1116500000 | 3.1e-15 | 4.4e-16 | 5.2467 | 5.2467 | 1.00 | 0.07355 | 0.91107 | ED | 8.89e-02 | 0.2672 | 39 | 21 | 137 | 134 | 159 | local | yes |
| 10 | 0.2411 | -3.9361281870 | -3.6170270856 | 2.7322000000 | 4.4e-15 | 1.7e-14 | 6.6683 | 6.6683 | 1.00 | 0.04785 | 0.86897 | ED | 1.31e-01 | 0.2223 | 51 | 27 | 172 | 168 | 203 | local | yes |
| 12 | 0.2411 | -4.7380770283 | -4.4660645236 | 3.3527500000 | 3.6e-14 | 5.0e-14 | 8.0908 | 8.0908 | 1.00 | 0.03362 | 0.82515 | ED | 1.75e-01 | 0.1897 | 63 | 33 | 216 | 208 | 247 | local | yes |
| 14 | 0.2411 | -5.5405306286 | -5.3034976901 | 3.9733000000 | 1.1e-13 | 7.7e-14 | 9.5138 | 9.5138 | 1.00 | 0.02491 | 0.78098 | ED | 2.19e-01 | 0.1651 | 75 | 39 | 252 | 245 | 291 | local | yes |
| 16 | 0.2411 | -6.3433164211 | -6.1332908175 | 4.5938500000 | n/a (no ED) | n/a (no ED) | 10.9372 | 10.9372 | 1.00 | 0.01920 | 0.73648 | DMRG-MPS overlap (uncertified) | 2.64e-01 | 0.1461 | 87 | 45 | 292 | 284 | 335 | local | yes |
| 18 | 0.2411 | -7.1463324487 | -6.9577882796 | 5.2144000000 | n/a (no ED) | n/a (no ED) | 12.3607 | 12.3607 | 1.00 | 0.01525 | 0.69439 | DMRG-MPS overlap (uncertified) | 3.06e-01 | 0.1310 | 99 | 51 | 332 | 320 | 379 | local | yes |
| 20 | 0.2411 | -7.9495146057 | -7.7784645690 | 5.8349500000 | n/a (no ED) | n/a (no ED) | 13.7845 | 13.7845 | 1.00 | 0.01241 | 0.65278 | DMRG-MPS overlap (uncertified) | 3.47e-01 | 0.1186 | 111 | 57 | 363 | 355 | 423 | local | yes |
| 24 | 0.2411 | -9.5562212540 | -9.4119433252 | 7.0760500000 | n/a (no ED) | n/a (no ED) | 16.6323 | 16.6323 | 1.00 | 0.00867 | 0.57568 | DMRG-MPS overlap (uncertified) | 4.24e-01 | 0.0998 | 135 | 69 | 443 | 431 | 511 | local | yes |
| 28 | 0.2411 | -11.1632298850 | -11.0384766852 | 8.3171500000 | n/a (no ED) | n/a (no ED) | 19.4804 | 19.4804 | 1.00 | 0.00640 | 0.50458 | DMRG-MPS overlap (uncertified) | 4.95e-01 | 0.0861 | 159 | 81 | 517 | 504 | 599 | local | yes |
| 32 | 0.2411 | -12.7704324927 | -12.6605489205 | 9.5582500000 | n/a (no ED) | n/a (no ED) | 22.3287 | 22.3287 | 1.00 | 0.00492 | 0.44173 | DMRG-MPS overlap (uncertified) | 5.58e-01 | 0.0757 | 183 | 93 | 599 | 576 | 687 | local | yes |
| 4 | 0.4 | -1.5067764363 | -0.8330951895 | 0.9500000000 | 6.7e-16 | 1.1e-16 | 2.4568 | 2.4568 | 1.00 | 0.27421 | 0.99594 | ED | 4.06e-03 | 0.4708 | 15 | 9 | 61 | 58 | 70 | local | yes |
| 6 | 0.4 | -2.2655348801 | -1.7633485573 | 1.6500000000 | 4.9e-15 | 2.2e-15 | 3.9155 | 3.9155 | 1.00 | 0.12825 | 0.98805 | ED | 1.20e-02 | 0.3957 | 27 | 15 | 98 | 96 | 111 | local | yes |
| 8 | 0.4 | -3.0251905645 | -2.6214818416 | 2.3500000000 | 1.8e-15 | 4.0e-15 | 5.3752 | 5.3752 | 1.00 | 0.07511 | 0.97704 | ED | 2.30e-02 | 0.3229 | 39 | 21 | 138 | 134 | 155 | local | yes |
| 10 | 0.4 | -3.7852888948 | -3.4451673343 | 3.0500000000 | 1.8e-14 | 1.3e-14 | 6.8353 | 6.8353 | 1.00 | 0.04976 | 0.96379 | ED | 3.62e-02 | 0.2695 | 51 | 27 | 180 | 173 | 203 | local | yes |
| 12 | 0.4 | -4.5456261627 | -4.2500350676 | 3.7500000000 | 1.2e-13 | 4.1e-14 | 8.2956 | 8.2956 | 1.00 | 0.03563 | 0.94899 | ED | 5.10e-02 | 0.2303 | 63 | 33 | 219 | 212 | 247 | local | yes |
| 14 | 0.4 | -5.3061027329 | -5.0435118006 | 4.4500000000 | 1.3e-13 | 9.2e-14 | 9.7561 | 9.7561 | 1.00 | 0.02692 | 0.93337 | ED | 6.66e-02 | 0.2008 | 75 | 39 | 254 | 247 | 291 | local | yes |
| 16 | 0.4 | -6.0666655720 | -5.8295569546 | 5.1500000000 | n/a (no ED) | n/a (no ED) | 11.2167 | 11.2167 | 1.00 | 0.02114 | 0.91702 | DMRG-MPS overlap (uncertified) | 8.30e-02 | 0.1778 | 87 | 45 | 294 | 285 | 335 | local | yes |
| 18 | 0.4 | -6.8272844181 | -6.6104747998 | 5.8500000000 | n/a (no ED) | n/a (no ED) | 12.6773 | 12.6773 | 1.00 | 0.01710 | 0.90000 | DMRG-MPS overlap (uncertified) | 1.00e-01 | 0.1595 | 99 | 51 | 333 | 324 | 379 | local | yes |
| 20 | 0.4 | -7.5879410154 | -7.3876995115 | 6.5500000000 | n/a (no ED) | n/a (no ED) | 14.1379 | 14.1379 | 1.00 | 0.01416 | 0.88328 | DMRG-MPS overlap (uncertified) | 1.17e-01 | 0.1446 | 111 | 57 | 366 | 354 | 423 | local | yes |
| 24 | 0.4 | -9.1093253675 | -8.9345357926 | 7.9500000000 | n/a (no ED) | n/a (no ED) | 17.0593 | 17.0593 | 1.00 | 0.01025 | 0.84890 | DMRG-MPS overlap (uncertified) | 1.51e-01 | 0.1217 | 135 | 69 | 444 | 432 | 511 | local | yes |
| 28 | 0.4 | -10.6307656965 | -10.4746379729 | 9.3500000000 | n/a (no ED) | n/a (no ED) | 19.9808 | 19.9808 | 1.00 | 0.00781 | 0.81538 | DMRG-MPS overlap (uncertified) | 1.85e-01 | 0.1051 | 159 | 81 | 523 | 506 | 599 | local | yes |
| 32 | 0.4 | -12.1522370500 | -12.0103923579 | 10.7500000000 | n/a (no ED) | n/a (no ED) | 22.9022 | 22.9022 | 1.00 | 0.00619 | 0.78219 | DMRG-MPS overlap (uncertified) | 2.18e-01 | 0.0924 | 183 | 93 | 597 | 577 | 686 | local | yes |

*Table 1.* The energies are DMRG values (10 digits); `|E0−E0(ED)|`, `|E1−E1(ED)|` are the DMRG errors against exact diagonalisation (N ≤ 14 here; N = 16–20 in Section 5; N ≥ 24 has no exact reference, see the run-to-run table in Section 5). **1−γ is not a DMRG error:** it is the infidelity of the low-depth (L = 1 brickwork) trial *circuit* with respect to the ground state, which is large by construction and is what the filter has to remove. *DMRG converged*: the ground-state, top-state and orthogonality-penalty excited-state DMRG runs each reached the sweep tolerance within 100 sweeps; a **NO** means that quantity (for the excited state, E1 and hence the gap Δ) is not converged and the design is only indicative. *E1 method*: `local` = projector penalty λ|ψ0⟩⟨ψ0| applied as a rank-1 term in every two-site local problem, giving the lowest state of any symmetry (checked against ED in Section 5). Trial circuit = L = 1 approximate brickwork circuit from the DMRG MPS. α is in scaled units: the exact commutator sum Σ_g‖[H_g, Σ_{g'>g}H_g']‖ for every N (computed on the local support of each term). W is the certified (Weyl) bandwidth, which coincides with E_top − E0 here, so Δ = (E1−E0)/W is the exact scaled gap.

**α: exact commutator sum against the triangle bound**

| N | J2 | α (exact commutator sum) | α triangle (diagnostic) | triangle / exact |
|---|---|---|---|---|
| 4 | 0 | 0.26795 | 0.26795 | 1.000 |
| 6 | 0 | 0.21407 | 0.21407 | 1.000 |
| 8 | 0 | 0.17133 | 0.17133 | 1.000 |
| 10 | 0 | 0.14166 | 0.14166 | 1.000 |
| 12 | 0 | 0.12041 | 0.12041 | 1.000 |
| 14 | 0 | 0.10458 | 0.10458 | 1.000 |
| 16 | 0 | 0.09237 | 0.09237 | 1.000 |
| 18 | 0 | 0.08268 | 0.08268 | 1.000 |
| 20 | 0 | 0.07482 | 0.07482 | 1.000 |
| 24 | 0 | 0.06284 | 0.06284 | 1.000 |
| 28 | 0 | 0.05416 | 0.05416 | 1.000 |
| 32 | 0 | 0.04757 | 0.04757 | 1.000 |
| 4 | 0.2411 | 0.39773 | 0.44532 | 1.120 |
| 6 | 0.2411 | 0.32946 | 0.38375 | 1.165 |
| 8 | 0.2411 | 0.26717 | 0.31433 | 1.177 |
| 10 | 0.2411 | 0.22230 | 0.26281 | 1.182 |
| 12 | 0.2411 | 0.18967 | 0.22487 | 1.186 |
| 14 | 0.2411 | 0.16513 | 0.19615 | 1.188 |
| 16 | 0.2411 | 0.14610 | 0.17378 | 1.189 |
| 18 | 0.2411 | 0.13095 | 0.15591 | 1.191 |
| 20 | 0.2411 | 0.11862 | 0.14134 | 1.192 |
| 24 | 0.2411 | 0.09977 | 0.11901 | 1.193 |
| 28 | 0.2411 | 0.08607 | 0.10275 | 1.194 |
| 32 | 0.2411 | 0.07566 | 0.09037 | 1.194 |
| 4 | 0.4 | 0.47080 | 0.54674 | 1.161 |
| 6 | 0.4 | 0.39573 | 0.48528 | 1.226 |
| 8 | 0.4 | 0.32286 | 0.40079 | 1.241 |
| 10 | 0.4 | 0.26946 | 0.33646 | 1.249 |
| 12 | 0.4 | 0.23033 | 0.28859 | 1.253 |
| 14 | 0.4 | 0.20080 | 0.25215 | 1.256 |
| 16 | 0.4 | 0.17783 | 0.22366 | 1.258 |
| 18 | 0.4 | 0.15951 | 0.20085 | 1.259 |
| 20 | 0.4 | 0.14457 | 0.18221 | 1.260 |
| 24 | 0.4 | 0.12170 | 0.15360 | 1.262 |
| 28 | 0.4 | 0.10505 | 0.13271 | 1.263 |
| 32 | 0.4 | 0.09240 | 0.11679 | 1.264 |

*Table 1b.* α is the exact value of Σ_g‖[H_g, Σ_{g'>g}H_g']‖ for every N (for Pauli terms it reduces to a small dense eigenvalue problem on the support of each term). The triangle-inequality value, which earlier versions used for N > 14, is shown for comparison: it is exact at J2 = 0 and up to about 1.25× larger otherwise.

## 2. Guaranteed designs, N ≤ 14, ε = 10⁻²

**Guaranteed: design and evolution time**

| N | J2 | pulses | n (steps) | kᵢ | T (scaled) | T (phys) | η cert. | ε bound | P_succ lb |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 4 | 85 | [1, 16, 25, 43] | 5.64 | 2.38 | 0.143 | 9.91e-03 | 0.535 |
| 6 | 0 | 4 | 591 | [64, 104, 161, 262] | 17.94 | 4.79 | 0.092 | 9.99e-03 | 0.716 |
| 8 | 0 | 6 | 1,819 | [15, 107, 168, 283, 464, 782] | 30.75 | 6.00 | 0.099 | 9.99e-03 | 0.631 |
| 10 | 0 | 6 | 3,905 | [134, 214, 364, 666, 1013, 1514] | 46.84 | 7.20 | 0.082 | 1.00e-02 | 0.503 |
| 12 | 0 | 8 | 8,800 | [55, 168, 292, 494, 875, 1469, 2112, 3335] | 88.25 | 11.18 | 0.050 | 9.99e-03 | 0.600 |
| 14 | 0 | 8 | 16,975 | [106, 307, 543, 917, 1606, 2891, 4199, 6406] | 118.36 | 12.76 | 0.062 | 1.00e-02 | 0.595 |
| 4 | 0.2411 | 4 | 94 | [1, 10, 36, 47] | 4.58 | 1.90 | 0.313 | 9.93e-03 | 0.631 |
| 6 | 0.2411 | 4 | 579 | [44, 83, 162, 290] | 12.32 | 3.22 | 0.189 | 9.98e-03 | 0.594 |
| 8 | 0.2411 | 6 | 2,398 | [20, 141, 221, 373, 612, 1031] | 32.04 | 6.11 | 0.099 | 1.00e-02 | 0.687 |
| 10 | 0.2411 | 6 | 5,777 | [127, 284, 494, 839, 1469, 2564] | 49.24 | 7.38 | 0.113 | 1.00e-02 | 0.691 |
| 12 | 0.2411 | 6 | 10,960 | [418, 756, 881, 1650, 2733, 4522] | 70.08 | 8.66 | 0.095 | 1.00e-02 | 0.574 |
| 14 | 0.2411 | 8 | 22,902 | [143, 413, 769, 1260, 2167, 3832, 5619, 8699] | 126.09 | 13.25 | 0.059 | 9.99e-03 | 0.694 |
| 4 | 0.4 | 4 | 76 | [1, 4, 33, 38] | 4.58 | 1.87 | 0.486 | 9.98e-03 | 0.891 |
| 6 | 0.4 | 4 | 394 | [5, 56, 136, 197] | 9.80 | 2.50 | 0.399 | 9.98e-03 | 0.739 |
| 8 | 0.4 | 4 | 1,058 | [58, 145, 326, 529] | 16.73 | 3.11 | 0.286 | 1.00e-02 | 0.578 |
| 10 | 0.4 | 4 | 2,891 | [173, 357, 915, 1446] | 25.25 | 3.69 | 0.321 | 1.00e-02 | 0.613 |
| 12 | 0.4 | 6 | 5,320 | [44, 156, 334, 718, 1408, 2660] | 44.08 | 5.31 | 0.189 | 1.00e-02 | 0.562 |
| 14 | 0.4 | 6 | 10,571 | [88, 306, 706, 1430, 2756, 5285] | 70.03 | 7.18 | 0.164 | 1.00e-02 | 0.689 |

**Guaranteed: gate counts (level-3 columns include the trial circuit; direct columns are filter only)**

| N | J2 | CX trial | CX filter | CX total | CX filter (direct) | RZ non-Cliff. (level 3) | RZ non-Cliff. (direct) | depth total | E[CX]/accept | CX total (line) |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 9 | 1,785 | 1,794 | 3,060 | 2,528 | 854 | 7,636 | 3,353 | 2,304 |
| 6 | 0 | 15 | 20,685 | 20,700 | 35,460 | 27,877 | 9,460 | 86,988 | 28,889 | 27,792 |
| 8 | 0 | 21 | 89,131 | 89,152 | 152,796 | 118,376 | 40,024 | 373,054 | 141,202 | 121,894 |
| 10 | 0 | 27 | 246,015 | 246,042 | 421,740 | 324,294 | 109,346 | 1,027,218 | 488,487 | 339,762 |
| 12 | 0 | 33 | 677,600 | 677,633 | 1,161,600 | 889,020 | 299,208 | 2,825,047 | 1,129,200 | 941,633 |
| 14 | 0 | 39 | 1,544,725 | 1,544,764 | 2,648,100 | 2,020,281 | 679,008 | 6,433,816 | 2,596,452 | 2,155,864 |
| 4 | 0.2411 | 9 | 3,290 | 3,299 | 5,640 | 4,480 | 1,508 | 13,885 | 5,217 | 5,398 |
| 6 | 0.2411 | 15 | 36,477 | 36,492 | 62,532 | 48,158 | 16,216 | 152,392 | 61,335 | 67,951 |
| 8 | 0.2411 | 21 | 218,218 | 218,239 | 374,088 | 285,502 | 95,926 | 909,001 | 317,694 | 393,293 |
| 10 | 0.2411 | 27 | 687,463 | 687,490 | 1,178,508 | 895,609 | 300,410 | 2,859,818 | 994,954 | 1,309,480 |
| 12 | 0.2411 | 33 | 1,611,120 | 1,611,153 | 2,761,920 | 2,093,574 | 701,446 | 6,696,807 | 2,803,988 | 3,098,060 |
| 14 | 0.2411 | 39 | 4,007,850 | 4,007,889 | 6,870,600 | 5,199,007 | 1,740,560 | 16,650,045 | 5,772,340 | 7,756,183 |
| 4 | 0.4 | 9 | 2,660 | 2,669 | 4,560 | 3,634 | 1,220 | 11,242 | 2,993 | 4,366 |
| 6 | 0.4 | 15 | 24,822 | 24,837 | 42,552 | 32,802 | 11,036 | 103,733 | 33,581 | 46,244 |
| 8 | 0.4 | 21 | 96,278 | 96,299 | 165,048 | 126,040 | 42,324 | 401,137 | 166,370 | 173,533 |
| 10 | 0.4 | 27 | 344,029 | 344,056 | 589,764 | 448,282 | 150,336 | 1,431,248 | 560,526 | 655,320 |
| 12 | 0.4 | 33 | 782,040 | 782,073 | 1,340,640 | 1,016,338 | 340,486 | 3,250,767 | 1,391,435 | 1,503,820 |
| 14 | 0.4 | 39 | 1,849,925 | 1,849,964 | 3,171,300 | 2,399,870 | 803,402 | 7,685,408 | 2,683,271 | 3,580,084 |

**Guaranteed: measured quality (noiseless)**

| N | J2 | measured 1−F | F (ED) | P_succ meas. | state err vs exact | bound holds | full-transpile CX (n ≤ 200) |
|---|---|---|---|---|---|---|---|
| 4 | 0 | 2.97e-04 | 0.99970 | 0.535 | 1.17e-04 | yes | 1,785 |
| 6 | 0 | 6.87e-04 | 0.99931 | 0.717 | 5.09e-04 | yes | – |
| 8 | 0 | 1.25e-03 | 0.99875 | 0.631 | 3.83e-04 | yes | – |
| 10 | 0 | 1.23e-03 | 0.99877 | 0.504 | 3.10e-04 | yes | – |
| 12 | 0 | 3.29e-04 | 0.99967 | 0.600 | 1.43e-04 | yes | – |
| 14 | 0 | 6.34e-04 | 0.99937 | 0.595 | 1.07e-04 | yes | – |
| 4 | 0.2411 | 9.79e-04 | 0.99902 | 0.632 | 1.63e-04 | yes | 3,290 |
| 6 | 0.2411 | 1.56e-03 | 0.99844 | 0.595 | 6.37e-04 | yes | – |
| 8 | 0.2411 | 2.42e-04 | 0.99976 | 0.687 | 3.52e-04 | yes | – |
| 10 | 0.2411 | 5.03e-04 | 0.99950 | 0.691 | 2.00e-04 | yes | – |
| 12 | 0.2411 | 4.09e-04 | 0.99959 | 0.575 | 1.59e-04 | yes | – |
| 14 | 0.2411 | 3.83e-04 | 0.99962 | 0.694 | 1.48e-04 | yes | – |
| 4 | 0.4 | 2.50e-04 | 0.99975 | 0.892 | 3.60e-03 | yes | 2,660 |
| 6 | 0.4 | 4.13e-04 | 0.99959 | 0.740 | 9.72e-04 | yes | – |
| 8 | 0.4 | 5.06e-04 | 0.99949 | 0.579 | 9.33e-04 | yes | – |
| 10 | 0.4 | 7.36e-04 | 0.99926 | 0.614 | 5.27e-04 | yes | – |
| 12 | 0.4 | 4.80e-04 | 0.99952 | 0.562 | 4.34e-04 | yes | – |
| 14 | 0.4 | 4.64e-04 | 0.99954 | 0.689 | 3.04e-04 | yes | – |

*Tables 2a–c.* Rigorous a-priori n: the bound on 1−F is ≤ ε. Measured 1−F is typically far below ε (the bound is conservative).

## 3. Empirical (oracle) designs, N ≤ 14, ε = 10⁻²

**Empirical: design and evolution time**

| N | J2 | pulses | n (steps) | kᵢ | T (scaled) | T (phys) | η cert. | ε bound | P_succ lb |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 2 | 2 | [1, 1] | 5.64 | 2.38 | 0.764 | 1.00e+00 | 0.857 |
| 6 | 0 | 4 | 5 | [1, 1, 1, 2] | 17.94 | 4.79 | 1.000 | 1.00e+00 | 0.897 |
| 8 | 0 | 4 | 11 | [1, 2, 3, 5] | 30.75 | 6.00 | 0.209 | 1.00e+00 | 0.683 |
| 10 | 0 | 4 | 10 | [1, 2, 3, 4] | 46.84 | 7.20 | 1.002 | 1.00e+00 | 0.773 |
| 12 | 0 | 6 | 29 | [1, 2, 3, 5, 7, 11] | 88.25 | 11.18 | 0.176 | 1.00e+00 | 0.248 |
| 14 | 0 | 7 | 39 | [1, 1, 2, 4, 6, 10, 15] | 118.36 | 12.76 | 0.091 | 1.00e+00 | 0.306 |
| 4 | 0.2411 | 2 | 3 | [1, 2] | 4.58 | 1.90 | 0.313 | 1.00e+00 | 0.820 |
| 6 | 0.2411 | 3 | 6 | [1, 2, 3] | 12.32 | 3.22 | 0.268 | 1.00e+00 | 0.719 |
| 8 | 0.2411 | 4 | 11 | [1, 2, 3, 5] | 32.04 | 6.11 | 0.267 | 1.00e+00 | 0.661 |
| 10 | 0.2411 | 5 | 17 | [1, 1, 3, 4, 8] | 49.24 | 7.38 | 0.163 | 1.00e+00 | 0.318 |
| 12 | 0.2411 | 6 | 25 | [1, 2, 2, 4, 6, 10] | 70.08 | 8.66 | 0.095 | 1.00e+00 | 0.565 |
| 14 | 0.2411 | 7 | 41 | [1, 1, 2, 4, 7, 10, 16] | 126.09 | 13.25 | 0.208 | 1.00e+00 | 0.396 |
| 4 | 0.4 | 2 | 5 | [2, 3] | 4.58 | 1.87 | 0.486 | 1.00e+00 | 0.945 |
| 6 | 0.4 | 3 | 6 | [1, 2, 3] | 9.80 | 2.50 | 0.399 | 1.00e+00 | 0.738 |
| 8 | 0.4 | 4 | 11 | [1, 2, 3, 5] | 16.73 | 3.11 | 0.286 | 1.00e+00 | 0.503 |
| 10 | 0.4 | 4 | 16 | [1, 2, 5, 8] | 25.25 | 3.69 | 0.320 | 1.00e+00 | 0.612 |
| 12 | 0.4 | 5 | 26 | [1, 2, 3, 7, 13] | 44.08 | 5.31 | 0.189 | 1.00e+00 | 0.537 |
| 14 | 0.4 | 5 | 38 | [1, 3, 5, 10, 19] | 70.03 | 7.18 | 0.189 | 1.00e+00 | 0.645 |

**Empirical: gate counts (level-3 columns include the trial circuit; direct columns are filter only)**

| N | J2 | CX trial | CX filter | CX total | CX filter (direct) | RZ non-Cliff. (level 3) | RZ non-Cliff. (direct) | depth total | E[CX]/accept | CX total (line) |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 9 | 42 | 51 | 72 | 119 | 22 | 249 | 57 | 63 |
| 6 | 0 | 15 | 175 | 190 | 300 | 335 | 84 | 846 | 223 | 250 |
| 8 | 0 | 21 | 539 | 560 | 924 | 854 | 246 | 2,414 | 848 | 758 |
| 10 | 0 | 27 | 630 | 657 | 1,080 | 1,007 | 284 | 2,833 | 981 | 897 |
| 12 | 0 | 33 | 2,233 | 2,266 | 3,828 | 3,147 | 992 | 9,556 | 9,423 | 3,136 |
| 14 | 0 | 39 | 3,549 | 3,588 | 6,084 | 4,896 | 1,567 | 15,072 | 13,035 | 4,992 |
| 4 | 0.2411 | 9 | 105 | 114 | 180 | 201 | 50 | 508 | 143 | 181 |
| 6 | 0.2411 | 15 | 378 | 393 | 648 | 598 | 171 | 1,693 | 570 | 719 |
| 8 | 0.2411 | 21 | 1,001 | 1,022 | 1,716 | 1,447 | 444 | 4,328 | 1,606 | 1,825 |
| 10 | 0.2411 | 27 | 2,023 | 2,050 | 3,468 | 2,808 | 889 | 8,618 | 7,171 | 3,880 |
| 12 | 0.2411 | 33 | 3,675 | 3,708 | 6,300 | 4,989 | 1,606 | 15,522 | 7,719 | 7,100 |
| 14 | 0.2411 | 39 | 7,175 | 7,214 | 12,300 | 9,559 | 3,123 | 30,098 | 18,659 | 13,924 |
| 4 | 0.4 | 9 | 175 | 184 | 300 | 295 | 82 | 805 | 196 | 296 |
| 6 | 0.4 | 15 | 378 | 393 | 648 | 597 | 171 | 1,689 | 564 | 719 |
| 8 | 0.4 | 21 | 1,001 | 1,022 | 1,716 | 1,447 | 444 | 4,324 | 2,146 | 1,825 |
| 10 | 0.4 | 27 | 1,904 | 1,931 | 3,264 | 2,657 | 836 | 8,123 | 3,286 | 3,654 |
| 12 | 0.4 | 33 | 3,822 | 3,855 | 6,552 | 5,183 | 1,669 | 16,133 | 7,688 | 7,382 |
| 14 | 0.4 | 39 | 6,650 | 6,689 | 11,400 | 8,878 | 2,893 | 27,917 | 10,971 | 12,908 |

**Empirical: measured quality (noiseless)**

| N | J2 | measured 1−F | F (ED) | P_succ meas. | state err vs exact | bound holds | full-transpile CX (n ≤ 200) |
|---|---|---|---|---|---|---|---|
| 4 | 0 | 8.49e-03 | 0.99151 | 0.889 | 5.83e-02 | yes | 42 |
| 6 | 0 | 6.96e-03 | 0.99304 | 0.850 | 8.84e-02 | yes | 175 |
| 8 | 0 | 5.25e-03 | 0.99475 | 0.661 | 7.56e-02 | yes | 539 |
| 10 | 0 | 7.85e-03 | 0.99215 | 0.670 | 1.04e-01 | yes | 630 |
| 12 | 0 | 5.75e-03 | 0.99425 | 0.240 | 6.28e-02 | yes | 2,233 |
| 14 | 0 | 3.45e-03 | 0.99655 | 0.275 | 6.12e-02 | yes | 3,549 |
| 4 | 0.2411 | 5.78e-03 | 0.99422 | 0.798 | 4.34e-02 | yes | 105 |
| 6 | 0.2411 | 9.50e-03 | 0.99050 | 0.689 | 6.13e-02 | yes | 378 |
| 8 | 0.2411 | 7.74e-03 | 0.99226 | 0.636 | 8.07e-02 | yes | 1,001 |
| 10 | 0.2411 | 7.05e-03 | 0.99295 | 0.286 | 7.54e-02 | yes | 2,023 |
| 12 | 0.2411 | 7.80e-03 | 0.99220 | 0.480 | 8.48e-02 | yes | 3,675 |
| 14 | 0.2411 | 7.83e-03 | 0.99217 | 0.387 | 9.72e-02 | yes | 7,175 |
| 4 | 0.4 | 3.58e-03 | 0.99642 | 0.938 | 6.24e-02 | yes | 175 |
| 6 | 0.4 | 6.01e-03 | 0.99399 | 0.697 | 7.98e-02 | yes | 378 |
| 8 | 0.4 | 8.69e-03 | 0.99131 | 0.476 | 9.39e-02 | yes | 1,001 |
| 10 | 0.4 | 9.44e-03 | 0.99056 | 0.588 | 9.53e-02 | yes | 1,904 |
| 12 | 0.4 | 9.35e-03 | 0.99065 | 0.501 | 8.95e-02 | yes | 3,822 |
| 14 | 0.4 | 9.75e-03 | 0.99025 | 0.610 | 8.73e-02 | yes | 6,650 |

*Tables 3a–c.* Smallest n on the sweep whose *measured* 1−F ≤ ε; an oracle (uses the exact ground state). At very small n the first-order Trotter error can be accidentally small, so the empirical n is a lower envelope, not a design a practitioner could trust.

## 4. Cost of one Trotter step

| N | J2 | bonds | Pauli terms | CX (level 3) | RZ non-Cliff. (level 3) | CX (direct) | RZ non-Cliff. (direct) | depth (level 3) | CX (line, 3-step avg) | line / a2a |
|---|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 3 | 9 | 21 | 29 | 36 | 10 | 89 | 27.0 | 1.29 |
| 6 | 0 | 5 | 15 | 35 | 47 | 60 | 16 | 147 | 47.0 | 1.34 |
| 8 | 0 | 7 | 21 | 49 | 65 | 84 | 22 | 205 | 67.0 | 1.37 |
| 10 | 0 | 9 | 27 | 63 | 83 | 108 | 28 | 263 | 87.0 | 1.38 |
| 12 | 0 | 11 | 33 | 77 | 101 | 132 | 34 | 321 | 107.0 | 1.39 |
| 14 | 0 | 13 | 39 | 91 | 119 | 156 | 40 | 379 | 127.0 | 1.40 |
| 16 | 0 | 15 | 45 | 105 | 137 | 180 | 46 | 437 | 147.0 | 1.40 |
| 18 | 0 | 17 | 51 | 119 | 155 | 204 | 52 | 495 | 167.0 | 1.40 |
| 20 | 0 | 19 | 57 | 133 | 173 | 228 | 58 | 553 | 187.0 | 1.41 |
| 24 | 0 | 23 | 69 | 161 | 209 | 276 | 70 | 669 | 232.3 | 1.44 |
| 28 | 0 | 27 | 81 | 189 | 245 | 324 | 82 | 785 | 622.3 | 3.29 |
| 32 | 0 | 31 | 93 | 217 | 281 | 372 | 94 | 901 | 742.0 | 3.42 |
| 4 | 0.2411 | 5 | 15 | 35 | 47 | 60 | 16 | 147 | 57.3 | 1.64 |
| 6 | 0.2411 | 9 | 27 | 63 | 83 | 108 | 28 | 263 | 117.3 | 1.86 |
| 8 | 0.2411 | 13 | 39 | 91 | 119 | 156 | 40 | 379 | 164.0 | 1.80 |
| 10 | 0.2411 | 17 | 51 | 119 | 155 | 204 | 52 | 495 | 226.7 | 1.90 |
| 12 | 0.2411 | 21 | 63 | 147 | 191 | 252 | 64 | 611 | 282.7 | 1.92 |
| 14 | 0.2411 | 25 | 75 | 175 | 227 | 300 | 76 | 727 | 338.7 | 1.94 |
| 16 | 0.2411 | 29 | 87 | 203 | 263 | 348 | 88 | 843 | 394.7 | 1.94 |
| 18 | 0.2411 | 33 | 99 | 231 | 299 | 396 | 100 | 959 | 450.7 | 1.95 |
| 20 | 0.2411 | 37 | 111 | 259 | 335 | 444 | 112 | 1,075 | 506.7 | 1.96 |
| 24 | 0.2411 | 45 | 135 | 315 | 407 | 540 | 136 | 1,307 | 1230.7 | 3.91 |
| 28 | 0.2411 | 53 | 159 | 371 | 479 | 636 | 160 | 1,539 | 1999.0 | 5.39 |
| 32 | 0.2411 | 61 | 183 | 427 | 551 | 732 | 184 | 1,771 | 2540.7 | 5.95 |
| 4 | 0.4 | 5 | 15 | 35 | 47 | 60 | 16 | 147 | 57.3 | 1.64 |
| 6 | 0.4 | 9 | 27 | 63 | 83 | 108 | 28 | 263 | 117.3 | 1.86 |
| 8 | 0.4 | 13 | 39 | 91 | 119 | 156 | 40 | 379 | 164.0 | 1.80 |
| 10 | 0.4 | 17 | 51 | 119 | 155 | 204 | 52 | 495 | 226.7 | 1.90 |
| 12 | 0.4 | 21 | 63 | 147 | 191 | 252 | 64 | 611 | 282.7 | 1.92 |
| 14 | 0.4 | 25 | 75 | 175 | 227 | 300 | 76 | 727 | 338.7 | 1.94 |
| 16 | 0.4 | 29 | 87 | 203 | 263 | 348 | 88 | 843 | 394.7 | 1.94 |
| 18 | 0.4 | 33 | 99 | 231 | 299 | 396 | 100 | 959 | 450.7 | 1.95 |
| 20 | 0.4 | 37 | 111 | 259 | 335 | 444 | 112 | 1,075 | 506.7 | 1.96 |
| 24 | 0.4 | 45 | 135 | 315 | 407 | 540 | 136 | 1,307 | 1230.7 | 3.91 |
| 28 | 0.4 | 53 | 159 | 371 | 479 | 636 | 160 | 1,539 | 1999.0 | 5.39 |
| 32 | 0.4 | 61 | 183 | 427 | 551 | 732 | 184 | 1,771 | 2540.7 | 5.95 |

*Table 4.* Transpiled counts of one Trotter step of one pulse (all bonds, ancilla-controlled), all-to-all, plus the line-connectivity cost. CX per step grows linearly with the number of bonds; J2 ≠ 0 roughly doubles it (J2 = 0.2411 and 0.4 have identical circuits per step, only the angles differ). The line-connectivity figure depends on the transpiler's routing, which is not smooth in N (e.g. N = 24, J2 = 0.2411); treat it as indicative.

## 5. Accuracy of the DMRG inputs

The certificate depends on E0, E1 (through Δ), the DMRG cross-check of E_top, and γ, all of which come from DMRG (`get_energies.py`). The settings, chosen so that two-site DMRG reaches roughly machine precision on these 1D chains:

| setting | value | role |
|---|---|---|
| bond dimensions | 8, 16, 32, 64, 128, 256 (one per sweep stage) | the cutoff, not this cap, limits the bond dimension here (χ reached: 40–100 for N = 12–24) |
| SVD cutoff | 1e-14 | discarded weight per bond ≈ energy error in 1D |
| energy tolerance | 1e-12 | stop when the energy changes by less between sweeps |
| inner Lanczos tolerance | 1e-10 (retry ×100, ×1e4 if ARPACK fails) | quimb's default is 1e-3 |
| E1 (all N) | H + 1.1·bandwidth·\|ψ0⟩⟨ψ0\|, projector applied as a rank-1 term in each two-site local solve | lowest state of any symmetry; avoids the χ²-bond projector MPO, which ran out of memory at N ≥ 16 with these tolerances |
| max sweeps | 100 | convergence flags are recorded in Table 1 |

**DMRG against exact diagonalisation** (N ≤ 14 from the main runs; N ≥ 16 from `scripts/validate_dmrg.py ed`, exact reference used for testing only)

| N | J2 | \|E0 − E0(ED)\| | \|E1 − E1(ED)\| | \|E_top − E_top(ED)\| | \|ΔE1\| / gap | 1 − \|⟨ψ0\|ψ0(ED)⟩\|² | E1 method | converged | max bond |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 0 | 2.7e-15 | 1.1e-16 | 1.0e-15 | 1.7e-16 | 2.2e-16 | local | yes | 4 |
| 6 | 0 | 4.4e-16 | 4.4e-16 | 4.6e-16 | 9.0e-16 | -1.8e-15 | local | yes | 8 |
| 8 | 0 | 7.6e-15 | 1.3e-15 | 2.5e-15 | 3.4e-15 | -8.9e-16 | local | yes | 16 |
| 10 | 0 | 1.2e-14 | 7.1e-15 | 1.4e-15 | 2.2e-14 | -1.3e-15 | local | yes | 32 |
| 12 | 0 | 2.0e-14 | 5.0e-14 | 3.6e-15 | 1.8e-13 | 4.2e-15 | local | yes | 44 |
| 14 | 0 | 1.0e-13 | 1.0e-13 | 3.6e-15 | 4.1e-13 | 3.4e-14 | local | yes | 49 |
| 4 | 0.2411 | 1.6e-15 | 7.8e-16 | 6.7e-16 | 1.2e-15 | 0.0e+00 | local | yes | 4 |
| 6 | 0.2411 | 1.1e-16 | 6.7e-16 | 4.0e-15 | 1.4e-15 | 0.0e+00 | local | yes | 8 |
| 8 | 0.2411 | 3.1e-15 | 4.4e-16 | 2.7e-15 | 1.2e-15 | -1.3e-15 | local | yes | 16 |
| 10 | 0.2411 | 4.4e-15 | 1.7e-14 | 6.7e-15 | 5.4e-14 | 0.0e+00 | local | yes | 32 |
| 12 | 0.2411 | 3.6e-14 | 5.0e-14 | 1.8e-15 | 1.8e-13 | 5.6e-15 | local | yes | 42 |
| 14 | 0.2411 | 1.1e-13 | 7.7e-14 | 1.0e-14 | 3.3e-13 | 2.8e-14 | local | yes | 53 |
| 4 | 0.4 | 6.7e-16 | 1.1e-16 | 1.8e-15 | 1.6e-16 | 0.0e+00 | local | yes | 4 |
| 6 | 0.4 | 4.9e-15 | 2.2e-15 | 1.1e-15 | 4.4e-15 | 2.2e-16 | local | yes | 8 |
| 8 | 0.4 | 1.8e-15 | 4.0e-15 | 8.9e-15 | 9.9e-15 | 1.1e-15 | local | yes | 16 |
| 10 | 0.4 | 1.8e-14 | 1.3e-14 | 1.8e-15 | 3.8e-14 | 1.3e-15 | local | yes | 29 |
| 12 | 0.4 | 1.2e-13 | 4.1e-14 | 1.2e-14 | 1.4e-13 | 2.0e-14 | local | yes | 46 |
| 14 | 0.4 | 1.3e-13 | 9.2e-14 | 2.5e-14 | 3.5e-13 | 3.6e-14 | local | yes | 60 |
| 16 | 0 | 2.4e-13 | 2.1e-13 | 1.4e-14 | 1.5e-13 | 4.4e-14 | local | yes | 62 |
| 18 | 0 | 1.9e-13 | 2.7e-13 | 8.9e-16 | 4.1e-13 | 6.0e-14 | local | yes | 66 |
| 20 | 0 | 3.5e-13 | 3.3e-13 | 2.7e-14 | 8.9e-14 | 8.1e-14 | local | yes | 77 |
| 16 | 0.2411 | 1.6e-13 | 2.4e-13 | 6.2e-15 | 3.6e-13 | 5.1e-14 | local | yes | 65 |
| 18 | 0.2411 | 2.2e-13 | 2.6e-13 | 1.3e-13 | 2.3e-13 | 5.4e-14 | local | yes | 68 |
| 20 | 0.2411 | 2.9e-13 | 3.4e-13 | 3.6e-15 | 3.2e-13 | 6.3e-14 | local | yes | 77 |
| 16 | 0.4 | 3.0e-13 | 2.8e-13 | 8.0e-15 | 7.9e-14 | 6.8e-14 | local | yes | 70 |
| 18 | 0.4 | 3.7e-13 | 2.4e-13 | 8.9e-16 | 6.0e-13 | 7.9e-14 | local | yes | 79 |
| 20 | 0.4 | 3.7e-13 | 2.9e-13 | 5.3e-15 | 4.3e-13 | 6.8e-14 | local | yes | 86 |

*Table 5a.* Absolute errors of the DMRG energies, the relative error of the gap (|ΔE1|/gap for N ≤ 14, where the E0 error is negligible next to it), and the overlap deficit of the ground state, against exact diagonalisation.

**No exact reference: two independent runs from different random starts (N ≥ 24)**

| N | J2 | E0 | E1 | run-to-run \|ΔE0\| | \|ΔE1\| | \|ΔE_top\| | 1 − overlap(ψ0, ψ0′) | variance of ψ0 (a / b) | variance of ψ1 (a / b) | converged | max bond |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 24 | 0 | -10.45378576 | -10.30083406 | 3.0e-14 | 3.6e-15 | 2.9e-14 | 4.7e-14 | 3.2e-12 / 2.6e-12 | 2.4e-12 / 2.4e-12 | yes | 91 |
| 28 | 0 | -12.22544055 | -12.09242500 | 1.4e-14 | 3.0e-14 | 8.8e-14 | 1.8e-14 | 3.2e-12 / 3.5e-12 | 2.9e-12 / 2.8e-12 | yes | 100 |
| 32 | 0 | -13.99731562 | -13.87957582 | 3.2e-14 | 1.6e-14 | 1.2e-13 | 2.4e-14 | 6.3e-12 / 5.6e-12 | 3.4e-12 / 3.7e-12 | yes | 115 |
| 24 | 0.2411 | -9.55622125 | -9.41194333 | 7.1e-15 | 2.3e-14 | 1.5e-14 | 5.2e-14 | 2.4e-12 / 2.4e-12 | 2.4e-12 / 2.3e-12 | yes | 88 |
| 28 | 0.2411 | -11.16322989 | -11.03847669 | 2.5e-14 | 4.4e-14 | 2.8e-14 | 3.5e-14 | 2.2e-12 / 2.6e-12 | 3.3e-12 / 3.8e-12 | yes | 101 |
| 32 | 0.2411 | -12.77043249 | -12.66054892 | 4.6e-14 | 1.8e-14 | 7.3e-14 | 3.0e-14 | 5.2e-12 / 4.1e-12 | 3.6e-12 / 2.9e-12 | yes | 117 |
| 24 | 0.4 | -9.10932537 | -8.93453579 | 2.0e-14 | 2.7e-14 | 3.6e-14 | 2.3e-14 | 3.0e-12 / 2.8e-12 | 2.3e-12 / 2.6e-12 | yes | 103 |
| 28 | 0.4 | -10.63076570 | -10.47463797 | 2.0e-14 | 7.1e-15 | 6.2e-14 | 4.3e-14 | 4.0e-12 / 3.4e-12 | 3.0e-12 / 2.9e-12 | yes | 113 |
| 32 | 0.4 | -12.15223705 | -12.01039236 | 5.7e-14 | 7.1e-15 | 9.2e-14 | 5.2e-14 | 5.5e-12 / 4.3e-12 | 3.3e-12 / 3.2e-12 | yes | 128 |

*Table 5b.* Agreement between two independent DMRG runs and the energy variance ⟨H²⟩ − ⟨H⟩² of each state. A state with variance s² has an exact eigenvalue within s of its energy; run-to-run agreement at 1e-13 or better together with variances of order 1e-12 is strong evidence of convergence, but it is not a proof that E1 is the second-lowest level.

**The ground-state weight γ used for N > 14 against the dense value**

| N | J2 | γ, MPS contraction (no ED, no statevector) | γ, dense overlap with exact ground state | difference |
|---|---|---|---|---|
| 16 | 0 | 0.606802310650 | 0.606802312131 | -1.5e-09 |
| 18 | 0 | 0.554637097041 | 0.554637099565 | -2.5e-09 |
| 16 | 0.2411 | 0.736784131818 | 0.736784132514 | -7.0e-10 |
| 18 | 0.2411 | 0.694543127658 | 0.694543127932 | -2.7e-10 |
| 16 | 0.4 | 0.916938666774 | 0.916938666766 | +7.8e-12 |
| 18 | 0.4 | 0.900035947685 | 0.900035947967 | -2.8e-10 |

*Table 5c.* γ is the overlap of the DMRG ground state with the trial circuit. The N > 14 path computes it by MPS contraction; the dense value is for checking only.

## 6. The certified bandwidth

The scaling H_s = (H − E0)/W needs W ≥ E_top − E0 so that the spectrum stays inside [0, 1]. A certified W must not rely on DMRG's E_top (DMRG on −H is variational, so E_top(DMRG) ≤ E_top and W(DMRG) ≤ W_exact: the unsafe direction). Weyl's inequality gives λmax(Σ_b h_b) ≤ Σ_b λmax(h_b); a bond J S·S has λmax = J/4 (J ≥ 0), so E_top ≤ (N−1)J1/4 + (N−2)J2/4, and the fully polarized state attains it. Hence W = (N−1)J1/4 + (N−2)J2/4 − E0 equals E_top − E0 (to rounding) with no DMRG input for E_top. An earlier version used the looser Pauli bound W = Σ|c_β| − E0 = 3(Σ_b λmax) − E0, which is 1.6–1.9× larger; Table 5 compares them.

| N | J2 | W (Weyl) = E_top − E0 | \|W − (E_top−E0)_DMRG\| | W (Pauli-ℓ1, earlier) | ℓ1 / Weyl |
|---|---|---|---|---|---|
| 4 | 0 | 2.366025 | 0.0e+00 | 3.8660 | 1.63 |
| 6 | 0 | 3.743577 | 0.0e+00 | 6.2436 | 1.67 |
| 8 | 0 | 5.124933 | 1.8e-15 | 8.6249 | 1.68 |
| 10 | 0 | 6.508035 | 0.0e+00 | 11.0080 | 1.69 |
| 12 | 0 | 7.892091 | 1.8e-15 | 13.3921 | 1.70 |
| 14 | 0 | 9.276725 | 3.6e-15 | 15.7767 | 1.70 |
| 16 | 0 | 10.661737 | 0.0e+00 | 18.1617 | 1.70 |
| 18 | 0 | 12.047011 | 1.6e-14 | 20.5470 | 1.71 |
| 20 | 0 | 13.432473 | 3.6e-15 | 22.9325 | 1.71 |
| 24 | 0 | 16.203786 | 5.3e-14 | 27.7038 | 1.71 |
| 28 | 0 | 18.975441 | 4.3e-14 | 32.4754 | 1.71 |
| 32 | 0 | 21.747316 | 1.8e-13 | 37.2473 | 1.71 |
| 4 | 0.2411 | 2.409291 | 0.0e+00 | 4.1504 | 1.72 |
| 6 | 0.2411 | 3.826447 | 0.0e+00 | 6.8086 | 1.78 |
| 8 | 0.2411 | 5.246651 | 1.8e-15 | 9.4700 | 1.80 |
| 10 | 0.2411 | 6.668328 | 5.3e-15 | 12.1327 | 1.82 |
| 12 | 0.2411 | 8.090827 | 0.0e+00 | 14.7963 | 1.83 |
| 14 | 0.2411 | 9.513831 | 7.1e-15 | 17.4604 | 1.84 |
| 16 | 0.2411 | 10.937166 | 0.0e+00 | 20.1249 | 1.84 |
| 18 | 0.2411 | 12.360732 | 0.0e+00 | 22.7895 | 1.84 |
| 20 | 0.2411 | 13.784465 | 0.0e+00 | 25.4544 | 1.85 |
| 24 | 0.2411 | 16.632271 | 7.1e-15 | 30.7844 | 1.85 |
| 28 | 0.2411 | 19.480380 | 8.5e-14 | 36.1147 | 1.85 |
| 32 | 0.2411 | 22.328682 | 1.7e-13 | 41.4452 | 1.86 |
| 4 | 0.4 | 2.456776 | 8.9e-16 | 4.3568 | 1.77 |
| 6 | 0.4 | 3.915535 | 0.0e+00 | 7.2155 | 1.84 |
| 8 | 0.4 | 5.375191 | 6.2e-15 | 10.0752 | 1.87 |
| 10 | 0.4 | 6.835289 | 0.0e+00 | 12.9353 | 1.89 |
| 12 | 0.4 | 8.295626 | 5.3e-15 | 15.7956 | 1.90 |
| 14 | 0.4 | 9.756103 | 0.0e+00 | 18.6561 | 1.91 |
| 16 | 0.4 | 11.216666 | 1.1e-14 | 21.5167 | 1.92 |
| 18 | 0.4 | 12.677284 | 0.0e+00 | 24.3773 | 1.92 |
| 20 | 0.4 | 14.137941 | 3.6e-15 | 27.2379 | 1.93 |
| 24 | 0.4 | 17.059325 | 7.1e-15 | 32.9593 | 1.93 |
| 28 | 0.4 | 19.980766 | 2.5e-14 | 38.6808 | 1.94 |
| 32 | 0.4 | 22.902237 | 4.6e-14 | 44.4022 | 1.94 |

*Table 5.* The Weyl bandwidth agrees with the DMRG E_top − E0 to rounding (≤ 1e-13), because E_top is the all-up energy. **W cancels in the step count:** α_scaled = α/W² and T_scaled = xπW/g (g the unscaled gap), so n ∝ α_scaled T² ∝ α/g² and T (phys) = T_scaled/W = xπ/g do not depend on W. W only changes the *window*: [Δ, 1] is wider than the occupied spectrum when W is loose, which can make the certified η slightly larger.

Measured effect of switching from the Pauli-ℓ1 bandwidth to the Weyl bandwidth, over 59 matching (N, J2, ε) guaranteed designs: n changed by a factor 0.63–1.14 (median 0.99), T (phys) by 0.67–1.50.

## 7. Larger J1–J2 chains (N = 14–32): resource costs

N ≥ 16 is **DMRG/MPS only**: γ is the DMRG-MPS overlap of the trial circuit, the gap is DMRG's, W is certified, α is the exact commutator sum, and there is no ED, no statevector and no dense operator anywhere. The designs are therefore conditional certificates: rigorous *given* γ, the gap and the shift (see `CERTIFICATION.tex`). Circuits are **not simulated**; their counts are n × the transpiled per-step cost plus the trial circuit, validated against full transpiles at n = 48 (Section 9). N = 14 (statevector-simulated, ED cross-checked) is included for continuity.

**Guaranteed designs, ε = 0.01**

| N | J2 | γ | Δ | W | α | pulses | n (steps) | T (scaled) | T (phys) | CX trial | CX filter | CX total | CX filter (direct) | RZ non-Cliff. (level 3) | RZ non-Cliff. (direct) | depth | P_succ lb | ε bound | CX total (line) | DMRG converged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 14 | 0 | 0.6609 | 0.0265 | 9.28 | 0.1046 | 8 | 16,975 | 118.4 | 12.76 | 39 | 1,544,725 | 1,544,764 | 2,648,100 | 2,020,281 | 679,008 | 6,433,816 | 0.595 | 1.00e-02 | 2,155,864 | yes |
| 16 | 0 | 0.6047 | 0.0206 | 10.66 | 0.0924 | 8 | 26,578 | 152.8 | 14.33 | 45 | 2,790,690 | 2,790,735 | 4,784,040 | 3,641,482 | 1,222,596 | 11,614,921 | 0.526 | 9.99e-03 | 3,907,011 | yes |
| 18 | 0 | 0.5547 | 0.0164 | 12.05 | 0.0827 | 8 | 42,194 | 191.4 | 15.89 | 51 | 5,021,086 | 5,021,137 | 8,607,576 | 6,540,404 | 2,194,096 | 20,886,409 | 0.454 | 9.48e-03 | 7,046,449 | yes |
| 20 | 0 | 0.5139 | 0.0134 | 13.43 | 0.0748 | 8 | 60,428 | 234.3 | 17.44 | 57 | 8,036,924 | 8,036,981 | 13,777,584 | 10,454,414 | 3,504,832 | 33,417,107 | 0.407 | 9.51e-03 | 11,300,093 | yes |
| 24 | 0 | 0.4323 | 0.0094 | 16.20 | 0.0628 | 8 | 168,243 | 332.8 | 20.54 | 69 | 27,087,123 | 27,087,192 | 46,435,068 | 35,163,231 | 11,777,018 | 112,555,078 | 0.332 | 9.67e-03 | 39,088,526 | yes |
| 14 | 0.2411 | 0.7810 | 0.0249 | 9.51 | 0.1651 | 8 | 22,902 | 126.1 | 13.25 | 39 | 4,007,850 | 4,007,889 | 6,870,600 | 5,199,007 | 1,740,560 | 16,650,045 | 0.694 | 9.99e-03 | 7,756,183 | yes |
| 16 | 0.2411 | 0.7365 | 0.0192 | 10.94 | 0.1461 | 8 | 37,598 | 163.6 | 14.96 | 45 | 7,632,394 | 7,632,439 | 13,084,104 | 9,888,566 | 3,308,632 | 31,695,449 | 0.629 | 9.36e-03 | 14,838,722 | yes |
| 18 | 0.2411 | 0.6944 | 0.0153 | 12.36 | 0.1310 | 8 | 56,486 | 206.0 | 16.66 | 51 | 13,048,266 | 13,048,317 | 22,368,456 | 16,889,642 | 5,648,608 | 54,170,453 | 0.562 | 9.36e-03 | 25,456,408 | yes |
| 20 | 0.2411 | 0.6528 | 0.0124 | 13.78 | 0.1186 | 8 | 89,033 | 253.2 | 18.37 | 57 | 23,059,547 | 23,059,604 | 39,530,652 | 29,826,418 | 9,971,704 | 95,710,898 | 0.583 | 1.00e-02 | 45,110,110 | yes |
| 24 | 0.2411 | 0.5757 | 0.0087 | 16.63 | 0.0998 | 8 | 174,959 | 271.6 | 16.33 | 69 | 55,112,085 | 55,112,154 | 94,477,860 | 71,208,752 | 23,794,432 | 228,671,924 | 0.313 | 9.98e-03 | 215,316,278 | yes |
| 14 | 0.4 | 0.9334 | 0.0269 | 9.76 | 0.2008 | 6 | 10,571 | 70.0 | 7.18 | 39 | 1,849,925 | 1,849,964 | 3,171,300 | 2,399,870 | 803,402 | 7,685,408 | 0.689 | 1.00e-02 | 3,580,084 | yes |
| 16 | 0.4 | 0.9170 | 0.0211 | 11.22 | 0.1778 | 6 | 17,333 | 89.2 | 7.95 | 45 | 3,518,599 | 3,518,644 | 6,031,884 | 4,558,870 | 1,525,310 | 14,612,054 | 0.587 | 9.53e-03 | 6,840,802 | yes |
| 18 | 0.4 | 0.9000 | 0.0171 | 12.68 | 0.1595 | 8 | 28,325 | 137.8 | 10.87 | 51 | 6,543,075 | 6,543,126 | 11,216,700 | 8,469,507 | 2,832,508 | 27,164,054 | 0.602 | 1.00e-02 | 12,765,184 | yes |
| 20 | 0.4 | 0.8833 | 0.0142 | 14.14 | 0.1446 | 8 | 42,215 | 166.4 | 11.77 | 57 | 10,933,685 | 10,933,742 | 18,743,460 | 14,142,387 | 4,728,088 | 45,381,548 | 0.717 | 9.99e-03 | 21,388,990 | yes |
| 24 | 0.4 | 0.8489 | 0.0102 | 17.06 | 0.1217 | 8 | 73,076 | 230.0 | 13.48 | 69 | 23,018,940 | 23,019,009 | 39,461,040 | 29,742,372 | 9,938,344 | 95,510,843 | 0.619 | 1.00e-02 | 89,932,266 | yes |
| 28 | 0.4 | 0.8154 | 0.0078 | 19.98 | 0.1051 | 8 | 115,832 | 301.5 | 15.09 | 81 | 42,973,672 | 42,973,753 | 73,669,152 | 55,484,042 | 18,533,128 | 178,266,047 | 0.544 | 9.99e-03 | 231,548,249 | yes |
| 32 | 0.4 | 0.7822 | 0.0062 | 22.90 | 0.0924 | 8 | 224,514 | 380.4 | 16.61 | 93 | 95,867,478 | 95,867,571 | 164,344,248 | 123,707,799 | 41,310,584 | 397,614,980 | 0.625 | 9.99e-03 | 570,415,329 | yes |

**Guaranteed designs, ε = 0.1**

| N | J2 | γ | Δ | W | α | pulses | n (steps) | T (scaled) | T (phys) | CX trial | CX filter | CX total | CX filter (direct) | RZ non-Cliff. (level 3) | RZ non-Cliff. (direct) | depth | P_succ lb | ε bound | CX total (line) | DMRG converged |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 0 | 0.6047 | 0.0206 | 10.66 | 0.0924 | 6 | 3,212 | 91.7 | 8.60 | 45 | 337,260 | 337,305 | 578,160 | 440,338 | 147,758 | 1,403,979 | 0.464 | 9.99e-02 | 472,209 | yes |
| 18 | 0 | 0.5547 | 0.0164 | 12.05 | 0.0827 | 6 | 4,946 | 114.8 | 9.53 | 51 | 588,574 | 588,625 | 1,008,984 | 766,962 | 257,198 | 2,448,649 | 0.388 | 9.99e-02 | 826,033 | yes |
| 20 | 0 | 0.5139 | 0.0134 | 13.43 | 0.0748 | 6 | 7,479 | 117.2 | 8.72 | 57 | 994,707 | 994,764 | 1,705,212 | 1,294,235 | 433,788 | 4,136,310 | 0.330 | 9.99e-02 | 1,398,630 | yes |
| 24 | 0 | 0.4323 | 0.0094 | 16.20 | 0.0628 | 8 | 18,710 | 249.6 | 15.40 | 69 | 3,012,310 | 3,012,379 | 5,163,960 | 3,910,834 | 1,309,708 | 12,517,501 | 0.353 | 9.95e-02 | 4,347,026 | yes |
| 28 | 0 | 0.3652 | 0.0070 | 18.98 | 0.0542 | 8 | 35,021 | 336.1 | 17.71 | 81 | 6,618,969 | 6,619,050 | 11,346,804 | 8,580,662 | 2,871,730 | 27,492,084 | 0.270 | 9.47e-02 | 21,794,817 | yes |
| 32 | 0 | 0.3081 | 0.0054 | 21.75 | 0.0476 | 8 | 56,303 | 435.2 | 20.01 | 93 | 12,217,751 | 12,217,844 | 20,944,716 | 15,821,729 | 5,292,490 | 50,729,690 | 0.205 | 1.00e-01 | 41,776,919 | yes |
| 16 | 0.2411 | 0.7365 | 0.0192 | 10.94 | 0.1461 | 6 | 3,769 | 65.4 | 5.98 | 45 | 765,107 | 765,152 | 1,311,612 | 991,537 | 331,678 | 3,177,602 | 0.479 | 1.00e-01 | 1,487,544 | yes |
| 18 | 0.2411 | 0.6944 | 0.0153 | 12.36 | 0.1310 | 6 | 5,803 | 103.0 | 8.33 | 51 | 1,340,493 | 1,340,544 | 2,297,988 | 1,735,423 | 580,306 | 5,565,456 | 0.456 | 9.99e-02 | 2,615,270 | yes |
| 20 | 0.2411 | 0.6528 | 0.0124 | 13.78 | 0.1186 | 6 | 8,594 | 126.6 | 9.18 | 57 | 2,225,846 | 2,225,903 | 3,815,736 | 2,879,351 | 962,534 | 9,238,973 | 0.390 | 9.99e-02 | 4,354,350 | yes |
| 24 | 0.2411 | 0.5757 | 0.0087 | 16.63 | 0.0998 | 8 | 20,722 | 217.3 | 13.06 | 69 | 6,527,430 | 6,527,499 | 11,189,880 | 8,434,293 | 2,818,200 | 27,084,165 | 0.413 | 9.99e-02 | 25,501,944 | yes |
| 28 | 0.2411 | 0.5046 | 0.0064 | 19.48 | 0.0861 | 8 | 38,646 | 245.3 | 12.59 | 81 | 14,337,666 | 14,337,747 | 24,578,856 | 18,511,946 | 6,183,368 | 59,476,793 | 0.315 | 1.00e-01 | 77,253,435 | yes |
| 32 | 0.2411 | 0.4417 | 0.0049 | 22.33 | 0.0757 | 8 | 80,733 | 383.0 | 17.15 | 93 | 34,472,991 | 34,473,084 | 59,096,556 | 44,484,467 | 14,854,880 | 142,978,830 | 0.342 | 9.90e-02 | 205,115,735 | yes |
| 16 | 0.4 | 0.9170 | 0.0211 | 11.22 | 0.1778 | 4 | 2,135 | 59.4 | 5.30 | 45 | 433,405 | 433,450 | 742,980 | 561,794 | 187,884 | 1,800,140 | 0.759 | 9.46e-02 | 842,658 | yes |
| 18 | 0.4 | 0.9000 | 0.0171 | 12.68 | 0.1595 | 6 | 2,925 | 73.5 | 5.80 | 51 | 675,675 | 675,726 | 1,158,300 | 874,905 | 292,506 | 2,805,454 | 0.688 | 1.00e-01 | 1,318,251 | yes |
| 20 | 0.4 | 0.8833 | 0.0142 | 14.14 | 0.1446 | 6 | 4,030 | 88.7 | 6.28 | 57 | 1,043,770 | 1,043,827 | 1,789,320 | 1,350,410 | 451,366 | 4,332,673 | 0.633 | 1.00e-01 | 2,041,924 | yes |
| 24 | 0.4 | 0.8489 | 0.0102 | 17.06 | 0.1217 | 6 | 6,966 | 122.6 | 7.19 | 69 | 2,194,290 | 2,194,359 | 3,761,640 | 2,835,600 | 947,382 | 9,105,073 | 0.549 | 9.99e-02 | 8,572,893 | yes |
| 28 | 0.4 | 0.8154 | 0.0078 | 19.98 | 0.1051 | 6 | 11,074 | 160.8 | 8.05 | 81 | 4,108,454 | 4,108,535 | 7,043,064 | 5,304,958 | 1,771,846 | 17,043,485 | 0.484 | 9.94e-02 | 22,137,007 | yes |
| 32 | 0.4 | 0.7822 | 0.0062 | 22.90 | 0.0924 | 6 | 21,424 | 202.9 | 8.86 | 93 | 9,148,048 | 9,148,141 | 15,682,368 | 11,805,207 | 3,942,022 | 37,942,590 | 0.555 | 9.96e-02 | 54,431,336 | yes |

**Pulse structure of the guaranteed designs (ε = 10⁻²)**

| N | J2 | kᵢ (guaranteed) | largest block kmax | CX in largest block |
|---|---|---|---|---|
| 14 | 0 | [106, 307, 543, 917, 1606, 2891, 4199, 6406] | 6406 | 582,946 |
| 16 | 0 | [166, 447, 895, 1477, 2579, 4480, 6459, 10075] | 10075 | 1,057,875 |
| 18 | 0 | [533, 1110, 1311, 2250, 3965, 6926, 10004, 16095] | 16095 | 1,915,305 |
| 20 | 0 | [575, 1120, 2068, 3348, 5802, 10072, 14594, 22849] | 22849 | 3,038,917 |
| 24 | 0 | [1492, 4232, 7046, 8144, 14111, 28322, 37845, 67051] | 67051 | 10,795,211 |
| 14 | 0.2411 | [143, 413, 769, 1260, 2167, 3832, 5619, 8699] | 8699 | 1,522,325 |
| 16 | 0.2411 | [304, 609, 1266, 2092, 3665, 6310, 9094, 14258] | 14258 | 2,894,374 |
| 18 | 0.2411 | [532, 1014, 1906, 3131, 5442, 9306, 13660, 21495] | 21495 | 4,965,345 |
| 20 | 0.2411 | [775, 1602, 2889, 4843, 8323, 14899, 22099, 33603] | 33603 | 8,703,177 |
| 24 | 0.2411 | [1375, 3978, 7547, 9236, 16024, 26438, 40961, 69400] | 69400 | 21,861,000 |
| 14 | 0.4 | [88, 306, 706, 1430, 2756, 5285] | 5285 | 924,875 |
| 16 | 0.4 | [358, 674, 1218, 2285, 4414, 8384] | 8384 | 1,701,952 |
| 18 | 0.4 | [177, 373, 882, 1509, 2456, 4381, 7535, 11012] | 11012 | 2,543,772 |
| 20 | 0.4 | [264, 314, 852, 1841, 3541, 5936, 10588, 18879] | 18879 | 4,889,661 |
| 24 | 0.4 | [457, 514, 1895, 3248, 6172, 10365, 17575, 32850] | 32850 | 10,347,750 |
| 28 | 0.4 | [724, 1136, 2436, 5581, 10251, 20717, 32528, 42459] | 42459 | 15,752,289 |
| 32 | 0.4 | [1403, 1867, 4587, 9837, 18761, 31824, 55755, 100480] | 100480 | 42,904,960 |

**Measured (statevector-simulated) points, N = 12 and 14**

| N | J2 | design | n | pulses | T (phys) | CX total | RZ non-Cliff. (direct) | measured 1−F | P_succ | ε bound (rigorous) |
|---|---|---|---|---|---|---|---|---|---|---|
| 12 | 0 | guaranteed | 8,800 | 8 | 11.18 | 677,633 | 299,208 | 3.29e-04 | 0.600 | 9.99e-03 |
| 12 | 0 | empirical | 29 | 6 | 11.18 | 2,266 | 992 | 5.75e-03 | 0.240 | 1.00e+00 |
| 14 | 0 | guaranteed | 16,975 | 8 | 12.76 | 1,544,764 | 679,008 | 6.34e-04 | 0.595 | 1.00e-02 |
| 14 | 0 | empirical | 39 | 7 | 12.76 | 3,588 | 1,567 | 3.45e-03 | 0.275 | 1.00e+00 |
| 12 | 0.2411 | guaranteed | 10,960 | 6 | 8.66 | 1,611,153 | 701,446 | 4.09e-04 | 0.575 | 1.00e-02 |
| 12 | 0.2411 | empirical | 25 | 6 | 8.66 | 3,708 | 1,606 | 7.80e-03 | 0.480 | 1.00e+00 |
| 14 | 0.2411 | guaranteed | 22,902 | 8 | 13.25 | 4,007,889 | 1,740,560 | 3.83e-04 | 0.694 | 9.99e-03 |
| 14 | 0.2411 | empirical | 41 | 7 | 13.25 | 7,214 | 3,123 | 7.83e-03 | 0.387 | 1.00e+00 |
| 12 | 0.4 | guaranteed | 5,320 | 6 | 5.31 | 782,073 | 340,486 | 4.80e-04 | 0.562 | 1.00e-02 |
| 12 | 0.4 | empirical | 26 | 5 | 5.31 | 3,855 | 1,669 | 9.35e-03 | 0.501 | 1.00e+00 |
| 14 | 0.4 | guaranteed | 10,571 | 6 | 7.18 | 1,849,964 | 803,402 | 4.64e-04 | 0.689 | 1.00e-02 |
| 14 | 0.4 | empirical | 38 | 5 | 7.18 | 6,689 | 2,893 | 9.75e-03 | 0.610 | 1.00e+00 |

*Empirical* rows are the oracle lower envelope; *guaranteed* rows show how far below ε the measured error sits.

### Scaling with N

| quantity | J2 | source | N range | p in y ~ N^p | R² | points |
|---|---|---|---|---|---|---|
| n (steps) | 0 | dmrg | 6–24 | 3.99 | 0.997 | 9 |
| T (scaled) | 0 | dmrg | 6–24 | 2.17 | 0.994 | 9 |
| CX per step | 0 | dmrg | 6–24 | 1.10 | 1.000 | 9 |
| CX total | 0 | dmrg | 6–24 | 5.08 | 0.999 | 9 |
| CX (direct) | 0 | dmrg | 6–24 | 5.08 | 0.999 | 9 |
| RZ non-Cliff. (level 3) | 0 | dmrg | 6–24 | 5.06 | 0.999 | 9 |
| RZ non-Cliff. (direct) | 0 | dmrg | 6–24 | 5.05 | 0.999 | 9 |
| depth total | 0 | dmrg | 6–24 | 5.08 | 0.999 | 9 |
| n (steps) | 0.2411 | dmrg | 6–24 | 4.07 | 0.998 | 9 |
| T (scaled) | 0.2411 | dmrg | 6–24 | 2.31 | 0.978 | 9 |
| CX per step | 0.2411 | dmrg | 6–24 | 1.16 | 0.999 | 9 |
| CX total | 0.2411 | dmrg | 6–24 | 5.22 | 0.998 | 9 |
| CX (direct) | 0.2411 | dmrg | 6–24 | 5.22 | 0.998 | 9 |
| RZ non-Cliff. (level 3) | 0.2411 | dmrg | 6–24 | 5.21 | 0.998 | 9 |
| RZ non-Cliff. (direct) | 0.2411 | dmrg | 6–24 | 5.20 | 0.998 | 9 |
| depth total | 0.2411 | dmrg | 6–24 | 5.22 | 0.998 | 9 |
| n (steps) | 0.4 | dmrg | 6–32 | 3.78 | 0.998 | 11 |
| T (scaled) | 0.4 | dmrg | 6–32 | 2.30 | 0.992 | 11 |
| CX per step | 0.4 | dmrg | 6–32 | 1.13 | 0.999 | 11 |
| CX total | 0.4 | dmrg | 6–32 | 4.92 | 0.999 | 11 |
| CX (direct) | 0.4 | dmrg | 6–32 | 4.92 | 0.999 | 11 |
| RZ non-Cliff. (level 3) | 0.4 | dmrg | 6–32 | 4.91 | 0.999 | 11 |
| RZ non-Cliff. (direct) | 0.4 | dmrg | 6–32 | 4.90 | 0.999 | 11 |
| depth total | 0.4 | dmrg | 6–32 | 4.91 | 0.999 | 11 |

*Table 7.* Power-law fits of guaranteed-design quantities (ε = 10⁻²) over points with converged DMRG and a certified design from the floor solver (fallback-builder points, marked † in Section 7, are excluded). Few points: indicative, not asymptotic.

## 8. Dependence on the target ε

**N = 6, J2 = 0**

| ε | design | n | pulses | T (phys) | CX total | RZ non-Cliff. (direct) | depth | measured 1−F | ε bound | P_succ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1e-01 | guaranteed | 67 | 4 | 2.56 | 2,360 | 1,076 | 9,960 | 7.23e-03 | 9.86e-02 | 0.701 |
| 1e-01 | empirical | 2 | 2 | 2.56 | 85 | 34 | 405 | 4.72e-02 | 1.00e+00 | 0.860 |
| 3e-02 | guaranteed | 237 | 4 | 3.83 | 8,310 | 3,796 | 34,950 | 4.86e-04 | 2.78e-02 | 0.678 |
| 3e-02 | empirical | 4 | 3 | 3.83 | 155 | 67 | 699 | 9.92e-03 | 1.00e+00 | 0.859 |
| 1e-02 | guaranteed | 591 | 4 | 4.79 | 20,700 | 9,460 | 86,988 | 6.87e-04 | 9.99e-03 | 0.717 |
| 1e-02 | empirical | 5 | 4 | 4.79 | 190 | 84 | 846 | 6.96e-03 | 1.00e+00 | 0.850 |
| 3e-03 | guaranteed | 1,809 | 6 | 6.39 | 63,330 | 28,950 | 266,034 | 6.99e-05 | 3.00e-03 | 0.804 |
| 3e-03 | empirical | 10 | 5 | 6.39 | 365 | 165 | 1,581 | 1.44e-03 | 1.00e+00 | 0.662 |
| 1e-03 | guaranteed | 4,332 | 6 | 7.99 | 151,635 | 69,318 | 636,915 | 2.96e-05 | 1.00e-03 | 0.801 |
| 1e-03 | empirical | 13 | 6 | 7.99 | 470 | 214 | 2,022 | 7.99e-04 | 1.00e+00 | 0.645 |

**N = 8, J2 = 0**

| ε | design | n | pulses | T (phys) | CX total | RZ non-Cliff. (direct) | depth | measured 1−F | ε bound | P_succ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1e-01 | guaranteed | 179 | 4 | 3.20 | 8,792 | 3,942 | 36,854 | 1.23e-02 | 9.92e-02 | 0.539 |
| 1e-01 | empirical | 3 | 2 | 3.20 | 168 | 68 | 774 | 2.79e-02 | 1.00e+00 | 0.749 |
| 3e-02 | guaranteed | 565 | 4 | 4.00 | 27,706 | 12,434 | 115,984 | 1.47e-03 | 2.83e-02 | 0.476 |
| 3e-02 | empirical | 6 | 3 | 4.00 | 315 | 135 | 1,389 | 4.71e-03 | 1.00e+00 | 0.813 |
| 1e-02 | guaranteed | 1,819 | 6 | 6.00 | 89,152 | 40,024 | 373,054 | 1.25e-03 | 9.99e-03 | 0.631 |
| 1e-02 | empirical | 11 | 4 | 6.00 | 560 | 246 | 2,414 | 5.25e-03 | 1.00e+00 | 0.661 |
| 3e-03 | guaranteed | 5,045 | 6 | 8.00 | 247,226 | 110,996 | 1,034,384 | 9.44e-05 | 2.80e-03 | 0.630 |
| 3e-03 | empirical | 15 | 6 | 8.00 | 756 | 336 | 3,234 | 1.38e-03 | 1.00e+00 | 0.548 |
| 1e-03 | guaranteed | 12,070 | 8 | 10.00 | 591,451 | 265,548 | 2,474,509 | 5.79e-05 | 1.00e-03 | 0.732 |
| 1e-03 | empirical | 24 | 7 | 10.00 | 1,197 | 535 | 5,079 | 9.01e-04 | 1.00e+00 | 0.652 |

**N = 6, J2 = 0.4**

| ε | design | n | pulses | T (phys) | CX total | RZ non-Cliff. (direct) | depth | measured 1−F | ε bound | P_succ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1e-01 | guaranteed | 79 | 4 | 2.50 | 4,992 | 2,216 | 20,888 | 1.19e-03 | 9.80e-02 | 0.976 |
| 1e-01 | empirical | 3 | 2 | 2.50 | 204 | 86 | 900 | 5.15e-02 | 1.00e+00 | 0.930 |
| 3e-02 | guaranteed | 197 | 4 | 2.50 | 12,426 | 5,520 | 51,922 | 1.75e-03 | 2.99e-02 | 0.990 |
| 3e-02 | empirical | 4 | 2 | 2.50 | 267 | 114 | 1,163 | 2.74e-02 | 1.00e+00 | 0.948 |
| 1e-02 | guaranteed | 394 | 4 | 2.50 | 24,837 | 11,036 | 103,733 | 4.13e-04 | 9.98e-03 | 0.740 |
| 1e-02 | empirical | 6 | 3 | 2.50 | 393 | 171 | 1,689 | 6.01e-03 | 1.00e+00 | 0.697 |
| 3e-03 | guaranteed | 1,322 | 4 | 3.75 | 83,301 | 37,020 | 347,797 | 1.99e-04 | 3.00e-03 | 0.733 |
| 3e-03 | empirical | 15 | 4 | 3.75 | 960 | 424 | 4,056 | 2.76e-03 | 1.00e+00 | 0.710 |
| 1e-03 | guaranteed | 3,033 | 4 | 3.75 | 191,094 | 84,928 | 797,790 | 1.39e-04 | 1.00e-03 | 0.630 |
| 1e-03 | empirical | 26 | 4 | 3.75 | 1,653 | 732 | 6,949 | 9.63e-04 | 1.00e+00 | 0.617 |

## 9. Checks and classical compute cost

**Per-step × n estimate against a full transpile of the real circuit (n = 48)**

| source | N | J2 | n | estimate CX (per-step × n) | full-transpile CX | CX error | estimate depth | full-transpile depth |
|---|---|---|---|---|---|---|---|---|
| dmrg | 4 | 0 | 48 | 1,008 | 1,008 | +0.0% | 4,272 | 4,196 |
| dmrg | 4 | 0.2411 | 48 | 1,680 | 1,680 | +0.0% | 7,056 | 6,980 |
| dmrg | 4 | 0.4 | 48 | 1,680 | 1,680 | +0.0% | 7,056 | 6,980 |
| dmrg | 6 | 0 | 48 | 1,680 | 1,680 | +0.0% | 7,056 | 6,980 |
| dmrg | 6 | 0.2411 | 48 | 3,024 | 3,024 | +0.0% | 12,624 | 12,548 |
| dmrg | 6 | 0.4 | 48 | 3,024 | 3,024 | +0.0% | 12,624 | 12,548 |
| dmrg | 8 | 0 | 48 | 2,352 | 2,352 | +0.0% | 9,840 | 9,769 |
| dmrg | 8 | 0.2411 | 48 | 4,368 | 4,368 | +0.0% | 18,192 | 18,121 |
| dmrg | 8 | 0.4 | 48 | 4,368 | 4,368 | +0.0% | 18,192 | 18,116 |
| dmrg | 10 | 0 | 48 | 3,024 | 3,024 | +0.0% | 12,624 | 12,558 |
| dmrg | 10 | 0.2411 | 48 | 5,712 | 5,712 | +0.0% | 23,760 | 23,694 |
| dmrg | 10 | 0.4 | 48 | 5,712 | 5,712 | +0.0% | 23,760 | 23,684 |
| dmrg | 12 | 0 | 48 | 3,696 | 3,696 | +0.0% | 15,408 | 15,347 |
| dmrg | 12 | 0.2411 | 48 | 7,056 | 7,056 | +0.0% | 29,328 | 29,262 |
| dmrg | 12 | 0.4 | 48 | 7,056 | 7,056 | +0.0% | 29,328 | 29,257 |
| dmrg | 14 | 0 | 48 | 4,368 | 4,368 | +0.0% | 18,192 | 18,131 |
| dmrg | 14 | 0.2411 | 48 | 8,400 | 8,400 | +0.0% | 34,896 | 34,835 |
| dmrg | 14 | 0.4 | 48 | 8,400 | 8,400 | +0.0% | 34,896 | 34,825 |
| dmrg | 16 | 0 | 48 | 5,040 | 5,040 | +0.0% | 20,976 | 20,915 |
| dmrg | 16 | 0.2411 | 48 | 9,744 | 9,744 | +0.0% | 40,464 | 40,403 |
| dmrg | 16 | 0.4 | 48 | 9,744 | 9,744 | +0.0% | 40,464 | 40,398 |
| dmrg | 18 | 0 | 48 | 5,712 | 5,712 | +0.0% | 23,760 | 23,704 |
| dmrg | 18 | 0.2411 | 48 | 11,088 | 11,088 | +0.0% | 46,032 | 45,971 |
| dmrg | 18 | 0.4 | 48 | 11,088 | 11,088 | +0.0% | 46,032 | 45,971 |
| dmrg | 20 | 0 | 48 | 6,384 | 6,384 | +0.0% | 26,544 | 26,483 |
| dmrg | 20 | 0.2411 | 48 | 12,432 | 12,432 | +0.0% | 51,600 | 51,539 |
| dmrg | 20 | 0.4 | 48 | 12,432 | 12,432 | +0.0% | 51,600 | 51,534 |
| dmrg | 24 | 0 | 48 | 7,728 | 7,728 | +0.0% | 32,112 | 32,056 |
| dmrg | 24 | 0.2411 | 48 | 15,120 | 15,120 | +0.0% | 62,736 | 62,675 |
| dmrg | 24 | 0.4 | 48 | 15,120 | 15,120 | +0.0% | 62,736 | 62,670 |
| dmrg | 28 | 0 | 48 | 9,072 | 9,072 | +0.0% | 37,680 | 37,619 |
| dmrg | 28 | 0.2411 | 48 | 17,808 | 17,808 | +0.0% | 73,872 | 73,811 |
| dmrg | 28 | 0.4 | 48 | 17,808 | 17,808 | +0.0% | 73,872 | 73,806 |
| dmrg | 32 | 0 | 48 | 10,416 | 10,416 | +0.0% | 43,248 | 43,187 |
| dmrg | 32 | 0.2411 | 48 | 20,496 | 20,496 | +0.0% | 85,008 | 84,947 |
| dmrg | 32 | 0.4 | 48 | 20,496 | 20,496 | +0.0% | 85,008 | 84,942 |

*Direct compilation:* per-step × n against a full transpile of a 48-step circuit (4 pulses of 12 steps) for 27 cases with N ≤ 12: maximum difference 0 CX and 0 non-Clifford Rz. For the level-3 compilation CX is exact and the per-step non-Clifford Rz estimate is up to about 4% too high (conservative). Per-step depth × n overestimates depth by about 1–2%.

**Classical compute (single core)**

| source | N | J2 | ε | context build (DMRG+ED+trial+per-step) s | filter design s | guaranteed search s | empirical search s |
|---|---|---|---|---|---|---|---|
| dmrg | 4 | 0 | 1e-02 | 0 | 33 | 1 | 0 |
| dmrg | 4 | 0.2411 | 1e-02 | 0 | 27 | 1 | 1 |
| dmrg | 4 | 0.4 | 1e-02 | 0 | 13 | 1 | 1 |
| dmrg | 6 | 0 | 1e-01 | 0 | 33 | 1 | 1 |
| dmrg | 6 | 0 | 3e-02 | 0 | 51 | 0 | 1 |
| dmrg | 6 | 0 | 1e-02 | 0 | 72 | 0 | 1 |
| dmrg | 6 | 0 | 3e-03 | 0 | 96 | 1 | 1 |
| dmrg | 6 | 0 | 1e-03 | 0 | 117 | 1 | 2 |
| dmrg | 6 | 0.2411 | 1e-02 | 0 | 59 | 0 | 1 |
| dmrg | 6 | 0.4 | 1e-01 | 0 | 8 | 2 | 1 |
| dmrg | 6 | 0.4 | 3e-02 | 0 | 21 | 5 | 1 |
| dmrg | 6 | 0.4 | 1e-02 | 0 | 33 | 0 | 1 |
| dmrg | 6 | 0.4 | 3e-03 | 0 | 52 | 1 | 2 |
| dmrg | 6 | 0.4 | 1e-03 | 0 | 75 | 1 | 3 |
| dmrg | 8 | 0 | 1e-01 | 0 | 53 | 3 | 1 |
| dmrg | 8 | 0 | 3e-02 | 0 | 77 | 0 | 1 |
| dmrg | 8 | 0 | 1e-02 | 0 | 111 | 1 | 2 |
| dmrg | 8 | 0 | 3e-03 | 0 | 119 | 2 | 3 |
| dmrg | 8 | 0 | 1e-03 | 0 | 130 | 4 | 5 |
| dmrg | 8 | 0.2411 | 1e-02 | 1 | 91 | 2 | 3 |
| dmrg | 8 | 0.4 | 1e-02 | 1 | 54 | 1 | 2 |
| dmrg | 10 | 0 | 1e-02 | 6 | 129 | 3 | 2 |
| dmrg | 10 | 0.2411 | 1e-02 | 8 | 106 | 7 | 4 |
| dmrg | 10 | 0.4 | 1e-02 | 9 | 81 | 3 | 3 |
| dmrg | 12 | 0 | 1e-02 | 334 | 143 | 18 | 11 |
| dmrg | 12 | 0.2411 | 1e-02 | 323 | 134 | 45 | 16 |
| dmrg | 12 | 0.4 | 1e-02 | 331 | 107 | 23 | 13 |
| dmrg | 14 | 0 | 1e-02 | 6 | 151 | 177 | 39 |
| dmrg | 14 | 0.2411 | 1e-02 | 8 | 145 | 455 | 100 |
| dmrg | 14 | 0.4 | 1e-02 | 9 | 134 | 212 | 56 |
| dmrg | 16 | 0 | 1e-01 | 7 | 125 | 0 | – |
| dmrg | 16 | 0 | 1e-02 | 7 | 122 | 0 | – |
| dmrg | 16 | 0.2411 | 1e-01 | 10 | 109 | 0 | – |
| dmrg | 16 | 0.2411 | 1e-02 | 10 | 142 | 1 | – |
| dmrg | 16 | 0.4 | 1e-01 | 15 | 76 | 0 | – |
| dmrg | 16 | 0.4 | 1e-02 | 15 | 139 | 1 | – |
| dmrg | 18 | 0 | 1e-01 | 8 | 158 | 0 | – |
| dmrg | 18 | 0 | 1e-02 | 8 | 112 | 1 | – |
| dmrg | 18 | 0.2411 | 1e-01 | 17 | 136 | 0 | – |
| dmrg | 18 | 0.2411 | 1e-02 | 17 | 126 | 1 | – |
| dmrg | 18 | 0.4 | 1e-01 | 15 | 82 | 0 | – |
| dmrg | 18 | 0.4 | 1e-02 | 15 | 151 | 1 | – |
| dmrg | 20 | 0 | 1e-01 | 38 | 152 | 0 | – |
| dmrg | 20 | 0 | 1e-02 | 38 | 104 | 2 | – |
| dmrg | 20 | 0.2411 | 1e-01 | 17 | 146 | 0 | – |
| dmrg | 20 | 0.2411 | 1e-02 | 17 | 113 | 0 | – |
| dmrg | 20 | 0.4 | 1e-01 | 20 | 92 | 0 | – |
| dmrg | 20 | 0.4 | 1e-02 | 20 | 142 | 0 | – |
| dmrg | 24 | 0 | 1e-01 | 44 | 139 | 0 | – |
| dmrg | 24 | 0 | 1e-02 | 44 | 53 | 2 | – |
| dmrg | 24 | 0.2411 | 1e-01 | 33 | 164 | 0 | – |
| dmrg | 24 | 0.2411 | 1e-02 | 33 | 81 | 0 | – |
| dmrg | 24 | 0.4 | 1e-01 | 48 | 122 | 0 | – |
| dmrg | 24 | 0.4 | 1e-02 | 48 | 145 | 1 | – |
| dmrg | 28 | 0 | 1e-01 | 64 | 150 | 1 | – |
| dmrg | 28 | 0 | 1e-02 | 64 | 118 | – | – |
| dmrg | 28 | 0.2411 | 1e-01 | 75 | 157 | 0 | – |
| dmrg | 28 | 0.2411 | 1e-02 | 75 | 122 | – | – |
| dmrg | 28 | 0.4 | 1e-01 | 77 | 152 | 0 | – |
| dmrg | 28 | 0.4 | 1e-02 | 77 | 122 | 0 | – |
| dmrg | 32 | 0 | 1e-01 | 64 | 112 | 0 | – |
| dmrg | 32 | 0 | 1e-02 | 64 | 121 | – | – |
| dmrg | 32 | 0.2411 | 1e-01 | 96 | 150 | 0 | – |
| dmrg | 32 | 0.2411 | 1e-02 | 96 | 148 | – | – |
| dmrg | 32 | 0.4 | 1e-01 | 173 | 172 | 0 | – |
| dmrg | 32 | 0.4 | 1e-02 | 173 | 149 | 0 | – |

Bound violations among simulated points: **0**.

**Points without a result**

* dmrg N=28 J2=0.2411 ε=0.01 guaranteed: not achievable (leakage alone exceeds the budget / not reached in sweep) [design source: builder; the floor solver (<= 8 pulses, t_max_frac 0.5) found no feasible filter for gamma = 0.505]
* dmrg N=28 J2=0 ε=0.01 guaranteed: not achievable (leakage alone exceeds the budget / not reached in sweep) [design source: builder; the floor solver (<= 8 pulses, t_max_frac 0.5) found no feasible filter for gamma = 0.365]
* dmrg N=32 J2=0.2411 ε=0.01 guaranteed: not achievable (leakage alone exceeds the budget / not reached in sweep) [design source: builder; the floor solver (<= 8 pulses, t_max_frac 0.5) found no feasible filter for gamma = 0.442]
* dmrg N=32 J2=0 ε=0.01 guaranteed: not achievable (leakage alone exceeds the budget / not reached in sweep) [design source: builder; the floor solver (<= 8 pulses, t_max_frac 0.5) found no feasible filter for gamma = 0.308]

## Figures

![Guaranteed designs vs N](figures/fig_scaling_guaranteed.png)

![Guaranteed vs empirical](figures/fig_guaranteed_vs_empirical.png)

![ε scan](figures/fig_eps_scan.png)

## Reproducing

```
export OMP_NUM_THREADS=1
for i in 0 1 2 3; do python scripts/run_results.py --worker $i --of 4 & done; wait   # resumable
python scripts/add_direct_counts.py   # direct-compilation counts (run_results.py now writes them itself)
python scripts/audit_results.py [--alpha]   # independent re-derivation of every stored bound and count; must print ALL CHECKS PASSED
python scripts/build_results.py
```

Environment: Python 3.12, qiskit, qiskit-aer, quimb, mps-to-circuit, scipy, numpy, matplotlib (see `results/data/*.json` for the exact settings of each point: order, bandwidth, trial layers, seed).