# Certified cosine-filter state preparation for the J1–J2 chain: certificates, gate counts, and where the certificate is loose

*Scope: only results that are produced by a script in this repository from stored data, each named below. Nothing here is a hardware result, and nothing is modelled or extrapolated. Anything not shown is listed in §8 as not done.*

## 1. Summary

We evaluate the Stetcu–Baroni–Carlson (SBC) cosine-filter state preparation for the spin-½ J1–J2 chain, H = Σ J1 **S**ᵢ·**S**ᵢ₊₁ + J2 **S**ᵢ·**S**ᵢ₊₂ (J1 = 1). A classical front end (DMRG or exact diagonalisation) supplies the squared trial overlap γ, the gap Δ and a certified bandwidth W. The filter F(E) = ∏ᵢ cos(E tᵢ + φᵢ) is Trotterised, and the prepared state carries a rigorous trace-distance certificate.

1. **Certificate and cost** (open chains, N = 4–32, J2 ∈ {0, 0.2411, 0.4}, ε = 10⁻²). With the first-order Trotter constant α computed exactly from Pauli commutators, the certified circuit for N ≤ 14 uses 76–24,027 Trotter steps and 1,539–3,604,089 RZZ gates. Certified RZZ scales as ≈ N⁵ and steps as ≈ N⁴ (`results/RESULTS.md`).
2. **The certificate is loose.** The smallest step count whose *measured* infidelity meets ε (the "oracle") is 15–586× smaller than the certified count (median 190×), N ≤ 14.
3. **Second-order Trotter formula.** Replacing the first-order bound by the second-order (Strang) commutator bound of Childs et al. (their Proposition 16) lowers the certified step count 20–59× (N = 6–12, 12 cases). Built and transpiled in the RZZ basis, the certified second-order circuit uses 3.0–8.6% of the first-order RZZ. The bound held in every simulation.
4. **Certified versus oracle at each order.** Certified/oracle RZZ is 66–438× at first order and 5–19× at second order (12 cases). At the oracle level the two orders cost comparable amounts (second-order oracle = 0.59–1.51× the first-order oracle).
5. **Why the certificate is loose.** At the certified step counts the first-order Trotter error is level mixing, not phase dephasing, and first-order perturbation theory reproduces the measured error to within 0.4% (9 cases, N ≤ 10).
6. **Spin symmetry.** Adding an S² projection stage raises the certified cost in all 29 feasible cases; a symmetric trial state without the stage lowers it in all 29 (§6).

## 2. Setup and certificate

*Circuit.* One pulse is W(t, φ) = H_anc · e^{−it H_s⊗Z_anc} · R_z(2φ) · H_anc on system plus one ancilla, with H_s = (H − E₀)/W (spectrum in [0, 1]). On an eigenstate of energy E it gives cos(tE + φ)|0⟩ − i sin(tE + φ)|1⟩ on the ancilla, so post-selecting the ancilla on 0 applies cos(tH_s + φ). The ancilla is measured and reset after every pulse and the run is aborted on outcome 1 (`core/circuits.py`). A product filter needs the reset: stacking pulses without it merges them into a single cosine.

*Filter quantities.* F(E) = ∏ᵢ cos(E tᵢ + φᵢ); η = sup over E ∈ [Δ, 1] of |F(E)|/|F(0)|; p_g = γ F(0)²; leakage ℓ = (1−γ)η²/(γ + (1−γ)η²) (Δ is the scaled gap, γ = |⟨E₀|ψ_trial⟩|²).

*Certificate.* The infidelity of the prepared state obeys 1 − 𝓕 ≤ (√ℓ + ε_T/√p_g)², where ε_T is a certified bound on the Trotter error of the filter. First order: ε_T = Σᵢ α tᵢ²/(2kᵢ) = αT²/(2n) for equal steps (kᵢ steps in pulse i, T = Σtᵢ, n = Σkᵢ), with α = Σ_g ‖[H_g, Σ_{g'>g} H_g']‖ evaluated exactly on local Pauli commutators (`core/trotter.py`). With target ε the step count is n = ⌈αT²/(2√p_g(√ε − √ℓ))⌉; the tabulated value is the smallest n for which the bound on the snapped pulse grid is ≤ ε (it agrees with the formula to within a few steps). Proofs: `CERTIFICATION.tex`, Appendix A of `paper.tex`.

*Design and grid.* The design (`floor.py`) fixes T and the pulses; the step size is dt = T/n and pulse i gets kᵢ = round(tᵢ/dt) steps (largest-remainder rule), after which the design is re-certified on the snapped grid and the phases re-optimised.

*Inputs.* γ and the gap are from ED for N ≤ 14. For N > 14 γ is a DMRG-based overlap and is not certified, so those rows are conditional on it (`results/RESULTS.md` §1).

*Gate counting.* RZZ(θ) = exp(−iθ ZZ/2); basis {rzz, rz, sx, x}; all-to-all; filter only (no trial circuit); success path. Level-3 compilation uses exactly 2 RZZ per Pauli term per step and direct compilation 4. Per-step × n equals the full transpile exactly for RZZ in a 48-step validation transpile of each of the 66 stored points (`scripts/audit_results.py`), and also in all 7 first-order full transpiles of the certified circuits run so far (level 3 and direct) compared with the formula 2 (4) × terms × n; non-Clifford counts from the per-step method are accurate to about ±20%.

## 3. Resource results (first order)

From `results/RESULTS.md` (66 (N, J2, ε) points; 0 bound violations among simulated points; `python scripts/audit_results.py` re-derives every stored quantity).

* N ≤ 14, guaranteed design, ε = 10⁻²: 76–24,027 steps, 1,539–3,604,089 RZZ, total filter time 1.9–13.3 in units of 1/J1.
* Guaranteed RZZ scales as N^5.1 (J2 = 0), N^5.2 (J2 = 0.2411), N^4.9 (J2 = 0.4); steps as N^4.0, N^4.0, N^3.7 (fits over N = 6–24/28/32).
* Guaranteed versus oracle (N ≤ 14): 15–586× in steps and 14–582× in RZZ, median 190×. The oracle needs the exact ground state, so it exists only where simulation is possible.
* N = 16–32 are counted, not simulated.

## 4. Second-order Trotter formula

### 4.1 The bound

For H = Σ_{a=1}^Γ H_a and S₂(t) = ∏_{a=Γ→1} e^{−iH_a t/2} ∏_{a=1→Γ} e^{−iH_a t/2}, with B_a = Σ_{b>a} H_b:

‖S₂(t) − e^{−itH}‖ ≤ t³B₂, B₂ = (1/12) Σ_a ‖[B_a,[B_a,H_a]]‖ + (1/24) Σ_a ‖[H_a,[H_a,B_a]]‖.

This is Proposition 16, Eq. (152) of Childs, Su, Tran, Wiebe, Zhu (arXiv:1912.08854v3), obtained there (Eqs. 146–151) by telescoping over a and the two-term bound of their Eq. (150); the first-order α above is their Proposition 15, Eq. (145). With k steps of length t/k the error is ≤ t³B₂/k². B₂ is evaluated exactly on local Pauli commutators, maximised over the forward and reversed term order (`trotter2/alpha2.py`); the constant for H⊗Z_anc equals that of H since Z² = 1.

Checks of my implementation (not of the formula): the two-term bound over 1,800 random Hermitian pairs, t from 0.05 to 30: maximum error/bound 0.998, never violated (`trotter2/validate_two_term.py`); the multi-term bound on J1–J2 chains with N = 4–8, open and periodic, t ∈ {3, 10}, k = 1–16: error ≤ 0.28 × bound in every case (`trotter2/validate_alpha2.py`).

### 4.2 Built circuits, step counts and RZZ

`trotter2/second_order_circuit.py` takes each stored guaranteed design (pulse times, phases, ε = 10⁻²), keeps ℓ and p_g, and chooses integer kᵢ = ⌈s·tᵢ⌉ with the smallest s whose bound (√ℓ + Σ B₂tᵢ³/kᵢ² /√p_g)² is ≤ ε. It builds the early-abort circuit with qiskit's second-order `SuzukiTrotter` (term order preserved), transpiles it fully to the RZZ basis (level 3 and direct), simulates the post-selected filter against ED with the stored design's trial state, and finds the oracle by bisection over the scaling s. The first-order oracle is the one in §3. Counts are totals for the whole filter circuit. First-order RZZ = 2 × (Pauli terms) × n, which equals the full transpile wherever both were run; second-order RZZ is a full transpile.

| N | J2 | n certified 1st | n certified 2nd | RZZ certified 1st | RZZ certified 2nd | ratio | RZZ oracle 1st (n) | RZZ oracle 2nd (n) | certified/oracle 1st | certified/oracle 2nd |
|--|--|--|--|--|--|--|--|--|--|--|
| 6 | 0 | 590 | 28 | 17,700 | 1,304 | 0.074 | 150 (5) | 224 (4) | 118× | 6× |
| 6 | 0.2411 | 578 | 26 | 31,212 | 2,358 | 0.076 | 324 (6) | 489 (5) | 96× | 5× |
| 6 | 0.4 | 394 | 20 | 21,276 | 1,824 | 0.086 | 324 (6) | 400 (4) | 66× | 5× |
| 8 | 0 | 1,817 | 56 | 76,314 | 3,730 | 0.049 | 462 (11) | 545 (7) | 165× | 7× |
| 8 | 0.2411 | 2,398 | 73 | 187,044 | 9,799 | 0.052 | 858 (11) | 1,021 (7) | 218× | 10× |
| 8 | 0.4 | 1,059 | 36 | 82,602 | 4,848 | 0.059 | 858 (11) | 725 (5) | 96× | 7× |
| 10 | 0 | 3,917 | 89 | 211,518 | 7,679 | 0.036 | 540 (10) | 709 (7) | 392× | 11× |
| 10 | 0.2411 | 5,777 | 125 | 589,254 | 22,239 | 0.038 | 1,734 (17) | 1,707 (9) | 340× | 13× |
| 10 | 0.4 | 2,890 | 66 | 294,780 | 11,758 | 0.040 | 1,632 (16) | 961 (5) | 181× | 12× |
| 12 | 0 | 8,797 | 168 | 580,602 | 17,824 | 0.031 | 1,914 (29) | 1,339 (11) | 303× | 13× |
| 12 | 0.2411 | 10,960 | 187 | 1,380,960 | 41,465 | 0.030 | 3,150 (25) | 2,127 (9) | 438× | 19× |
| 12 | 0.4 | 5,318 | 107 | 670,068 | 23,785 | 0.035 | 3,276 (26) | 2,348 (10) | 205× | 10× |

(Same table: `trotter2/circuit_table.md`, `results/RESULTS.md` §10.)

* Certified step counts: n₁/n₂ = 20–59, median 38, increasing with N for each J2.
* Certified RZZ: second order uses 3.0–8.6% of first order at level 3 (median 4.4%) and 3.4–10.2% with direct compilation (median 5.3%).
* A second-order step costs more than a first-order step: 1.55–1.73× the RZZ in the cases with both circuits transpiled in full (`results/RESULTS.md` §10, currently 7 cases).
* Certified/oracle is 66–438× at first order and 5–19× at second order. It is not monotone in N for every J2 (e.g. J2 = 0 first order: 118, 165, 392, 303 for N = 6, 8, 10, 12).
* Second-order oracle / first-order oracle RZZ = 0.59–1.51 (above 1 in 6 of the 12 cases). The oracle is a practical value from a bisection (the measured infidelity is not strictly monotone in n), not a guaranteed minimum, and it is 4–29 steps, so it is coarsely quantised.
* In all 12 simulations the measured infidelity (3×10⁻⁴–1.6×10⁻³, against ε = 10⁻²) and the state error were below their certified bounds; the measured state error is at most 1.0% of its bound. The qiskit circuit and an independent simulation agree to ≤ 9×10⁻¹⁶ in 1 − |⟨·|·⟩|² for N ≤ 8 (for N ≥ 10 only the simulation was run).

## 5. Why the certificate is loose

Detailed in `trotter2/DEPHASING.md`; ED, N = 6–10 (N = 12 for the level counts).

* **Level weights** (`trotter2/dephasing.py`). In the eigenbasis of H_s the trial state has weight on 7–79 levels above 10⁻⁴ (N = 6–12, 10 cases; a level is a distinct (energy, total spin) pair), reduced to 1–5 after the filter; 99.99% of the weight needs 9–291 levels before and 3–9 after (`trotter2/dephasing_table.md`).
* **No O(τ) phase error.** For a Lie–Trotter step e^{−iτH_eff}, H_eff = H + τX₁ + …, X₁ = −(i/2) Σ_{a>b}[H_a, H_b]. All Pauli terms XX, YY, ZZ are real, so X₁ is imaginary Hermitian with zero diagonal in a real eigenbasis. Checked: max |diag X₁| ≤ 3.8×10⁻¹⁸ over the 9 cases. Measured quasi-energy shifts of one step scale as τ²: the ground-state shift has fitted exponent 2.00 for both orders in all 9 cases, and the largest shift over the whole spectrum has exponent 1.76–2.00.
* **The O(τ) error is level mixing.** First-order perturbation theory in τX₁ gives the amplitude generated from level n to level m in a pulse of length t as δV_mn = −τ X_mn (e^{−iEₙt} − e^{−iE_m t})/(E_m − Eₙ), bounded by 2τ|X_mn|/|E_m − Eₙ| independent of t for non-resonant levels (the commutator bound αt²/(2k) = αtτ/2 grows with t). Propagated through all pulses it reproduces the measured first-order trace distance to within 0.4% (ratio 0.996–1.003, 9 cases; `trotter2/mixing_theory.py`, `trotter2/mixing_table.md`).
* **Which states.** 99% of the predicted first-order error comes from 2–27 source eigenstates (of 64–1024), and 84–98% of it is generated in the last two (longest) pulses.
* **Scaling** (N = 6, J2 = 0.4 only, `trotter2/scaling_check.py`). First-order trace distance ∝ 1/n (fitted exponent −1.00 over n = 98–1576). Second order falls by a factor of 4 per doubling of n at n = 20, 40, 80 (5.7×10⁻⁴, 1.4×10⁻⁴, 3.6×10⁻⁵), as expected from the time-symmetry of the Strang step; the fit over all five points is −2.3 because of the first two points.
* At the certified first-order step counts the measured Trotter trace distance is 2×10⁻⁴–10⁻³ against a certified budget of about 0.056.

## 6. Spin symmetry

`s2proj/` (exact classical inputs for open and periodic chains, N = 8–16, J2 ∈ {0, 0.2411, 0.4}, ε = 10⁻²; costs from the first-order certificate with designs from the pipeline's floor solver and RZZ per step from the model 2 × (Pauli terms) for the energy stage and a transpiled count for the S² stage; filter RZZ per run, no trial-circuit cost).

* **S² projection stage.** S² = const + Σ_{i<j} P_ij and the Jucys–Murphy blocks X_k commute, so exp(−itS²) = ∏_k exp(−itX_k) exactly. Zero-phase pulses at t = π/4, π/8 remove all J > 0 components (checked on the open chain N = 8, J2 = 0.2411 with exact evolution: J > 0 weight ≤ 8×10⁻¹⁸ after the S² stage, `s2proj/validate_s2.py`). A: filter without projection (gap = min of triplet and singlet gaps); B: with the S² stage. In all 29 feasible cases B costs more than A: B/A = 1.30–10.1, median 3.8 (`python s2proj/analyze.py`). One case (periodic N = 16, J2 = 0.4) is infeasible.
* **Trotterised circuit preserves total spin.** Every bond's XX, YY, ZZ factors commute with each other and with the ancilla Z. Checked numerically: starting from the ground state, two pulses of the first- and second-order circuits leave weight of magnitude ≤ 7×10⁻¹⁷ outside J = 0 in both ancilla branches (N = 6, 8; `s2proj/check_trotter_spin.py`). So an S = 0 trial state stays S = 0 and the filter only has to resolve the singlet-sector gap.
* **Symmetric trial.** A trial of singlet pairs plus two layers of SU(2)-invariant gates, with the filter resolving the singlet gap and no S² stage, costs 0.05–0.99× A (median 0.24×, 29 cases; `s2proj/analyze.py`). In the 17 cases with an attribution run (`s2proj/sym_attrib.py`) the gain A/SYM has median 4.2 (1.2–18.6), which splits into a better overlap, median 1.6 (1.2–2.8), times the larger singlet gap, median 2.7 (1.0–8.0).

## 7. Level structure and filter-prepared periodic states

* **Open chain** (`python scripts/summarize_low_levels.py`, ED, J2 grid of 74 points in [0, 1]): for N = 8, 12, 16, 20 the lowest excitation is a triplet at every grid point. The smallest value of E(first excited singlet) − E(lowest triplet) over the grid is 0.0112 (J2 = 0.925), 0.0059 (0.625), 0.0020 (0.595), 0.0009 (0.575) respectively, so the closest approach shrinks with N; at J2 = 0.2411 the difference is 0.398, 0.277, 0.212, 0.172. A crossing between grid points cannot be excluded, and nothing is claimed for larger N.
* **Periodic chain** (same script): the first excited singlet and the lowest triplet cross at J2 = 0.2463, 0.2435, 0.2425, 0.2420 for N = 8, 12, 16, 20 (linear interpolation on a 0.025 grid), approaching the reference value J2c ≈ 0.2411 (S. Eggert, Phys. Rev. B 54, R9612) from above.
* **Filter-prepared periodic states** (`python pbc_corr/summarize_filter_states.py`; ideal exact filter with the certified design, *not* the Trotterised circuit): in 19 cases (N = 8–20, J2 ∈ {0, 0.2411, 0.4, 0.6 at N = 8}, ε = 10⁻²) the energy is within 1.9×10⁻⁴ relative error of the exact ground state, the largest relative error over the stored correlators is 4.4% and over the block entropies 5.7%. Nine cases have no feasible design (J2 = 0.6 for N ≥ 10; J2 = 0.4 for N ≥ 16).

## 8. What is not done or not shown

* No hardware run, no noise model, no noisy simulation.
* The first-order full transpiles exist for 7 of the 12 second-order cases (per-step × n is exact for RZZ, so the RZZ totals do not depend on them; depth and non-Clifford counts of the first-order circuits do). The N = 12 dephasing analysis is missing for J2 = 0.2411 and 0.4.
* Second-order circuits and the mixing analysis are computed only for N ≤ 12 and ≤ 10 respectively; the first-order perturbation theory has no second-order counterpart here, and its second-order scaling is measured only at N = 6, J2 = 0.4.
* Pulse times and phases for the second-order circuits are the first-order designs, not re-optimised for the second-order budget.
* The trial circuit is excluded from all gate counts; for N > 14 the trial overlap is a DMRG-based number, not certified.
* Level-3 non-Clifford counts are ±20% estimates.
* The oracle is only available where the exact ground state is (N ≤ 14).
* The spin-symmetry costs (§6) use ED inputs and the cost model stated there (energy-stage RZZ per step 2 × terms; first-order certificate), not full transpiles of the combined circuits.
* `paper.tex` is an earlier draft in the CX basis and does not contain §4–§7; `results/RESULTS.tex` does not yet contain §10 of `results/RESULTS.md`.

## 9. Where things are

| Item | Location |
|--|--|
| Pipeline, certificate, designs | `pipeline.py`, `floor.py`, `core/` |
| First-order results (RZZ basis) and audit | `results/RESULTS.md`, `results/data/`, `scripts/audit_results.py` |
| Proofs | `CERTIFICATION.tex`, `paper.tex` App. A–B |
| Second-order bound, circuits, oracle | `trotter2/` (`alpha2.py`, `second_order_circuit.py`, `first_order_full.py`, `validate_*.py`, `circuit_table.md`) |
| Mixing / dephasing analysis | `trotter2/DEPHASING.md` (`dephasing.py`, `mixing_theory.py`, `scaling_check.py`) |
| Spin symmetry | `s2proj/` |
| Level structure, periodic filter states | `scripts/summarize_low_levels.py`, `pbc_corr/summarize_filter_states.py` |
| Local reproduction | `LOCAL_VERIFY.md`, `requirements.txt` |

## References

* A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, Theory of Trotter error with commutator scaling, Phys. Rev. X 11, 011020 (2021), arXiv:1912.08854 (Propositions 15 and 16 read in the arXiv v3 text).
* The SBC filter paper supplied with the project (`a7585736-FIlter.pdf`).
* S. Eggert, Phys. Rev. B 54, R9612 (1996), for the reference value J2c ≈ 0.2411.
