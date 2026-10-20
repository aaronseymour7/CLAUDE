# A certified cosine-filter state-preparation workflow for the J1–J2 Heisenberg chain: what it costs, what tightens it, and what does not help

*Status: technical note summarising the work in this directory. Every number is produced by a script in this repository (named below). Claims that rest on modelled inputs or on literature I could only read at abstract level are marked **[modelled]** or **[abstract-level]**.*

## 1. Summary

We study the Stetcu–Baroni–Carlson (SBC) cosine-filter state preparation for the spin-1/2 J1–J2 chain, H = Σ J1 S_i·S_{i+1} + J2 S_i·S_{i+2}. A classical front end (DMRG/ED) supplies the trial overlap γ, the gap Δ and a certified bandwidth W. From these we design a filter F(E) = ∏ cos(E t_i + φ_i), and Trotterise it with a rigorous trace-distance certificate on the prepared state.

Findings:

1. **Certificate and cost (N = 4–32, J2 ∈ {0, 0.2411, 0.4}, ε = 10⁻²).** With a first-order Trotter error constant α computed *exactly* from the Pauli commutators, the certified circuit needs 76–24,027 steps for N ≤ 14. Counted in an RZZ gate basis, certified RZZ count scales as ~N^5 and steps as ~N^4 (`results/RESULTS.md`, §Summary).
2. **The certificate is loose.** The smallest step count that *empirically* meets ε (oracle, measured against the exact ground state) is 15–586× smaller than the certified count (median 190×).
3. **A second-order (Strang) certificate closes most of that gap, and built circuits confirm it.** Replacing the first-order Trotter error αT²/(2n) by the commutator-scaling Strang bound T³B₂/n² reduces the certified step count by 20–59× at ε = 10⁻² for N = 6–12 (median 38×), growing with N over that range. I then built the second-order circuits, transpiled them to the RZZ basis and simulated them against ED for all 12 cases (N = 6–12, J2 ∈ {0, 0.2411, 0.4}): certified filter RZZ falls to 3.0–8.6% of the first-order count (median 4.4%, level-3 compilation) and 3.4–10.2% with direct compilation. The bound held in every simulation, and the qiskit circuit matches the simulation to 10⁻¹⁵ (N ≤ 8). Against a second-order oracle (smallest n whose measured 1−F ≤ ε) the certified count is 5–21× too large, versus 66–438× at first order. (`trotter2/`, §4.3.)
4. **Negative results.** (a) Adding an exact total-spin (S²) projection stage to the filter is not cost-effective in 1D. (b) An SU(2)-symmetric trial state gains about ×1.4 from the ansatz and a median ×2.7 from the enlarged gap, but this is very likely known in the literature. (c) In 1D, MPS-to-circuit compilation is the stronger competitor.
5. **Benchmark and outlook.** Filter-prepared periodic-chain states reproduce exact SzSz correlators, energies and entanglement entropy to within 0.1–3% at sizes the filter can reach. A 2D crossover estimate against MPS circuits is **[modelled]** and is not a result.

## 2. Setup and certificate

*Filter.* An ancilla block applies cos(H_s t_i + φ_i) with H_s the scaled Hamiltonian. F(E) = ∏_i cos(E t_i + φ_i). Let η = sup_{[Δ,1]} |F| / |F(0)|, p_g = γ F(0)² the ideal success probability, and leakage ℓ = (1−γ)η² / (γ + (1−γ)η²).

*Certificate.* The prepared state's infidelity obeys 1 − F ≤ (√ℓ + ε_T/√p_g)², where ε_T is the operator-norm error of the Trotterised filter. For first-order Trotter, ε_T = Σ α t_i²/(2k_i) = αT²/(2n) and α = Σ_g ‖[H_g, Σ_{g'>g} H_g']‖, computed exactly from local Pauli commutators (`core/trotter.py`, `alpha_commutator`). The step count is n = ⌈αT²/(2√p_g(√ε − √ℓ))⌉. The bandwidth is the Weyl bound W = Σ_bonds λ_max − E0, which equals E_top − E0 for these chains. Proofs are in `CERTIFICATION.tex` and Appendix A of `paper.tex`.

*Inputs.* γ and Δ come from DMRG (and ED for N ≤ 14). For N > 14 the overlap γ is a DMRG-MPS overlap and is **uncertified**; the certificate is then conditional on that input (`results/RESULTS.md`, §1).

*Gate counting.* Two-qubit gate is RZZ(θ) = exp(−iθ ZZ/2); the basis is {rzz, rz, sx, x}. Level-3 compilation uses exactly 2 RZZ per Pauli term per Trotter step; direct compilation uses 4. Per-step counts times n reproduce full transpiled counts exactly for RZZ; the non-Clifford estimate is accurate to about ±20% (`scripts/audit_results.py`).

## 3. Resource results

All from `results/RESULTS.md` (66 (N, J2, ε) points, N = 4–32, 0 bound violations among simulated points).

* N ≤ 14, guaranteed design, ε = 10⁻²: 76–24,027 steps, 1,539–3,604,089 RZZ, T(phys) 1.9–13.3 in 1/J1.
* Scaling of guaranteed RZZ with N: ~N^5.1 (J2 = 0), N^5.2 (J2 = 0.2411), N^4.9 (J2 = 0.4); steps ~N^4.0, N^4.0, N^3.7; RZZ per step ~N^1.1. Fits over N = 6–24/28/32.
* Guaranteed versus empirical: 15–586× steps, 14–582× RZZ, median 190×.
* N = 16–32 are counted, not simulated; no ED cross-check exists there.

## 4. Tightening the Trotter certificate

### 4.1 The bound

For H = Σ_{a=1}^Γ H_a, let S₂(t) = ∏_{a=Γ→1} e^{−iH_a t/2} ∏_{a=1→Γ} e^{−iH_a t/2} and B_a = Σ_{b>a} H_b. Then

‖S₂(t) − e^{−itH}‖ ≤ t³ B₂,  B₂ = (1/12) Σ_a ‖[B_a,[B_a,H_a]]‖ + (1/24) Σ_a ‖[H_a,[H_a,B_a]]‖,

so k steps of t/k give error ≤ t³B₂/k². This is the second-order case of the commutator-scaling theory of Childs–Su–Tran–Wiebe–Zhu, PRX 11, 011020 (2021) [arXiv:1912.08854].

*Verification status.* Checked against the paper (arXiv:1912.08854v3, supplied by the user): **Proposition 16, Eq. (152)** states exactly this bound, ‖S₂(t) − e^{−itH}‖ ≤ (t³/12) Σ_{γ₁}‖[Σ_{γ₃>γ₁}H_{γ₃}, [Σ_{γ₂>γ₁}H_{γ₂}, H_{γ₁}]]‖ + (t³/24) Σ_{γ₁}‖[H_{γ₁}, [H_{γ₁}, Σ_{γ₂>γ₁}H_{γ₂}]]‖, for S₂(t) = ∏_{γ=Γ}^{1} e^{−i(t/2)H_γ} ∏_{γ=1}^{Γ} e^{−i(t/2)H_γ}, derived (their Eqs. 146–151) by the same telescoping over γ₁ and the two-term bound (their Eq. 150). **Proposition 15, Eq. (145)**, ‖S₁(t) − e^{−itH}‖ ≤ (t²/2) Σ_{γ₁}‖[Σ_{γ₂>γ₁}H_{γ₂}, H_{γ₁}]‖, is the first-order α used throughout this note. The paper states the bound for the term ordering in its definition; the code takes the maximum over the forward and reversed order, which is a superset. The numerical checks below are therefore a check of my implementation, not of the formula:

* **Telescoping argument (as in the paper).** Write S₂ for terms a..Γ as e^{−iH_a t/2} S₂^{(a+1..Γ)} e^{−iH_a t/2}. Unitarity of the outer factors gives ‖S₂ − e^{−itH}‖ ≤ ‖S₂^{(a+1..Γ)} − e^{−itB_a}‖ + ‖e^{−iH_a t/2}e^{−itB_a}e^{−iH_a t/2} − e^{−it(H_a+B_a)}‖. Iterating, the multi-term bound is the sum over a of the **two-term** Strang bound with A = H_a, B = B_a. The two-term bound ‖e^{−itA/2}e^{−itB}e^{−itA/2} − e^{−it(A+B)}‖ ≤ t³(‖[B,[B,A]]‖/12 + ‖[A,[A,B]]‖/24) is the paper's Eq. (150).
* **Numerical check of the two-term bound** (`trotter2/validate_two_term.py`): 1,800 random Hermitian pairs, t from 0.05 to 30. Maximum error/bound = **0.998**. The bound was never violated, and a ratio this close to 1 means the constants cannot be made smaller (smaller constants would be violated). This is evidence, not a proof.
* **Numerical check of the multi-term bound on the actual Hamiltonians** (`trotter2/validate_alpha2.py`): J1–J2 chains with N = 4–8, open and periodic, t ∈ {3, 10}, k = 1–16; measured error ≤ 0.28 × bound in every case. Bound code: `trotter2/alpha2.py` (exact local evaluation of the nested commutators; maximum over forward and reversed term order).

The paper also notes that these two-term-based bounds are tight up to the triangle inequality, which is consistent with the 0.998 ratio above.

### 4.2 Effect on the certified step count

Only ε_T changes in the certificate; ℓ, p_g and the trace-distance form are unchanged. With the allowed Trotter budget d_T√p_lb taken from the stored guaranteed design, n₁ = αT²/(2·allowed) and n₂ = √(B₂T³/allowed). The recomputed n₁ equals the stored certified n to within +1 in every case (a consistency check of the pipeline). Results at ε = 10⁻² (`trotter2/table.md`):

| N | J2 | n₁ (first order) | n₂ (second order) | n₁/n₂ |
|--|--|--|--|--|
| 6 | 0.2411 | 579 | 26 | 22 |
| 8 | 0.2411 | 2,399 | 73 | 33 |
| 10 | 0.2411 | 5,777 | 125 | 46 |
| 12 | 0.2411 | 10,960 | 186 | 59 |
| 12 | 0.4 | 5,319 | 107 | 50 |

Over all twelve (N, J2) cases n₁/n₂ = 20–59, median 38. (§4.3 replaces the two-sweeps accounting estimate by built and transpiled circuits.)

### 4.3 Built second-order circuits: RZZ counts and simulation

`trotter2/second_order_circuit.py` takes each stored guaranteed design (pulse times t_i, phases, ε = 10⁻²), keeps the leakage and p_g terms of the certificate, and chooses integer per-pulse step counts k_i = ⌈s t_i⌉ with the smallest s whose total bound (√ℓ + Σ B₂ t_i³/k_i² / √p_g)² ≤ ε. The B₂ of H⊗Z_anc equals that of H since Z² = 1. It then (i) builds the early-abort filter circuit with qiskit's second-order `SuzukiTrotter` (term order preserved), (ii) transpiles the full circuit to {rzz, rz, sx, x} at optimisation level 3 and level 1 ("direct"), (iii) simulates the post-selected second-order filter against ED using the stored design's trial state, and (iv) bisects for the smallest k-scaling whose measured 1 − F ≤ ε (second-order oracle). Counts are for the filter only (no trial circuit), as in the first-order tables. First-order RZZ is the stored per-step count × n (2 per Pauli term at level 3, 4 direct), which the pipeline had verified against full transpiles.

Result (`trotter2/circuit_table.md`, `trotter2/circuit_N*_J2*.json`):

| N | J2 | n₁ cert | n₂ cert | RZZ first (L3) | RZZ second (L3) | ratio | RZZ first (direct) | RZZ second (direct) | ratio | n₂ oracle | cert/oracle: 1st → 2nd order |
|--|--|--|--|--|--|--|--|--|--|--|--|
| 6 | 0 | 590 | 28 | 17,700 | 1,304 | 0.074 | 35,400 | 3,360 | 0.095 | 4 | 118× → 7× |
| 6 | 0.2411 | 578 | 26 | 31,212 | 2,358 | 0.076 | 62,424 | 5,616 | 0.090 | 5 | 96× → 5× |
| 6 | 0.4 | 394 | 20 | 21,276 | 1,824 | 0.086 | 42,552 | 4,320 | 0.102 | 4 | 66× → 5× |
| 8 | 0 | 1,817 | 56 | 76,314 | 3,730 | 0.049 | 152,628 | 9,408 | 0.062 | 7 | 165× → 8× |
| 8 | 0.2411 | 2,398 | 73 | 187,044 | 9,799 | 0.052 | 374,088 | 22,776 | 0.061 | 7 | 218× → 10× |
| 8 | 0.4 | 1,059 | 36 | 82,602 | 4,848 | 0.059 | 165,204 | 11,232 | 0.068 | 5 | 96× → 7× |
| 10 | 0 | 3,917 | 89 | 211,518 | 7,679 | 0.036 | 423,036 | 19,224 | 0.045 | 7 | 392× → 13× |
| 10 | 0.2411 | 5,777 | 125 | 589,254 | 22,239 | 0.038 | 1,178,508 | 51,000 | 0.043 | 9 | 340× → 14× |
| 10 | 0.4 | 2,890 | 66 | 294,780 | 11,758 | 0.040 | 589,560 | 26,928 | 0.046 | 5 | 181× → 13× |
| 12 | 0 | 8,797 | 168 | 580,602 | 17,824 | 0.031 | 1,161,204 | 44,352 | 0.038 | 11 | 303× → 15× |
| 12 | 0.2411 | 10,960 | 187 | 1,380,960 | 41,465 | 0.030 | 2,761,920 | 94,248 | 0.034 | 9 | 438× → 21× |
| 12 | 0.4 | 5,318 | 107 | 670,068 | 23,785 | 0.035 | 1,340,136 | 53,928 | 0.040 | 10 | 205× → 11× |

* **RZZ.** The built second-order circuits use 3.0–8.6% of the first-order certified RZZ at level 3 (median 4.4%) and 3.4–10.2% with direct compilation (median 5.3%). The transpiler merges the adjacent half-steps of consecutive steps, so second-order RZZ per step is 0.78–0.88 of the naive two sweeps; my earlier "two sweeps" accounting was mildly pessimistic.
* **Bound and circuit.** In all 12 simulations the measured 1 − F (3×10⁻⁴ to 1.6×10⁻³, against ε = 10⁻²) and the state error were below their certified bounds. The measured state error is at most 1.0% of its bound, so the Trotter part of the second-order certificate is itself still about two orders of magnitude loose. The qiskit circuit and my independent simulation agree to ≤ 9×10⁻¹⁶ in 1 − |⟨·|·⟩|² (checked for N ≤ 8; at N ≥ 10 only the simulation was run).
* **Remaining gap.** At equal order the certified n is 5–21× the oracle at second order versus 66–438× at first order; the ratio grows with N in both. In absolute RZZ at N = 12, J2 = 0.4: certified second order 23,785; second-order oracle 2,348; first-order oracle 3,276 (= 2 × 63 terms × 26 steps); certified first order 670,068.

*Limits.* (i) N ≤ 12 only: the commutator evaluation is slow in Python. (ii) The pulse times and phases are the stored first-order designs; re-optimising the pulse design for a second-order Trotter budget was not attempted. (iii) The trial circuit is not included in these counts. (iv) The bound's constants are verified numerically only (§4.1).

## 5. What does not help

**S² projection (1D).** Total spin satisfies S² = const + Σ_{i<j} P_ij, and the Jucys–Murphy blocks X_k commute, so exp(−itS²) = ∏_k exp(−itX_k) exactly. A block-aware commutator constant (exact for N ≤ 12: 114, 231, 408 at N = 8, 10, 12; ~N³) and zero-phase pulses at t = π/4, π/8, (π/16) remove all J > 0 components. The stage costs more than the gap it opens in 1D, so the unprojected filter was better (`s2proj/`, `results/` S² files). Because each bond's XX, YY, ZZ factors commute with each other and with the ancilla Z, the filter circuit itself (including its Trotter error) preserves total spin: an exactly S = 0 trial stays S = 0, and the filter then resolves only the singlet-sector gap.

**Symmetric trial.** An SU(2)-symmetric trial (singlet pairs plus swap-power gates) gains about ×1.4 from the ansatz and a median ×2.7 from the larger singlet gap (`s2proj/sym_*`, `attrib_*`). The prior-art search found related work **[abstract-level]**: the QPF paper itself notes that projection raises the gap (arXiv:2211.10545); total-spin filtering (arXiv:2106.10867); spin-adapted circuits (arXiv:2412.14797); a crossover study (arXiv:2210.14109). I treat this component as likely known and make no novelty claim.

**MPS-to-circuit.** In 1D the MPS-circuit competitor is cheaper than the filter at the sizes studied (`s2proj/results/mps_circuit_vs_layers.json`).

**Low-energy Trotter bounds from the literature.** Şahinoğlu–Somma (arXiv:2006.02660), arXiv:2402.10362 and Mizuta–Kuwahara (PRL 2025) **[abstract-level]** were examined as a route to a low-energy-restricted certificate. At N ≤ 32 their constants were vacuous in my earlier analysis; that conclusion is from the earlier session and is not re-run here.

## 6. Physics benchmark: first excitation and emergent field theory

* **First excitation is a triplet** for open chains N ≤ 16, J2 ∈ [0, 1] (`scripts/check_first_excitation.py`, `results/low_levels*.json`).
* **Periodic chains** show a singlet–triplet level crossing whose location converges toward J2 ≈ 0.2411 with N; the open chain has no such crossing, as expected since it carries only the vacuum tower (triplet at level 1, singlet at level 2, ratio → 2) whereas on a ring the SU(2)_L × SU(2)_R multiplet makes them degenerate at the critical point (`scripts/plot_low_levels.py --pbc`, `pbc_corr/true_singlet_gap.py`).
* **Filter-prepared periodic states** reproduce exact SzSz correlators, energies and entanglement entropy to within 0.1–3% at the sizes where the filter can reach a useful purification (`pbc_corr/`). Extrapolating the observable beyond those sizes is not done, and the pipeline is not a way to *derive* the field theory; it reproduces ED observables that ED already supplies.

## 7. 2D outlook **[modelled]**

A cost model (calibrated on the 1D designs: T = xπ/Δ with x ≈ 1.0–1.25, allowed Trotter budget ≈ 0.056√p_g, RZZ per step = 2 × Pauli terms) is combined with ED gaps on small clusters (up to 4×6), a fitted MPS bond-dimension growth and an assumed singlet-gap scaling (`crossover2d/`).

* The exact Trotter α per term grows slowly with cluster size (0.427, 0.472, 0.493, 0.505, 0.513 for 4×4 … 12×12).
* J2 = 0: first-order filter crosses MPS-circuit cost at L ≈ 15–16 for the hybrid variant A and L ≈ 11–12 for variant S (the two cost variants printed by `crossover2d/crossover.py`).
* J2 = 0.5, width 4, first-order filter, `inv_L` singlet model: filter/MPS cost ratio 562, 678, 591, 249, 66 at L = 4, 5, 6, 8, 10 for variant A, and 562, 439, 266, 63, 10.7 for variant S. The run was still going when this note was written, so larger L is not available.

These are projections from fitted inputs. None of the 2D numbers has been run through the pipeline, and the singlet-gap extrapolation is the dominant uncertainty. The second-order certificate of §4 has **not** been applied to the 2D estimate, although it would lower the filter side by roughly an order of magnitude on the 1D evidence.

## 8. Limits and open items

* The second-order bound is Proposition 16, Eq. (152) of Childs et al. (checked against the paper, §4.1); my implementation of it is validated numerically (two-term bound ratio ≤ 0.998 over 1,800 random cases; multi-term bound ≤ 0.28 on chains; all 12 built-circuit simulations within the bound).
* Second-order circuits are built and counted for N ≤ 12 only (12 cases); N > 12 is not computed.
* For N > 14 the trial overlap γ is a DMRG-MPS overlap and is not certified; those rows are conditional.
* Level-3 non-Clifford counts are accurate to about ±20%, not exact.
* `paper.tex` still reports CX counts; `results/RESULTS.*` is in the RZZ basis. Periodic and S² results are not yet in `results/RESULTS.*`.
* The 2D section is a model; the S² and symmetric-trial gains are measured in 1D only.

## 9. Where things are

| Item | Location |
|--|--|
| Pipeline, certificate, designs | `pipeline.py`, `floor.py`, `core/` |
| Results tables (RZZ basis) | `results/RESULTS.md`, `results/data/` |
| Proofs | `CERTIFICATION.tex`, `paper.tex` App. A–B |
| Second-order bound | `trotter2/` (`alpha2.py`, `compare.py`, `validate_*.py`, `table.md`, `README.md`) |
| S² projection and symmetric trial | `s2proj/` |
| Periodic chain, field-theory checks | `pbc_corr/` |
| 2D crossover model | `crossover2d/` |
| NetKet RBM energies | `rbm_netket/` |

## References (as cited above)

* A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, Theory of Trotter error with commutator scaling, Phys. Rev. X 11, 011020 (2021), arXiv:1912.08854 (equation not read; see §4.1).
* Stetcu–Baroni–Carlson filter paper (supplied as `a7585736-FIlter.pdf`) and arXiv:2211.10545, 2106.10867, 2412.14797, 2210.14109, 2006.02660, 2402.10362; Mizuta–Kuwahara, PRL 2025 — all abstract-level only.
* S. Eggert, Phys. Rev. B 54, R9612 (1996), SzSz observable definition (abstract as supplied).
