# Dephasing, mixing and "how many states" in the Trotterised filter

Scope: open J1–J2 chains, N = 6–10 (N = 12 for the level counts, pending), J2 ∈ {0, 0.2411, 0.4}, ε = 10⁻² guaranteed designs, exact diagonalisation. Everything below is computed by the scripts in this directory (`dephasing.py`, `mixing_theory.py`, `scaling_check.py`, `analyze_dephasing.py`); energies are in the scaled units of the filter (spectrum in [0, 1]).

## 1. Set-up in the eigenbasis

Write the trial state as ψ = Σ cₙ|n⟩ in eigenstates of H_s, with levels Eₙ. The exact filter is diagonal: aₙ = cₙ F(Eₙ), F(E) = ∏ᵢ cos(E tᵢ + φᵢ). The exact filter is diagonal in this basis: it reweights amplitudes and mixes nothing. All Trotter error is therefore a statement about how the Trotterised pulse differs from cos(H_s t + φ) *in this basis*.

## 2. How many states are there? (`dephasing_table.md`, `dephasing.png`)

Counting **levels** (distinct energy and total spin S; a spin-S level has 2S+1 degenerate Sz components):

| N | J2 | levels | γ | levels with >10⁻⁴ of the trial weight | after the filter | levels for 99.99% of the weight (trial → filtered) |
|--|--|--|--|--|--|--|
| 6 | 0 / 0.2411 / 0.4 | 20 | 0.90 / 0.95 / 0.99 | 14 / 14 / 7 | 3 / 3 / 4 | 15→3 / 14→6 / 9→4 |
| 8 | 0 / 0.2411 / 0.4 | 70 | 0.84 / 0.91 / 0.98 | 27 / 18 / 12 | 5 / 1 / 3 | 40→5 / 31→4 / 17→5 |
| 10 | 0 / 0.2411 / 0.4 | 252 | 0.78 / 0.87 / 0.96 | 49 / 43 / 18 | 5 / 3 / 3 | 106→6 / 127→5 / 47→9 |

So out of 2ᴺ basis states (64–1024) and 20–252 distinct levels, the trial state occupies tens of levels at the 10⁻⁴ level, with a weight spread over many decades (see the figure). The filter reduces that to a handful. The weight above half the gap is 1–22% (larger for smaller J2).

## 3. What a Trotter step does to those states

A Lie–Trotter step (later terms applied on the left) is exactly e^{−iτ H_eff} with, by BCH,

H_eff = H + τ X₁ + τ² X₂ + …,  X₁ = −(i/2) Σ_{a>b} [H_a, H_b].

All Pauli terms XX, YY, ZZ are real matrices, so every commutator is real antisymmetric and **X₁ is imaginary Hermitian**. In a real eigenbasis it has zero diagonal. Checked: max |diag X₁| = 4×10⁻¹⁸ (N = 6, J2 = 0.4), and X₁ equals the numerically extracted (H_eff − H)/τ to 0.1% at τ = 0.02, 0.01.

Consequences:

* **No O(τ) dephasing.** Energy (quasi-energy) shifts start at O(τ²). Measured from the eigenvalues of one step: ground-state shift exponent 2.00 for **both** first and second order, for all N ≤ 10 and J2. Magnitudes: N = 10, J2 = 0.4, τ = 0.05: δE₀ = 1.8×10⁻⁶ (first order) and 6.7×10⁻⁷ (second order). At the certified step sizes (τ ≈ 0.009 for first order, ≈ 0.38 for second order at this case, T ≈ 25), extrapolating the measured τ² law, the ground-state phase T·δE₀ is ~10⁻⁶ rad (first order) and ~10⁻³ rad (second order, which uses a much larger τ). For first order this is negligible next to the measured error (2×10⁻⁴–10⁻³ in trace distance). For second order I did not separate the phase contribution from the mixing contribution, so I make no claim about which dominates there.
* **The O(τ) error is level mixing.** Treating τX₁ as a perturbation (first order, non-degenerate levels), the amplitude of level m generated from level n over a pulse of length t is
  δV₀_mn = −τ X_mn (e^{−iEₙt} − e^{−iE_m t}) / (E_m − Eₙ),
  with the ancilla-1 branch V₁ = V₀* giving the mirror term, and the pulse block Ã = (e^{−iφ}V₀ + e^{iφ}V₁)/2. Its magnitude is bounded by 2τ|X_mn| / |E_m − Eₙ|, **independent of t** for non-resonant levels. The commutator bound αt²/(2k) = αtτ/2 instead grows linearly in t and uses a worst-case operator norm over all states.
* **Test.** Propagating these first-order errors through all pulses (exact cosines after each) and removing the overall post-selection amplitude, the predicted trace distance to the exact filter state matches the measured first-order Trotter error to within 0.4% in all 9 cases (`mixing_table.md`, ratio 0.996–1.003). The first-order error is therefore, to this accuracy, pure first-order mixing.
* **Scaling.** First-order measured error ∝ 1/n exactly (fitted exponent −1.00 over n = 98…1576, N = 6, J2 = 0.4). Second order: the argument is that the Strang step is time-symmetric, S(−τ) = S(τ)⁻¹, so H_eff has only even powers of τ, there is no τX₁ term, the mixing is O(τ²) and the error should go as 1/n². Measured (N = 6, J2 = 0.4 only): factors of 4 per doubling, 5.7×10⁻⁴ → 1.4×10⁻⁴ → 3.6×10⁻⁵ at n = 20, 40, 80; the fit over all five points is −2.3 because of the pre-asymptotic n = 6, 11 points).

## 4. Which states carry the error

From the perturbation theory (`mixing_N*_J2*.json`):

| N | J2 | measured trace distance | predicted | eigenstates carrying 99% of the error | share of the error from the last two pulses |
|--|--|--|--|--|--|
| 6 | 0 / 0.2411 / 0.4 | 5.1e-4 / 6.4e-4 / 9.7e-4 | ratio 1.000 / 1.003 / 0.996 | 11 / 4 / 2 | 0.97 / 0.92 / 0.97 |
| 8 | 0 / 0.2411 / 0.4 | 3.8e-4 / 3.5e-4 / 9.3e-4 | 1.000 / 1.000 / 1.000 | 10 / 4 / 2 | 0.98 / 0.88 / 0.97 |
| 10 | 0 / 0.2411 / 0.4 | 3.1e-4 / 2.0e-4 / 5.3e-4 | 1.000 / 0.999 / 1.000 | 27 / 14 / 4 | 0.96 / 0.84 / 0.93 |

(Eigenstates here are individual eigenvectors, i.e. Sz components of a multiplet count separately.) So, although the trial state has tens of levels with weight, the first-order Trotter error is generated by only 2–27 source eigenstates (of 64–1024), and 84–98% of it is generated in the **last two pulses**, the longest ones. The number of source eigenstates is larger for smaller trial overlap γ (J2 = 0: 11, 10, 27 for N = 6, 8, 10; J2 = 0.2411: 4, 4, 14; J2 = 0.4: 2, 2, 4). The error resolved by level is the right-hand panel of `dephasing.png` (N = 10, J2 = 0.4 only).

Where does the second-order error sit? In the 9 cases, 8–96% of the state error² is in levels that the exact filter suppresses (`dephasing_table.md`, "share in suppressed levels"), i.e. leakage into excited levels rather than error on the ground-state amplitude. The ground level itself carries at most 8×10⁻¹² (first order) and 8×10⁻¹¹ (second order) of the state error², against totals of 4×10⁻⁸ to 9×10⁻⁷ (first order) and 8×10⁻⁹ to 3×10⁻⁷ (second order), i.e. under 10⁻³ of the error (the global post-selection amplitude is removed before this split).

## 5. What this says about the certificate

* The certified Trotter trace-distance budget in the designs is d_T ≈ 0.056; the measured first-order Trotter error is 2×10⁻⁴ to 1×10⁻³, i.e. 60–280× smaller, which is of the order of the certified-versus-oracle step ratios (66–438× at first order).
* The commutator bound αt²/(2k) grows linearly with the pulse length t and is a worst case over all states; the first-order perturbative expression of §3, which matches the measured error, has no growth with t for non-resonant levels.

## 6. Limits and open points

* ED only (N ≤ 10 for the mixing theory and Floquet check; N = 12 level counts pending). The number of states contributing to the error is **not** shown to scale in any particular way with N: it depends on J2 and γ and is non-monotonic over N = 6–10.
* First-order perturbation theory reproduces the measured error (ratio 0.996–1.003) in the 9 cases tested without any secular (t-linear) term. X₁ commutes with S² and S^z_tot to ≤ 3×10⁻¹⁶ (N = 6, 8; `check_x1_symmetry.py`), so it connects only states with the same S and S^z.
* One design per case (ε = 10⁻², guaranteed), step counts from the certificate; second-order pulses were not re-optimised.
* The second-order mixing theory (O(τ²)) was not implemented; for second order only the scaling (1/n², measured at N = 6) and the level-resolved error are shown.
