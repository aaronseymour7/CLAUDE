# Response to the second set of review comments

Every change is on branch `claude/project-thread-6thtun` of the CLAUDE repo (PR #7). Line numbers below refer to `j1j2_filter/sbpipe_study/paper/paper.tex` as it now stands. After the changes, all 96 simulated circuits still respect the bound and the 10 tests still pass.

## 1. Cost per success uses p_g, but the theorem gives (√p_g − ε_U)²: **valid, fixed**

The comment is right. Dividing by p_g underestimated the certified expected cost.

- `pipeline._plan_for` now divides the cost of one attempt by the certified bound (√p_g − ε_U)². This value is also stored as `p_succ_bound`. When no filter is used, the success probability is 1, because nothing is post-selected.
- Sec. II E states the new denominator. The pipeline spends the whole error budget, ε_U = (√ε − √ℓ)√p_g, so the denominator equals p_g(1 − √ε + √ℓ)², which is at least 0.81 p_g at ε = 10⁻². Costs therefore rise by a factor between 1.06 and 1.23.
- I reran every study. Plan selection uses the new cost, so the chosen plans can change.

Numbers that changed:

| Quantity | Before | After |
|--|--|--|
| N = 48, χ = 2, CNOTs | 4.8×10⁸ | 5.4×10⁸ |
| N = 48, χ = 4, CNOTs | 8.8×10⁶ | 9.7×10⁶ |
| N = 48, χ = 4, T gates | 1.9×10⁸ | 2.1×10⁸ |
| Median looseness of the bound | 17× | 14× |
| Dyadic median leak ratio | 56 | 54 |

The medians of 92× (second order against first order) and 2.9× (optimized against dyadic), the 96/96 validation result and every qualitative conclusion are unchanged.

## 2. Π₀ may be a degenerate projector, and 1 − F ≤ D² for mixed states: **partly valid, fixed by restating**

The proof already measured distance to the pure state |g⟩ = Π₀|ψ⟩/‖Π₀|ψ⟩‖, and in a noiseless circuit the output is pure. The notation D(ρ̃, Π₀) hid both facts.

Theorem 1 now:
- defines |g⟩ explicitly, and notes that it is pure even when the ground space is degenerate;
- states that, without noise, the post-selected output is pure;
- states 1 − Tr(Π₀ρ̃) ≤ 1 − |⟨g|ψ̃_out⟩|² = D².

The proof adds that last step, which holds because |g⟩ lies in the ground space. It also says that a mixed output, for example under noise, is outside the theorem. The paper defines 1 − F as 1 − |⟨g|ψ̃_out⟩|². For the chains studied here (even N, open boundaries) the ground state is nondegenerate, which ED confirms for N ≤ 14, so this is the usual infidelity.

## 3. W scaling of the Trotter constants: **not an error, convention now stated**

The code evaluates α₁ and β₂ directly on the terms of H_s, with coefficients c_a/W (`pipeline.trotter_constant(Hs, …)`). It uses the scaled times t_i throughout. I checked numerically that α₁(H_s)·W² / α₁(H) = 1 and β₂(H_s)·W³ / β₂(H) = 1, to 10⁻¹⁵. So no factor of W is missing, and the resource estimates are unaffected.

Sec. II C now says which convention the code uses and gives the equivalent physical form: α₁(H)/W² and β₂(H)/W³, with evolution time t_i/W.

## 4. The gap is not rigorous, and the relevant level is the one the trial state can reach: **agreed, text strengthened**

This limitation was already acknowledged. The text now spells it out:

- **Sec. II A:** E₁ is defined as the lowest level the trial state can overlap with. The DMRG trial states are not symmetry-adapted, so this is the global first excited level. Restricting to the trial state's symmetry sector could only raise it, so using the global level is conservative. The certificate applies to the actual Hamiltonian only if the assumptions hold for it. With an overestimated E₁, a real excitation can sit in the filter's transition region, and the bound then says nothing.
- **Sec. III A:** the lowest triplet is used as the global first excitation, which ED confirms for every case with N ≤ 14. DMRG energies are variational upper bounds, so a residual DMRG error biases the gap upward, as does a penalty run that converges to a higher state. That is the dangerous direction. The bond-dimension error bars do not protect against this bias, so for N > 14 the certificate is conditional on the gap.

A rigorous lower bound on the gap is outside this paper. Two items from the first review remain open. Retitling the paper "conditionally certified" goes further in this direction. Widening the gap-overestimation tests to more instances is still to do. Both are paused until you say to go ahead.

## 5. Dirichlet identity and interval bound: **valid, fixed**

Sec. II B now separates the two steps:
1. The identity is exact for every E. It follows from repeated use of sin 2x = 2 sin x cos x.
2. The suppression bound uses the window: for E ∈ [Δ, 1] and τ = π/(1+Δ), Eτ ∈ [θ, π−θ]. It depends on the scaling H_s and the choice of τ, and does not hold outside [Δ, 1].

The uniform-schedule paragraph states the same range Et ∈ [θ, π−θ] for its bound on cos.

## 6. Optimality of k_i = ⌈|t_i|/δt⌉: **valid, assumptions now stated**

Sec. II C now explains that the error constants are equal across pulses by construction:
- every pulse evolves under the same H_s ⊗ Z_a;
- the phase enters only the ancilla rotation;
- the bound is symmetric in the sign of t.

With equal constants, the continuous relaxation is solved by k_i ∝ |t_i| for both orders. Rounding up keeps the bound and costs at most m extra steps, but it need not be the integer optimum. Any budget freed by rounding goes to rotation synthesis.

## 7. Report synthesis tolerance and rotation count with the T count: **done**

Tables II and IV now have two new columns:
- **rot.:** arbitrary rotations per attempt (trial plus filter).
- **δ:** the synthesis accuracy per rotation.

Sec. II E says this explicitly. In Table II, δ ranges from about 10⁻¹¹ to 10⁻⁷.

## Not done (paused)

From the first review:
- rigorous norms in place of the Lanczos estimate for 9–12-qubit supports (you declined that edit);
- the energy–gap (Eckart) certificate baseline;
- a QSVT or qubitization cost comparison;
- the isometry costing of MPS circuits;
- a ladder model;
- broader sensitivity tests;
- the title and abstract reframing.

## Follow-up comments

1. **"≥ (1−ε)² p_g" in Sec. II E.** The source already reads `\ge(1-\sqrt\epsilon)^2p_g`, which is (1−√ε)² p_g, and that matches the 0.81 figure. The radical probably got lost when the PDF was read as text. No change is needed.
2. **The overlap condition deserves equal prominence.** Agreed, and done.
   - The abstract now says the certificate is conditional on the DMRG energies and on an overlap estimated from a reference MPS.
   - Sec. III A states that γ estimates the overlap with the exact ground state but is not a lower bound on it.
   - The Discussion names all three conditions: the E₀ interval, the gap and the overlap.
   - The Discussion also notes that the Eckart inequality 1−γ ≤ (⟨H⟩−E₀)/(E₁−E₀) would turn energy bounds into a lower bound on γ.
3. **Title.** Changed to "From classical energy estimates to a conditionally certified ground state: …".
