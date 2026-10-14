# Certification of the SBC filter pipeline — every chosen value, derived step by step

This document walks through every number the J1–J2 pipeline chooses or certifies, from the Hamiltonian
to the final Trotter step count `n`. Each step states (i) what is chosen, (ii) the derivation, (iii) what
it relies on, and (iv) where it lives in the code. Section 12 lists which inputs are rigorous and which are
assumptions, and Section 13 records what was simplified and what could still be dropped.

**Conventions.** `‖·‖` is the operator norm for matrices and the 2-norm for vectors. `|g⟩` is the exact
ground state, `|ψ⟩` the trial state, `γ = |⟨g|ψ⟩|²`.

**The claim being certified.** Given a trial state with ground-state weight `γ`, a spectral gap `Δ`, a
bandwidth `W`, and a Trotter step budget `n`, the post-selected output state `ψ̃` of the filter circuit
satisfies

```
1 − |⟨g|ψ̃⟩|²  ≤  ( √(2ℓ) + d_T )²                                           (★)
```

where `ℓ` is the *leakage* of the exact filter (Section 5) and `d_T` is the *Trotter distance*
(Section 8). Everything below is the derivation of (★) and of how its ingredients are chosen.

---

## 0. The pipeline in one picture

```
H (Pauli sum) ──► scale: H_s = (H − E0)/W            spec(H_s) ⊂ {0} ∪ [Δ, 1]            §1
trial state ψ, γ = |<g|ψ>|²                                                              §4
ancilla pulse  H·Rz(2φ)·exp(−i t H_s⊗Z)·H   ──►  block cos(H_s t + φ)                    §2
filter F(E) = Π_i cos(E t_i + φ_i)           P_succ ≥ γ F(0)²                           §3
η = sup_{[Δ,1]} |F| / |F(0)|  (certified)    leakage ℓ ≤ (1−γ)η²/(γ+(1−γ)η²)           §4–5
design (T = xπ/Δ, m pulses, t_i, φ_i)        maximize F(0)² s.t. η ≤ η_target            §6
Lie–Trotter, k_i steps/pulse                 ‖W̃_i − W_i‖ ≤ α t_i²/(2k_i), α = Σ‖[H_g, Σ_{g'>g}H_g']‖   §7–8
combine                                       1−F ≤ (√(2ℓ) + 2ε_T/√p_g)²,  ε_T = αT²/(2n)             §9
choose n and the leakage/Trotter split       minimize n / p_g                              §10
cost: CX, RZ, depth                           per-step transpile × n                      §11
```

---

## 1. The Hamiltonian, the shift and the bandwidth

**Model.** Open chain, `H = J1 Σ_i S_i·S_{i+1} + J2 Σ_i S_i·S_{i+2}`, `S = σ/2`. Writing each bond as
`¼(XX+YY+ZZ)`, `H = Σ_β c_β P_β` over Pauli strings `P_β` (plus an identity coefficient `c_I`, zero for plain Heisenberg bonds). There are
`3(N−1) + 3(N−2)` Pauli terms. The gate-level Hamiltonian lives in `hamiltonians.j1j2_hamiltonian`.

**Shift.** `E0` is the DMRG ground-state energy (Section 12 discusses its accuracy). We use
`H_s = (H − E0)/W`, so that the ground state sits at `0`.

**Chosen value — the bandwidth `W`.**

```
W = c_I + Σ_β |c_β| − E0
```

*Derivation.* For any Pauli sum, `‖H − c_I‖ ≤ Σ_β |c_β|‖P_β‖ = Σ_β |c_β|` (triangle inequality, Pauli strings have
norm 1). Hence `E_top ≤ c_I + Σ_β|c_β|`, so `W ≥ E_top − E0`, and therefore `spec(H_s) ⊂ [0, 1]` *without* needing
the DMRG value of `E_top`. The price is that `W` is 1.6–1.9× larger than `E_top − E0` on the J1–J2 chain
(measured), which shrinks the scaled gap by the same factor.

*Code:* `core/spectrum.py: get_spectrum` (`W_cert = (cI + l1) − E0`, `W = max(W_cert, Etop − E0)`).

**Chosen value — the scaled gap.** `Δ = (E1 − E0)/W`, with `E1` the lowest energy of any symmetry sector other
than the ground state (Section 12 explains why this matters near the BKT point `J2 ≈ 0.2411` [3,4]).

**What we need from this step.** The spectrum of `H_s` is contained in `{0} ∪ [Δ, 1]`. The upper end `1` is
rigorous (above). The lower structure (`0` is the non-degenerate ground state, `Δ` is the true gap) is an
*input*; see Section 12.

---

## 2. The pulse: why `H·Rz(2φ)·e^{−itH⊗Z}·H` applies `cos(Ht + φ)`

Take the ancilla in `|0⟩`. The pulse is the product of four operations, read right to left:

1. Hadamard: `|0⟩ → |+⟩ = (|0⟩+|1⟩)/√2`.
2. Controlled evolution `e^{−itH⊗Z}`, which equals `e^{−itH}` on the ancilla-`|0⟩` branch and `e^{+itH}` on the
   ancilla-`|1⟩` branch (because `Z|0⟩ = |0⟩`, `Z|1⟩ = −|1⟩`).
3. `Rz(2φ) = diag(e^{−iφ}, e^{+iφ})`, a phase `e^{∓iφ}` on the two branches.
4. Hadamard again.

Project the final ancilla on `⟨0|`: since `⟨0|H = ⟨+|`,

```
⟨0|_anc · (pulse) · |0⟩_anc = ½ ( e^{−iφ} e^{−iHt} + e^{+iφ} e^{+iHt} ) = cos(Ht + φ).
```

So reading `0` on the ancilla applies the operator `cos(Ht + φ)` to the system. The `1` outcome applies
`−i sin(Ht + φ)` and is discarded. `pulse convention` is checked numerically in
`core/circuits.py: verify_pulse_convention`, including that flipping `φ → −φ` is rejected.

*Citation:* this is the structure of the SBC projection filter [1] (and, with random times, of the rodeo
algorithm [9]).

---

## 3. The filter, the success path and the success probability

Chaining `m` pulses with the ancilla reset between them (or, in the early-abort circuit, stopping as soon as
the ancilla reads `1`), the all-zeros outcome applies

```
F(H_s) = Π_{i=1}^{m} cos(H_s t_i + φ_i),     F(E) = Π_i cos(E t_i + φ_i).
```

The *unnormalized* post-selected state is `a = F(H_s)ψ`, with `‖a‖² = P_succ`. Expanding `ψ = Σ_k ψ_k |k⟩` in
eigenstates `|k⟩` of `H_s` (`E_0 = 0`),

```
P_succ = Σ_k |ψ_k|² F(E_k)²  ≥  γ F(0)².
```

so the ground-state term alone gives a lower bound on the success probability. The expected number of
repetitions is at most `1/P_succ`.

**Code:** `core/filter_design.py: filter_values`, `core/simulate.py`.

---

## 4. Leakage: fidelity of the exact filter

Define the *worst-case contraction*

```
η = sup_{E ∈ [Δ, 1]} |F(E)| / |F(0)|.
```

If `spec(H_s) ⊂ {0} ∪ [Δ, 1]`, then the excited components satisfy `Σ_{k≥1} |ψ_k|² F(E_k)² ≤ (1−γ) η² F(0)²`.
The fidelity of the exact (un-Trotterized) filtered state with `|g⟩` is

```
|⟨g|ã⟩|² = γ F(0)² / (γ F(0)² + Σ_{k≥1}|ψ_k|² F(E_k)²)  ≥  γ / (γ + (1−γ)η²).
```

The *leakage* is the complementary quantity

```
ℓ = 1 − fidelity ≤ (1−γ)η² / (γ + (1−γ)η²).                                 (L)
```

**Inverting for the design target.** Requiring `ℓ ≤ ε` is equivalent to `(1−γ)η²(1−ε) ≤ γε`, i.e.

```
η_target(ε, γ) = √( γ ε / ((1−ε)(1−γ)) ).                                   (η*)
```

*Code:* `floor.py: eta_target`; `core/trotter.py: total_error_bound` (the `leak` line).

**Why a pure-state trace distance would not help here.** For pure states the trace distance is `sin θ` with
`cos²θ = fidelity`, i.e. `D = √(1 − F)`. That is exactly the relation we use in (★); Pinsker's inequality
relates *relative entropy* to trace distance and is vacuous for a pure target (relative entropy is infinite).
We therefore do not use it.

---

## 5. Certifying `η`: a rigorous bound on `sup|F|` over `[Δ, 1]`

`η` cannot be taken from a grid maximum alone, because the maximum between grid points may exceed it. We bound
it with the second derivative.

**Lemma.** `|F''(E)| ≤ L²` with `L = Σ_i |t_i|`.

*Proof.* Write `f_i(E) = cos(E t_i + φ_i)`. Then `|f_i| ≤ 1`, `|f_i'| ≤ |t_i|`, `|f_i''| ≤ t_i²`. By the product
rule,

```
F'' = Σ_i f_i'' Π_{j≠i} f_j  +  Σ_{i≠j} f_i' f_j' Π_{l≠i,j} f_l,
|F''| ≤ Σ_i t_i² + Σ_{i≠j} |t_i||t_j| = (Σ_i |t_i|)² = L².   ∎
```

The same argument gives `|F'| ≤ L`.

**Interval bound.** On an interval with midpoint `m` and radius `d`, Taylor's theorem with the remainder gives

```
|F(E)| ≤ |F(m)| + |F'(m)| d + ½ L² d².
```

Starting from a coarse grid, intervals whose upper bound is below `(1+rtol) × (current max)` are *pruned
(keeping that upper bound)*, and the rest are bisected. The returned supremum is the maximum of the pruned
upper bounds and the best grid value, so it is a rigorous upper bound on `sup_{[Δ,1]} |F|` (looser only if the
interval cap is hit, which is flagged by `converged = False`).

**Denominator.** `F(0) = Π cos φ_i` is exact for the DMRG reference energy. If the DMRG `E0` sits `δ` above the
true ground energy, the true ground state is at `−δ/W` in scaled units and `|F(−δ/W)| ≥ F(0) − L δ/W` (Lipschitz).
This is the `e0_slack` term: `f0_lb = max(F(0) − L·e0_slack, 0)` and `η = sup / f0_lb`. How `e0_slack` is
chosen is in Section 12.

*Code:* `floor.py: certify` (used during design), `core/filter_design.py: certify_filter` (grid max `+ L h/2`,
used on the snapped grid). `pipeline.make_grid` takes the smaller of the two certified values, both of which are
valid.

---

## 6. Choosing the filter: `T`, `m`, `t_i`, `φ_i`

**What is chosen.** The total time `T = Σ t_i = xπ/Δ`, the number of pulses `m`, the individual times `t_i`
and the phases `φ_i`.

**Optimization problem.** For given `(x, m)` and target `η_target`:

```
maximize  F(0)² = Π cos² φ_i
subject to  certified η ≤ η_target,   Σ t_i = T,   t_min ≤ t_i ≤ t_max,   |φ_i| < π/2.
```

Maximizing `F(0)²` is the right objective because `F(0)²` multiplies the ground-state probability: it is the
factor in `P_succ ≥ γF(0)²`. It is solved with SLSQP using analytic Jacobians, in two stages (find a feasible
point, then maximize), with cutting-plane refinement of the energy grid, followed by the certificate of Section 5.

**Why `T ∝ 1/Δ`.** A filter that is `1` (in relative terms) at `E = 0` and `≤ η` at `E = Δ` has to change by
`(1−η)F(0)` over an energy interval of length `Δ`. Since `|F'| ≤ L = T`,

```
(1−η) F(0) ≤ |F(0)| − |F(Δ)| ≤ T Δ        ⟹        T Δ ≥ (1 − η) F(0).
```

So `x = TΔ/π` cannot be taken arbitrarily small. The search covers `x ∈ {0.4, 0.5, 0.6, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5}`.

**Choices made by hand (not derived).**
- `m ∈ {4, 6, 8}` pulses. More pulses help only slightly at fixed `T` and add Trotter overhead per pulse.
- `t_min = 5%·T/m` so no pulse has zero duration, `t_max = T/2` (`FLOOR_KW`).
- The cost proxy for ranking `(x, m)` is `x²/F(0)²`: Trotter steps scale as `T²` (Section 8) and the
  repetition count as `1/F(0)²`. Ties within 2% go to fewer pulses.

**Grid snapping.** The pipeline runs in whole Trotter steps: pulse `i` gets an integer `k_i` of steps of equal
length `dt = T/n`, `t_i = k_i dt`, `Σ k_i = n` (largest-remainder rounding). Pulses with `k_i = 0` are dropped
(a zero-time pulse reduces to the scalar `cos φ_i` after post-selection, so `η` and the state are unchanged).
The phases are then re-optimized for the snapped times under the same constraint. A *uniform* `dt` is the
cost-optimal allocation: minimizing `Σ k_i` subject to `Σ α t_i²/(2k_i) ≤ ε` is a Cauchy–Schwarz problem
whose solution is `k_i ∝ t_i`.

*Code:* `floor.py: solve_floor, floor_vs_precision, grid_design`.

---

## 7. Trotter error of one controlled evolution

Each pulse needs `e^{−itH_s⊗Z}`, implemented with `k` first-order Lie–Trotter steps. Write `H_s = Σ_g H_g` with
`H_g = c_g P_g`.

**Lemma (first-order commutator bound).** For terms ordered `g = 1…Γ` and `S(τ) = Π e^{−iτH_g}` (in the order
that the circuit applies them),

```
‖e^{−iτH} − S(τ)‖ ≤ (τ²/2) · Σ_g ‖[H_g, Σ_{g'>g} H_{g'}]‖ =: (τ²/2) α.
```

*Proof for two terms.* Let `V(s) = e^{−isA}e^{−isB}`. Then `V' = −i(A+B)V − i[e^{−isA}, B]e^{−isB}`, so by
Duhamel `e^{−iτ(A+B)} − V(τ) = i∫₀^τ e^{−i(τ−s)(A+B)} [e^{−isA}, B] e^{−isB} ds`. Also
`[e^{−isA}, B] e^{isA} = e^{−isA}Be^{isA} − B`, whose norm is at most `s‖[A,B]‖` (differentiate in `s`).
Hence `‖e^{−iτ(A+B)} − V(τ)‖ ≤ ∫₀^τ s‖[A,B]‖ ds = (τ²/2)‖[A,B]‖`. For `Γ` terms, peel off `H_1` against the
sum of the others and induct, which produces `Σ_g ‖[H_g, Σ_{g'>g} H_{g'}]‖`. ∎ (This is the standard result;
see Childs, Su, Tran, Wiebe, Zhu [5], and the survey [6].)

**`k` steps.** With `τ = t/k` and telescoping over the `k` steps (all factors are unitary),

```
‖e^{−itH} − S(t/k)^k‖ ≤ k · (t/k)² α / 2 = α t² / (2k).
```

**The controlled version has the same `α`.** The circuit uses `H⊗Z`, whose terms are `c_g P_g ⊗ Z`. Because
`Z² = 1`, `[P_g⊗Z, P_{g'}⊗Z] = [P_g, P_{g'}] ⊗ 1`, so the commutators, and hence `α`, are those of `H_s`
itself. Therefore

```
‖W̃_i − W_i‖ ≤ ε_i := α t_i² / (2 k_i),
```

where `W_i = e^{−it_iH_s⊗Z}` and `W̃_i` is its Trotterized version.

**Choosing `α` — and the simplification.** The definition is just

```
α = Σ_g ‖[H_g, Σ_{g'>g} H_{g'}]‖.                                           (α)
```

No further inequality is needed: for Pauli strings, `[P_g, P_{g'}] = 2P_gP_{g'}` if they anticommute and `0`
otherwise, and `P_g` is unitary, so

```
‖[H_g, Σ_{g'>g} H_{g'}]‖ = 2 |c_g| · ‖ Σ_{g'>g, anticommuting} c_{g'} P_{g'} ‖.
```

The inner sum acts on the few qubits around term `g` (a handful of qubits for a J1–J2 bond), so its norm is a
tiny dense eigenvalue problem. The cost is `O(#terms)`, independent of `N`, and the value is *exact*. It agrees
with the exact sparse computation to `1e-14` for `N = 4, 6, 8` and `J2 ∈ {0, 0.2411, 0.4}` in both term orders,
and runs in 16 s at `N = 64`. The previous `N > 14` fallback, the triangle bound `Σ_{pairs} 2|c||c'|`,
overestimates `α` by up to 1.25× (measured `N=16..64`: 28.1 vs 22.4, 61.3 vs 48.5, 127.5 vs 100.6 at `J2 = 0.4`)
and is no longer needed. `α` grows linearly with `N` (22.4, 48.5, 100.6 for `N = 16, 32, 64`).

**Term order.** The bound holds for the order the circuit realizes. The pipeline takes the larger of the two
possible orders (`trotter_alpha`), and `lie_trotter_order` checks numerically which one qiskit's
`LieTrotter(preserve_order=True)` realizes.

*Code:* `core/trotter.py: alpha_commutator` (new, exact, local), `alpha_comm(tight=True)` (dense/sparse check),
`trotter_alpha`.

---

## 8. From pulse error to state error

Let `A_i = ⟨0|_anc W'_i |0⟩_anc` where `W'_i` is the full pulse of Section 2 and `Ã_i` the Trotterized
counterpart. They differ only in the `e^{−itH⊗Z}` factor, and the Hadamards and `Rz` are unitary, so
`‖A_i − Ã_i‖ ≤ ‖W̃_i − W_i‖ ≤ ε_i`. (Compressing a unitary to the `⟨0|·|0⟩` block cannot increase a norm.)

**Telescoping.** All `A_i` and `Ã_i` are blocks of unitaries, so `‖A_i‖, ‖Ã_i‖ ≤ 1`. Hence for the unnormalized
exact vector `a = A_m⋯A_1ψ` and the Trotterized vector `b = Ã_m⋯Ã_1ψ`,

```
‖a − b‖ ≤ Σ_i ‖A_i − Ã_i‖ ≤ Σ_i α t_i²/(2k_i) =: ε_T.
```

On the uniform grid, `t_i = k_i dt`, `Σ k_i = n`, `T = n dt`:

```
ε_T = Σ_i α k_i dt² / 2 = α dt² n / 2 = α T² / (2n).                          (ε_T)
```

**Normalization.** The pipeline compares the *normalized* states `a/‖a‖` and `b/‖b‖`:

```
‖a/‖a‖ − b/‖b‖‖ ≤ ‖a − b‖/‖a‖ + | ‖b‖ − ‖a‖ |/‖a‖ ≤ 2 ‖a − b‖ / ‖a‖.
```

With `‖a‖² = P_succ ≥ p_g := γ f0_lb²` (Section 3 and the denominator of Section 5),

```
d_T := ‖ψ̃ − ψ_f‖ ≤ 2 ε_T / √p_g = α T² / (n √p_g).                           (d_T)
```

where `ψ_f = a/‖a‖` is the exact filtered state and `ψ̃ = b/‖b‖` the Trotterized one (the global phase is
irrelevant to fidelity). Solving for the step count: for a target `d_T`,

```
n = ⌈ α T² / ( √p_g · d_T ) ⌉.                                               (n)
```

*Code:* `core/trotter.py: trotter_bounds, steps_needed`.

---

## 9. Combining leakage and Trotter error

Let `g` be the exact ground state and `θ` a phase aligning `ψ_f` with `g`.

**Leakage as a distance.** `‖ψ_f − e^{iθ}g‖² = 2 − 2|⟨g|ψ_f⟩|`. Since `|⟨g|ψ_f⟩|² ≥ 1 − ℓ` and
`1 − √(1−ℓ) ≤ ℓ`, this is at most `2ℓ`. So `‖ψ_f − e^{iθ}g‖ ≤ √(2ℓ)`.

**Triangle inequality.** `‖ψ̃ − e^{iθ}g‖ ≤ ‖ψ̃ − ψ_f‖ + ‖ψ_f − e^{iθ}g‖ ≤ d_T + √(2ℓ) =: r`.

**Back to infidelity.** If `‖ψ̃ − e^{iθ}g‖ ≤ r` then `|⟨g|ψ̃⟩| ≥ 1 − r²/2`, hence
`1 − |⟨g|ψ̃⟩|² ≤ r² − r⁴/4 ≤ r²`. Therefore

```
1 − F_total ≤ ( √(2ℓ) + d_T )².                                              (★)
```

This is the whole bound: one triangle inequality and two squarings. Nothing here depends on the Hamiltonian
beyond `α`, `Δ`, `W`, `γ`.

*Code:* `core/trotter.py: total_error_bound` (`total = min(1, (√(2·leak) + d_T)²)`).

---

## 10. Choosing the step count `n` and the leakage/Trotter split

Given a total target `ε`, require `(√(2ℓ) + d_T)² ≤ ε`, i.e. `d_T ≤ √ε − √(2ℓ)`. Substituting (n):

```
n(ε, ℓ) = ⌈ α T² / ( √p_g · (√ε − √(2ℓ)) ) ⌉.
```

A looser leakage target (large `ℓ`) lets the filter have a larger `F(0)²` and a shorter design, but leaves less room
for Trotter error, so the optimum is interior. The pipeline therefore tries the leakage budgets
`ℓ-budget = f·ε`, `f ∈ {0.4, 0.2, 0.1, 0.05, 0.02, 0.01}` (`FLOOR_EPS_FRACS`), builds the best design for each, and picks
the one minimizing `n/p_g`. This is the *expected* cost, since the expected number of repetitions is `1/P_succ`
and the circuit length is proportional to `n`.

The final guaranteed `n` is then found on the *snapped grid* by increasing `n` in 5% steps until (★) holds
for the actual `(k_i, φ_i)`, because snapping changes `η` and the `t_i`.

*Code:* `pipeline.py: choose_floor_design, find_guaranteed`.

**Empirical (oracle) `n`.** Where the exact ground state is available (`N ≤ 14`), the same circuit is
simulated and the smallest `n` whose measured infidelity is `≤ ε` is also reported. That number is not
guaranteed (the measured error is not monotone in `n`) and is reported only as the practical cost.

---

## 11. Resource costs

**Gate counts.** One Trotter step of one pulse is transpiled once (opt-level 3, basis `{h, s, cx, rz}`) and its
CX count and depth multiplied by `n`; the `m` ancilla `Rz(2φ_i)` and the trial-state circuit are added. This was
validated against a full transpile of the entire circuit (0.0% CX difference at `n = 48`).

**Non-Clifford `RZ` and `T`.** An `rz` whose angle is a multiple of `π/2` is Clifford and free; the rest are
counted. The `T`-count of synthesizing them is approximated by `1.15 log₂(1/ε_syn) + 9.2` per rotation; the
`1.15` slope is from the repeat-until-success synthesis analysis [10] and the additive `9.2` is **not verified**
(see Section 12). Ross–Selinger [11] gives `≈ 3 log₂(1/ε)` ancilla-free.

**Expected cost.** Dividing by the success probability `p_g` gives the expected CX count including repeats.

*Code:* `core/resources.py`, `pipeline.per_step_costs, costs, measured_cost`.

---

## 12. What is rigorous and what is an assumption

| Ingredient | Status | Notes |
|---|---|---|
| `spec(H_s) ≤ 1` | **Rigorous** | triangle inequality on the Pauli sum (§1) |
| Pulse implements `cos(Ht+φ)` | **Rigorous** | algebra (§2), also checked numerically |
| `P_succ ≥ γF(0)²`, leakage (L) | **Rigorous** given `γ`, `Δ`, ground state non-degenerate | §3–4 |
| `η` upper bound | **Rigorous** given the window `[Δ, 1]` | second-derivative branch and bound (§5) |
| Trotter bound `α t²/(2k)` | **Rigorous** | commutator bound, exact `α` (§7) |
| Combination (★) | **Rigorous** | §9 |
| `γ` | **Assumption** | DMRG/MPS overlap `|⟨ψ_DMRG|ψ⟩|²`, not a lower bound on `|⟨g|ψ⟩|²`. |
| `Δ` (the gap) | **Assumption** | DMRG `E1 − E0`. A variational excited energy is an *upper* bound on `E1`, which is the unsafe direction for the certificate; no rigorous lower bound on `E1` is available. |
| `E0` shift (`e0_slack`) | **Heuristic, conservative in practice** | see below |
| `T`-count constants | **Unverified** | `+9.2` |

**`γ` and `Δ` are the two inputs the guarantee depends on.** The claim (★) is therefore a
*conditional* certificate: it holds exactly when the DMRG gap is not larger than the true gap and the trial
overlap is not larger than the true overlap. With DMRG energies accurate to `~1e-13` (validated against exact
diagonalization for `N = 8–20` including the dimerized phase) these assumptions are very likely to hold in the
cases we run, but they are assumptions. Practical mitigations: use `Δ_used = (1−ρ)Δ_DMRG` and `γ_used = (1−ρ)γ_DMRG`
with a safety margin `ρ`; the sensitivity of `n` to wrong `Δ`, `γ` is in the sensitivity study.

**The `E1` sector.** The first excitation of the J1–J2 chain is a triplet for `J2 < J2c ≈ 0.2411` [3,4] but can be
a singlet in the dimerized phase, so `E1` is taken as the lowest state of any symmetry orthogonal to `ψ0`
(penalty `λ|ψ0⟩⟨ψ0|` applied as a rank-1 local term, `get_energies.py: OrthoDMRG2`), never from an `S^z`-sector
shortcut.

**`e0_slack`.** The DMRG `E0` is variational, so `E0_DMRG ≥ E0_true`; the true ground state is slightly below
`0` in scaled units. The code uses `e0_slack = √Var(H)/W` with `Var = ⟨H²⟩ − ⟨H⟩²`. A derivable bound is the
Temple inequality [12]: for a trial energy `⟨H⟩` below the first excited energy,
`E0 ≥ ⟨H⟩ − Var/(E1 − ⟨H⟩)`, giving `δ ≤ Var/Δ_raw` (in scaled units `Var/(W² Δ)`). Whenever `√Var < Δ_raw`
the code's `√Var/W` is *larger* than the Temple value, so the code is conservative, but note that
Temple needs the same gap input `E1` as everything else. At DMRG accuracies of order `1e-13` in energy the variance is tiny (the actual `e0_slack` of each run is stored in the
result data), so `L·e0_slack` is negligible against `F(0)`; this was not tabulated here.

---

## 13. Simplifications made and still possible

**Done.**
1. **`α`.** Replaced the triangle-inequality bound (and its fallback logic) by the definition (α), evaluated exactly on
   the local support of each term: one code path for all `N`, exact, linear cost (§7). *Not yet wired into
   `pipeline.build_ctx`*: the large-`N` result tables currently in `RESULTS.md` still use the triangle value
   (a factor up to 1.25 larger in `n`, still valid); switching it reduces the guaranteed `n` at `N > 14` by that factor.
2. **Pinsker/trace-distance machinery** is not used; pure-state distance `√(1−F)` is enough (§4, §9).
3. **One chain of three inequalities** (leakage, telescoping, triangle) instead of separate lemmas for each effect.

**Could still be dropped or merged (not changed).**
- The two interval certificates (`floor.certify` and `filter_design.certify_filter`) are both valid; keeping only the
  branch-and-bound one removes a code path (the min of the two is used today).
- The 2-norm form (★) carries a factor 2 from normalization (§8). A tighter statement uses the *trace distance*
  `D ≤ √ℓ + ε_T/√p_g`; (★) is simpler and, for the cases where `d_T` is small, the same order. Switching would shrink
  `n` by up to 2× at the cost of one more line of algebra.
- `e0_slack` can be replaced by the Temple bound or simply dropped (it is `~1e-6` of the leading terms).
- The hand-chosen design grid (`x`, `m`, `t_min`, `t_max`) could be replaced by a single scalar search over `T`.

---

## References

[1] I. Stetcu, A. Baroni, J. Carlson, *Projection algorithm for state preparation on quantum computers*,
Phys. Rev. C **108**, L031306 (2023), arXiv:2211.10545.

[3] T. Okamoto, K. Nomura, *Fluid-dimer critical point in S=1/2 antiferromagnetic Heisenberg chain with next
nearest neighbor interactions*, Phys. Lett. A **169**, 433 (1992). *(From memory.)*

[4] S. Eggert, *Numerical evidence for multiplicative logarithmic corrections from marginal operators*,
Phys. Rev. B **54**, R9612 (1996), arXiv:cond-mat/9602026 (reports `J2c = 0.241167 ± 0.000005`).

[5] A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, *Theory of Trotter error with commutator scaling*,
Phys. Rev. X **11**, 011020 (2021), arXiv:1912.08854.

[6] S. Lloyd, *Universal quantum simulators*, Science **273**, 1073 (1996); H. F. Trotter, *On the product of
semi-groups of operators*, Proc. Amer. Math. Soc. **10**, 545 (1959).

[7] S. R. White, *Density matrix formulation for quantum renormalization groups*, Phys. Rev. Lett. **69**, 2863 (1992);
U. Schollwöck, *The density-matrix renormalization group in the age of matrix product states*, Ann. Phys. **326**, 96 (2011).

[8] J. Gray, *quimb: a python package for quantum information and many-body calculations*, J. Open Source Softw. **3**, 819 (2018).

[9] K. Choi, D. Lee, J. Bonitati, Z. Qian, J. Watkins, *Rodeo algorithm for quantum computing*,
Phys. Rev. Lett. **127**, 040505 (2021). *(From memory.)*

[10] A. Bocharov, M. Roetteler, K. M. Svore, *Efficient synthesis of universal repeat-until-success quantum circuits*,
Phys. Rev. Lett. **114**, 080502 (2015).

[11] N. J. Ross, P. Selinger, *Optimal ancilla-free Clifford+T approximation of z-rotations*,
Quantum Inf. Comput. **16**, 901 (2016), arXiv:1403.2975.

[12] G. Temple, *The theory of Rayleigh's principle as applied to continuous systems*,
Proc. R. Soc. A **119**, 276 (1928).

[2] M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information*, Cambridge Univ. Press (2000) — fidelity
and trace distance of pure states.

Items marked "from memory" were not checked against the source in this session; [1], [4] and [5] were confirmed by
search. Note that [12] and the J2c value attribution to [3] should be verified before submission.
