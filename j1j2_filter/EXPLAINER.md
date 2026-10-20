# Explainer: level crossings, movers, towers, and what the filter evolves

Written for talking through the J1–J2 chain results. Numbers marked (ours) are from this repository; the field-theory statements are standard SU(2)₁ results, not new.

## 1. Why periodic chains cross and open chains do not

**Terminology.** In open boundary conditions (OBC) there is no avoided crossing: the levels never cross. The first excitation is a triplet (S = 1) and the competing level is a singlet (S = 0). Total spin is conserved, so the two sectors cannot hybridize, and a crossing between them would be a true crossing. The claim is: periodic chains (PBC) have a true singlet–triplet crossing near J2 ≈ 0.2411, and open chains do not.

**Evidence (ours).** First excitation is a triplet for open chains N ≤ 16 and J2 ∈ [0, 1] (`scripts/check_first_excitation.py`, `results/low_levels*.json`). The PBC crossing converges toward 0.2411 with N (`results/low_levels_ed_pbc.json`, `pbc_corr/true_singlet_gap.py`).

**Short version.** At the critical point J2c ≈ 0.2411 the chain is described by an SU(2)₁ conformal field theory. On a ring, left- and right-moving sectors both carry spin, so the lowest triplet and lowest singlet come from the same multiplet and are degenerate in the thermodynamic limit. A marginal operator, coupling J_L·J_R, splits them, and its coupling changes sign at J2c, so the levels cross on a ring. An open chain has a boundary that ties left- and right-movers into one chiral sector. There is only one tower of states, with the triplet at level 1 and the singlet at level 2, so the singlet/triplet ratio tends to 2, not 1. There is nothing degenerate to cross.

**Details.**
- Ring: lowest triplet and singlet sit in the (½, ½) primary multiplet of SU(2)_L × SU(2)_R; at the fixed point the gap ratio → 1.
- Marginal operator J_L·J_R splits the multiplet with coupling ~ 1/ln L; the finite-size crossing J2c(L) drifts to 0.2411 with 1/ln L corrections.
- Open chain: the boundary identifies left- and right-movers, so only the vacuum Kac–Moody tower survives. Level 1 is J₋₁ᵃ|0⟩ (a triplet); the singlet first appears at level 2; energies go as E ∝ (πv/L)·level, so the ratio is 2. The marginal operator shifts levels by ~1/ln L and cannot reorder them.

**Caveat.** For J2 well above J2c (dimerized phase) open chains have edge spin-½ states and the singlet and triplet become close in energy. We checked only N ≤ 16, so no claim is made for larger N there.

## 2. What "movers" are

At low energy the critical chain is gapless and its excitations are waves with E = v|k|, travelling right (right-movers) or left (left-movers). In the field theory the spin current splits into J_R(x − vt) and J_L(x + vt), each with its own SU(2).

- **Ring:** the two are independent (a wave goes around and returns still moving the same way). The symmetry is SU(2)_L × SU(2)_R, states are labelled (spin_L, spin_R), and the lowest nontrivial state (½, ½) contains 2 × 2 = 4 states = 1 singlet + 3 triplet, exactly degenerate at the critical point.
- **Open chain:** a wave that hits an end reflects, turning a right-mover into a left-mover with the same spin. The boundary forces J_L = J_R. Unfolding the chain gives one chiral field on a circle of length 2L, with a single SU(2) and a single tower.

## 3. What towers are

A tower is the ladder of evenly spaced levels obtained from one starting state by repeatedly adding quanta, like 0, ħω, 2ħω, … for a harmonic oscillator. In the chain the spacing is πv/L (open) or 2πv/L (ring), with v the wave speed.

Vacuum tower of SU(2)₁ (open chain):

| Level | States | Spin content |
|--|--|--|
| 0 | \|0⟩ | singlet |
| 1 | J₋₁ᵃ\|0⟩ | one triplet |
| 2 | J₋₂ᵃ\|0⟩, J₋₁ᵃJ₋₁ᵇ\|0⟩ | a triplet and a singlet |

So the first excitation is a triplet at 1×(πv/L) and the first excited singlet is at 2×(πv/L): ratio 2.

On a ring there are several towers (left and right each build one; the spectrum is their product). The (½, ½) state sits at the bottom of its own tower, degenerate with the singlet there; the marginal operator splits it, which is why levels can cross.

(In 2D Néel order the "Anderson tower" is a different ladder: E(S) − E₀ ∝ S(S+1)/N, the spectrum of a rigid rotor, spaced ~1/N rather than ~1/L.)

## 4. What Hamiltonian is time-evolved, and what is the initial object

The filter evolves under the **target Hamiltonian itself**; it is not adiabatic evolution and there is no separate initial Hamiltonian.

- H = Σ J1 S_i·S_{i+1} + J2 S_i·S_{i+2}, J1 = 1, open or periodic.
- In the circuit it is shifted and rescaled, H_s = (H − E0)/W, with E0 the ground-state energy and W the certified bandwidth (Weyl bound), so the spectrum lies in [0, 1] with the ground state at 0.
- It is written as Pauli terms (XX, YY, ZZ per bond), evolved controlled on an ancilla: e^{−i t_i H_s⊗Z}. Sequence: H on the ancilla, phase rotation by φ_i, controlled evolution, H on the ancilla. Reading the ancilla as 0 applies cos(H_s t_i + φ_i) to the system. The product over pulses is F(E) = ∏ cos(E t_i + φ_i).
- The Trotter formula splits the controlled evolution into the Pauli terms (first and second order are compared in `results/RESULTS.md` §10 and `WRITEUP.md` §4).

The initial object is a **state**: a shallow-circuit trial state (MPS-to-circuit compilation of the DMRG ground state of the same H) with overlap γ with the true ground state.

Variants: the S² study evolves under the total-spin operator S² with its own pulses; the periodic benchmark uses the same H with periodic bonds.

## 5. How the expectation values differ

Write the trial state as |ψ⟩ = Σ cₙ|n⟩ in eigenstates of H.

**Plain time evolution** e^{−iHt} multiplies each cₙ by a phase. The weights |cₙ|² are unchanged, so ⟨H⟩ = Σ|cₙ|²Eₙ is conserved and observables only oscillate. The ground-state weight γ = |c₀|² never grows.

**The filter** applies cos(H_s t + φ) and post-selects the ancilla on 0. That operation is not unitary on the system: cₙ → cₙ·F(Eₙ), then the state is renormalized. The weights become |cₙ|²F(Eₙ)² / Σ|cₘ|²F(Eₘ)². The design makes |F| ≤ η·|F(0)| on [Δ, 1], so the ground-state weight rises and everything above the gap is suppressed.

| Quantity | Trial state | After the filter |
|--|--|--|
| Ground-state weight | γ = \|c₀\|² | ≥ 1 − ε (certified) |
| Energy ⟨H⟩ | Σ\|cₙ\|²Eₙ, above E0 | → E0 |
| Observables (e.g. SzSz) | trial-state values | → ground-state values |

Example (ours, N = 12, J2 = 0.4, ε = 10⁻²): trial overlap γ = 0.949; after the filter 1 − F ≈ 4×10⁻⁴ in the ED simulation (certified bound 10⁻²).

**Cost of post-selection.** The filter works only on the ancilla-0 branch, with probability p = Σ|cₙ|²F(Eₙ)² (about 0.56 in the same case); repeating until success costs about 1/p runs.

So plain evolution preserves expectation values, and the filter works because post-selection makes the evolution non-unitary and reweights the spectrum toward E0.
