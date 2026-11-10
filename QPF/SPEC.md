# QPF as defined by Stetcu, Baroni and Carlson

Source: I. Stetcu, A. Baroni, J. Carlson, "Projection algorithm for state preparation on quantum computers",
Phys. Rev. C 108, L031306 (2023), [arXiv:2211.10545](https://arxiv.org/abs/2211.10545) (main text and Supplemental
Material). Every rule below is the paper's unless it is marked **choice**. Choices only cover points the paper
leaves open.

## The step (`qpf/core.py`)

| Paper | Code |
|--|--|
| Ancilla a in \|0⟩, unitary exp(−i(t_i O + δ_i) ⊗ Y_a), measure a (Eq. 1–2) | `ancilla_step(O, psi, t, delta)` returns both branches; checked against the dense unitary in tests |
| Keep reading 0: cos(t_i O + δ_i)\|ψ⟩, probability ⟨ψ\|cos²(t_i O + δ_i)\|ψ⟩ | `run(O, psi, times, phases)` |
| All steps succeed: N ∏_i cos(t_i O + δ_i)\|ψ(0)⟩, success = product of step probabilities | `Run.p_total`, `Run.state` |
| Filter function f = ∏_i cos(t_i o_α + δ_i) on the eigenvalues o_α | `filter_function(o, times, phases)` |
| Converged success P = f(0)² \|⟨Ψ0\|ψ(0)⟩\|², f(0) = ∏_i cos δ_i | `success_probability_limit` |
| H shifted and scaled to [0, 1] with the target at 0 | `scale(H, E0, Emax)` |
| **choice**: exact evolution (no Trotter), through `expm_multiply` | `Operator`; `Diagonal` when O is given by its eigenvalues |

## Quantum-number projection (`qpf/symmetry.py`)

| Paper | Code |
|--|--|
| All δ_i = 0 | `run(..., phases=None)` |
| t1 sets the resolution, t_{i+1} = t_i / 2; total time → 2 t1 | `halving_times`, `ladder_for(values, t1)` (shortest ladder that zeroes every unwanted value) |
| J² with target J = 0: t1 = π/4 removes J = 1, 2, 5, 6, ...; three steps for J ≤ 8 | `j2_ladder(0, 8)` → π/4, π/8, π/16 |
| Other targets by shifting the operator (J² − 12 for J = 3) | `j2_values(J_target, J_max)`, `j2_ladder(J_target, ...)` |

## Energy projection design (`qpf/design.py`)

| Paper | Code |
|--|--|
| Inputs: target near 0, gap estimate Δe, overlap estimate αe | `optimize(gap, alpha, m, T, phases)` |
| Supplement: 1000 levels, equal amplitudes, fraction αe at zero energy (50 of 1000), the rest on a grid from the gap to 1, "randomly displaced within a grid spacing" | `pseudo_spectrum(gap, alpha, 1000)`; **choice**: uniform displacement inside each grid cell |
| Supplement: minimise the mean squared difference between the filtered state and the target (1 at zero energy, 0 from the gap to 1) | `mse` = mean_k (f(E_k) − target_k)² |
| Total time limited to π/Δe (Fig. 1, Heisenberg) | `T = π/gap`; **choice**: the total is used in full, Σ t_i = T, t_i ≥ 0 |
| Six (t_i, δ_i) for the Heisenberg case | `m = 6` |
| Steps applied shortest time first | `optimize` returns times sorted ascending |
| **choice**: L-BFGS-B with 30 random restarts; phases in (−π/2, π/2) | `optimize(..., restarts=30)` |
| Constant times, "similar to Ge et al.", number of steps chosen for the lowest energy at fixed T | `constant(T, m)`; the example scans m = 1..30 |
| Gaussian times (Rodeo), 7 draws from a normal distribution scaled to T, averaged over draws; negative times allowed, total = Σ\|t_i\| | `gaussian(T, 7, seed)`; **choice**: N(0, 1), 20 draws |

With these definitions the optimiser reproduces what the paper reports about its phases: small, positive, larger
for the longer times, about 0 for the shortest one, and f(0) close to 1 (0.96 for Δe = 0.15, αe = 0.05).

## Reproduction of the paper's Heisenberg example (`examples/heisenberg_4x4.py`)

4 × 4 periodic lattice, H = Σ σ_i·σ_j, h = 0, J_z = 0 sector (12 870 states), Néel initial state.

| Quantity | Paper | This code |
|--|--|--|
| Néel scaled energy | ≈ 0.17 | 0.168 |
| Gap to the lowest relevant state, before / after J² projection | ≈ 0.03 / 0.15 | 0.030 / 0.155 |
| J² projection steps | 3 | 3 (π/4, π/8, π/16) |
| J² projection success | ≈ 11% | 11.1% (exactly 1/9) |
| Energy after J² projection | ≈ 0.06 | 0.064 |
| Ordering after energy projection at T = π/Δ | constant worst, QPF(t_i, δ_i) best | constant 1.7e−3, Gaussian 1.0e−3, QPF(t_i) 5.5e−5, QPF(t_i, δ_i) 1.5e−5 |
| Effect of phases on success | slightly lower | 0.0826 → 0.0766 |

The paper only gives the energies after the energy projection as a plot, so this comparison is by ordering and
order of magnitude.

## Known ambiguities in the paper

- The supplement says the Fig. 1 filter assumes a gap of 0.2, while the main text says 0.15. We use 0.15.
- The main text reads the optimisation as "maximize the squared overlap of the final state with the ground
  state", while the supplement specifies the mean-squared-error objective. We use the supplement's version,
  because maximising the normalised overlap alone drives f(0)² down to about 0.1, which contradicts the paper's
  "small non-zero phases" and "typically close to one".
