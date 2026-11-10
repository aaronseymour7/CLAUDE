# QPF as defined by Stetcu, Baroni and Carlson

Source: I. Stetcu, A. Baroni, J. Carlson, "Projection algorithm for state preparation on quantum computers",
Phys. Rev. C 108, L031306 (2023), [arXiv:2211.10545](https://arxiv.org/abs/2211.10545) (main text and Supplemental
Material). Every rule below is the paper's (or, for the design step, the authors' notebook) unless it is marked **choice**. Choices only cover points the paper
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

Reference: the authors' design notebook (A. Baroni, `reference/optimize_times_phases.ipynb`). `design.optimize`
is a line-by-line port and reproduces the notebook's example output exactly (cost 1.319e−01,
fdiff 2.498529e−03; see the tests).

| Notebook | Code |
|--|--|
| Design spectrum: E = [0, linspace(gap, 1, N − 1)], N = 1000, one ground level, uniform trial | `baroni_spectrum(gap, 1000)` |
| Objective v3: C = \|1 − ∏cos δ_i / sqrt(Σ_n ∏cos²(t_i E_n + δ_i))\| (1 − ground amplitude of the normalised filtered uniform state) | `baroni_cost` |
| Σ t_i = T = π/gap (equality); 0 ≤ t_i ≤ T/3; −π/2 ≤ δ_i ≤ π/2 | `optimize(gap, m, T)` |
| v3: ∏cos²δ_i ≥ 0.9 (keeps the success within 10% of the zero-phase value); v3b: ∏cos²(E_0 t_i + δ_i) ≥ 0.9; v4: no phases | `method="v3" / "v3b" / "v4"`, `p_min=0.9` |
| Start: halving times scaled to T, zero phases; SLSQP, maxiter 5000, ftol 1e−12, tol 1e−13 | same |
| Evaluation: trial with amplitude γ on level 0, the rest uniform; filtered, normalised; fdiff = Σ(gs − f)² | `evaluate`, `notebook_trial` |

From the paper (not in the notebook): steps applied shortest time first (used in the example's trajectories);
comparison schedules **constant** (m equal times, m chosen for the lowest energy at fixed T) and **gaussian**
(7 normal draws scaled to Σ\|t_i\| = T, averaged; **choice**: N(0, 1), 20 draws).

Where the notebook differs from the paper's text: no random displacement of the levels and no αe in the design
(the overlap only enters the evaluation), one ground level instead of a fraction αe of 1000, and the 0.9
constraint on f(0)², which the text does not mention. That constraint is why the phases are small. The text's
reading is kept as `optimize_supplement` for comparison only.

## Reproduction of the paper's Heisenberg example (`examples/heisenberg_4x4.py`)

4 × 4 periodic lattice, H = Σ σ_i·σ_j, h = 0, J_z = 0 sector (12 870 states), Néel initial state.

| Quantity | Paper | This code |
|--|--|--|
| Néel scaled energy | ≈ 0.17 | 0.168 |
| Gap to the lowest relevant state, before / after J² projection | ≈ 0.03 / 0.15 | 0.030 / 0.155 |
| J² projection steps | 3 | 3 (π/4, π/8, π/16) |
| J² projection success | ≈ 11% | 11.1% (exactly 1/9) |
| Energy after J² projection | ≈ 0.06 | 0.064 |
| Ordering after energy projection at T = π/Δ | constant worst, QPF(t_i, δ_i) best | constant 1.7e−3, Gaussian 1.0e−3, QPF(t_i) 3.6e−4, QPF(t_i, δ_i) 1.4e−5 |
| Effect of phases on success | slightly lower | 0.0828 → 0.0743 |
| Phases | small, positive, larger for longer times | 0.06 … 0.23 rad, increasing with t_i |

The paper only gives the energies after the energy projection as a plot, so this comparison is by ordering and
order of magnitude.

## Known ambiguities in the paper

- The supplement says the Fig. 1 filter assumes a gap of 0.2, while the main text says 0.15. We use 0.15.
- The objective and constraints are only fully defined by the authors' notebook (above), not by the text.
