# Total spin and the filter cost (S² stage, singlet knob)

Why: Stetcu et al. (QPF) first project the trial onto the target total spin J, then filter in energy. Our study
skipped that step, and its staggered-Rz overlap knob puts 90-95% of the missing weight into S > 0, so the energy
filter had to resolve the triplet gap. See README, "Total spin (S²) options", for the code.

All numbers: exact trial base, ε = 1e-2, RUS, CNOT per success including the knob layer and the S² stage.

## 1. Original grid (triplet gap 0.005-1.7 J1), with and without the S² stage

`results/map_N*.json` (A) vs `results/su2_D_triplet_grid_s2stage/` (D = same Hamiltonians and trials + S² stage).

| N | triplet gap | singlet gap | γ = 0.5: D/A | γ = 0.1: D/A | γ = 0.01: D/A |
|--|--|--|--|--|--|
| 12 | 0.005 | 0.53 | 0.006 | 0.002 | 0.003 |
| 12 | 0.030 | 0.38 | 0.029 | 0.052 | 0.054 |
| 12 | 0.115 | 0.38 | 0.32 | 0.33 | 0.53 |
| 12 | 0.443 | 0.79 | 0.90 | 1.37 | 1.05 |
| 12 | 1.70 | 3.31 | 1.66 | 1.39 | 1.09 |

N = 8 and 10 look the same (D/A 0.002-0.007 at the smallest gap, 0.7-2.3 at the largest). The whole small-gap end
of the original maps (δ < 0, where the gap closes) is the edge-spin triplet: one S² stage removes it and the
energy filter then works at a singlet gap of 0.35-0.75 J1. Where the triplet gap is not small, the S² stage
(all N(N-1)/2 pairs) costs more than it saves, up to about 1.7-2.3x.

## 2. Singlet-gap grid (0.3-3.5 J1, right of the singlet dip): three variants on the same Hamiltonians

`results/su2_compare.csv`, `figs/su2_cuts_N*.png`; ratios to A at the same (δ, γ), median (min-max) over the map:

| N | B/A (Rz knob + S² stage) | C/A (singlet knob) |
|--|--|--|
| 8 | 1.03 (0.16-8.4) | 0.53 (0.19-1.14) |
| 10 | 1.11 (0.38-11.6) | 0.53 (0.22-1.14) |
| 12 | 1.28 (0.21-18) | 0.53 (0.18-1.14) |
| 14 | 1.60 (0.43-27) | 0.52 (0.22-1.13) |

Fitted exponents, ln C = c + a ln(W/Δ_singlet) + b ln(1/γ) (RUS, ε = 1e-2, N = 14):
A a = 1.17, b = 1.15; B a = 0.41, b = 1.10; C a = 0.88, b = 1.11.

- B's overhead grows with N (the stage is O(N²) pairs); it pays only where the triplet gap is much smaller than the
  singlet gap.
- C, the singlet knob, is about half the cost of A at every N with a near-linear gap exponent (0.8-1.0). With an SU(2)
  trial no S² stage is needed, so this is the cleaner benchmark of the energy filter itself.
- With the singlet gap as the variable, the reachable range at J2 = 0.2411 is only about one decade (0.3-3.5 J1);
  the dimerisation knob has a singlet dip near δ ≈ -0.05 to -0.15. J2 = 0.5 would give about two decades
  (minimum singlet gap 0.03 J1 at N = 14).
