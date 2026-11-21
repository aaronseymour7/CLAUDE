# QPF energy-pulse design: authors' v4 vs a minimax option

**This is an option reported next to the published QPF. It does not replace it.** The circuit is the published one:
the J² halving ladder, then Π_i cos(t_i H + δ_i) with one ancilla and one measurement per pulse. Only the pulse times
(and optionally the phases) change. The manuscript's published-QPF numbers are untouched.

Code: `aaronseymour7/CLAUDE`, branch `claude/project-thread-a8kce9`. The design is in `ORNL/qpf/design_plus.py`, the
pipeline driver in `ORNL/j1j2/qpf_design_compare.py`, and the analysis, fits and figures in
`ORNL/j1j2/qpf_design_compare_{analyze,fits,fig}.py` and `qpf_design_step_fig.py`. Raw per-point JSON is in `raw/` here.

## The authors' v4, and what limits it
- **Objective:** |1 − 1/√(Σ_n Π_i cos²(t_i E_n))| on E = {0} ∪ linspace(Δ, 1, 999). This is a monotone function of
  the **mean** of f(E)² over [Δ, 1], i.e. a trial spread evenly over the band.
- **Constraints:** Σ t_i = T = π/Δ, 0 ≤ t_i ≤ T/3. The optimizer is SLSQP from a halving start, with m = 4–12 scanned.
  v4 passes the phases to SLSQP but ignores them.
- **The optimizer is not the bottleneck.** Multistart SLSQP on the same objective and constraints lowers the v4 cost by
  at most about 1%.
- **The average hides the gap edge.** The mean of f² is about 1e-3, but f² at E = Δ is 0.03–0.08. On the real chain,
  the lowest excited states sit right there, so one block suppresses them by only 12–30x. The pipeline then repeats the
  block: the median is R = 2.
- **The cap binds at T = π/Δ.** A zero of f at E = Δ needs one pulse of length π/(2Δ) = T/2 > T/3.
- **The design is fragile to round-off.** A 1e-14 change in the scaled gap can flip the chosen m between a design with a
  zero-length pulse and one without. A zero-length pulse still costs a Trotter step, an ancilla and a measurement. It can
  also land in a different local minimum (N = 12, singlet gap 0.55: v4 cost 0.289 vs the manuscript's 0.274). Because of
  this, 6% of the manuscript points (131 of 2224) differ by more than 1% when the same code is rerun, mostly by
  5–30% and at worst 0.68x to 2.0x (median ratio 1.00).

## The option
- **Objective (minimax):** C(t) = max over E ∈ [Δ, 1] of Π_i cos²(t_i E), the worst leaked level. With phases, it is
  max f(E)²/f(0)² subject to f(0)² = Π cos²δ_i ≥ p₀.
  - For any trial with ground weight γ, one block gives 1 − F ≤ (1 − γ)/γ · C.
  - Minimax is therefore optimal against the worst trial with that γ. It needs no trial or spectrum beyond Δ and W.
- **Time knob:** T = c·π/Δ with c ∈ {1, 1.25, 1.5, 2}. The pipeline chooses c per point, as it already chooses R, the
  Trotter steps and RUS/AA. The authors' T/3 cap is kept. One extra schedule drops the cap at c = 1.
- **Optimizer:** a soft-max (log-sum-exp of log f², β = 2 → 2048) with analytic gradients. It uses SLSQP with 24 starts
  per m and m = 4–12, and drops pulses shorter than 1e-3·T.
- **Phases:** p₀ = 0.9 with c ∈ {1, 1.25, 1.5}, and p₀ = 0.7 with c ∈ {1, 1.25}. On the design band, phases give 4–8x
  more suppression at T = π/Δ for p₀ = 0.9. On the real chain they add almost nothing once c is free (see below).

`fig_filter_step.png` shows f(E)² for v4, minimax at T = π/Δ, and minimax at T = 1.5π/Δ, at scaled gaps 0.1 and 0.3.

## Result on the Section A grid
The grid matches the manuscript exactly: N = 8–14, 6 singlet gaps × 6 overlaps, both knobs, ε = 1e-2 and 1e-3, R ≤ 4,
RUS and AA, the analytic 6-CNOT model, and Qiskit 2.6.2 for the AA reflections. My v4 rerun matches the manuscript's v4
(median ratio 1.00).

CNOT per success, as the ratio of the option to the authors' v4 (median, with 10–90% in brackets):

| route | protocol | ε | minimax, c = 1 only | minimax, c chosen per point | + phases (all options) |
|--|--|--|--|--|--|
| J² → H (published) | RUS | 1e-2 | 1.00 | **0.85** (0.70–1.00) | 0.85 |
| J² → H | RUS | 1e-3 | 1.00 | **0.81** (0.66–0.96) | 0.81 |
| J² → H | AA | 1e-2 | 1.00 | **0.84** (0.68–1.00) | 0.82 |
| J² → H | AA | 1e-3 | 1.00 | **0.79** (0.64–0.94) | 0.77 |
| H only | RUS | 1e-2 | 1.00 | **0.68** (0.46–0.95) | 0.62 |
| H only | RUS | 1e-3 | 1.00 | **0.56** (0.44–0.80) | 0.56 |

- The option is cheaper at 90–98% of points. Against the manuscript's "best of v3 and v4", the J² → H ratio is 0.83–0.89.
- **Mechanism:** a slightly longer, sharper single block replaces repeated blocks. For J² → H RUS, the median R drops
  from 2 to 1, the energy-stage time from 44 to 32 (scaled units), Trotter steps from 14 to 10, and CNOT per attempt from
  6.3k to 5.0k. Success p is unchanged.
- Minimax at the authors' T = π/Δ alone gives nothing (median 1.00): with the T/3 cap it cannot null the gap edge.
- The gain is flat in N for H only (about 0.6). For J² → H it shrinks slightly at N = 14 (0.80 at N = 8 to 0.86 at N = 14),
  because the J² stage is a larger share there.

## Impact on the manuscript (flag, nothing changed)
- **Section A exponents barely move** (`design_compare_fits.csv`). For J² → H, staggered knob, RUS, ε = 1e-2, the fit
  goes from a = 0.38, b = 1.15, α = 2.28 (manuscript "best", reproduced) to a = 0.37, b = 1.10, α = 2.33. Other cells move
  by ≤ 0.05 in a and b, and by up to 0.2 in α. The qualitative story holds, with a prefactor about 15% lower.
- **Section B (QPF vs QETU):** published QPF is already cheaper on most grid points. The option would widen that by about
  15–20% for QPF only. If it is shown at all, it should be a separate, clearly labelled row ("QPF, minimax design (ours)"),
  not the published baseline.
- **Reproducibility note for the methods section:** the v4 design is sensitive to round-off (see above). Pinning the design
  times, as the raw JSON already does, is what makes the manuscript numbers reproducible.

## Files
- `design_compare_points.csv`: one row per point × schedule (analytic model). Columns include R, Trotter steps, times, p
  and CNOT per success.
- `design_compare_wide.csv`: CNOT per success per schedule side by side, plus the manuscript's `ms v3`, `ms v4` and
  `ms best`.
- `design_compare_summary.csv`: ratio statistics per route, protocol, ε, option and baseline.
- `design_compare_fits.csv`: Section A fits for v4, the manuscript's best, and the options.
- `design_compare_designs.csv`: every design (times, phases, m, design cost).
- `fig_filter_step.png`, `fig_cost_ratio.png`.
