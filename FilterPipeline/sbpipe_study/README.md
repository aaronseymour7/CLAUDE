# J1–J2 application study and paper for `sbpipe`

DMRG energy estimates and matrix-product trial states for the open spin-½ J1–J2 chain are fed to the
generic pipeline in `../sbpipe/`; the study validates the certificate by simulation (N ≤ 14) and counts
resources up to N = 48.

| file | purpose |
|--|--|
| `example_j1j2.py` | one case end to end: `python sbpipe_study/example_j1j2.py 12 0.2411 2 1e-2` |
| `run_study.py` | all experiments; `python sbpipe_study/run_study.py all` (or `inputs`, `validate`, `scaling`, `sensitivity`, `ansatz`, `oracle`, comma-separated) |
| `make_figures.py` | figures, tables and the numeric macros the paper uses (`paper/generated/`) |
| `results/*.json` | raw results; `results/cache/*.json` are the DMRG inputs (the `.npz` statevectors and design cache are not committed) |
| `paper/paper.tex`, `paper/paper.pdf` | the paper (revtex4-2; `pdflatex paper.tex` twice) |

Run from the repository root. Requirements: numpy, scipy, matplotlib, quimb; qiskit and pytest only for the
tests. `JOBS` sets the number of worker processes (default 4). Approximate run times on 4 cores: inputs
5 min, validate 15–35 min, scaling 30 min (the optimized designs at N = 48 dominate; they are cached in
`results/cache/designs/`), the rest a few minutes each.

## Headline results (from `paper/generated/macros.tex`)

* The certificate held in all 96 simulated circuits (N = 8–14, three couplings, ε = 10⁻², 10⁻³, three
  schedules, both Trotter orders); median looseness 14× in infidelity, mostly from the Trotter bound
  (median 185×). Certified step counts are 9–39× the smallest step count that actually meets ε.
* Second-order Trotter is a median 92× cheaper than first order; the optimized filter a median 2.9×
  cheaper than the closed-form dyadic schedule, but the optimizer fails for small overlaps (14 of 144
  scaling plans, all N ≥ 32), where the dyadic schedule still works.
* Overestimating the gap by ≥ 50% breaks the target; underestimating it by 50% costs 2.9×.
* N = 48, J2 = 0.2411, ε = 10⁻²: 9.7×10⁶ CNOTs per success from a χ = 4 trial. A χ = 8 MPS circuit alone
  already meets ε = 10⁻² (4,500 CNOTs), so in 1D at this target direct MPS preparation wins.
