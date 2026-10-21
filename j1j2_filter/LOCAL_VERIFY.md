# Running and verifying this locally (macOS, Apple silicon)

I produced everything on a Linux x86-64 container with Python 3.12; I have **not** run it on a Mac. Expect the setup below to work, but if a package fails to install on arm64, that is the place to look first.

Tested: the pinned install into a fresh Python 3.12 venv took ~1 min here (Linux), and checks 1–4 below then gave exactly the numbers listed. Not tested on macOS.

## Setup

```bash
brew install python@3.12 git                       # skip what you already have
git clone https://github.com/aaronseymour7/CLAUDE.git
cd CLAUDE && git checkout claude/elegant-bell-ypgph4
python3.12 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip
pip install -r j1j2_filter/requirements.txt         # pinned versions used for the results
cd j1j2_filter
export OMP_NUM_THREADS=1                            # avoids thread oversubscription; the code is mostly single-core
```

All commands below are run from `CLAUDE/j1j2_filter`. The stored results are in `results/data/*.json`, so the checks below do not need DMRG unless stated.

## Checks, fastest first

| # | Command | Time (my machine) | What to expect |
|--|--|--|--|
| 1 | `python scripts/audit_results.py` | ~5 s | `66 data files, 93 designs, 3830 checks` and `ALL CHECKS PASSED`. Independently re-derives every stored bound, step count, gate count and time from the stored inputs (clean re-implementation of the formulas in `CERTIFICATION.tex`; no call into the pipeline's bound code). `--alpha` also recomputes the commutator constant α for every point (~2 min). |
| 2 | `python trotter2/validate_two_term.py` | ~2 s | `max error/bound over 1800 random cases: 0.9981`. Random-matrix check of the two-term Strang bound (Childs et al., Prop. 16 / Eq. 150). |
| 3 | `python trotter2/first_order_full.py 6 0.4` | ~20 s | Full transpile of the first-order certified circuit: `rzz: 21276` (= 2 × 27 Pauli terms × 394 steps), level 3. Writes `trotter2/first_6_0.4.json`. |
| 4 | `python trotter2/second_order_circuit.py 6 0.4` | ~20 s–2 min (18 s in a clean venv here) | Rebuilds the trial state (runs DMRG + ED), builds the second-order circuit, transpiles it, simulates it. Expect `n1=394, n2=20`, second-order level-3 `rzz: 1824`, `bound_holds: true`, circuit-vs-simulation deviation ~1e-16, oracle steps 6 (1st order) and 4 (2nd order). |
| 5 | `python trotter2/dephasing.py 6 0.4` then `python trotter2/mixing_theory.py 6 0.4` | ~3–4 min | `mixing_theory.py` prints `predicts trace distance 9.683e-04, measured 9.720e-04 (ratio 0.996)` and `max|diag X1| ~ 4e-18`. The first-order perturbation theory of the Trotter error against the measured error. |
| 6 | `python trotter2/scaling_check.py 6 0.4` | ~1–2 min | First-order error exponent −1.00 in n; second-order −2.3 over the whole range (−2.0 asymptotically; the table printed shows factors of 4 per doubling at n = 20, 40, 80). |
| 7 | `python trotter2/compare.py` needs `NONLY=<N>`; e.g. `NONLY=6 python trotter2/compare.py` | ~1 min | Certified step counts first vs second order from the stored designs. |
| 8 | `python scripts/summarize_low_levels.py` | ~1 s | Open chains N = 8, 12, 16, 20: triplet is the lowest excitation at every grid point, smallest singlet–triplet difference 0.0112, 0.0059, 0.0020, 0.0009; periodic chains: crossing at J2 = 0.2463, 0.2435, 0.2425, 0.2420. (`WRITEUP.md` §7) |
| 9 | `python pbc_corr/summarize_filter_states.py` | ~1 s | 19 periodic cases, ideal filter: max relative correlator error 0.044, entropy 0.057, energy 1.9e-4; 9 cases without a feasible design. (`WRITEUP.md` §7) |
| 10 | `python s2proj/analyze.py` | ~1 s | Certified filter RZZ without projection (A), with an S² stage (B) and with a symmetric trial (SYM): `B/A over 29 feasible cases: min 1.30, median 3.82, max 10.07`; `SYM/A ... min 0.05, median 0.24, max 0.99`. (`WRITEUP.md` §6) |
| 11 | `python s2proj/check_trotter_spin.py` | ~10 s | Weight outside J = 0 after Trotterised pulses of both orders: magnitude ≤ 7e-17 (N = 6, 8). |
| 12 | `python trotter2/check_x1_symmetry.py` | ~10 s | The first-order Trotter correction X₁ commutes with S² and S^z: 6e-17 / 3e-16 and 0 (N = 6, 8). |

Cases are `N J2` with N ∈ {6, 8, 10, 12} and J2 ∈ {0, 0.2411, 0.4}. Stay at N ≤ 10 on a laptop: N = 12 needs several GB of memory (dense 4096 × 4096 matrices) and hours for the first-order full transpile. The Mac Air has no fan, so long runs will throttle.

## Re-deriving a stored point from scratch

```bash
RESULTS_DATA=/tmp/check ONLY_N=6 python scripts/run_results.py     # reruns DMRG + ED + design + search for N=6 (all J2); resumable
```
then compare `/tmp/check/dmrg_N6_J2*_eps0.01.json` with `results/data/` (certified n, designs, counts). The design search involves a random seed fixed in the settings, so the pulse times/phases should reproduce but I have not verified bit-for-bit reproduction on another machine.

## Rebuilding the documents

```bash
python scripts/build_results.py          # results/RESULTS.md (+ csv tables, figures; HTML if pandoc is installed)
python trotter2/analyze_dephasing.py     # trotter2/dephasing_table.md and dephasing.png
```

## Where the claims live

| Claim | File to read | Check |
|--|--|--|
| Certificate and its proof | `CERTIFICATION.tex`, `paper.tex` App. A–B | #1 |
| RZZ gate counts, both orders, certified vs oracle | `results/RESULTS.md` §4, §9, §10; `trotter2/circuit_table.md` | #1, #3, #4 |
| Second-order bound (Childs et al. Prop. 16, Eq. 152) | `WRITEUP.md` §4.1; `trotter2/alpha2.py` | #2, `trotter2/validate_alpha2.py` (slow: ~10+ min) |
| Dephasing / mixing analysis | `trotter2/DEPHASING.md` | #5, #6 |
| Circuit and operators (W_i, ancilla, Rz) | `core/circuits.py` (`apply_filter_pulse`, `build_early_abort_circuit`), `EXPLAINER.md` | `python -c "from core.circuits import verify_pulse_convention as v; v()"` |

## Known limits of the verification

* `WRITEUP.md` §8 lists what is not shown (no hardware or noisy simulation; second-order and mixing results only for N ≤ 12 / ≤ 10; the oracle only where the exact ground state is available). The 2D crossover code in `crossover2d/` is a cost model, is not used in `WRITEUP.md`, and should not be cited as a result.
* N > 14 results use a DMRG-based trial overlap that is not certified, and are counted, not simulated.
* The first-order full-transpile runs exist for 7 of the 12 cases (`trotter2/first_*.json`); the second-order circuits and oracle values exist for all 12, and the first-order RZZ totals do not depend on the missing transpiles (per-step × n matched the full transpile in all 7).
