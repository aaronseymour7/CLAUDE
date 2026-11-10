# ORNL: filter costs for J1-J2 ground states

Code and data for the manuscript draft "Filter costs for J1-J2 ground states". We benchmark the fault-tolerant
hardware costs (CNOT and non-Clifford Rz) of ground-state preparation with two cosine filters, the
Stetcu-Baroni-Carlson (SBC) projection algorithm with a dyadic scheduler and QETU, each with and without
amplitude amplification. The testbed is the open 1D J1-J2 chain up to N = 16, in an ideal noiseless simulation
with second-order Trotter evolution.

```
pip install -r requirements.txt
pip install --ignore-requires-python git+https://github.com/qiskit-community/mps-to-circuit
python -m pytest tests -q
```

## What is in the paper and where it comes from

| Paper section | File | Output |
|--|--|--|
| 2. Filters (worked example of SBC, AA, QETU, and the Qiskit circuits that set the cost model) | `filter_walkthrough.ipynb` | printed tables |
| 5. Filter-only scaling: SBC, SBC + AA, QETU, QETU + AA vs N, η = Δ/W and γ; trial population before the filter | `scaling_study.ipynb` | `results/scaling_study.csv`, `results/scaling_study_population.csv` |
| 5. Gap-vs-overlap maps and the range dependence of the η exponent | `run_maps.py`, `analyze.py` | `results/map_N*.json`, `fits.csv`, `all_points.csv`, `figs/` |
| Total spin: current knob vs + S² stage (Stetcu et al. step 1) vs singlet knob, same Hamiltonians | `run_maps.py --config configs/su2_*.json`, `compare_su2.py` | `results/su2_*/`, `results/su2_compare.csv`, `su2_fits.csv`, `figs/su2_cuts_N*.png` |
| 3. Gap knob: dimerisation vs J2 | `plot_gap_vs_delta.py` | `figs/gap_vs_delta.png` |
| 6. Full workflow in Qiskit Aer: DMRG spectral estimates, default mps-to-circuit trial, SBC S^z projection, filter circuits | `workflow_study.ipynb` | `results/workflow_study.csv`, `workflow_prep.csv`, `workflow_study_seed_stats.csv`, `workflow_qetu_check.csv` (5 seeded DMRG draws per N) |
| All of it in one place (standalone, no `ornl` import): both Hamiltonians, ED + DMRG, gap vs overlap study, SBC + S^z projection, AA and QETU in Aer, summary | `final_filters.ipynb` | printed tables |

The notebooks run top to bottom. `scaling_study.ipynb` takes about 25 minutes for N = 6..16;
`workflow_study.ipynb` takes about 20 minutes (5 DMRG seeds).

## The `ornl` package

| File | What it does |
|--|--|
| `ornl/config.py` | `StudyConfig`: every setting for the maps, plus `--config` / `--set` parsing |
| `ornl/model.py` | `Chain`: dimerised J1-J2 Hamiltonian in the S^z = 0 sector, spectrum, the two knobs |
| `ornl/trotter.py` | Strang steps, the exact eigenphase form of a Trotterised pulse, the dyadic filter |
| `ornl/costs.py` | CNOT per attempt and per reflection, RUS and AA cost per success |
| `ornl/search.py` | Measured search for the cheapest (K, r, n) per protocol and target |
| `ornl/symmetry.py` | S² in the sector, singlet gap, the SBC S² projection stage (Trotterised, with its CNOT count) |
| `ornl/prep.py` | S^z-conserving bond-2 truncation and one-layer circuit (the `trial="mps1"` option of the maps) |
| `tests/test_engine.py` | Eigenphase engine vs direct gate-by-gate simulation and `expm`, AA toy model, cost edge cases |

## Model and knobs

H = J1 Σ (1 + (−1)^i δ) S_i·S_{i+1} + J2 Σ S_i·S_{i+2}, open chain, J1 = 1, J2 = 0.2411, S^z = 0 sector.

- Gap: dimerisation δ keeps SU(2) and sets the scaled gap η = Δ/W, W = Emax − E0.
- Overlap: exp(−iε Σ(−1)^i Z_i) on the exact ground state; ε is solved for each target γ. It commutes with
  S^z, so the trial stays in the sector.

## Cost model

Each ancilla-controlled bond is a Clifford Bell change plus a 4-CNOT parity network: 6 CNOT, 3 Rz. Strang steps
merge their end bonds, so m steps take (2B − 2)m + 1 blocks, B = 2N − 3. Reflections for AA are counted from a
Qiskit transpile of `mcx` to {CX, Rz, SX, X}. RUS costs C_attempt / p; AA costs
((2k + 1) C_A + k C_R) / sin²((2k + 1)θ), sin θ = √p, and is used only where it beats RUS. A filter counts only if
it reaches 1 − F ≤ ε at n Strang steps per τ and also at the next finer n.

## Speed

A Trotterised pulse is never time-stepped. One Strang step S is complex symmetric and unitary, so
S = O e^{−iφ} O^T with O real orthogonal, and a pulse of m steps is O cos(mφ) O^T. Every (K, r) at a given n
then costs one vector product. At N = 16 the S^z = 0 sector has 12870 states and one decomposition takes about
2 minutes.

## Total spin (S²) options

Stetcu et al. first project the trial onto the target total spin and only then filter in energy. The default study
skips that step, and its overlap knob (staggered Rz) puts 90-95% of the missing weight in S > 0, so the energy
filter has to resolve the triplet gap. Three `StudyConfig` fields change this:

- `knob = "singlet_bond"`: the trial is one brick layer of SU(2)-invariant gates on the base,
  e^{-iε SWAP} on even bonds and e^{+iε SWAP} on odd bonds (3 CNOT each, counted per attempt). The trial stays a
  singlet, all missing weight sits on excited singlets, and γ reaches ~0 at ε = π/2.
- `s2_stage = True`: the SBC S² projection before the energy filter: pulses cos(t S²), t = π/4, π/8, ...
  (two pulses up to N = 8 and three up to N = 16 remove every J > 0), each Trotterised over all N(N-1)/2 pairs with
  the same 6-CNOT controlled block. Its Strang step count is searched (`s2_n_grid`) together with the energy filter,
  and p, 1 - F and CNOTs are those of the whole two-stage circuit.
- `gap_kind = "singlet"`: the `gaps` grid targets the first excited singlet. The dimerisation knob has a singlet
  dip at small negative δ (two dimer patterns cross), so only the branch to its right is used; at J2 = 0.2411 the
  reachable singlet gap is about 0.31-3.5 J1 (N = 14) or 0.5-3.6 J1 (N = 8).

The energy filter always resolves the gap the trial actually sees: the singlet gap with the singlet knob or the
S² stage, the triplet gap otherwise (`ch.gap`, `ch.gap_s`; both gaps are stored in the maps).
