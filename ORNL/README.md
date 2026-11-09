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
| 3. Gap knob: dimerisation vs J2 | `plot_gap_vs_delta.py` | `figs/gap_vs_delta.png` |
| 6. Full workflow in Qiskit Aer: DMRG spectral estimates, default mps-to-circuit trial, filter circuits | `workflow_study.ipynb` | `results/workflow_study.csv`, `workflow_prep.csv`, `workflow_qetu_check.csv` |

The notebooks run top to bottom. `scaling_study.ipynb` takes about 25 minutes for N = 6..16;
`workflow_study.ipynb` takes about 5 minutes.

## The `ornl` package

| File | What it does |
|--|--|
| `ornl/config.py` | `StudyConfig`: every setting for the maps, plus `--config` / `--set` parsing |
| `ornl/model.py` | `Chain`: dimerised J1-J2 Hamiltonian in the S^z = 0 sector, spectrum, the two knobs |
| `ornl/trotter.py` | Strang steps, the exact eigenphase form of a Trotterised pulse, the dyadic filter |
| `ornl/costs.py` | CNOT per attempt and per reflection, RUS and AA cost per success |
| `ornl/search.py` | Measured search for the cheapest (K, r, n) per protocol and target |
| `ornl/prep.py` | S^z-conserving bond-dimension-2 truncation used for the mps-to-circuit trial |
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
