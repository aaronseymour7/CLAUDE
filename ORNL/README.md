# QPF: the Stetcu-Baroni-Carlson projection filter

The base engine, built to the paper's definitions (main text and Supplemental Material). Our own Hamiltonians and
experiments will be added on top of it later. `SPEC.md` maps every rule of the paper to code, flags the few
choices the paper leaves open, and compares our numbers with the paper's 4 × 4 Heisenberg example.

```
pip install -r requirements.txt
pip install --ignore-requires-python git+https://github.com/qiskit-community/mps-to-circuit   # j1j2 only
python -m pytest tests -q                 # 9 tests, ~15 s
python examples/heisenberg_4x4.py         # paper's example, ~20 s -> figs/, results/heisenberg_4x4.json
python j1j2/run_pipeline.py               # J1-J2 chain, N = 8..16 -> figs/j1j2_pipeline.png, results/j1j2_pipeline.json
python j1j2/trotter_costs.py              # same pipeline with 1st/2nd-order Trotter: CNOT and non-Clifford Rz counts -> results/j1j2_trotter.{json,md}
```

| File | Content |
|--|--|
| `qpf/core.py` | the ancilla step exp(−i(tO + δ) ⊗ Y), post-selected runs, filter function, success probability, scaling of H |
| `qpf/symmetry.py` | quantum-number projection: zero phases, halving times from t1 (J² ladder, shifted targets) |
| `qpf/design.py` | energy projection: port of the authors' design notebook (`reference/`), optimised (t_i) and (t_i, δ_i); constant and Gaussian schedules |
| `reference/optimize_times_phases.ipynb` | the authors' design notebook (shared by A. Baroni), the reference for `design.optimize` |
| `qpf/models.py` | the paper's test system: Heisenberg model on a periodic square lattice, J_z sectors, J², Néel state |
| `j1j2/inputs.py` | J1-J2 chain: operators, DMRG spectral estimates (E0, Emax, first excited singlet), mps-to-circuit trial, ED reference |
| `j1j2/run_pipeline.py` | the paper's pipeline on the open J1-J2 chain (J² projection, designed energy filter iterated to convergence) |
| `j1j2/trotter.py` | QPF step as a Qiskit circuit (ancilla Y via S·H, LieTrotter / SuzukiTrotter(2) PauliEvolutionGate), level-3 transpile counts, and an exact simulator of the same Trotter product |
| `j1j2/trotter_costs.py` | Trotter step sizes per stage for 1−F ≤ 1e-2, 1e-3 and the CNOT / non-Clifford Rz cost per attempt and per success |
| `examples/heisenberg_4x4.py` | J² projection then energy projection, as in the paper's Figs. 1-2 |
