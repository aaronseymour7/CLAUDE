# ORNL: measured cost of a ground-state filter vs. gap and overlap

A small, standalone study of what a Trotterised cosine filter costs, and how the spectral gap and the
trial-state overlap drive that cost, with and without amplitude amplification. The example is the open 1D
J1-J2 chain. Everything is measured by exact simulation of the circuit; there are no bounds or certificates.
Small cases are cross-checked against real Qiskit circuits in noiseless state-vector simulation.

The engine needs only numpy, scipy and matplotlib. Qiskit is optional and only used for verification.

```
pip install -r requirements.txt
python -m pytest tests -q        # 10 checks (Qiskit ones are skipped if qiskit is missing)
python run_maps.py               # N = 8..16, eps = 1e-2 1e-3  -> results/map_N*.json
python analyze.py                # fits -> results/fits.csv, figures -> figs/
python verify_qiskit.py          # engine vs Qiskit circuits at N = 4, 6, 8 -> results/verify_qiskit.json
```

Open `playground.ipynb` to play with single points interactively (runs in about 30 s).

## Changing variables

Every study setting lives in one place, `ornl/config.py` (`StudyConfig`): couplings, system sizes, gap and
overlap grids, target infidelities, Trotter-step grid, K and r limits, protocols, and the CNOT cost model.
Every script takes the same two options:

```
python run_maps.py --set Ns=[8,10] --set eps=[0.01]          # override any field, values are JSON
python run_maps.py --set J2=0.3 --set out_dir='"results_J2_0.3"'
python run_maps.py --config my_study.json                     # a JSON file with any subset of fields
python analyze.py --config my_study.json                      # reads the same out_dir, writes fig_dir
python verify_qiskit.py --N 6 --gaps 0.3 1.0 --gammas 0.1 0.5 --compile pauli
```

Unknown field names are rejected. `run_maps.py` writes its JSON after every gap value and resumes where it
stopped; it refuses to mix a saved file with a different grid.

From Python, each piece can be used alone:

```python
from ornl import StudyConfig, chain_for_gap, measure_point, costs, trotter
cfg = StudyConfig(J2=0.3)
ch = chain_for_gap(10, gap=0.2, cfg=cfg)                 # gap knob: dimerisation
psi = ch.trial(ch.eps_for_overlap(0.1))                  # overlap knob: staggered Rz
v, p, infid = trotter.dyadic_filter(ch, psi, K=6, r=2, n=2)
c = costs.cnot_per_attempt(10, 6, 2, 2)
print(p, infid, costs.cost(c, p, "rus"), costs.cost(c, p, "aa", costs.reflection_cnots(10, 6, 2)))
rows, eps_knob = measure_point(ch, [0.1, 0.5], [1e-2, 1e-3], cfg)   # cheapest filter per protocol
```

## Files

| File | What it does |
|--|--|
| `ornl/config.py` | `StudyConfig`: every variable, plus `--config` / `--set` parsing |
| `ornl/model.py` | `Chain`: Hamiltonian in the S^z = 0 sector, spectrum, the two knobs, map to full 2^N vectors |
| `ornl/trotter.py` | Strang steps, filter pulses, the exact eigenphase form, the dyadic filter |
| `ornl/costs.py` | CNOT count per attempt and per reflection, RUS and AA cost |
| `ornl/search.py` | Measured search for the cheapest (K, r, n) per protocol and target |
| `ornl/qiskit_check.py` | Real Qiskit circuits: KAK bond blocks, the attempt A, the AA circuit, post-selection |
| `run_maps.py` | Runs the (gap, overlap) grid for each N and eps |
| `analyze.py` | Power-law fits per N and over all N, and the figures |
| `verify_qiskit.py` | Engine vs Qiskit, point by point |
| `playground.ipynb` | Interactive notebook; regenerate with `make_notebook.py` |
| `tests/` | `test_engine.py` (direct simulation, `expm`, Grover toy model), `test_qiskit.py` (circuits) |
| `results/` | One JSON per N, `fits.csv`, `verify_qiskit.json` |
| `figs/` | `fig_knobs`, `fig_maps_eps*`, `fig_cuts_eps*`, `fig_exponents`, `fig_local_exponents` |

## Qiskit verification

`ornl/qiskit_check.py` builds the circuit the engine costs, gate for gate:

- Register: system qubits 0..N−1 (site i = qubit i), then one ancilla per pulse, so one attempt is a single
  unitary A and post-selection on "all ancillas 0" happens at the end. This is the form AA needs.
- Trial: an exact `StatePreparation` of the trial vector (not costed, same as the engine).
- Controlled bond exp(−iφ (XX+YY+ZZ) ⊗ Z_a), two compile modes:
  - `bond` (default): a Bell-basis change on the pair, synthesised by Qiskit's KAK
    `TwoQubitBasisDecomposer(CXGate())` (1 CNOT each side), which makes XX+YY+ZZ diagonal; then a 4-CNOT
    parity network onto the ancilla. 6 CNOTs per block, which is what `ornl.costs` counts.
  - `pauli`: one `PauliEvolutionGate` per bond, left to Qiskit's level-3 transpiler (which re-synthesises
    two-qubit blocks with KAK). This shows what a generic compile gives.
- AA: k rounds of R_good (sign flip on "ancillas all 0"), A†, R_0 (sign flip on |0…0⟩ of every qubit), A,
  with Qiskit's `mcx` for the multi-controlled Z.

Results of `python verify_qiskit.py` (27 points, N = 4, 6, 8, up to 17 qubits, in `results/verify_qiskit.json`):

- Success probability, 1 − F and AA success after k rounds agree with the engine to 1e−12 or better.
- The CNOT count of the built filter equals the formula exactly.
- The `pauli` compile costs about 5% more CNOTs (1262 vs 1200 at N = 6); transpiling the `bond` circuit at
  level 3 removes about 7% (1118).
- Qiskit's `mcx` without helper qubits costs 1.5 to 3 times the 12n CNOTs per reflection that the cost
  model assumes. Constructions that use spare
  qubits get closer to linear; set `mcz_cnots_per_qubit` to see how the AA costs move.

## What is computed

**Model.** H = Σ J1 (1 + (−1)^i δ) S_i·S_{i+1} + J2 Σ S_i·S_{i+2}, open chain, J1 = 1, J2 = 0.2411.

**Gap knob: dimerisation δ.** It keeps SU(2), so the ground state is a singlet and E1 a triplet, both in
S^z = 0. δ is solved for each target gap. The reachable range is from about 0 (δ < 0, where weak end bonds
leave nearly free edge spins) up to about 1.8 J1 (δ → 1, isolated dimers). The grid uses 14 gaps from
0.005 to 1.7 J1, the same physical gaps at every N.

**Overlap knob: one layer of staggered Rz.** The trial is exp(−iε Σ(−1)^i Z_i)|E0⟩, with ε solved for each
target γ = |⟨E0|trial⟩|². It conserves S^z and moves the missing weight onto the low triplets, E1 included,
which is the hard case for a filter. The grid uses 14 overlaps from 0.005 to 0.99.

**Filter.** One ancilla. Pulse t applies cos(t H_s) on success, with H_s = (H − E0)/(Emax − E0). The dyadic
schedule uses t = τ 2^k for k < K, the ladder repeated r times, τ = π/(1 + Δ_s). Time evolution is
second-order Strang with n steps per τ.

**Measured cost.** For each (gap, γ, ε) the code finds the cheapest (K, r, n) whose simulated 1 − F ≤ ε,
where F is the fidelity with the exact ground state. A choice only counts if the same (K, r) still meets ε
at the next finer n, so a lucky Trotter cancellation is not counted. n is raised until finer steps stop
helping.

**CNOTs.** Each ancilla-controlled Heisenberg bond is one 6-CNOT block, and Strang merges the outer half
steps: a pulse of k steps costs 6(2k(B − 1) + 1) CNOTs, with B = 2N − 3 bonds.

**Protocols.** Both are optimised separately, so each picks its own cheapest filter.
- RUS, repeat until success: cost = CNOT per attempt / p.
- AA, amplitude amplification: k Grover rounds, (2k + 1) uses of the attempt, success sin²((2k+1)θ) with
  sin θ = √p. It needs a coherent attempt (no mid-circuit measurement, so one ancilla per pulse) and an
  invertible trial circuit. Each round also pays for two reflections: R_good, a sign flip of "all K·r
  ancillas = 0", and R_0, a sign flip of |0…0⟩ on all N + K·r qubits. Each is an n-qubit controlled Z,
  counted as 12n CNOTs. Cost = ((2k+1)·C_attempt + k·C_reflections) / sin²((2k+1)θ). AA is used only where
  it beats RUS.

`attempt` in the fits is the CNOT count of one attempt of the RUS-optimal filter.

## How it is fast enough for N = 16

The Trotterised filter is never time-stepped. One Strang step S is complex symmetric and unitary, so
S = O e^{−iφ} O^T with O real orthogonal. A Trotterised pulse of m steps is then exactly O cos(mφ) O^T,
and every (K, r) at a given n costs one vector product. O comes from a real symmetric eigensolver. Spin flip
commutes with S, which splits it into two blocks of half the size. At N = 16 the S^z = 0 sector has 12870
states, each block 6435, and one decomposition takes about 2 minutes. `tests/test_engine.py` checks this against a
direct gate-by-gate simulation to 1e−12, and `verify_qiskit.py` checks it against real Qiskit circuits.

Rough run times on 4 cores: N ≤ 12 under a minute each, N = 14 a few minutes, N = 16 a few hours.
Memory peaks at about 4 GB at N = 16.
