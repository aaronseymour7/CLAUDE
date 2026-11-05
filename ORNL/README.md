# ORNL: measured cost of a ground-state filter vs. gap and overlap

A small, standalone study of what a Trotterised cosine filter costs, and how the spectral gap and the
trial-state overlap drive that cost. The example is the open 1D J1-J2 chain. Everything is measured by exact
simulation of the circuit; there are no bounds or certificates.

It needs only numpy, scipy and matplotlib, and nothing outside this folder.

```
pip install -r requirements.txt
python test_filtercost.py        # 6 checks, a few seconds
python run_maps.py               # N = 8 10 12 14 16, eps = 1e-2 1e-3  -> results/map_N*.json
python analyze.py                # fits -> results/fits.csv, figures -> figs/
```

`python run_maps.py --N 8 10 12 --eps 1e-2` runs a subset. A run writes its JSON after every gap value and
resumes where it stopped.

## Files

| File | What it does |
|--|--|
| `filtercost.py` | The whole engine: model, the two knobs, Trotter simulation, measured search, CNOT count, RUS and AA cost |
| `run_maps.py` | Runs the (gap, overlap) grid for each N and eps |
| `analyze.py` | Power-law fits per N and over all N, and the figures |
| `test_filtercost.py` | Checks the engine against direct simulation, `expm` and a Grover toy model |
| `results/` | One JSON per N, and `fits.csv` |
| `figs/` | `fig_knobs`, `fig_maps_eps*`, `fig_cuts_eps*`, `fig_exponents`, `fig_local_exponents` |

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
  sin θ = √p. It needs a coherent attempt (no mid-circuit measurement) and an invertible trial circuit.
  It is used only where it beats RUS.

`attempt` in the fits is the CNOT count of one attempt of the RUS-optimal filter.

## How it is fast enough for N = 16

The Trotterised filter is never time-stepped. One Strang step S is complex symmetric and unitary, so
S = O e^{−iφ} O^T with O real orthogonal. A Trotterised pulse of m steps is then exactly O cos(mφ) O^T,
and every (K, r) at a given n costs one vector product. O comes from a real symmetric eigensolver. Spin flip
commutes with S, which splits it into two blocks of half the size. At N = 16 the S^z = 0 sector has 12870
states, each block 6435, and one decomposition takes about 2 minutes. The test file checks this against a
direct gate-by-gate simulation to 1e−12.

Rough run times on 4 cores: N ≤ 12 under a minute each, N = 14 a few minutes, N = 16 a few hours.
Memory peaks at about 4 GB at N = 16.
