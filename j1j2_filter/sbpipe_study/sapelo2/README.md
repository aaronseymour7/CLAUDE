# Large-N studies on Sapelo2

Two studies from the review of the measured draft, too large for a laptop. Each case is one SLURM
array task writing one JSON file, so partial runs are still usable.

| study | what it answers | systems | tasks |
| --- | --- | --- | --- |
| `chain` | Where SBC and QETU cross over as N grows | open J1–J2 chain, N = 16, 18, 20, 22, 24, gap held at 0.40 by tuning J2; χ = 16 MPS rotated to γ = 0.1, 0.5; ε = 1e-2, 1e-3, 1e-4; dyadic, optimized, QETU at τ = π/2, 3π/4, π/(1+δ) | 150 (30 per N) |
| `cyl` | Filter vs MPS alone beyond 16 sites, and a non-snake trial | 4×5 and 4×6 cylinders, J2 = 0, 0.3; snake MPS χ = 2–64, plus products of exact 4-site ring and 8-site block ground states; ε = 1e-2 … 1e-5; dyadic, optimized, QETU at π/(1+δ) | 544 (272 per length) |

## Setup (once, from the repository root)

    bash j1j2_filter/sbpipe_study/sapelo2/setup.sh

Edit the `module load Python/...` line in `setup.sh` and `task.sbatch` if that module name is not
available (`module spider Python`).

## Submit

From the repository root. Each call submits a `prep` array (ground state, E0, E1, J2 tuning, trial
states; one task per system) and a `case` array that waits for it.

    CHAIN_N=16,18,20 bash j1j2_filter/sbpipe_study/sapelo2/submit.sh chain 6G 2-00:00:00
    CHAIN_N=22       bash j1j2_filter/sbpipe_study/sapelo2/submit.sh chain 12G 5-00:00:00
    CHAIN_N=24       bash j1j2_filter/sbpipe_study/sapelo2/submit.sh chain 32G 7-00:00:00
    CYL_LX=5         bash j1j2_filter/sbpipe_study/sapelo2/submit.sh cyl 6G 2-00:00:00
    CYL_LX=6         bash j1j2_filter/sbpipe_study/sapelo2/submit.sh cyl 32G 7-00:00:00

`CHAIN_N` and `CYL_LX` fix the task numbering, so use the same value when merging. Rough cost per task
on one core: minutes at 16 qubits, about an hour at 20, many hours at 22, and a day or more at 24
(a Strang step at 24 qubits takes about 15 s). Start with N ≤ 20 and 4×5.

## Collect

    cd j1j2_filter/sbpipe_study
    CHAIN_N=16,18,20 PYTHONPATH=.. python run_cluster.py merge chain   # -> results/cluster_chain.jsonl
    CYL_LX=5         PYTHONPATH=.. python run_cluster.py merge cyl     # -> results/cluster_cyl.jsonl

Send back the `results/cluster_*.jsonl` files (or the `results/cluster/` folder). A merge with a
different `CHAIN_N` overwrites the jsonl, so rename it between merges, for example
`cluster_chain_N16-20.jsonl`. The ground-state cache (`results/cache_cluster/`, up to 0.5 GB per
24-qubit system) need not be sent.

## What each record holds

`gamma`, `gap`, `delta` (scaled gap), `trial_cnot`, `mps_alone` (trial already meets ε), and
`result`: CNOTs per success, success probability, Strang steps, coarsest passing step (`dt_phys`, in
1/J1), filter length, plus `trotter_log`, every (step, total steps, infidelity) evaluated, which also
checks that passing is monotone in the step.

## Check before submitting

    PYTHONPATH=j1j2_filter python j1j2_filter/sbpipe_study/bigsim.py

compares the matrix-free H|v>, Strang step and exact evolution with sbpipe on 10 and 8 qubits (all
differences should be ~1e-14).
