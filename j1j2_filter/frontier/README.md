# Larger cylinders on Frontier

These scripts run the paper's size-transfer test on 4×6 (24 qubits) and 4×7 (28 qubits) J1–J2 cylinders. They plan the SBC and QETU filters from DMRG inputs and simulate the Trotterized circuits exactly on Frontier's GPUs. Each circuit runs at the step taken from the next-smaller cylinder. The output is the actual behavior of each circuit: infidelity, success probability, post-selected energy and CNOT count. It is the same measurement as the 4×5 check in the paper (`sbpipe_study/run_transfer20.py`), one or two sizes further.

## What is here

| File | Purpose |
|---|---|
| `gpusim.py` | Statevector Strang product formula on the GPU (CuPy) or CPU (NumPy). Each Heisenberg bond is applied as one exact SWAP-power, giving the same circuit as `sbpipe.trotter.product_formula`. |
| `run_large.py` | The driver, with three stages: `inputs` (DMRG, CPU), `circuits` (GPU) and `selftest`. Run `python3 frontier/run_large.py -h` for options. |
| `env.sh` | Modules and a Python venv with numpy, scipy, quimb and CuPy built for ROCm/MI250X. |
| `inputs.slurm` | Stage 1: one node per J2. |
| `circuits.slurm` | Stage 2: one case per GCD, 28 cases on 4 nodes. |

## Steps

1. Clone the repository and check out this branch. On a login node, run `REPO=$PWD/j1j2_filter source j1j2_filter/frontier/env.sh` once; this builds the venv and CuPy.
2. Run `python3 frontier/run_large.py selftest`, first on a login node and then inside a short GPU allocation. Every number it prints should be about 1e-13 or smaller, and the backend should say `GPU (CuPy)` on a compute node. If the fused kernel does not compile under your ROCm version, set `SBPIPE_NOFUSE=1`; it is slower but uses only standard CuPy operations.
3. In both `.slurm` files, replace `PROJID` with your allocation and set `REPO` and `OUT`. Then run `sbatch frontier/inputs.slurm`.
4. When the inputs are done, run `sbatch frontier/circuits.slurm`. The results go to `$OUT/circuits_cyl6.jsonl`, one row per (J2, eps, schedule).
5. For 4×7, run both jobs again with `LX=7` and `PREV=$OUT/circuits_cyl6.jsonl`, for example `LX=7 PREV=... sbatch frontier/inputs.slurm`.

## What the rows mean

The step rule is the paper's: `dt_prev` is the coarsest passing step of the next-smaller cylinder. Each case is simulated at `dt_prev`, then `dt_prev/2`, and finer only if needed. The rule passes when the `dt_prev/2` point has `inf <= eps`. A pass at `dt_prev` itself means the step did not need to shrink.

The trial χ is the one with the cheapest certified plan among χ = 8, 16, 32, 64 that do not meet ε on their own. `--chi` fixes it instead.

`cert_cnot` is the a-priori certified cost, and `cnot` is the expected cost at each simulated step, both including the trial circuit.

## Limits

- **Size.** Memory is 16·2^N bytes per state vector: 256 MiB at 24 qubits and 4 GiB at 28. QETU keeps four vectors, which is well inside one 64 GB GCD at 4×7. A 4×8 cylinder (32 qubits) needs 64 GiB per vector, which no single GCD holds. That size needs a distributed statevector simulator, which these scripts do not provide.
- **Reference state.** At these sizes there is no exact ground state, so infidelity is measured against the DMRG ground state at the last rung of `--ladder` (default χ = 1024). The inputs json records `ref_err`, the infidelity between the last two rungs. Infidelities near `ref_err` are not resolved, so keep it well below the smallest ε, and extend the ladder if it is not.
- **Inputs.** The gap and overlap come from DMRG, as in the paper's planned route. They are not certified.
- **Run time.** These numbers are rough guesses, not measurements. DMRG to χ = 1024 on a 4×7 cylinder may take several hours per J2 on a CPU node, so raise `-t` in `inputs.slurm` if needed. The circuits are memory-bandwidth bound and should take minutes per case at 24 qubits and tens of minutes at 28.
