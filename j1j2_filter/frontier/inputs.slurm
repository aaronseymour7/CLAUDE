#!/bin/bash
# Stage 1: DMRG inputs on CPU, one J2 per node. Edit the account and paths, then:  sbatch frontier/inputs.slurm
#SBATCH -A PROJID
#SBATCH -J sbpipe-inputs
#SBATCH -p batch
#SBATCH -N 2
#SBATCH -t 12:00:00
#SBATCH -o sbpipe-inputs-%j.out

LX=${LX:-6}                                  # 4 x LX cylinder: 6 -> 24 qubits, 7 -> 28 qubits
OUT=${OUT:-$MEMBERWORK/PROJID/sbpipe}        # dense vectors are 2^N * 16 bytes each (4 GiB at 28 qubits)
REPO=${REPO:-$HOME/claude/j1j2_filter}

source "$REPO/frontier/env.sh"
export OMP_NUM_THREADS=56

# one task per node, each takes one J2 (SLURM_PROCID picks it)
srun -N 2 -n 2 -c 56 python3 "$REPO/frontier/run_large.py" inputs --Lx "$LX" --J2 0 0.3 --out "$OUT"
