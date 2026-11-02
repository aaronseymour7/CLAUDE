#!/bin/bash
# Stage 2: planned filter circuits on GPU, one case per GCD (8 per node). Run after inputs.slurm:
#   sbatch frontier/circuits.slurm
# 2 couplings x 7 targets x 2 schedules = 28 cases; 4 nodes = 32 GCDs runs them all at once.
#SBATCH -A PROJID
#SBATCH -J sbpipe-circuits
#SBATCH -p batch
#SBATCH -N 4
#SBATCH -t 06:00:00
#SBATCH -o sbpipe-circuits-%j.out

LX=${LX:-6}
OUT=${OUT:-$MEMBERWORK/PROJID/sbpipe}
REPO=${REPO:-$HOME/claude/j1j2_filter}
# curves of the next-smaller cylinder: for LX=6 the 4x5 results in the repo, for LX=7 the LX=6 output of this job
PREV=${PREV:-"$REPO/sbpipe_study/results/cross20.jsonl $REPO/sbpipe_study/results/transfer20.jsonl"}

source "$REPO/frontier/env.sh"
export OMP_NUM_THREADS=7

srun -N 4 -n 32 -c 7 --gpus-per-task=1 --gpu-bind=closest \
    python3 "$REPO/frontier/run_large.py" circuits --Lx "$LX" --J2 0 0.3 --prev $PREV --out "$OUT"

cat "$OUT"/circuits_cyl${LX}.rank*.jsonl > "$OUT/circuits_cyl${LX}.jsonl"
