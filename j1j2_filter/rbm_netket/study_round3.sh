#!/bin/bash
# Round 3 at N=16: push the learning rate / SR regularisation (round 2 best: lr 0.03, diag_shift 1e-3), and check the S^z=1 sector.
PY=${PY:-/tmp/nkenv/bin/python}; cd "$(dirname "$0")"; N=16
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1
C="--model symm --samples 2048 --iters 1500"
{
echo "0.6 r3_lr06ds1e3  $C --sz 0 --lr 0.06 --diag-shift 0.001"
echo "0.6 r3_lr10ds1e3  $C --sz 0 --lr 0.1 --diag-shift 0.001"
echo "0.6 r3_lr03ds3e4  $C --sz 0 --lr 0.03 --diag-shift 0.0003"
echo "0.6 r3_lr06ds3e4  $C --sz 0 --lr 0.06 --diag-shift 0.0003"
echo "0.6 r3_sz1_lr03ds $C --sz 1 --lr 0.03 --diag-shift 0.001"
echo "0.2411 r3_sz1_lr03ds $C --sz 1 --lr 0.03 --diag-shift 0.001"
} | xargs -P 4 -L 1 bash -c "\"$PY\" -W ignore run_rbm.py --N $N --J2 \$0 --tag _\$1 \${@:2} --outdir study_N16 2>&1 | grep -E '^ *_|Error|error' "
