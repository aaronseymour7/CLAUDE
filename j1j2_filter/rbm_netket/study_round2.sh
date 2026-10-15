#!/bin/bash
# Round 2 at N=16: optimiser / symmetry variants of the symmetrised RBM (2048 samples, 1500 iterations).
PY=${PY:-/tmp/nkenv/bin/python}; cd "$(dirname "$0")"; N=16
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1
C="--model symm --samples 2048 --iters 1500 --sz 0"
{
echo "0.6 r2_ref     $C"
echo "0.6 r2_lr03    $C --lr 0.03"
echo "0.6 r2_ds1e3   $C --diag-shift 0.001"
echo "0.6 r2_lr03ds  $C --lr 0.03 --diag-shift 0.001"
echo "0.6 r2_flip    $C --flip 1"
echo "0.6 r2_nomar   $C --marshall 0"
echo "0.2411 r2_ref  $C"
echo "0.2411 r2_lr03ds $C --lr 0.03 --diag-shift 0.001"
echo "0.2411 r2_flip $C --flip 1"
} | xargs -P 4 -L 1 bash -c "\"$PY\" -W ignore run_rbm.py --N $N --J2 \$0 --tag _\$1 \${@:2} --outdir study_N16 2>&1 | grep -E '^ *_|Error|error' "
