#!/bin/bash
# Setup study at N=16 (exact references available): one run per (config, J2); 4 in parallel, 1 thread each.
PY=${PY:-/tmp/nkenv/bin/python}; cd "$(dirname "$0")"; N=${N:-16}
export XLA_FLAGS="--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1" OMP_NUM_THREADS=1
cfgs=("base:" "nomarshall:--marshall 0" "alpha1:--alpha 1" "alpha4:--alpha 4" "samples4096:--samples 4096" "cos:--sched cos" "symm:--model symm" "iters2000:--iters 2000")
for J2 in ${J2S:-0.2411 0.6}; do for c in "${cfgs[@]}"; do echo "${c%%:*} $J2 ${c#*:}"; done; done |
 xargs -P 4 -L 1 bash -c "\"$PY\" -W ignore run_rbm.py --N $N --J2 \$1 --sz 0 --tag _\$0 \${@:2} --outdir study_N$N 2>&1 | grep \"^ *_\"" 
