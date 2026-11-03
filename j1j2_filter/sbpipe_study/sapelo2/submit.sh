#!/bin/bash
# Usage (from the repository root):  bash j1j2_filter/sbpipe_study/sapelo2/submit.sh STUDY MEM TIME
#   STUDY = chain | cyl ; MEM e.g. 8G ; TIME e.g. 3-00:00:00
# Optional: CHAIN_N=16,18,20 (chain sizes), CYL_LX=5 (cylinder lengths). Use the same values for every
# submit and for `merge`, since they fix the task numbering.
set -e
STUDY=$1; MEM=${2:-8G}; TIME=${3:-2-00:00:00}
D=j1j2_filter/sbpipe_study
export PYTHONPATH=$PWD/j1j2_filter
mkdir -p logs
NP=$(python $D/run_cluster.py count $STUDY prep)
NC=$(python $D/run_cluster.py count $STUDY case)
P=$(sbatch --parsable --array=0-$((NP-1)) --mem=$MEM --time=$TIME \
    --export=ALL,STAGE=prep,STUDY=$STUDY $D/sapelo2/task.sbatch)
C=$(sbatch --parsable --array=0-$((NC-1))%100 --mem=$MEM --time=$TIME --dependency=afterok:$P \
    --export=ALL,STAGE=case,STUDY=$STUDY $D/sapelo2/task.sbatch)
echo "prep $P ($NP tasks), cases $C ($NC tasks)"
