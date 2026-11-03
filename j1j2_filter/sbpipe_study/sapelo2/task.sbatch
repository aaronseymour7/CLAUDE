#!/bin/bash
#SBATCH --job-name=sbpipe
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%A_%a.out
# Submitted by submit.sh with --array, --mem, --time and the variables STAGE, STUDY (and CHAIN_N / CYL_LX).
module load Python/3.11.3-GCCcore-12.3.0
source $HOME/venv-sbpipe/bin/activate
export OMP_NUM_THREADS=1 PYTHONPATH=$SLURM_SUBMIT_DIR/j1j2_filter
cd $SLURM_SUBMIT_DIR/j1j2_filter/sbpipe_study
python run_cluster.py $STAGE $STUDY $SLURM_ARRAY_TASK_ID
