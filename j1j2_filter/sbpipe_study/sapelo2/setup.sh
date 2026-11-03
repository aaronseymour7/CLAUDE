#!/bin/bash
# One-time environment on Sapelo2 (run on a login or interactive node from the repository root).
module load Python/3.11.3-GCCcore-12.3.0          # any Python >= 3.10 module works; check `module spider Python`
python -m venv $HOME/venv-sbpipe
source $HOME/venv-sbpipe/bin/activate
pip install --upgrade pip
pip install numpy scipy quimb
python j1j2_filter/sbpipe_study/bigsim.py           # self-check: three numbers ~1e-14 per line
