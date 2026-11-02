# Environment for frontier/*.slurm. Module names follow OLCF's Frontier documentation; check `module avail`
# if a version has moved. Creates the virtual environment on first use (run once on a login node first).
module load PrgEnv-gnu rocm cray-python
VENV=${VENV:-$HOME/venv-sbpipe}
if [ ! -d "$VENV" ]; then
    python3 -m venv "$VENV"
    source "$VENV/bin/activate"
    pip install numpy scipy quimb cotengra autoray
    # CuPy for ROCm, built against the loaded rocm module for the MI250X (gfx90a)
    CUPY_INSTALL_USE_HIP=1 ROCM_HOME=$ROCM_PATH HCC_AMDGPU_TARGET=gfx90a pip install cupy
else
    source "$VENV/bin/activate"
fi
export PYTHONPATH="${REPO:-$HOME/claude/j1j2_filter}:$PYTHONPATH"
