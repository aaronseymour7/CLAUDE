"""mps_circuit_vs_layers.py -- the competitor: infidelity of the pipeline's own MPS -> circuit (mps_to_circuit, approximate, L layers) versus L, with RZZ-basis gate counts.
usage: python s2proj/mps_circuit_vs_layers.py --bc obc pbc --N 16 --J2 0.2411 --layers 1 2 3 4 6 8"""
import argparse, json, os, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
import quimb.tensor as qtn
from mps_to_circuit import mps_to_circuit
from qiskit.quantum_info import Statevector
from hamiltonians import reorder_axes
from pbc_ctx import ed_inputs
from core.resources import resource_costs, RZZ_BASIS

ap = argparse.ArgumentParser(); ap.add_argument("--bc", nargs="+", default=["obc"]); ap.add_argument("--N", type=int, nargs="+", default=[16]); ap.add_argument("--J2", type=float, nargs="+", default=[0.2411])
ap.add_argument("--layers", type=int, nargs="+", default=[1, 2, 3, 4, 6, 8]); ap.add_argument("--chi", type=int, default=96); a = ap.parse_args(); out = {}
for bc in a.bc:
    for N in a.N:
        for J2 in a.J2:
            psi = ed_inputs(N, 1.0, J2, bc == "pbc")["psi"]
            mps = qtn.MatrixProductState.from_dense(reorder_axes(psi), dims=[2] * N); mps.compress(cutoff=1e-12, max_bond=a.chi)
            for L in a.layers:
                t0 = time.time(); qc = mps_to_circuit(mps.arrays, method="approximate", shape="lpr", num_layers=L)
                g = float(abs(np.vdot(psi, Statevector(qc).data)) ** 2); c = resource_costs(qc, verbose=False, basis=RZZ_BASIS)
                out[f"{bc}_N{N}_J2{J2:g}_L{L}"] = dict(gamma=g, infidelity=1 - g, rzz=c["rzz"], nc=c["nc_total"], depth=c["depth"], seconds=time.time() - t0)
                print(f"{bc} N={N} J2={J2:g} L={L}: 1-F = {1 - g:.3e}  RZZ = {c['rzz']:.0f}  non-Clifford rot = {c['nc_total']:.0f}  depth {c['depth']:.0f}  ({time.time() - t0:.0f}s)", flush=True)
(HERE / "results" / "mps_circuit_vs_layers.json").write_text(json.dumps(out))
