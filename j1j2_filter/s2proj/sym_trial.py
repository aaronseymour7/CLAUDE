"""sym_trial.py -- SU(2)-symmetric trial states: nearest-neighbour singlet pairs followed by layers of SU(2)-invariant two-qubit gates
G(theta) = exp(-i theta sigma.sigma / 4) = e^{i theta/4} (cos(theta/2) - i sin(theta/2) SWAP).  Total spin stays J = 0 exactly (p_S = 1), and so does every
Lie-Trotter circuit of the filter (each bond's XX, YY, ZZ factors commute and commute with the ancilla Z), so the energy filter only has to resolve the SINGLET gap
with no S^2 pulses.  Angles are optimised classically (statevector) to maximise the overlap with the exact ground state.
usage: python s2proj/sym_trial.py --bc obc pbc --N 8 12 16 --J2 0 0.2411 --layers 1 2 3"""
import argparse, json, os, sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / "pbc_corr")); sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT.parent)); os.chdir(ROOT)
from pbc_ctx import ed_inputs


def singlet_pairs(N):
    v = np.array([1.0 + 0j])
    for _ in range(N // 2):                       # qubits (2p, 2p+1): (|01> - |10>)/sqrt2 ; qiskit ordering: later qubits more significant
        v = np.kron(np.array([0, 1, -1, 0], complex) / np.sqrt(2), v)
    return v


def swap_bonds(vec, N, i, j):
    """SWAP_ij |psi>."""
    s = np.arange(1 << N); bi, bj = (s >> i) & 1, (s >> j) & 1
    t = np.where(bi != bj, s ^ ((1 << i) | (1 << j)), s)
    out = np.empty_like(vec); out[t] = vec; return out


def layer_bonds(N, parity, pbc):
    bl = [(i, i + 1) for i in range(parity, N - 1, 2)]
    if pbc and parity == 1:
        bl.append((N - 1, 0))
    return bl


def build_state(theta, N, layers, pbc):
    v = singlet_pairs(N); k = 0
    for l in range(layers):
        for (i, j) in layer_bonds(N, (l + 1) % 2, pbc):               # first layer acts on the bonds BETWEEN the singlet pairs
            th = theta[k]; k += 1
            v = np.exp(1j * th / 4) * (np.cos(th / 2) * v - 1j * np.sin(th / 2) * swap_bonds(v, N, i, j))
    return v


def n_params(N, layers, pbc):
    return sum(len(layer_bonds(N, (l + 1) % 2, pbc)) for l in range(layers))


def optimise(psi, N, layers, pbc, restarts=6, seed=0):
    rng = np.random.default_rng(seed); npar = n_params(N, layers, pbc); best = (0.0, None)
    f = lambda th: -abs(np.vdot(psi, build_state(th, N, layers, pbc))) ** 2
    for r in range(restarts):
        x0 = rng.uniform(-0.8, 0.8, npar) if r else np.full(npar, 0.4)
        res = minimize(f, x0, method="L-BFGS-B", options=dict(maxiter=400))
        if -res.fun > best[0]:
            best = (-res.fun, res.x)
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--bc", nargs="+", default=["obc"]); ap.add_argument("--N", type=int, nargs="+", default=[8]); ap.add_argument("--J2", type=float, nargs="+", default=[0.0])
    ap.add_argument("--layers", type=int, nargs="+", default=[1, 2, 3]); a = ap.parse_args(); out = {}
    for bc in a.bc:
        for N in a.N:
            for J2 in a.J2:
                psi = ed_inputs(N, 1.0, J2, bc == "pbc")["psi"]
                for L in a.layers:
                    g, th = optimise(psi, N, L, bc == "pbc"); out[f"{bc}_N{N}_J2{J2:g}_L{L}"] = dict(gamma=g, theta=th.tolist())
                    print(f"{bc} N={N} J2={J2:g} layers={L}: gamma_sym = {g:.4f}  ({n_params(N, L, bc == 'pbc')} angles)", flush=True)
    (HERE / "results" / "sym_trial.json").write_text(json.dumps(out))
