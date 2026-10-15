"""run_rbm.py -- variational (RBM-type) energies of the OPEN J1-J2 chain with NetKet.

Follows the NetKet tutorial https://netket.readthedocs.io/en/stable/tutorials/gs-j1j2.html (source: docs/tutorials/gs-j1j2.ipynb):
custom graph with colour-1 (J1) and colour-2 (J2) bonds, GraphOperator, Spin hilbert space with fixed total S^z, the tutorial's
"FFNN" ansatz (complex Dense layer with 2N hidden units -> log cosh -> sum, i.e. an RBM with alpha=2), MetropolisExchange sampler
(d_max=2), SGD(lr=0.01) + stochastic reconfiguration (VMC_SR, diag_shift=0.01).
Differences from the tutorial (all deliberate):
  * OPEN boundary conditions (no wrap-around bonds), as in the rest of this repository;
  * Hamiltonian normalised as sum J S.S (S=sigma/2), i.e. the tutorial's Pauli form divided by 4, so energies match the DMRG/ED files;
  * two sectors:  S^z=0  -> lowest level of the S^z=0 sector = E0 (singlet ground state);
                  S^z=1  -> lowest level of the S^z=1 sector = lowest S=1 (triplet) level;
  * Marshall sign rotation on the nearest-neighbour exchange (as in the tutorial) switchable with --marshall 0/1.
The RBM energy is a variational UPPER bound on the lowest level of its sector, so (E_rbm - E_exact) >= 0 up to sampling error.

usage: python run_rbm.py --N 12 --J2 0 0.2411 0.5 --sz 0 1 [--iters 600 --samples 1008 --marshall 1 --ed]
output: results/rbm_N{N}_J2{J2}_sz{sz}.json  (one file per run)
"""
import os, json, time, argparse
os.environ.setdefault("JAX_PLATFORM_NAME", "cpu")
from pathlib import Path
import numpy as np
import netket as nk
import netket.nn as nknn
import flax.linen as nn
import jax.numpy as jnp

OUT = Path(__file__).resolve().parent / "results"


class RBM(nn.Module):                       # the tutorial's "FFNN" (== RBM with hidden density alpha)
    alpha: int = 2

    @nn.compact
    def __call__(self, x):
        x = nn.Dense(features=self.alpha * x.shape[-1], use_bias=True, param_dtype=np.complex128,
                     kernel_init=nn.initializers.normal(stddev=0.01), bias_init=nn.initializers.normal(stddev=0.01))(x)
        return jnp.sum(nknn.log_cosh(x), axis=-1)


class RBMSymm(nn.Module):                   # same RBM, but symmetrised over the chain's automorphisms (reflection for an open chain)
    perms: object = None
    alpha: int = 2
    flip: bool = False                          # also symmetrise over the global spin flip x -> -x (valid in the S^z=0 sector)

    @nn.compact
    def __call__(self, x):
        P = jnp.asarray(self.perms)             # (n_sym, N)
        xs = x[..., P]                          # (..., n_sym, N)
        if self.flip:
            xs = jnp.concatenate([xs, -xs], axis=-2)
        h = nn.Dense(features=self.alpha * x.shape[-1], use_bias=True, param_dtype=np.complex128,
                     kernel_init=nn.initializers.normal(stddev=0.01), bias_init=nn.initializers.normal(stddev=0.01))(xs)
        from jax.scipy.special import logsumexp
        return logsumexp(jnp.sum(nknn.log_cosh(h), axis=-1), axis=-1) - jnp.log(xs.shape[-2])


def ref_energy(N, J2, sz):
    """Exact reference from results/low_levels_ed.json (S^z=0 -> E0, S^z=1 -> lowest S=1 level) when the grid contains J2."""
    f = Path(__file__).resolve().parents[1] / "results" / "low_levels_ed.json"
    if f.exists():
        for J, e0, s1, t in json.load(open(f)).get(str(N), []):
            if abs(J - J2) < 1e-9:
                return e0 if sz == 0 else (t if sz == 1 else None)
    return None


def build(N, J2, sz, marshall=True):
    J = [1.0, J2]
    edges = []
    for i in range(N - 1):
        edges.append([i, i + 1, 1])
    for i in range(N - 2):
        edges.append([i, i + 2, 2])
    g = nk.graph.Graph(edges=edges)
    mszsz = np.kron([[1, 0], [0, -1]], [[1, 0], [0, -1]])
    exch = np.asarray([[0, 0, 0, 0], [0, 0, 2, 0], [0, 2, 0, 0], [0, 0, 0, 0]])
    s = -1.0 if marshall else 1.0               # Marshall rotation flips the sign of the J1 exchange only (J2 joins equal sublattices)
    ops = [(J[0] * mszsz / 4).tolist(), (J[1] * mszsz / 4).tolist(), (s * J[0] * exch / 4).tolist(), (J[1] * exch / 4).tolist()]
    hi = nk.hilbert.Spin(s=0.5, total_sz=sz, N=N)
    return g, hi, nk.operator.GraphOperator(hi, graph=g, bond_ops=ops, bond_ops_colors=[1, 2, 1, 2])


def run(N, J2, sz, iters, samples, marshall, alpha, lr, diag_shift, seed, do_ed, model='rbm', sched='const', chains=16, discard=None, flip=False):
    g, hi, H = build(N, J2, sz, marshall)
    sa = nk.sampler.MetropolisExchange(hilbert=hi, graph=g, d_max=2, n_chains=chains)
    if model == "symm":
        ma = RBMSymm(perms=tuple(map(tuple, np.array(g.automorphisms().to_array()).tolist())), alpha=alpha, flip=flip)
    else:
        ma = RBM(alpha=alpha)
    vs = nk.vqs.MCState(sa, ma, n_samples=samples, n_discard_per_chain=discard, seed=seed)
    import optax
    rate = lr if sched == "const" else optax.cosine_decay_schedule(lr, iters, alpha=0.05)
    gs = nk.driver.VMC_SR(hamiltonian=H, optimizer=nk.optimizer.Sgd(learning_rate=rate), variational_state=vs, diag_shift=diag_shift)
    log = nk.logging.RuntimeLog()
    t0 = time.time(); gs.run(n_iter=iters, out=log, show_progress=False); dt = time.time() - t0
    e = np.real(np.array(log.data["Energy"]["Mean"]))
    tail = e[-max(1, iters // 10):]
    vs.n_samples = max(samples * 8, 8192)       # final, higher-statistics evaluation of the trained state
    fin = vs.expect(H)
    r = dict(N=N, J2=J2, sz=sz, marshall=marshall, alpha=alpha, iters=iters, samples=samples, n_params=int(vs.n_parameters),
             E_tail_mean=float(tail.mean()), E_tail_std=float(tail.std()), E_final=float(np.real(fin.mean)), E_final_err=float(fin.error_of_mean),
             E_final_var=float(fin.variance), R_hat=float(fin.R_hat), tau_corr=float(fin.tau_corr), model=model, sched=sched, chains=chains, lr=lr, diag_shift=diag_shift, seconds=dt, energy_curve=e.tolist())
    ref = ref_energy(N, J2, int(sz))
    if ref is not None:
        r["E_exact"] = ref; r["rel_err"] = (r["E_final"] - ref) / abs(ref)
    if do_ed and "E_exact" not in r:
        r["E_exact"] = float(nk.exact.lanczos_ed(H, k=1)[0])
        r["rel_err"] = (r["E_final"] - r["E_exact"]) / abs(r["E_exact"])
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, nargs="+", default=[12]); ap.add_argument("--J2", type=float, nargs="+", default=[0.0])
    ap.add_argument("--sz", type=float, nargs="+", default=[0.0, 1.0])
    ap.add_argument("--iters", type=int, default=600); ap.add_argument("--samples", type=int, default=1008)
    ap.add_argument("--marshall", type=int, default=1); ap.add_argument("--alpha", type=int, default=2)
    ap.add_argument("--lr", type=float, default=0.01); ap.add_argument("--diag-shift", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=0); ap.add_argument("--model", default="rbm", choices=["rbm", "symm"]); ap.add_argument("--flip", type=int, default=0)
    ap.add_argument("--sched", default="const", choices=["const", "cos"]); ap.add_argument("--chains", type=int, default=16)
    ap.add_argument("--discard", type=int, default=None); ap.add_argument("--tag", default=""); ap.add_argument("--outdir", default="results")
    ap.add_argument("--ed", action="store_true", help="NetKet Lanczos ED reference (N<=20)")
    a = ap.parse_args(); OUT = Path(__file__).resolve().parent / a.outdir; OUT.mkdir(exist_ok=True)
    for N in a.N:
        for J2 in a.J2:
            for sz in a.sz:
                r = run(N, J2, sz, a.iters, a.samples, bool(a.marshall), a.alpha, a.lr, a.diag_shift, a.seed, a.ed, a.model, a.sched, a.chains, a.discard, bool(a.flip))
                json.dump(r, open(OUT / f"rbm_N{N}_J2{J2:g}_sz{sz:g}{a.tag}.json", "w"))
                ex = f"  exact {r['E_exact']:.6f}  rel.err {r['rel_err']:.2e}" if "E_exact" in r else ""
                print(f"{a.tag} N={N} J2={J2:g} Sz={sz:g}: E={r['E_final']:.6f}+-{r['E_final_err']:.1e}{ex}  R^={r['R_hat']:.3f} tau={r['tau_corr']:.1f} ({r['seconds']:.0f}s)", flush=True)
