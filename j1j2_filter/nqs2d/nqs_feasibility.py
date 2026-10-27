"""nqs_feasibility.py -- can a neural-network quantum state (NQS) serve the filter pipeline on the 2D J1-J2 model?

Open W x L square lattice, H = J1 sum_<ij> S_i.S_j + J2 sum_<<ij>> S_i.S_j (same geometry and qubit order as ../crossover2d/geometry.py:
q = x*W + y, columns of height W, so the snake MPS cut crosses W sites).  Everything is checked against exact diagonalisation.

  train : complex RBM (log psi = a.s + sum_j log cosh(b_j + W_j.s), times the Marshall sign) in the S^z = 0 sector, optimised by
          stochastic reconfiguration with exact local energies (full sum over the sector when it has <= 20000 states, otherwise
          i.i.d. samples drawn exactly from |psi|^2).  The same RBM in the S^z = 1 sector gives a variational triplet energy.
          Also: best MPS of bond dimension chi (TT-SVD truncation of the exact ground state, snake order) and its overlap.
  cost  : feed each trial's overlap into sbpipe.refine (eps = 1e-2, second order, cheaper of dyadic/optimized schedule) with exact
          energies, and compare the expected CNOTs per success.  The NQS has no known compact circuit, so its row gives the filter
          cost alone and the loading budget at which it would tie the best MPS-trial plan.
  refcheck : NQS as the classical reference instead of the trial: overlap of each MPS trial with the NQS vs with the exact state,
          and NQS energy and gap vs exact.
usage (from j1j2_filter/):  python nqs2d/nqs_feasibility.py train 4x4 --j2 0.5 --alpha 1 2 4 ;  python nqs2d/nqs_feasibility.py cost
"""
import argparse, glob, json, os, sys, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
RES = HERE / "results"


# ---------------------------------------------------------------- lattice and exact diagonalisation
def bonds(W, L):
    idx = lambda x, y: x * W + y
    b1, b2 = set(), set()
    for x in range(L):
        for y in range(W):
            if y + 1 < W: b1.add((idx(x, y), idx(x, y + 1)))
            if x + 1 < L: b1.add((idx(x, y), idx(x + 1, y)))
            if x + 1 < L and y + 1 < W: b2.add(tuple(sorted((idx(x, y), idx(x + 1, y + 1)))))
            if x + 1 < L and y - 1 >= 0: b2.add(tuple(sorted((idx(x, y), idx(x + 1, y - 1)))))
    return sorted(b1), sorted(b2)


def sector(W, L, nup, j2):
    """Sparse H in the sector with nup up spins; basis = ints, bit q = spin q (1 = up)."""
    N = W * L; b1, b2 = bonds(W, L)
    s = np.arange(1 << N, dtype=np.int64)
    pc = np.zeros(1 << N, dtype=np.int8)
    for q in range(N):
        pc += ((s >> q) & 1).astype(np.int8)
    st = s[pc == nup]; del s, pc
    d = len(st); diag = np.zeros(d); R, C, V = [], [], []
    for J, bl in ((1.0, b1), (j2, b2)):
        if J == 0:
            continue
        for i, j in bl:
            bi, bj = (st >> i) & 1, (st >> j) & 1
            diag += J * np.where(bi == bj, 0.25, -0.25)
            src = np.where(bi != bj)[0]; tgt = np.searchsorted(st, st[src] ^ ((1 << i) | (1 << j)))
            R.append(tgt); C.append(src); V.append(np.full(len(src), 0.5 * J))
    H = sp.diags(diag) + sp.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(d, d))
    spins = (2 * ((st[:, None] >> np.arange(N)) & 1) - 1).astype(np.float64)
    return H.tocsr(), st, spins


def marshall(W, L, spins):
    """(-1)^(number of up spins on sublattice A), A = (x + y) even."""
    N = W * L; A = np.array([((q // W) + (q % W)) % 2 == 0 for q in range(N)])
    n_up_A = ((spins[:, A] + 1) / 2).sum(1)
    return np.where(n_up_A % 2 == 0, 1.0, -1.0)


# ---------------------------------------------------------------- RBM
def logcosh(z):
    z = np.where(z.real < 0, -z, z)
    return z + np.log1p(np.exp(-2 * z)) - np.log(2.0)


class RBM:
    def __init__(self, N, alpha, seed=0):
        r = np.random.default_rng(seed); M = alpha * N
        self.N, self.M = N, M
        self.a = 1e-3 * (r.normal(size=N) + 1j * r.normal(size=N))
        self.b = 1e-2 * (r.normal(size=M) + 1j * r.normal(size=M))
        self.W = 1e-2 * (r.normal(size=(M, N)) + 1j * r.normal(size=(M, N)))

    @property
    def n_params(self):                 # real parameters
        return 2 * (self.N + self.M + self.M * self.N)

    def theta(self, s):
        return s @ self.W.T + self.b

    def log_psi(self, s):
        return s @ self.a + logcosh(self.theta(s)).sum(1)

    def O(self, s):
        t = np.tanh(self.theta(s))
        return np.concatenate([s, t, (t[:, :, None] * s[:, None, :]).reshape(len(s), -1)], axis=1)

    def update(self, dp):
        N, M = self.N, self.M
        self.a = self.a - dp[:N]; self.b = self.b - dp[N:N + M]; self.W = self.W - dp[N + M:].reshape(M, N)


def amplitudes(rbm, spins, sign):
    lp = rbm.log_psi(spins); lp = lp - lp.real.max()
    psi = sign * np.exp(lp); return psi / np.linalg.norm(psi)


def train_rbm(H, spins, sign, alpha, iters, lr, seed, n_samples=10000, log_every=50, label=""):
    rbm = RBM(spins.shape[1], alpha, seed); rng = np.random.default_rng(seed + 1)
    d = len(spins); full = d <= 20000; hist = []
    E_best, best = np.inf, None
    for it in range(iters):
        psi = amplitudes(rbm, spins, sign); Hpsi = H @ psi
        E = float(np.vdot(psi, Hpsi).real)
        if it % log_every == 0 or it == iters - 1:
            hist.append((it, E))
            print(f"  {label} it {it:5d}  E = {E:.8f}", flush=True)
        if E < E_best:
            E_best, best = E, (rbm.a.copy(), rbm.b.copy(), rbm.W.copy())
        p = np.abs(psi) ** 2
        if full:
            idx, w = np.arange(d), p
        else:
            idx = rng.choice(d, size=n_samples, p=p / p.sum()); w = np.full(n_samples, 1.0 / n_samples)
        Eloc = Hpsi[idx] / psi[idx]
        O = rbm.O(spins[idx])
        Om = w @ O; Oc = O - Om; Ec = Eloc - w @ Eloc
        S = (Oc.conj().T * w) @ Oc; F = (Oc.conj().T * w) @ Ec
        lam = max(1e-4, 1e-2 * 0.995 ** it)
        S[np.diag_indices_from(S)] += lam * (1 + S.diagonal().real) if it < 200 else lam
        dp = np.linalg.solve(S, F)
        rbm.update(lr * dp)
    rbm.a, rbm.b, rbm.W = best
    psi = amplitudes(rbm, spins, sign)
    return rbm, psi, hist


# ---------------------------------------------------------------- MPS truncation of the exact state
def mps_truncate(vec, N, chi):
    """Left-to-right TT-SVD truncated to bond dimension chi; returns the normalised dense state."""
    cores, rest, r = [], vec.reshape(1, -1), 1
    for k in range(N - 1):
        m = rest.reshape(r * 2, -1)
        U, s, Vh = np.linalg.svd(m, full_matrices=False)
        keep = min(chi, len(s)); U, s, Vh = U[:, :keep], s[:keep], Vh[:keep]
        cores.append(U.reshape(r, 2, keep)); rest = s[:, None] * Vh; r = keep
    cores.append(rest.reshape(r, 2, 1))
    out = cores[0].reshape(2, -1)
    for c in cores[1:]:
        out = (out @ c.reshape(c.shape[0], -1)).reshape(-1, c.shape[2])
    out = out.reshape(-1); return out / np.linalg.norm(out)


def embed(st, v, N):
    full = np.zeros(1 << N, dtype=complex); full[st] = v
    # bit q = spin q, so the big-endian reshape of `full` lists sites N-1 ... 0: the snake order reversed (same chain)
    return full


# ---------------------------------------------------------------- commands
def cmd_train(a):
    W, L = (int(v) for v in a.lattice.split("x")); N = W * L; j2 = a.j2; t0 = time.time()
    H0, st0, s0 = sector(W, L, N // 2, j2); H1, st1, s1 = sector(W, L, N // 2 + 1, j2)
    w0, v0 = sla.eigsh(H0, k=4, which="SA", tol=1e-12); o = np.argsort(w0); w0, v0 = w0[o], v0[:, o]
    w1 = np.sort(sla.eigsh(H1, k=2, which="SA", tol=1e-12, return_eigenvectors=False))
    g = v0[:, 0]; sg0, sg1 = marshall(W, L, s0), marshall(W, L, s1)
    out = dict(W=W, L=L, N=N, J2=j2, dim_sz0=len(st0), E0=float(w0[0]), E1_sz0=float(w0[1]), E_triplet=float(w1[0]),
               marshall_weight=float(np.sum(np.abs(g) ** 2 * (np.sign(g * sg0) == np.sign(np.sum(g * sg0))))),
               rbm=[], mps=[])
    print(f"{W}x{L} J2={j2}: E0={w0[0]:.8f}  E1(Sz=0)={w0[1]:.8f}  triplet={w1[0]:.8f}  Marshall-sign weight {out['marshall_weight']:.4f}")
    gfull = embed(st0, g, N)
    for chi in a.chis:
        m = mps_truncate(gfull, N, chi)
        out["mps"].append(dict(chi=chi, gamma=float(abs(np.vdot(gfull, m)) ** 2), n_params=int(min(chi, 2 ** (N // 2)) ** 2 * 2 * N),
                               energy=float(np.vdot(m[st0], H0 @ m[st0]).real) if np.allclose(np.linalg.norm(m[st0]), 1, atol=1e-6) else None,
                               psi=None))
        out["mps"][-1].pop("psi"); np.save(RES / f"mps_{W}x{L}_J2{j2:g}_chi{chi}.npy", m)
        print(f"  MPS chi={chi:3d}: gamma = {out['mps'][-1]['gamma']:.6f}")
    np.save(RES / f"gs_{W}x{L}_J2{j2:g}.npy", gfull)
    for alpha in a.alpha:
        rbm, psi, hist = train_rbm(H0, s0, sg0, alpha, a.iters, a.lr, a.seed, label=f"RBM a={alpha} Sz=0")
        E = float(np.vdot(psi, H0 @ psi).real); var = float(np.vdot(H0 @ psi, H0 @ psi).real - E ** 2)
        gam = float(abs(np.vdot(g, psi)) ** 2)
        rbm1, psi1, hist1 = train_rbm(H1, s1, sg1, alpha, a.iters, a.lr, a.seed, label=f"RBM a={alpha} Sz=1")
        E1 = float(np.vdot(psi1, H1 @ psi1).real)
        np.save(RES / f"rbm_{W}x{L}_J2{j2:g}_a{alpha}.npy", embed(st0, psi, N))
        out["rbm"].append(dict(alpha=alpha, n_params=rbm.n_params, E=E, rel_err=(E - out["E0"]) / abs(out["E0"]), variance=var, gamma=gam,
                               E_triplet=E1, gap=E1 - E, gap_exact=out["E_triplet"] - out["E0"], history=hist, history_sz1=hist1))
        print(f"  RBM alpha={alpha} ({rbm.n_params} params): E={E:.6f} rel.err {out['rbm'][-1]['rel_err']:.2e}  var {var:.2e}  gamma {gam:.5f}  "
              f"gap {E1 - E:.4f} (exact {out['E_triplet'] - out['E0']:.4f})", flush=True)
    out["seconds"] = time.time() - t0
    (RES / f"train_{W}x{L}_J2{j2:g}.json").write_text(json.dumps(out, indent=1))


def cmd_cost(a):
    sys.path.insert(0, str(ROOT))
    import sbpipe
    from sbpipe.frontends.dmrg_quimb import mps_circuit_cost
    rows = []
    for f in sorted(glob.glob(str(RES / "train_*.json"))):
        r = json.load(open(f)); W, L, N, j2 = r["W"], r["L"], r["N"], r["J2"]
        if a.only and f"{W}x{L}_J2{j2:g}" not in a.only:
            continue
        b1, b2 = bonds(W, L); labels, coeffs, groups = [], [], []
        for gi, (i, j, J) in enumerate([(i, j, 1.0) for i, j in b1] + [(i, j, j2) for i, j in b2 if j2]):
            for P in "XYZ":
                lab = ["I"] * N; lab[i] = lab[j] = P; labels.append("".join(lab)); coeffs.append(0.25 * J); groups.append(gi)
        H = sbpipe.PauliSum(labels, np.array(coeffs), groups)
        est = sbpipe.EnergyEstimates(E0=r["E0"], E0_err=1e-6, E1=r["E1_sz0"], E1_err=1e-6, source="ED")
        cfg = sbpipe.Config(eps=a.eps, schedule="auto", order=2)

        def plan(gamma, cost, name):
            p = sbpipe.refine(H, est, sbpipe.Ansatz(gamma=gamma, cost=cost, name=name), cfg); s = p.summary()
            return dict(trial=name, gamma=gamma, trial_cnot=cost.get("cnot", 0), m=s["m"], schedule=s["schedule"],
                        filter_cnot=s["res_filter_cnot"], p_succ=s["res_p_succ"], expected_cnot=s["res_expected_cnot"],
                        expected_t=s["res_expected_t"])
        out = dict(W=W, L=L, N=N, J2=j2, eps=a.eps, delta=None, plans=[])
        for m in r["mps"]:
            if m["gamma"] < 1e-3:
                continue
            out["plans"].append(plan(m["gamma"], mps_circuit_cost(N, m["chi"]), f"MPS chi={m['chi']}"))
        best = min(out["plans"], key=lambda q: q["expected_cnot"])
        for q in r["rbm"]:
            pl = plan(q["gamma"], {}, f"RBM alpha={q['alpha']} (loading not counted)")
            pl["breakeven_loading_cnot"] = best["expected_cnot"] * pl["p_succ"] - pl["filter_cnot"]
            out["plans"].append(pl)
        # how the filtered cost depends on the overlap alone, for a free trial
        out["gamma_scan"] = [plan(gm, {}, f"free trial gamma={gm}") for gm in (0.05, 0.1, 0.2, 0.5, 0.8, 0.95)]
        out["best_mps"] = best["trial"]
        print(f"\n{W}x{L} J2={j2:g}  (eps={a.eps})")
        for q in out["plans"] + out["gamma_scan"]:
            extra = f"  break-even loading {q['breakeven_loading_cnot']:.3g} CNOT" if "breakeven_loading_cnot" in q else ""
            print(f"  {q['trial']:38s} gamma={q['gamma']:.4f}  trial {q['trial_cnot']:9.3g}  filter/attempt {q['filter_cnot']:9.3g}  "
                  f"p_succ {q['p_succ']:.3f}  per success {q['expected_cnot']:9.3g}{extra}")
        (RES / f"cost_{W}x{L}_J2{j2:g}.json").write_text(json.dumps(out, indent=1))
        rows.append(out)


def cmd_refcheck(a):
    for f in sorted(glob.glob(str(RES / "train_*.json"))):
        r = json.load(open(f)); W, L, j2 = r["W"], r["L"], r["J2"]
        g = np.load(RES / f"gs_{W}x{L}_J2{j2:g}.npy"); out = []
        for q in r["rbm"]:
            nqs = np.load(RES / f"rbm_{W}x{L}_J2{j2:g}_a{q['alpha']}.npy")
            for m in r["mps"]:
                t = np.load(RES / f"mps_{W}x{L}_J2{j2:g}_chi{m['chi']}.npy")
                out.append(dict(alpha=q["alpha"], chi=m["chi"], gamma_exact=float(abs(np.vdot(g, t)) ** 2),
                                gamma_vs_nqs=float(abs(np.vdot(nqs, t)) ** 2), nqs_fidelity=q["gamma"]))
        r["refcheck"] = out; Path(f).write_text(json.dumps(r, indent=1))
        print(f"{W}x{L} J2={j2:g}")
        for o in out:
            print(f"  RBM a={o['alpha']} (F={o['nqs_fidelity']:.4f})  MPS chi={o['chi']:3d}: gamma exact {o['gamma_exact']:.4f}  vs NQS {o['gamma_vs_nqs']:.4f}")


if __name__ == "__main__":
    P = argparse.ArgumentParser(); sub = P.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train"); t.add_argument("lattice"); t.add_argument("--j2", type=float, default=0.5)
    t.add_argument("--alpha", type=int, nargs="+", default=[1, 2, 4]); t.add_argument("--chis", type=int, nargs="+", default=[1, 2, 4, 8, 16, 32])
    t.add_argument("--iters", type=int, default=1500); t.add_argument("--lr", type=float, default=0.05); t.add_argument("--seed", type=int, default=0)
    c = sub.add_parser("cost"); c.add_argument("--eps", type=float, default=1e-2); c.add_argument("--only", nargs="*", help="e.g. 4x5_J20")
    sub.add_parser("refcheck")
    a = P.parse_args(); {"train": cmd_train, "cost": cmd_cost, "refcheck": cmd_refcheck}[a.cmd](a)
