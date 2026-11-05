"""Matrix-product-state version of the study, for chains beyond exact diagonalisation (needs quimb).

Reference (cfg.mps_reference): "ed" takes E0, E1, Emax and |E0> from exact diagonalisation (small N, turned
into an MPS), "dmrg" uses no exact spectrum at all, "auto" (default) picks ED up to cfg.mps_ed_max_N:
  DMRG         E0 and |E0> from DMRG; Emax from DMRG on -H; E1 from DMRG in S^z = 1 (penalty lam (S^z_tot - 1)^2),
               which is the lowest triplet. With SU(2) and a singlet ground state this is the gap the engine uses.
  trial        exp(-i eps sum_i (-1)^i Z_i) |E0>, one layer of single-site gates (bond dimension unchanged).
  overlap      gamma = |<E0|trial>|^2, and eps solved for a target gamma, as in the exact model.
  filter       the same Trotterised circuit. One pulse on ancilla a, post-selected on a = 0, leaves the system in
               (S^m + S^-m)/2 |psi> (H on the ancilla, every bond gate controlled by Z_a, H, project on 0). So each
               pulse runs the merged Strang gate sequence forward (+dt) and backward (-dt) on the system MPS, adds
               the two and compresses. No ancilla and no long-range gates are needed; the J2 bonds use quimb's
               swap+split. The norm of the result is the success probability, so p is exact up to truncation.
  fidelity     F = |<E0|v>|^2 / p, against the same reference |E0>.

All pulses are functions of the same Strang step, so they commute and the order of the ladder does not matter.
Truncation is set by cfg.mps_max_bond and cfg.mps_cutoff; `truncation_check` reruns a point at a larger bond.
"""
import numpy as np
import quimb.operator as qop
import quimb.tensor as qtn
from scipy.optimize import brentq

from .costs import aa_cost, cnot_per_attempt, cost, reflection_cnots
from .model import Chain, chain_bonds
from .model import delta_for_gap as ed_delta_for_gap
from .trotter import strang_sequence

SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)   # S_i.S_j = (2 SWAP - 1)/4


def _builder(N, bonds):
    J = {(i, j): x for i, j, x in bonds}
    return qop.heisenberg_from_edges(list(J), j=lambda i, j: J[(i, j)] if (i, j) in J else J[(j, i)])


def _dmrg(mpo, bond_dims, tol):
    d = qtn.DMRG2(mpo, bond_dims=list(bond_dims), cutoffs=1e-12)
    d.solve(tol=tol, verbosity=0, max_sweeps=40)
    return float(np.real(d.energy)), d.state


class MPSChain:
    """Same model and knobs as `model.Chain`, from DMRG. Attributes: N, delta, bonds, E0, E1, Emax, gap, W, gap_s,
    g0 (ground-state MPS), max_bond, cutoff."""

    def __init__(self, N, delta, J1=1.0, J2=0.2411, max_bond=128, cutoff=1e-12, reference="auto",
                 ed_max_N=16, dmrg_bonds=(16, 32, 64, 128), dmrg_tol=1e-10, penalty=2.0):
        assert N % 2 == 0, "even N only"
        self.N, self.delta, self.J1, self.J2 = N, delta, J1, J2
        self.max_bond, self.cutoff = max_bond, cutoff
        self.bonds = chain_bonds(N, delta, J1, J2)
        if reference == "auto":
            reference = "ed" if N <= ed_max_N else "dmrg"
        self.reference = reference
        if reference == "ed":
            self._reference_ed()
        elif reference == "dmrg":
            self._reference_dmrg(dmrg_bonds, dmrg_tol, penalty)
        else:
            raise ValueError("reference must be 'auto', 'ed' or 'dmrg'")
        self.gap, self.W = self.E1 - self.E0, self.Emax - self.E0
        self.gap_s = self.gap / self.W

    def _reference_ed(self):
        """Exact E0, E1, Emax and |E0> (as an MPS) from `model.Chain`."""
        ch = Chain(self.N, self.delta, self.J1, self.J2)
        self.E0, self.E1, self.Emax = ch.E0, ch.E1, ch.Emax
        full = np.zeros(1 << self.N, dtype=complex)
        full[ch.states] = ch.g0                               # same order as quimb: site 0 = most significant bit
        g0 = qtn.MatrixProductState.from_dense(full, dims=[2] * self.N, cutoff=1e-14)
        self.g0 = g0 / g0.norm()

    def _reference_dmrg(self, dmrg_bonds, dmrg_tol, penalty):
        N = self.N
        H = _builder(N, self.bonds)
        self.E0, g0 = _dmrg(H.build_mpo(), dmrg_bonds, dmrg_tol)
        self.g0 = g0 / g0.norm()
        neg = _builder(N, [(i, j, -x) for i, j, x in self.bonds])
        self.Emax = -_dmrg(neg.build_mpo(), dmrg_bonds, dmrg_tol)[0]
        # H + penalty (S^z_tot - 1)^2. The builder mishandles constants and same-site products, so the constant
        # part, penalty (N/4 + 1) from (S^z_i)^2 = 1/4 and the +1, is added by hand.
        trip = _builder(N, self.bonds)
        for i in range(N):
            for j in range(N):
                if i != j:
                    trip.add_term(penalty, ("sz", i), ("sz", j))
            trip.add_term(-2 * penalty, ("sz", i))
        self.E1 = _dmrg(trip.build_mpo(), dmrg_bonds, dmrg_tol)[0] + penalty * (N / 4 + 1)

    # ---- overlap knob
    def trial(self, eps):
        v = self.g0.copy()
        for i in range(self.N):
            s = eps * (-1) ** i
            v.gate_(np.diag([np.exp(-1j * s), np.exp(1j * s)]), i, contract=True)
        return v

    def overlap(self, eps):
        return abs(self.g0.overlap(self.trial(eps))) ** 2

    def eps_for_overlap(self, gamma):
        """Smallest eps with overlap = gamma (coarse scan to the first minimum, then brentq)."""
        es = np.linspace(0, 1.5, 151)
        g = np.array([self.overlap(e) for e in es])
        rise = np.flatnonzero(np.diff(g) > 0)
        stop = rise[0] if len(rise) else len(es) - 1
        if not g[stop] <= gamma <= 1.0:
            return None
        return brentq(lambda e: self.overlap(e) - gamma, 0.0, es[stop], xtol=1e-12)

    # ---- filter
    def _strang(self, v, dt, m):
        v = v.copy()
        info, gates = {}, {}                                   # info: quimb tracks the orthogonality centre
        for b, w in strang_sequence(len(self.bonds), m):
            i, j, J = self.bonds[b]
            th = w * (dt / 2) * J / self.W                     # exp(-i th S.S) = e^{i th/4}(cos(th/2) - i sin(th/2) SWAP)
            if th not in gates:
                gates[th] = np.exp(1j * th / 4) * (np.cos(th / 2) * np.eye(4) - 1j * np.sin(th / 2) * SWAP)
            v.gate_(gates[th], (i, j), contract="swap+split", max_bond=self.max_bond, cutoff=self.cutoff,
                    info=info)
        return v * np.exp(1j * m * dt * self.E0 / self.W)          # the -E0/W part of H_s, a global phase

    def pulse(self, v, dt, m):
        """(S^m + S^-m)/2 v: what one post-selected pulse leaves on the system."""
        out = (self._strang(v, dt, m) + self._strang(v, -dt, m)) * 0.5
        out.compress(max_bond=self.max_bond, cutoff=self.cutoff)
        return out

    def filter(self, psi, K, r, n, out=None):
        """Dyadic filter (pulses t = tau 2^k, k < K, ladder repeated r times, n Strang steps per tau).
        Returns (v, p, 1 - F). With `out` a list, (r, p, 1 - F) is appended after every ladder."""
        tau = np.pi / (1 + self.gap_s)
        v = psi
        for rr in range(r):
            for k in range(K):
                v = self.pulse(v, tau / n, n * 2 ** k)
            if out is not None:
                out.append((rr + 1,) + self.measure(v))
        return (v,) + self.measure(v)

    def measure(self, v):
        p = float(np.real(v.overlap(v)))
        amp = abs(self.g0.overlap(v)) ** 2
        return p, max(0.0, 1 - amp / p)


def delta_for_gap(N, gap, cfg, **kw):
    """Dimerisation for the target gap, from DMRG gaps (a few DMRG runs per brentq step)."""
    f = lambda d: MPSChain(N, d, cfg.J1, cfg.J2, **kw).gap - gap        # noqa: E731
    return brentq(f, *cfg.delta_range, xtol=1e-6)


def mps_chain(N, gap, cfg, delta=None):
    """MPSChain at the target gap. cfg.mps_reference: 'ed' (exact, small N), 'dmrg', or 'auto' (ED up to
    cfg.mps_ed_max_N). The dimerisation is solved with the same reference."""
    kw = dict(max_bond=cfg.mps_max_bond, cutoff=cfg.mps_cutoff, reference=cfg.mps_reference,
              ed_max_N=cfg.mps_ed_max_N, dmrg_bonds=cfg.mps_dmrg_bonds)
    use_ed = cfg.mps_reference == "ed" or (cfg.mps_reference == "auto" and N <= cfg.mps_ed_max_N)
    if delta is None:
        delta = (ed_delta_for_gap(N, gap, cfg.J1, cfg.J2, *cfg.delta_range) if use_ed
                 else delta_for_gap(N, gap, cfg, **kw))
    return MPSChain(N, delta, cfg.J1, cfg.J2, **kw)


def _table(mch, psi, n, cfg, e_min, bound):
    """(K, r) -> (p, 1 - F) at n steps per tau. Pulses commute, so at fixed r the filter for K + 1 is the one for
    K followed by r copies of pulse K: one pass over K per r. A pass stops once eps_min is met (more K only costs
    more), once 1 - F stops falling (saturated by Trotter error or too few repetitions), or once one attempt alone
    costs more than `bound`. The exact engine scans every K instead; `run_mps.py check` compares the two."""
    tau = np.pi / (1 + mch.gap_s)
    tab = {}
    for r in range(1, cfg.r_max + 1):
        if cnot_per_attempt(mch.N, 1, r, n, cfg.bond_cnots) > bound:
            break
        v = psi
        for K in range(1, cfg.K_max + 1):
            if cnot_per_attempt(mch.N, K, r, n, cfg.bond_cnots) > bound:
                break
            for _ in range(r):
                v = mch.pulse(v, tau / n, n * 2 ** (K - 1))
            tab[(K, r)] = mch.measure(v)
            if tab[(K, r)][1] <= e_min:
                break
            if K > 2 and tab[(K, r)][1] > 0.98 * tab[(K - 1, r)][1]:
                break                                   # 1 - F has saturated at this r and n: more K cannot help
    return tab


def measure_point(mch, gammas, eps_list, cfg, verbose=False):
    """Cheapest (K, r, n) per (gamma, eps, protocol) from MPS runs, with the rules of `search.measure_point`:
    for each K the fewest r that meets eps; a choice counts only if the same (K, r) also meets eps at the next n
    in cfg.n_grid; n is raised until a level brings no cheaper choice. Returns (rows, eps_knob), same format."""
    eps_knob = [mch.eps_for_overlap(g) for g in gammas]
    rows = [{e: {pr: None for pr in cfg.protocols} for e in eps_list} for _ in gammas]
    for jj, ek in enumerate(eps_knob):
        if ek is None:
            continue
        psi = mch.trial(ek)
        best = {}
        # Every passing filter has p <= gamma / (1 - eps) (the ground state is never amplified), so its RUS cost is
        # >= C (1 - eps)/gamma and its AA cost >= aa_cost(C, gamma / (1 - eps)): attempts above `bound` cannot win.
        p_up = min(1.0, gammas[jj] / (1 - max(eps_list)))
        for i, n in enumerate(cfg.n_grid[:-1]):
            bound = max([b["cost"] * p_up for (e, pr), b in best.items() if pr == "rus"] +
                        [b["cost"] / aa_cost(1.0, p_up) for (e, pr), b in best.items() if pr == "aa"],
                        default=np.inf) if len(best) == len(eps_list) * len(cfg.protocols) else np.inf
            tab = _table(mch, psi, n, cfg, min(eps_list), bound)
            improved = False
            for e in eps_list:
                cands = []
                for K in sorted({k for k, _ in tab}):
                    rs = sorted(r for k, r in tab if k == K and tab[(k, r)][1] <= e)
                    if rs:
                        r = rs[0]
                        p, inf = tab[(K, r)]
                        c = cnot_per_attempt(mch.N, K, r, n, cfg.bond_cnots)
                        refl = reflection_cnots(mch.N, K, r, cfg.mcz_cnots_per_qubit)
                        cands.append((K, r, p, inf, c, refl))
                for pr in cfg.protocols:
                    for K, r, p, inf, c, refl in sorted(cands, key=lambda x: cost(x[4], x[2], pr, x[5])):
                        cc = float(cost(c, p, pr, refl))
                        if (e, pr) in best and cc >= best[(e, pr)]["cost"]:
                            break
                        if mch.filter(psi, K, r, cfg.n_grid[i + 1])[2] > e:   # must hold at the next n too
                            continue
                        best[(e, pr)] = dict(K=K, r=r, n=n, cnot=c, p=p, infid=inf, cost=cc,
                                             trotter_steps=int(n * r * (2 ** K - 1)))
                        improved = True
                        break
            if verbose:
                print(f"    gamma {gammas[jj]}: n={n} " +
                      " ".join(f"{e:g}/{pr}:K{b['K']}r{b['r']}n{b['n']}" for (e, pr), b in best.items()), flush=True)
            if i > 0 and not improved and len(best) == len(eps_list) * len(cfg.protocols):
                break
        for (e, pr), b in best.items():
            rows[jj][e][pr] = b
    return rows, eps_knob


def truncation_check(mch, psi, K, r, n, factor=2):
    """Rerun one filter with factor x the bond dimension and a 100x smaller cutoff; returns both (p, 1 - F)."""
    a = mch.filter(psi, K, r, n)[1:]
    mb, co = mch.max_bond, mch.cutoff
    mch.max_bond, mch.cutoff = mb * factor, co / 100
    try:
        b = mch.filter(psi, K, r, n)[1:]
    finally:
        mch.max_bond, mch.cutoff = mb, co
    return a, b
