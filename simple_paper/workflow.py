"""
workflow.py -- minimal, self-contained version of the DMRG-trial-state + Stetcu-Baroni filter workflow.

Inputs  : a Hamiltonian (open J1-J2 chain, grouped into bond terms), energies E0, E1, Etop,
          a trial state |psi>, and a filter F(E) = prod_i cos(E t_i + phi_i).
Outputs : an a-priori bound on the trace distance D(final state, |E0>) and the measured value.

Only numpy/scipy are needed.  Filter design and certification are taken from ../j1j2_filter/floor.py.
Everything is dense-vector simulation (N <= 12).
"""
import sys, pathlib
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "j1j2_filter"))
import floor as FL  # noqa: E402  (filter design + rigorous eta certificate)

# ----------------------------------------------------------------------------- Hamiltonian
X = sp.csr_matrix([[0, 1], [1, 0]], dtype=complex)
Y = sp.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
Z = sp.csr_matrix([[1, 0], [0, -1]], dtype=complex)
I2 = sp.identity(2, format="csr", dtype=complex)


def bonds(N, J1=1.0, J2=0.0):
    """Bond list [(i, j, J)] in Trotter order: nn even, nn odd, nnn even, nnn odd (i = left site)."""
    out = []
    for par in (0, 1):
        out += [(i, i + 1, J1) for i in range(par, N - 1, 2)]
    if J2 != 0.0:
        for par in (0, 1):
            out += [(i, i + 2, J2) for i in range(par, N - 2, 2)]
    return out


def _two_site(P, i, j, N):
    ops = [I2] * N
    ops[i] = ops[j] = P
    M = ops[0]
    for o in ops[1:]:
        M = sp.kron(M, o, format="csr")
    return M


def bond_matrix(i, j, N):
    """S_i . S_j = (XX+YY+ZZ)/4 as a sparse 2^N matrix."""
    return 0.25 * (_two_site(X, i, j, N) + _two_site(Y, i, j, N) + _two_site(Z, i, j, N))


def spectrum_and_ground_state(Hs):
    """Exact E0, E1, Etop and ground state (N <= 12)."""
    d = Hs.shape[0]
    if d <= 1024:
        w, v = np.linalg.eigh(Hs.toarray())
        return w[0], w[np.searchsorted(w, w[0] + 1e-9)], w[-1], v[:, 0], w
    v0 = np.random.default_rng(0).standard_normal(d)       # fixed start: reproducible gamma
    w, v = sla.eigsh(Hs, k=3, which="SA", v0=v0)
    wt = sla.eigsh(Hs, k=1, which="LA", return_eigenvectors=False, v0=v0)[0]
    o = np.argsort(w)
    w, v = w[o], v[:, o]
    return w[0], w[np.searchsorted(w, w[0] + 1e-9)], wt, v[:, 0], None


class Chain:
    """Scaled Hamiltonian Hs = (H - E0)/W = sum_g H_g + c0 with H_g = (J_g/W) S_i.S_j and c0 = -E0/W."""

    def __init__(self, N, J2=0.0, J1=1.0):
        self.N, self.J2 = N, J2
        self.bonds = bonds(N, J1, J2)
        H = sum(Jg * bond_matrix(i, j, N) for i, j, Jg in self.bonds)
        self.E0, self.E1, self.Etop, g, self.allE = spectrum_and_ground_state(H.tocsr())
        g = g / np.linalg.norm(g)
        self.g = g
        self.W = self.Etop - self.E0
        self.gap = (self.E1 - self.E0) / self.W            # scaled gap Delta
        self.c0 = -self.E0 / self.W
        self.coef = [Jg / self.W for _, _, Jg in self.bonds]  # H_g = coef_g * S.S
        self.Hg = [c * bond_matrix(i, j, N) for (i, j, _), c in zip(self.bonds, self.coef)]
        self.Hs = sum(self.Hg) + self.c0 * sp.identity(2 ** N, format="csr")
        self.H = H.tocsr()
        self._alpha = None

    # --- Trotter commutator constant (first order), max over forward and reversed order
    @property
    def alpha(self):
        if self._alpha is None:
            def opnorm(A):
                if A.shape[0] <= 1024:
                    return np.linalg.norm(A.toarray(), 2)
                if A.nnz == 0 or abs(A).max() < 1e-14:
                    return 0.0                      # commuting pair
                return sla.svds(A, k=1, return_singular_vectors=False,
                                v0=np.random.default_rng(0).standard_normal(A.shape[0]))[0]
            a = []
            for order in (self.Hg, self.Hg[::-1]):
                tot, tail = 0.0, None
                # sum_g || [H_g, sum_{g'>g} H_g'] ||
                suffix = [None] * len(order)
                acc = None
                for idx in range(len(order) - 1, -1, -1):
                    suffix[idx] = acc
                    acc = order[idx] if acc is None else acc + order[idx]
                for idx, Hg in enumerate(order):
                    if suffix[idx] is None:
                        continue
                    C = Hg @ suffix[idx] - suffix[idx] @ Hg
                    tot += opnorm(C.tocsr())
                a.append(tot)
            self._alpha = max(a)
        return self._alpha


# ----------------------------------------------------------------------------- trial state
def mps_truncate(psi, N, chi):
    """TT-SVD truncation of a state vector to bond dimension chi (stand-in for a bond-dimension-chi DMRG state)."""
    M, left, tensors = psi.reshape(1, -1), 1, []
    for _ in range(N - 1):
        U, s, Vh = np.linalg.svd(M.reshape(left * 2, -1), full_matrices=False)
        r = min(chi, len(s))
        tensors.append(U[:, :r].reshape(left, 2, r))
        M, left = s[:r, None] * Vh[:r], r
    tensors.append(M.reshape(left, 2, 1))
    out = tensors[0].reshape(2, -1)
    for T_ in tensors[1:]:
        out = (out @ T_.reshape(T_.shape[0], -1)).reshape(-1, T_.shape[2])
    out = out.reshape(-1)
    return out / np.linalg.norm(out)


# ----------------------------------------------------------------------------- pulses
def _swap_gate(theta_c):
    """exp(-i * theta_c * S.S) on two qubits:  S.S = 1/4 on the triplet, -3/4 on the singlet."""
    P_s = np.zeros((4, 4), complex)
    P_s[0, 0] = P_s[3, 3] = 1.0
    P_s[1:3, 1:3] = 0.5
    P_a = np.eye(4) - P_s
    return np.exp(-1j * theta_c / 4) * P_s + np.exp(1j * 3 * theta_c / 4) * P_a


def apply_bond(psi, gate, i, j, N):
    t = psi.reshape((2,) * N)
    t = np.moveaxis(t, (i, j), (0, 1)).reshape(4, -1)
    t = (gate @ t).reshape((2, 2) + (2,) * (N - 2))
    return np.moveaxis(t, (0, 1), (i, j)).reshape(-1)


def trotter_evolve(ch, psi, tau, k, sign=+1):
    """k first-order steps of exp(-i sign*tau*sum_g H_g), step dt = tau/k (the constant c0 is applied separately)."""
    if k == 0:
        return psi
    dt = tau / k
    gates = [_swap_gate(sign * dt * c) for c in ch.coef]
    for _ in range(k):
        for (i, j, _J), g in zip(ch.bonds, gates):
            psi = apply_bond(psi, g, i, j, ch.N)
    return psi


def pulse_trotter(ch, psi, t, phi, k):
    """Post-selected (ancilla = 0) block  <0|W|0>  with k Trotter steps: 1/2 (e^{-i(phi+t c0)} V+ + e^{+i(phi+t c0)} V-)."""
    vp = trotter_evolve(ch, psi, t, k, +1)
    vm = trotter_evolve(ch, psi, t, k, -1)
    return 0.5 * (np.exp(-1j * (phi + t * ch.c0)) * vp + np.exp(1j * (phi + t * ch.c0)) * vm)


def pulse_exact(ch, psi, t, phi):
    A = sla.expm_multiply(-1j * t * ch.Hs, psi)
    B = sla.expm_multiply(+1j * t * ch.Hs, psi)
    return 0.5 * (np.exp(-1j * phi) * A + np.exp(1j * phi) * B)


def run_filter_exact(ch, psi, times, phases):
    for t, p in zip(times, phases):
        psi = pulse_exact(ch, psi, t, p)
    return psi


def run_filter_trotter(ch, psi, ks, dt, phases):
    for k, p in zip(ks, phases):
        psi = pulse_trotter(ch, psi, k * dt, p, int(k))
    return psi


def dist(a, b):
    """Trace distance between the pure states |a>, |b> (= sin of the angle) = sqrt(1 - |<a|b>|^2)."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.sqrt(max(0.0, 1.0 - abs(np.vdot(a, b)) ** 2)))


# ----------------------------------------------------------------------------- bound
def leakage_bound(gamma, eta):
    """l <= (1-gamma) eta^2 / (gamma + (1-gamma) eta^2)."""
    return (1 - gamma) * eta ** 2 / (gamma + (1 - gamma) * eta ** 2)


def trace_distance_bound(gamma, eta, f0sq, alpha, T, n):
    """D(final, E0) <= sqrt(l) + eps_T / sqrt(p_g),  eps_T = alpha T^2 / (2 n),  p_g = gamma F(0)^2."""
    dl = float(np.sqrt(leakage_bound(gamma, eta)))
    pg = gamma * f0sq
    dt = min(1.0, alpha * T ** 2 / (2 * n) / np.sqrt(pg))
    return min(1.0, dl + dt), dl, dt
