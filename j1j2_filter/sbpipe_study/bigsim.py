"""Matrix-free state-vector tools for Heisenberg-type J1-J2 Hamiltonians at 16-24 qubits.

Every bond of the model is J/4 (XX + YY + ZZ) = J/4 (2 SWAP - 1), so H|v> needs only axis swaps of the
state tensor (no sparse matrix, no index arrays), and the exact bond unitary is
exp(-i th (XX + YY + ZZ)) = e^{i th} (cos 2th - i sin 2th SWAP). The Strang product applies the bonds in
the PauliSum's stored order and then in reverse, so it equals sbpipe.trotter.product_formula(H, ., ., ., 2)
and fastsim.Strang to rounding (bigsim.check()). Exact evolution uses a Chebyshev expansion on the known
spectral interval.
"""
import math

import numpy as np
from scipy.special import jv


class Heis:
    def __init__(self, H):
        L, c = list(H.labels), np.asarray(H.coeffs, float)
        self.N = len(L[0])
        self.bonds, self.const = [], 0.0
        i = 0
        while i < len(L):
            lab = L[i]
            if set(lab) == {"I"}:
                self.const += c[i]
                i += 1
                continue
            sup = [q for q, p in enumerate(lab) if p != "I"]
            assert (len(sup) == 2 and L[i + 1] == lab.replace("X", "Y") and L[i + 2] == lab.replace("X", "Z")
                    and np.allclose(c[i:i + 3], c[i])), "not a Heisenberg bond triple"
            self.bonds.append((sup[0], sup[1], float(c[i])))      # coefficient of XX (= of YY, ZZ)
            i += 3

    def _view(self, v, i, j):
        """v as (2^i, 2, 2^(j-i-1), 2, 2^(N-j-1)) for qubits i < j (qubit 0 is the most significant bit)."""
        if i > j:
            i, j = j, i
        return v.reshape(2 ** i, 2, 2 ** (j - i - 1), 2, 2 ** (self.N - j - 1))

    def swap(self, v, i, j):
        out = np.array(v, copy=True)
        a, b = self._view(out, i, j), self._view(v, i, j)
        a[:, 0, :, 1, :] = b[:, 1, :, 0, :]
        a[:, 1, :, 0, :] = b[:, 0, :, 1, :]
        return out

    def matvec(self, v):
        """sum_b 2c_b SWAP_b + (const - sum_b c_b), with SWAP acting only on the 01/10 blocks."""
        out = (self.const + sum(c for _, _, c in self.bonds)) * v     # 2c SWAP - c = c on 00/11 ...
        o = out
        for i, j, c in self.bonds:
            a, b = self._view(o, i, j), self._view(v, i, j)
            a[:, 0, :, 1, :] += (2 * c) * b[:, 1, :, 0, :] - (2 * c) * b[:, 0, :, 1, :]
            a[:, 1, :, 0, :] += (2 * c) * b[:, 0, :, 1, :] - (2 * c) * b[:, 1, :, 0, :]
        return out

    def bond(self, v, i, j, th):
        """exp(-i th (XX + YY + ZZ)) on qubits i, j, in place: e^{-i th} on |00>, |11>, and
        e^{i th} (cos 2th - i sin 2th SWAP) on the |01>, |10> pair."""
        x = self._view(v, i, j)
        e = np.exp(-1j * th)
        x[:, 0, :, 0, :] *= e
        x[:, 1, :, 1, :] *= e
        a = x[:, 0, :, 1, :].copy()
        b = x[:, 1, :, 0, :]
        cc, ss = np.exp(1j * th) * math.cos(2 * th), -1j * np.exp(1j * th) * math.sin(2 * th)
        x[:, 0, :, 1, :] = cc * a + ss * b
        x[:, 1, :, 0, :] = cc * b + ss * a
        return v

    def strang(self, v, t, k):
        dt = t / k
        v = np.array(v, dtype=complex, copy=True)
        seq = self.bonds + self.bonds[::-1]
        ph = np.exp(-1j * self.const * dt)
        for _ in range(int(k)):
            for i, j, c in seq:
                v = self.bond(v, i, j, c * dt / 2)
            v *= ph
        return v

    def expm(self, v, t, lo, hi, tol=1e-13):
        """exp(-i t H) v by a Chebyshev series on [lo, hi] (must contain the spectrum)."""
        if t == 0:
            return np.array(v, dtype=complex, copy=True)
        a, b = 0.5 * (hi - lo), 0.5 * (hi + lo)
        z = a * abs(t)
        n = int(1.3 * z + 40)
        while abs(jv(n, z)) > tol:
            n += 10
        s = -1j * np.sign(t)
        mv = lambda x: (self.matvec(x) - b * x) / a
        T0 = np.array(v, dtype=complex, copy=True)
        T1 = mv(T0)
        out = jv(0, z) * T0 + 2 * s * jv(1, z) * T1
        for m in range(2, n + 1):
            T0, T1 = T1, 2 * mv(T1) - T0
            out += 2 * (s ** m) * jv(m, z) * T1
        return np.exp(-1j * t * b) * out


def rotate(psi, N, theta):
    """exp(-i theta Y/2) on the even qubits (lowers the singlet overlap smoothly)."""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    v = np.array(psi, copy=True).reshape([2] * N)
    for q in range(0, N, 2):
        v = np.moveaxis(v, q, 0)
        a, b = v[0].copy(), v[1].copy()
        v = np.moveaxis(np.stack([c * a - s * b, s * a + c * b]), 0, q)
    return np.ascontiguousarray(v).reshape(-1)


def lowest(Hh, k=6, tol=1e-11):
    """Lowest k eigenpairs of the (real) Hamiltonian, matrix-free Lanczos."""
    from scipy.sparse.linalg import LinearOperator, eigsh
    n = 2 ** Hh.N
    op = LinearOperator((n, n), matvec=lambda x: Hh.matvec(x.astype(float)).real, dtype=float)
    rng = np.random.default_rng(0)
    w, V = eigsh(op, k=k, which="SA", tol=tol, v0=rng.normal(size=n))
    o = np.argsort(w)
    return w[o], V[:, o]


def check():
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from sbpipe import j1j2_chain
    from sbpipe.models import j1j2_square
    from sbpipe.trotter import product_formula
    import scipy.sparse.linalg as sla
    rng = np.random.default_rng(1)
    for H in (j1j2_chain(10, 1.0, 0.4).shifted_scaled(-5.0, 7.0), j1j2_square(2, 4, 1.0, 0.3).shifted_scaled(-5, 7)):
        h = Heis(H)
        v = rng.normal(size=2 ** h.N) + 1j * rng.normal(size=2 ** h.N)
        Hsp = H.to_sparse()
        print("matvec", np.abs(h.matvec(v) - Hsp @ v).max(),
              "strang", np.abs(h.strang(v, 1.3, 3) - product_formula(H, v, 1.3, 3, 2)).max(),
              "expm", np.abs(h.expm(v, 7.0, -1.0, 1.5) - sla.expm_multiply(-1j * 7.0 * Hsp, v)).max())


if __name__ == "__main__":
    check()
