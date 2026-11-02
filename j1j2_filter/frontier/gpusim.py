"""Statevector Strang product formula for Heisenberg-bond Hamiltonians, on a GPU (CuPy, ROCm or CUDA) or on the
CPU (NumPy), for the SBC / QETU filter circuits of sbpipe at 24-30 qubits.

Every term of the J1-J2 Hamiltonian is a bond c (XX + YY + ZZ) = c (2 SWAP - 1), and the three Paulis of a bond
are consecutive in sbpipe's term (= Trotter) order and commute, so applying them one after another equals the exact
bond unitary

    exp(-i theta (XX + YY + ZZ)) = e^{-i theta} [ 1 on |00>, |11> ;  e^{2 i theta} (cos 2theta - i sin 2theta SWAP)
                                                  on the |01>, |10> block ]

This is the same circuit as sbpipe.trotter.product_formula(H, psi, t, k, order=2), to rounding (selftest()).
No index arrays: each bond is two strided views of the state, so memory is the state vectors only.

Memory per state vector: 16 * 2^N bytes (complex128): 4 GiB at N = 28, 16 GiB at N = 30. filter_trotter keeps
2 vectors plus a quarter-size temporary; qetu_trotter keeps 4. One Frontier GCD has 64 GB of HBM.
"""
import os

import numpy as np

_FORCE_CPU = os.environ.get("SBPIPE_CPU", "0") == "1"
try:
    if _FORCE_CPU:
        raise ImportError
    import cupy as xp
    GPU = True
except ImportError:                                     # CPU fallback (tests, small N)
    xp = np
    GPU = False

_mix = None
if GPU:
    _mix = xp.ElementwiseKernel(
        "complex128 x, complex128 y", "complex128 a, complex128 b",
        "complex<double> ta = a; complex<double> tb = b; a = x * ta + y * tb; b = x * tb + y * ta;",
        "sbpipe_bond_mix")


def to_host(v):
    return xp.asnumpy(v) if GPU else np.asarray(v)


def to_dev(v):
    return xp.asarray(np.asarray(v, dtype=np.complex128))


class BondStrang:
    """Strang splitting of H (a sbpipe PauliSum, possibly shifted and scaled) in its stored term order."""

    def __init__(self, H, fused=True):
        L, c = list(H.labels), np.asarray(H.coeffs)
        self.N = H.num_qubits
        self.ops, self.const = [], 0.0          # ops: (i, j, coeff) per bond; const: identity coefficient
        i = 0
        while i < len(L):
            lab = L[i]
            if set(lab) == {"I"}:
                self.const += float(np.real(c[i]))
                i += 1
                continue
            sup = [q for q, p in enumerate(lab) if p != "I"]
            ok = (i + 2 < len(L) and len(sup) == 2 and set(lab) == {"I", "X"}
                  and L[i + 1] == lab.replace("X", "Y") and L[i + 2] == lab.replace("X", "Z")
                  and np.allclose(c[i:i + 3], c[i]))
            if not ok:
                raise ValueError(f"term {lab} is not part of an XX+YY+ZZ bond; only Heisenberg bonds are supported")
            self.ops.append((sup[0], sup[1], float(np.real(c[i]))))
            i += 3
        self.fused = fused and GPU and os.environ.get("SBPIPE_NOFUSE", "0") != "1"

    def _views(self, v, i, j):
        N = self.N
        V = v.reshape(2 ** i, 2, 2 ** (j - i - 1), 2, 2 ** (N - j - 1))
        return V[:, 0, :, 1, :], V[:, 1, :, 0, :]

    def _bond(self, v, i, j, x, y):
        a, b = self._views(v, i, j)
        if self.fused:
            _mix(x, y, a, b)
        else:
            t = a.copy()
            a *= x
            a += y * b
            b *= x
            b += y * t

    def evolve_(self, v, t, k):
        """In place: v <- S(t/k)^k v, S the second-order (Strang) step."""
        k = int(k)
        dt = t / k
        seq = self.ops + self.ops[::-1]
        coef, phase = [], 0.0
        for i, j, cf in seq:
            th = cf * dt / 2
            phase += -th                                   # e^{-i theta} on the whole state
            e2 = np.exp(2j * th)
            coef.append((i, j, complex(e2 * np.cos(2 * th)), complex(-1j * e2 * np.sin(2 * th))))
        phase += -self.const * dt                          # identity term, both halves
        for _ in range(k):
            for i, j, x, y in coef:
                self._bond(v, i, j, x, y)
        v *= np.exp(1j * phase * k)
        return v

    def energy(self, v):
        """<v|H|v> / <v|v> for the Hamiltonian this object was built from."""
        n = float(xp.vdot(v, v).real)
        e = self.const * n
        for i, j, cf in self.ops:
            a, b = self._views(v, i, j)
            sw = n - float(xp.vdot(a, a).real) - float(xp.vdot(b, b).real) + 2 * float(xp.vdot(a, b).real)
            e += cf * (2 * sw - n)
        return e / n


def filter_trotter(S, psi, times, phases, steps):
    """SBC filter, success branch: prod_i [e^{-i phi_i} U(t_i) + e^{i phi_i} U(-t_i)] / 2 on psi (unnormalised)."""
    v = xp.array(psi, dtype=xp.complex128, copy=True)
    for t, ph, k in zip(times, phases, steps):
        a = S.evolve_(v.copy(), t, int(k))
        S.evolve_(v, -t, int(k))                           # v is now U(-t) psi
        v *= 0.5 * np.exp(1j * ph)
        v += (0.5 * np.exp(-1j * ph)) * a
        del a
    return v


def qetu_trotter(S, psi, d, tau, x_s, steps):
    """Ancilla-0 block of the Trotterised QETU circuit: T_d(X/x_s) psi / T_d(1/x_s), X = (S(tau) + S(-tau)) / 2."""
    from sbpipe.design import _cheb
    it = iter(np.asarray(steps, int))

    def X(v):
        k = int(next(it))
        a = S.evolve_(v.copy(), tau, k)
        b = S.evolve_(v.copy(), -tau, k)
        a += b
        a *= 0.5
        return a

    v0 = xp.array(psi, dtype=xp.complex128, copy=True)
    if d == 0:
        return v0
    v1 = X(v0) / x_s
    for _ in range(d - 1):
        w = X(v1)
        w *= 2.0 / x_s
        w -= v0
        v0, v1 = v1, w
    return v1 / _cheb(d, np.array([1.0 / x_s]))[0]


def selftest():
    """Compare with sbpipe's term-by-term product formula on a 4x3 cylinder (12 qubits), J2 = 0.3."""
    from sbpipe.models import j1j2_square
    from sbpipe.trotter import product_formula
    H = j1j2_square(3, 4, 1.0, 0.3).shifted_scaled(-5.0, 7.0)
    rng = np.random.default_rng(1)
    v = rng.normal(size=2 ** H.num_qubits) + 1j * rng.normal(size=2 ** H.num_qubits)
    a = product_formula(H, v, 1.3, 3, 2)
    out = {}
    for fused in ((False, True) if GPU else (False,)):
        S = BondStrang(H, fused=fused)
        b = to_host(S.evolve_(to_dev(v), 1.3, 3))
        out[f"evolve fused={fused}"] = float(np.abs(a - b).max())
    S = BondStrang(H)
    E_ref = float(np.vdot(v, H.to_sparse() @ v).real / np.vdot(v, v).real)
    out["energy"] = abs(S.energy(to_dev(v)) - E_ref)
    out["backend"] = "GPU (CuPy)" if GPU else "CPU (NumPy)"
    return out


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    print(selftest())
