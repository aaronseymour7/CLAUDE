"""Faster Strang product formula for Heisenberg-type Hamiltonians, for the 20-qubit runs.

The stored term order puts the XX, YY and ZZ terms of each bond next to each other. Those three
commute, so applying them one after another is the same as applying the exact bond unitary
exp(-i theta (XX + YY + ZZ)) = e^{i theta} (cos 2 theta - i sin 2 theta SWAP) when the three
coefficients are equal. One SWAP permutation per bond replaces three Pauli rotations, and the result
equals sbpipe.trotter.product_formula(H, psi, t, k, 2) to rounding (checked in check()).
Bonds whose three coefficients differ, and any other term, fall back to the Pauli rotation.
"""
import numpy as np

from sbpipe.pauli import _action, apply_rotation


class Strang:
    def __init__(self, H):
        L, c = list(H.labels), np.asarray(H.coeffs)
        self.ops = []                      # ("bond", src, coeff) | ("pauli", label, coeff) | ("id", coeff)
        i = 0
        while i < len(L):
            lab = L[i]
            if set(lab) == {"I"}:
                self.ops.append(("id", None, complex(c[i])))
                i += 1
                continue
            sup = [q for q, p in enumerate(lab) if p != "I"]
            if (i + 2 < len(L) and len(sup) == 2 and set(lab) == {"I", "X"}
                    and L[i + 1] == lab.replace("X", "Y") and L[i + 2] == lab.replace("X", "Z")
                    and np.allclose(c[i:i + 3], c[i])):
                sx, _ = _action(lab)
                _, pz = _action(L[i + 2])
                idx = np.arange(sx.size)
                src = np.where(np.real(pz) > 0, idx, sx).astype(np.int64)
                self.ops.append(("bond", src, float(np.real(c[i]))))
                i += 3
            else:
                self.ops.append(("pauli", lab, complex(c[i])))
                i += 1

    def _apply(self, op, theta_scale, v):
        kind, a, cf = op
        th = cf * theta_scale
        if kind == "id":
            return np.exp(-1j * th) * v
        if kind == "pauli":
            return apply_rotation(a, th, v)
        return np.exp(1j * th) * (np.cos(2 * th) * v - 1j * np.sin(2 * th) * v[a])

    def evolve(self, psi, t, k):
        dt = t / k
        v = np.array(psi, dtype=complex, copy=True)
        buf = np.empty_like(v)
        seq = self.ops + self.ops[::-1]
        coef = []                                   # per op: (kind, src/label, a, b) precomputed
        for kind, a, cf in seq:
            th = cf * dt / 2
            if kind == "bond":
                coef.append((kind, a, np.exp(1j * th) * np.cos(2 * th), -1j * np.exp(1j * th) * np.sin(2 * th)))
            else:
                coef.append((kind, a, th, None))
        for _ in range(int(k)):
            for kind, a, x, y in coef:
                if kind == "bond":
                    np.take(v, a, out=buf)
                    buf *= y
                    v *= x
                    v += buf
                elif kind == "id":
                    v *= np.exp(-1j * x)
                else:
                    v = apply_rotation(a, x, v)
        return v


def filter_trotter(S, psi, times, phases, steps):
    v = np.array(psi, dtype=complex, copy=True)
    for t, ph, k in zip(times, phases, steps):
        a = S.evolve(v, t, int(k))
        b = S.evolve(v, -t, int(k))
        v = 0.5 * (np.exp(-1j * ph) * a + np.exp(1j * ph) * b)
    return v


def qetu_trotter(S, psi, d, tau, x_s, steps):
    from sbpipe.simulate import _chebyshev_apply
    it = iter(np.asarray(steps, int))

    def X(v):
        k = int(next(it))
        return 0.5 * (S.evolve(v, tau, k) + S.evolve(v, -tau, k))
    return _chebyshev_apply(X, psi, d, x_s)


def check():
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from run_extended import model
    from sbpipe.trotter import product_formula
    rng = np.random.default_rng(1)
    for kind, size in (("chain", 10), ("cyl", 2)):
        H = model(kind, size, 0.3).shifted_scaled(-5.0, 7.0)
        v = rng.normal(size=2 ** H.num_qubits) + 1j * rng.normal(size=2 ** H.num_qubits)
        a = product_formula(H, v, 1.3, 3, 2)
        b = Strang(H).evolve(v, 1.3, 3)
        print(kind, size, np.abs(a - b).max(), sum(o[0] == "bond" for o in Strang(H).ops))


if __name__ == "__main__":
    check()
