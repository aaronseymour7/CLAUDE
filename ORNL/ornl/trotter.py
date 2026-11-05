"""Trotterised evolution under H_s = (H - E0) / W, second-order Strang, in the S^z = 0 sector.

Filter pulse (one ancilla, post-selected on 0): (S^m + S^-m) / 2, with S one Strang step of length dt.

Two equivalent ways to apply it:
  strang(ch, V, dt, m)   gate by gate: each bond factor is exp(-i b SWAP) = cos b - i sin b SWAP, one gather.
  step_spectrum(ch, dt)  exact eigenbasis of S, so that (S^m + S^-m)/2 = O cos(m phi) O^T for every m at once.
                         S = A_1..A_B A_B..A_1 is complex symmetric (S = P P^T) and unitary, so O is real orthogonal;
                         it comes from a real symmetric eigensolver on Re S + alpha Im S. Spin flip commutes with
                         every SWAP, which splits S into two blocks of half the size.
"""
import numpy as np


def _half_factors(ch, dt):
    """(cos b, -i sin b, perm) per bond for a half step dt/2: exp(-i (dt/2)(J/W) S.S) = phase * exp(-i b SWAP),
    b = J dt / (4 W)."""
    return [(np.cos(0.5 * J * dt / (2 * ch.W)), -1j * np.sin(0.5 * J * dt / (2 * ch.W)), p)
            for (_, _, J), p in zip(ch.bonds, ch.perm)]


def _phase(ch, dt):
    """The identity part of H_s per step, collected into one global phase."""
    return np.exp(-1j * (ch.c_id - ch.E0) / ch.W * dt)


def strang(ch, V, dt, steps):
    """`steps` Strang steps of e^{-i dt H_s} applied to V (dim or dim x m), gate by gate."""
    seq = _half_factors(ch, dt)
    seq = seq + seq[::-1]
    V = np.array(V, dtype=complex)
    for _ in range(steps):
        for a, b, p in seq:
            V = a * V + b * V[p]
    return _phase(ch, dt) ** steps * V


def pulse(ch, V, dt, m):
    """One Trotterised filter pulse, (S^m + S^-m) / 2, gate by gate."""
    return 0.5 * (strang(ch, V, dt, m) + strang(ch, V, -dt, m))


def step_spectrum(ch, dt):
    """Eigenphases phi of one Strang step S(dt) = O e^{-i phi} O^T, and a function writing vectors in that basis.
    Returns (phi, project): project(v) = O^T v for v of shape (dim,) or (dim, m)."""
    d, N = ch.dim, ch.N
    flip = ch.index_of(((1 << N) - 1) ^ ch.states)
    a = np.flatnonzero(np.arange(d) < flip)                       # one state of each spin-flip pair
    b = flip[a]
    half = _half_factors(ch, dt)
    phis, Os = [], []
    for sign in (1.0, -1.0):
        G = np.zeros((d, len(a)), dtype=complex)                   # columns: spin-flip (anti)symmetric states
        G[a, np.arange(len(a))] = 1 / np.sqrt(2)
        G[b, np.arange(len(a))] = sign / np.sqrt(2)
        for x, y, p in half:                                       # G = P^T Q
            Gp = G[p]
            G *= x
            Gp *= y
            G += Gp
            del Gp
        S = (G.T @ G) * _phase(ch, dt)                             # S restricted to the sector
        del G
        _, O = np.linalg.eigh(S.real + 0.1234 * S.imag)
        lam = np.einsum("ij,ij->j", O, S.real @ O) + 1j * np.einsum("ij,ij->j", O, S.imag @ O)
        phis.append(-np.angle(lam))
        Os.append((sign, O))

    def project(v):
        return np.concatenate([O.T @ ((v[a] + sign * v[b]) / np.sqrt(2)) for sign, O in Os], axis=0)
    return np.concatenate(phis), project


def dyadic_filter(ch, psi, K, r, n):
    """Apply the dyadic filter (pulses tau 2^k, k < K, ladder repeated r times, n Strang steps per tau) gate by gate.
    Returns the unnormalised post-selected state, its success probability p and infidelity 1 - F."""
    tau = np.pi / (1 + ch.gap_s)
    v = np.array(psi, dtype=complex)
    for _ in range(r):
        for k in range(K):
            v = pulse(ch, v, tau / n, n * 2 ** k)
    p = float(np.vdot(v, v).real)
    return v, p, 1 - abs(np.vdot(ch.g0, v)) ** 2 / p
