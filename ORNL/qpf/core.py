"""The QPF step and its repetition, exactly as in Stetcu et al.

One step: attach an ancilla a in |0>, apply exp(-i (t O + delta) (x) Y_a), measure a.

    exp(-i theta Y)|0> = cos(theta)|0> + sin(theta)|1>, so
    |psi>|0>  ->  cos(t O + delta)|psi> |0>  +  sin(t O + delta)|psi> |1>          (paper Eq. 1)

Reading 0 keeps cos(t O + delta)|psi> (unnormalised), with probability <psi| cos^2(t O + delta) |psi>.
A sequence of steps that all read 0 leaves N prod_i cos(t_i O + delta_i) |psi(0)>, with success probability equal
to the product of the step probabilities. The filter function is f(o) = prod_i cos(t_i o + delta_i) on the
eigenvalues o of O.

O is any Hermitian operator: a scipy sparse matrix, a dense array, or a `Diagonal` (O given in its eigenbasis,
as the paper does for design: "all states are written as superposition of H~ eigenstates").
Evolution is exact (no Trotter): cos(A) v = (e^{-iA} v + e^{iA} v) / 2 via scipy's expm_multiply.
"""
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply


@dataclass
class Diagonal:
    """An operator given by its eigenvalues; vectors are amplitudes in its eigenbasis."""
    eigvals: np.ndarray

    def cos(self, t, delta, v):
        return np.cos(t * self.eigvals + delta) * v

    def sin(self, t, delta, v):
        return np.sin(t * self.eigvals + delta) * v

    def expect(self, v):
        return float(np.real(np.vdot(v, self.eigvals * v)) / np.real(np.vdot(v, v)))


class Operator:
    """A Hermitian operator O (sparse or dense matrix) with exact cos/sin(t O + delta) actions on vectors."""

    def __init__(self, O):
        self.O = sp.csr_matrix(O) if not sp.issparse(O) else O.tocsr()

    def _exp(self, s, v):                                   # e^{-i s O} v
        if s == 0:
            return np.array(v, dtype=complex)
        return expm_multiply(-1j * s * self.O, np.asarray(v, dtype=complex))

    def cos(self, t, delta, v):
        a, b = self._exp(t, v), self._exp(-t, v)
        return 0.5 * (np.exp(-1j * delta) * a + np.exp(1j * delta) * b)

    def sin(self, t, delta, v):
        a, b = self._exp(t, v), self._exp(-t, v)
        return (np.exp(1j * delta) * b - np.exp(-1j * delta) * a) / 2j

    def expect(self, v):
        return float(np.real(np.vdot(v, self.O @ v)) / np.real(np.vdot(v, v)))


def as_operator(O):
    return O if isinstance(O, (Diagonal, Operator)) else Operator(O)


def ancilla_step(O, psi, t, delta=0.0):
    """Joint state after exp(-i (t O + delta) (x) Y_a) on |psi>|0>: returns (branch a=0, branch a=1)."""
    O = as_operator(O)
    return O.cos(t, delta, psi), O.sin(t, delta, psi)


@dataclass
class Run:
    """Result of a sequence of post-selected steps."""
    state: np.ndarray                 # normalised state after the last successful step
    p_steps: list                     # probability of reading 0 at each step (given all earlier steps read 0)
    p_total: float                    # product of p_steps: success probability of the whole sequence
    times: list = field(default_factory=list)
    phases: list = field(default_factory=list)
    trajectory: list = field(default_factory=list)   # normalised state after each step (if requested)


def run(O, psi, times, phases=None, keep=False):
    """Apply the steps (t_i, delta_i) in order, post-selecting 0 each time (paper: Eq. 1 iterated)."""
    O = as_operator(O)
    phases = np.zeros(len(times)) if phases is None else np.asarray(phases, float)
    v = np.asarray(psi, dtype=complex)
    v = v / np.linalg.norm(v)
    ps, traj = [], []
    for t, d in zip(times, phases):
        w = O.cos(t, d, v)
        p = float(np.real(np.vdot(w, w)))
        ps.append(p)
        if p <= 0:
            raise FloatingPointError("post-selection probability is zero")
        v = w / np.sqrt(p)
        if keep:
            traj.append(v.copy())
    return Run(v, ps, float(np.prod(ps)), list(times), list(phases), traj)


def filter_function(o, times, phases=None):
    """f(o) = prod_i cos(t_i o + delta_i) (the paper writes the product over steps, evaluated at eigenvalues o)."""
    o = np.asarray(o, float)
    phases = np.zeros(len(times)) if phases is None else np.asarray(phases, float)
    f = np.ones_like(o)
    for t, d in zip(times, phases):
        f = f * np.cos(t * o + d)
    return f


def success_probability_limit(overlap_sq, phases):
    """Converged success probability P = f(0)^2 |<Psi0|psi(0)>|^2, f(0) = prod_i cos(delta_i)."""
    return float(np.prod(np.cos(phases)) ** 2 * overlap_sq)


def scale(H, E0, Emax):
    """The paper's dimensionless Hamiltonian: shifted and scaled so the spectrum runs from 0 (target) to 1."""
    H = sp.csr_matrix(H) if not sp.issparse(H) else H
    return (H - E0 * sp.identity(H.shape[0], format="csr")) / (Emax - E0)
