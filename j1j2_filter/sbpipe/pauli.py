"""Minimal Pauli-sum algebra and statevector action (numpy only).

Conventions
-----------
* A Pauli string is a label of length N over {I, X, Y, Z}; label[q] acts on qubit q.
* Statevectors are big-endian: qubit q is bit (N-1-q) of the basis index, so qubit 0
  is the most significant bit. This matches quimb's ``MatrixProductState.to_dense()``
  with site q = qubit q.
* A ``PauliSum`` keeps its terms in a fixed order; that order is the Trotter order.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import scipy.sparse as sp

_MUL = {  # single-qubit products a*b = phase * c
    ("I", "I"): (1, "I"), ("I", "X"): (1, "X"), ("I", "Y"): (1, "Y"), ("I", "Z"): (1, "Z"),
    ("X", "I"): (1, "X"), ("X", "X"): (1, "I"), ("X", "Y"): (1j, "Z"), ("X", "Z"): (-1j, "Y"),
    ("Y", "I"): (1, "Y"), ("Y", "X"): (-1j, "Z"), ("Y", "Y"): (1, "I"), ("Y", "Z"): (1j, "X"),
    ("Z", "I"): (1, "Z"), ("Z", "X"): (1j, "Y"), ("Z", "Y"): (-1j, "X"), ("Z", "Z"): (1, "I"),
}


def pauli_mul(a: str, b: str) -> Tuple[complex, str]:
    """Product of two Pauli strings: a*b = phase * label."""
    ph, out = 1 + 0j, []
    for x, y in zip(a, b):
        p, c = _MUL[(x, y)]
        ph *= p
        out.append(c)
    return ph, "".join(out)


def commutes(a: str, b: str) -> bool:
    n = sum(1 for x, y in zip(a, b) if x != "I" and y != "I" and x != y)
    return n % 2 == 0


def support(label: str) -> List[int]:
    return [q for q, c in enumerate(label) if c != "I"]


def weight(label: str) -> int:
    return len(label) - label.count("I")


@dataclass
class PauliSum:
    """H = sum_k coeffs[k] * P_k with real coefficients (Hermitian)."""

    labels: List[str]
    coeffs: np.ndarray
    groups: List[int] = field(default=None)   # optional group id per term (e.g. bond), for Weyl bounds

    def __post_init__(self):
        self.coeffs = np.asarray(self.coeffs, dtype=float)
        if len(self.labels) != len(self.coeffs):
            raise ValueError("labels and coeffs differ in length")
        if len({len(l) for l in self.labels}) > 1:
            raise ValueError("all labels must have the same length")
        if self.groups is None:
            self.groups = list(range(len(self.labels)))

    @property
    def num_qubits(self) -> int:
        return len(self.labels[0])

    def __len__(self):
        return len(self.labels)

    def identity_coeff(self) -> float:
        idn = "I" * self.num_qubits
        return float(sum(c for l, c in zip(self.labels, self.coeffs) if l == idn))

    def without_identity(self) -> "PauliSum":
        idn = "I" * self.num_qubits
        keep = [k for k, l in enumerate(self.labels) if l != idn]
        return PauliSum([self.labels[k] for k in keep], self.coeffs[keep],
                        [self.groups[k] for k in keep])

    def shifted_scaled(self, shift: float, scale: float) -> "PauliSum":
        """(H - shift) / scale, with the identity term appended last."""
        h = self.without_identity()
        c0 = (self.identity_coeff() - shift) / scale
        return PauliSum(h.labels + ["I" * self.num_qubits],
                        np.concatenate([h.coeffs / scale, [c0]]),
                        h.groups + [-1])

    def l1_norm(self) -> float:
        return float(np.sum(np.abs(self.without_identity().coeffs)))

    def to_sparse(self) -> sp.csr_matrix:
        N = self.num_qubits
        M = sp.csr_matrix((2 ** N, 2 ** N), dtype=complex)
        for l, c in zip(self.labels, self.coeffs):
            M = M + c * pauli_sparse(l)
        return M.tocsr()


def pauli_sparse(label: str) -> sp.csr_matrix:
    x, z, ny = _masks(label)
    N = len(label)
    idx = np.arange(2 ** N)
    sgn = _parity(idx & z)
    data = (1j) ** ny * np.where(sgn, -1.0, 1.0)
    # (P psi)[c] = data[c^x] * psi[c^x]  ->  column c^x, row c
    return sp.csr_matrix((data, (idx ^ x, idx)), shape=(2 ** N, 2 ** N))


def local_matrix(terms: Dict[str, complex], qubits: Sequence[int], sparse: bool = False):
    """Matrix of sum_P c_P P restricted to `qubits` (P must be identity elsewhere)."""
    M = None
    for lab, c in terms.items():
        m = c * pauli_sparse("".join(lab[q] for q in qubits))
        M = m if M is None else M + m
    return M.tocsr() if sparse else M.toarray()


def operator_norm(terms: Dict[str, complex], max_exact_qubits: int = 12) -> Tuple[float, bool]:
    """Spectral norm of a Hermitian Pauli sum (real coefficients). Exact on its support when the
    support has at most `max_exact_qubits` qubits (dense eigvalsh up to 8 qubits, Lanczos above),
    otherwise the triangle-inequality bound sum |c|. Returns (norm, exact)."""
    terms = {l: c for l, c in terms.items() if abs(c) > 1e-14}
    if not terms:
        return 0.0, True
    qs = sorted({q for l in terms for q in support(l)})
    tri = float(sum(abs(c) for c in terms.values()))
    if len(qs) > max_exact_qubits:
        return tri, False
    if len(qs) <= 8:
        return float(np.max(np.abs(np.linalg.eigvalsh(local_matrix(terms, qs))))), True
    from scipy.sparse.linalg import eigsh
    M = local_matrix(terms, qs, sparse=True)
    lam = abs(eigsh(M, k=1, which="LM", return_eigenvectors=False, tol=1e-12)[0])
    return float(min(tri, lam * (1 + 1e-9))), True


# ----------------------------------------------------------------------------- statevector action
def _masks(label: str) -> Tuple[int, int, int]:
    N = len(label)
    x = z = ny = 0
    for q, c in enumerate(label):
        bit = 1 << (N - 1 - q)
        if c in "XY":
            x |= bit
        if c in "ZY":
            z |= bit
        if c == "Y":
            ny += 1
    return x, z, ny


def _parity(a: np.ndarray) -> np.ndarray:
    a = a.copy()
    p = np.zeros_like(a)
    while np.any(a):
        p ^= a & 1
        a >>= 1
    return p.astype(bool)


@lru_cache(maxsize=4096)
def _action(label: str):
    """(perm, phase) with (P psi) = phase * psi[perm]."""
    x, z, ny = _masks(label)
    idx = np.arange(2 ** len(label))
    src = idx ^ x
    phase = ((1j) ** ny) * np.where(_parity(src & z), -1.0, 1.0)
    return src, phase


def apply_pauli(label: str, psi: np.ndarray) -> np.ndarray:
    src, phase = _action(label)
    return phase * psi[src]


def apply_rotation(label: str, theta: float, psi: np.ndarray) -> np.ndarray:
    """exp(-i theta P) psi."""
    if set(label) == {"I"}:
        return np.exp(-1j * theta) * psi
    return np.cos(theta) * psi - 1j * np.sin(theta) * apply_pauli(label, psi)
