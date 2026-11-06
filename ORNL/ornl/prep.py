"""One-layer MPS-to-circuit trial preparation.

The reference ground state (exact or DMRG) is truncated to bond dimension 2 by sequential SVD from the right.
A bond-2 MPS in right-canonical form is exactly a staircase of N - 1 two-qubit gates (one layer):

    G_i on qubits (i, i+1):  |a>_i |0>_{i+1}  ->  sum_{s, b} B_i^s[a, b] |s>_i |b>_{i+1}

so qubit i+1 carries the bond into the next gate, and the last gate also writes the last site. Each G_i is the
isometry completed to a unitary; the KAK decomposition compiles it to at most 3 CNOTs.

S^z is conserved: the SVD is done block by block in the S^z of the right part, so each kept Schmidt vector has a
definite S^z. When singular values are degenerate across sectors (the triplet multiplet at even cuts of a singlet)
the one with the smallest |S^z| is kept, so the truncated state stays in S^z = 0, the sector the filter works in.

Basis: engine order, site 0 = most significant bit, bit 0 = spin up (S^z = +1/2).
"""
import numpy as np

TIE = 1e-8


def _sz(s):
    return 0.5 - s                                     # bit 0 -> +1/2, bit 1 -> -1/2


def truncate_bond2(full, N, chi=2):
    """Right-canonical bond-`chi` MPS of the 2^N vector `full` (engine order). Returns B[0..N-1] with
    B[i] of shape (chi_left, 2, chi_right); B[0] carries the norm."""
    B = [None] * N
    psi = np.asarray(full, dtype=complex).reshape(1 << (N - 1), 2)        # rows: sites 0..N-2; cols: site N-1
    q_right = np.array([0.0])                                             # S^z of the right bond states
    for i in range(N - 1, 0, -1):
        chi_r = len(q_right)
        cols_q = np.array([_sz(s) + q for s in (0, 1) for q in q_right])  # column index = s * chi_r + b
        svals, us, vs, qs = [], [], [], []
        for q in np.unique(cols_q):
            idx = np.flatnonzero(np.isclose(cols_q, q))
            U, S, Vh = np.linalg.svd(psi[:, idx], full_matrices=False)
            for k in range(len(S)):
                svals.append(S[k]); us.append(U[:, k]); qs.append(q)
                v = np.zeros(2 * chi_r, dtype=complex); v[idx] = Vh[k]; vs.append(v)
        svals = np.array(svals); qs = np.array(qs)
        order = sorted(range(len(svals)), key=lambda k: (-np.round(svals[k] / svals.max() / TIE), abs(qs[k])))
        keep = order[:chi]
        Vk = np.array([vs[k] for k in keep])                              # (chi, 2 chi_r)
        B[i] = Vk.reshape(len(keep), 2, chi_r)
        left = np.stack([us[k] * svals[k] for k in keep], axis=1)         # (2^i, chi)
        q_right = qs[keep]
        if i > 1:
            psi = left.reshape(1 << (i - 1), 2 * len(keep))               # rows: sites 0..i-2; cols: (s_{i-1}, a)
        else:
            psi = left.reshape(1, 2 * len(keep))
    B[0] = psi.reshape(1, 2, -1)
    return B


def mps_to_vector(B):
    v = B[0].reshape(2, -1)
    for T in B[1:]:
        v = np.einsum("xa,asb->xsb", v, T).reshape(-1, T.shape[2])
    return v.reshape(-1)


def _complete(iso):
    """4 x k isometry (orthonormal columns) -> 4 x 4 unitary whose first k columns are `iso`."""
    k = iso.shape[1]
    Q, _ = np.linalg.qr(np.hstack([iso, np.random.default_rng(0).normal(size=(4, 4 - k))]))
    Q[:, :k] = iso
    return Q


def staircase(B):
    """Two-qubit unitaries G[i] on (i, i+1), i = 0..N-2, with |psi> = G[N-2] ... G[0] |0...0> (engine order:
    G[i] row/column index = 2 s_i + s_{i+1})."""
    N = len(B)
    Bs = [T.copy() for T in B]
    Bs[0] = Bs[0] / np.linalg.norm(Bs[0])
    pad = lambda T: np.pad(T, ((0, 2 - T.shape[0]), (0, 0), (0, 2 - T.shape[2])))   # noqa: E731  bond dim -> 2
    Bs = [pad(T) for T in Bs]
    gates = []
    for i in range(N - 1):
        T = Bs[i]                                                     # (a, s, b)
        if i == N - 2:                                                # merge the last site: b -> s_{N-1}
            T = np.einsum("asb,tb->ast", T, Bs[N - 1][:, :, 0].T)
        iso = T.reshape(2, 4).T                                       # column a: sum_{s,b} T[a,s,b] |s b>
        if i == 0:
            iso = iso[:, :1]                                          # input |0>|0>
        else:
            norms = np.linalg.norm(iso, axis=0)
            iso = iso[:, norms > 1e-12]
        cols = iso
        if i > 0 and cols.shape[1] < 2:                               # bond dim 1 here: input a = 0 only
            pass
        U = _complete(cols)
        # gate acts on |a>_i |0>_{i+1}: input index 2 a + 0 -> place the isometry columns at inputs 0 and 2
        G = np.zeros((4, 4), dtype=complex)
        ins = [0, 2][:cols.shape[1]]
        rest = [x for x in range(4) if x not in ins]
        G[:, ins] = U[:, :cols.shape[1]]
        G[:, rest] = U[:, cols.shape[1]:]
        gates.append(G)
    return gates


def apply_staircase(gates, N):
    v = np.zeros(1 << N, dtype=complex)
    v[0] = 1
    for i, G in enumerate(gates):
        t = v.reshape(1 << i, 4, 1 << (N - i - 2))
        v = np.einsum("xy,ayb->axb", G, t).reshape(-1)
    return v


def kak_cnots(G):
    """CNOTs of one two-qubit gate under Qiskit's KAK decomposer (3 if Qiskit is not installed)."""
    try:
        from qiskit.circuit.library import CXGate
        from qiskit.synthesis import TwoQubitBasisDecomposer
    except ImportError:
        return 3
    return int(TwoQubitBasisDecomposer(CXGate()).num_basis_gates(G))


def prepare(ch):
    """One-layer MPS-to-circuit approximation of ch.g0. Returns (sector vector of the prepared state, gates,
    CNOTs, gamma_0 = |<E0|prepared>|^2)."""
    N = ch.N
    full = np.zeros(1 << N, dtype=complex)
    full[ch.states] = ch.g0
    gates = staircase(truncate_bond2(full, N))
    out = apply_staircase(gates, N)
    leak = 1 - np.sum(abs(out[ch.states]) ** 2)
    assert leak < 1e-10, f"prepared state left the S^z = 0 sector (weight {leak:.1e})"
    v = out[ch.states]
    return v, gates, sum(kak_cnots(G) for G in gates), float(abs(np.vdot(ch.g0, v)) ** 2)
