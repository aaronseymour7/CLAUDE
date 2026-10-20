"""validate_two_term.py -- random-matrix check of the two-term Strang bound used (by telescoping) in alpha2.py:
   || e^{-itA/2} e^{-itB} e^{-itA/2} - e^{-it(A+B)} ||  <=  t^3 ( ||[B,[B,A]]||/12 + ||[A,[A,B]]||/24 )   for Hermitian A, B.
Prints the max ratio error/bound over random cases (a ratio near 1 means the constants cannot be lowered; > 1 would refute them)."""
import numpy as np, scipy.linalg as sl
rng = np.random.default_rng(0); worst = 0.0
herm = lambda n: (lambda M: (M + M.conj().T) / 2)(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
c = lambda X, Y: X @ Y - Y @ X
nrm = lambda X: np.linalg.norm(X, 2)
for _ in range(300):
    n = int(rng.integers(2, 7)); A, B = herm(n) * rng.uniform(.2, 2), herm(n) * rng.uniform(.2, 2)
    for t in (0.05, 0.3, 1, 3, 10, 30):
        S = sl.expm(-1j * t * A / 2) @ sl.expm(-1j * t * B) @ sl.expm(-1j * t * A / 2)
        worst = max(worst, nrm(S - sl.expm(-1j * t * (A + B))) / (t ** 3 * (nrm(c(B, c(B, A))) / 12 + nrm(c(A, c(A, B))) / 24)))
print("max error/bound over 1800 random cases:", round(worst, 4))
