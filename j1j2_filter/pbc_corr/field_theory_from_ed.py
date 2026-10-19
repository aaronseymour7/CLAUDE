"""field_theory_from_ed.py -- what field-theory data can be read off exact periodic J1-J2 ground states at the sizes a filter can reach (N = 8..20)?
For each J2 and N: E0, lowest triplet gap, singlet gap, half-chain entanglement entropy S(N/2), translation-averaged SzSz correlators.
Field-theory quantities (c = 1, SU(2)_1 at J2 <= J2c) extracted from them:
  v     from the triplet gap    :  gap = 2 pi v x / N with x = 1/2          ->  v = N gap / pi
  c     from the Casimir energy :  E0 = N e_inf - pi c v / (6 N)            (fit of e_inf and c*v, divided by v)
  c     from entanglement       :  S(N/2; N) = (c/3) ln N + const           (slope between consecutive sizes)
  eta   from correlations       :  N |C(N/2)| ~ N^{1-eta} (times logs)      (local log-log slope between same-family sizes)
  ratio singlet/triplet gap     :  1 (periodic: one SU(2)xSU(2) multiplet at the fixed point)
usage: python pbc_corr/field_theory_from_ed.py [J2 ...]"""
import sys, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from pbc_ctx import ed_inputs, correlations

J2S = [float(a) for a in sys.argv[1:]] or [0.0, 0.2411, 0.4]
NS = (8, 10, 12, 14, 16, 18, 20)


def entropy(psi, N, l):
    s = np.linalg.svd(psi.reshape(1 << (N - l), 1 << l), compute_uv=False) ** 2; s = s[s > 1e-16]
    return float(-(s * np.log(s)).sum())


out = {}
for J2 in J2S:
    rows = []
    for N in NS:
        ed = ed_inputs(N, 1.0, J2); psi = ed["psi"]
        rows.append(dict(N=N, E0=ed["E0"], gap_t=ed["E1_triplet"] - ed["E0"], gap_s=ed["E1_singlet"] - ed["E0"], S_half=entropy(psi, N, N // 2),
                         S_l={l: entropy(psi, N, l) for l in range(1, N // 2 + 1)}, C=correlations(psi, N)))
        print("J2", J2, "N", N, flush=True)
    out[J2] = rows
json.dump({str(k): v for k, v in out.items()}, open(HERE / "results" / "field_theory_ed.json", "w"), default=float)
