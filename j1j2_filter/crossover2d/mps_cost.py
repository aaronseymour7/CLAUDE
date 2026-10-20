"""mps_cost.py -- bond dimension and gate cost of the exact-ish sequential MPS state-preparation circuit from the Schmidt spectra saved by ed_measure.py.
chi_k(eps): smallest bond dimension at cut k whose discarded weight is <= eps/(N-1) (first-order infidelity accounting: the sum of per-cut discarded weights).
Cost model (RZZ gates): the isometry at site k acts on ceil(log2 max(chi_{k-1}, chi_k)) + 1 qubits; a generic n-qubit unitary costs (23/48) 4^n CNOT (quantum Shannon
decomposition), i.e. the same number of RZZ.  This is an optimistic estimate for the competitor (isometries are cheaper than unitaries, lower-depth routes exist)."""
import json, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent


def chi_per_cut(spectra, eps):
    N = len(spectra) + 1; out = []
    for s in spectra:
        s = np.sort(np.array(s))[::-1]; tail = s.sum() - np.cumsum(s)         # tail[i] = weight discarded when keeping i+1 values
        k = int(np.searchsorted(-(tail <= eps / (N - 1)).astype(int), 0)) if False else int(np.argmax(tail <= eps / (N - 1))) + 1
        out.append(min(k, len(s)))
    return out


def mps_rzz(chis):
    N = len(chis) + 1; tot = 0.0; prev = 1
    for c in chis + [1]:
        n = int(np.ceil(np.log2(max(prev, c)))) + 1
        tot += 23.0 / 48.0 * 4.0 ** n; prev = c
    return tot


if __name__ == "__main__":
    print(f"{'WxL':>6s} {'N':>3s} {'chi_max(1e-2)':>13s} {'chi_max(1e-3)':>13s} {'RZZ_MPS(1e-2)':>14s} {'S_mid (nats)':>12s}")
    for f in sorted(glob.glob(str(HERE / "results" / "ed_*_J20.json")), key=lambda f: (json.load(open(f))["W"], json.load(open(f))["N"])):
        r = json.load(open(f)); sp = r["spectra"]; mid = sp[len(sp) // 2]; p = np.array(mid); p = p[p > 1e-16]
        c2, c3 = chi_per_cut(sp, 1e-2), chi_per_cut(sp, 1e-3)
        print(f"{r['W']}x{r['L']:<4d} {r['N']:>3d} {max(c2):13d} {max(c3):13d} {mps_rzz(c2):14.3e} {float(-(p * np.log(p)).sum()):12.3f}")
