"""Check the fast engine against real Qiskit circuits (noiseless state vector), at small N.

For each (N, gap, overlap), the engine picks its cheapest filter; Qiskit then builds that exact circuit
(trial StatePreparation, one ancilla per pulse, controlled bonds via the KAK-synthesised Bell change), simulates
it, and runs amplitude amplification on it. Compared: success probability, 1 - F, CNOTs, AA success after k
rounds, and the CNOTs of the AA reflections.

    python verify_qiskit.py                                   # N = 4 6 8, a few gaps and overlaps
    python verify_qiskit.py --N 6 --gaps 0.3 1.0 --gammas 0.1 0.5 --compile pauli
    python verify_qiskit.py --max-qubits 18 --set J2=0.3      # model settings come from StudyConfig
"""
import argparse
import json
import os

from ornl import StudyConfig, chain_for_gap, measure_point
from ornl import qiskit_check as qk
from ornl.config import add_cli


def main():
    ap = add_cli(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter))
    ap.add_argument("--N", type=int, nargs="+", default=[4, 6, 8])
    ap.add_argument("--gaps", type=float, nargs="+", default=[0.3, 0.8, 1.5])
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.1, 0.5, 0.9])
    ap.add_argument("--eps", type=float, default=1e-2)
    ap.add_argument("--protocol", default="rus", choices=["rus", "aa"], help="whose cheapest filter to build")
    ap.add_argument("--compile", default="bond", choices=["bond", "pauli"])
    ap.add_argument("--max-qubits", type=int, default=20, help="skip points needing more (N + one ancilla per pulse)")
    ap.add_argument("--out", default=None, help="JSON file for the rows (default: <out_dir>/verify_qiskit.json)")
    a = ap.parse_args()
    cfg = StudyConfig.from_cli(a.config, a.set)

    rows = []
    hdr = (f"{'N':>2} {'gap':>5} {'gam':>5} {'K':>2} {'r':>2} {'n':>2} {'qb':>3} {'p engine':>10} {'p qiskit':>10} "
           f"{'1-F engine':>10} {'1-F qiskit':>10} {'CX model':>8} {'CX circ':>8} {'k':>2} {'AA sim':>7} {'AA form':>7} "
           f"{'refl circ/model':>15}")
    print(hdr)
    for N in a.N:
        for g in a.gaps:
            ch = chain_for_gap(N, g, cfg)
            best, eps_knob = measure_point(ch, a.gammas, [a.eps], cfg)
            for gam, b, e in zip(a.gammas, best, eps_knob):
                b = b[a.eps][a.protocol]
                if b is None:
                    continue
                nq = N + b["K"] * b["r"]
                if nq > a.max_qubits:
                    print(f"{N:2d} {g:5.2f} {gam:5.2f}  skipped: {nq} qubits > --max-qubits {a.max_qubits}")
                    continue
                v = qk.verify(ch, ch.trial(e), b["K"], b["r"], b["n"], compile=a.compile, eps=e)
                v.update(gap=g, gamma=gam, p_engine=b["p"], infid_engine=b["infid"])
                v = {k: (float(x) if hasattr(x, "dtype") else x) for k, x in v.items()}
                rows.append(v)
                print(f"{N:2d} {g:5.2f} {gam:5.2f} {b['K']:2d} {b['r']:2d} {b['n']:2d} {nq:3d} {b['p']:10.6f} "
                      f"{v['p']:10.6f} {b['infid']:10.3e} {v['infid']:10.3e} {v['cnot_formula']:8d} "
                      f"{v['cnot_circuit']:8d} {v['aa_k']:2d} {v['aa_success']:7.4f} {v['aa_success_formula']:7.4f} "
                      f"{v['reflect_cnot_circuit']:7d}/{v['reflect_cnot_model']:<7d}", flush=True)
    if rows:
        worst_p = max(abs(r["p"] - r["p_engine"]) for r in rows)
        worst_f = max(abs(r["infid"] - r["infid_engine"]) for r in rows)
        worst_aa = max(abs(r["aa_success"] - r["aa_success_formula"]) for r in rows)
        print(f"\nlargest differences: p {worst_p:.1e}, 1-F {worst_f:.1e}, AA success {worst_aa:.1e}")
        out = a.out or os.path.join(cfg.out_dir, "verify_qiskit.json")
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        json.dump(rows, open(out, "w"), indent=1)
        print(f"rows written to {out}")


if __name__ == "__main__":
    main()
