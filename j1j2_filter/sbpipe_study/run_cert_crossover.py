"""Certified (implementable) version of the measured crossovers: for every (system, eps, trial chi) whose
trial circuit is cheaper than the cheapest MPS that meets eps alone, plan the filter with the pipeline
(Theorem 1, exact energies, schedules 'auto' and 'qetu') -> results/cert_crossover.jsonl."""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import run_extended as RX

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "cert_crossover.jsonl")


def _records():
    A = [json.loads(l) for l in open(os.path.join(RES, "emp_actual1d.jsonl"))]
    A += json.load(open(os.path.join(RES, "emp_actual.json")))
    return A


def cases():
    A = _records()
    out, seen = [], set()
    for r in A:
        if r["mps_alone"]:
            continue
        sys_ = (r["kind"], r["size"], r["J2"], r["eps"])
        m = [q["trial_cnot"] for q in A if q["mps_alone"] and (q["kind"], q["size"], q["J2"], q["eps"]) == sys_]
        if m and r["trial_cnot"] >= min(m):
            continue
        for s in ("auto", "qetu"):
            k = sys_ + (r["chi"], s)
            if k not in seen:
                seen.add(k)
                out.append((r["kind"], r["size"], r["J2"], r["eps"], r["chi"], s))
    return out


def main():
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l)
            done.add((r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"]))
    todo = [c for c in cases() if c not in done]
    todo.sort(key=lambda c: c[5] != "qetu")             # the fast QETU plans first
    print(len(todo), "plans", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 4))) as ex, open(OUT, "a") as fh:
        for r in ex.map(RX._regime_case, todo):
            pl = r["plan"]
            fh.write(json.dumps(r) + "\n"); fh.flush()
            print(r["kind"], r["size"], r["J2"], r["eps"], r["chi"], r["schedule"],
                  f"{pl['res_expected_cnot']:.3g}" if pl else r.get("error"), flush=True)


if __name__ == "__main__":
    main()
