"""Filter-only resources for every gap/overlap case (results/gapov.jsonl -> results/gapov_filter.json).

Per shot: CNOTs, arbitrary Rz rotations, Strang steps, Hamiltonian evolution time (physical units 1/J1) and
pulse/query count. Per success: the same divided by the success probability. Trial-state cost excluded.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from sbpipe import j1j2_chain
from sbpipe.resources import per_block_counts

import run_study as RS


def counts(H, steps, pulses, anc):
    """filter_counts() for a second-order schedule given only the totals (it is linear in them)."""
    cw, rw = per_block_counts(H)
    cnot = steps * (2 * cw[1:-1].sum() + cw[0] + cw[-1]) + pulses * cw[0]
    rot = steps * (2 * rw[1:-1].sum() + rw[0] + rw[-1]) + pulses * rw[0] + anc * pulses
    return int(cnot), int(rot)


def _key(r):
    return (r["N"], r["gap_target"], r["gamma_target"], r["eps"])


def _jsonl(name):
    p = os.path.join(RS.RES, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def optimized_map():
    """(N, gap, gamma, eps) -> {meas_optimized, cert_optimized}: the single-seed run, replaced by the retry
    (more optimizer seeds / leak fractions) wherever the retry is cheaper."""
    out = {}
    for r in _jsonl("gapov_opt.jsonl"):
        out[_key(r)] = {k: r[k] for k in ("meas_optimized", "cert_optimized")}
    for r in _jsonl("gapov_opt_fix.jsonl"):
        o = out.setdefault(_key(r), {"meas_optimized": None, "cert_optimized": None})
        for k in ("meas_optimized", "cert_optimized"):
            if k in r and r[k] and (o[k] is None or r[k]["cnot"] < o[k]["cnot"]):
                o[k] = r[k]
    return out


def load_cases(files=("gapov.jsonl",)):
    """gap/overlap cases with the optimized results attached (None where not run)."""
    opt = optimized_map()
    rows = []
    for f in files:
        for r in _jsonl(f):
            r.update(opt.get(_key(r), {"meas_optimized": None, "cert_optimized": None}))
            rows.append(r)
    return rows


def main():
    rows = load_cases()
    out, cache = [], {}
    for r in rows:
        key = (r["N"], r["J2"])
        if key not in cache:
            cache[key] = j1j2_chain(r["N"], 1.0, r["J2"])
        H = cache[key]
        W = r["W"]
        o = {k: r[k] for k in ("study", "N", "gap", "gap_target", "J2", "gamma", "gamma_target", "eps", "W", "delta")}
        md, mq = r["meas_dyadic"], r["meas_qetu"]
        cd, cq = r["cert_dyadic"], r["cert_qetu"]
        ent = {}
        c, z = counts(H, md["steps"], md["m"], 1)
        ent["meas_dyadic"] = dict(cnot=c, rz=z, steps=md["steps"], pulses=md["m"], time=md["T"] / W, p=md["p"],
                                  dt=md["dt"] / W)
        tq = np.pi / 2
        c, z = counts(H, mq["d"] * mq["k"], mq["d"], 2)
        ent["meas_qetu"] = dict(cnot=c, rz=z, steps=mq["d"] * mq["k"], pulses=mq["d"], time=mq["d"] * tq / W,
                                p=mq["p"], dt=tq / mq["k"] / W)
        for name, cc, anc in (("cert_dyadic", cd, 1), ("cert_qetu", cq, 2)):
            c, z = counts(H, cc["steps"], cc["m"], anc)
            # certified success probability: the bound's expected CNOTs = (filter + trial) / p
            p = (c + r["trial_cnot"]) / cc["cnot"]
            ent[name] = dict(cnot=c, rz=z, steps=cc["steps"], pulses=cc["m"], time=cc["T"] / W, p=p,
                             dt=cc["T"] / cc["steps"] / W)
        mo, co = r["meas_optimized"], r["cert_optimized"]
        ent["meas_optimized"] = ent["cert_optimized"] = None
        if mo:
            c, z = counts(H, mo["steps"], mo["m"], 1)
            assert c == mo["cnot_filter"], (c, mo["cnot_filter"])
            ent["meas_optimized"] = dict(cnot=c, rz=z, steps=mo["steps"], pulses=mo["m"], time=mo["T"] / W,
                                         p=mo["p"], dt=mo["dt"] / W)
        if co:
            c, z = counts(H, co["steps"], co["m"], 1)
            ent["cert_optimized"] = dict(cnot=c, rz=z, steps=co["steps"], pulses=co["m"], time=co["T"] / W,
                                         p=(c + r["trial_cnot"]) / co["cnot"], dt=co["T"] / co["steps"] / W)
        o.update(ent)
        # consistency: reconstructed measured CNOTs must reproduce the stored per-success total
        for k, m in (("meas_dyadic", md), ("meas_qetu", mq)):
            tot = (ent[k]["cnot"] + r["trial_cnot"]) / m["p"]
            assert abs(tot / m["cnot"] - 1) < 1e-6, (k, tot, m["cnot"])
        out.append(o)
    RS.dump("gapov_filter.json", out)
    print(len(out), "cases")


if __name__ == "__main__":
    main()
