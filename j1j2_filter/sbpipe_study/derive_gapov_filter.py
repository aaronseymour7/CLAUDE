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
from sbpipe.resources import per_step_counts

import run_study as RS


def counts(H, steps, pulses, anc):
    """filter_counts() for a second-order schedule given only the totals (it is linear in them)."""
    s = per_step_counts(H)
    G, cw = s["terms"], s["cnot_per_term"]
    inner = cw[1:-1].sum()
    cnot = steps * (2 * inner + cw[0] + cw[-1]) + pulses * cw[0]
    rot = steps * (2 * G - 2) + pulses + anc * pulses
    return int(cnot), int(rot)


def main():
    rows = [json.loads(l) for l in open(os.path.join(RS.RES, "gapov.jsonl"))]
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
