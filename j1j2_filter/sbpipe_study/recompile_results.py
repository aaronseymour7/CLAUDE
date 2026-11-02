"""Re-count every stored filter CNOT / Rz under the bond compile (sbpipe.resources, compile="bond").

The study was run with the ladder rule (each Pauli term (x) Z_anc a separate rotation, 4 CNOTs for a
weight-2 term). The bond compile applies the three commuting terms of a Heisenberg bond as one exact
6-CNOT block, so the circuits, Trotter errors, infidelities and success probabilities are unchanged and
only the counts move. Counts are linear in the totals, for B bonds, S Strang steps and m pulses/queries:

    order 2   ladder  CNOT = S (24B - 8) + 4m     Rz = S (6B - 2) + m
              bond    CNOT = S (12B - 12) + 6m    Rz = S (6B - 6) + 3m
    order 1   ladder  CNOT = 12 B S               bond  CNOT = 6 B S      (Rz = 3 B S in both)

(+ the ancilla rotations per pulse, unchanged). Every stored filter count is checked against the ladder
formula before it is replaced, and the stored success probability (measured, or the certified bound
recovered from the stored per-success cost) is kept. Run once: results/COMPILE records the conversion.
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
RES = os.path.join(HERE, "results")

from sbpipe import j1j2_chain                    # noqa: E402
from sbpipe.models import j1j2_square            # noqa: E402

LY = 4
_B = {}
stats = {}


def bonds(rec):
    kind = rec.get("kind", "chain")
    if kind not in ("chain", "cyl"):
        kind = "chain"
    size = rec.get("size", rec["N"]) if kind == "cyl" else rec["N"]
    key = (kind, size, round(rec["J2"], 9))
    if key not in _B:
        H = j1j2_chain(size, 1.0, rec["J2"]) if kind == "chain" else j1j2_square(size, LY, 1.0, rec["J2"])
        h = H.without_identity()
        assert len(h) % 3 == 0 and len(set(h.groups)) == len(h) // 3
        _B[key] = len(h) // 3
    return _B[key]


def old(S, m, B, order=2):
    return S * (24 * B - 8) + 4 * m if order == 2 else 12 * B * S


def new(S, m, B, order=2):
    return S * (12 * B - 12) + 6 * m if order == 2 else 6 * B * S


def rz_shift(S, m, B, order=2):
    """bond minus ladder rotation count."""
    return (-4 * S + 2 * m) if order == 2 else 0


def check(f_old, S, m, B, order, tag, tol=1e-6):
    want = old(S, m, B, order)
    if abs(f_old - want) > tol * max(1.0, want) + 0.6:
        raise ValueError(f"{tag}: stored filter CNOT {f_old} != ladder formula {want} (S={S}, m={m}, B={B})")


def bump(tag):
    stats[tag] = stats.get(tag, 0) + 1


def meas(d, S, m, B, trial, tag, pkey="p", ckey="cnot", order=2):
    """Measured result: per-success cost (filter + trial) / p with p stored."""
    if d is None or ckey not in d or d.get(pkey) is None:
        return
    p = d[pkey]
    if trial is None:
        trial = d[ckey] * p - d["cnot_filter"]
    f_old = d.get("cnot_filter", d[ckey] * p - trial)
    if B is None:                                   # infer it from the stored ladder count
        B = round(((f_old - 4 * m) / S + 8) / 24) if order == 2 else round(f_old / (12 * S))
    check(f_old, S, m, B, order, tag, tol=1e-6)
    f = new(S, m, B, order)
    d[ckey] = (f + trial) / p
    if "cnot_filter" in d:
        d["cnot_filter"] = f
    bump(tag)


def cert(d, S, m, B, trial, tag, ckey="cnot", order=2):
    """Certified plan: per-success cost with the (unstored) certified p recovered from the old cost."""
    if d is None or d.get(ckey) is None or S is None:
        return
    p = (old(S, m, B, order) + trial) / d[ckey]
    d[ckey] = (new(S, m, B, order) + trial) / p
    bump(tag)


def plan(P, B, tag):
    if "res_filter_cnot" not in P or not P.get("res_steps"):
        return
    order, S, m = P.get("order", 2), P["res_steps"], P["m"]
    check(P["res_filter_cnot"], S, m, B, order, tag)
    f, dr = new(S, m, B, order), rz_shift(S, m, B, order)
    P["res_filter_cnot"] = f
    P["res_attempt_cnot"] = f + P["res_trial_cnot"]
    P["res_filter_rotations"] += dr
    P["res_attempt_rotations"] += dr
    P["res_attempt_t"] = P["res_attempt_rotations"] * P["res_t_per_rotation"]
    P["res_expected_cnot"] = P["res_attempt_cnot"] / P["res_p_succ"]
    P["res_expected_t"] = P["res_attempt_t"] / P["res_p_succ"]
    bump(tag)


def candidates(C, B, trial, order, tag):
    for c in C or []:
        if c.get("expected_cnot") and c.get("steps"):
            cert(c, c["steps"], c["m"], B, trial, tag, ckey="expected_cnot", order=order)


def walk_plans(rec, B, tag):
    if isinstance(rec, dict):
        if "res_filter_cnot" in rec:
            plan(rec, B, tag)
        for v in rec.values():
            walk_plans(v, B, tag)
    elif isinstance(rec, list):
        for v in rec:
            walk_plans(v, B, tag)


TRIAL = {}


GTRIAL = {}


def _trials():
    for l in open(os.path.join(RES, "gapov.jsonl")):
        r = json.loads(l)
        assert GTRIAL.setdefault(r["N"], r["trial_cnot"]) == r["trial_cnot"]
    for f in ("emp_actual1d.jsonl", "emp_actual.json"):
        p = os.path.join(RES, f)
        recs = [json.loads(l) for l in open(p)] if f.endswith(".jsonl") else json.load(open(p))
        for r in recs:
            TRIAL[(r["kind"], r["size"], round(r["J2"], 9), r["chi"])] = r["trial_cnot"]


def do_record(f, r):
    B = bonds(dict(r, N=r.get("N", r.get("size")))) if r.get("J2") is not None and ("N" in r or "size" in r) else None
    tr = r.get("trial_cnot")
    if tr is None and f.startswith("emp_actual") and "chi" in r:
        tr = TRIAL[(r["kind"], r["size"], round(r["J2"], 9), r["chi"])]
    if f.startswith("gapov") and f != "gapov_tune.json":
        if B is None:
            B = 2 * r["N"] - 3                      # J2 is tuned (non-zero) in the gap/overlap study
        if tr is None:
            tr = GTRIAL[r["N"]]
        meas(r.get("meas_dyadic"), *(lambda d: (d["steps"], d["m"]) if d else (0, 0))(r.get("meas_dyadic")), B, tr, f + ":meas_dyadic")
        q = r.get("meas_qetu")
        if q:
            meas(q, q["d"] * q["k"], q["d"], B, tr, f + ":meas_qetu")
        o = r.get("meas_optimized")
        if o:
            meas(o, o["steps"], o["m"], B, tr, f + ":meas_opt")
        for k in ("cert_dyadic", "cert_qetu", "cert_optimized"):
            c = r.get(k)
            if c:
                cert(c, c["steps"], c["m"], B, tr, f + ":" + k)
    elif f.startswith("emp_actual"):
        d = r.get("dyadic")
        if d:
            meas(d, d["steps"], d["K"] * d["r"], B, tr, f + ":dyadic")
        o = r.get("optimized")
        if o:
            meas(o, o["steps"], o["m"], B, tr, f + ":optimized")
        q = r.get("qetu")
        if q:
            meas(q, q["d"] * q["k"], q["d"], B, tr, f + ":qetu")
    elif f == "emp_sbc.json":
        if r.get("cnot_per_success") is None:
            return
        trial = r["cnot_per_success"] * r["p_succ"] - r["cnot_filter"]
        meas(r, r["steps"], r["emp_m"], B, trial, f + ":meas", pkey="p_succ", ckey="cnot_per_success")
        if r.get("cert_cnot_per_success"):
            cert(r, r["cert_steps"], r["cert_m"], B, trial, f + ":cert", ckey="cert_cnot_per_success")
    elif f == "emp_trotter.json":
        if r.get("cnot_per_success") is None:
            return
        trial = r["cnot_per_success"] * r["p_succ_trotter"] - r["cnot_filter"]
        meas(r, r["steps_min"], r["d_trot"], B, trial, f + ":meas", pkey="p_succ_trotter",
             ckey="cnot_per_success")
        if r.get("cert_cnot_per_success"):
            p = r["cert_p_succ"]
            S, m = r["cert_steps"], r["cert_d"]
            r["cert_cnot_per_success"] = (new(S, m, B) + trial) / p
            bump(f + ":cert")
    elif f == "emp_frontier.json":
        g = r["gamma"]
        pts = []
        for eps, d, c in r["points"]:
            if d:
                check(c * g - tr, d, d, B, 2, f)
                c = (new(d, d, B) + tr) / g
                bump(f)
            pts.append([eps, d, c])
        r["points"] = pts
    else:
        walk_plans(r, B, f)
        if "candidates" in r and "plan" in r:
            candidates(r["candidates"], B, r["plan"].get("res_trial_cnot", 0), r["plan"].get("order", 2),
                       f + ":candidates")


def main():
    flag = os.path.join(RES, "COMPILE")
    if os.path.exists(flag):
        print("already converted:", open(flag).read().strip())
        return
    _trials()
    out = []
    for p in sorted(glob.glob(os.path.join(RES, "*.json")) + glob.glob(os.path.join(RES, "*.jsonl"))):
        f = os.path.basename(p)
        if f in ("gapov_filter.json", "gapov_tune.json"):           # derived; regenerate with derive_gapov_filter.py
            continue
        if f.endswith(".jsonl"):
            recs = [json.loads(l) for l in open(p) if l.strip()]
        else:
            recs = json.load(open(p))
        lst = recs if isinstance(recs, list) else [recs]
        for r in lst:
            if isinstance(r, dict):
                try:
                    do_record(f, r)
                except Exception as e:
                    raise RuntimeError(f"{f}: {e!r} in {str(r)[:200]}")
        out.append((p, f, recs, lst))
    for p, f, recs, lst in out:                     # write only once everything converted
        with open(p, "w") as fh:
            if f.endswith(".jsonl"):
                fh.write("".join(json.dumps(r) + "\n" for r in lst))
            else:
                json.dump(recs, fh, indent=1)
    with open(flag, "w") as fh:
        fh.write("bond6: filter counts use the 6-CNOT controlled-bond block (recompile_results.py)\n")
    for k in sorted(stats):
        print(f"{stats[k]:6d}  {k}")


if __name__ == "__main__":
    main()
