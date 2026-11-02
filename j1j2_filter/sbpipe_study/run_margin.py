"""Gap margin from DMRG convergence data, tested against exact gaps.

The certified plans are planned from a DMRG gap Delta(chi) = E1(chi) - E0(chi) (E1 = lowest S^z = 1 state), and
an overestimated gap is the input error that breaks them. A rule that uses only DMRG data:

    Delta_plan = Delta(chi) - c * |Delta(chi) - Delta(chi/2)|        (c = 1, 2, 4)

i.e. a margin set by the last change of the gap along the bond-dimension ladder. For every system with an exact
gap (chains N = 12, 14, 16 and cylinders 4x3, 4x4, 4x5, all J2) we run the DMRG ladder chi = 4 ... 128 for the
ground state and the S^z = 1 state, and record at each final chi whether Delta_plan <= Delta_exact (the margin
covers the error) and f = Delta_plan / Delta_exact (how much gap the plan gives away; f < 1 is safe).
-> results/margin.json

    python sbpipe_study/run_margin.py
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import json
import time

import numpy as np
from concurrent.futures import ProcessPoolExecutor

import run_extended as RX
import run_transfer20 as RT

RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "margin.json")
LADDER = (4, 8, 16, 32, 64, 128)


def systems():
    S = []
    for r in map(json.loads, open(os.path.join(RES, "aposteriori.jsonl"))):
        k = (r["kind"], r["size"], r["J2"])
        if k not in [s[:3] for s in S]:
            S.append(k + (r["gap"],))
    for J2 in (0.0, 0.3):
        E = json.load(open(os.path.join(RT.CACHE, RT._key("cyl", 5, J2) + ".json")))["E01"]
        S.append(("cyl", 5, J2, E[1] - E[0]))
    return S


def _one(s):
    from sbpipe.frontends import dmrg_quimb as DQ
    kind, size, J2, gap = s
    t0 = time.time()
    H = RX.model(kind, size, J2)
    mpo = DQ.pauli_to_mpo(H)
    g = DQ.ground_state(mpo, bond_dims=LADDER)
    e = DQ.sector_ground_state(mpo, DQ.total_sz(H.num_qubits), 1.0, 2.0, bond_dims=LADDER)
    return dict(kind=kind, size=size, J2=J2, N=H.num_qubits, gap_exact=gap,
                E0=[h[1] for h in g["history"]], E1=[h[1] for h in e["history"]], chis=list(LADDER),
                time=time.time() - t0)


def main():
    res = json.load(open(OUT)) if os.path.exists(OUT) else []
    have = {(r["kind"], r["size"], r["J2"]) for r in res}
    todo = [s for s in systems() if s[:3] not in have]
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 2))) as ex:
        for r in ex.map(_one, todo):
            res.append(r)
            print(r["kind"], r["size"], r["J2"], f"{r['time']:.0f}s", flush=True)
            json.dump(res, open(OUT, "w"), indent=1)


# ----------------------------------------------------------------------------- what the filter does with each margin
BOUT = os.path.join(RES, "margin_behave.jsonl")


def _filter_values(des, Es):
    if des.schedule == "qetu":
        from sbpipe.design import _cheb
        dm = des.meta
        x = np.cos(dm["tau"] * Es)
        return _cheb(dm["degree"], x / dm["x_s"]) / _cheb(dm["degree"], np.array([1.0 / dm["x_s"]]))[0]
    f = np.ones_like(Es)
    for t, ph in zip(des.times, des.phases):
        f = f * np.cos(Es * t + ph)
    return f


def _behave(args):
    """Plan from DMRG energies at one rung of the ladder with error bars c * (last change), and apply the exact
    filter function to the exact spectrum and spectral weights of the trial (no Trotter error)."""
    import run_empirical as RE
    from sbpipe import Ansatz, EnergyEstimates
    from sbpipe.pipeline import Config, refine
    import sbpipe.pipeline as PL
    r, i, c, chi, base = args
    sys_ = (r["kind"], r["size"], r["J2"])
    c0 = float(RX.model(*sys_).without_identity().coeffs[0])
    # beta2 is homogeneous of degree 3 in the coefficients (identity excluded): rescale the unshifted value
    PL.trotter_constant = lambda Hs, order: base * (float(Hs.without_identity().coeffs[0]) / c0) ** 3
    z = np.load(os.path.join(RE.CACHEE, RE._act_key(*sys_) + ".npz"))
    E, w = z["E"], z[f"w{chi}"]
    H = RX.model(*sys_)
    d0 = abs(r["E0"][i] - r["E0"][i - 1]); d1 = abs(r["E1"][i] - r["E1"][i - 1])
    est = EnergyEstimates(E0=r["E0"][i], E1=r["E1"][i], E0_err=max(c * d0, 1e-10), E1_err=max(c * d1, 1e-10))
    gam = float(w[0])
    out = []
    for eps in (1e-2, 1e-3, 1e-4):
        if 1 - gam <= eps:
            continue
        for s in ("dyadic", "optimized", "qetu"):
            row = dict(kind=r["kind"], size=r["size"], J2=r["J2"], chi_dmrg=r["chis"][i], c=c, chi=chi, eps=eps,
                       schedule=s, gap_plan=(est.E1 - est.E1_err) - (est.E0 - est.E0_err), gap_exact=r["gap_exact"])
            try:
                pl = refine(H, est, Ansatz(gamma=gam, statevector=None, cost=dict(cnot=0), name="t", meta={}),
                            Config(eps=eps, schedule=s, order=2))
            except (RuntimeError, ValueError) as e:
                row["error"] = str(e)[:80]
                out.append(row)
                continue
            Es = (E - pl.window["shift"]) / pl.window["W"]
            f = _filter_values(pl.design, Es)
            q = w * f ** 2
            row.update(inf=float(1 - q[0] / q.sum()), cert_cnot=pl.resources["expected_cnot"], m=pl.design.m)
            out.append(row)
    return out


def _beta(sys_):
    from sbpipe.pipeline import beta2
    return list(sys_), beta2(RX.model(*sys_).shifted_scaled(0.0, 1.0))


def behave():
    import run_empirical as RE
    R = [r for r in json.load(open(OUT)) if (r["kind"], r["size"], r["J2"]) in RE.ACT_SYSTEMS]
    bp = os.path.join(RES, "margin_beta2.json")
    B = {tuple(k): v for k, v in json.load(open(bp))} if os.path.exists(bp) else {}
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 2))) as ex:
        for k, v in ex.map(_beta, [s for s in RE.ACT_SYSTEMS if s not in B]):
            B[tuple(k)] = v
            print("beta2", k, v, flush=True)
            json.dump([[list(k), v] for k, v in B.items()], open(bp, "w"))
    done = set()
    if os.path.exists(BOUT):
        done = {(x["kind"], x["size"], x["J2"], x["chi_dmrg"], x["c"], x["chi"]) for x in map(json.loads, open(BOUT))}
    todo = [(r, i, c, chi, B[(r["kind"], r["size"], r["J2"])]) for r in R for i in range(1, len(r["chis"]))
            for c in (0, 1, 2) for chi in (16, 32)
            if (r["kind"], r["size"], r["J2"], r["chis"][i], c, chi) not in done]
    print(len(todo), "tasks", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("JOBS", 2))) as ex, open(BOUT, "a") as fh:
        for rows in ex.map(_behave, todo):
            fh.write("".join(json.dumps(x) + "\n" for x in rows)); fh.flush()
            if rows:
                print(rows[0]["kind"], rows[0]["size"], rows[0]["J2"], rows[0]["chi_dmrg"], rows[0]["c"], len(rows), flush=True)


if __name__ == "__main__":
    behave() if sys.argv[1:] == ["behave"] else main()
