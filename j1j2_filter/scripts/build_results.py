"""
build_results.py -- turn results/data/*.json (written by run_results.py) into results/RESULTS.md,
CSV tables and figures.   usage: python scripts/build_results.py        (html/docx via pandoc if available)
"""
import csv
import glob
import json
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG, TAB = RES / "figures", RES / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

R = [json.load(open(f)) for f in sorted(glob.glob(str(RES / "data" / "*.json")))]
ERR = [json.load(open(f)) for f in sorted(glob.glob(str(RES / "data" / "*.err")))]
EPS_MAIN = 1e-2
ST = json.load(open(RES / "dmrg_status.json")) if (RES / "dmrg_status.json").exists() else {}


def conv_str(r):
    """DMRG convergence of the three runs (ground / top / excited-by-penalty), from the data file or dmrg_status.json."""
    runs = r["inputs"].get("dmrg_runs") or ST.get(f"N{r['N']}_J2{r['J2']:g}", [])
    if not runs:
        return "–"
    bad = [x["role"] for x in runs if not x["converged"]]
    return "yes" if not bad else "NO (" + ", ".join(bad) + ")"


def sel(source=None, N=None, J2=None, eps=None):
    return [r for r in R if (source is None or r["source"] == source) and (N is None or r["N"] == N)
            and (J2 is None or abs(r["J2"] - J2) < 1e-9) and (eps is None or abs(r["eps"] - eps) < 1e-12)]


def one(source, N, J2, eps):
    s = sel(source, N, J2, eps)
    return s[0] if s else None


def D(r, which):
    return None if r is None else r["designs"].get(which)


def ok(d):
    return d is not None and "n" in d


def from_floor(r):
    """True if the guaranteed design came from the precision-targeted floor solver (not the fallback FilterBuilder)."""
    return r["design"]["source"] in ("floor", "none")


# ------------------------------------------------------------------ formatting
def f(x, spec=".3g", dash="–"):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return dash
    if spec == "int":
        return f"{int(round(x)):,}"
    return format(x, spec)


def _cell(x):
    """Escape pipes so |x| inside a header/cell does not split the markdown table column."""
    return str(x).replace("|", "\\|")


def md(headers, rows):
    out = ["| " + " | ".join(_cell(h) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(_cell(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def ks(k):
    return "[" + ", ".join(str(x) for x in k) + "]"


def csv_dump(name, headers, rows):
    with open(TAB / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(headers)
        w.writerows(rows)


J2_ORDER = (0.0, 0.2411, 0.4)
jl = lambda j: f"{j:g}"
tot = lambda d, topo="a2a": d[f"counts_rzz_{topo}"]["total"]
fil = lambda d, topo="a2a": d[f"counts_rzz_{topo}"]["filter"]
drc = lambda d: d["counts_rzz_direct"]["filter"]      # direct (level-1) compilation, filter only
cxl3 = lambda d: d["counts_a2a"]["filter"]            # legacy CX-basis counts (level 3), kept for the CX -> RZZ comparison table
cxdr = lambda d: d["counts_direct"]["filter"]         # legacy CX-basis counts (direct)

L = []
A = L.append

# ------------------------------------------------------------------ 0 header
Ns_large = sorted({r["N"] for r in R if r["N"] > 14})
A("# Results: resource costs of the cosine-filter pipeline for the open J1–J2 chain")
A("")
A("*Produced by `scripts/run_results.py` + `scripts/build_results.py` from the pipeline in this directory "
  "(`pipeline.py`, `floor.py`, `core/*`). Every number below is read from `results/data/*.json`; "
  "raw tables are in `results/tables/*.csv`.*")
A("")
A("## 0. How to read this document")
A("")
A("* **Pipeline.** DMRG (`get_energies`) → spectrum (`get_spectrum`) → L-layer approximate trial circuit "
  "(`mps-to-circuit`, L = 1) → filter design (`floor.py`, certified η) → snapped Trotter grid → counts. "
  "Settings are those in `pipeline.py`: first-order Lie–Trotter, ε targets as given, `SWEEP_L = 1`, `BANDWIDTH = \"certified\"` "
  "for the DMRG path, `SYNTH_EPS = 1e-3`.")
A("* **Two designs per case.** *Guaranteed* = smallest n whose a-priori rigorous bound on 1−F is ≤ ε. "
  "*Empirical* = smallest n whose **measured** 1−F (against the exact ground state) is ≤ ε. The empirical design "
  "uses the exact ground state, so it is an oracle lower envelope, not something a practitioner can know in advance; "
  "it needs simulation and is therefore only available for N ≤ 14.")
A("* **Everything uses the DMRG spectrum.** E0, E1 (orthogonality-penalty DMRG: lowest state of any symmetry orthogonal to ψ0) and E_top from DMRG; the bandwidth is the **certified** "
  "W = Σ_bonds λmax(bond) − E0 (Weyl's inequality: λmax of a sum ≤ sum of λmax; a bond J S·S has λmax = J/4), which equals E_top − E0 for J1, J2 ≥ 0 "
  "because the fully polarized state attains it (Section 6). The DMRG gap, shift and γ are not certified: every guarantee below is a *conditional certificate* (exact only if the DMRG reference is exact; see `CERTIFICATION.tex`, Section 1 and Section 12).")
A("* **N ≤ 14 versus N > 14.** For N ≤ 14 an exact-diagonalisation cross-check and statevector simulation of the filter are available: "
  "γ is the overlap with the ED ground state, and measured 1−F, P_succ and the state error are reported. "
  "For **N > 14 there is no ED, no statevector and no dense or sparse 2^N operator**: γ is the overlap of the DMRG MPS with the "
  "trial circuit simulated as an MPS (Aer matrix-product-state), the energy variance comes from MPO contractions, α is the "
  "same exact commutator sum as for small N (Table 1b), and circuits are **counted, not simulated** (no measured 1−F or P_succ; only the rigorous bound).")
A("* **Time.** `T (scaled)` is the total evolution time Σtᵢ in units of 1/W; `T (phys)` = T/W in units of 1/J1 "
  "(H = Σ J S·S, J1 = 1). Pulses with kᵢ = 0 on the grid are dropped.")
A("* **Gate basis and counting.** The two-qubit gate is **RZZ(θ) = exp(−iθ ZZ/2)**; CX is *not* in the basis (qiskit basis `rzz, rz, sx, x`; RZZ(π/2) is Clifford). "
  "One Trotter step of one pulse (every Pauli term of H, controlled by the ancilla) is transpiled once and multiplied by n; per pulse one Rz(2φ) is added. "
  "Two compilations of the same circuit are reported: **direct** (qiskit `optimization_level=1`: every Pauli rotation exp(−iθ P⊗P⊗Z) becomes a ladder of four Clifford RZZ(π/2) around one non-Clifford Rz): "
  "RZZ = 4 × (#Pauli terms) per step, non-Clifford rotations = #Pauli terms + 1 per step; "
  "**level 3** (block re-synthesis): two RZZ per Pauli term, one of them a non-Clifford RZZ(θ), plus extra non-Clifford single-qubit Rz, i.e. half the RZZ of the direct compilation but more non-Clifford rotations. "
  "`Non-Cliff. rot.` = Rz gates plus RZZ gates whose angle is not a multiple of π/2 (the cost metric for fault-tolerant hardware; no T-count is reported). Clifford single-qubit gates (Hadamards, sx, x) are not tracked. "
  "Per-step × n was checked against a full transpile of a 48-step circuit in the RZZ basis (Section 9). "
  "**Level-3 totals include the trial circuit** (L = 1 brickwork, N−1 two-qubit gates, 3 non-Clifford RZZ each); direct counts are for the filter only. Connectivity is all-to-all unless marked *line* "
  "(ancilla at the end of a line: pessimistic). Counts in the previous CX basis are kept in the data files and compared in Section 4b.")
A("* **Success path.** All counts are for the all-zeros (success) path. The expected cost per accepted run is "
  "total ÷ P_succ (early abort makes it lower).")
A("")
A(f"Data present: {len(R)} (N, J2, ε) points, N from {min(r['N'] for r in R)} to {max(r['N'] for r in R)}; {len(ERR)} failed point(s).")
A("")

# ------------------------------------------------------------------ 1 inputs
A("## 1. Case inputs (DMRG spectrum, certified bandwidth)")
A("")
rows, hd = [], ["N", "J2", "E0 (DMRG)", "E1", "E_top", "|E0−E0(ED)|", "|E1−E1(ED)|", "E_top−E0", "W (certified)", "W / (E_top−E0)", "Δ (scaled gap)", "γ (trial)", "γ source",
                "trial-circuit infidelity 1−γ", "α used", "Pauli terms", "trial RZZ", "trial Rz", "trial non-Cliff. rot.", "trial depth", "E1 method", "DMRG converged"]
for J2 in J2_ORDER:
    for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
        r = one("dmrg", N, J2, EPS_MAIN)
        if r is None:
            continue
        i = r["inputs"]
        tc = i["trial_cost_rzz"]
        ed_ = i.get("dmrg_vs_ed") or {}
        rows.append([N, jl(J2), f(i["E0"], ".10f"), f(i["E1"], ".10f"), f(i["Etop"], ".10f"),
                     f(ed_.get("E0_diff"), ".1e") if ed_ else "n/a (no ED)", f(ed_.get("E1_diff"), ".1e") if ed_ else "n/a (no ED)",
                     f(i["W_exact_dmrg"], ".4f"), f(i["W"], ".4f"),
                     f(i["W"] / i["W_exact_dmrg"], ".2f"), f(i["gap_scaled"], ".5f"), f(i["gamma"], ".5f"), i["gamma_src"],
                     f(1 - i["gamma"], ".2e"), f(i["alpha"], ".4f"), i["n_pauli_terms"], f(tc["rzz"], "int"), f(tc["rz"], "int"),
                     f(tc["nc_total"], "int"), f(tc["depth"], "int"),
                     i.get("e1_method") or "–", conv_str(r)])
A(md(hd, rows))
csv_dump("inputs.csv", hd, rows)
A("")
A("*Table 1.* The energies are DMRG values (10 digits); `|E0−E0(ED)|`, `|E1−E1(ED)|` are the DMRG errors against exact diagonalisation (N ≤ 14 here; N = 16–20 in Section 5; N ≥ 24 has no exact reference, see the run-to-run table in Section 5). "
  "**1−γ is not a DMRG error:** it is the infidelity of the low-depth (L = 1 brickwork) trial *circuit* with respect to the ground state, which is large by construction and is what the filter has to remove. "
  "*DMRG converged*: the ground-state, top-state and orthogonality-penalty excited-state DMRG runs each reached the sweep tolerance within 100 sweeps; "
  "a **NO** means that quantity (for the excited state, E1 and hence the gap Δ) is not converged and the design is only indicative. "
  "*E1 method*: `local` = projector penalty λ|ψ0⟩⟨ψ0| applied as a rank-1 term in every two-site local problem, giving the lowest state of any symmetry (checked against ED in Section 5). "
  "Trial circuit = L = 1 approximate brickwork circuit from the DMRG MPS. α is in scaled units: the exact commutator sum Σ_g‖[H_g, Σ_{g'>g}H_g']‖ for every N (computed on the local support of each term). W is the certified (Weyl) bandwidth, which coincides with E_top − E0 here, so Δ = (E1−E0)/W is the exact scaled gap.")
A("")
rows, hd = [], ["N", "J2", "α (exact commutator sum)", "α triangle (diagnostic)", "triangle / exact"]
for J2 in J2_ORDER:
    for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
        r = one("dmrg", N, J2, EPS_MAIN)
        if r is not None and "alpha_tri" in r["inputs"]:
            i = r["inputs"]
            rows.append([N, jl(J2), f(i["alpha"], ".5f"), f(i["alpha_tri"], ".5f"), f(i["alpha_tri"] / i["alpha"], ".3f")])
if rows:
    A("**α: exact commutator sum against the triangle bound**")
    A("")
    A(md(hd, rows))
    csv_dump("alpha_check.csv", hd, rows)
    A("")
    A("*Table 1b.* α is the exact value of Σ_g‖[H_g, Σ_{g'>g}H_g']‖ for every N (for Pauli terms it reduces to a small dense eigenvalue problem on the support of each term). "
      "The triangle-inequality value, which earlier versions used for N > 14, is shown for comparison: it is exact at J2 = 0 and up to about 1.25× larger otherwise.")
    A("")


def design_tables(source, which, N_list, title_prefix, eps=EPS_MAIN, tag=""):
    """Three tables: design/time, resource counts, measured quality."""
    out = []
    rows1, rows2, rows3 = [], [], []
    h1 = ["N", "J2", "pulses", "n (steps)", "kᵢ", "T (scaled)", "T (phys)", "η cert.", "ε bound", "P_succ lb"]
    h2 = ["N", "J2", "RZZ trial", "RZZ filter", "RZZ total", "RZZ filter (direct)", "Non-Cliff. rot. (level 3)", "Non-Cliff. rot. (direct)", "depth total",
          "E[RZZ]/accept", "RZZ total (line)"]
    h3 = ["N", "J2", "measured 1−F", "F (ED)", "P_succ meas.", "state err vs exact", "bound holds", "full-transpile RZZ (n ≤ 200)"]
    for J2 in J2_ORDER:
        for N in N_list:
            r = one(source, N, J2, eps)
            d = D(r, which)
            if not ok(d):
                continue
            c = tot(d)
            ps = d["p_succ"] if d.get("simulated") else d["p_succ_lb"]
            rows1.append([N, jl(J2), d["n_pulses_used"], f(d["n"], "int"), ks(d["k"]), f(d["T_scaled"], ".2f"),
                          f(d["T_phys"], ".2f"), f(d["eta"], ".3f"), f(d["eps_bound"], ".2e"), f(d["p_succ_lb"], ".3f")])
            rows2.append([N, jl(J2), f(r["inputs"]["trial_cost_rzz"]["rzz"], "int"), f(fil(d)["rzz"], "int"), f(c["rzz"], "int"),
                          f(drc(d)["rzz"], "int"), f(c["nc_total"], "int"), f(drc(d)["nc_total"], "int"), f(c["depth"], "int"),
                          f(c["rzz"] / ps, "int"), f(tot(d, "line")["rzz"], "int")])
            ft = d.get("full_transpile_check")
            rows3.append([N, jl(J2), f(d.get("eps_actual"), ".2e"), f(d.get("F_ed"), ".5f"), f(d.get("p_succ"), ".3f"),
                          f(d.get("state_err"), ".2e"),
                          ("yes" if d["eps_actual"] <= d["eps_bound"] + 1e-9 else "NO") if d.get("simulated") else "not simulated",
                          f(ft["rzz"], "int") if ft and "rzz" in ft else "–"])
    if rows1:
        out += [f"**{title_prefix}: design and evolution time**", "", md(h1, rows1), "",
                f"**{title_prefix}: gate counts (level-3 columns include the trial circuit; direct columns are filter only)**", "", md(h2, rows2), "",
                f"**{title_prefix}: measured quality (noiseless)**", "", md(h3, rows3), ""]
        csv_dump(f"{source}_{which}_{tag}eps{eps:g}_design.csv", h1, rows1)
        csv_dump(f"{source}_{which}_{tag}eps{eps:g}_counts.csv", h2, rows2)
        csv_dump(f"{source}_{which}_{tag}eps{eps:g}_measured.csv", h3, rows3)
    return out


A("## 2. Guaranteed designs, N ≤ 14, ε = 10⁻²")
A("")
L += design_tables("dmrg", "guaranteed", (4, 6, 8, 10, 12, 14), "Guaranteed")
A("*Tables 2a–c.* Rigorous a-priori n: the bound on 1−F is ≤ ε. Measured 1−F is typically far below ε (the bound is conservative).")
A("")

A("## 3. Empirical (oracle) designs, N ≤ 14, ε = 10⁻²")
A("")
L += design_tables("dmrg", "empirical", (4, 6, 8, 10, 12, 14), "Empirical")
A("*Tables 3a–c.* Smallest n on the sweep whose *measured* 1−F ≤ ε; an oracle (uses the exact ground state). "
  "At very small n the first-order Trotter error can be accidentally small, so the empirical n is a lower envelope, "
  "not a design a practitioner could trust.")
A("")

# ------------------------------------------------------------------ 4 per step
A("## 4. Cost of one Trotter step")
A("")
rows, hd = [], ["N", "J2", "bonds", "Pauli terms", "RZZ (level 3)", "of which non-Cliff.", "Non-Cliff. rot. (level 3)", "RZZ (direct)", "Non-Cliff. rot. (direct)", "depth (level 3)", "RZZ (line, 3-step avg)", "line / a2a"]
for J2 in J2_ORDER:
    for N in sorted({r["N"] for r in R}):
        r = (sel("dmrg", N, J2, EPS_MAIN) or sel("dmrg", N, J2) or [None])[0]
        if r is None:
            continue
        s1, sl = r["per_step"]["rzz"]["a2a"]["one"], r["per_step"]["rzz"]["line"]["avg3"]
        rows.append([N, jl(J2), (N - 1) + ((N - 2) if J2 else 0), r["inputs"]["n_pauli_terms"], f(s1["rzz"], "int"),
                     f(s1["rzz_nonclifford"], "int"), f(s1["nc_total"], "int"), f(r["per_step"]["rzz"]["direct"]["rzz"], "int"), f(r["per_step"]["rzz"]["direct"]["nc_total"], "int"), f(s1["depth"], "int"),
                     f(sl["rzz"], ".1f"), f(sl["rzz"] / s1["rzz"], ".2f")])
A(md(hd, rows))
csv_dump("per_step.csv", hd, rows)
A("")
A("*Table 4.* Transpiled counts of one Trotter step of one pulse (all bonds, ancilla-controlled), all-to-all, plus the "
  "line-connectivity cost. RZZ per step grows linearly with the number of bonds; J2 ≠ 0 roughly doubles it (J2 = 0.2411 and 0.4 have identical circuits per step, only the angles differ). "
  "The line-connectivity figure depends on the transpiler's routing, which is not smooth in N (e.g. N = 24, J2 = 0.2411); treat it as indicative.")
A("")

# ------------------------------------------------------------------ 4b CX vs RZZ
A("### 4b. Previous CX basis versus the RZZ basis (guaranteed designs, ε = 10⁻², totals including the trial circuit for level 3)")
A("")
rows, hd = [], ["N", "J2", "n (steps)", "CX total (level 3)", "RZZ total (level 3)", "RZZ / CX (level 3)", "CX filter (direct)", "RZZ filter (direct)", "RZZ / CX (direct)",
                "Rz non-Cliff. (level 3, CX basis)", "Non-Cliff. rot. (level 3, RZZ basis)"]
for J2 in J2_ORDER:
    for N in sorted({r["N"] for r in R}):
        r = one("dmrg", N, J2, EPS_MAIN)
        d = D(r, "guaranteed")
        if not ok(d):
            continue
        c3 = d["counts_a2a"]["total"]; z3 = tot(d)
        rows.append([N, jl(J2), f(d["n"], "int"), f(c3["cx"], "int"), f(z3["rzz"], "int"), f(z3["rzz"] / c3["cx"], ".2f"),
                     f(cxdr(d)["cx"], "int"), f(drc(d)["rzz"], "int"), f(drc(d)["rzz"] / cxdr(d)["cx"], ".2f"),
                     f(c3["rz_nonclifford"], "int"), f(z3["nc_total"], "int")])
A(md(hd, rows))
csv_dump("cx_vs_rzz.csv", hd, rows)
A("")
A("*Table 4b.* The direct compilation maps each CX one-to-one onto a Clifford RZZ(π/2) (ratio 1.00). At level 3 qiskit had already re-synthesised the CX circuit to about 2.3 CX per ancilla-controlled Pauli rotation; "
  "in the RZZ basis each such rotation costs exactly two RZZ (one Clifford, one non-Clifford RZZ(θ)), so the RZZ count is about 0.86× the level-3 CX count (and half of the direct count). "
  "In the RZZ basis the non-Clifford count includes the RZZ(θ) rotations, so `Non-Cliff. rot.` is not directly comparable with the Rz-only count of the CX basis.")
A("")

# ------------------------------------------------------------------ 5 certified bandwidth
A("## 5. Accuracy of the DMRG inputs")
A("")
A("The certificate depends on E0, E1 (through Δ), the DMRG cross-check of E_top, and γ, all of which come from DMRG (`get_energies.py`). "
  "The settings, chosen so that two-site DMRG reaches roughly machine precision on these 1D chains:")
A("")
A(md(["setting", "value", "role"], [
    ["bond dimensions", "8, 16, 32, 64, 128, 256 (one per sweep stage)", "the cutoff, not this cap, limits the bond dimension here (χ reached: 40–100 for N = 12–24)"],
    ["SVD cutoff", "1e-14", "discarded weight per bond ≈ energy error in 1D"],
    ["energy tolerance", "1e-12", "stop when the energy changes by less between sweeps"],
    ["inner Lanczos tolerance", "1e-10 (retry ×100, ×1e4 if ARPACK fails)", "quimb's default is 1e-3"],
    ["E1 (all N)", "H + 1.1·bandwidth·|ψ0⟩⟨ψ0|, projector applied as a rank-1 term in each two-site local solve", "lowest state of any symmetry; avoids the χ²-bond projector MPO, which ran out of memory at N ≥ 16 with these tolerances"],
    ["max sweeps", "100", "convergence flags are recorded in Table 1"]]))
A("")
rows, hd = [], ["N", "J2", "|E0 − E0(ED)|", "|E1 − E1(ED)|", "|E_top − E_top(ED)|", "|ΔE1| / gap", "1 − |⟨ψ0|ψ0(ED)⟩|²", "E1 method", "converged", "max bond"]
VAL = json.load(open(RES / "dmrg_validation.json")) if (RES / "dmrg_validation.json").exists() else {}
for J2 in J2_ORDER:
    for N in (4, 6, 8, 10, 12, 14):
        r = one("dmrg", N, J2, EPS_MAIN)
        if r is None or not r["inputs"].get("dmrg_vs_ed"):
            continue
        i = r["inputs"]; e = i["dmrg_vs_ed"]
        rows.append([N, jl(J2), f(e.get("E0_diff"), ".1e"), f(e.get("E1_diff"), ".1e"), f(e.get("E_top_diff"), ".1e"),
                     f(e.get("E1_diff", 0) / i["gap_raw"], ".1e"), f(1 - i["dmrg_ed_overlap"], ".1e"), i.get("e1_method") or "–",
                     conv_str(r), max(x["max_bond"] for x in i["dmrg_runs"]) if i.get("dmrg_runs") else "–"])
for key in sorted((k for k in VAL if k.startswith("ed_")), key=lambda k: (VAL[k]["J2"], VAL[k]["N"])):
    v = VAL[key]
    rows.append([v["N"], jl(v["J2"]), f(v["dE0"], ".1e"), f(v["dE1"], ".1e"), f(v["dEtop"], ".1e"), f(v["gap_rel"], ".1e"), f(v["one_minus_overlap"], ".1e"),
                 v["e1_method"], "yes" if v["converged"] else "NO", max(v["max_bond"])])
if rows:
    A("**DMRG against exact diagonalisation** (N ≤ 14 from the main runs; N ≥ 16 from `scripts/validate_dmrg.py ed`, exact reference used for testing only)")
    A("")
    A(md(hd, rows))
    csv_dump("dmrg_vs_ed.csv", hd, rows)
    A("")
    A("*Table 5a.* Absolute errors of the DMRG energies, the relative error of the gap (|ΔE1|/gap for N ≤ 14, where the E0 error is negligible next to it), "
      "and the overlap deficit of the ground state, against exact diagonalisation.")
    A("")
cons = sorted((k for k in VAL if k.startswith("consist_")), key=lambda k: (VAL[k]["J2"], VAL[k]["N"]))
if cons:
    rows, hd = [], ["N", "J2", "E0", "E1", "run-to-run |ΔE0|", "|ΔE1|", "|ΔE_top|", "1 − overlap(ψ0, ψ0′)", "variance of ψ0 (a / b)", "variance of ψ1 (a / b)", "converged", "max bond"]
    for key in cons:
        v = VAL[key]
        rows.append([v["N"], jl(v["J2"]), f(v["E0"], ".8f"), f(v["E1"], ".8f"), f(v["dE0"], ".1e"), f(v["dE1"], ".1e"), f(v["dEtop"], ".1e"),
                     f(v["one_minus_overlap"], ".1e"), f"{v['var0']:.1e} / {v['var0_b']:.1e}", f"{v['var1']:.1e} / {v['var1_b']:.1e}",
                     "yes" if v["converged"] else "NO", max(v["max_bond"])])
    A("**No exact reference: two independent runs from different random starts (N ≥ 24)**")
    A("")
    A(md(hd, rows))
    csv_dump("dmrg_consistency.csv", hd, rows)
    A("")
    A("*Table 5b.* Agreement between two independent DMRG runs and the energy variance ⟨H²⟩ − ⟨H⟩² of each state. A state with variance s² has an exact eigenvalue within s "
      "of its energy; run-to-run agreement at 1e-13 or better together with variances of order 1e-12 is strong evidence of convergence, but it is not a proof that E1 is the second-lowest level.")
    A("")
gm = sorted((k for k in VAL if k.startswith("gamma_")), key=lambda k: (VAL[k]["J2"], VAL[k]["N"]))
if gm:
    rows, hd = [], ["N", "J2", "γ, MPS contraction (no ED, no statevector)", "γ, dense overlap with exact ground state", "difference"]
    for key in gm:
        v = VAL[key]
        rows.append([v["N"], jl(v["J2"]), f(v["gamma_mps_contraction"], ".12f"), f(v["gamma_ed_dense"], ".12f"), f"{v['diff']:+.1e}"])
    A("**The ground-state weight γ used for N > 14 against the dense value**")
    A("")
    A(md(hd, rows))
    csv_dump("gamma_validation.csv", hd, rows)
    A("")
    A("*Table 5c.* γ is the overlap of the DMRG ground state with the trial circuit. The N > 14 path computes it by MPS contraction; the dense value is for checking only.")
    A("")

A("## 6. The certified bandwidth")
A("")
A("The scaling H_s = (H − E0)/W needs W ≥ E_top − E0 so that the spectrum stays inside [0, 1]. A certified W must not rely on DMRG's E_top "
  "(DMRG on −H is variational, so E_top(DMRG) ≤ E_top and W(DMRG) ≤ W_exact: the unsafe direction). Weyl's inequality gives "
  "λmax(Σ_b h_b) ≤ Σ_b λmax(h_b); a bond J S·S has λmax = J/4 (J ≥ 0), so E_top ≤ (N−1)J1/4 + (N−2)J2/4, and the fully polarized state attains it. "
  "Hence W = (N−1)J1/4 + (N−2)J2/4 − E0 equals E_top − E0 (to rounding) with no DMRG input for E_top. "
  "An earlier version used the looser Pauli bound W = Σ|c_β| − E0 = 3(Σ_b λmax) − E0, which is 1.6–1.9× larger; Table 5 compares them.")
A("")
rows, hd = [], ["N", "J2", "W (Weyl) = E_top − E0", "|W − (E_top−E0)_DMRG|", "W (Pauli-ℓ1, earlier)", "ℓ1 / Weyl"]
for J2 in J2_ORDER:
    for N in (4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32):
        r = one("dmrg", N, J2, EPS_MAIN)
        if r is None:
            continue
        i = r["inputs"]
        wl1 = 3.0 * (i["W"] + i["E0"]) - i["E0"]
        rows.append([N, jl(J2), f(i["W"], ".6f"), f(abs(i["W"] - i["W_exact_dmrg"]), ".1e"), f(wl1, ".4f"), f(wl1 / i["W"], ".2f")])
A(md(hd, rows))
csv_dump("bandwidth.csv", hd, rows)
A("")
BWC = json.load(open(RES / "bandwidth_comparison.json")) if (RES / "bandwidth_comparison.json").exists() else None
A("*Table 5.* The Weyl bandwidth agrees with the DMRG E_top − E0 to rounding (≤ 1e-13), because E_top is the all-up energy. "
  "**W cancels in the step count:** α_scaled = α/W² and T_scaled = xπW/g (g the unscaled gap), so n ∝ α_scaled T² ∝ α/g² and T (phys) = T_scaled/W = xπ/g do not depend on W. "
  "W only changes the *window*: [Δ, 1] is wider than the occupied spectrum when W is loose, which can make the certified η slightly larger.")
if BWC:
    A("")
    A(f"Measured effect of switching from the Pauli-ℓ1 bandwidth to the Weyl bandwidth, over {BWC['count']} matching (N, J2, ε) guaranteed designs: "
      f"n changed by a factor {BWC['n_ratio_min']:.2f}–{BWC['n_ratio_max']:.2f} (median {BWC['n_ratio_median']:.2f}), T (phys) by {BWC['T_ratio_min']:.2f}–{BWC['T_ratio_max']:.2f}.")
A("")
# ------------------------------------------------------------------ 6 larger systems
A(f"## 7. Larger J1–J2 chains (N = 14–{max(Ns_large) if Ns_large else 14}): resource costs")
A("")
A("N ≥ 16 is **DMRG/MPS only**: γ is the DMRG-MPS overlap of the trial circuit, the gap is DMRG's, W is certified, "
  "α is the exact commutator sum, and there is no ED, no statevector and no dense operator anywhere. The designs are therefore conditional certificates: rigorous *given* γ, the gap "
  "and the shift (see `CERTIFICATION.tex`). Circuits are **not simulated**; their counts are n × the transpiled per-step cost plus the trial circuit, validated "
  "against full transpiles at n = 48 (Section 9). N = 14 (statevector-simulated, ED cross-checked) is included for continuity.")
A("")
Nl = [N for N in (14, 16, 18, 20, 24, 28, 32) if sel("dmrg", N)]
for eps in (EPS_MAIN, 1e-1):
    rows, hd = [], ["N", "J2", "γ", "Δ", "W", "α", "pulses", "n (steps)", "T (scaled)", "T (phys)", "RZZ trial", "RZZ filter", "RZZ total",
                    "RZZ filter (direct)", "Non-Cliff. rot. (level 3)", "Non-Cliff. rot. (direct)", "depth", "P_succ lb", "ε bound", "RZZ total (line)", "DMRG converged"]
    for J2 in J2_ORDER:
        for N in Nl:
            r = one("dmrg", N, J2, eps)
            d = D(r, "guaranteed")
            if not ok(d):
                continue
            i, c = r["inputs"], tot(d)
            rows.append([N, jl(J2), f(i["gamma"], ".4f"), f(i["gap_scaled"], ".4f"), f(i["W"], ".2f"), f(i["alpha"], ".4f"),
                         d["n_pulses_used"], f(d["n"], "int") + ("" if from_floor(r) else " †"), f(d["T_scaled"], ".1f"), f(d["T_phys"], ".2f"),
                         f(i["trial_cost_rzz"]["rzz"], "int"), f(fil(d)["rzz"], "int"), f(c["rzz"], "int"), f(drc(d)["rzz"], "int"),
                         f(c["nc_total"], "int"), f(drc(d)["nc_total"], "int"), f(c["depth"], "int"),
                         f(d["p_succ_lb"], ".3f"), f(d["eps_bound"], ".2e"), f(tot(d, "line")["rzz"], "int"), conv_str(r)])
    if rows:
        A(f"**Guaranteed designs, ε = {eps:g}**")
        A("")
        A(md(hd, rows))
        csv_dump(f"large_guaranteed_eps{eps:g}.csv", hd, rows)
        A("")
        if any(" †" in str(row[7]) for row in rows):
            A("† The floor solver (≤ 8 pulses, t_max = T/2) found no feasible filter for this γ and ε, so the design fell back to the generic minimax FilterBuilder "
              "(up to 15 pulses, x not optimized). It is a valid certified design but far from the trend; do not read scaling from it.")
            A("")
rows, hd = [], ["N", "J2", "kᵢ (guaranteed)", "largest block kmax", "RZZ in largest block"]
for J2 in J2_ORDER:
    for N in Nl:
        r = one("dmrg", N, J2, EPS_MAIN)
        d = D(r, "guaranteed")
        if ok(d):
            rows.append([N, jl(J2), ks(d["k"]), d["kmax"], f(d["kmax"] * r["per_step"]["rzz"]["a2a"]["one"]["rzz"], "int")])
if rows:
    A("**Pulse structure of the guaranteed designs (ε = 10⁻²)**")
    A("")
    A(md(hd, rows))
    A("")

# empirical at N = 14
rows, hd = [], ["N", "J2", "design", "n", "pulses", "T (phys)", "RZZ total", "Non-Cliff. rot. (direct)", "measured 1−F", "P_succ", "ε bound (rigorous)"]
for J2 in J2_ORDER:
    for N in (12, 14):
        for src in ("dmrg",):
            r = one(src, N, J2, EPS_MAIN)
            for which in ("guaranteed", "empirical"):
                d = D(r, which)
                if ok(d) and d.get("simulated"):
                    c = tot(d)
                    rows.append([N, jl(J2), which, f(d["n"], "int"), d["n_pulses_used"], f(d["T_phys"], ".2f"), f(c["rzz"], "int"),
                                 f(drc(d)["nc_total"], "int"), f(d["eps_actual"], ".2e"), f(d["p_succ"], ".3f"), f(d["eps_bound"], ".2e")])
if rows:
    A("**Measured (statevector-simulated) points, N = 12 and 14**")
    A("")
    A(md(hd, rows))
    A("")
    A("*Empirical* rows are the oracle lower envelope; *guaranteed* rows show how far below ε the measured error sits.")
    A("")

# scaling fits
A("### Scaling with N")
A("")


def fit(xs, ys):
    xs, ys = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    if len(xs) < 3:
        return None
    p, c = np.polyfit(xs, ys, 1)
    pred = p * xs + c
    ss = float(np.sum((ys - ys.mean()) ** 2))
    return p, (1 - float(np.sum((ys - pred) ** 2)) / ss) if ss > 0 else 1.0


rows, hd = [], ["quantity", "J2", "source", "N range", "p in y ~ N^p", "R²", "points"]
FITS = {}
metrics = [("n (steps)", lambda r, d: d["n"]), ("T (scaled)", lambda r, d: d["T_scaled"]), ("RZZ per step", lambda r, d: r["per_step"]["rzz"]["a2a"]["one"]["rzz"]),
           ("RZZ total", lambda r, d: tot(d)["rzz"]), ("RZZ (direct)", lambda r, d: drc(d)["rzz"]),
           ("Non-Cliff. rot. (level 3)", lambda r, d: tot(d)["nc_total"]), ("Non-Cliff. rot. (direct)", lambda r, d: drc(d)["nc_total"]),
           ("depth total", lambda r, d: tot(d)["depth"])]
for src, Nrange in (("dmrg", (6, 32)),):
    for J2 in J2_ORDER:
        pts = [(N, one(src, N, J2, EPS_MAIN)) for N in range(Nrange[0], Nrange[1] + 1, 2)]
        pts = [(N, r, D(r, "guaranteed")) for N, r in pts if r is not None and ok(D(r, "guaranteed")) and r["inputs"].get("dmrg_converged", True) and from_floor(r)]
        for name, fn in metrics:
            if len(pts) >= 3:
                fr = fit([N for N, _, _ in pts], [fn(r, d) for N, r, d in pts])
                FITS[(name, J2, src)] = (fr, pts)
                rows.append([name, jl(J2), src, f"{pts[0][0]}–{pts[-1][0]}", f(fr[0], ".2f"), f(fr[1], ".3f"), len(pts)])
if rows:
    A(md(hd, rows))
    csv_dump("scaling_fits.csv", hd, rows)
    A("")
    A("*Table 7.* Power-law fits of guaranteed-design quantities (ε = 10⁻²) over points with converged DMRG and a certified design from the floor solver (fallback-builder points, marked † in Section 7, are excluded). Few points: indicative, not asymptotic.")
    A("")

# ------------------------------------------------------------------ 7 eps scan
A("## 8. Dependence on the target ε")
A("")
for (N, J2) in ((6, 0.0), (8, 0.0), (6, 0.4)):
    rows, hd = [], ["ε", "design", "n", "pulses", "T (phys)", "RZZ total", "Non-Cliff. rot. (direct)", "depth", "measured 1−F", "ε bound", "P_succ"]
    for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
        r = one("dmrg", N, J2, eps)
        for which in ("guaranteed", "empirical"):
            d = D(r, which)
            if ok(d):
                c = tot(d)
                rows.append([f(eps, ".0e"), which, f(d["n"], "int"), d["n_pulses_used"], f(d["T_phys"], ".2f"), f(c["rzz"], "int"),
                             f(drc(d)["nc_total"], "int"), f(c["depth"], "int"), f(d.get("eps_actual"), ".2e"), f(d["eps_bound"], ".2e"),
                             f(d.get("p_succ", d["p_succ_lb"]), ".3f")])
    if rows:
        A(f"**N = {N}, J2 = {J2:g}**")
        A("")
        A(md(hd, rows))
        csv_dump(f"eps_scan_N{N}_J2{J2:g}.csv", hd, rows)
        A("")

# ------------------------------------------------------------------ 8 validation + timing
A("## 9. Checks and classical compute cost")
A("")
rows, hd = [], ["source", "N", "J2", "n", "compilation", "estimate RZZ (per-step × n)", "full-transpile RZZ", "estimate non-Cliff. rot.", "full-transpile non-Cliff. rot.", "estimate depth", "full-transpile depth"]
for r in sorted(R, key=lambda r: (r["source"], r["N"], r["J2"])):
    v = r.get("validation_rzz")
    if v and "l3" in v and abs(r["eps"] - EPS_MAIN) < 1e-12:
        for tag, name in (("l3", "level 3"), ("direct", "direct")):
            x = v[tag]
            rows.append([r["source"], r["N"], jl(r["J2"]), v["n"], name, f(x["est_rzz"], "int"), f(x["full_rzz"], "int"),
                         f(x["est_nc"], "int"), f(x["full_nc"], "int"), f(x["est_depth"], "int"), f(x["full_depth"], "int")])
if rows:
    A("**Per-step × n estimate against a full transpile of the real circuit in the RZZ basis (4 pulses × 12 steps = 48 steps, all-to-all)**")
    A("")
    A(md(hd, rows))
    csv_dump("validation.csv", hd, rows)
    A("")
    vr = [r["validation_rzz"] for r in R if r.get("validation_rzz") and "l3" in r["validation_rzz"]]
    if vr:
        dz3 = max(abs(v["l3"]["est_rzz"] - v["l3"]["full_rzz"]) for v in vr); dzd = max(abs(v["direct"]["est_rzz"] - v["direct"]["full_rzz"]) for v in vr)
        dnd = max(abs(v["direct"]["est_nc"] - v["direct"]["full_nc"]) for v in vr)
        o3 = [(v["l3"]["est_nc"] - v["l3"]["full_nc"]) / v["l3"]["full_nc"] for v in vr]
        d3 = [(v["l3"]["est_depth"] - v["l3"]["full_depth"]) / v["l3"]["full_depth"] for v in vr]
        A(f"*Summary ({len(vr)} cases):* maximum |estimate − full transpile| = {dz3:.0f} RZZ (level 3) and {dzd:.0f} RZZ (direct); direct non-Clifford rotations differ by at most {dnd:.0f}; "
          f"the level-3 non-Clifford-rotation estimate is off by {100 * min(o3):+.0f}% to {100 * max(o3):+.0f}% (not conservative: block re-synthesis in a full circuit differs from the single-step transpile, so treat level-3 non-Clifford counts as ±20% estimates; the direct counts are exact), "
          f"and per-step depth × n is off by {100 * min(d3):+.0f}% to {100 * max(d3):+.0f}%.")
        A("")
rows, hd = [], ["source", "N", "J2", "ε", "context build (DMRG+ED+trial+per-step) s", "filter design s", "guaranteed search s", "empirical search s"]
for r in sorted(R, key=lambda r: (r["source"], r["N"], r["J2"], -r["eps"])):
    rows.append([r["source"], r["N"], jl(r["J2"]), f(r["eps"], ".0e"), f(r["inputs"]["build_s"], ".0f"), f(r["design"]["design_s"], ".0f"),
                 f((D(r, "guaranteed") or {}).get("search_s"), ".0f"), f((D(r, "empirical") or {}).get("search_s"), ".0f")])
if rows:
    A("**Classical compute (single core)**")
    A("")
    A(md(hd, rows))
    csv_dump("timing.csv", hd, rows)
    A("")
viol = [(r["source"], r["N"], r["J2"], r["eps"], w) for r in R for w, d in r["designs"].items()
        if ok(d) and d.get("simulated") and d["eps_actual"] > d["eps_bound"] + 1e-9]
A(f"Bound violations among simulated points: **{len(viol)}**" + (f" — {viol}" if viol else "") + ".")
A("")
errs = [(r["source"], r["N"], r["J2"], r["eps"], w, d["error"] + (f" [design source: {r['design']['source']}; the floor solver (<= 8 pulses, t_max_frac 0.5) found no feasible filter for gamma = {r['inputs']['gamma']:.3f}]" if r["design"]["source"] == "builder" else "")) for r in R for w, d in r["designs"].items() if not ok(d)]
if errs or ERR:
    A("**Points without a result**")
    A("")
    for e in errs:
        A(f"* {e[0]} N={e[1]} J2={e[2]:g} ε={e[3]:g} {e[4]}: {e[5]}")
    for e in ERR:
        A(f"* {e['source']} N={e['N']} J2={e['J2']:g} ε={e['eps']:g}: {e['error']}")
    A("")

# ------------------------------------------------------------------ 10 first vs second order Trotter
SO = {}
for fpath in sorted(glob.glob(str(ROOT / "trotter2" / "circuit_N*_J2*.json"))):
    c2 = json.load(open(fpath))
    c1p = ROOT / "trotter2" / f"first_{c2['N']}_{c2['J2']:g}.json"
    if c1p.exists():
        SO[(c2["N"], round(c2["J2"], 4))] = (json.load(open(c1p)), c2)
if SO:
    A("## 10. First-order versus second-order (Strang) Trotter formula, ε = 10⁻²")
    A("")
    A("Same classical front end, same filter design (pulse times and phases of the *guaranteed* design in Sections 2–3), same trial state and the same certificate "
      "structure (leakage and p_g unchanged); only the product formula and the matching Trotter-error bound change. Both bounds are from Childs–Su–Tran–Wiebe–Zhu, "
      "arXiv:1912.08854v3: first order = Proposition 15, Eq. (145), error ≤ (t²/2)·Σ_γ‖[Σ_{γ'>γ}H_γ', H_γ]‖ = αt²/2 (the α used everywhere in this document); "
      "second order = Proposition 16, Eq. (152), error ≤ t³·B₂ with B₂ = (1/12)Σ_γ‖[B_γ,[B_γ,H_γ]]‖ + (1/24)Σ_γ‖[H_γ,[H_γ,B_γ]]‖, B_γ = Σ_{γ'>γ}H_γ' "
      "(exact local evaluation, maximum over the forward and reversed term order; `trotter2/alpha2.py`). With k steps per pulse the Trotter error of a pulse of length t is "
      "αt²/(2k) (first order) or B₂t³/k² (second order). **Per step the second-order circuit is more expensive** (palindromic: every Pauli rotation twice, minus merging of adjacent "
      "half-steps), **but it needs far fewer steps**; the tables give both and the net effect, with all counts from a full transpile of the actual circuit "
      "(RZZ basis {rzz, rz, sx, x}, filter only, all-to-all, early-abort form, no trial circuit).")
    A("")
    rows, hd = [], ["N", "J2", "α (1st)", "B₂ (2nd)", "n 1st (certified)", "n 2nd (certified)", "n₁/n₂", "n 1st (oracle)", "n 2nd (oracle)",
                    "certified/oracle 1st", "certified/oracle 2nd", "bound held (2nd, simulated)"]
    for (N, J2), (c1, c2) in sorted(SO.items()):
        rr = one("dmrg", N, J2, EPS_MAIN); a1 = rr["inputs"]["alpha"]
        rows.append([N, jl(J2), f(a1, ".4f"), f(c2["B2"], ".5f"), f(c2["n1"], "int"), f(c2["n2"], "int"), f(c2["n1"] / c2["n2"], ".0f"),
                     f(c2["empirical1"]["n"], "int"), f(c2["empirical2"]["n"], "int"), f(c2["n1"] / c2["empirical1"]["n"], ".0f") + "×",
                     f(c2["n2"] / c2["empirical2"]["n"], ".0f") + "×", str(c2["bound_holds"])])
    A("**Steps and certificate** (oracle = smallest uniform k-scaling whose *measured* 1−F ≤ ε, found by bisection against the exact ground state; "
      "the first-order oracle is the one of Section 3)")
    A("")
    A(md(hd, rows)); csv_dump("first_vs_second_steps.csv", hd, rows); A("")
    rows, hd = [], ["N", "J2", "RZZ/step 1st", "RZZ/step 2nd", "per-step ratio", "n 1st", "n 2nd", "total RZZ 1st (L3)", "total RZZ 2nd (L3)", "2nd/1st",
                    "total RZZ 1st (direct)", "total RZZ 2nd (direct)", "2nd/1st (direct)", "non-Cliff. rot. 1st (L3)", "non-Cliff. rot. 2nd (L3)",
                    "depth 1st (L3)", "depth 2nd (L3)", "depth 2nd/1st"]
    for (N, J2), (c1, c2) in sorted(SO.items()):
        s1, s2 = c1["l3"]["rzz"] / c1["n"], c2["second_l3"]["rzz"] / c2["n2"]
        rows.append([N, jl(J2), f(s1, ".1f"), f(s2, ".1f"), f(s2 / s1, ".2f"), f(c1["n"], "int"), f(c2["n2"], "int"),
                     f(c1["l3"]["rzz"], "int"), f(c2["second_l3"]["rzz"], "int"), f(c2["second_l3"]["rzz"] / c1["l3"]["rzz"], ".3f"),
                     f(c1["direct"]["rzz"], "int"), f(c2["second_direct"]["rzz"], "int"), f(c2["second_direct"]["rzz"] / c1["direct"]["rzz"], ".3f"),
                     f(c1["l3"]["nc_total"], "int"), f(c2["second_l3"]["nc_total"], "int"),
                     f(c1["l3"]["depth"], "int"), f(c2["second_l3"]["depth"], "int"), f(c2["second_l3"]["depth"] / c1["l3"]["depth"], ".3f")])
    A("**Gate counts, certified designs (full transpile of both circuits)**")
    A("")
    A(md(hd, rows)); csv_dump("first_vs_second_gates.csv", hd, rows); A("")
    rr_ = [c2["second_l3"]["rzz"] / c1["l3"]["rzz"] for c1, c2 in SO.values()]; rd_ = [c2["second_direct"]["rzz"] / c1["direct"]["rzz"] for c1, c2 in SO.values()]
    ps_ = [(c2["second_l3"]["rzz"] / c2["n2"]) / (c1["l3"]["rzz"] / c1["n"]) for c1, c2 in SO.values()]
    dp_ = [c2["second_l3"]["depth"] / c1["l3"]["depth"] for c1, c2 in SO.values()]
    nc_ = [c2["second_l3"]["nc_total"] / c1["l3"]["nc_total"] for c1, c2 in SO.values()]
    A(f"*Summary ({len(SO)} cases, N = {min(k[0] for k in SO)}–{max(k[0] for k in SO)}):* a second-order step costs {f(min(ps_), '.2f')}–{f(max(ps_), '.2f')}× the RZZ of a first-order step, "
      f"but the certified step count is {f(min(c2['n1'] / c2['n2'] for _, c2 in SO.values()), '.0f')}–{f(max(c2['n1'] / c2['n2'] for _, c2 in SO.values()), '.0f')}× smaller, so the certified circuit uses "
      f"{f(min(rr_), '.3f')}–{f(max(rr_), '.3f')}× the first-order RZZ (level 3; median {f(float(np.median(rr_)), '.3f')}) and "
      f"{f(min(rd_), '.3f')}–{f(max(rd_), '.3f')}× with direct compilation, {f(min(nc_), '.3f')}–{f(max(nc_), '.3f')}× the non-Clifford rotations and "
      f"{f(min(dp_), '.3f')}–{f(max(dp_), '.3f')}× the depth. The second-order formula is *not* cheaper at fixed step count; it wins because the certificate's Trotter budget is so tight "
      "(first-order error ∝ T²/n, second-order ∝ T³/n²) that the certified filter lives deep in the asymptotic regime. The advantage grows with N in this range "
      "(longer T), and nothing was computed beyond N = 12. Against each order's own oracle the certified count is "
      f"{f(min(c2['n1'] / c2['empirical1']['n'] for _, c2 in SO.values()), '.0f')}–{f(max(c2['n1'] / c2['empirical1']['n'] for _, c2 in SO.values()), '.0f')}× too large at first order and "
      f"{f(min(c2['n2'] / c2['empirical2']['n'] for _, c2 in SO.values()), '.0f')}–{f(max(c2['n2'] / c2['empirical2']['n'] for _, c2 in SO.values()), '.0f')}× at second order.")
    A("")
    A("*Checks:* in every case the second-order circuit was simulated against ED (same trial state): measured 1−F and the state error are below their certified bounds "
      f"(bound held: {all(c2['bound_holds'] for _, c2 in SO.values())}); for N ≤ 8 the qiskit `SuzukiTrotter(order=2)` circuit matches the independent simulation to ≤ "
      f"{max(abs(c2.get('circ_vs_sim_dev', 0)) for _, c2 in SO.values()):.0e} in 1−|⟨ψ|φ⟩|²; first-order RZZ from the full transpile equals the per-step × n estimate of Section 4. "
      "The pulse times and phases are the first-order designs (not re-optimised for second order); the second-order step counts use integer k_i = ⌈s·t_i⌉ with the smallest s meeting ε.")
    A("")

# ------------------------------------------------------------------ figures
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
COL = {0.0: "#1f4e79", 0.2411: "#2e8b57", 0.4: "#c0392b"}


def series(src, J2, which, fn, Nmin=4, Nmax=34):
    xs, ys = [], []
    for N in range(Nmin, Nmax + 1, 2):
        r = one(src, N, J2, EPS_MAIN)
        d = D(r, which)
        if ok(d) and r["inputs"].get("dmrg_converged", True):          # unconverged excited-state DMRG: not plotted
            xs.append(N); ys.append(fn(r, d))
    return xs, ys


panels = [("n (Trotter steps)", lambda r, d: d["n"]), ("T (physical, 1/J1)", lambda r, d: d["T_phys"]),
          ("total RZZ", lambda r, d: tot(d)["rzz"]), ("non-Clifford rotations (direct)", lambda r, d: drc(d)["nc_total"])]
fig, axes = plt.subplots(1, 4, figsize=(14, 3.3))
for ax, (title, fn) in zip(axes, panels):
    for J2 in J2_ORDER:
        xs, ys = series("dmrg", J2, "guaranteed", fn)
        if xs:
            ax.semilogy(xs, ys, "-", marker="o", color=COL[J2], ms=4, label=f"J2={J2:g}")
    ax.set_xlabel("N"); ax.set_title(title, fontsize=9)
axes[0].legend(fontsize=6, frameon=False)
fig.suptitle("Guaranteed designs, ε = 10⁻² (DMRG spectrum, certified bandwidth; N > 14 counted, not simulated; converged-DMRG points only)", fontsize=9)
fig.tight_layout(); fig.savefig(FIG / "fig_scaling_guaranteed.png", dpi=170); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
for J2 in J2_ORDER:
    for which, ls in (("guaranteed", "-"), ("empirical", ":")):
        xs, ys = series("dmrg", J2, which, lambda r, d: tot(d)["rzz"])
        if xs:
            ax[0].semilogy(xs, ys, ls, marker="o", ms=4, color=COL[J2], label=f"J2={J2:g} {which}")
        xs, ys = series("dmrg", J2, which, lambda r, d: d["n"])
        if xs:
            ax[1].semilogy(xs, ys, ls, marker="o", ms=4, color=COL[J2])
ax[0].set_title("total RZZ (N ≤ 14)", fontsize=9); ax[1].set_title("Trotter steps n (N ≤ 14)", fontsize=9)
for a in ax:
    a.set_xlabel("N")
ax[0].legend(fontsize=6, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig_guaranteed_vs_empirical.png", dpi=170); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.6, 3.3))
for (N, J2), mk in zip(((6, 0.0), (8, 0.0), (6, 0.4)), "osd"):
    for which, ls in (("guaranteed", "-"), ("empirical", ":")):
        xs, ys = [], []
        for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
            d = D(one("dmrg", N, J2, eps), which)
            if ok(d):
                xs.append(eps); ys.append(tot(d)["rzz"])
        if xs:
            ax.loglog(xs, ys, ls, marker=mk, ms=4, label=f"N={N}, J2={J2:g}, {which}")
ax.invert_xaxis(); ax.set_xlabel("target ε"); ax.set_ylabel("total RZZ"); ax.legend(fontsize=6, frameon=False)
fig.tight_layout(); fig.savefig(FIG / "fig_eps_scan.png", dpi=170); plt.close(fig)

if SO:
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    for J2 in J2_ORDER:
        pts = sorted((N, c1, c2) for (N, j), (c1, c2) in SO.items() if abs(j - J2) < 1e-9)
        xs = [p[0] for p in pts]
        ax[0].semilogy(xs, [p[1]["l3"]["rzz"] for p in pts], "-", marker="o", ms=4, color=COL[J2], label=f"J2={J2:g} first order")
        ax[0].semilogy(xs, [p[2]["second_l3"]["rzz"] for p in pts], "--", marker="s", ms=4, color=COL[J2], label=f"J2={J2:g} second order")
        ax[1].semilogy(xs, [p[1]["n"] for p in pts], "-", marker="o", ms=4, color=COL[J2])
        ax[1].semilogy(xs, [p[2]["n2"] for p in pts], "--", marker="s", ms=4, color=COL[J2])
    ax[0].set_title("certified filter, total RZZ (level 3)", fontsize=9); ax[1].set_title("certified Trotter steps n", fontsize=9)
    for a in ax:
        a.set_xlabel("N")
    ax[0].legend(fontsize=5.5, frameon=False)
    fig.tight_layout(); fig.savefig(FIG / "fig_first_vs_second_order.png", dpi=170); plt.close(fig)

A("## Figures")
A("")
if SO:
    A("![First vs second order](figures/fig_first_vs_second_order.png)")
    A("")
A("![Guaranteed designs vs N](figures/fig_scaling_guaranteed.png)")
A("")
A("![Guaranteed vs empirical](figures/fig_guaranteed_vs_empirical.png)")
A("")
A("![ε scan](figures/fig_eps_scan.png)")
A("")
A("## Reproducing")
A("")
A("```")
A("export OMP_NUM_THREADS=1")
A("for i in 0 1 2 3; do python scripts/run_results.py --worker $i --of 4 & done; wait   # resumable")
A("python scripts/add_direct_counts.py   # direct-compilation counts (run_results.py now writes them itself)")
A("python scripts/audit_results.py [--alpha]   # independent re-derivation of every stored bound and count; must print ALL CHECKS PASSED")
A("python scripts/build_results.py")
A("```")
A("")
A("Environment: Python 3.12, qiskit, qiskit-aer, quimb, mps-to-circuit, scipy, numpy, matplotlib (see `results/data/*.json` for the "
  "exact settings of each point: order, bandwidth, trial layers, seed).")

# ------------------------------------------------------------------ summary (inserted after section 0)
S = ["## Summary", ""]
g_small = [(r, D(r, "guaranteed"), D(r, "empirical")) for r in sel("dmrg", eps=EPS_MAIN) if r["N"] <= 14]
g_small = [(r, g, e) for r, g, e in g_small if ok(g)]
g_large_all = [(r, D(r, "guaranteed")) for r in sel("dmrg", eps=EPS_MAIN) if r["N"] > 14 and ok(D(r, "guaranteed"))]
g_large = [(r, g) for r, g in g_large_all if r["inputs"].get("dmrg_converged", True) and from_floor(r)]
n_fallback = sum(1 for r, _ in g_large_all if not from_floor(r))
if g_small:
    ns = [g["n"] for _, g, _ in g_small]
    cxs = [tot(g)["rzz"] for _, g, _ in g_small]
    S.append(f"* **N ≤ 14, guaranteed designs, ε = 10⁻²:** {f(min(ns), 'int')}–{f(max(ns), 'int')} Trotter steps, "
             f"{f(min(cxs), 'int')}–{f(max(cxs), 'int')} total RZZ, "
             f"T (phys) {f(min(g['T_phys'] for _, g, _ in g_small), '.1f')}–{f(max(g['T_phys'] for _, g, _ in g_small), '.1f')} (1/J1).")
    pairs = [(g, e) for _, g, e in g_small if ok(e)]
    if pairs:
        rn = [g["n"] / e["n"] for g, e in pairs]
        rc = [tot(g)["rzz"] / tot(e)["rzz"] for g, e in pairs]
        S.append(f"* **Guaranteed versus empirical (oracle):** the guaranteed design uses {f(min(rn), '.0f')}–{f(max(rn), '.0f')}× the steps and "
                 f"{f(min(rc), '.0f')}–{f(max(rc), '.0f')}× the RZZ of the smallest n that empirically meets ε (median {f(float(np.median(rc)), '.0f')}×).")
if g_large:
    cxs = [tot(g)["rzz"] for _, g in g_large]
    nunc = sum(1 for r, _ in g_large_all if not r["inputs"].get("dmrg_converged", True))
    S.append(f"* **N = {min(r['N'] for r, _ in g_large)}–{max(r['N'] for r, _ in g_large)} (DMRG/MPS only, counted not simulated), ε = 10⁻²:** {f(min(g['n'] for _, g in g_large), 'int')}–"
             f"{f(max(g['n'] for _, g in g_large), 'int')} steps, {min(cxs):.2e}–{max(cxs):.2e} total RZZ (level 3; direct compilation {min(drc(g)['rzz'] for _, g in g_large):.2e}–{max(drc(g)['rzz'] for _, g in g_large):.2e}), "
             f"{min(drc(g)['nc_total'] for _, g in g_large):.2e}–{max(drc(g)['nc_total'] for _, g in g_large):.2e} non-Clifford rotations (direct; level 3: {min(tot(g)['nc_total'] for _, g in g_large):.2e}–{max(tot(g)['nc_total'] for _, g in g_large):.2e}), "
             f"T (phys) {f(min(g['T_phys'] for _, g in g_large), '.1f')}–{f(max(g['T_phys'] for _, g in g_large), '.1f')} "
             f"(converged-DMRG points only; {nunc} further point(s) have an unconverged excited-state DMRG and are flagged in Section 7; {n_fallback} point(s) designed by the fallback FilterBuilder, marked † there, are excluded from these ranges and from the fits).")
for J2 in J2_ORDER:
    k = ("RZZ total", J2, "dmrg")
    if k in FITS:
        S.append(f"* **Scaling, J2 = {jl(J2)}:** guaranteed total RZZ ~ N^{FITS[k][0][0]:.1f} (R² {FITS[k][0][1]:.3f}), "
                 f"steps ~ N^{FITS[('n (steps)', J2, 'dmrg')][0][0]:.1f}, RZZ per step ~ N^{FITS[('RZZ per step', J2, 'dmrg')][0][0]:.1f} "
                 f"(fits over N = {FITS[k][1][0][0]}–{FITS[k][1][-1][0]}).")
ratios = [r["inputs"]["W"] / r["inputs"]["W_exact_dmrg"] for r in R]
S.append(f"* **Certified bandwidth:** the Weyl bound W = Σ_bonds λmax − E0 equals E_top − E0 (ratio {min(ratios):.3f}–{max(ratios):.3f}); W cancels in n and in T (phys), see Section 6.")
S.append(f"* **Bound violations among simulated points:** {len(viol)}.")
if SO:
    _rr = [c2["second_l3"]["rzz"] / c1["l3"]["rzz"] for c1, c2 in SO.values()]
    _ps = [(c2["second_l3"]["rzz"] / c2["n2"]) / (c1["l3"]["rzz"] / c1["n"]) for c1, c2 in SO.values()]
    _dp = [c2["second_l3"]["depth"] / c1["l3"]["depth"] for c1, c2 in SO.values()]
    S.append(f"* **First versus second order (Section 10, N = 6–12, ε = 10⁻², full transpiles of both circuits):** a second-order step costs {min(_ps):.2f}–{max(_ps):.2f}× a first-order step, "
             f"but the certified step count is {min(c2['n1'] / c2['n2'] for _, c2 in SO.values()):.0f}–{max(c2['n1'] / c2['n2'] for _, c2 in SO.values()):.0f}× smaller: the certified second-order circuit uses "
             f"{min(_rr):.3f}–{max(_rr):.3f}× the first-order RZZ (level 3) and {min(_dp):.3f}–{max(_dp):.3f}× the depth.")
S.append("")
i1 = next(i for i, l in enumerate(L) if l.startswith("## 1. Case inputs"))
L[i1:i1] = S

text = "\n".join(L)
(RES / "RESULTS.md").write_text(text)
print(f"wrote {RES / 'RESULTS.md'} ({len(text.split())} words), {len(list(TAB.glob('*.csv')))} csv tables, {len(list(FIG.glob('*.png')))} figures")
(RES / "style.css").write_text("body{max-width:none!important;margin:0 auto;padding:1.5em 2em;font:15px/1.45 Georgia,serif}"
                               "table{display:block;overflow-x:auto;border-collapse:collapse;font-size:.74em;margin:1em 0}"
                               "td,th{border-bottom:1px solid #ccc;padding:2px 7px;white-space:nowrap}img{max-width:100%}")
if shutil.which("pandoc"):
    subprocess.run(["pandoc", "RESULTS.md", "-o", "RESULTS.html", "--standalone", "--embed-resources", "--metadata", "title=Filter pipeline: resource results",
                    "--from=markdown+pipe_tables", "--css=style.css"], cwd=RES, check=False)
