---
title: "From DMRG to a Certified Ground State: A Workflow for Stetcu–Baroni Filtering of Matrix-Product Trial States, with Its Accuracy and Gate Cost"
author: "Draft — authors to be added"
date: "8 October 2026"
---

## Abstract

Density-matrix renormalisation group (DMRG) calculations give a matrix-product state and the low-lying energies of a spin chain classically, but a compact trial circuit built from them keeps some excited-state weight. The projection filter of Stetcu, Baroni and Carlson removes it on a quantum computer using one ancilla, controlled time evolution and post-selection. We give a complete workflow around this filter, from the DMRG front end to a verified state, and evaluate it on the open $J_1$–$J_2$ spin-½ chain for $N=4$–$16$ with real DMRG trial states compiled to circuits. The workflow contains a certificate: an a-priori bound on the trace distance $D$ to the exact ground state, $D\le\sqrt{\bar\ell}+\epsilon_T/\sqrt{p_g}$, a leakage term plus a Trotter term joined by the triangle inequality. We find that (i) DMRG supplies the inputs to $10^{-8}$ relative accuracy; (ii) the certificate is never violated ({{npts}} simulated points) and is within {{tot_lo}}–{{tot_hi}}$\times$ of the measured distance, because its leakage term is nearly tight, while its Trotter term is {{trot_lo}}–{{trot_hi}}$\times$ too large; (iii) a non-oracle *hybrid* rule (certified leakage plus a measured $n$-versus-$2n$ Trotter estimate) met every target ({{hyb_met}} of {{hyb_cases}} cases) with {{nbh_lo}}–{{nbh_hi}}$\times$ fewer Trotter steps than the certificate; (iv) the certificate fails predictably if the classical inputs are optimistic (a gap over-estimated by $25\%$ already breaks it); and (v) in gates, the hybrid filter costs $3\times10^{2}$–$4\times10^{4}$ CX including the trial circuit, growing as roughly $N^{2.5}$–$N^{3.1}$, which is {{dir_hyb_lo}}–{{dir_hyb_hi}}$\times$ the cost of preparing the state directly from the MPS at the same fidelity. The workflow is therefore a certified, reproducible way to use the filter, and a calibrated account of its cost, but not a resource advantage on these chains.

## 1. Introduction

A standard hybrid route to a ground state is to run DMRG, compile the resulting matrix-product state (MPS) into a shallow circuit, and repair the circuit's residual error with a quantum projection step. The Stetcu–Baroni–Carlson (SBC) projection [1], related to the rodeo algorithm [2], multiplies the trial state's energy components by a filter $F(E)=\prod_i\cos(Et_i+\phi_i)$ realised with pulses of controlled time evolution on a single ancilla. It sits among a family of ground-state projectors that use polynomial or Fourier filters of the Hamiltonian: phase-estimation-based preparation [3], quantum-eigenvalue-transformation filters [4,5] and, more generally, quantum singular value transformation. The SBC filter trades the optimal $\log(1/\varepsilon)$ scaling of those methods for a measurement-based, single-ancilla circuit whose filter can be designed and certified classically.

Someone deciding whether to use it needs to know **what the full workflow is, how much of its guarantee is rigorous and how much is assumed, and what it costs in gates**. This paper answers those questions for one concrete, fully reproducible workflow:

1. Section 2 specifies the workflow stage by stage, with the classical inputs it consumes and the certificate it produces.
2. Section 3 describes the experiments. Section 4 reports accuracy: the quality of the DMRG inputs, how tight the certificate is, how the step count should be chosen, how the trial state matters, what happens when inputs are wrong, and how sensitive the result is to design choices.
3. Section 5 reports cost in CX gates, including the trial circuit, how it scales, and how it compares with preparing the state directly from the MPS.
4. Section 6 states what the workflow does and does not establish.

We do not claim a resource advantage; Section 5 shows there is none at these sizes.

## 2. The workflow

![Figure 1. The workflow.](figs/fig_workflow.png){width=100%}

**Stage 1 — classical front end (DMRG).** From a Hamiltonian $H=\sum_gH_g$ (here the bonds of the chain) run DMRG for the ground-state energy $E_0$ and a high-bond-dimension reference MPS, DMRG on $-H$ for the top energy $E_{\mathrm{top}}$, and DMRG on $H+\lambda|\psi_0\rangle\langle\psi_0|$ with $\lambda=1.1(E_{\mathrm{top}}-E_0)$ for the first excitation $E_1$ (penalty method). Set $W=E_{\mathrm{top}}-E_0$, $H_s=(H-E_0)/W$, so that $\operatorname{spec}H_s\subset\{0\}\cup[\Delta,1]$ with $\Delta=(E_1-E_0)/W$.

**Stage 2 — trial state.** Run DMRG at a small maximum bond dimension $\chi$ and compile the MPS to a circuit (here exactly, by sequential isometries with `mps-to-circuit`; an approximate brickwork compilation is also tested). Estimate the ground-state weight $\gamma=|\langle E_0|\psi\rangle|^2$ as the overlap with the reference MPS. If $1-\gamma\le\varepsilon$ for the target $\varepsilon$, stop: the trial state already suffices.

**Stage 3 — filter.** One pulse $(t,\phi)$ applies $\mathsf H_{\mathrm{a}}\,e^{-itH_s\otimes Z_{\mathrm{a}}}\,\mathrm{Rz}_{\mathrm{a}}(2\phi)\,\mathsf H_{\mathrm{a}}$; its ancilla-$|0\rangle$ block is $\cos(H_st+\phi)$. Keeping outcome $0$ after $m$ pulses applies $F(E)=\prod_{i=1}^m\cos(Et_i+\phi_i)$, with total time $T=\sum_it_i$ and $F(0)=\prod_i\cos\phi_i$. Two numbers characterise the filter: $\eta=\sup_{E\in[\Delta,1]}|F(E)|/|F(0)|$ (suppression of excited states relative to the ground state) and $p_g=\gamma F(0)^2$ (a lower bound on the success probability). We design $(t_i,\phi_i)$ with `floor.py`, which maximises $F(0)^2$ subject to a *certified* bound on $\eta$ (interval bisection using $|F'|\le\sum|t_i|$ and $|F''|\le(\sum|t_i|)^2$), for a leakage budget $\varepsilon_\ell=\varepsilon/4$.

**Stage 4 — Trotterisation and the choice of $n$.** Each controlled evolution is $k_i$ first-order Lie–Trotter steps over the bonds, $n=\sum_ik_i$, $t_i=k_i\,dt$. The Trotter error is governed by the commutator constant $\alpha=\sum_g\big\|[H_g,\sum_{g'>g}H_{g'}]\big\|$ (the larger of the two product orders; the ancilla does not change it). Three rules choose $n$:

* **Certified.** The smallest $n$ with $\bar D(n)\le\sqrt\varepsilon$, where the *certificate* is
$$\boxed{\;D\le\bar D:=\sqrt{\bar\ell}+\frac{\alpha T^2}{2n\sqrt{\gamma F(0)^2}},\qquad\bar\ell=\frac{(1-\gamma)\eta^2}{\gamma+(1-\gamma)\eta^2}\;}$$
$D$ is the trace distance between the implemented and the exact ground state. The first term is **leakage** (what an exact filter would leave) and depends only on $\gamma$ and $\eta$; the second is the **Trotter** term. Equivalently, infidelity $\le\bar D^2$ and energy error $\langle H\rangle-E_0\le W\bar D^2$. The proof (Appendix A) is the triangle inequality for the trace distance plus a leakage lemma and a first-order Trotter lemma.
* **Hybrid.** Keep the rigorous leakage term but replace the worst-case Trotter term by a measured one. For each candidate $n$ on a geometric grid, run the same filter at $n$ and $2n$ steps; because the first-order error scales as $1/n$, the Trotter distance is $\approx2D(\tilde\varphi_n,\tilde\varphi_{2n})$. Accept the smallest $n$ such that $\sqrt{\bar\ell}+2D(\tilde\varphi_n,\tilde\varphi_{2n})\le\sqrt\varepsilon$ at that grid point and the next. This uses no ground state and is *not* a theorem (it assumes the $1/n$ regime), but it is what a practitioner can actually evaluate.
* **Measured (oracle).** The smallest $n$ whose measured $D$ meets the target. It needs the exact ground state, so it is a reference, not a rule.

**Stage 5 — run and check.** Execute the filter; check the post-selection rate, which is at least $p_g$, and, where the energy is measurable, $\langle H\rangle-E_0\le W\bar D^2$. If the device is noisy and prepares a mixed state $\tilde\rho$, the triangle inequality adds a noise term: $D(\tilde\rho,g)\le D(\tilde\rho,\tilde\varphi)+\bar D$.

**What is rigorous and what is assumed.**

| Quantity | Source | Status |
|---|---|---|
| $E_0$ | DMRG | variational upper bound; accurate to $\sim10^{-9}$ here |
| $E_{\mathrm{top}}$, $W$ | DMRG (or $W=\tfrac14\sum_gJ_g-E_0$ from the bond norms, valid for $J_g>0$) | the bond bound is rigorous up to the error in $E_0$; DMRG is accurate but uncertified |
| $E_1$, $\Delta$ | DMRG penalty method | an *upper* bound on $E_1$: the wrong side for a gap lower bound |
| $\gamma$ | overlap with a high-$\chi$ reference MPS | an overlap with an approximation of $|E_0\rangle$ |
| $\eta$, $F(0)$ | interval certificate | rigorous given $\Delta$ |
| $\alpha$ | exact operator norms of commutators | rigorous |
| leakage and Trotter terms | Lemmas 1 and 2 | rigorous *given* the inputs above |

The certificate is therefore as trustworthy as $\Delta$ and $\gamma$; Section 4.5 quantifies what happens when they are wrong.

## 3. Experiments

* **System.** Open spin-½ $J_1$–$J_2$ chain, $J_1=1$, $J_2\in\{0,0.4\}$, $N\in\{4,6,\dots,16\}$. Bonds are ordered: nearest-neighbour bonds with even then odd left site, then next-nearest-neighbour bonds likewise.
* **Front end.** Two-site DMRG (`quimb`) with the repository's `get_energies.py`; reference MPS up to bond dimension 64. The workflow sees only DMRG quantities. Exact diagonalisation is used afterwards to validate (true ground state, gap, overlap) and to measure distances.
* **Trial state.** DMRG at maximum bond dimension $\chi=2$ (main experiments), compiled exactly. $\chi\in\{1,2,4,8\}$ and $L$-layer approximate circuits are compared at $N=8$.
* **Targets.** $\varepsilon=10^{-1},10^{-2}$ ($N\le16$) and $10^{-3}$ ($N\le8$), i.e. $D\le0.32,0.1,0.032$.
* **Simulation.** Exact statevector simulation of the Trotterised filter, the ideal filter, and the ground state. A run is simulated only if $n\times\text{bonds}\times2\times2^N\le1.5\times10^{10}$; the $N=16$, $\varepsilon=10^{-2}$ certified circuits are not simulated, and their step counts and CX are reported from the certificate.
* **Cost.** One Trotter step of one pulse (for each bond, the $XX,YY,ZZ$ rotations tensored with $Z$ on the ancilla) is transpiled once with qiskit (level 3) to $\{\mathrm{CX},\mathrm{Rz},\mathrm{H},\mathrm{S}\}$; total CX is trial CX plus $n\times$ CX per step. All-to-all connectivity unless stated. The trial circuit is the exact `mps-to-circuit` compilation, transpiled the same way.
* **Direct-preparation baseline.** Preparing the state directly from the MPS: `mps-to-circuit` sequential (exact) circuits for the reference MPS truncated to $\chi\in\{1,2,4,8,16\}$, $L$-layer approximate circuits ($L\le8$), and generic state preparation of the exact ground state ($N\le12$). Each is measured against the exact ground state; the baseline for a target $\varepsilon$ is the cheapest with infidelity $\le\varepsilon$.

## 4. Accuracy

### 4.1 The classical front end supplies the inputs

{{inputs}}

*Table 1.* DMRG-derived inputs against exact diagonalisation. $\chi_{\mathrm{ref}}$ is the bond dimension the reference MPS used.

DMRG reproduces $E_0$, $E_1$ and the gap to better than {{max_E}} (absolute) and {{max_gap}} (relative), and the estimated $\gamma$ differs from the exact overlap by at most {{max_gam}}. At these sizes the certificate computed from DMRG inputs is therefore indistinguishable from the one computed from exact inputs. This is an empirical statement about well-behaved chains, not a guarantee: DMRG does not certify these numbers (Section 2) and the penalty $E_1$ is variationally on the unsafe side. Section 4.5 shows what an error of that sign does.

### 4.2 The certificate is accurate, mostly because leakage is tight

{{main}}

*Table 2.* $\varepsilon=10^{-2}$, $\chi=2$. $^\dagger$$1-\gamma\le\varepsilon$: no filter needed, and the entry is the trial state's own distance. $D$ is at $n_{\mathrm{bound}}$. The last two columns compare the certificate with the measurement, term by term.

* **Never violated.** Across {{npts}} simulated points (every case at $n_{\mathrm{bound}}$ plus the $n$-sweeps) the total certificate, its leakage term, its Trotter term, the triangle inequality, the energy bound and $P_{\mathrm{succ}}\ge\gamma F(0)^2$ each held {{npts}} out of {{npts}} times. This checks the derivation and the code; the proofs are in Appendix A.
* **Close in total.** At $n_{\mathrm{bound}}$ the certificate exceeds the measured $D$ by {{tot_lo}}–{{tot_hi}}$\times$ (median {{tot_med}}$\times$).
* **Leakage is nearly tight; Trotter is not.** The leakage term is within {{leak_lo}}–{{leak_hi}}$\times$; the Trotter term is {{trot_lo}}–{{trot_hi}}$\times$ too large and the discrepancy grows with $N$ (from $9\times$ at $N=6$ to $52$–$77\times$ at $N=14$), because $\alpha$ is a worst case over all states. The measured $D$ at $n_{\mathrm{bound}}$ is dominated by leakage ({{leak_meas_lo}}–{{leak_meas_hi}}, against a Trotter part of at most {{trot_meas_max}}).
* **Scaling in $n$.** The Trotter certificate and the measured Trotter distance both fall as $1/n$ (fitted log–log slopes {{slope_b_lo}} to {{slope_b_hi}} and {{slope_m_lo}} to {{slope_m_hi}}), so the certificate has the right scaling in $n$ and is off by a roughly constant factor (Figure 2). The measured total saturates at the leakage floor.
* **Other observables.** The energy bound $W\bar D^2$ holds but is {{en_lo}}–{{en_hi}}$\times$ above the measured energy error. The success rate matches the lower bound $\gamma F(0)^2$ to {{psucc_max}}% or better, so $1/(\gamma F(0)^2)$ is an accurate estimate of the repetition overhead.

![Figure 2. Dependence on $n$ at fixed design ($\varepsilon=10^{-2}$, $J_2=0$, $N=4,6,8,10$ as solid, dashed, dotted, dash-dotted). (a) Certificate (blue) and measured distance (red). (b) Trotter part only: certificate (blue) and measurement (green).](figs/fig_sweep.png){width=100%}

### 4.3 Choosing $n$: certified, hybrid, measured

{{hybrid}}

*Table 3.* The hybrid rule against the certificate and the oracle. $D$ is measured at $n_{\mathrm{hyb}}$ against the exact ground state, which the rule itself never uses.

* The hybrid rule met the target in {{hyb_met}} of {{hyb_cases}} cases, with a margin $\sqrt\varepsilon/D$ of {{hmarg_lo}}–{{hmarg_hi}}$\times$.
* It uses {{nbh_lo}}–{{nbh_hi}}$\times$ (median {{nbh_med}}$\times$) fewer steps than the certificate and only {{nhm_lo}}–{{nhm_hi}}$\times$ (median {{nhm_med}}$\times$) more than the oracle. The gap to the oracle is mostly the coarse geometric grid and the requirement of confirmation at the next grid point.
* The n-versus-$2n$ estimate of the Trotter distance matched the measured one to within [{{rich_lo}}, {{rich_hi}}] over {{rich_n}} sweep points with $n\ge n_{\mathrm{meas}}$.
* The hybrid rule is not a theorem. It assumes the $1/n$ regime, which held in all runs here, and inherits the assumptions on $\gamma$ and $\Delta$ through its leakage term.

Step counts scale with the target roughly as the first-order theory predicts for the certificate ($n\propto\varepsilon^{-1/2}$ at fixed filter; observed ratios $n(10^{-3})/n(10^{-2})$ of $4$–$7$ include the tightening of the filter itself):

{{eps}}

*Table 4.* Certified step counts at three targets, and hybrid and measured counts at $10^{-2}$ and $10^{-3}$.

### 4.4 The choice of trial state

At $N=8$, $\varepsilon=10^{-2}$:

{{trials}}

*Table 5.* Trial state against total CX (trial plus filter). "Approx." is the $L$-layer `mps-to-circuit` brickwork compilation of the reference MPS.

A better trial state needs a much shorter filter: going from $\chi=1$ ($\gamma=0.13$ at $J_2=0$) to $\chi=2$ cuts the certified step count by $7.7\times$ ($31\times$ at $J_2=0.4$, where $\chi=1$ gives $\gamma=0.04$). For $\chi=4$ the DMRG state already meets the target and no filter is needed. A trial state produced by approximate compilation behaves like a DMRG state of the same $\gamma$: the filter only sees $\gamma$.

### 4.5 When the classical inputs are wrong

The certificate assumes $\Delta$ and $\gamma$ are known. We design and certify with deliberately wrong values and measure the true distance ($N=8$ and $12$, $J_2=0$, $\varepsilon=10^{-2}$):

{{sens}}

*Table 6.* "Claimed" is the certificate computed from the (wrong) inputs; "with true inputs" recomputes it with the true gap and overlap; "measured" is the actual distance.

* **Under-estimating the gap is safe** (the filter is designed for a wider window and costs more steps).
* **Over-estimating the gap is not.** At $+10\%$ the claimed bound still holds, but the guaranteed distance with true inputs is already $30$–$45\%$ above the claim; at $+25\%$ the claim is violated and the target missed at both sizes.
* **Over-estimating $\gamma$** behaves the same way, at a threshold that depends on how much of the excited weight is hidden.
* This is why the penalty-method $E_1$, which is an *upper* bound, must be treated with care: in this study it was accurate to $10^{-8}$, but the workflow should either certify $\Delta$ independently or design with a deliberately *reduced* gap: a $20\%$ reduction was safe and cost $1.55\times$ more steps at both sizes.

### 4.6 Sensitivity to design choices

At $N=8$, $J_2=0$, $\varepsilon=10^{-2}$, certified $n$ against the leakage share $\varepsilon_\ell/\varepsilon$ and the number of pulses $m$ (the solver searched $x=T\Delta/\pi\in\{0.6,1,1.5,2\}$):

{{design}}

*Table 7.* Certified step count against design choices. Four pulses are infeasible for this trial state.

The certified $n$ varies by $3\times$ ({{des_lo}}–{{des_hi}}) over the shares tried, with a minimum near the $1/4$ used elsewhere; six and eight pulses are equivalent. The leakage share is the one design parameter worth tuning.

## 5. Cost

{{cost}}

*Table 8.* $\varepsilon=10^{-2}$, $\chi=2$, all-to-all connectivity. Filter CX is $n\times$ CX per step; totals add the trial circuit. The last column is the cheapest direct preparation of the state at infidelity $\le\varepsilon$ (Section 3).

![Figure 3. Total CX against $N$ for $\varepsilon=10^{-2}$.](figs/fig_cx.png){width=65%}

* **Per step**, CX grows about linearly with the number of bonds ($N^{{{fit0_step}}}$ at $J_2=0$); $J_2=0.4$ costs about $1.9\times$ more because it has about twice as many bonds. A line with the ancilla at the end costs {{line_lo}}–{{line_hi}}$\times$ more than all-to-all.
* **Scaling.** Over $N={{fit_Nlo}}$–${{fit_Nhi}}$, filter CX grows as $N^{{{fit0_filt_bound}}}$ (certified), $N^{{{fit0_filt_hyb}}}$ (hybrid) and $N^{{{fit0_filt_meas}}}$ (measured) at $J_2=0$, and as $N^{{{fit4_filt_bound}}}$, $N^{{{fit4_filt_hyb}}}$, $N^{{{fit4_filt_meas}}}$ at $J_2=0.4$. These are fits to at most six points. The steps themselves scale as $N^{{{fit0_n_bound}}}$ (certified) and $N^{{{fit0_n_meas}}}$ (measured), driven by $T^2\propto\Delta^{-2}$: the scaled gap falls from $0.28$ to $0.02$ over $N=4$–$16$.
* **Certified versus hybrid.** The certificate costs {{cxbh_lo}}–{{cxbh_hi}}$\times$ (median {{cxbh_med}}$\times$) the hybrid circuit; the hybrid circuit is {{cxhm_lo}}–{{cxhm_hi}}$\times$ the oracle.
* **Trial circuit.** For $\chi=2$ it costs $3(N-1)$ CX, negligible next to the filter.
* **Against direct preparation.** Directly compiling the MPS reaches the same infidelity with $11$–$250$ CX at every size except $N=16$, $J_2=0$ (a sequential $\chi=4$ circuit suffices at most sizes; there the required $\chi$ rises to $8$ and the cost to $1{,}089$ CX). The certified filter costs {{dir_bound_lo}}–{{dir_bound_hi}}$\times$ more, the hybrid filter {{dir_hyb_lo}}–{{dir_hyb_hi}}$\times$ more, and even the oracle filter {{dir_meas_lo}}–{{dir_meas_hi}}$\times$ more. The ratio for the hybrid filter drops where the direct method must double its bond dimension (from about $125\times$ at $N=14$ to $35\times$ at $N=16$, $J_2=0$), because the filter's cost depends on $1/\Delta^2$ and $\alpha$ and not on entanglement; but at these sizes it never crosses.

The direct-preparation baselines for each size:

{{baseline}}

*Table 9.* Direct preparation CX (infidelity to the exact ground state in parentheses), at $\chi=2,4,8$, $L=1,3$ approximate layers, and generic preparation ($N\le12$).

## 6. Discussion

**What the workflow gives.** A single procedure from DMRG output to a Trotter step count, with a certificate that held on every point tested, an explicit statement of which inputs it trusts, and a practical hybrid rule that is much cheaper than the certificate and met all targets. The certificate's leakage half is tight; its Trotter half is the price of a worst-case proof.

**What it costs.** At $\varepsilon=10^{-2}$ the hybrid circuit costs roughly $10^3$–$4\times10^4$ CX at $N=6$–$16$ against $30$–$1{,}100$ for direct MPS preparation. The filter is not competitive on these chains.

**When it could matter.** The filter's cost is set by the gap and the commutator constant, not by the state's entanglement; direct preparation's cost steps up each time the bond dimension must double. The data show the ratio narrowing at one such step ($N=16$, $J_2=0$) but not crossing. A genuine crossover would need systems where the required MPS bond dimension grows quickly with size (critical systems at larger $N$, or two dimensions) and where the filter can be run without paying for dense simulation. Neither is tested here.

**Limitations.**

* *Inputs.* The certificate needs $\gamma$ and $\Delta$. Here DMRG supplies them to $10^{-8}$ and exact diagonalisation validates them ($N\le16$). At larger $N$ there is no validation, and Section 4.5 shows the certificate breaks if the gap is optimistic. Certifying $\Delta$ at large $N$ is the main open point.
* *Models.* Two couplings of one chain.
* *Trotter order.* First order only; higher-order formulas would shrink the Trotter term but need a different constant.
* *Noise.* Not simulated; the triangle inequality gives an additive noise term.
* *Cost model.* CX only, all-to-all (a line costs about $2$–$3\times$ more), no T-count, no synthesis error; the CX per step is for a Heisenberg bond decomposition into three Pauli rotations and is not the cheapest possible.
* *Simulation.* Dense statevectors, so the large-$N$ reach of DMRG is not exercised, and the largest certified circuits are not simulated.
* *Hybrid rule.* An empirical rule that assumed the $1/n$ regime and had 25 test cases.

**Reproducibility.** All code is in `simple_paper/`: `dmrg_inputs.py` (DMRG front end and trial circuits), `workflow.py` (Hamiltonian, pulses, certificate), `run_experiments.py` (main grid), `run_extras.py` (trial states, wrong inputs, design choices, direct-preparation baseline), `resources.py` (CX counts), `pinsker_check.py`, `make_report.py`, `build_paper.py`; filter design is `j1j2_filter/floor.py`. They run in a Python 3.12 environment with `numpy`, `scipy`, `quimb`, `qiskit` and `mps-to-circuit`; set `OMP_NUM_THREADS=1` when running processes in parallel. The full set (all tables here) runs in one to two hours on four cores. Seeds are fixed, including the Lanczos start vector used for the exact ground state at $N\ge12$ (a $\chi=2$ cut passes through degenerate Schmidt values, so $\gamma$ is otherwise only reproducible to the third digit).

## References

1. I. Stetcu, A. Baroni, J. Carlson, "Projection algorithm for state preparation on quantum computers," *Phys. Rev. C* **105**, 064308 (2022).
2. K. Choi, D. Lee, J. Bonitati, Z. Qian, J. Watkins, "Rodeo algorithm for quantum computing," *Phys. Rev. Lett.* **127**, 040505 (2021).
3. D. Poulin, P. Wocjan, "Preparing ground states of quantum many-body systems on a quantum computer," *Phys. Rev. Lett.* **102**, 130503 (2009).
4. L. Lin, Y. Tong, "Near-optimal ground state preparation," *Quantum* **4**, 372 (2020).
5. Y. Dong, L. Lin, Y. Tong, "Ground-state preparation and energy estimation on early fault-tolerant quantum computers via quantum eigenvalue transformation of unitary matrices," *PRX Quantum* **3**, 040305 (2022).
6. A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, "Theory of Trotter error with commutator scaling," *Phys. Rev. X* **11**, 011020 (2021).
7. C. Schön, E. Solano, F. Verstraete, J. I. Cirac, M. M. Wolf, "Sequential generation of entangled multiqubit states," *Phys. Rev. Lett.* **95**, 110503 (2005).
8. S.-J. Ran, "Encoding of matrix product states into quantum circuits of one- and two-qubit gates," *Phys. Rev. A* **101**, 032310 (2020).
9. S. R. White, "Density matrix formulation for quantum renormalization groups," *Phys. Rev. Lett.* **69**, 2863 (1992); U. Schollwöck, "The density-matrix renormalization group in the age of matrix product states," *Ann. Phys.* **326**, 96 (2011).
10. M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* (Cambridge University Press, 2000), Ch. 9.
11. T. M. Cover, J. A. Thomas, *Elements of Information Theory* (Wiley, 2006), Lemma 11.6.1.

## Appendix A. Proof of the certificate

Write $D(\rho,\sigma)=\tfrac12\|\rho-\sigma\|_1$. For pure states $D=\sqrt{1-|\langle a|b\rangle|^2}=\sin\theta_{ab}$, with $\theta_{ab}$ the angle between the rays. $D$ obeys the triangle inequality. For any state $\rho$ and pure $g$, $1-\langle g|\rho|g\rangle\le D(\rho,g)\le\sqrt{1-\langle g|\rho|g\rangle}$; for pure $\rho$ the right side is an equality, so infidelity is $D^2$. Let $\varphi=F(H_s)\psi/\|F(H_s)\psi\|$ be the ideal filtered state and $\tilde\varphi$ the implemented one. The certificate follows from $D(\tilde\varphi,g)\le D(\tilde\varphi,\varphi)+D(\varphi,g)$ and two lemmas.

**Lemma 1 (leakage).** $D(\varphi,g)^2=\ell\le\bar\ell$.
*Proof.* The ground state sits at $E=0$, so $|\langle g|\varphi\rangle|^2=\gamma F(0)^2/(\gamma F(0)^2+y)$ with $y=\sum_{k\ge1}|c_k|^2F(E_k)^2$. All $E_k$, $k\ge1$, lie in $[\Delta,1]$, so $F(E_k)^2\le\eta^2F(0)^2$ and $y\le(1-\gamma)\eta^2F(0)^2$. Since $\ell=y/(\gamma F(0)^2+y)$ is increasing in $y$, the bound follows. $\square$

**Lemma 2 (Trotter).** $D(\tilde\varphi,\varphi)\le\min\{1,\epsilon_T/\sqrt{p_g}\}$ with $\epsilon_T=\sum_i\alpha t_i^2/(2k_i)=\alpha T^2/(2n)$.
*Proof.* Let $A_i=\langle0_{\mathrm{a}}|\mathcal W_i|0_{\mathrm{a}}\rangle$ and $\tilde A_i$ the Trotterised block; $a=A_m\cdots A_1\psi$, $b=\tilde A_m\cdots\tilde A_1\psi$. The first-order product-formula bound (below) gives $\|\mathcal W_i-\widetilde{\mathcal W}_i\|\le\alpha t_i^2/(2k_i)$; the blocks are compressions, so $\|A_i-\tilde A_i\|$ is no larger, and all blocks have norm $\le1$. Telescoping, $a-b=\sum_i(A_m\cdots A_{i+1})(A_i-\tilde A_i)(\tilde A_{i-1}\cdots\tilde A_1\psi)$, so $\|a-b\|\le\epsilon_T$. The rays of $a,b$ are $\varphi,\tilde\varphi$, so $D=\sin\angle(a,b)=\mathrm{dist}(a,\mathrm{span}\,b)/\|a\|\le\|a-b\|/\|a\|$, and $\|a\|^2=P_{\mathrm{succ}}\ge p_g$. $\square$

**Corollaries.** Infidelity $=D^2\le\bar D^2$. Energy: with $\tilde c_k$ the components of $\tilde\varphi$, $\langle H\rangle-E_0=\sum_{k\ge1}|\tilde c_k|^2(E_k-E_0)\le W(1-|\tilde c_0|^2)=WD^2$. Step count: for $\bar\ell\le\varepsilon/4$, $\bar D\le\sqrt\varepsilon$ holds if $n\ge\alpha T^2/\sqrt{\gamma F(0)^2\varepsilon}$; with $T=x\pi/\Delta$ the certified cost scales as $\alpha/\Delta^2$.

**First-order Trotter bound.** For Hermitian $A,B$, $U(s)=e^{-is(A+B)}$, $V(s)=e^{-isA}e^{-isB}$: $\frac{d}{ds}[U(\delta-s)V(s)]=iU(\delta-s)(B-e^{-isA}Be^{isA})V(s)$, so $U(\delta)-V(\delta)=-i\int_0^\delta U(\delta-s)(B-e^{-isA}Be^{isA})V(s)\,ds$. Since $\|B-e^{-isA}Be^{isA}\|\le s\|[A,B]\|$ and $U,V$ are unitary, $\|U(\delta)-V(\delta)\|\le\tfrac{\delta^2}2\|[A,B]\|$. For many terms apply this with $A=H_1$, $B=\sum_{g>1}H_g$ and recurse, giving $\frac{\delta^2}{2}\alpha$ per step; with $k$ steps of size $t/k$ the error is $\alpha t^2/(2k)$ [6]. The reversed product order gives the reversed sum, hence the maximum over orders. For $H_s\otimes Z$: $[H_gZ,H_{g'}Z]=[H_g,H_{g'}]\otimes\mathbb 1$, so $\alpha$ is unchanged; the constant shift commutes with everything and is applied exactly.

## Appendix B. Where Pinsker's inequality applies

For the target $|g\rangle\langle g|$ the relative entropy $S(\rho\|\,|g\rangle\langle g|)$ is infinite unless $\rho=|g\rangle\langle g|$, so $D\le\sqrt{S/2}$ gives nothing here, which is why the certificate uses the triangle inequality and fidelity instead. Pinsker does apply to the energy-resolved distributions $p_k=|c_k|^2$ (trial) and $q_k=p_kF(E_k)^2/P_{\mathrm{succ}}$ (filtered): $\mathrm{KL}(q\|p)=\sum_kq_k\ln(F(E_k)^2/P_{\mathrm{succ}})\le\ln(1/P_{\mathrm{succ}})$, and $\mathrm{TV}(p,q)\ge q_0-p_0=1-\ell-\gamma$, so with $\mathrm{TV}\le\sqrt{\mathrm{KL}/2}$,
$$P_{\mathrm{succ}}\le\exp\!\big(-2(1-\ell-\gamma)^2\big).$$
This is a necessary condition on the post-selection cost, not an error bound. We checked it on the {{pn}} $\varepsilon=10^{-2}$ cases with $N\le8$ that need a filter: every inequality held, but the cap ({{pk_lo}}–{{pk_hi}}) is far above the actual $P_{\mathrm{succ}}$ ({{pP_lo}}–{{pP_hi}}), so it is uninformative for these trial states.
