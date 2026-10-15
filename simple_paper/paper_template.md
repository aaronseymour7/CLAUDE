---
title: "From DMRG to a Certified Ground State: A Workflow for Stetcu–Baroni Filtering of Matrix-Product Trial States, with Its Accuracy and Gate Cost"
author: "Draft — authors to be added"
date: "15 October 2026"
---


## Abstract

Density-matrix renormalisation group (DMRG) gives a matrix-product state and the low-lying energies of a spin chain classically, but a compact trial circuit built from them keeps some excited-state weight. The projection filter of Stetcu, Baroni and Carlson removes it on a quantum computer using one ancilla, controlled time evolution and post-selection. We give a complete workflow around this filter, from the DMRG front end to a verified state, and evaluate it on the open $J_1$–$J_2$ spin-½ chain: $N=4$–$16$ with exact statevector simulation and real DMRG trial states compiled to circuits, and $N$ up to {{large_Nmax}} in a DMRG/MPS-only regime where circuits are counted. The front end is accurate to {{ed_max}} against exact diagonalisation ($N\le20$) and reproducible to {{cons_max}} between independent runs ($N\le32$); the excited state is the lowest state of any symmetry, as required across the Berezinskii–Kosterlitz–Thouless point of the chain, and the bandwidth is certified by Weyl's inequality without a DMRG value for the top of the spectrum. The paper is organised around three questions. **(1) Can the filter be certified without knowing the exact ground state?** Yes, conditionally: an a-priori bound on the trace distance to the ground state, $D\le\sqrt{\bar\ell}+\epsilon_T/\sqrt{p_g}$, needs only DMRG-derived inputs and held on all {{npts}} simulated points. It is a *conditional certificate*: it assumes the DMRG shift, gap and ground-state weight are right, which DMRG cannot certify, and we state exactly how it changes when they are off by given amounts. It fails if the weight $\gamma$ or the gap $\Delta$ is over-estimated (a gap $25\%$ too large already breaks it). **(2) How conservative is the certificate?** Its leakage term is reasonably tight ({{leak_lo}}–{{leak_hi}}$\times$); its worst-case Trotter term is extremely conservative ({{trot_lo}}–{{trot_hi}}$\times$, growing with $N$), so the certified step count is {{nbm_lo}}–{{nbm_hi}}$\times$ what is empirically sufficient. **(3) Can the practical cost be recovered without the exact ground state?** A hybrid rule, certified leakage plus a measured $n$-versus-$2n$ Trotter estimate, met every target ({{hyb_met}} of {{hyb_cases}}) with {{nbh_lo}}–{{nbh_hi}}$\times$ fewer steps than the certificate, but it is empirical, not rigorous. In gates, the hybrid filter costs {{hyb_cx_lo}}–{{hyb_cx_hi}} CX including the trial circuit for $N\le16$, growing as roughly $N^{{{hybfit_lo}}}$–$N^{{{hybfit_hi}}}$, which is {{dir_hyb_lo}}–{{dir_hyb_hi}}$\times$ the cost of preparing the state directly from the MPS at the same fidelity; the certified filter for $N=16$–${{large_Nmax}}$ needs {{lg_cx_lo}}–{{lg_cx_hi}} CX and {{lg_rz_lo}}–{{lg_rz_hi}} non-Clifford rotations. The workflow is a conditionally certified, reproducible way to use the filter and a calibrated account of its cost, not a resource advantage on these chains.

## 1. Introduction

A standard hybrid route to a ground state is to run DMRG, compile the resulting matrix-product state (MPS) into a shallow circuit, and repair the circuit's residual error with a quantum projection step. The Stetcu–Baroni–Carlson (SBC) projection [1], related to the rodeo algorithm [2], multiplies the trial state's energy components by a filter $F(E)=\prod_i\cos(Et_i+\phi_i)$ realised with pulses of controlled time evolution on a single ancilla. It sits among a family of ground-state projectors that use polynomial or Fourier filters of the Hamiltonian: phase-estimation-based preparation [3], quantum-eigenvalue-transformation filters [4,5] and, more generally, quantum singular value transformation. The SBC filter trades the optimal $\log(1/\varepsilon)$ scaling of those methods for a measurement-based, single-ancilla circuit whose filter can be designed and certified classically.

The practical difficulty is that the guarantee one would like, a bound on how far the prepared state is from the ground state, seems to need the ground state. This paper builds a workflow that avoids that, and asks three questions about it:

| | Question | Short answer | Where |
|---|---|---|---|
| **Q1** | Can SBC filtering be certified without knowing the exact ground state? | **Yes, conditionally** on conservative bounds for the ground-state weight $\gamma$ and the gap $\Delta$. A certificate built from DMRG-derived inputs was never violated, and it fails if $\gamma$ or $\Delta$ is over-estimated. | Section 4 |
| **Q2** | How conservative is that certificate? | **Leakage is reasonably tight; worst-case Trotter error is extremely conservative** (and more so as $N$ grows). | Section 5 |
| **Q3** | Can we recover the practical cost without the exact ground state? | **Largely, with the $n$-versus-$2n$ hybrid rule**, which met every target, but it is empirical rather than rigorous. | Section 6 |

Section 2 specifies the workflow and the certificate (including exactly which inputs are rigorous and which are assumed), Section 3 the experiments, and Section 7 the remaining practical questions: which trial state to use, how sensitive the result is to design choices, what the whole thing costs in CX gates, including against preparing the state directly from the MPS, and what it costs for chains up to $N={{large_Nmax}}$ where only DMRG and counted circuits are possible. Section 8 states the limitations. We do not claim a resource advantage; Section 7 shows there is none at these sizes.

## 2. The workflow

![Figure 1. The workflow.](figs/fig_workflow.png){width=100%}

**Stage 1 — classical front end (DMRG).** From a Hamiltonian $H=\sum_gH_g$ (here the bonds of the chain) run two-site DMRG for the ground-state energy $E_0$ and a high-bond-dimension reference MPS $\psi_0$, and again for the first excitation $E_1$ with the penalty $H+\lambda|\psi_0\rangle\langle\psi_0|$, $\lambda=1.1W$. The projector is applied as a rank-one term $\lambda|w\rangle\langle w|$ in every two-site local problem, where $w$ is $\psi_0$ contracted with the current environments; this is the same penalty as adding the projector as an MPO (whose bond dimension is $\chi^2$ and which ran out of memory beyond $N=14$ at the tolerances used here) at $O(\chi^3)$ cost. It returns the lowest state of *any* symmetry orthogonal to $\psi_0$, which matters because the first excitation of the $J_1$–$J_2$ chain is a triplet below the Berezinskii–Kosterlitz–Thouless point $J_2^c\approx0.2411$ [12,13] but can be a singlet in the dimerised phase. The bandwidth is certified without DMRG: by Weyl's inequality $\lambda_{\max}(\sum_bh_b)\le\sum_b\lambda_{\max}(h_b)$, and a bond $J\,\mathbf S\!\cdot\!\mathbf S$ has $\lambda_{\max}=J/4$ for $J\ge0$, so
$$W=\sum_b\lambda_{\max}(h_b)-E_0=\tfrac{N-1}{4}J_1+\tfrac{N-2}{4}J_2-E_0\ \ (J_1,J_2\ge0)$$
is an upper bound on $E_{\mathrm{top}}-E_0$, attained by the fully polarised state, so it equals $E_{\mathrm{top}}-E_0$. DMRG on $-H$ would be variational from the wrong side ($E_{\mathrm{top}}^{\rm DMRG}\le E_{\mathrm{top}}$) and is used only as a cross-check. Set $H_s=(H-E_0)/W$, so that $\operatorname{spec}H_s\subset\{0\}\cup[\Delta,1]$ with $\Delta=(E_1-E_0)/W$. Because $\alpha_{\rm scaled}=\alpha/W^2$ and $T_{\rm scaled}=x\pi W/g$ ($g=E_1-E_0$), the step count $n\propto\alpha T^2$ and the physical time $x\pi/g$ do not depend on $W$; a looser $W$ only widens the suppression window $[\Delta,1]$.

**Stage 2 — trial state.** Run DMRG at a small maximum bond dimension $\chi$ and compile the MPS to a circuit (here exactly, by sequential isometries with `mps-to-circuit`; an approximate brickwork compilation is also tested). Estimate the ground-state weight $\gamma=\|\Pi_0\psi\|^2$ (the weight in the whole ground space; a single nondegenerate vector in the cases studied) as the overlap with the reference MPS, by contraction when $N$ is too large for a statevector. If $1-\gamma\le\varepsilon$ for the target $\varepsilon$, stop: the trial state already suffices.

**Stage 3 — filter.** One pulse $(t,\phi)$ applies $\mathsf H_{\mathrm{a}}\,e^{-itH_s\otimes Z_{\mathrm{a}}}\,\mathrm{Rz}_{\mathrm{a}}(2\phi)\,\mathsf H_{\mathrm{a}}$; its ancilla-$|0\rangle$ block is $\cos(H_st+\phi)$. Keeping outcome $0$ after $m$ pulses applies $F(E)=\prod_{i=1}^m\cos(Et_i+\phi_i)$, with total time $T=\sum_it_i$ and $F(0)=\prod_i\cos\phi_i$. Two numbers characterise the filter: $\eta=\sup_{E\in[\Delta,1]}|F(E)|/|F(0)|$ (suppression of excited states relative to the ground state) and $p_g=\gamma F(0)^2$ (a lower bound on the success probability). We design $(t_i,\phi_i)$ with `floor.py`, which maximises $F(0)^2$ subject to a *certified* bound on $\eta$ (interval bisection using $|F'|\le\sum|t_i|$ and $|F''|\le(\sum|t_i|)^2$), for a leakage budget $\varepsilon_\ell=\varepsilon/4$.

**Stage 4 — Trotterisation and the choice of $n$.** Each controlled evolution is $k_i$ first-order Lie–Trotter steps over the bonds, $n=\sum_ik_i$, $t_i=k_i\,dt$. The Trotter error is governed by the commutator constant $\alpha=\sum_g\big\|[H_g,\sum_{g'>g}H_{g'}]\big\|$ (the larger of the two product orders; the ancilla does not change it). It is evaluated exactly, with no further inequality: for Pauli strings $[P_g,P_{g'}]=2P_gP_{g'}$ when they anticommute and $0$ otherwise, so each term is $2|c_g|\,\|\sum_{g'>g,\ \text{anticomm.}}c_{g'}P_{g'}\|$, a tiny dense eigenvalue problem on the few qubits near term $g$, with cost linear in the number of terms for any $N$. Three rules choose $n$:

* **Certified.** The smallest $n$ with $\bar D(n)\le\sqrt\varepsilon$, where the *certificate* is
$$\boxed{\;D\le\bar D:=\sqrt{\bar\ell}+\frac{\alpha T^2}{2n\sqrt{\gamma F(0)^2}},\qquad\bar\ell=\frac{(1-\gamma)\eta^2}{\gamma+(1-\gamma)\eta^2}\;}$$
$D$ is the trace distance between the implemented and the exact ground state. The first term is **leakage** (what an exact filter would leave) and depends only on $\gamma$ and $\eta$; the second is the **Trotter** term. Equivalently, infidelity $\le\bar D^2$ and energy error $\langle H\rangle-E_0\le W\bar D^2$. The proof (Appendix A) is the triangle inequality for the trace distance plus a leakage lemma and a first-order Trotter lemma.
* **Hybrid.** Keep the rigorous leakage term but replace the worst-case Trotter term by a measured one. For each candidate $n$ on a geometric grid, run the same filter at $n$ and $2n$ steps; because the first-order error scales as $1/n$, the Trotter distance is $\approx2D(\tilde\varphi_n,\tilde\varphi_{2n})$. Accept the smallest $n$ such that $\sqrt{\bar\ell}+2D(\tilde\varphi_n,\tilde\varphi_{2n})\le\sqrt\varepsilon$ at that grid point and the next. This uses no ground state and is *not* a theorem (it assumes the $1/n$ regime), but it is what a practitioner can actually evaluate.
* **Measured (oracle).** The smallest $n$ whose measured $D$ meets the target. It needs the exact ground state, so it is a reference, not a rule.

**The certificate is conditional on the DMRG inputs, and this can be made explicit.** Let $E_s=E_0^{\rm DMRG}$ be the shift, $\varepsilon_0=E_s-E_0\ge0$ (DMRG is variational), $\mu_0=\min_{\varphi\perp\psi_0}\langle\varphi|H|\varphi\rangle$ and $\varepsilon_1=E_1^{\rm DMRG}-\mu_0$, the optimisation error of the excited-state solve. Assume user-supplied bounds $\varepsilon_0\le\bar\varepsilon_0$, $\varepsilon_1\le\bar\varepsilon_1$ (neither is computable from DMRG alone). Courant–Fischer gives $E_0\le\mu_0\le E_1$ for *every* reference $\psi_0$ (Appendix A.2), so with $\delta_i=\bar\varepsilon_i/W$ the scaled spectrum obeys: ground space in $[-\delta_0,0]$, everything else in $[\Delta_*,1]$ with $\Delta_*=\Delta_{\rm D}-\delta_1$. The weight of $\psi_0$ outside the ground space is at most $s^2=\bar\varepsilon_0/(W\Delta_*)$, and with $\theta_1=\arcsin s$, $\theta_2=\arccos|\langle\psi_0|\psi\rangle|$, the trial state has ground-space weight $\gamma\ge\gamma_*=\cos^2(\theta_1+\theta_2)$. The theorem is then: *designing with $\eta_*=\sup_{[\Delta_*,1]}|F|/(F(0)-L\delta_0)$ ($L=\sum|t_i|$), $p_g=\gamma_*(F(0)-L\delta_0)^2$ and $\gamma_*$ in place of $\gamma$, the certificate above bounds the weight of the output outside the ground space.* The experiments set $\bar\varepsilon_0=\bar\varepsilon_1=0$ and take $\gamma$ as the overlap with the DMRG state, so every certified number below is the theorem evaluated at zero reference error, i.e. **conditional on the DMRG shift, gap and reference state being exact**, not unconditional.

**Stage 5 — run and check.** Execute the filter; check the post-selection rate, which is at least $p_g$, and, where the energy is measurable, $\langle H\rangle-E_0\le W\bar D^2$. If the device is noisy and prepares a mixed state $\tilde\rho$, the triangle inequality adds a noise term: $D(\tilde\rho,g)\le D(\tilde\rho,\tilde\varphi)+\bar D$.

**What is rigorous and what is assumed.**

| Quantity | Source | Status |
|---|---|---|
| $E_0$ | DMRG | variational ($\varepsilon_0\ge0$); accurate to $\sim10^{-13}$ here; enters the theorem through $\bar\varepsilon_0$ |
| $W$ | Weyl bound on the bonds | **rigorous**; equals $E_{\rm top}-E_0$ for $J_1,J_2\ge0$; no DMRG value of $E_{\rm top}$ needed |
| $E_1$, $\Delta$ | DMRG, orthogonality penalty | $\Delta_{\rm true}\ge\Delta_{\rm D}-\varepsilon_1/W$ by interlacing, so conditional on $\bar\varepsilon_1$ (zero in the runs) |
| $\gamma$ | overlap with a high-$\chi$ reference MPS | $\gamma\ge\gamma_*$ by the purity lemma, conditional on $\bar\varepsilon_0,\bar\varepsilon_1$ (zero in the runs) |
| $\eta$, $F(0)$ | interval certificate | rigorous in real arithmetic given $\Delta$; evaluated in double precision, **no interval arithmetic** (rounding $\lesssim10^{-14}\ll$ the pruning tolerance) |
| $\alpha$ | exact commutator sum | rigorous |
| ground-state degeneracy | — | the target is the ground space $\Pi_0$ and $\gamma$ its total weight; the chains studied have a nondegenerate singlet ground state (a degenerate space would show up as $\Delta\approx0$) |
| leakage and Trotter terms | Lemmas 1 and 2 | rigorous *given* the inputs above |

The certificate is therefore as trustworthy as $\Delta$ and $\gamma$; Section 4.3 quantifies what happens when they are wrong.

## 3. Experiments

* **System.** Open spin-½ $J_1$–$J_2$ chain, $J_1=1$, $J_2\in\{0,0.4\}$, $N\in\{4,6,\dots,16\}$. Bonds are ordered: nearest-neighbour bonds with even then odd left site, then next-nearest-neighbour bonds likewise.
* **Front end.** Two-site DMRG (`quimb` [14]) with the repository's `get_energies.py`: SVD cutoff $10^{-14}$, energy tolerance $10^{-12}$, inner Lanczos tolerance $10^{-10}$, bond dimensions up to 256 (the cutoff, not the cap, limits the bond dimension reached: $\chi\le114$ at $N\le32$), $E_1$ by the rank-one orthogonality penalty of Section 2. The workflow sees only DMRG quantities; the bandwidth is the Weyl bound. Exact diagonalisation is used afterwards to validate (true ground state, gap, overlap) and to measure distances.
* **Trial state.** DMRG at maximum bond dimension $\chi=2$ (main experiments), compiled exactly. $\chi\in\{1,2,4,8\}$ and $L$-layer approximate circuits are compared at $N=8$.
* **Targets.** $\varepsilon=10^{-1},10^{-2}$ ($N\le16$) and $10^{-3}$ ($N\le8$), i.e. $D\le0.32,0.1,0.032$.
* **Simulation.** Exact statevector simulation of the Trotterised filter, the ideal filter, and the ground state. A run is simulated only if $n\times\text{bonds}\times2\times2^N\le1.5\times10^{10}$; the $N=16$, $\varepsilon=10^{-2}$ certified circuits are not simulated, and their step counts and CX are reported from the certificate.
* **Larger chains ($N\le32$, Section 7.4).** The pipeline in `j1j2_filter/` runs DMRG and the circuits without any statevector or dense operator above $N=14$: $\gamma$ is the overlap of the DMRG MPS with the trial circuit simulated as an MPS (Aer), the trial state is an $L=1$ brickwork `mps-to-circuit` compilation, $J_2\in\{0,0.2411,0.4\}$, $\varepsilon\in\{10^{-1},10^{-2}\}$, and the circuits are *counted, not simulated* (per-step transpile times $n$, validated against a full transpile at $n=48$). Trotter terms are the individual Pauli rotations ($\alpha$ at Pauli level), so its step counts are not directly comparable with the bond-gate counts of Sections 4–7.3.
* **Cost.** One Trotter step of one pulse (for each bond, the $XX,YY,ZZ$ rotations tensored with $Z$ on the ancilla) is transpiled once with qiskit (level 3) to $\{\mathrm{CX},\mathrm{Rz},\mathrm{H},\mathrm{S}\}$; total CX is trial CX plus $n\times$ CX per step. The fault-tolerant metric reported for the larger chains is the number of *non-Clifford* $R_z$ rotations (angles not a multiple of $\pi/2$); no $T$-count is reported, because it depends on the synthesis method and its constants. All-to-all connectivity unless stated. The trial circuit is the exact `mps-to-circuit` compilation, transpiled the same way.
* **Direct-preparation baseline.** Preparing the state directly from the MPS: `mps-to-circuit` sequential (exact) circuits for the reference MPS truncated to $\chi\in\{1,2,4,8,16\}$, $L$-layer approximate circuits ($L\le8$), and generic state preparation of the exact ground state ($N\le12$). Each is measured against the exact ground state; the baseline for a target $\varepsilon$ is the cheapest with infidelity $\le\varepsilon$.

## 4. Q1: Can the filter be certified without knowing the exact ground state?

**Short answer: yes, conditionally.** Every quantity in the certificate, $\gamma$, $\Delta$, $W$, $\eta$, $F(0)$, $\alpha$ and $T$, can be computed from DMRG output and the circuit, with no access to the exact ground state. The certificate is then a theorem *given* $\gamma$ and $\Delta$ (Appendix A), and it is only as safe as those two inputs.

### 4.1 The classical front end supplies the inputs

{{inputs}}

*Table 1.* DMRG-derived inputs against exact diagonalisation. $\chi_{\mathrm{ref}}$ is the bond dimension the reference MPS used.

DMRG reproduces $E_0$, $E_1$ and the gap to better than {{max_E}} (absolute) and {{max_gap}} (relative), and the estimated $\gamma$ differs from the exact overlap by at most {{max_gam}} ($\gamma$ is first order in the state error, so it is the least accurate of the inputs). At these sizes the certificate computed from DMRG inputs is therefore indistinguishable from the one computed from exact inputs. The accuracy depends on the DMRG settings: with a looser SVD cutoff ($10^{-10}$) and sweep tolerance ($10^{-10}$), and quimb's default inner-solver tolerance ($10^{-3}$), the same calculation is limited to about $10^{-9}$, which propagates into $\Delta$ and $\gamma$ and so into the uncertified columns; the settings of Section 3 remove that.

Beyond $N=16$ the same front end was validated in two ways. Against exact diagonalisation (sparse Lanczos, used for testing only) for $N=16$–$20$:

{{dmrg_ed}}

*Table 1b.* DMRG against exact diagonalisation, $N=16$–$20$. All errors are at most {{ed_max}}.

and, with no exact reference for $N=24$–$32$, by two independent runs from different random starts, plus the energy variance $\langle H^2\rangle-\langle H\rangle^2$ of each state (a state with variance $s^2$ has an eigenvalue within $s$ of its energy):

{{dmrg_cons}}

*Table 1c.* Run-to-run agreement of DMRG at $N=24$–$32$. Energies agree to {{cons_max}} or better and the variances are at most {{cons_var_max}}. The trial-state overlap $\gamma$ computed by MPS contraction agrees with the dense overlap with the exact ground state to {{gamma_diff_max}} ($N=16$, $18$).

Run-to-run agreement shows that the solves are converged and reproducible; it does not prove that they sit at the true minimum, which is why Section 2 states the certificate as conditional. This is an empirical statement about well-behaved chains: DMRG does not certify these numbers.

### 4.2 The certificate computed from DMRG inputs was never violated

Across {{npts}} simulated points (every case at the certified step count plus the $n$-sweeps of Section 5), with the certificate computed from DMRG inputs only, the total certificate, its leakage term, its Trotter term, the triangle inequality, the energy bound $\langle H\rangle-E_0\le W\bar D^2$ and the success-rate bound $P_{\mathrm{succ}}\ge\gamma F(0)^2$ each held {{npts}} out of {{npts}} times. The exact ground state was used only afterwards, to measure $D$. This checks the derivation and the code; the proofs are in Appendix A.

### 4.3 What happens when the inputs are wrong

The condition in the short answer is not a formality. We design and certify with deliberately wrong values and measure the true distance ($N=8$ and $12$, $J_2=0$, $\varepsilon=10^{-2}$):

{{sens}}

*Table 2.* "Claimed" is the certificate computed from the (wrong) inputs; "with true inputs" recomputes it with the true gap and overlap; "measured" is the actual distance.

* **Under-estimating the gap is safe**: the filter is designed for a wider window and costs more steps.
* **Over-estimating the gap is not.** At $+10\%$ the claimed bound still holds, but the bound recomputed with the true gap is already $30$–$45\%$ above the claim; at $+25\%$ the claim is violated and the target missed at both sizes.
* **Over-estimating $\gamma$** behaves the same way, with a threshold that depends on how much excited weight is hidden.
* **Practical rule.** By interlacing (Section 2), the penalised solve targets $\mu_0\le E_1$, so the DMRG gap can overestimate the true gap only by the solve's own optimisation error $\varepsilon_1$; DMRG's $E_0$ is above the true one, which also works in the safe direction. $\varepsilon_1$ is not computable, but here it was at most {{max_E}}, so a gap error of $10$–$25\%$ is far outside what the front end produces. The workflow should still design with a deliberately *reduced* gap $\Delta_*=\Delta_{\rm D}-\bar\varepsilon_1/W$, or with a relative margin: a $20\%$ reduction was safe and cost {{gapm20_lo}}–{{gapm20_hi}}$\times$ more steps at both sizes.

**Answer to Q1.** Yes: the certificate needs no ground state, and it held on every point when built from DMRG inputs. It is conditional on $\gamma$ and $\Delta$ not being optimistic. The workflow bounds the error of those two numbers in terms of the DMRG errors $\bar\varepsilon_0,\bar\varepsilon_1$ (Section 2) but cannot compute those errors; the experiments quantify the cost of getting $\gamma$ and $\Delta$ wrong.

## 5. Q2: How conservative is the certificate?

**Short answer: the leakage term is reasonably tight, the worst-case Trotter term is extremely conservative.**

{{main}}

*Table 3.* $\varepsilon=10^{-2}$, $\chi=2$. $^\dagger$$1-\gamma\le\varepsilon$: no filter needed, and the entry is the trial state's own distance. $D$ is the measured distance at $n_{\mathrm{bound}}$. The last two columns compare the certificate with the measurement, term by term. $n_{\mathrm{hyb}}$ is discussed in Section 6.

* **Close in total.** At $n_{\mathrm{bound}}$ the certificate exceeds the measured $D$ by {{tot_lo}}–{{tot_hi}}$\times$ (median {{tot_med}}$\times$).
* **Leakage is reasonably tight.** The leakage term is within {{leak_lo}}–{{leak_hi}}$\times$ of the measured leakage distance. The only slack is that every excited component is charged the worst-case filter value $\eta F(0)$.
* **Trotter is extremely conservative.** The Trotter term is {{trot_lo}}–{{trot_hi}}$\times$ too large and the discrepancy grows with $N$ (from $9\times$ at $N=6$ to $52$–$77\times$ at $N=14$), because $\alpha$ is a worst case over all states. The measured $D$ at $n_{\mathrm{bound}}$ is dominated by leakage ({{leak_meas_lo}}–{{leak_meas_hi}}, against a Trotter part of at most {{trot_meas_max}}).
* **Consequence for step counts.** The certified step count is {{nbm_lo}}–{{nbm_hi}}$\times$ (median {{nbm_med}}$\times$) the smallest count that empirically meets the target.
* **Scaling in $n$ is right, the constant is not.** The Trotter certificate and the measured Trotter distance both fall as $1/n$ (fitted log–log slopes {{slope_b_lo}} to {{slope_b_hi}} and {{slope_m_lo}} to {{slope_m_hi}}), so the certificate has the right scaling and is off by a roughly constant factor at each $N$ (Figure 2). The measured total saturates at the leakage floor.
* **Other observables.** The energy bound $W\bar D^2$ holds but is {{en_lo}}–{{en_hi}}$\times$ above the measured energy error. The success rate matches the lower bound $\gamma F(0)^2$ to {{psucc_max}}% or better, so $1/(\gamma F(0)^2)$ is an accurate estimate of the repetition overhead.

![Figure 2. Dependence on $n$ at fixed design ($\varepsilon=10^{-2}$, $J_2=0$, $N=4,6,8,10$ as solid, dashed, dotted, dash-dotted). (a) Certificate (blue) and measured distance (red). (b) Trotter part only: certificate (blue) and measurement (green).](figs/fig_sweep.png){width=100%}

Step counts also scale with the target roughly as first-order theory predicts for the certificate ($n\propto\varepsilon^{-1/2}$ at a fixed filter; the observed ratios $n(10^{-3})/n(10^{-2})$ of $4$–$7$ include the tightening of the filter itself):

{{eps}}

*Table 4.* Certified step counts at three targets, and hybrid and measured counts at $10^{-2}$ and $10^{-3}$.

**Answer to Q2.** The certificate is within a factor of about {{tot_med}} in distance, because leakage dominates the error and is bounded tightly. The price of a worst-case proof is in the Trotter term: it is one to two orders of magnitude too large, increasingly so with $N$, which turns directly into an order-of-magnitude excess in circuit length (Section 7).

## 6. Q3: Can the practical cost be recovered without the exact ground state?

**Short answer: largely, with the $n$-versus-$2n$ hybrid rule, which is empirical rather than rigorous.** The hybrid rule keeps the rigorous leakage term and replaces the worst-case Trotter term by a measured one (Section 2, Stage 4). It never uses the ground state.

{{hybrid}}

*Table 5.* The hybrid rule against the certificate and the oracle. $D$ is measured at $n_{\mathrm{hyb}}$ against the exact ground state, which the rule itself never uses.

* **It met the target** in {{hyb_met}} of {{hyb_cases}} cases, with a margin $\sqrt\varepsilon/D$ of {{hmarg_lo}}–{{hmarg_hi}}$\times$.
* **It is close to the oracle.** It uses {{nbh_lo}}–{{nbh_hi}}$\times$ (median {{nbh_med}}$\times$) fewer steps than the certificate and only {{nhm_lo}}–{{nhm_hi}}$\times$ (median {{nhm_med}}$\times$) more than the oracle. The gap to the oracle is mostly the coarse geometric grid and the requirement of confirmation at the next grid point.
* **The Trotter estimate is accurate.** The $n$-versus-$2n$ estimate of the Trotter distance matched the measured one to within [{{rich_lo}}, {{rich_hi}}] over {{rich_n}} sweep points with $n\ge n_{\mathrm{meas}}$.
* **It is not a theorem.** It assumes the $1/n$ error regime, which held in every run here (Section 5), and it inherits the assumptions on $\gamma$ and $\Delta$ through its leakage term, so Section 4.3 applies to it unchanged. The evidence is 25 cases on one chain.

**Answer to Q3.** Yes, to a good approximation: the hybrid rule recovers all but a factor of about {{nhm_med}} of the empirically sufficient cost with no ground state, and met every target. What it gives up is the proof: its Trotter part is a measurement under an assumed scaling, not a bound.

## 7. Using the workflow: trial state, design choices and cost

### 7.1 The choice of trial state

At $N=8$, $\varepsilon=10^{-2}$:

{{trials}}

*Table 6.* Trial state against total CX (trial plus filter). "Approx." is the $L$-layer `mps-to-circuit` brickwork compilation of the reference MPS.

A better trial state needs a much shorter filter: going from $\chi=1$ ($\gamma=0.13$ at $J_2=0$) to $\chi=2$ cuts the certified step count by $7.7\times$ ($31\times$ at $J_2=0.4$, where $\chi=1$ gives $\gamma=0.04$). For $\chi=4$ the DMRG state already meets the target and no filter is needed. A trial state produced by approximate compilation behaves like a DMRG state of the same $\gamma$: the filter only sees $\gamma$.

### 7.2 Sensitivity to design choices

At $N=8$, $J_2=0$, $\varepsilon=10^{-2}$, certified $n$ against the leakage share $\varepsilon_\ell/\varepsilon$ and the number of pulses $m$ (the solver searched $x=T\Delta/\pi\in\{0.6,1,1.5,2\}$):

{{design}}

*Table 7.* Certified step count against design choices. Four pulses are infeasible for this trial state.

The certified $n$ varies by $3\times$ ({{des_lo}}–{{des_hi}}) over the shares tried, with a minimum near the $1/4$ used elsewhere; six and eight pulses are equivalent. The leakage share is the one design parameter worth tuning.

### 7.3 Gate cost

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

### 7.4 Larger chains: the DMRG/MPS-only regime ($N\le32$)

Above $N=14$ the workflow needs no statevector: DMRG supplies $E_0,E_1$ and the reference MPS, the trial circuit is simulated as an MPS to get $\gamma$, the filter is designed and certified from $(\gamma,\Delta,\alpha)$, and the circuits are counted. Guaranteed designs at $\varepsilon=10^{-2}$ (pipeline of `j1j2_filter/`, $L=1$ trial circuit, Pauli-level Trotter terms; $J_2=0.2411$ is the BKT point):

{{large}}

*Table 10.* Guaranteed designs, $\varepsilon=10^{-2}$, all-to-all connectivity; $T$ in units of $1/J_1$; totals include the trial circuit and are counted from the per-step transpile. $^\ddagger$: the floor solver found no feasible filter and the generic minimax builder was used. "Not achievable": no certified design within at most 8 pulses.

* **Reach.** The guaranteed circuits for $N=16$–${{large_Nmax}}$ need {{lg_n_lo}}–{{lg_n_hi}} Trotter steps, {{lg_cx_lo}}–{{lg_cx_hi}} CX and {{lg_rz_lo}}–{{lg_rz_hi}} non-Clifford $R_z$, with a physical evolution time of {{lg_T_lo}}–{{lg_T_hi}} (units of $1/J_1$).
* **Scaling.** Total CX grows as $N^{{{lgfit0_cx}}}$ ($J_2=0$, fit over $N={{lgfit0_N}}$), $N^{{{lgfit1_cx}}}$ ($J_2=0.2411$) and $N^{{{lgfit2_cx}}}$ ($J_2=0.4$, $N={{lgfit2_N}}$); the steps themselves as $N^{{{lgfit0_n}}}$–$N^{{{lgfit1_n}}}$. This is the certified count with a worst-case Trotter term; Section 5 shows that term is one to two orders of magnitude pessimistic.
* **Where the certificate runs out.** For $N=28$ and $32$ at $J_2=0$ and $0.2411$ the $L=1$ trial circuit has only $\gamma={{lg_gmin}}$–${{lg_gmax}}$, and no filter with at most 8 pulses reaches $\varepsilon=10^{-2}$ ($\eta_{\rm target}\approx0.07$ would be needed): at $\varepsilon=10^{-1}$ the same cases are feasible. This is a limit of the trial circuit, not of the filter; a deeper brickwork circuit raises $\gamma$.
* **The bandwidth does not matter for cost.** Replacing the earlier Pauli-$\ell_1$ bandwidth (1.6–1.9$\times$ larger than $E_{\rm top}-E_0$) by the Weyl bandwidth changed $n$ by a factor of {{bw_lo}}–{{bw_hi}} (median {{bw_med}}) and the physical time by {{bwT_lo}}–{{bwT_hi}} over {{bw_count}} matching designs; the spread is the discreteness of the design search, not a systematic effect, as expected from $n\propto\alpha/g^2$.

## 8. Discussion

**The three answers together.** The filter can be certified without the exact ground state, as long as $\gamma$ and $\Delta$ are not optimistic (Q1); the certificate is accurate where it matters for the distance, because leakage dominates and is tight, but its Trotter half is one to two orders of magnitude too large (Q2); and a measured $n$-versus-$2n$ Trotter estimate recovers nearly all of the practical cost without the ground state, at the price of being empirical (Q3). The recommended practice is therefore: certify $\gamma$, $\Delta$ and the leakage rigorously, choose $n$ by the hybrid rule, and report the certified count as the conservative fallback.

**What it costs.** At $\varepsilon=10^{-2}$ the hybrid circuit costs {{hyb_cx_lo}}–{{hyb_cx_hi}} CX at $N=6$–$16$ against $30$–$1{,}100$ for direct MPS preparation, and the certified circuit for $N=16$–${{large_Nmax}}$ costs {{lg_cx_lo}}–{{lg_cx_hi}} CX. The filter is not competitive on these chains.

**When it could matter.** The filter's cost is set by the gap and the commutator constant, not by the state's entanglement; direct preparation's cost steps up each time the bond dimension must double. The data show the ratio narrowing at one such step ($N=16$, $J_2=0$) but not crossing. A genuine crossover would need systems where the required MPS bond dimension grows quickly with size (critical systems at larger $N$, or two dimensions) and where the filter can be run without paying for dense simulation. Neither is tested here.

**Limitations.**

* *Inputs.* The certificate needs $\gamma$ and $\Delta$, which come from DMRG. They are validated against exact diagonalisation to $N=20$ (errors $\lesssim4\times10^{-13}$) and by run-to-run agreement to $N=32$, but nothing certifies them: a rigorous lower bound on the gap is not available from DMRG (a variational run cannot bound $E_1$ from below; interlacing bounds the *direction* of the error, not its size), and Section 4.3 shows the certificate breaks if the gap is optimistic. The results are therefore *conditional certificates* (Section 2). Obtaining a lower bound on $\mu_0$ or on the gap independently is the main open point.
* *Models.* Three couplings of one chain ($J_2=0$, $0.2411$, $0.4$) at larger $N$, two at $N\le16$. The ground state is a nondegenerate singlet; the theorem is stated for a ground space, but a degenerate target is not tested.
* *Interval arithmetic.* The bound on $\eta$ is computed in double precision without outward rounding; its rounding error ($\lesssim10^{-14}$) is far below the pruning tolerance, but it is not a machine-checked proof.
* *Trotter order.* First order only; higher-order formulas would shrink the Trotter term but need a different constant.
* *Noise.* Not simulated; the triangle inequality gives an additive noise term.
* *Cost model.* Gate counts of an ideal circuit: CX (all-to-all; a line costs about $2$–$3\times$ more) and, for the larger chains, non-Clifford $R_z$; no $T$-count and no synthesis error. The expected repetition count $1/p_g$ assumes independent attempts. Ancilla measurement and reset, classical feed-forward for the early abort, the native implementation of controlled Pauli rotations on a real connectivity, and any error-mitigation or error-correction overhead are not included. CX counts come from a heuristic compiler and are not minimal.
* *Simulation.* Dense statevectors to $N=16$; for $N=18$–$32$ the circuits are counted, not simulated, so only the certified bound is reported there and the hybrid and measured step counts are not available.
* *Hybrid rule.* An empirical rule that assumed the $1/n$ regime and had 25 test cases.


**Reproducibility.** All code is in `simple_paper/`: `dmrg_inputs.py` (DMRG front end and trial circuits), `workflow.py` (Hamiltonian, pulses, certificate), `run_experiments.py` (main grid), `run_extras.py` (trial states, wrong inputs, design choices, direct-preparation baseline), `resources.py` (CX counts), `pinsker_check.py`, `make_report.py`, `build_paper.py`; filter design is `j1j2_filter/floor.py`. They run in a Python 3.12 environment with `numpy`, `scipy`, `quimb`, `qiskit` and `mps-to-circuit`; set `OMP_NUM_THREADS=1` when running processes in parallel. The full set (all tables here) runs in one to two hours on four cores. Seeds are fixed, including the Lanczos start vector used for the exact ground state at $N\ge12$ (a $\chi=2$ cut passes through degenerate Schmidt values, so $\gamma$ is otherwise only reproducible to the third digit).

## References

1. I. Stetcu, A. Baroni, J. Carlson, "Projection algorithm for state preparation on quantum computers," *Phys. Rev. C* **108**, L031306 (2023), arXiv:2211.10545.
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
12. K. Okamoto, K. Nomura, "Fluid-dimer critical point in $S=\tfrac12$ antiferromagnetic Heisenberg chain with next nearest neighbor interactions," *Phys. Lett. A* **169**, 433 (1992).
13. S. Eggert, "Numerical evidence for multiplicative logarithmic corrections from marginal operators," *Phys. Rev. B* **54**, R9612 (1996), arXiv:cond-mat/9602026 (reports $J_2^c=0.241167\pm0.000005$).
14. J. Gray, "quimb: a python package for quantum information and many-body calculations," *J. Open Source Softw.* **3**, 819 (2018).
15. R. A. Horn, C. R. Johnson, *Matrix Analysis*, 2nd ed. (Cambridge University Press, 2013), Thm. 4.2.6 (Courant–Fischer) and Thm. 4.3.1 (Weyl).

## Appendix A. Proof of the certificate

### A.1 The certificate

Write $D(\rho,\sigma)=\tfrac12\|\rho-\sigma\|_1$. For pure states $D=\sqrt{1-|\langle a|b\rangle|^2}=\sin\theta_{ab}$, with $\theta_{ab}$ the angle between the rays. $D$ obeys the triangle inequality. For any state $\rho$ and pure $g$, $1-\langle g|\rho|g\rangle\le D(\rho,g)\le\sqrt{1-\langle g|\rho|g\rangle}$; for pure $\rho$ the right side is an equality, so infidelity is $D^2$. Let $\varphi=F(H_s)\psi/\|F(H_s)\psi\|$ be the ideal filtered state and $\tilde\varphi$ the implemented one. The certificate follows from $D(\tilde\varphi,g)\le D(\tilde\varphi,\varphi)+D(\varphi,g)$ and two lemmas.

**Lemma 1 (leakage).** $D(\varphi,g)^2=\ell\le\bar\ell$.
*Proof.* The ground state sits at $E=0$, so $|\langle g|\varphi\rangle|^2=\gamma F(0)^2/(\gamma F(0)^2+y)$ with $y=\sum_{k\ge1}|c_k|^2F(E_k)^2$. All $E_k$, $k\ge1$, lie in $[\Delta,1]$, so $F(E_k)^2\le\eta^2F(0)^2$ and $y\le(1-\gamma)\eta^2F(0)^2$. Since $\ell=y/(\gamma F(0)^2+y)$ is increasing in $y$, the bound follows. $\square$

**Lemma 2 (Trotter).** $D(\tilde\varphi,\varphi)\le\min\{1,\epsilon_T/\sqrt{p_g}\}$ with $\epsilon_T=\sum_i\alpha t_i^2/(2k_i)=\alpha T^2/(2n)$.
*Proof.* Let $A_i=\langle0_{\mathrm{a}}|\mathcal W_i|0_{\mathrm{a}}\rangle$ and $\tilde A_i$ the Trotterised block; $a=A_m\cdots A_1\psi$, $b=\tilde A_m\cdots\tilde A_1\psi$. The first-order product-formula bound (below) gives $\|\mathcal W_i-\widetilde{\mathcal W}_i\|\le\alpha t_i^2/(2k_i)$; the blocks are compressions, so $\|A_i-\tilde A_i\|$ is no larger, and all blocks have norm $\le1$. Telescoping, $a-b=\sum_i(A_m\cdots A_{i+1})(A_i-\tilde A_i)(\tilde A_{i-1}\cdots\tilde A_1\psi)$, so $\|a-b\|\le\epsilon_T$. The rays of $a,b$ are $\varphi,\tilde\varphi$, so $D=\sin\angle(a,b)=\mathrm{dist}(a,\mathrm{span}\,b)/\|a\|\le\|a-b\|/\|a\|$, and $\|a\|^2=P_{\mathrm{succ}}\ge p_g$. $\square$

**Corollaries.** Infidelity $=D^2\le\bar D^2$. Energy: with $\tilde c_k$ the components of $\tilde\varphi$, $\langle H\rangle-E_0=\sum_{k\ge1}|\tilde c_k|^2(E_k-E_0)\le W(1-|\tilde c_0|^2)=WD^2$. Step count: for $\bar\ell\le\varepsilon/4$, $\bar D\le\sqrt\varepsilon$ holds if $n\ge\alpha T^2/\sqrt{\gamma F(0)^2\varepsilon}$; with $T=x\pi/\Delta$ the certified cost scales as $\alpha/\Delta^2$.

**First-order Trotter bound.** For Hermitian $A,B$, $U(s)=e^{-is(A+B)}$, $V(s)=e^{-isA}e^{-isB}$: $\frac{d}{ds}[U(\delta-s)V(s)]=iU(\delta-s)(B-e^{-isA}Be^{isA})V(s)$, so $U(\delta)-V(\delta)=-i\int_0^\delta U(\delta-s)(B-e^{-isA}Be^{isA})V(s)\,ds$. Since $\|B-e^{-isA}Be^{isA}\|\le s\|[A,B]\|$ and $U,V$ are unitary, $\|U(\delta)-V(\delta)\|\le\tfrac{\delta^2}2\|[A,B]\|$. For many terms apply this with $A=H_1$, $B=\sum_{g>1}H_g$ and recurse, giving $\frac{\delta^2}{2}\alpha$ per step; with $k$ steps of size $t/k$ the error is $\alpha t^2/(2k)$ [6]. The reversed product order gives the reversed sum, hence the maximum over orders. For $H_s\otimes Z$: $[H_gZ,H_{g'}Z]=[H_g,H_{g'}]\otimes\mathbb 1$, so $\alpha$ is unchanged; the constant shift commutes with everything and is applied exactly.

### A.2 The bandwidth, the spectral window and the purity of the reference

**Bandwidth (Weyl).** For Hermitian $A,B$ and any unit $v$, $\langle v|(A+B)|v\rangle=\langle v|A|v\rangle+\langle v|B|v\rangle\le\lambda_{\max}(A)+\lambda_{\max}(B)$; maximising over $v$ gives $\lambda_{\max}(A+B)\le\lambda_{\max}(A)+\lambda_{\max}(B)$ [15]. A bond $J\,\mathbf S\!\cdot\!\mathbf S$ is $J/4$ on the triplet and $-3J/4$ on the singlet, so for $J\ge0$ its largest eigenvalue is $J/4$; the fully polarised state has $\mathbf S_i\!\cdot\!\mathbf S_j=1/4$ on every pair and attains the sum.

**Interlacing.** For any unit $\psi_0$ let $\mu_0=\min_{\varphi\perp\psi_0}\langle\varphi|H|\varphi\rangle$. Then $E_0\le\mu_0$ trivially and $E_1=\max_v\min_{\varphi\perp v}\langle\varphi|H|\varphi\rangle\ge\mu_0$ (Courant–Fischer), so $E_0\le\mu_0\le E_1$ for every $\psi_0$, accurate or not. The solve returns $E_1^{\rm DMRG}=\mu_0+\varepsilon_1$, hence $E_1-E_0\ge\mu_0-E_0=(E_1^{\rm DMRG}-E_s)-\varepsilon_1+\varepsilon_0\ge(E_1^{\rm DMRG}-E_s)-\varepsilon_1$, with $E_s=E_0^{\rm DMRG}=E_0+\varepsilon_0$. (The penalty is finite and the returned state has a residual overlap with $\psi_0$ of order $10^{-14}$; that difference is absorbed in $\varepsilon_1$, which is only required to be bounded above.)

**Window lemma.** With $H_s=(H-E_s)/W$ and $W=\sum_b\lambda_{\max}(h_b)-E_s$: the ground space sits at $-\varepsilon_0/W\in[-\delta_0,0]$; all other eigenvalues are $\ge(\mu_0-E_s)/W\ge\Delta_{\rm D}-\varepsilon_1/W\ge\Delta_*$ and $\le(E_{\rm top}-E_s)/W\le1$.

**Purity lemma.** Let $w_k$ be the weight of $\psi_0$ on level $E_k$. Then $\varepsilon_0=\sum_kw_k(E_k-E_0)\ge(E_1-E_0)(1-w_0)\ge W\Delta_*(1-w_0)$, so $\sin^2\theta(\psi_0,\Pi_0)\le s^2=\bar\varepsilon_0/(W\Delta_*)$. The angle $\theta(a,b)=\arccos|\langle a|b\rangle|$ is a metric on rays (Fubini–Study) and $\theta(\psi,\Pi_0)=\min_{v\in\Pi_0}\theta(\psi,v)$, so $\theta(\psi,\Pi_0)\le\theta(\psi,\psi_0)+\theta(\psi_0,\Pi_0)\le\theta_2+\theta_1$ and $\gamma=\cos^2\theta(\psi,\Pi_0)\ge\cos^2(\theta_1+\theta_2)=\gamma_*$ when $\theta_1+\theta_2\le\pi/2$.

**Theorem.** The proof of Appendix A.1 goes through with the window of the lemma: $|F(e_g)|\ge F(0)-L\delta_0$ for $e_g\in[-\delta_0,0]$ by $|F'|\le L$, every excited component obeys $|F|\le\sup_{[\Delta_*,1]}|F|$, and $\gamma\ge\gamma_*$, each entering the bounds monotonically. For a degenerate ground space, $g$ in Lemma 1 is the normalised $\Pi_0\varphi$ and the weight in $\Pi_0$ replaces the overlap with one vector; all ground components share the eigenvalue $e_g$.

**Floating point.** The interval bound $|F(E)|\le|F(\mu)|+|F'(\mu)|d+\tfrac12L^2d^2$ holds exactly in real arithmetic. In double precision each evaluation of $F$ and $F'$ carries a relative rounding error of order $m\,u$ ($u\approx1.1\times10^{-16}$, $\lesssim10^{-14}$ for $m\le12$), far below the pruning tolerance ($0.5$–$2\%$); a formal claim would need interval arithmetic with outward rounding, which is not implemented.

## Appendix B. Where Pinsker's inequality applies

For the target $|g\rangle\langle g|$ the relative entropy $S(\rho\|\,|g\rangle\langle g|)$ is infinite unless $\rho=|g\rangle\langle g|$, so $D\le\sqrt{S/2}$ gives nothing here, which is why the certificate uses the triangle inequality and fidelity instead. Pinsker does apply to the energy-resolved distributions $p_k=|c_k|^2$ (trial) and $q_k=p_kF(E_k)^2/P_{\mathrm{succ}}$ (filtered): $\mathrm{KL}(q\|p)=\sum_kq_k\ln(F(E_k)^2/P_{\mathrm{succ}})\le\ln(1/P_{\mathrm{succ}})$, and $\mathrm{TV}(p,q)\ge q_0-p_0=1-\ell-\gamma$, so with $\mathrm{TV}\le\sqrt{\mathrm{KL}/2}$,
$$P_{\mathrm{succ}}\le\exp\!\big(-2(1-\ell-\gamma)^2\big).$$
This is a necessary condition on the post-selection cost, not an error bound. We checked it on the {{pn}} $\varepsilon=10^{-2}$ cases with $N\le8$ that need a filter: every inequality held, but the cap ({{pk_lo}}–{{pk_hi}}) is far above the actual $P_{\mathrm{succ}}$ ({{pP_lo}}–{{pP_hi}}), so it is uninformative for these trial states.
