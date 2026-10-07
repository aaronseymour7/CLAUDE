---
title: "Filtering a Matrix-Product Trial State to the Ground State with the Stetcu–Baroni Projection: A Certified Workflow, Its Accuracy and Its Cost"
author: "Draft — authors to be added"
date: "7 October 2026"
---

## Abstract

DMRG gives a matrix-product trial state and the low-lying energies of a spin chain classically, but the trial state keeps some excited-state weight. The projection filter of Stetcu, Baroni and Carlson removes it on a quantum computer using one ancilla, controlled time evolution and post-selection. We describe the complete workflow, from trial state and energies to a filter and a Trotter step count, and give an a-priori guarantee on the trace distance $D$ between the prepared state and the exact ground state: $D\le\sqrt{\bar\ell}+\epsilon_T/\sqrt{p_g}$, a leakage term plus a Trotter term joined by the triangle inequality. We evaluate it on the open $J_1$–$J_2$ spin-½ chain ($N=4$–$12$, bond-dimension-2 trial states, exact simulation, 88 simulated points). The guarantee is never violated and lands within $2.7$–$4.9\times$ of the measured distance, because the leakage term is nearly tight ($1.3$–$2.5\times$). The Trotter term is conservative ($10$–$39\times$), so the guaranteed step count is $11$–$60\times$ the step count that is empirically sufficient. In gates, the guaranteed filter costs $2\times10^3$–$6\times10^5$ CX and the empirically sufficient one $2\times10^2$–$1.4\times10^4$, growing polynomially in $N$ ($\sim N^{4.6}$–$N^{5.2}$ and $\sim N^{2.8}$–$N^{3.9}$). The guaranteed filter costs more than exact generic state preparation at every size tested. The paper therefore documents a workflow that certifiably works and what it costs, not a resource advantage on these chains.

## 1. Introduction

A standard hybrid route to a ground state is to run DMRG, compile the resulting matrix-product state (MPS) into a shallow circuit, and then repair the circuit's residual error with a quantum projection step. The Stetcu–Baroni–Carlson (SBC) projection [1], related to the rodeo algorithm [2], does this by multiplying the trial state's energy components by a filter $F(E)=\prod_i\cos(Et_i+\phi_i)$, realised with pulses of controlled time evolution on one ancilla.

For someone deciding whether to use this, three questions matter. **Does it work?** **How accurate is the guarantee, and what does it assume?** **What does it cost?** This paper answers them for one concrete workflow:

1. Section 2 states the workflow and the guarantee. Everything the practitioner needs is one formula and a recipe.
2. Section 3 tests it. The guarantee is never violated; we show which part of it is tight and which is loose.
3. Section 4 reports the cost in Trotter steps and CX gates, how it scales with $N$, and how it compares with exact state preparation.
4. Section 5 states what the workflow does not do.

Proofs are short and sit in the appendix. We do not claim a resource advantage; Section 4 shows there is none for the sizes studied.

## 2. The workflow and its guarantee

![Figure 1. The workflow.](figs/fig_workflow.png){width=100%}

**Inputs.** A Hamiltonian written as a sum of terms $H=\sum_gH_g$ (here the bonds of the chain); the energies $E_0$, $E_1$, $E_{\mathrm{top}}$; and a trial state $|\psi\rangle=\sum_kc_k|E_k\rangle$. From these:

* the scaled Hamiltonian $H_s=(H-E_0)/W$ with $W=E_{\mathrm{top}}-E_0$, whose spectrum lies in $\{0\}\cup[\Delta,1]$ with $\Delta=(E_1-E_0)/W$;
* the trial-state ground-state weight $\gamma=|c_0|^2$.

**Filter.** One pulse $(t,\phi)$ applies $\mathsf H_{\mathrm{a}}\,e^{-itH_s\otimes Z_{\mathrm{a}}}\,\mathrm{Rz}_{\mathrm{a}}(2\phi)\,\mathsf H_{\mathrm{a}}$ to the system and an ancilla in $|0\rangle$; its ancilla-$|0\rangle$ block is $\cos(H_st+\phi)$. Keeping outcome $0$ after $m$ pulses applies, up to normalisation, $F(E)=\prod_{i=1}^m\cos(Et_i+\phi_i)$ with total time $T=\sum_it_i$ and $F(0)=\prod_i\cos\phi_i$. Two numbers characterise a filter: $\eta=\sup_{E\in[\Delta,1]}|F(E)|/|F(0)|$ (how strongly excited states are suppressed relative to the ground state) and $p_g=\gamma F(0)^2$, a lower bound on the success probability. We design $(t_i,\phi_i)$ with `floor.py`, which maximises $F(0)^2$ subject to a *certified* upper bound on $\eta$ (interval bisection using $|F'|\le\sum|t_i|$, $|F''|\le(\sum|t_i|)^2$).

**Implementation.** Each controlled evolution is replaced by $k_i$ first-order Lie–Trotter steps over the terms $H_g$, with $n=\sum_ik_i$ and $t_i=k_i\,dt$. The Trotter error is controlled by the commutator constant $\alpha=\sum_g\big\|[H_g,\sum_{g'>g}H_{g'}]\big\|$ (the larger of the two product orders; the ancilla does not change it).

**Guarantee.** Let $|\tilde\varphi\rangle$ be the implemented state, $g$ the ground state, and $D$ the trace distance. Then

$$\boxed{\;D(\tilde\varphi,g)\;\le\;\bar D:=\sqrt{\bar\ell}\;+\;\frac{\alpha T^2}{2n\sqrt{\gamma F(0)^2}},\qquad \bar\ell=\frac{(1-\gamma)\eta^2}{\gamma+(1-\gamma)\eta^2}\;}$$

The first term is **leakage**: what an exact filter would leave behind. It depends only on $\gamma$ and $\eta$. The second is the **Trotter** term and depends on $n$. Equivalent statements: infidelity $\le\bar D^2$ and energy error $\langle H\rangle-E_0\le W\bar D^2$. To reach a target $D\le\sqrt\varepsilon$ with leakage budget $\bar\ell\le\varepsilon/4$, it suffices that
$$n\;\ge\;\frac{\alpha T^2}{\sqrt{\gamma F(0)^2\,\varepsilon}},\qquad T=\frac{x\pi}{\Delta}\ (x=O(1)\text{ from the design}),$$
so the guaranteed cost scales as $\alpha/\Delta^2$. If the device is noisy and prepares a mixed state $\tilde\rho$, the triangle inequality adds a noise term: $D(\tilde\rho,g)\le D(\tilde\rho,\tilde\varphi)+\bar D$.

**Recipe.**

1. Compute $\gamma$ and $\Delta$. If $1-\gamma\le\varepsilon$ the trial state already meets the target and no filter is needed.
2. Design $(t_i,\phi_i)$ for a leakage budget $\varepsilon/4$ and certify $\eta$.
3. Compute $\alpha$ for the chosen term order.
4. For each candidate $n$, snap the pulses to integer step counts, re-optimise the phases, re-certify $\eta$, and evaluate $\bar D(n)$. Take the smallest $n$ with $\bar D(n)\le\sqrt\varepsilon$.
5. Run the filter and check the observables: infidelity, energy error, and the post-selection rate.

**Where the other tools in the toolbox fit.** The guarantee uses the trace distance, its triangle inequality, and the fact that for pure states $D=\sqrt{1-\text{fidelity}}$. Pinsker's inequality is not needed: relative entropy to a pure target is infinite, so it cannot bound $D$ here (Appendix B says where it does apply, and why it is weak).

## 3. Accuracy: what the guarantee delivers

### 3.1 Setup

* **System.** Open spin-½ $J_1$–$J_2$ chain, $J_1=1$, $J_2\in\{0,0.4\}$, $N\in\{4,\dots,12\}$. Terms are the bonds, ordered: nearest-neighbour bonds with even then odd left site, then next-nearest-neighbour bonds likewise. Exact $E_0,E_1,E_{\mathrm{top}}$ and the ground state come from diagonalisation.
* **Trial state.** The exact ground state truncated by successive SVDs to bond dimension $\chi=2$, as a stand-in for a bond-dimension-2 DMRG state ($\gamma=0.72$–$0.996$).
* **Filter.** $m\in\{4,6\}$ pulses, designed for $\varepsilon=10^{-2}$ ($D\le0.1$) with leakage budget $\varepsilon/4$, snapped to $n$ steps.
* **Measurement.** Exact statevector simulation of the Trotterised filter, of the exact filter, and of the ground state. We record $D$ and its two components.
* **Two step counts.** $n_{\mathrm{bound}}$ is the smallest $n$ whose guarantee is $\bar D\le0.1$. $n_{\mathrm{meas}}$ is the smallest $n$ whose *measured* $D$ is $\le0.1$. The latter uses the exact ground state, so it is an oracle: it shows what is sufficient, not something a practitioner can know in advance.

### 3.2 Results

| $N$ | $J_2$ | $\gamma$ | $\Delta$ | $n_{\mathrm{bound}}$ | $n_{\mathrm{meas}}$ | bound $\bar D$ | measured $D$ | leakage: bound / meas. | Trotter: bound / meas. |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 0.0 | 0.955 | 0.278 | 88 | 8 | 0.100 | 0.020 | 2.5$\times$ | 29$\times$ |
| 4 | 0.4 | 0.996 | 0.274 | – | – | – | 0.064$^\dagger$ | – | – |
| 6 | 0.0 | 0.898 | 0.131 | 375 | 18 | 0.100 | 0.037 | 1.3$\times$ | 20$\times$ |
| 6 | 0.4 | 0.988 | 0.128 | 405 | 29 | 0.100 | 0.026 | 2.2$\times$ | 10$\times$ |
| 8 | 0.0 | 0.837 | 0.077 | 992 | 37 | 0.100 | 0.033 | 1.5$\times$ | 16$\times$ |
| 8 | 0.4 | 0.977 | 0.075 | 1,075 | 47 | 0.100 | 0.024 | 2.1$\times$ | 19$\times$ |
| 10 | 0.0 | 0.776 | 0.050 | 2,183 | 74 | 0.099 | 0.033 | 1.5$\times$ | 20$\times$ |
| 10 | 0.4 | 0.964 | 0.050 | 2,256 | 59 | 0.100 | 0.021 | 2.5$\times$ | 29$\times$ |
| 12 | 0.0 | 0.717 | 0.036 | 6,967 | 117 | 0.099 | 0.023 | 2.2$\times$ | 39$\times$ |
| 12 | 0.4 | 0.949 | 0.036 | 4,060 | 93 | 0.100 | 0.022 | 2.3$\times$ | 31$\times$ |

*Table 1.* $\varepsilon=10^{-2}$. $^\dagger$$1-\gamma\le\varepsilon$: no filter needed, and the entry is the trial state's own distance. The last two columns are the ratio of the guarantee to the measurement for each term.

* **The guarantee holds.** Across all 88 simulated points (every case, at $n_{\mathrm{bound}}$ and along $n$-sweeps) the total guarantee, each of its two terms, the triangle inequality, the energy bound, and $P_{\mathrm{succ}}\ge\gamma F(0)^2$ were each satisfied 88 out of 88 times. This checks the derivation and the code; the proofs are in the appendix.
* **It is close in total.** At $n_{\mathrm{bound}}$ the guarantee exceeds the measured $D$ by $2.7$–$4.9\times$ (median $3.5\times$).
* **Leakage is nearly tight, Trotter is not.** The leakage term is within $1.3$–$2.5\times$. The Trotter term is $10$–$39\times$ too large, because $\alpha$ is a worst case over all states. Measured $D$ at $n_{\mathrm{bound}}$ is dominated by leakage ($0.02$–$0.04$ against a Trotter part of at most $0.005$).
* **Dependence on $n$.** Both the Trotter guarantee and the measured Trotter distance fall as $1/n$ (fitted slopes $-0.98$ to $-1.01$ and $-0.99$ to $-1.09$), so the guarantee has the right scaling and is off by a roughly constant factor (Figure 2). The measured total saturates at the leakage floor.
* **Other observables.** The energy bound $W\bar D^2$ holds but is $16$–$217\times$ above the measured energy error. The post-selection rate matches the lower bound $\gamma F(0)^2$ to $0.6$% or better, so $1/(\gamma F(0)^2)$ is an accurate estimate of the repetition overhead.
* **A non-oracle estimate of the Trotter part.** Comparing the same filter at $n$ and $2n$ steps, $\hat D_{\mathrm{Trot}}=2D(\tilde\varphi_n,\tilde\varphi_{2n})$ reproduces the measured Trotter distance to within $[0.84,1.07]$ over all 60 sweep points with $n\ge n_{\mathrm{meas}}$. It is an estimate (it assumes $1/n$ scaling), not a bound, but it needs no ground state. Rigorous leakage plus this estimate is the practical alternative to the full guarantee.

![Figure 2. Dependence on $n$ at fixed design ($\varepsilon=10^{-2}$, $J_2=0$, $N=4,6,8,10$ as solid, dashed, dotted, dash-dotted). (a) Total guarantee (blue) and measured distance (red). (b) Trotter part only: guarantee (blue) and measurement (green).](figs/fig_sweep.png){width=100%}

### 3.3 The trial state matters

At $N=8$, $J_2=0$, $\varepsilon=10^{-2}$:

| $\chi$ | $\gamma$ | $\sqrt{1-\gamma}$ (trial) | $T$ | $\eta$ | $n_{\mathrm{bound}}$ | $n_{\mathrm{emp}}$ | bound $\bar D$ | measured $D$ |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.1336 | 0.931 | 41.0 | 0.019 | 8519 | 93 | 0.0993 | 0.0330 |
| 2 | 0.8371 | 0.404 | 24.6 | 0.111 | 992 | 37 | 0.0996 | 0.0326 |
| 3 | 0.9110 | 0.298 | 24.6 | 0.157 | 877 | 37 | 0.0995 | 0.0214 |
| 4 | 0.9988 | 0.035 | – | – | no filter needed | – | – | 0.035 |

*Table 2.* Better trial states need shorter, cheaper filters; at $\chi=4$ ($\gamma=0.9988$) no filter is needed at this target. A poor trial state ($\chi=1$, $\gamma=0.13$) needs $\approx8.6\times$ the guaranteed steps of $\chi=2$.

## 4. Cost

One Trotter step of one pulse is, for each bond, the three Pauli rotations $XX,YY,ZZ$ tensored with $Z$ on the ancilla. We transpiled it with qiskit (level 3) to $\{\mathrm{CX},\mathrm{Rz},\mathrm{H},\mathrm{S}\}$. Ancilla Hadamards and phase rotations add no CX, so total CX is $n\times(\text{CX per step})$. The count was the same for 1 and 3 consecutive steps with all-to-all connectivity. Counts are CX only, noiseless, for the filter alone (not the trial circuit).

| $N$ | $J_2$ | CX / step | $n_{\mathrm{bound}}$ | CX at $n_{\mathrm{bound}}$ | $n_{\mathrm{meas}}$ | CX at $n_{\mathrm{meas}}$ | exact generic prep |
|---|---|---|---|---|---|---|---|
| 4 | 0.0 | 21 | 88 | 1,848 | 8 | 168 | 11 |
| 6 | 0.0 | 35 | 375 | 13,125 | 18 | 630 | 57 |
| 8 | 0.0 | 49 | 992 | 48,608 | 37 | 1,813 | 247 |
| 10 | 0.0 | 63 | 2,183 | 137,529 | 74 | 4,662 | 1,013 |
| 12 | 0.0 | 77 | 6,967 | 536,459 | 117 | 9,009 | 4,083 |
| 6 | 0.4 | 63 | 405 | 25,515 | 29 | 1,827 | 57 |
| 8 | 0.4 | 91 | 1,075 | 97,825 | 47 | 4,277 | 247 |
| 10 | 0.4 | 119 | 2,256 | 268,464 | 59 | 7,021 | 1,013 |
| 12 | 0.4 | 147 | 4,060 | 596,820 | 93 | 13,671 | 4,083 |

*Table 3.* $\varepsilon=10^{-2}$, $\chi=2$, all-to-all connectivity. Last column: qiskit's exact preparation of a generic $N$-qubit state, which uses no structure.

![Figure 3. CX count against $N$.](figs/fig_cx.png){width=60%}

* **Per step**, CX grows about linearly with the number of bonds ($N^{1.1}$ at $J_2=0$); $J_2=0.4$ costs about $1.9\times$ more because it has about twice as many bonds. A line with the ancilla at the end costs $1.4$–$2.6\times$ more than all-to-all.
* **Steps** grow as $N^{4.1}$ at the guarantee and $N^{2.7}$ in the measured case ($J_2=0$), driven by $T^2\propto\Delta^{-2}$: the scaled gap falls from $0.28$ to $0.036$ over $N=4$–$12$.
* **Total CX** at the guarantee grows as $N^{5.2}$ ($N^{4.6}$ at $J_2=0.4$) and at the measured step count as $N^{3.9}$ ($N^{2.8}$). The guarantee costs 11–60$\times$ the measured count.
* **Against exact generic preparation**, which grows $4\times$ per two qubits, the guaranteed filter is 131–448$\times$ more expensive at every size. Even the oracle measured filter is $15\times$ the generic cost at $N=4$ and $2.2\times$ at $N=12$ ($J_2=0$). The trend is toward a crossover just above $N=12$ *if* these power laws continued and the oracle step count were available, but they are fits to five points, and exact generic preparation is a weak baseline: structured MPS preparation is cheaper still for these chains, and a companion study in this repository (`evaluation/`) finds no advantage for the filter at $N\le12$.

**Summary of capability and cost.** The workflow produces a filter and a step count with a proof that $D\le\bar D$, and the proof is accurate to a factor of about 3–5 in $D$. What the proof costs is the Trotter term: the certified circuit is 11–60× longer than needed. The filter's cost is set by $1/\Delta^2$ and the commutator constant, *not* by how entangled the state is; it can therefore only compete where state preparation itself becomes expensive, which none of these chains are.

## 5. Limitations

* **Inputs are assumed known.** The guarantee requires $\gamma$, $\Delta$ and $E_0$. Here they come from exact diagonalisation, so everything is limited to $N\le12$. At larger $N$ DMRG supplies estimates, not bounds: the penalty-method first excitation is an *upper* bound on $E_1$, the wrong side for a gap lower bound, and an overlap with a DMRG state is an overlap with an approximation of $|E_0\rangle$. Certifying $\Delta$ and $\gamma$ at large $N$ is the main open step.
* **Trial states are a stand-in.** SVD truncations of the exact ground state, not compiled circuits; compilation lowers $\gamma$ further. The guarantee applies unchanged given the true $\gamma$.
* **First-order Trotter only,** noiseless, no gate-synthesis error.
* **Cost scope.** CX only, filter only, from fits over $N\le12$ and two couplings.
* **Dense simulation.** All checks use statevectors, so the study does not exercise DMRG's large-$N$ reach.

Reproducibility: all code is in `simple_paper/` (`workflow.py`, `run_experiments.py`, `resources.py` [needs qiskit], `pinsker_check.py`, `make_report.py`, `build_paper.py`); filter design is `j1j2_filter/floor.py`. Seeds are fixed. The $N=12$ ground state uses a fixed Lanczos start vector because a $\chi=2$ cut passes through degenerate Schmidt values and $\gamma$ otherwise varies in the third digit between runs. The full grid runs in a few minutes on one core.

## References

1. I. Stetcu, A. Baroni, J. Carlson, "Projection algorithm for state preparation on quantum computers," *Phys. Rev. C* **105**, 064308 (2022).
2. K. Choi, D. Lee, J. Bonitati, Z. Qian, J. Watkins, "Rodeo algorithm for quantum computing," *Phys. Rev. Lett.* **127**, 040505 (2021).
3. A. M. Childs, Y. Su, M. C. Tran, N. Wiebe, S. Zhu, "Theory of Trotter error with commutator scaling," *Phys. Rev. X* **11**, 011020 (2021).
4. M. A. Nielsen, I. L. Chuang, *Quantum Computation and Quantum Information* (Cambridge University Press, 2000), Ch. 9.
5. T. M. Cover, J. A. Thomas, *Elements of Information Theory* (Wiley, 2006), Lemma 11.6.1.
6. S. R. White, *Phys. Rev. Lett.* **69**, 2863 (1992); U. Schollwöck, *Ann. Phys.* **326**, 96 (2011).

## Appendix A. Proof of the guarantee

Write $D(\rho,\sigma)=\tfrac12\|\rho-\sigma\|_1$. For pure states $D=\sqrt{1-|\langle a|b\rangle|^2}=\sin\theta_{ab}$, with $\theta_{ab}$ the angle between the rays. $D$ obeys the triangle inequality. For any state $\rho$ and pure $g$, $1-\langle g|\rho|g\rangle\le D(\rho,g)\le\sqrt{1-\langle g|\rho|g\rangle}$; for pure $\rho$ the right side is an equality, so infidelity is $D^2$. Let $\varphi=F(H_s)\psi/\|F(H_s)\psi\|$ be the ideal filtered state. The theorem follows from $D(\tilde\varphi,g)\le D(\tilde\varphi,\varphi)+D(\varphi,g)$ and two lemmas.

**Lemma 1 (leakage).** $D(\varphi,g)^2=\ell\le\bar\ell$.
*Proof.* The ground state sits at $E=0$, so $|\langle g|\varphi\rangle|^2=\gamma F(0)^2/(\gamma F(0)^2+y)$ with $y=\sum_{k\ge1}|c_k|^2F(E_k)^2$. All $E_k$, $k\ge1$, lie in $[\Delta,1]$, so $F(E_k)^2\le\eta^2F(0)^2$ and $y\le(1-\gamma)\eta^2F(0)^2$. Since $\ell=y/(\gamma F(0)^2+y)$ is increasing in $y$, the bound follows. $\square$

**Lemma 2 (Trotter).** $D(\tilde\varphi,\varphi)\le\min\{1,\epsilon_T/\sqrt{p_g}\}$ with $\epsilon_T=\sum_i\alpha t_i^2/(2k_i)=\alpha T^2/(2n)$.
*Proof.* Let $A_i=\langle0_{\mathrm{a}}|\mathcal W_i|0_{\mathrm{a}}\rangle$ and $\tilde A_i$ the Trotterised block; $a=A_m\cdots A_1\psi$, $b=\tilde A_m\cdots\tilde A_1\psi$. The first-order product-formula bound (below) gives $\|\mathcal W_i-\widetilde{\mathcal W}_i\|\le\alpha t_i^2/(2k_i)$; the blocks are compressions, so $\|A_i-\tilde A_i\|$ is no larger, and all blocks have norm $\le1$. Telescoping, $a-b=\sum_i(A_m\cdots A_{i+1})(A_i-\tilde A_i)(\tilde A_{i-1}\cdots\tilde A_1\psi)$, so $\|a-b\|\le\epsilon_T$. The rays of $a,b$ are $\varphi,\tilde\varphi$, so $D=\sin\angle(a,b)=\mathrm{dist}(a,\mathrm{span}\,b)/\|a\|\le\|a-b\|/\|a\|$, and $\|a\|^2=P_{\mathrm{succ}}\ge p_g$. $\square$

**Corollaries.** Infidelity $=D^2\le\bar D^2$. Energy: with $\tilde c_k$ the components of $\tilde\varphi$, $\langle H\rangle-E_0=\sum_{k\ge1}|\tilde c_k|^2(E_k-E_0)\le W(1-|\tilde c_0|^2)=WD^2$.

**First-order Trotter bound.** For Hermitian $A,B$, $U(s)=e^{-is(A+B)}$, $V(s)=e^{-isA}e^{-isB}$: $\frac{d}{ds}[U(\delta-s)V(s)]=iU(\delta-s)(B-e^{-isA}Be^{isA})V(s)$, so $U(\delta)-V(\delta)=-i\int_0^\delta U(\delta-s)(B-e^{-isA}Be^{isA})V(s)\,ds$. Since $\|B-e^{-isA}Be^{isA}\|\le s\|[A,B]\|$ and $U,V$ are unitary, $\|U(\delta)-V(\delta)\|\le\tfrac{\delta^2}2\|[A,B]\|$. For many terms apply this with $A=H_1$, $B=\sum_{g>1}H_g$ and recurse, giving $\frac{\delta^2}{2}\alpha$ per step; with $k$ steps of size $t/k$ the error is $\alpha t^2/(2k)$ [3]. The reversed product order gives the reversed sum, hence the maximum over orders. For $H_s\otimes Z$: $[H_gZ,H_{g'}Z]=[H_g,H_{g'}]\otimes\mathbb 1$, so $\alpha$ is unchanged; the constant shift commutes with everything and is applied exactly.

## Appendix B. Where Pinsker's inequality applies

For the target $|g\rangle\langle g|$ the relative entropy $S(\rho\|\,|g\rangle\langle g|)$ is infinite unless $\rho=|g\rangle\langle g|$, so $D\le\sqrt{S/2}$ gives nothing here. It does apply to the energy-resolved distributions $p_k=|c_k|^2$ (trial) and $q_k=p_kF(E_k)^2/P_{\mathrm{succ}}$ (filtered): $\mathrm{KL}(q\|p)=\sum_kq_k\ln(F(E_k)^2/P_{\mathrm{succ}})\le\ln(1/P_{\mathrm{succ}})$, and $\mathrm{TV}(p,q)\ge q_0-p_0=1-\ell-\gamma$, so with $\mathrm{TV}\le\sqrt{\mathrm{KL}/2}$,
$$P_{\mathrm{succ}}\le\exp\!\big(-2(1-\ell-\gamma)^2\big).$$
This is a necessary condition on the post-selection cost, not an error bound. We checked it on the five $\varepsilon=10^{-2}$ cases with $N\le8$: all inequalities held, but the cap ($0.95$–$1.00$) is far above the actual $P_{\mathrm{succ}}$ ($0.36$–$0.94$), so it is uninformative for these trial states.
