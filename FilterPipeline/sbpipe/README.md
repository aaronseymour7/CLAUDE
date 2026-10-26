# sbpipe: a modular, certified Stetcu–Baroni filter pipeline

`sbpipe` takes **classical energy estimates** and a **trial state (ansatz)** for a Hamiltonian written as a
sum of Pauli strings, and returns a certified plan for refining the trial state with the
Stetcu–Baroni–Carlson cosine filter (one ancilla, controlled time evolution, post-selection), together
with a simple resource count. For small systems it also simulates the planned circuit and checks the
certificate.

```
EnergyEstimates (E0 ± dE0, E1 - dE1, E_max)      Ansatz (gamma, statevector?, circuit cost)
                 \                                  /
                  refine(H, estimates, ansatz, Config)            H : PauliSum
                    1. window   H_s = (H - shift)/W, pass band [0,w0], stop band [delta,1]
                    2. design   F(E) = prod cos(E t_i + phi_i), certified eta   (uniform | dyadic | optimized | auto)
                    3. Trotter  first/second-order commutator bounds -> per-pulse steps k_i
                    4. certificate  1 - F <= (sqrt(ell) + eps_U / sqrt(p_g))^2 <= eps
                    5. resources  qubits, CNOT, rotations, T count, expected cost / success
                  -> Plan  --validate()-->  simulated infidelity, success probability, bound check
```

Nothing in the pipeline is specific to spin chains or to DMRG. A front end only has to fill the two
input containers; `sbpipe/frontends/dmrg_quimb.py` is the one used in the J1–J2 study.

## Quick start

```python
import sbpipe
from sbpipe.frontends.dmrg_quimb import dmrg_inputs

H = sbpipe.j1j2_chain(12, J1=1.0, J2=0.2411)          # any PauliSum works
est, ansatz, _ = dmrg_inputs(H, chi_trial=2)          # or build EnergyEstimates / Ansatz yourself
plan = sbpipe.refine(H, est, ansatz, sbpipe.Config(eps=1e-2, schedule="optimized", order=2))
print(plan.summary())                                  # design, certificate, resources
print(sbpipe.validate(plan, H, ansatz))                # N <= ~14: statevector check
```

Bring your own inputs:

```python
est = sbpipe.EnergyEstimates(E0=-4.5456, E0_err=1e-6, E1=-4.2500, E1_err=1e-3)  # E_max: Weyl bound if omitted
ans = sbpipe.Ansatz(gamma=0.9, cost={"cnot": 120, "rotations": 300}, statevector=None)
```

## Modules

| module | content |
|--|--|
| `pauli.py` | `PauliSum`, Pauli products, local exact operator norms, statevector action |
| `models.py` | `j1j2_chain`, rigorous Weyl spectral bounds for any grouped `PauliSum` |
| `inputs.py` | `EnergyEstimates`, `Ansatz` |
| `design.py` | filter schedules: `uniform` and `dyadic` (closed-form certificates), `optimized` (SLSQP design + interval certificate from `floor.py`); `Config(schedule="auto")` keeps the cheaper of dyadic and optimized |
| `trotter.py` | exact local commutator constants α₁ (Childs et al. Prop. 15) and β₂ (Prop. 16), product formulas |
| `pipeline.py` | `Config`, `refine`, `validate`, `certificate` |
| `resources.py` | gate-count model (`GateModel`) |
| `simulate.py` | exact and Trotterized post-selected filter on a statevector |
| `frontends/dmrg_quimb.py` | DMRG E0, E1 (symmetry-sector penalty), reference MPS, bond-dimension-χ trial |

If the trial state already satisfies 1 − γ ≤ ε, `refine` also considers the empty filter and keeps it when
it is cheaper. Synthesis error of the trial circuit's rotations is included in the certificate.

## What is certified and what is assumed

The certificate is a theorem **given the inputs**: E0 lies in [E0 − dE0, E0 + dE0], every other level the
trial state touches is ≥ E1 − dE1, the spectrum is ≤ E_max (automatic with the Weyl bound) and the trial
state's ground-space weight is ≥ gamma. A DMRG front end makes all of these plausible but does not prove
them (DMRG E1 is an upper bound, and gamma is an overlap with a DMRG reference). The study in
`../sbpipe_study/` measures what happens when they are wrong.

## Tests

```
python -m pytest sbpipe/tests -q
```

Checks: Pauli algebra against Kronecker products; Weyl bound against exact spectra; α₁ and β₂ against the
exact Trotter error of dense matrices; dyadic block against the Dirichlet-kernel identity; filter
certificates against a 2×10⁶-point grid; Strang gate counts against an explicit merged sequence; CNOT
counts against qiskit; an end-to-end certificate check.

Performance note: `sbpipe` limits BLAS to one thread on import (`SBPIPE_BLAS_THREADS` overrides). The
many tiny eigenproblems in `trotter.py` and the SLSQP designer are 10–50× slower with threaded BLAS.
