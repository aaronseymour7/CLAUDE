# QPF: the Stetcu-Baroni-Carlson projection filter

The base engine, built to the paper's definitions (main text and Supplemental Material). Our own Hamiltonians and
experiments will be added on top of it later. `SPEC.md` maps every rule of the paper to code, flags the few
choices the paper leaves open, and compares our numbers with the paper's 4 × 4 Heisenberg example.

```
pip install -r requirements.txt
python -m pytest tests -q                 # 7 tests, ~5 s
python examples/heisenberg_4x4.py         # paper's example, ~20 s -> figs/, results/heisenberg_4x4.json
```

| File | Content |
|--|--|
| `qpf/core.py` | the ancilla step exp(−i(tO + δ) ⊗ Y), post-selected runs, filter function, success probability, scaling of H |
| `qpf/symmetry.py` | quantum-number projection: zero phases, halving times from t1 (J² ladder, shifted targets) |
| `qpf/design.py` | energy projection: the supplement's pseudo-spectrum and objective, optimised (t_i) and (t_i, δ_i); constant and Gaussian schedules |
| `qpf/models.py` | the paper's test system: Heisenberg model on a periodic square lattice, J_z sectors, J², Néel state |
| `examples/heisenberg_4x4.py` | J² projection then energy projection, as in the paper's Figs. 1-2 |
