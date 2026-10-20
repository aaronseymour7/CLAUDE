# rbm_netket -- RBM variational energies of the open J1-J2 chain (NetKet)

`run_rbm.py` follows the NetKet tutorial <https://netket.readthedocs.io/en/stable/tutorials/gs-j1j2.html>
(source read from `docs/tutorials/gs-j1j2.ipynb` in the netket repository): custom graph with J1/J2 bond colours, `GraphOperator`,
fixed-`S^z` Spin space, complex RBM-type ansatz (Dense 2N -> log cosh -> sum), `MetropolisExchange`, `Sgd(0.01)` + `VMC_SR(diag_shift=0.01)`.

Differences from the tutorial: open boundaries; `sum J S.S` normalisation (tutorial's Pauli form / 4); `S^z=0` gives E0, `S^z=1` gives the
lowest triplet level; `--marshall 0/1` toggles the Marshall sign rotation of the J1 exchange. An RBM energy is a variational upper bound on the
lowest level of its sector.

Environment: NetKet 3.22.4 fails to import on Python < 3.14 (upstream annotation bug), so use Python 3.14:
`uv venv --python 3.14 nkenv && uv pip install --python nkenv/bin/python netket scipy matplotlib`.

    nkenv/bin/python rbm_netket/run_rbm.py --N 12 --J2 0 0.2411 0.5 --sz 0 1 --iters 600 --ed

Results: `results/rbm_N*_J2*_sz*.json` (final energy +- error, variance, learning curve, exact reference when `--ed`).
Check (N=8, J2=0.2411, 300 iterations): E0 = -3.13451(5) vs exact -3.135001 (rel. err 1.6e-4); S=1 min = -2.74881(2) vs -2.749123 (1.1e-4).
