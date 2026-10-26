# sbpipe: certified Stetcu–Baroni filter pipeline

Code and data for the paper in `sbpipe_study/paper/`.

- `sbpipe/`: the pipeline package (energy estimates + ansatz → certified filter plan, resources); see `sbpipe/README.md`.
- `sbpipe_study/`: DMRG / J1–J2 application, all data, figures and the paper; see `sbpipe_study/README.md`.

To reproduce:

    pip install -r requirements.txt
    python -m pytest sbpipe/tests -q            # 10 tests
    python sbpipe_study/example_j1j2.py         # one case end to end
    python sbpipe_study/make_figures.py         # figures and tables from sbpipe_study/results/*.json
    cd sbpipe_study/paper && pdflatex paper && pdflatex paper
