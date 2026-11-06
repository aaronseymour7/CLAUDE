"""Every knob of the study in one place.

    cfg = StudyConfig()                                  # defaults
    cfg = StudyConfig(Ns=(8, 10), eps=(1e-2,))           # override in code
    cfg = StudyConfig.load("my.json")                    # or from a JSON file with any subset of the fields
    cfg = cfg.updated(J2=0.3, gaps=[0.05, 0.5])          # copy with changes

Scripts accept `--config file.json` and `--set key=value` (value parsed as JSON), e.g.
    python run_maps.py --set Ns=[8,10] --set eps=[0.01]
"""
import json
from dataclasses import asdict, dataclass, field, fields, replace

import numpy as np


@dataclass(frozen=True)
class StudyConfig:
    # model: open J1-J2 chain, gap tuned by dimerisation, overlap by a staggered-Rz layer on the ground state
    J1: float = 1.0
    J2: float = 0.2411
    Ns: tuple = (8, 10, 12, 14, 16)
    gaps: tuple = tuple(np.round(np.geomspace(0.005, 1.7, 14), 10))     # physical gap E1 - E0 (units of J1)
    gammas: tuple = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)
    eps: tuple = (1e-2, 1e-3)                                             # target infidelities 1 - F
    delta_range: tuple = (-0.9, 0.95)                                     # dimerisation search interval
    trial: str = "exact"              # trial base: "exact" ground state (study) or "mps1" (one-layer MPS-to-circuit, application)

    # measured search over the dyadic filter: K ladder rungs, r repetitions, n Strang steps per tau
    n_grid: tuple = (1, 2, 3, 4, 6, 8, 12, 16)
    K_max: int = 20
    r_max: int = 64
    protocols: tuple = ("rus", "aa")

    # gate-cost model
    bond_cnots: int = 6               # one ancilla-controlled Heisenberg bond block (Bell change + parity network)
    mcz_cnots_per_qubit: int = 12     # n-qubit controlled Z ~ 2n Toffolis x 6 CNOTs (AA reflections)

    # MPS mode (ornl.mps, run_mps.py; needs quimb): DMRG reference, filter on the system MPS
    mps_Ns: tuple = (16, 20, 24, 32)
    mps_reference: str = "auto"       # "ed", "dmrg", or "auto" = ED up to mps_ed_max_N, DMRG above
    mps_ed_max_N: int = 16
    mps_max_bond: int = 128
    mps_cutoff: float = 1e-12
    mps_dmrg_bonds: tuple = (16, 32, 64, 128)

    # output
    out_dir: str = "results"
    fig_dir: str = "figs"

    def updated(self, **kw):
        kw = {k: tuple(v) if isinstance(v, list) else v for k, v in kw.items()}
        return replace(self, **kw)

    @classmethod
    def load(cls, path):
        return cls().updated(**json.load(open(path)))

    def save(self, path):
        json.dump(asdict(self), open(path, "w"), indent=1)

    @classmethod
    def from_cli(cls, config=None, sets=()):
        """Build from an optional JSON file and a list of 'key=value' strings (value in JSON)."""
        cfg = cls.load(config) if config else cls()
        names = {f.name for f in fields(cls)}
        kw = {}
        for s in sets or ():
            k, v = s.split("=", 1)
            if k not in names:
                raise SystemExit(f"unknown setting {k!r}; choose from {sorted(names)}")
            kw[k] = json.loads(v)
        return cfg.updated(**kw)


def add_cli(ap):
    """Add --config / --set to an argparse parser."""
    ap.add_argument("--config", help="JSON file with any subset of StudyConfig fields")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="override one field, value in JSON, e.g. --set Ns=[8,10]")
    return ap
