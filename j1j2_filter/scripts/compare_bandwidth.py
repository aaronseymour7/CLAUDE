"""compare_bandwidth.py OLD_DATA_DIR -- compare guaranteed designs computed with the Pauli-l1 bandwidth (old data) against the
Weyl bandwidth (results/data); writes results/bandwidth_comparison.json (read by build_results.py)."""
import json, sys, glob, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
old_dir = Path(sys.argv[1])
rn, rt = [], []
for f in sorted(glob.glob(str(ROOT / "results" / "data" / "dmrg_N*_eps*.json"))):
    g = old_dir / os.path.basename(f)
    if not g.exists():
        continue
    new, old = json.load(open(f)), json.load(open(g))
    dn, do = new["designs"].get("guaranteed"), old["designs"].get("guaranteed")
    if not dn or not do or "n" not in dn or "n" not in do:
        continue
    if new["design"]["source"] != "floor" or old["design"]["source"] != "floor":
        continue
    rn.append(dn["n"] / do["n"]); rt.append(dn["T_phys"] / do["T_phys"])
out = dict(count=len(rn), n_ratio_min=min(rn), n_ratio_max=max(rn), n_ratio_median=float(np.median(rn)),
           T_ratio_min=min(rt), T_ratio_max=max(rt))
json.dump(out, open(ROOT / "results" / "bandwidth_comparison.json", "w"), indent=1)
print(out)
