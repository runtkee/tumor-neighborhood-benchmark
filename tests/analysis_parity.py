"""Run the Python kit and the browser recipe (via Node) on the same CSV + manifest; compare every per-patient number.
Usage: uv run python analysis/parity_check.py <manifest.yaml> [--tol 1e-9]"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from bio.kit import manifest as mf
from bio.kit.compute import compute
from bio.kit.importers import load_cells

WEB = Path(__file__).resolve().parents[1] / "web"


def compare(py, js, tol):
    py = py.set_index("patient").sort_index()
    js = js.set_index("patient").sort_index()
    report = {"patients_py": len(py), "patients_js": len(js), "missing_in_js": sorted(set(py.index) - set(js.index)),
              "extra_in_js": sorted(set(js.index) - set(py.index))}
    cols = sorted(set(py.columns) | set(js.columns) - {"dataset", "recipe"})
    cols = [c for c in cols if c not in ("dataset", "recipe")]
    worst, mism = 0.0, []
    for c in cols:
        if c not in py or c not in js:
            mism.append((c, "column missing in " + ("js" if c not in js else "py")))
            continue
        a = pd.to_numeric(py[c], errors="coerce").reindex(js.index).to_numpy(float)
        b = pd.to_numeric(js[c], errors="coerce").to_numpy(float)
        nan_mismatch = np.isnan(a) != np.isnan(b)
        both = ~np.isnan(a) & ~np.isnan(b)
        d = np.abs(a[both] - b[both])
        if nan_mismatch.any():
            mism.append((c, f"{nan_mismatch.sum()} NaN mismatches"))
        if d.size and d.max() > tol:
            mism.append((c, f"max diff {d.max():.3g} in {(d > tol).sum()} patients"))
        if d.size:
            worst = max(worst, d.max())
    report.update(n_values=int(sum(pd.to_numeric(py[c], errors="coerce").notna().sum() for c in cols if c in py)),
                  worst_abs_diff=worst, mismatches=mism)
    return report


def main():
    path, tol = sys.argv[1], float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-9
    m = mf.load(path)
    py = compute(load_cells(m), m)
    cfg = {k: v for k, v in m.items() if not k.startswith("_")}
    with tempfile.TemporaryDirectory() as tmp:
        cfg_p, out_p = Path(tmp) / "cfg.json", Path(tmp) / "out.json"
        cfg_p.write_text(json.dumps(cfg))
        subprocess.run(["node", str(WEB / "node_run.js"), str(cfg_p), str(out_p), m["path"]], check=True)
        js = pd.DataFrame(json.loads(out_p.read_text(), parse_constant=lambda x: float("nan"))["rows"])
    rep = compare(py, js, tol)
    print(json.dumps(rep, indent=1, default=str))
    sys.exit(1 if rep["mismatches"] or rep["missing_in_js"] or rep["extra_in_js"] else 0)


if __name__ == "__main__":
    main()
