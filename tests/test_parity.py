"""Python kit and browser recipe (web/recipe.js via Node) must give identical per-patient results."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from tests.test_kit import _write
from bio.kit.compute import compute
from bio.kit.importers import load_cells


from tests.analysis_parity import WEB, compare


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_python_and_browser_agree(tmp_path):
    m = _write(tmp_path, px=0.5)
    py = compute(load_cells(m), m)
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({k: v for k, v in m.items() if not k.startswith("_")}))
    out = tmp_path / "out.json"
    subprocess.run(["node", str(WEB / "node_run.js"), str(cfg), str(out), m["path"]], check=True)
    js = pd.DataFrame(json.loads(out.read_text())["rows"])
    rep = compare(py, js, 1e-9)
    assert not rep["mismatches"] and not rep["missing_in_js"] and not rep["extra_in_js"], rep
