"""Python kit and browser recipe (web/recipe.js via Node) must give identical per-patient results (recipe v3p.1)."""
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from bio.kit import manifest
from bio.kit.compute import compute
from bio.kit.importers import load_cells
from tests.analysis_parity import WEB, compare
from tests.test_kit import _write
from tests.test_v3p import _write as _write_v3p, _write_segregated

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "codex_subset_3pt.csv"
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def js_rows(m, tmp_path):
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({k: v for k, v in m.items() if not k.startswith("_")}))
    out = tmp_path / "out.json"
    subprocess.run(["node", str(WEB / "node_run.js"), str(cfg), str(out), m["path"]], check=True)
    return pd.DataFrame(json.loads(out.read_text())["rows"])


def assert_same(m, tmp_path, v3p_columns=()):
    py = compute(load_cells(m), m)
    js = js_rows(m, tmp_path)
    rep = compare(py, js, 1e-9)
    assert not rep["mismatches"] and not rep["missing_in_js"] and not rep["extra_in_js"], rep
    for c in v3p_columns:  # the new columns exist and hold numbers, so the comparison above is not vacuous
        assert c in py.columns and py[c].notna().any(), c
    return py, rep


def test_python_and_browser_agree(tmp_path):
    m = _write(tmp_path, px=0.5)
    assert_same(m, tmp_path, ["comp_nest_pct_T", "nest_T_near_mac_log2", "grad_mac_0-5_log2", "nest_macTouch_prolif_diff"])


def test_agree_on_fixture_with_codex_manifest(tmp_path):
    """3 public CODEX patients with the reviewed manifest: T subsets from the labels, tumor compartments from geometry."""
    m = manifest.load(ROOT / "kit_manifests" / "CODEX_CRC_Schurch2020.yaml")
    m.update(path=str(FIX), dataset="T")
    py, rep = assert_same(m, tmp_path, ["comp_nest_pct_CD8", "comp_nest_pct_Treg", "nest_CD8_near_neut_log2", "nest_T_near_mac_log2",
                                        "grad_neut_0-5_log2", "grad_mac_20-30_log2", "nest_macTouch_prolif_n"])
    assert rep["n_values"] > 400  # many columns x 3 patients really were compared (nest contact diffs are NaN here: too few nest T cells)


@pytest.mark.parametrize("variant", ["subsets_macstates_tumorfree", "annotated"])
def test_agree_on_planted_synthetic(tmp_path, variant):
    """Subsets, macrophage states, an image without tumor cells (compartments unknown), or an explicit compartment column."""
    m = _write_v3p(tmp_path, tumor_free=variant != "annotated", annotate=variant == "annotated")
    py, rep = assert_same(m, tmp_path, ["comp_nest_pct_CD8", "nest_Treg_near_mac_log2", "TmacTouch_PDL1_diff", "grad_mac_0-5_prolif",
                                        "nest_macTouch_CD8_prolif_diff", "nest_macTouch_prolif_near_sum", "nest_macTouch_prolif_far_sum",
                                        "nest_macTouch_CD8_prolif_near_sum", "TmacTouch_PDL1_near_sum", "TmacTouch_PDL1_far_sum",
                                        "nest_TmacTouch_PDL1_near_sum", "nest_TmacTouch_PDL1_far_n", "nest_TmacTouch_PDL1_diff"])
    assert rep["n_values"] > 500


def test_agree_on_segregated_nest_with_tumor_free_pool(tmp_path):
    """Planted case from tests/test_v3p.py: the nest null pool without tumor cells gives the same numbers in both languages."""
    m = _write_segregated(tmp_path)
    py, rep = assert_same(m, tmp_path, ["nest_T_near_mac_log2", "nest_B_near_neut_log2", "grad_mac_0-5_log2", "nest_neutTouch_prolif_near_sum"])
    assert abs(py["nest_T_near_mac_log2"]).max() < 0.3
