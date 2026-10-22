"""The browser loader must read HALO-style and QuPath-style exports of the same cells and reproduce the results
of the tidy CSV exactly. Fixture: 3 public CODEX patients (Schürch et al. 2020, CC BY 4.0)."""
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

from bio.kit import manifest as mf
from bio.kit.compute import compute
from bio.kit.importers import load_cells

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "codex_subset_3pt.csv"
MAN = yaml.safe_load((ROOT / "kit_manifests" / "CODEX_CRC_Schurch2020.yaml").read_text())
C = MAN["columns"]
ROLE_DYE = {"T": "Dye 1", "neut": "Dye 2", "mac": "Dye 3", "B": "Dye 4", "stroma": "Dye 5", "tumor": "Dye 6", "drop": "Dye 7"}
STATE_COL = {"prolif": "Dye 8 Nucleus Intensity", "PD1": "Dye 9 Membrane Intensity", "GZMB": "Dye 10 Cytoplasm Intensity",
             "LAG3": "Dye 11 Membrane Intensity", "ICOS": "Dye 12 Membrane Intensity"}
BLEED_COL = {"neut": "Dye 2 Membrane Intensity", "mac": "Dye 3 Cytoplasm Intensity", "stroma": "Dye 5 Cytoplasm Intensity"}


def _js(cfg, files, tmp):
    cfg_p, out = tmp / "cfg.json", tmp / "out.json"
    cfg_p.write_text(json.dumps(cfg))
    subprocess.run(["node", str(ROOT / "web" / "node_run.js"), str(cfg_p), str(out)] + [str(f) for f in files], check=True)
    return pd.DataFrame(json.loads(out.read_text())["rows"])


def _python_reference(tmp, subsets=True):
    m = {k: v for k, v in MAN.items() if k != "mac_states"}  # dye-numbered exports carry no macrophage-state columns
    m.update(path=str(FIX), dataset="T", _manifest_path="x")
    if not subsets:  # dye-numbered flags can't name CD8 / CD4 / Treg, so the reference reads every T cell as plain T
        m["roles"] = {k: v.split(":")[0] for k, v in MAN["roles"].items()}
    return compute(load_cells(m), m)


def _same(py, js):
    from tests.analysis_parity import compare
    rep = compare(py, js, 1e-9)
    assert not rep["mismatches"] and not rep["missing_in_js"] and not rep["extra_in_js"], rep


@pytest.fixture(scope="module")
def src():
    return pd.read_csv(FIX, dtype={C[k]: str for k in ("patient", "field", "type")}, float_precision="round_trip")


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_halo_style_export(tmp_path, src):
    roles = {k: v.split(":")[0] for k, v in MAN["roles"].items()}  # dye flags carry no T subset
    files = []
    for field, g in src.groupby(C["field"]):
        pat = g[C["patient"]].iloc[0]
        out = pd.DataFrame({"Image Location": field, "Analysis Region": "TISSUE", "Object Id": range(len(g)),
                            "XMin": g[C["x"]] - 6, "XMax": g[C["x"]] + 6, "YMin": g[C["y"]] - 7, "YMax": g[C["y"]] + 7})
        role = g[C["type"]].map(roles).fillna("other")
        for r, dye in ROLE_DYE.items():
            out[f"{dye} Positive"] = (role == r).astype(int)
        for s, col in STATE_COL.items():
            out[col] = g[MAN["states"][s]].values
        for r, col in BLEED_COL.items():
            out[col] = g[MAN["bleed"][r]].values
        f = tmp_path / f"P{pat}_CRC_{field}.csv"
        f.write_bytes(b"\xef\xbb\xbf" + out.to_csv(index=False, lineterminator="\r\n").encode())  # BOM + CRLF, like HALO
        files.append(f)
    cfg = {"dataset": "T", "pixel_size_um": MAN["pixel_size_um"], "roles": {},
           "columns": {"patient": "__file__", "field": "Image Location", "x": {"min": "XMin", "max": "XMax"},
                       "y": {"min": "YMin", "max": "YMax"}, "type": None},
           "patient_rule": {"mode": "regex", "regex": "^P(\\d+)_"},
           "phenotype_rules": [{"role": r, "positive": [f"{d} Positive"]} for r, d in ROLE_DYE.items()],
           "states": STATE_COL, "bleed": BLEED_COL}
    js = _js(cfg, files, tmp_path)
    py = _python_reference(tmp_path, subsets=False)
    # field names differ (file names vs CODEX field ids) but each field maps 1:1, so per-patient results must match
    _same(py, js)


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_qupath_style_export(tmp_path, src):
    px = MAN["pixel_size_um"]
    out = pd.DataFrame({"Image": "pt" + src[C["patient"]] + " - " + src[C["field"]] + ".ome.tif", "Name": "",
                        "Class": src[C["type"]], "Centroid X µm": src[C["x"]] * px, "Centroid Y µm": src[C["y"]] * px})
    for s, col in MAN["states"].items():
        out[f"Cell: {s} mean"] = src[col].values
    for r, col in MAN["bleed"].items():
        out[f"Cell: {r} marker mean"] = src[col].values
    f = tmp_path / "measurements.tsv"
    out.to_csv(f, sep="\t", index=False, float_format=None)
    cfg = {"dataset": "T", "pixel_size_um": 1.0, "roles": MAN["roles"],
           "columns": {"patient": "Image", "field": "Image", "x": "Centroid X µm", "y": "Centroid Y µm", "type": "Class"},
           "patient_rule": {"mode": "regex", "regex": "^pt(\\d+)"}, "field_rule": {"mode": "regex", "regex": " - (.+)\\.ome\\.tif$"},
           "states": {s: f"Cell: {s} mean" for s in MAN["states"]}, "bleed": {r: f"Cell: {r} marker mean" for r in MAN["bleed"]}}
    js = _js(cfg, [f], tmp_path)
    py = _python_reference(tmp_path)
    _same(py, js)


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_real_halo_sample_parses(tmp_path):
    """SPIAT's tiny HALO export (Artistic-2.0): BOM, μm² headers, box coordinates, Dye flags."""
    cfg = {"dataset": "T", "pixel_size_um": 1.0, "roles": {},
           "columns": {"patient": "Analysis Region", "field": "Analysis Region", "x": {"min": "XMin", "max": "XMax"},
                       "y": {"min": "YMin", "max": "YMax"}, "type": None},
           "phenotype_rules": [{"role": "T", "positive": ["Dye 3 Positive"]}, {"role": "tumor", "positive": ["Dye 2 Positive"]}],
           "states": {}, "bleed": {}}
    cfg_p, out = tmp_path / "c.json", tmp_path / "o.json"
    cfg_p.write_text(json.dumps(cfg))
    subprocess.run(["node", str(ROOT / "web" / "node_run.js"), str(cfg_p), str(out), str(ROOT / "tests/fixtures/tiny_halo.csv")], check=True)
    assert json.loads(out.read_text())["n_cells"] == 10
