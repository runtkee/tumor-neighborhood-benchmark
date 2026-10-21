"""Automatic mapping (web/automap.js, run through web/node_auto.js) must read common export styles of the same
cells and reproduce the hand-reviewed Python results exactly, and must refuse files it can't interpret."""
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

from bio.kit.compute import compute
from bio.kit.importers import load_cells
from tests.analysis_parity import compare

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "codex_subset_3pt.csv"
MAN = yaml.safe_load((ROOT / "kit_manifests" / "CODEX_CRC_Schurch2020.yaml").read_text())
C = MAN["columns"]
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


def auto(tmp, *files):
    out = tmp / "auto.json"
    subprocess.run(["node", str(ROOT / "web" / "node_auto.js"), str(out)] + [str(f) for f in files], check=True, capture_output=True)
    return json.loads(out.read_text())


@pytest.fixture(scope="module")
def reference():
    m = dict(MAN, path=str(FIX), dataset="T", _manifest_path="x")
    return compute(load_cells(m), m)


@pytest.fixture(scope="module")
def src():
    return pd.read_csv(FIX, dtype={C[k]: str for k in ("patient", "field", "type")}, float_precision="round_trip")


def same(py, out):
    js = pd.DataFrame(out["rows"])
    js["patient"] = js["patient"].astype(str)
    rep = compare(py, js, 1e-9)
    assert not rep["mismatches"] and not rep["missing_in_js"] and not rep["extra_in_js"], rep


def test_codex_export_maps_exactly(tmp_path, reference):
    out = auto(tmp_path, FIX)
    assert out["config"]["pixel_size_um"] == 0.377  # from the CODEX column-naming fingerprint
    assert out["config"]["roles"] == {k: v for k, v in MAN["roles"].items() if k in out["config"]["roles"]}
    same(reference, out)


def test_sorted_file_still_finds_cell_types(tmp_path, reference, src):
    f = tmp_path / "sorted.csv"
    src.sort_values(C["type"]).to_csv(f, index=False)  # start of file holds a single cell type
    out = auto(tmp_path, f)
    assert out["config"]["columns"]["type"] == C["type"]


def test_halo_style_marker_flags(tmp_path, reference, src):
    roles = src[C["type"]].map(MAN["roles"]).fillna("other")
    marker = {"T": "CD3", "neut": "CD66b", "mac": "CD68", "B": "CD20", "stroma": "aSMA", "tumor": "PanCK"}
    files = []
    for field, g in src.groupby(C["field"]):
        r = roles.loc[g.index]
        keep = r != "drop"
        g, r = g[keep], r[keep]
        out = pd.DataFrame({"Image Location": field, "Object Id": range(len(g)),
                            "XMin": g[C["x"]] - 6, "XMax": g[C["x"]] + 6, "YMin": g[C["y"]] - 7, "YMax": g[C["y"]] + 7})
        for role, mk in marker.items():
            out[f"{mk} Positive Classification"] = (r == role).astype(int).values
        names = {"prolif": "Ki67 Nucleus Intensity", "PD1": "PD-1 Membrane Intensity", "GZMB": "Granzyme B Cytoplasm Intensity",
                 "LAG3": "LAG-3 Membrane Intensity", "ICOS": "ICOS Membrane Intensity"}
        for s, col in names.items():
            out[col] = g[MAN["states"][s]].values
        out["CD15 Membrane Intensity"] = g[MAN["bleed"]["neut"]].values
        out["CD68 Cytoplasm Intensity"] = g[MAN["bleed"]["mac"]].values
        out["aSMA Cytoplasm Intensity"] = g[MAN["bleed"]["stroma"]].values
        pat = g[C["patient"]].iloc[0]
        fp = tmp_path / f"P{pat}_{field}.csv"
        fp.write_bytes(b"\xef\xbb\xbf" + out.to_csv(index=False, lineterminator="\r\n").encode())
        files.append(fp)
    res = auto(tmp_path, *files)
    assert res["config"]["phenotype_rules"], res.get("report")
    assert res["config"]["columns"]["field"] == "Image Location"
    # pixel size is not stated in this export, so it is estimated: compare with the reference recomputed at that size
    m = dict(MAN, path=str(FIX), dataset="T", _manifest_path="x", pixel_size_um=res["config"]["pixel_size_um"])
    df = load_cells(m)
    df["patient"] = "P" + df["patient"].astype(str)  # same IDs as the file names, so shuffle seeds match
    same(compute(df, m), res)


def test_qupath_style_microns(tmp_path, reference, src):
    px = MAN["pixel_size_um"]
    out = pd.DataFrame({"Image": "pt" + src[C["patient"]] + "_" + src[C["field"]] + ".ome.tif", "Class": src[C["type"]],
                        "Centroid X µm": src[C["x"]] * px, "Centroid Y µm": src[C["y"]] * px})
    for s, col in MAN["states"].items():
        out[{"prolif": "Nucleus: Ki67: Mean", "PD1": "Membrane: PD-1: Mean", "GZMB": "Cytoplasm: Granzyme B: Mean",
             "LAG3": "Membrane: LAG-3: Mean", "ICOS": "Membrane: ICOS: Mean"}[s]] = src[col].values
    out["Membrane: CD15: Mean"] = src[MAN["bleed"]["neut"]].values
    out["Cell: CD68: Mean"] = src[MAN["bleed"]["mac"]].values
    out["Cell: aSMA: Mean"] = src[MAN["bleed"]["stroma"]].values
    f = tmp_path / "measurements.tsv"
    out.to_csv(f, sep="\t", index=False)
    res = auto(tmp_path, f)
    assert res["config"]["pixel_size_um"] == 1 and not res["fails"], res.get("fails")
    assert res["config"]["patient_rule"] == {"mode": "before", "sep": "_"}
    ref = reference.copy()
    ref["patient"] = "pt" + ref["patient"].astype(str)
    js = pd.DataFrame(res["rows"])
    # image names differ from CODEX field names, so shuffle seeds differ: compare only the seed-independent contact metrics
    cols = ["patient"] + [c for c in ref.columns if c.endswith("_diff")]
    rep = compare(ref[cols], js[[c for c in cols if c in js]], 1e-9)
    assert not rep["mismatches"], rep


def test_dye_numbered_halo_fails_with_reason(tmp_path):
    res = auto(tmp_path, ROOT / "tests" / "fixtures" / "tiny_halo.csv")
    assert any("antibody" in f["why"] for f in res["fails"])


def test_excel_file_is_refused(tmp_path):
    f = tmp_path / "x.xlsx"
    f.write_bytes(b"PK\x03\x04garbage")
    res = auto(tmp_path, f)
    assert "Excel" in res["fails"][0]["title"]
