"""Kit smoke tests on synthetic data: wrong pixel size must fail; a planted effect must be found."""
import numpy as np
import pandas as pd
import yaml

from bio.kit import checks, manifest
from bio.kit.compute import compute
from bio.kit.importers import load_cells


def _write(tmp_path, px, plant=True):
    rng = np.random.default_rng(1)
    rows = []
    for p in range(6):
        n = 2000
        xy = rng.uniform(0, 600, (n, 2))  # µm; ~8 µm spacing at this density
        t = rng.choice(["Tcell", "Neu", "Mac", "Bcell", "Tum"], n, p=[.2, .1, .1, .1, .5])
        ki = rng.normal(0, 1, n)
        if plant:  # T cells near macrophages get higher Ki67
            mac = xy[t == "Mac"]
            d = np.min(np.linalg.norm(xy[:, None, :] - mac[None, :, :], axis=2), axis=1)
            ki[(t == "Tcell") & (d < 10)] += 2
        mk = rng.normal(0, 1, n); mk[t == "Mac"] += 3
        mn = rng.normal(0, 1, n); mn[t == "Neu"] += 3
        rows.append(pd.DataFrame({"pt": f"P{p}", "img": f"P{p}_1", "cx": xy[:, 0] / px, "cy": xy[:, 1] / px, "lab": t,
                                  "Ki67": ki, "CD68": mk, "MPO": mn}))
    pd.concat(rows).to_csv(tmp_path / "cells.csv", index=False)
    m = {"dataset": "SYN", "modality": "protein", "path": str(tmp_path / "cells.csv"),
         "columns": {"patient": "pt", "field": "img", "x": "cx", "y": "cy", "type": "lab"},
         "pixel_size_um": px, "roles": {"Tcell": "T", "Neu": "neut", "Mac": "mac", "Bcell": "B", "Tum": "tumor"},
         "states": {"prolif": "Ki67"}, "bleed": {"mac": "CD68", "neut": "MPO"}}
    (tmp_path / "m.yaml").write_text(yaml.safe_dump(m))
    return manifest.load(tmp_path / "m.yaml")


def test_wrong_pixel_size_fails(tmp_path):
    m = _write(tmp_path, px=0.5, plant=False)
    m["pixel_size_um"] = 5.0  # 10x too large
    res = checks.run_checks(load_cells(m), m)
    assert any(n.startswith("cell spacing") and s == "fail" for n, s, _ in res)


def test_planted_contact_effect_found(tmp_path):
    m = _write(tmp_path, px=0.5)
    assert manifest.problems(m) == []
    df = load_cells(m)
    assert not checks.blocking(checks.run_checks(df, m))
    rows = compute(df, m)
    assert (rows["macTouch_prolif_diff"] > 10).all()          # planted effect, every patient
    assert rows["neutTouch_prolif_diff"].abs().median() < 10  # no effect where none was planted


def test_cell_spacing_lower_bound(tmp_path):
    m = _write(tmp_path, px=0.5, plant=False)
    df = load_cells(m)
    nn = np.median([np.median(__import__("scipy.spatial", fromlist=["cKDTree"]).cKDTree(g[["x_um", "y_um"]].to_numpy()).query(
        g[["x_um", "y_um"]].to_numpy(), k=2)[0][:, 1]) for _, g in df.groupby("field")])
    for target, want in ((4.5, "pass"), (3.5, "fail")):
        m2 = dict(m, pixel_size_um=0.5 * target / nn)
        res = [r for r in checks.run_checks(load_cells(m2), m2) if r[0].startswith("cell spacing")][0]
        assert res[1] == want and "4–15" in res[2], res


def test_patient_rule_derives_ids_like_the_browser():
    from bio.kit.importers import derive
    rule = {"mode": "regex", "regex": "crc_([a-z]+)_"}
    assert derive("run1_crc_d_c", rule) == "d" and derive("run2_crc_d_m", rule) == "d"
    assert derive("no_match", rule) == "no_match"                      # unmatched values stay whole
    assert derive("P01_core", {"mode": "before", "sep": "_"}) == "P01"
    assert derive("abc", {"mode": "regex", "regex": "b"}) == "b"       # no group -> whole match
    assert derive("abc", {"mode": "regex", "regex": "("}) == "abc"     # bad pattern -> value
    assert derive("abc", None) == "abc"


def test_patient_rule_is_applied_to_the_cell_table(tmp_path):
    m = _write(tmp_path, px=0.5)
    m["patient_rule"] = {"mode": "before", "sep": "_"}   # "P0" .. "P5" have no "_": unchanged
    m["field_rule"] = {"mode": "regex", "regex": "^(P\\d)"}  # "P3_1" -> "P3"
    df = load_cells(m)
    assert set(df["patient"]) == {f"P{i}" for i in range(6)} and set(df["field"]) == {f"P{i}" for i in range(6)}
