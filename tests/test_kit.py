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
