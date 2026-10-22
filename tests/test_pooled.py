"""Pooled estimate for the sparse contact tests (recipe v3p.1): bio.pooled in Python, Findings.pooled in the browser file."""
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from bio.pooled import N_BOOT, pooled

ROOT = Path(__file__).resolve().parents[1]
KEY = "nest_neutTouch_PD1"


def fixture_rows(n=60, seed=3):
    """Per-patient rows as the recipe writes them: some patients without near cells, some with NaN columns, one without the columns."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        near = int(rng.poisson(3)) if i % 5 else 0
        far = int(rng.integers(0, 400))
        rows.append({"patient": f"P{i}", f"{KEY}_n": near, f"{KEY}_near_sum": float(near * rng.uniform(35, 60)),
                     f"{KEY}_far_n": far, f"{KEY}_far_sum": float(far * rng.uniform(45, 55)), f"{KEY}_diff": float("nan")})
    rows[7][f"{KEY}_near_sum"] = float("nan")   # a patient with a missing value is skipped
    del rows[9][f"{KEY}_far_n"]                 # so is one without the column
    return rows


def test_pooled_formula_and_counts():
    rows = fixture_rows()
    p = pooled(rows, KEY)
    use = [r for i, r in enumerate(rows) if i not in (7, 9) and (r[f"{KEY}_n"] > 0 or r[f"{KEY}_far_n"] > 0)]
    sn, ss = sum(r[f"{KEY}_n"] for r in use), sum(r[f"{KEY}_near_sum"] for r in use)
    fn, fs = sum(r[f"{KEY}_far_n"] for r in use), sum(r[f"{KEY}_far_sum"] for r in use)
    assert p["diff"] == pytest.approx(ss / sn - fs / fn, abs=1e-12)
    assert (p["n_near"], p["n_far"], p["used"]) == (sn, fn, len(use))
    assert p["patients"] == sum(1 for r in use if r[f"{KEY}_n"] > 0)
    assert p["lo"] < p["diff"] < p["hi"] and p["hi"] - p["lo"] < 15
    assert pooled(pd.DataFrame(rows), KEY) == p      # DataFrame and list of dicts agree
    assert pooled(rows, KEY) == p                    # deterministic (seeded by the key)
    assert pooled(rows, "no_such_test") is None


def test_pooled_bootstrap_is_the_documented_one():
    """200 resamples of patients with replacement, draws from Mulberry32(fnv1a("pooled|" + key)), 2.5 / 97.5 percentiles."""
    from bio.prng import Mulberry32, fnv1a
    rows = [r for i, r in enumerate(fixture_rows()) if i not in (7, 9) and (r[f"{KEY}_n"] > 0 or r[f"{KEY}_far_n"] > 0)]
    a = np.array([[r[f"{KEY}_n"], r[f"{KEY}_near_sum"], r[f"{KEY}_far_n"], r[f"{KEY}_far_sum"]] for r in rows])
    rng, est = Mulberry32(fnv1a("pooled|" + KEY)), []
    for _ in range(N_BOOT):
        s = a[[int(rng.random() * len(a)) for _ in range(len(a))]].sum(axis=0)
        if s[0] > 0 and s[2] > 0:
            est.append(s[1] / s[0] - s[3] / s[2])
    lo, hi = np.percentile(est, [2.5, 97.5])
    p = pooled(rows, KEY)
    assert N_BOOT == 200 and p["lo"] == pytest.approx(lo, abs=1e-9) and p["hi"] == pytest.approx(hi, abs=1e-9)


def test_pooled_empty_arm_gives_nan():
    p = pooled([{"patient": "A", f"{KEY}_n": 4, f"{KEY}_near_sum": 200.0, f"{KEY}_far_n": 0, f"{KEY}_far_sum": 0.0}], KEY)
    assert math.isnan(p["diff"]) and math.isnan(p["lo"]) and p["n_near"] == 4 and p["patients"] == 1


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_python_and_browser_pooled_agree(tmp_path):
    f = tmp_path / "rows.json"
    f.write_text(json.dumps([{k: (None if isinstance(v, float) and math.isnan(v) else v) for k, v in r.items()} for r in fixture_rows()]))
    out = json.loads(subprocess.run(["node", str(ROOT / "tests" / "pooled_parity.js"), str(f), KEY], capture_output=True, text=True, check=True).stdout)
    want = pooled(fixture_rows(), KEY)
    assert set(out) == set(want)
    for k, v in want.items():
        assert out[k] == pytest.approx(v, abs=1e-9), k
