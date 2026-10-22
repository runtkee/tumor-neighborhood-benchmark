"""Headless-Chrome smoke test of the built page (skipped without Chrome): the example renders, a bad file is
rejected with an explanation, and there are no console errors."""
import glob
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _find_chrome():
    if os.environ.get("CHROME"):
        return os.environ["CHROME"]
    for pattern in ("/opt/pw-browsers/chromium-*/chrome-linux/chrome",):  # cloud sandbox (Playwright's Chromium)
        found = sorted(glob.glob(pattern))
        if found:
            return found[-1]
    return "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


CHROME = _find_chrome()


@pytest.mark.skipif(not (Path(CHROME).exists() and shutil.which("node")), reason="needs Chrome and node")
def test_page_end_to_end(tmp_path):
    xlsx = tmp_path / "report.xlsx"
    xlsx.write_bytes(b"PK\x03\x04not really excel")
    r = subprocess.run(["node", str(ROOT / "tests" / "browser_smoke.js"), str(ROOT / "docs" / "index.html"), str(xlsx)],
                       capture_output=True, text=True, timeout=300, env={**os.environ, "CHROME": CHROME})
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "Excel workbook" in r.stdout


@pytest.mark.skipif(not (Path(CHROME).exists() and shutil.which("node")), reason="needs Chrome and node")
def test_kit_results_table_loads_as_dataset():
    """A per-patient table written by the kit is shown as its own dataset, with findings, no pixel-size question and no computation."""
    ref = ROOT / "src" / "bio" / "kit" / "reference" / "CODEX_CRC_Schurch2020.csv"
    r = subprocess.run(["node", str(ROOT / "tests" / "browser_smoke.js"), str(ROOT / "docs" / "index.html"), str(ref)],
                       capture_output=True, text=True, timeout=300, env={**os.environ, "CHROME": CHROME})
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    state = json.loads(next(line for line in r.stdout.splitlines() if line.startswith("FINAL_STATE "))[len("FINAL_STATE "):])
    extra = [d for d in state if "kit results" in d.get("meta", "")]
    assert len(extra) == 1, state
    d = extra[0]
    assert d["name"] == "CODEX_CRC_Schurch2020"
    assert "35 patients" in d["meta"] and "ready" in d["chips"].lower()
    assert d["findings"] >= 1 and not d["px"]
    last = state[-1]
    detail = next(x for x in last["details"] if x.startswith("CODEX_CRC_Schurch2020"))
    assert "How this file was read" in detail and "Checks (4)" in detail and "Share a summary" not in detail
    assert "Runs in your browser. Loaded files are never uploaded." in last["privacy"]
    assert "pixel" not in " ".join(x.get("chips", "") for x in state).lower()
