"""Headless-Chrome smoke test of the built page (skipped without Chrome): the example renders, a bad file is
rejected with an explanation, and there are no console errors."""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CHROME = os.environ.get("CHROME", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")


@pytest.mark.skipif(not (Path(CHROME).exists() and shutil.which("node")), reason="needs Chrome and node")
def test_page_end_to_end(tmp_path):
    xlsx = tmp_path / "report.xlsx"
    xlsx.write_bytes(b"PK\x03\x04not really excel")
    r = subprocess.run(["node", str(ROOT / "tests" / "browser_smoke.js"), str(ROOT / "docs" / "index.html"), str(xlsx)],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
    assert "Excel workbook" in r.stdout
