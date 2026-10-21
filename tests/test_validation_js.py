"""Runs the browser-side validation tests (tests/test_validation.js) as part of pytest."""
import shutil
import subprocess
from pathlib import Path

import pytest


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_browser_validation():
    r = subprocess.run(["node", str(Path(__file__).with_name("test_validation.js"))], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
