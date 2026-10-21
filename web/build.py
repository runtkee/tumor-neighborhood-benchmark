"""Build the single-file browser app.
Outputs:
  web/dist/tumor-neighborhood-benchmark.html   standalone: open from disk, works offline; a CSP blocks all network access
  web/dist/artifact.html                       same app as page content for hosting as a private link
  docs/index.html                              copy of the standalone file, served by GitHub Pages
Usage: uv run python web/build.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
WEB, DIST, REF = ROOT / "web", ROOT / "web" / "dist", ROOT / "src" / "bio" / "kit" / "reference"
EXAMPLE_PATIENTS = 6
LABELS = {  # id -> (label, short label for plot columns)
    "IMMUcan_IMC_panel1": ("IMMUcan · IMC · 5 cancers", "IMMUcan IMC"),
    "CODEX_CRC_Schurch2020": ("CODEX · colorectal", "CODEX colon"),
    "SPATCH_COAD_xenium5k": ("SPATCH colon · Xenium", "SPATCH colon-X"),
    "SPATCH_COAD_cosmx6k": ("SPATCH colon · CosMx", "SPATCH colon-C"),
    "SPATCH_OV_xenium5k": ("SPATCH ovary · Xenium", "SPATCH ovary"),
    "SPATCH_HCC_xenium5k": ("SPATCH liver · Xenium", "SPATCH liver"),
}


def reference():
    rows, datasets, recipe = [], [], None
    for f in sorted(REF.glob("*.csv")):
        df = pd.read_csv(f)
        recipe = df["recipe"].iloc[0]
        keep = ["dataset", "patient"] + [c for c in df.columns if c.endswith("_log2") or c.endswith("_diff")]
        df = df[keep]
        did = df["dataset"].iloc[0]
        label, short = LABELS.get(did, (did, did))
        datasets.append({"id": did, "label": label, "short": short, "n": int(df["patient"].nunique())})
        for r in df.to_dict("records"):
            rows.append({k: (None if isinstance(v, float) and np.isnan(v) else (round(v, 6) if isinstance(v, float) else str(v) if k == "patient" else v))
                         for k, v in r.items()})
    order = list(LABELS)
    datasets.sort(key=lambda d: order.index(d["id"]) if d["id"] in order else 99)
    return {"recipe": recipe, "datasets": datasets, "rows": rows}


def example():
    m = yaml.safe_load((ROOT / "kit_manifests" / "CODEX_CRC_Schurch2020.yaml").read_text())
    c = m["columns"]
    need = [c[k] for k in ("patient", "field", "x", "y", "type")] + list(m["states"].values()) + list(m["bleed"].values())
    src = pd.read_csv(ROOT / m["path"], usecols=need, dtype={c[k]: str for k in ("patient", "field", "type")},
                      float_precision="round_trip")
    counts = src.groupby(c["patient"])[c["type"]].apply(lambda s: min((s == "granulocytes").sum(), s.str.contains("T cells").sum()))
    pats = counts.sort_values(ascending=False).index[:EXAMPLE_PATIENTS]
    sub = src[src[c["patient"]].isin(pats)][need]
    cfg = {"dataset": "EXAMPLE_CODEX_subset", "columns": c, "pixel_size_um": m["pixel_size_um"], "roles": m["roles"],
           "states": m["states"], "bleed": m["bleed"]}
    return sub.to_csv(index=False), cfg


def main():
    DIST.mkdir(parents=True, exist_ok=True)
    ref = reference()
    ex_csv, ex_cfg = example()
    css = (WEB / "app.css").read_text()
    body = (WEB / "body.html").read_text()
    blocks = (f'<script type="application/json" id="ref-data">{json.dumps(ref, separators=(",", ":"))}</script>\n'
              f'<script type="application/json" id="example-cfg">{json.dumps(ex_cfg)}</script>\n'
              f'<script type="text/plain" id="example-data">{ex_csv}</script>\n'
              f'<script type="text/plain" id="recipe-src">{(WEB / "recipe.js").read_text()}</script>\n'
              f'<script>{(WEB / "app.js").read_text()}</script>\n')
    for s in [ex_csv, (WEB / "recipe.js").read_text()]:
        assert "</script" not in s.lower()
    title = "<title>Tumor Neighborhood Benchmark</title>"
    fonts = ('<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
             '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans+Condensed:wght@600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">')
    # hosted page: content only (the host adds the document skeleton and its own sandbox)
    (DIST / "artifact.html").write_text(f'{title}\n<meta name="claude:viewport-fit-sides" content="cover">\n{fonts}\n<style>{css}</style>\n{body}\n'
                                        f'<script>window.__HOSTED = true;</script>\n{blocks}')
    # standalone: full document; CSP forbids every network request (fonts fall back to system faces)
    csp = ("default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data: blob:; "
           "worker-src blob:; connect-src 'none'; font-src 'none'; form-action 'none'; base-uri 'none'")
    (DIST / "tumor-neighborhood-benchmark.html").write_text(
        f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
        f'<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">'
        f'<meta http-equiv="Content-Security-Policy" content="{csp}">{title}<style>{css}</style></head><body>{body}\n{blocks}</body></html>')
    # GitHub Pages serves docs/ on main: the standalone file is the public site
    (ROOT / "docs").mkdir(exist_ok=True)
    (ROOT / "docs" / "index.html").write_text((DIST / "tumor-neighborhood-benchmark.html").read_text())
    (ROOT / "docs" / ".nojekyll").write_text("")
    for f in DIST.iterdir():
        print(f"{f.name}: {f.stat().st_size / 1e6:.1f} MB")
    print(f"reference: {sum(d['n'] for d in ref['datasets'])} patients, recipe {ref['recipe']}; example: {ex_csv.count(chr(10)) - 1} cells")


if __name__ == "__main__":
    main()
