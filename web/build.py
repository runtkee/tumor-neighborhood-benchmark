"""Build the single-file browser app.
Outputs:
  web/dist/tumor-neighborhood-benchmark.html   standalone: open from disk, works offline; a CSP blocks all network access
  web/dist/artifact.html                       same app as page content for hosting as a private link
  docs/index.html                              copy of the standalone file, served by GitHub Pages
Usage: uv run python web/build.py
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
WEB, DIST, REF = ROOT / "web", ROOT / "web" / "dist", ROOT / "src" / "bio" / "kit" / "reference"
LABELS = {  # id -> (label, short label for plot columns)
    "IMMUcan_IMC_panel1": ("IMMUcan · IMC · 5 cancers", "IMMUcan IMC"),
    "CODEX_CRC_Schurch2020": ("CODEX · colorectal", "CODEX colon"),
    "SPATCH_COAD_xenium5k": ("SPATCH colon · Xenium", "SPATCH colon-X"),
    "SPATCH_COAD_cosmx6k": ("SPATCH colon · CosMx", "SPATCH colon-C"),
    "SPATCH_OV_xenium5k": ("SPATCH ovary · Xenium", "SPATCH ovary"),
    "SPATCH_HCC_xenium5k": ("SPATCH liver · Xenium", "SPATCH liver"),
    "MARTEAU_CRC_IMC_core": ("Marteau colorectal · IMC · tumor core", "Marteau core"),
    "MARTEAU_CRC_IMC_margin": ("Marteau colorectal · IMC · invasive margin", "Marteau margin"),
    "CERVILLA_COLORECTAL_xenium5k": ("Cervilla colorectal · Xenium", "Cervilla colon"),
    "CERVILLA_LUNG_xenium5k": ("Cervilla lung · Xenium", "Cervilla lung"),
    "CERVILLA_BREAST_xenium5k": ("Cervilla breast · Xenium", "Cervilla breast"),
    "CERVILLA_OVARIAN_xenium5k": ("Cervilla ovary · Xenium", "Cervilla ovary"),
    "MARTEAU_XENIUM_core": ("Marteau colorectal · Xenium · tumor core", "Marteau X core"),
    "MARTEAU_XENIUM_margin": ("Marteau colorectal · Xenium · invasive margin", "Marteau X margin"),
    "MARTEAU_XENIUM_normal": ("Marteau colorectal · Xenium · adjacent normal", "Marteau X normal"),
}

# Where each dataset comes from, by id prefix: the sidebar blurb "Public reference: N patients from K published studies ..." is
# generated from this at load (web/app.js refNote). source = study name, detail = one-line description, unit = what n counts,
# part = what this dataset is within the study (Marteau core vs margin).
SOURCE = [
    ("IMMUcan", {"source": "IMMUcan IMC", "detail": "5 cancer types"}),
    ("CODEX", {"source": "Sch\u00fcrch et al. 2020 CODEX", "detail": "colorectal"}),
    ("SPATCH", {"source": "SPATCH benchmark", "detail": ""}),
    ("MARTEAU_CRC_IMC_core", {"source": "Marteau et al. 2025 IMC", "detail": "colorectal", "unit": "region", "part": "tumor-core"}),
    ("MARTEAU_CRC_IMC_margin", {"source": "Marteau et al. 2025 IMC", "detail": "colorectal", "unit": "region", "part": "margin"}),
    ("MARTEAU_XENIUM_core", {"source": "Marteau et al. 2025 Xenium", "detail": "colorectal", "part": "tumor-core"}),
    ("MARTEAU_XENIUM_margin", {"source": "Marteau et al. 2025 Xenium", "detail": "colorectal", "part": "margin"}),
    ("MARTEAU_XENIUM_normal", {"source": "Marteau et al. 2025 Xenium", "detail": "colorectal", "part": "adjacent-normal"}),
    ("CERVILLA", {"source": "Cervilla et al. 2025 Xenium 5K", "detail": ""}),
]


def source_of(did):
    for prefix, info in SOURCE:
        if did.startswith(prefix):
            return info
    return {"source": did, "detail": ""}


# Reference columns the page reads (web/findings.js and web/app.js). Everything else in the CSVs is left out of the page:
#   *_log2, *_diff                  per-patient results; the page finds its charts and sentences by scanning these names
#                                   (T_near_*_log2, nest_*_near_*_log2, grad_<src>_<bin>_log2, *Touch_*_diff, TmacTouch_*_diff)
#   comp_<nest|margin|far>_cells / _pct_<group>   composition table and the "Tumor nest" sentence
#   <test>_n / _near_sum / _far_n / _far_sum      arm sizes and sums behind the pooled estimate, only for the tests the page pools
#                                   (nest_*Touch_*, TmacTouch_*, nest_TmacTouch_*; see pooledLine in app.js, nestContact / macState in findings.js)
# Not read, so dropped: grad_<src>_<bin>_<state> state means, comp_*_immshare_*, per-test _n for tests the page does not pool, recipe.
READ_RE = re.compile(r"_(log2|diff)$|^comp_(nest|margin|far)_(cells|pct_(T|CD8|CD4|Treg|neut|mac|B))$")
POOLED_RE = re.compile(r"^(nest_\w+Touch|TmacTouch|nest_TmacTouch)_")
POOLED_SUFFIX = ("_n", "_near_sum", "_far_n", "_far_sum")


def keep_columns(columns):
    cols = list(columns)
    pooled = {c[:-9] for c in cols if c.endswith("_near_sum") and POOLED_RE.match(c)}
    return ["patient"] + [c for c in cols if c != "patient" and (READ_RE.search(c) or any(c == b + x for b in pooled for x in POOLED_SUFFIX))]


def _cell(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    if isinstance(v, (np.integer, int)):
        return int(v)
    if isinstance(v, (np.floating, float)):
        v = round(float(v), 6)
        return int(v) if v.is_integer() else v  # 1.0 -> 1: same number in JavaScript, two bytes shorter
    return str(v)


def reference():
    """{recipe, datasets: [{id, label, short, n, cols, data}]}: columns once per dataset, one array per patient
    (the page rebuilds the row objects; key names repeated on every row were three quarters of the size)."""
    datasets, recipe = [], None
    for f in sorted(REF.glob("*.csv")):
        df = pd.read_csv(f)
        recipe = df["recipe"].iloc[0]
        did = df["dataset"].iloc[0]
        cols = keep_columns(df.columns)
        label, short = LABELS.get(did, (did, did))
        data = [[str(r[0]) if i == 0 else _cell(r[i]) for i in range(len(cols))] for r in df[cols].itertuples(index=False, name=None)]
        datasets.append({"id": did, "label": label, "short": short, **source_of(did), "n": int(df["patient"].nunique()), "cols": cols, "data": data})
    order = list(LABELS)
    datasets.sort(key=lambda d: order.index(d["id"]) if d["id"] in order else 99)
    return {"recipe": recipe, "datasets": datasets}


def main():
    DIST.mkdir(parents=True, exist_ok=True)
    ref = reference()
    css = (WEB / "app.css").read_text()
    body = (WEB / "body.html").read_text()
    blocks = (f'<script type="application/json" id="ref-data">{json.dumps(ref, separators=(",", ":"))}</script>\n'
              f'<script type="text/plain" id="recipe-src">{(WEB / "recipe.js").read_text()}</script>\n'
              f'<script type="text/plain" id="automap-src">{(WEB / "automap.js").read_text()}</script>\n'
              f'<script>{(WEB / "findings.js").read_text()}</script>\n'
              f'<script>{(WEB / "app.js").read_text()}</script>\n')
    for s in [(WEB / n).read_text() for n in ("recipe.js", "automap.js", "findings.js", "app.js")]:
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
    print(f"reference: {sum(d['n'] for d in ref['datasets'])} patients, {sum(len(d['cols']) for d in ref['datasets'])} columns in total, recipe {ref['recipe']}")


if __name__ == "__main__":
    main()
