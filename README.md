# Tumor Neighborhood Benchmark

Compare a tumor-imaging cell table with 218 public patients on a fixed set of neighborhood tests, without uploading anything:

- Do T cells keep away from neutrophils, and do B cells (the control) keep away just as much?
- Do T cells gather near macrophages?
- Are T cells touching neutrophils or macrophages different (dividing, PD-1, granzyme B, LAG-3, ICOS, dying)? Each test reports how much marker signal bleeds between touching cells.

**Open it:** https://runtkee.github.io/tumor-neighborhood-benchmark/

The page opens with a public example loaded. Drop in your own CSV or tab-separated files, check the guessed column mapping, run the checks, and compare.

## Your data stays on your computer
- The page computes everything inside your browser tab. It has no server component.
- Its Content-Security-Policy (first lines of `docs/index.html`) forbids all network connections, so the browser itself blocks any attempt to send data out.
- For offline use, save `docs/index.html` and open it from disk.
- The optional summary export holds per-patient effect sizes with random patient codes. The key linking codes to patients stays with you.

## Input
One row per cell, from one or several files (e.g. one per slide). Supported:
- **Cell type:** a label column, or HALO-style marker positivity flags with ordered rules (e.g. T cell = CD3 positive; neutrophil subset = CD66b and LOX-1 positive).
- **Position:** centroid columns, or bounding boxes (`XMin`/`XMax`).
- **Patient / image IDs:** taken from a column or the file name, either whole, the text before a separator, or a regular-expression match.
- **Formats:** comma, tab or semicolon delimited, with or without a byte-order mark; `NA` / `#N/A` treated as missing.
- **Private-only roles and readouts:** analysed separately and shown next to the closest public role, never pooled.

Tested with HALO-style, QuPath-style and inForm exports (`tests/test_formats.py`).

## Method (locked recipe v2p.1)
- Touching = under 10 µm; nearby = under 30 µm.
- Every proximity test is compared against labels shuffled within each image (20 shuffles), so tissue layout is held fixed.
- Marker levels are within-image percentiles, so brightness differences between images cancel.
- Each patient contributes one value.
- Automatic checks must pass first: cell spacing (catches a wrong pixel size), cells per role, each cell type's own marker being high in it, and number of patients.

Details: [docs/recipe.md](docs/recipe.md), [docs/kit.md](docs/kit.md).

## Reference data
218 patients: IMMUcan imaging mass cytometry (179 patients, 5 cancers), CODEX colorectal cancer (35), SPATCH platform benchmark (3 patients, 4 runs). Only per-patient summary statistics are bundled. Sources and licenses: [NOTICE.md](NOTICE.md).

## Python kit
The same recipe in Python, for RNA data (AnnData `.h5ad`), batch runs, and rebuilding the reference:
```
uv run python -m bio.kit suggest cells.csv --dataset MY_COHORT > my.yaml   # draft a mapping, then review it
uv run python -m bio.kit run my.yaml --out reports/MY_COHORT
```

## Development
```
uv run pytest -q              # includes Python-vs-browser parity and file-format tests (needs Node)
uv run python web/build.py    # rebuild docs/index.html (needs the public datasets locally, see kit_manifests/)
```
The browser recipe (`web/recipe.js`) and the Python recipe (`src/bio/recipe_protein.py`) share a random-number generator and must agree to 1e-9 on every per-patient number. On CODEX (1,311 values) and IMMUcan (5,655 values) they agree to 1e-12.

## Limits
- Protein-imaging recipe only in the browser. RNA data needs the Python kit.
- Browser memory limits apply to multi-gigabyte whole-slide exports. Export only the needed columns.
- These are associations in tissue snapshots. They do not show that one cell type causes a change in another.

Code: MIT. Data: see NOTICE.md.
