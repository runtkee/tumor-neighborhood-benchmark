# Tumor Neighborhood Benchmark

Compare a tumor-imaging cell table with 218 public patients on a fixed set of neighborhood tests, without uploading anything:

- Do T cells keep away from neutrophils, and do B cells (the control) keep away just as much?
- Do T cells gather near macrophages?
- Are T cells touching neutrophils or macrophages different (dividing, PD-1, granzyme B, LAG-3, ICOS, dying)? Each test reports how much marker signal bleeds between touching cells.

**Open it:** https://runtkee.github.io/tumor-neighborhood-benchmark/

The page opens with the public datasets in the left panel, each with its summary findings. Add your own data with **Add data**: either a cell table (one or several files), which is read automatically, checked, analysed and listed next to them, or a per-patient results table written by the Python kit (`share_summary.csv`, `patients_full.csv`, or a `src/bio/kit/reference/*.csv` file), which is shown as it is, without computing anything. A results table is recognised from its header (`dataset`, `patient`, `recipe` and at least one `_log2` or `_diff` column) and is listed with a "kit results" tag, one entry per `dataset` value. This is how RNA cohorts (Xenium, CosMx) run in the kit appear in the page. "How this file was read" lists the columns used and the recipe-version check; results from another recipe version are shown but flagged as possibly not comparable. Tick datasets to compare them side by side.

## Loaded data stays on your computer
- The page computes everything inside your browser tab. It has no server component.
- Its Content-Security-Policy (first lines of `docs/index.html`) forbids all network connections, so the browser itself blocks any attempt to send data out.
- For offline use, save `docs/index.html` and open it from disk.
- The optional summary export holds per-patient effect sizes with random patient codes. The key linking codes to patients stays with you.

## How files are read
There are no settings, apart from the pixel size question below. The app decides how to read a file from its column names and from sample values taken across the whole file. Files are often sorted, so it doesn't rely on the first rows. "How this file was read" shows every decision, its evidence and a confidence level.

| Decision | How |
|---|---|
| Cell position | Numeric columns named like coordinates (X, Centroid X µm, Pos_X, Cell X Position), or bounding boxes (XMin/XMax) |
| Units | Stated in the column name (µm); recognised from the export format (CODEX processor naming, inForm pixels); otherwise estimated from cell spacing |
| Patient, image | Columns named like patient / donor / case and image / core / FOV; checked that every image belongs to one patient; otherwise taken from file or image names |
| Cell types | The label column whose values are recognised cell types ("CD8+ T cells", "CD68+CD163+ macrophages", "Neutrophil", "MacCD163") **and** whose labelled cells actually carry that type's marker. Region columns such as "neighborhood name" are ruled out this way. If there are no labels, yes/no columns named after lineage markers (CD3, CD68, CD66b, …) are turned into rules |
| T-cell subsets | Read from the cell-type label: Treg / FOXP3 → Treg, CD8 → CD8, CD4 → CD4. A label naming both CD4 and CD8, or neither, stays plain T. With marker yes/no columns, rules are applied in the order Treg, CD8, CD4, then plain T |
| Readouts | Columns named after Ki-67, PD-1, granzyme B, LAG-3, ICOS and cleaved PARP; mean intensity preferred |
| Macrophage markers | Columns named after PD-L1, IDO1, CD163 or HLA-DR (mean intensity preferred); compared between macrophages touching T cells and macrophages far from them |
| Tumor compartments | Computed from cell positions: nest = tumor cells are at least half of the cells within 50 µm; margin = within 100 µm of a nest; far = the rest. No tumor label means whole-image results only. Each decision appears in "How this file was read" |

**What it can't read:**
- **HALO exports that name markers only by dye** ("Dye 3 Positive"). The file never says which antibody Dye 3 is, so the app explains this and asks for an export with marker names or a phenotype column.
- **Tables with intensities but no cell types.** The app doesn't assign cell types from raw intensities.

**Pixel size is the one thing the app may ask for.** When a file doesn't state units and the format isn't recognised, the app asks how many µm one pixel is, because the 10 and 30 µm thresholds depend on it. It pre-fills an estimate from cell spacing and offers matching instrument values. If the user keeps the estimate, the analysis is also rerun at −20% and +25%, and any finding that changes is flagged. Entered values still go through the cell-spacing check.

Tested: automatic reading of the raw CODEX (35 patients) and IMMUcan (179 patients) files reproduces the hand-reviewed reference exactly. HALO-style, QuPath-style and inForm exports, sorted files, and refusals are covered in `tests/test_automap.py`.

## When data doesn't fit
Nothing is compared until the data passes its checks. Every problem names the file and says how to fix it.
- **On load:** Excel, zip, gzip, PDF and HDF5 (`.h5ad`) files are recognised and refused with instructions. The page also catches empty files, files without a usable delimiter, repeated column names, and files whose columns differ.
- **Before comparing:** the page blocks the comparison when it finds any of these:
  - rows with missing columns;
  - unusable coordinates;
  - empty patient or image IDs;
  - duplicated cells (the same file loaded twice);
  - wrong pixel size (cell spacing outside 4–15 µm, with the implied correct value);
  - implausible image sizes;
  - no T cells;
  - a marker rule built on an intensity column instead of a yes/no flag;
  - a cell type whose own marker isn't high in it (a wrong label or a wrong marker).
- **Notes, which don't block:** constant or text-only readouts (left out), yes/no flags used as readouts, small images skipped, unmapped labels, and patients with too few T cells.

The tests in `tests/test_validation.js` cover each case, and `tests/browser_smoke.js` runs the page in headless Chrome.

## How big a dataset fits
Files are read in chunks, and only the mapped columns are kept. Measured on a 10-core laptop with 26 GB RAM:

| Cells | File size | Time | Peak memory |
|---|---|---|---|
| 50,000 (example) | 6 MB | under 1 s | small |
| 2 million | 223 MB | about 15 s | 1.3 GB |
| 6 million | 671 MB | about 35 s | 2.9 GB |

Expect a typical 8 GB laptop to handle about 2–3 million cells and a 16 GB laptop about 6–8 million. Above that, run groups of patients separately. Results are per patient, so the groups combine exactly. The page warns when a load looks too large.

The Python kit's RNA path reads whole AnnData files into memory. The SPATCH files used for the reference (1–3 GB each) need a 16 GB+ machine.

## Method (locked recipe v3p.1)
- Touching = under 10 µm; nearby = under 30 µm.
- Every proximity test is compared against labels shuffled within each image (20 shuffles), so tissue layout is held fixed.
- Marker levels are within-image percentiles, so brightness differences between images cancel.
- Each patient contributes one value.
- v3p adds tumor compartments (nest, margin, far), the immune composition of each, proximity tests inside the nest (labels shuffled within each compartment), T-cell subsets (CD8, CD4, Treg), distance gradients (0–5 to 100+ µm), and macrophage markers next to T cells. The v2p.1 numbers are unchanged.
- v3p.1 changes two things inside the nest. The shuffled labels are drawn from non-tumor cells only, so immune and stromal cells that share the non-epithelial strands of a nest do not read as enrichment. And the sparse contact tests (nest T cells touching neutrophils or macrophages, macrophage markers next to T cells) also report per-patient sums, so a pooled estimate with a patient-bootstrap 95% interval can be shown when too few patients have 20 cells per arm. A pooled estimate is not a per-patient result, and the page says so. `TmacTouch_*` now covers all macrophages in the image; `nest_TmacTouch_*` is the nest-only version.
- The public reference was built with v3p. The tumor-free null pool and the pooled lines appear once it is rebuilt with v3p.1; the page shows pooled lines only where the sums exist.
- Automatic checks must pass first: cell spacing (catches a wrong pixel size), cells per role, each cell type's own marker being high in it, and number of patients.

Details: [docs/recipe.md](docs/recipe.md), [docs/kit.md](docs/kit.md).

## Reference data: where it comes from
218 patients: IMMUcan imaging mass cytometry (179 patients, 5 cancers), CODEX colorectal cancer (35), SPATCH platform benchmark (3 patients, 4 runs).
- **Not in this repository:** the raw public datasets (about 12 GB). They were downloaded from their public archives (Zenodo, Mendeley Data, the SPATCH portal).
- **In this repository:** they were run once through the Python kit, and only the per-patient summary numbers are committed (`src/bio/kit/reference/*.csv`, about 100 KB) and embedded in the page.
- **Cell-level data in git:** just a 3-patient CODEX subset (`tests/fixtures/`, CC BY 4.0) and two 10-cell SPIAT examples (Artistic-2.0), used by the tests.

Sources and licenses: [NOTICE.md](NOTICE.md).

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
- Very large cohorts (more than about 8 million cells) may exceed browser memory; split them by patient.
- These are associations in tissue snapshots. They do not show that one cell type causes a change in another.
- Research use only; not for diagnosis or treatment decisions.

Code: MIT. Data: see NOTICE.md.
