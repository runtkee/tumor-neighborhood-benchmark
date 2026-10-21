# Benchmark kit: compare private data against the public reference, locally

**What it does:** runs the same locked recipe on your dataset that was run on 250+ public patients, then writes one HTML page that puts your patients next to the public ones, claim by claim, each with its control. Nothing leaves the machine. There are no network calls, no AI, and no telemetry.

## Use
```
python -m bio.kit suggest my_cells.csv --dataset MY_COHORT > my.yaml   # 1. draft a manifest
#    2. open my.yaml and fix every line marked "# CHECK" (columns, pixel size, label -> role, markers)
python -m bio.kit check my.yaml                                         # 3. automatic checks
python -m bio.kit run my.yaml --out reports/MY_COHORT                   # 4. report
```
Outputs in `--out`:
- `report.html`: open in any browser; works offline.
- `patients_full.csv`: per-patient results with your real patient IDs. **Keep it local.**
- `share_summary.csv`: the same numbers with each patient replaced by a random code, and no cell-level data. Share it only if you choose to.
- `patient_key_KEEP_PRIVATE.csv`: which code is which patient. **Keep it local.**

## Input
One row per cell with these columns: patient, image/field, x, y, cell-type label, plus marker intensities (protein) or an AnnData `.h5ad` with counts (RNA). The manifest maps your labels onto roles: `T, neut, mac, B, stroma, tumor, other, drop`.
- **Private-only roles:** declare them with `extra_roles: {pmn_mdsc: neut}`. The report shows each one beside its nearest public proxy and never pools them.
- **Private-only T states:** add them under `extra_states: [lmp]`, with a matching `states:` entry pointing to the column.

## Checks (a failed check blocks the run unless `--force`, and the report flags forced runs)
- Cell spacing 5–15 µm after applying `pixel_size_um`. This catches pixel/µm mistakes.
- Enough cells per role (≥300).
- Each role's own marker is clearly higher in that role than in T cells. This catches wrong labels or wrong markers; it already caught vimentin as a bad stroma marker in CODEX.
- Which T-state markers are available, and the number of patients.

## Recipe
Locked v2p (`src/bio/recipe_protein.py`, `docs/recipe.md`): 10 µm = touching, 30 µm = nearby, shuffle null within each image, protein levels as within-image percentiles, one value per patient.

## Reference (`src/bio/kit/reference/`)
Per-patient summary tables from IMMUcan (179 patients), CODEX CRC (35) and SPATCH (3 patients, 4 runs). Sources and licenses are in `SOURCES.html`. Rebuild with `python -m bio.kit build-reference kit_manifests/*.yaml`.

## Not yet supported
- HALO object-data export importer (needs a sample file).
- Steinbock/Hyperion folders (convert to CSV or h5ad first).
- An interactive mapping screen in the Python kit (the browser version below has one).

## Browser version (no setup): Tumor Neighborhood Benchmark
- **Standalone file:** `web/dist/tumor-neighborhood-benchmark.html` (6.6 MB). Opens from disk in any browser and works offline. Its Content-Security-Policy blocks every network request, so the browser itself prevents data from leaving.
- **Hosted:** https://runtkee.github.io/tumor-neighborhood-benchmark/ (the same standalone file, served by GitHub Pages from `docs/index.html`). Computation stays local to the viewer's browser.
- Source: `web/recipe.js` (port of recipe v2p.1), `web/app.js`, `web/app.css`, `web/body.html`. Build with `uv run python web/build.py`; it writes `docs/index.html`.
- **Parity:** `analysis/parity_check.py <manifest>` runs Python and JS on the same CSV. CODEX (1,311 values) and IMMUcan (5,655 values) match to ≤1e-12. `tests/test_parity.py` keeps them in sync.
- CSV input and the protein recipe only. RNA data stays with the Python kit.

## Validation and size (browser version)
- `runChecks` in `web/recipe.js` returns `[check, pass|warn|fail, what we found, how to fix]`. Any fail blocks the comparison unless the user forces it, and forced runs are flagged on the results.
- The loader streams files (`importFiles`), stores coordinates and readouts in typed arrays, and interns repeated strings.
- An image is identified by patient plus image name, so names reused across patients ("Core 1") stay separate. The Python kit applies the same rule.
- Readouts with fewer than two distinct values are left out in both versions.
- Tests: `node tests/test_validation.js` (parser chunking, streaming vs whole-file, 13 bad-data cases), `node tests/browser_smoke.js docs/index.html [files]`.
