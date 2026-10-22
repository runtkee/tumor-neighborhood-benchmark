# Locked recipe (`src/bio/recipe.py`)

Fixed parameters: touching <10 µm, nearby <30 µm, activated neutrophil = top 10% of the activation ratio, ≥3 identity counts to trust a ratio, ≥300 cells per type to run a test, 50 shuffles.

Metrics per sample:
- **`T_near_X_log2`**: share of T cells within 30 µm of cell type X, compared with X's labels shuffled among non-T cells. Below 0 means T cells avoid X.
- **`XTouch_<state>_pct`**: % change in T-cell state (state genes ÷ T-identity genes) for T cells under 10 µm from X vs those ≥30 µm away. Adjusted for region, local T density, CD8 vs CD4, and RNA total.
- **`actN_T30_*`**: T cells within 30 µm per SD of neutrophil activation, raw and adjusted (negative binomial model). Log2 scale.
- **Controls:**
  - `CTRL_B_near_neut`: do B/plasma cells avoid neutrophils too? If yes, a T-cell result isn't specific.
  - `CTRL_fibTouch`: does touching a fibroblast give the same "effect"?

Versions:
- **v1:** scored T state per 1000 counts. Neighbor spillover diluted the T cells' own genes, which faked contact effects (the fibroblast control matched the macrophage result).
- **v2 (current):** T state scored as a log ratio to T-identity genes. Frozen 2026-10-20.

Output: `results/recipe_v2.csv`, one row per sample.

Cell labels: the recipe reads CD4T, CD8T, Tprolif, Treg (T cells), Neutrophil, Macrophage, B, Plasma, Fibroblast and Epithelial. `Treg` was added for the kit's Xenium labels (SPATCH has no Treg, so SPATCH rows are unchanged); the method is the same, so the version stays v2.

Running it on Xenium data: a runner script (research-only, not part of the kit) runs the recipe on every section of the h5ad files that `python -m bio.kit xenium` writes (region term dropped when the file has no `spatial_cluster`), together with the Project 2 senescence test. It also adds `actN_spearman` and `actN_top33_*` / `actN_top10_*` columns (share of T cells within 30 µm of high-ratio neutrophils vs a shuffled-neutrophil null, the analysis/06 test).

# Protein recipe v3p (adds to v2p.1; `src/bio/recipe_protein.py`, `web/recipe.js`)

v2p.1 columns keep their names and exact values. v3p adds, per patient:
- **Compartments** from cell positions and the tumor-cell label, fixed thresholds: *nest* = tumor cells are at least
  half of the cells within 50 µm; *margin* = within 100 µm of a nest cell; *far* = the rest. Images without tumor-cell
  labels get no compartment results.
- **Composition** per compartment: percent T, CD8, CD4, Treg, neutrophil, macrophage and B cells, of all cells and of
  immune cells (`comp_<nest|margin|far>_pct_<role>`, `comp_<c>_immshare_<role>`).
- **Nest-restricted proximity**: the v2p.1 tests for T cells inside the nest, and per T subset (CD8, CD4, Treg), with
  the chance level from labels shuffled within each compartment of the image (`nest_<sub>_near_<src>_log2`).
- **Distance gradient**: share of nest T cells, and their marker percentiles, in 0–5, 5–10, 10–20, 20–30, 30–50, 50–100
  and 100+ µm from the nearest neutrophil or macrophage, each against chance (`grad_<src>_<bin>_*`).
- **Nest contact tests** per subset (`nest_<src>Touch_<sub>_<state>_diff`), with the bleed columns as before.
- **Macrophage state next to T cells**: PD-L1, IDO1, CD163, HLA-DR percentiles in macrophages under 10 µm from a T cell
  minus those 30 µm or more away (`TmacTouch_<marker>_diff`).

Why: the questions are about T cells inside the cancerous part of the tissue, about which T cells, and about how the
relationship changes with distance; whole-image numbers answer a different question. Limits: the nest is defined by
cell labels, so on whole sections it may include normal epithelium unless the labels separate malignant cells.
