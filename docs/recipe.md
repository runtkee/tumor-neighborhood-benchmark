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
