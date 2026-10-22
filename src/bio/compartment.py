"""Tumor compartments from cell positions (recipe v3p; see docs/recipe-v3.md). Ported 1:1 to web/recipe.js.

Per image:
  f50    = share of tumor-role cells among all cells within 50 um of the cell (self included, distance <= 50)
  nest   = f50 >= 0.5
  margin = not nest, and a nest cell within 100 um (distance <= 100)
  far    = everything else
Images with no tumor-role cell are "unknown" (compartments_for) and are skipped by the compartment metrics.
Thresholds are fixed in the recipe, not tuned per dataset.
"""
import numpy as np
from scipy.spatial import cKDTree

F50_RADIUS_UM = 50.0
NEST_MIN_FRAC = 0.5
MARGIN_RADIUS_UM = 100.0
COMPARTMENTS = ["nest", "margin", "far"]
ORDER = ["nest", "margin", "far", "unknown"]  # fixed order used when drawing stratified shuffles


def assign(x_um, y_um, is_tumor):
    """One image. Returns a numpy array of str, each "nest", "margin" or "far"."""
    xy = np.column_stack([np.asarray(x_um, float), np.asarray(y_um, float)])
    tumor = np.asarray(is_tumor, bool)
    n = len(xy)
    out = np.full(n, "far", dtype="<U6")
    if n == 0:
        return out
    # counts within 50 um (distance <= 50, the cell itself included), so f50 = n_tumor / n_all
    n_all = cKDTree(xy).query_ball_point(xy, r=F50_RADIUS_UM, return_length=True)
    n_tum = cKDTree(xy[tumor]).query_ball_point(xy, r=F50_RADIUS_UM, return_length=True) if tumor.any() else np.zeros(n, int)
    nest = n_tum / n_all >= NEST_MIN_FRAC
    if nest.any():
        d_nest = cKDTree(xy[nest]).query(xy, k=1)[0]
        out[(~nest) & (d_nest <= MARGIN_RADIUS_UM)] = "margin"
    out[nest] = "nest"
    return out


def compartments_for(df):
    """Fill df["comp"] per image (patient + field) from geometry and the "tumor" type. Images without tumor cells get
    "unknown". If df already has a "comp" column (explicit annotation from the importer) it is kept as is.
    Adds the column in place and returns df."""
    if "comp" in df.columns:
        return df
    comp = np.full(len(df), "unknown", dtype="<U7")
    tumor = (df["type"] == "tumor").to_numpy()
    x = df["x_um"].to_numpy(float)
    y = df["y_um"].to_numpy(float)
    for _, ix in df.groupby(["patient", "field"]).indices.items():  # positional indices
        if not tumor[ix].any():
            continue
        comp[ix] = assign(x[ix], y[ix], tumor[ix])
    df["comp"] = comp
    return df
