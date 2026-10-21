"""Turn a dataset + manifest into the kit's cell table: one row per cell with
patient, field, x_um, y_um, raw_type, type (role), one column per state, one bleed_<role> column per bleed marker."""
import numpy as np
import pandas as pd


def load_cells(m):
    if m["modality"] == "protein":
        return _protein(m)
    return _rna(m)


def _protein(m):
    path = m["path"]
    if str(path).endswith(".h5ad"):
        import anndata as ad
        a = ad.read_h5ad(path)
        layer = m.get("layer")
        X = pd.DataFrame(a.layers[layer] if layer else _dense(a.X), columns=a.var_names)
        src = pd.concat([a.obs.reset_index(drop=True), X], axis=1)
    else:
        c = m["columns"]
        # keys as raw strings and correctly-rounded floats, so the browser version parses identically
        src = pd.read_csv(path, dtype={c[k]: str for k in ("patient", "field", "type")}, float_precision="round_trip")
    for col, val in m.get("filters", {}).items():
        src = src[src[col].astype(str) == str(val)]
    src = src.reset_index(drop=True)
    c = m["columns"]
    px = float(m["pixel_size_um"])
    df = pd.DataFrame({
        "patient": src[c["patient"]].astype(str).values, "field": src[c["field"]].astype(str).values,
        "x_um": src[c["x"]].values * px, "y_um": src[c["y"]].values * px, "raw_type": src[c["type"]].astype(str).values,
    })
    for s, col in m.get("states", {}).items():
        df[s] = src[col].values
    for r, col in m.get("bleed", {}).items():
        df[f"bleed_{r}"] = src[col].values
    return _roles(df, m)


def _rna(m):
    import anndata as ad
    a = ad.read_h5ad(m["path"])
    for col, val in m.get("filters", {}).items():
        a = a[a.obs[col].astype(str) == str(val)].copy()
    X = a.X.tocsc()
    tot = np.maximum(np.asarray(X.sum(1)).ravel(), 1)
    def cnt(genes):
        idx = [a.var_names.get_loc(g) for g in genes if g in a.var_names]
        return np.asarray(X[:, idx].sum(1)).ravel() if idx else np.full(a.n_obs, np.nan)
    c = m["columns"]
    def get(spec):
        if str(spec).startswith("obsm:"):
            _, key, i = spec.split(":")
            return np.asarray(a.obsm[key])[:, int(i)]
        if spec in a.obs:
            return a.obs[spec].astype(str).values
        return np.full(a.n_obs, spec)  # a constant, e.g. one patient per file
    px = float(m["pixel_size_um"])
    df = pd.DataFrame({"patient": get(c["patient"]), "field": get(c["field"]),
                       "x_um": get(c["x"]).astype(float) * px, "y_um": get(c["y"]).astype(float) * px,
                       "raw_type": a.obs[c["type"]].astype(str).values})
    ident = cnt(m["identity_genes"])
    for s, genes in m.get("states", {}).items():
        df[s] = np.log((cnt(genes) + 0.5) / (ident + 0.5))  # relative to T identity: cancels dilution
    for r, genes in m.get("bleed", {}).items():
        df[f"bleed_{r}"] = 1000 * cnt(genes) / tot
    return _roles(df, m)


def _roles(df, m):
    df["type"] = df["raw_type"].map(m["roles"]).fillna("unmapped")
    df = df[df["type"] != "drop"]
    df = df[np.isfinite(df["x_um"].astype(float)) & np.isfinite(df["y_um"].astype(float))].reset_index(drop=True)  # no position, no distance
    min_f = int(m.get("min_cells_per_field", 0) or 0)
    if min_f > 0:  # drop tiny images (same rule as web/recipe.js); an image is patient + image name
        n = df.groupby(["patient", "field"])["field"].transform("size")
        df = df[n >= min_f].reset_index(drop=True)
    return df


def _dense(X):
    return X.toarray() if hasattr(X, "toarray") else np.asarray(X)
