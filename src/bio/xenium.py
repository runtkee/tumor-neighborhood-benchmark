"""Read a 10x Xenium "outs" bundle into AnnData (coordinates in µm, raw counts, negative controls kept apart).

Public Xenium datasets (10x demo data, Zenodo cohorts) ship as a folder with
  cell_feature_matrix.h5   10x HDF5 matrix: group `matrix` with data / indices / indptr / shape / barcodes / features/*
  cells.csv.gz or cells.parquet   one row per cell: cell_id, x_centroid, y_centroid (µm), transcript_counts, cell_area, ...
and no cell-type labels. `read_bundle` returns:
  X                    CSR, raw counts, cells x genes (feature type "Gene Expression" only)
  obs                  every column of the cells table (except cell_id, which becomes the index), plus
                       neg_control_counts / neg_control_n_features (summed negative-control probes and codewords)
  obsm["spatial"]      [x_centroid, y_centroid] in µm
If the bundle has analysis/clustering/gene_expression_graphclust/clusters.csv (10x graph-based clusters, columns
Barcode,Cluster), the cluster of each cell is added as obs["cluster_10x"] (string, NaN for cells the file lacks);
bio.annotate then labels clusters instead of single cells. Labels are added afterwards by bio.annotate.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

MATRIX = "cell_feature_matrix.h5"
GENE_TYPE = "Gene Expression"


def _find(root, names):
    """First existing file among `names` in root, root/outs, then anywhere up to two levels down."""
    for base in (root, root / "outs"):
        for n in names:
            if (base / n).is_file():
                return base / n
    for n in names:
        hits = sorted(root.glob(f"*/{n}")) + sorted(root.glob(f"*/*/{n}"))
        if hits:
            return hits[0]
    return None


def _decode(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in np.asarray(a).ravel()], dtype=object)


def _read_matrix(h5_path):
    import h5py
    with h5py.File(h5_path, "r") as f:
        if "matrix" not in f:
            raise ValueError(f"{h5_path}: no 'matrix' group; this is not a 10x HDF5 matrix (groups: {list(f)})")
        g = f["matrix"]
        need = ["data", "indices", "indptr", "shape", "barcodes"]
        missing = [k for k in need if k not in g] + [k for k in ("name", "feature_type") if k not in g.get("features", {})]
        if missing:
            raise ValueError(f"{h5_path}: matrix is missing {missing}")
        n_feat, n_cell = (int(v) for v in g["shape"][:])
        # 10x stores features x cells in CSC order (columns = cells)
        m = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=(n_feat, n_cell))
        return (m.T.tocsr(), _decode(g["barcodes"][:]), _decode(g["features/name"][:]),
                _decode(g["features/feature_type"][:]),
                _decode(g["features/id"][:]) if "id" in g["features"] else None)


def _read_cells(root):
    csv = _find(root, ["cells.csv.gz", "cells.csv"])
    pq = _find(root, ["cells.parquet"])
    if csv is not None:
        return pd.read_csv(csv), csv
    if pq is not None:
        try:
            import pyarrow  # noqa: F401
        except ImportError as e:
            raise ImportError(f"found {pq} but pyarrow is not installed, and there is no cells.csv.gz next to it. "
                              "Install pyarrow (uv add pyarrow) or export/download cells.csv.gz.") from e
        return pd.read_parquet(pq), pq
    raise FileNotFoundError(f"no cells.csv.gz / cells.csv / cells.parquet under {root}")


CLUSTER_COL = "cluster_10x"
GRAPHCLUST = "analysis/clustering/gene_expression_graphclust/clusters.csv"


def _skip_clusters(clusters):
    return clusters is None or clusters is False or str(clusters).lower() == "none"


def read_clusters(csv_path, barcodes):
    """clusters.csv (Barcode,Cluster) -> (Categorical of cluster ids as strings aligned to barcodes, n cells without one).
    Raises ValueError when the file has no Barcode / Cluster columns or shares no barcode with the cells."""
    df = pd.read_csv(csv_path, dtype=str)
    if not {"Barcode", "Cluster"} <= set(df.columns):
        raise ValueError(f"{Path(csv_path).name}: expected the columns Barcode,Cluster (found {list(df.columns)[:6]})")
    s = df.drop_duplicates("Barcode").set_index("Barcode")["Cluster"]
    got = s.reindex(pd.Index(barcodes))
    n_missing = int(got.isna().sum())
    if n_missing == len(got):
        raise ValueError(f"{Path(csv_path).name}: none of its {len(s)} barcodes match the {len(got)} cells "
                         f"(e.g. {list(s.index[:2])} vs {list(barcodes[:2])}); it is from a different run")
    return pd.Categorical(got.to_numpy()), n_missing


def read_bundle(path, clusters="auto"):
    """Xenium outs directory (or its parent), or an .h5ad path (returned as is) -> AnnData.

    clusters: "auto" loads analysis/clustering/gene_expression_graphclust/clusters.csv when the bundle has it (a file
    that does not match the cells is skipped with a UserWarning); a path loads that file (errors are raised);
    None / "none" does not load clusters (and drops obs["cluster_10x"] from an .h5ad)."""
    import warnings

    import anndata as ad
    p = Path(path)
    skip = _skip_clusters(clusters)
    explicit = None if skip or str(clusters) == "auto" else Path(clusters)
    if explicit is not None and not explicit.is_file():
        raise FileNotFoundError(f"clusters file {explicit} not found")
    if p.suffix == ".h5ad":
        a = ad.read_h5ad(p)
        if skip and CLUSTER_COL in a.obs:
            del a.obs[CLUSTER_COL]
        if explicit is not None:
            a.obs[CLUSTER_COL], _ = read_clusters(explicit, a.obs_names)
        return a
    if not p.is_dir():
        raise FileNotFoundError(f"{p} is neither a directory nor an .h5ad file")
    h5 = _find(p, [MATRIX])
    if h5 is None:
        zarr = sorted(p.rglob("cell_feature_matrix.zarr*"))
        hint = (f" Found {zarr[0].name} but the HDF5 file is needed: download 'Cell feature matrix (HDF5)'."
                if zarr else "")
        raise FileNotFoundError(f"{MATRIX} not found under {p}.{hint}")
    X, barcodes, names, ftypes, ids = _read_matrix(h5)
    cells, cells_path = _read_cells(h5.parent if _find(h5.parent, ["cells.csv.gz", "cells.csv", "cells.parquet"]) else p)

    for c in ("cell_id", "x_centroid", "y_centroid"):
        if c not in cells.columns:
            raise ValueError(f"{cells_path.name}: column {c!r} missing (columns: {list(cells.columns)[:12]})")
    cells = cells.copy()
    cells["cell_id"] = _decode(cells["cell_id"].to_numpy())
    cells = cells.drop_duplicates("cell_id").set_index("cell_id")
    # hash lookup (np.setdiff1d on object arrays is quadratic: hours for 1M cells); sorted + unique like setdiff1d
    missing = np.sort(pd.unique(barcodes[~pd.Index(barcodes).isin(cells.index)]).astype(object)) if len(barcodes) else barcodes[:0]
    if len(missing):
        raise ValueError(f"{len(missing)} of {len(barcodes)} cells in {MATRIX} have no row in {cells_path.name} "
                         f"(e.g. {missing[:3].tolist()}); the two files are from different runs.")
    cells = cells.loc[barcodes]
    for col in cells.columns:  # bytes in parquet -> str, so the h5ad can be written
        if cells[col].dtype == object:
            cells[col] = _decode(cells[col].to_numpy())

    gene = ftypes == GENE_TYPE
    if not gene.any():
        raise ValueError(f"{MATRIX}: no features of type {GENE_TYPE!r} (types: {sorted(set(ftypes))})")
    neg = np.array([("negative control" in t.lower()) for t in ftypes])
    obs = cells.copy()
    obs["neg_control_counts"] = np.asarray(X[:, np.where(neg)[0]].sum(axis=1)).ravel().astype(np.float32) if neg.any() \
        else np.zeros(len(obs), np.float32)
    obs["neg_control_n_features"] = int(neg.sum())
    obs.index = pd.Index(barcodes, name=None)

    var = pd.DataFrame(index=pd.Index(names[gene]))
    if ids is not None:
        var["gene_id"] = ids[gene]
    Xg = X[:, np.where(gene)[0]].tocsr().astype(np.float32)
    a = ad.AnnData(X=Xg, obs=obs, var=var)
    if not a.var_names.is_unique:
        a.var_names_make_unique()
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=float)
    info = {"source": str(h5.parent), "cells_table": cells_path.name,
            "feature_types": {t: int((ftypes == t).sum()) for t in sorted(set(ftypes))}}
    cl_path = None
    if not skip:
        cl_path = explicit or _find(h5.parent, [GRAPHCLUST]) or _find(p, [GRAPHCLUST])
    if cl_path is not None:
        try:
            a.obs[CLUSTER_COL], n_missing = read_clusters(cl_path, a.obs_names)
            info.update(clusters_file=str(cl_path), clusters_cells_without_cluster=n_missing)
        except ValueError as e:
            if explicit is not None:
                raise
            warnings.warn(f"{e}; clusters ignored, cells are labelled one by one", UserWarning, stacklevel=2)
    a.uns["xenium_bundle"] = json.dumps(info)
    return a


def concat_sections(adatas, names):
    """Stack per-section AnnData (already annotated) into one cohort object.

    Genes: outer join, a gene missing from a section's panel is 0 there (reported in uns['panel_differences'] and a
    UserWarning). Cell ids are made unique by appending _<name>. Coordinates stay per section (distinct images)."""
    import warnings

    import anndata as ad
    panels = [set(a.var_names) for a in adatas]
    union = set.union(*panels)
    missing = {n: sorted(union - p) for n, p in zip(names, panels) if union - p}
    if missing:
        warnings.warn("panels differ between bundles; missing genes are treated as 0 in those sections: "
                      + "; ".join(f"{n}: {', '.join(g[:6])}{' ...' if len(g) > 6 else ''}" for n, g in missing.items()),
                      UserWarning, stacklevel=2)
    c = ad.concat(dict(zip(names, adatas)), join="outer", fill_value=0, index_unique="_")
    c.uns["panel_differences"] = json.dumps(missing)
    c.uns["xenium_bundles"] = json.dumps([json.loads(a.uns["xenium_bundle"]) if "xenium_bundle" in a.uns else {}
                                          for a in adatas])
    return c


# ---------------------------------------------------------------------------------------------------------------
# Sample areas: cut a slide that carries several patients into one section per selection polygon.
# ---------------------------------------------------------------------------------------------------------------
def read_areas(path):
    """Selection polygons from a `<region>_coordinates.csv` (10x Xenium Explorer selection export).

    Format: comment lines starting with `#` (selection names, areas), a header `Selection,X,Y`, then one row per
    polygon vertex in µm (the frame of x_centroid / y_centroid). Consecutive rows with the same Selection form one
    polygon; a name that comes back later in the file is a further polygon of the same selection.
    Returns a list of (name, vertices[n,2]) in file order."""
    import csv
    p = Path(path)
    polys, header_seen, cur, verts = [], False, None, []
    with open(p, newline="", encoding="utf-8-sig") as fh:
        for ln, row in enumerate(csv.reader(fh), 1):
            if not row or not "".join(row).strip() or row[0].lstrip().startswith("#"):
                continue
            if not header_seen:
                if [c.strip().lower() for c in row[:3]] != ["selection", "x", "y"]:
                    raise ValueError(f"{p.name} line {ln}: expected the header 'Selection,X,Y', found {row[:3]}")
                header_seen = True
                continue
            if len(row) < 3:
                raise ValueError(f"{p.name} line {ln}: expected Selection,X,Y, found {row}")
            try:
                xy = (float(row[1]), float(row[2]))
            except ValueError:
                raise ValueError(f"{p.name} line {ln}: X / Y are not numbers: {row[1:3]}") from None
            name = row[0].strip()
            if name != cur:
                if cur is not None:
                    polys.append((cur, np.array(verts, float)))
                cur, verts = name, []
            verts.append(xy)
    if cur is not None:
        polys.append((cur, np.array(verts, float)))
    if not header_seen or not polys:
        raise ValueError(f"{p.name}: no 'Selection,X,Y' header and polygon rows found")
    short = [n for n, v in polys if len(v) < 3]
    if short:
        raise ValueError(f"{p.name}: selection(s) with fewer than 3 vertices: {sorted(set(short))}")
    return polys


_PIP_CELLS = 100_000     # points x edges evaluated at once (small enough for the CPU cache)


def _contains_numpy(xy, verts):
    """Even-odd ray casting; xy [n,2], verts [m,2]. Points exactly on an edge may fall either way.
    Fully vectorised: edges x points in chunks, no Python loop over points (or over edges)."""
    inside = np.zeros(len(xy), bool)
    box = ((xy[:, 0] >= verts[:, 0].min()) & (xy[:, 0] <= verts[:, 0].max())
           & (xy[:, 1] >= verts[:, 1].min()) & (xy[:, 1] <= verts[:, 1].max()))
    if not box.any():
        return inside
    x1, y1 = verts[:, 0], verts[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
    keep = y1 != y2                                  # horizontal edges never cross a horizontal ray
    x1, y1, x2, y2 = x1[keep], y1[keep], x2[keep], y2[keep]
    slope = (x2 - x1) / (y2 - y1)
    pts = np.flatnonzero(box)
    step = max(1, _PIP_CELLS // max(len(x1), 1))
    for s in range(0, len(pts), step):
        sel = pts[s:s + step]
        px, py = xy[sel, 0][:, None], xy[sel, 1][:, None]
        cross = (y1[None, :] > py) != (y2[None, :] > py)
        xint = x1[None, :] + (py - y1[None, :]) * slope[None, :]
        inside[sel] = ((cross & (px < xint)).sum(axis=1) & 1).astype(bool)
    return inside


def points_in_polygon(xy, verts, backend="auto"):
    """Boolean mask: which points of xy [n,2] lie inside the polygon verts [m,2]. backend: auto (matplotlib if it
    is installed, else numpy ray casting), matplotlib, numpy."""
    xy = np.asarray(xy, float)
    if backend in ("auto", "matplotlib"):
        try:
            from matplotlib.path import Path as MplPath
        except ImportError:
            if backend == "matplotlib":
                raise
        else:
            inside = np.zeros(len(xy), bool)       # only the points in the polygon's bounding box need the exact test
            box = ((xy[:, 0] >= verts[:, 0].min()) & (xy[:, 0] <= verts[:, 0].max())
                   & (xy[:, 1] >= verts[:, 1].min()) & (xy[:, 1] <= verts[:, 1].max()))
            if box.any():
                inside[box] = MplPath(verts).contains_points(xy[box])
            return inside
    return _contains_numpy(xy, verts)


def assign_areas(xy, polys, backend="auto"):
    """For each point the index of the selection (in order of first appearance of the names) whose polygon holds
    it, or -1 outside every polygon. A point inside several polygons goes to the first selection in file order
    (the count of such points is returned). Returns (index[n], names, n_overlap)."""
    names = list(dict.fromkeys(n for n, _ in polys))
    which = np.full(len(xy), -1, np.int32)
    hits = np.zeros(len(xy), np.int16)
    for name, verts in polys:
        m = points_in_polygon(xy, verts, backend)
        hits += m
        which[m & (which < 0)] = names.index(name)
    return which, names, int((hits > 1).sum())


def cut_by_areas(adata, polys, backend="auto"):
    """Split one slide into sections by selection polygons -> (list of (name, AnnData), info).
    info: n_cells, n_outside (dropped), n_overlap, counts {name: cells}. Selections with no cell are left out and
    listed in info['empty']."""
    which, names, n_overlap = assign_areas(adata.obsm["spatial"], polys, backend)
    out, counts, empty = [], {}, []
    for i, n in enumerate(names):
        m = which == i
        counts[n] = int(m.sum())
        if counts[n] == 0:
            empty.append(n)
            continue
        out.append((n, adata[m].copy()))
    info = {"n_cells": int(adata.n_obs), "n_outside": int((which < 0).sum()), "n_overlap": n_overlap,
            "counts": counts, "empty": empty}
    return out, info
