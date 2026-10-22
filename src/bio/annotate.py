"""Deterministic marker-based lineage labels for imaging-transcriptomics cells (Xenium, CosMx). No clustering.

Per cell and lineage: score = (sum of the lineage's marker counts present in the panel) per 1000 counts.
Background (when the file has negative-control probes) = mean negative-probe rate per feature x number of markers,
on the same per-1000 scale. A lineage "fires" if score > max(3 x background, 1 per 1000), at least 2 of its markers
are detected (count >= 1), and, for T, at least one CD3 gene is detected (CD2/CD247 alone also mark NK cells).
The label is the lineage with the highest score / threshold ratio. Ambiguous = nothing fired, or the runner-up is
within 0.75 of the best (runner-ups inside the same tie group do not count: Macrophage/Monocyte/DC and T/NK
go to the higher ratio). LowCount = fewer than `min_counts` counts. T cells are split into CD8T / CD4T / Treg
(per cell with the cell method; with the cluster method the cluster decides "T" and each cell decides its subset, "T"
= no subset evidence).

Two methods. "cell" is the per-cell scoring above. "cluster" (the default when obs has `cluster_10x`, the 10x
graph-based clusters) labels whole clusters: for each cluster and lineage, m = mean per-1000 expression of the
lineage's marker set over the cluster's cells; scaled = m / (the largest m of that lineage over all clusters), so
every present lineage has one cluster at 1.0 however common the lineage is (the largest m is taken over the
clusters with at least max(CL_REF_MIN_CELLS, CL_REF_MIN_FRAC x the section's non-LowCount cells) cells, so a handful-cell
outlier cluster cannot become the reference; smaller clusters are still labelled and may scale above 1.0). A cluster is called lineage L when L has
the highest scaled value (>= CL_MIN_SCALED), m is above the lineage floor, enough of L's markers (>= MIN_MARKERS, or
the single SOLO_OK marker) are each >= CL_MARKER_FRAC of that marker's largest cluster mean (CL_MARKER_FRAC_2 when
the panel has exactly 2 markers of L), and the runner-up's scaled value is < CL_AMBIG x the best; otherwise Ambiguous.
A would-be tie is broken by the per-cell firing share: the best lineage is called when the share of the cluster's cells
in which it fires is >= CL_TIEBREAK x the runner-up's (note "tiebreak: fired"). A lineage whose largest cluster mean is below
CL_ABSENT_MULT x its floor is "absent" from the section and is never called. Cells inherit their cluster's label;
LowCount still applies per cell, and cells the clusters file lacks are "Unclustered" (no per-cell fallback).
Cluster labelling suits sparse panels, where single cells carry too few counts for per-cell scoring.
Inside a T cluster each cell gets its own subset (T_SUBSET_RULE): Treg if FOXP3 count >= 2 and no CD8A/CD8B counts (panels
without FOXP3: IL2RA >= 2 and CTLA4 >= 2, same CD8 exclusion), else CD8T if CD8A+CD8B > CD4, else CD4T if CD4 > CD8A+CD8B, else "T" (subset unknown).
A cluster whose T scaled value is >= CL_NK_T_RATIO x its NK value is T; NK also needs CL_NK_MIN_SPECIFIC enriched NK
markers besides GNLY.

These are marker calls, not a tumor-vs-normal call. `annotate_report` gives a marker-leakage table: how often each
lineage's markers fire in cells labelled as something else. High values mean the labels (or the signal) are not
reliable (the CosMx-6K neutrophil pitfall). With method="cell" a lineage firing in many other labels is flagged
unreliable; with method="cluster" a lineage is flagged when its clusters are not clearly separated from the runner-up lineage
(median of best minus runner-up scaled value below CL_UNRELIABLE_GAP).
"""
import json

import numpy as np
import pandas as pd
from scipy import sparse

LABELS = ["Epithelial", "Fibroblast", "Endothelial", "SMC", "CD8T", "CD4T", "Treg", "T", "B", "Plasma", "Macrophage",
          "Monocyte", "DC", "Neutrophil", "Mast", "NK", "Unclustered", "Ambiguous", "LowCount"]
LINEAGE_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7", "KRT17"],
    "Fibroblast": ["DCN", "LUM", "PDGFRA", "PDGFRB", "FAP", "THY1", "COL5A1", "COL6A1", "COL1A1", "COL1A2", "COL3A1"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"],
    "SMC": ["ACTA2", "MYH11", "TAGLN", "DES"],
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "CD247"],
    "B": ["MS4A1", "CD79A", "CD19", "BANK1"],
    # XBP1 and SDC1 are left out (epithelium and many other cells express them); immunoglobulin genes too, see AMBIENT_IG.
    # Both markers must be detected per cell (MIN_MARKERS = 2), and both must be enriched in a cluster.
    "Plasma": ["JCHAIN", "MZB1"],
    "Macrophage": ["CD68", "CD163", "C1QA", "C1QB", "C1QC", "MRC1", "MS4A7", "CSF1R"],
    "Monocyte": ["FCN1", "VCAN", "S100A12", "LYZ"],
    "DC": ["CLEC9A", "CD1C", "LAMP3", "CCR7", "FLT3"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CXCR2", "CXCR1", "FPR1", "CMTM2", "S100A8", "S100A9", "G0S2"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "NK": ["NCAM1", "KLRD1", "KLRF1", "GNLY", "NCR1"],
}
LINEAGES = list(LINEAGE_MARKERS)
# Immunoglobulin transcripts float free in tumour tissue and land in nearly every cell, so they are not used for Plasma.
AMBIENT_IG = ["IGKC", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "IGLC1", "IGLC2", "IGLC3"]
# A lineage with a single marker in the panel is still called when that marker is specific for the lineage.
SOLO_OK = {"EPCAM", "PECAM1", "MYH11", "ACTA2", "MS4A1", "CPA3"}
SOLO_MIN_COUNTS = 2     # counts of the single marker needed in a cell (one count is too easily background)
LABEL_LINEAGE = {"CD8T": "T", "CD4T": "T", "Treg": "T"}  # every other real label is its own lineage
CD3 = ["CD3D", "CD3E", "CD3G"]
TIE_GROUPS = [{"Macrophage", "Monocyte", "DC"}, {"T", "NK"}]
EXTRA_GENES = ["CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4"]

MIN_MARKERS = 2        # markers detected for a lineage to fire (and needed in the panel)
BG_MULT = 3.0          # score must exceed this x background
FLOOR = 1.0            # ... and this many counts per 1000 (the only bar when there are no negative controls)
AMBIG_RATIO = 0.75     # runner-up within this fraction of the best -> Ambiguous
UNRELIABLE_FRAC = 0.10  # a lineage fires in more than this share of ...
UNRELIABLE_LABELS = 3   # ... this many other labels -> unreliable
MIN_GROUP = 50          # other labels with fewer cells are shown but not counted
CHUNK = 200_000
CL_MIN_SCALED = 0.4     # cluster method: the best lineage's cluster mean / its largest cluster mean at least this
CL_MARKER_FRAC = 0.5    # ... with at least MIN_MARKERS (or the single SOLO_OK marker) markers at this fraction of their largest cluster mean
CL_MARKER_FRAC_2 = 0.3  # ... but a lineage with exactly 2 markers in the panel (e.g. Plasma) needs both at this fraction
CL_TIEBREAK = 1.5       # a would-be tie is called for the best lineage when its per-cell firing share is >= this x the runner-up's
CL_ABSENT_MULT = 2.0    # a lineage whose largest cluster mean is below this x its floor is absent from the section
CL_AMBIG = 0.75         # ... and the runner-up scaled value below this fraction of the best
CL_REF_MIN_CELLS = 50   # the per-lineage reference maximum is taken over clusters with at least this many usable cells ...
CL_REF_MIN_FRAC = 0.005  # ... and at least this share of the section's non-LowCount cells (whichever is larger)
T_SUBSET_RULE = ("per cell inside T clusters: Treg if FOXP3 count >= 2 and CD8A+CD8B = 0 (no FOXP3 in panel: IL2RA >= 2 and CTLA4 >= 2, same CD8 exclusion); "
                 "else CD8T if CD8A+CD8B > CD4; else CD4T if CD4 > CD8A+CD8B; else T (no subset evidence)")
CL_NK_T_RATIO = 0.5     # a cluster whose T scaled value is >= this x its NK scaled value is T, not NK
CL_NK_MIN_SPECIFIC = 2  # NK needs this many enriched markers other than GNLY (GNLY is also made by cytotoxic CD8 T cells)
NK_NOT_COUNTED = {"GNLY"}
CL_UNRELIABLE_GAP = 0.25  # a lineage whose clusters have a median (best - runner-up scaled) below this is unreliable
METHODS = ("auto", "cluster", "cell")
CLUSTER_COL = "cluster_10x"
MAX_CLUSTER_ROWS = 60   # report_text prints the Clusters table up to this many rows


def _tie_group(lin):
    for g in TIE_GROUPS:
        if lin in g:
            return g
    return {lin}


def _check_counts(X):
    d = X.data[:100_000] if sparse.issparse(X) else np.asarray(X).ravel()[:100_000]
    if len(d) and (np.any(d < 0) or not np.allclose(d, np.round(d))):
        raise ValueError("annotate needs raw integer counts in X (found normalized or negative values)")


def _nat_key(x):
    x = str(x)
    return (0, int(x), x) if x.isdigit() else (1, 0, x)


def _label_clusters(gsum, cnt, cluster_names, used, solo, idx, pos, thr, extras, gdet=None, fshare=None, n_usable=None):
    """Cluster-level calls from per-cluster sums of per-1000 marker rates (gsum [K, n_cols], cnt [K] cells counted;
    gdet [K, n_cols] = cells with count >= 1 of each column, among the same cells;
    fshare [K, n_lineages] = share of those cells in which the per-cell scoring fires each lineage, for the tie-break;
    n_usable = the section's non-LowCount cells, default cnt.sum()).
    The reference maximum of each lineage (and of each marker) is taken over the clusters with at least
    max(CL_REF_MIN_CELLS, CL_REF_MIN_FRAC x n_usable) cells (all clusters if none is that large; likewise, a lineage that
    is below its absent floor in every such cluster takes its reference from all clusters, so a lineage that only forms
    small clusters is not lost); smaller clusters are scaled against the reference too, so their scaled value can exceed 1.0.
    Returns (list of dicts, one per cluster in the order of cluster_names; list of absent lineages; reference info).
    T clusters are returned with label "T"; the subset is decided per cell by the caller."""
    K = len(cluster_names)
    has = cnt > 0
    mean_c = gsum / np.maximum(cnt, 1)[:, None]
    if gdet is None:
        gdet = np.zeros_like(gsum)
    if fshare is None:
        fshare = np.zeros((K, len(LINEAGES)))
    n_usable = float(cnt.sum() if n_usable is None else n_usable)
    ref_min = max(CL_REF_MIN_CELLS, CL_REF_MIN_FRAC * n_usable)
    ref = has & (cnt >= ref_min)
    fallback = bool(has.any() and not ref.any())
    if fallback:
        ref = has
    mean_c_max = mean_c[ref].max(axis=0) if ref.any() else np.zeros(mean_c.shape[1])
    m = np.zeros((K, len(LINEAGES)))
    for li, l in enumerate(LINEAGES):
        j = [pos[idx[g]] for g in used[l]]
        if len(j):
            m[:, li] = mean_c[:, j].sum(axis=1)
    m[~has] = 0.0
    m_max = m[ref].max(axis=0) if ref.any() else np.zeros(len(LINEAGES))
    # A lineage that is below its absent floor in every large cluster but present in a small one (a rare cell type that
    # only forms small clusters) takes its reference from all clusters, as the old rule did; otherwise it would be absent.
    lin_fallback = []
    if not fallback and has.any():
        m_all = m[has].max(axis=0)
        for li, l in enumerate(LINEAGES):
            if m_max[li] < CL_ABSENT_MULT * thr[li] <= m_all[li]:
                m_max[li] = m_all[li]
                lin_fallback.append(l)
                for g in used[l]:
                    mean_c_max[pos[idx[g]]] = mean_c[has, pos[idx[g]]].max()
    absent = [l for li, l in enumerate(LINEAGES)
              if not (len(used[l]) < MIN_MARKERS and l not in solo) and m_max[li] < CL_ABSENT_MULT * thr[li]]
    callable_ = np.array([not (len(used[l]) < MIN_MARKERS and l not in solo) and l not in absent for l in LINEAGES])
    with np.errstate(divide="ignore", invalid="ignore"):
        scaled = np.where(callable_[None, :] & (m_max[None, :] > 0), m / m_max[None, :], 0.0)

    rows = []
    for k in range(K):
        r = scaled[k]
        order = np.argsort(-r, kind="stable")
        top3 = [[LINEAGES[j], round(float(r[j]), 3)] for j in order[:3]]
        row = {"cluster": cluster_names[k], "n": 0, "n_used": int(cnt[k]), "label": "Ambiguous", "note": "",
               "best_scaled": 0.0, "second_scaled": 0.0, "lineage": "", "top3": top3}
        if cnt[k] == 0:
            row["note"] = "no usable cells (all below min_counts)"
            rows.append(row)
            continue
        b = int(order[0])
        L = LINEAGES[b]
        if L == "NK" and "T" in LINEAGES:   # a cytotoxic T cluster also scores high on GNLY: T unless clearly less T than NK
            tb = LINEAGES.index("T")
            if r[tb] >= CL_NK_T_RATIO * r[b]:
                b, L = tb, "T"
        best = float(r[b])
        others = [j for j, o in enumerate(LINEAGES) if o not in _tie_group(L)]
        second = max([float(r[j]) for j in others], default=0.0)
        run = max(others, key=lambda j: (r[j], -j)) if others else None      # the runner-up lineage
        frac = CL_MARKER_FRAC_2 if len(used[L]) == 2 else CL_MARKER_FRAC
        n_ok = sum(mean_c[k, pos[idx[g]]] >= frac * mean_c_max[pos[idx[g]]] > 0
                   for g in used[L] if not (L == "NK" and g in NK_NOT_COUNTED))
        need = 1 if L in solo else (CL_NK_MIN_SPECIFIC if L == "NK" else MIN_MARKERS)
        row["best_scaled"], row["second_scaled"] = best, second
        if best <= 0 or best < CL_MIN_SCALED or m[k, b] < thr[b] or n_ok < need:
            row["note"] = "no lineage fired"
        elif second >= CL_AMBIG * best and not (
                run is not None and fshare[k, b] > 0 and fshare[k, b] >= CL_TIEBREAK * fshare[k, run]):
            row["note"] = "tie"
        else:
            if second >= CL_AMBIG * best:
                row["note"] = "tiebreak: fired"
            row["lineage"] = L
            row["label"] = L
        rows.append(row)

    ref_info = {"min_n": float(ref_min), "n_reference_clusters": int(ref.sum()), "n_clusters": int(has.sum()),
                "fallback_all": fallback, "fallback_lineages": lin_fallback}
    return rows, absent, ref_info


def _t_subset(ex, n):
    """Per-cell T subset from raw counts (ex = {gene: counts} for the panel's CD4, CD8A, CD8B, FOXP3, IL2RA, CTLA4):
    T_SUBSET_RULE. Returns an object array of CD8T / CD4T / Treg / T."""
    z = np.zeros(n, np.float32)
    g = lambda name: ex.get(name, z)
    cd8, cd4 = g("CD8A") + g("CD8B"), g("CD4")
    treg = (g("FOXP3") >= 2 if "FOXP3" in ex else (g("IL2RA") >= 2) & (g("CTLA4") >= 2)) & (cd8 == 0)
    return np.where(treg, "Treg", np.where(cd8 > cd4, "CD8T", np.where(cd4 > cd8, "CD4T", "T"))).astype(object)


def annotate(adata, min_counts=20, method="auto"):
    """Add obs['annotation'] (+ annotation_ratio, annotation_note), obsm['lineage_ratio'/'lineage_fired'] and
    uns['annotation_params'] (JSON). Returns adata.

    method: "cell" (per-cell scoring), "cluster" (label each obs['cluster_10x'] cluster, cells inherit) or "auto"
    (cluster when obs has cluster_10x, else cell). The per-cell scores are computed either way (the leakage table
    uses them)."""
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    if method == "auto":
        method = "cluster" if CLUSTER_COL in adata.obs else "cell"
    if method == "cluster" and CLUSTER_COL not in adata.obs:
        raise ValueError(f"method='cluster' needs obs[{CLUSTER_COL!r}] (10x graph-based clusters)")
    X = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))
    X = X.tocsr()
    _check_counts(X)
    n = adata.n_obs
    idx = {g: i for i, g in enumerate(adata.var_names)}
    used = {l: [g for g in m if g in idx] for l, m in LINEAGE_MARKERS.items()}
    solo = {l: used[l][0] for l in LINEAGES if len(used[l]) == 1 and used[l][0] in SOLO_OK}
    unavailable = [l for l in LINEAGES if len(used[l]) < MIN_MARKERS and l not in solo]
    ig_excluded = [g for g in AMBIENT_IG if g in idx]
    cd3 = [g for g in CD3 if g in idx]
    extras = [g for g in EXTRA_GENES if g in idx]

    cols = sorted({idx[g] for gs in used.values() for g in gs} | {idx[g] for g in extras})
    pos = {c: j for j, c in enumerate(cols)}
    tot = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    low = tot < min_counts

    codes, cl_names = None, []
    if method == "cluster":
        codes, uniq = pd.factorize(adata.obs[CLUSTER_COL].astype(object))
        codes = np.asarray(codes)
        cl_names = [str(u) for u in uniq]
        gsum = np.zeros((len(cl_names), len(cols)), np.float64)
        gdet = np.zeros((len(cl_names), len(cols)), np.float64)
        cnt = np.zeros(len(cl_names), np.int64)
        n_in_cluster = np.bincount(codes[codes >= 0], minlength=len(cl_names))

    # Per-cell lineage sums / detected-marker counts are sparse products (counts x marker-indicator matrix); nothing
    # loops over cells. Chunks bound the memory of the temporaries (and fix the order cluster sums are accumulated in).
    Xc = X[:, cols].tocsr()
    ind = np.zeros((len(cols), len(LINEAGES)), np.float32)
    for li, l in enumerate(LINEAGES):
        if l not in unavailable:
            ind[[pos[idx[g]] for g in used[l]], li] = 1.0
    ind = sparse.csr_matrix(ind)
    cd3_ind = np.zeros((len(cols), 1), np.float32)
    cd3_ind[[pos[idx[g]] for g in cd3], 0] = 1.0
    cd3_ind = sparse.csr_matrix(cd3_ind)
    sums = np.zeros((n, len(LINEAGES)), np.float64)
    ndet = np.zeros((n, len(LINEAGES)), np.int16)
    cd3_det = np.zeros(n, bool)
    ex = {g: np.asarray(X[:, idx[g]].toarray()).ravel().astype(np.float32) for g in extras}
    for s in range(0, n, CHUNK):
        e = min(s + CHUNK, n)
        blk = Xc[s:e]
        sums[s:e] = (blk @ ind).toarray()          # counts are integers: float32 sums are exact
        det = blk.copy()
        det.data = (det.data > 0).astype(np.float32)
        ndet[s:e] = (det @ ind).toarray()
        if cd3:
            cd3_det[s:e] = (det @ cd3_ind).toarray().ravel() > 0
        if codes is not None:  # per-cluster sums of per-1000 rates, over cells with enough counts
            ok = (codes[s:e] >= 0) & ~low[s:e]
            if ok.any():
                c_ok = codes[s:e][ok]
                scale = sparse.diags(1000.0 / np.maximum(tot[s:e][ok], 1.0))
                rate = scale @ blk[ok].astype(np.float64)
                oh = sparse.csr_matrix((np.ones(len(c_ok)), (c_ok, np.arange(len(c_ok)))), shape=(len(cl_names), len(c_ok)))
                gsum += (oh @ rate).toarray()
                dd = blk[ok].copy()
                dd.data = (dd.data > 0).astype(np.float64)
                gdet += (oh @ dd).toarray()
                cnt += np.bincount(c_ok, minlength=len(cl_names))
    score = 1000.0 * sums / np.maximum(tot, 1.0)[:, None]

    # background from negative-control probes: mean rate per feature per 1000 counts, times number of markers
    bg_rate, bg_source = 0.0, f"fixed floor of {FLOOR:g} per 1000 (no negative-control probes in the file)"
    if "neg_control_counts" in adata.obs and "neg_control_n_features" in adata.obs:
        nf = float(np.max(adata.obs["neg_control_n_features"].to_numpy()))
        if nf > 0 and tot.sum() > 0:
            bg_rate = 1000.0 * float(adata.obs["neg_control_counts"].sum()) / (nf * tot.sum())
            bg_source = (f"{int(nf)} negative-control features, mean {bg_rate:.4f} counts per feature per 1000 counts; "
                         f"threshold = max({BG_MULT:g} x background, {FLOOR:g}) per 1000")
    bg = np.array([bg_rate * len(used[l]) for l in LINEAGES])
    thr = np.maximum(BG_MULT * bg, FLOOR)

    fired = (score > thr[None, :]) & (ndet >= MIN_MARKERS)
    for l in solo:  # one marker: at least SOLO_MIN_COUNTS counts of it (sums = that marker's counts)
        li = LINEAGES.index(l)
        fired[:, li] = (score[:, li] > thr[li]) & (sums[:, li] >= SOLO_MIN_COUNTS)
    for l in unavailable:
        fired[:, LINEAGES.index(l)] = False
    if cd3:
        fired[:, LINEAGES.index("T")] &= cd3_det
    ratio = score / thr[None, :]
    masked = np.where(fired, ratio, 0.0)

    best = masked.argmax(axis=1)
    best_r = masked[np.arange(n), best]
    # runner-up: best ratio among lineages outside the winner's tie group
    second_r = np.zeros(n)
    for li, l in enumerate(LINEAGES):
        rows = best == li
        if not rows.any():
            continue
        other = [LINEAGES.index(o) for o in LINEAGES if o not in _tie_group(l)]
        second_r[rows] = masked[np.ix_(rows, other)].max(axis=1)
    lin = np.array(LINEAGES, dtype=object)[best]
    none_fired = best_r == 0
    tie = (~none_fired) & (second_r >= AMBIG_RATIO * best_r)

    t_default, t_no_subset, n_t_cluster = 0, 0, 0
    clusters, absent, ref_info = [], [], None
    if method == "cell":
        label = lin.copy()
        # T subsets
        is_t = (lin == "T") & ~none_fired & ~tie
        cd8 = ex.get("CD8A", 0) + ex.get("CD8B", 0)
        cd4 = ex.get("CD4", np.zeros(n, np.float32))
        cd8 = cd8 if isinstance(cd8, np.ndarray) else np.zeros(n, np.float32)
        foxp3 = ex.get("FOXP3", np.zeros(n, np.float32)) > 0
        il2ra_ctla4 = (ex.get("IL2RA", np.zeros(n, np.float32)) > 0) & (ex.get("CTLA4", np.zeros(n, np.float32)) > 0)
        sub = np.where(cd8 > cd4, "CD8T", np.where(foxp3 | il2ra_ctla4, "Treg", "CD4T"))
        label[is_t] = sub[is_t]
        t_default = int((is_t & (cd8 == 0) & (cd4 == 0) & ~foxp3 & ~il2ra_ctla4).sum())
        n_t = int(is_t.sum())
        label[none_fired | tie] = "Ambiguous"
        note = np.where(none_fired, "no lineage fired", np.where(tie, "tie", ""))
        out_ratio = best_r
    else:
        okc = (codes >= 0) & ~low           # the cells behind gsum / cnt
        oh = sparse.csr_matrix((np.ones(int(okc.sum())), (codes[okc], np.arange(int(okc.sum())))), shape=(len(cl_names), int(okc.sum())))
        fshare = (oh @ sparse.csr_matrix(fired[okc].astype(np.float64))).toarray() / np.maximum(cnt, 1)[:, None]
        clusters, absent, ref_info = _label_clusters(gsum, cnt, cl_names, used, solo, idx, pos, thr, extras, gdet, fshare,
                                                     n_usable=int((~low).sum()))
        k_lab = np.array([c["label"] for c in clusters] + ["Unclustered"], dtype=object)   # last = no cluster id
        k_note = np.array([c["note"] for c in clusters] + ["no cluster id"], dtype=object)
        k_ratio = np.array([c["best_scaled"] for c in clusters] + [0.0])
        c_idx = np.where(codes >= 0, codes, len(clusters))
        label, note, out_ratio = k_lab[c_idx].copy(), k_note[c_idx].copy(), k_ratio[c_idx]
        is_t = (label == "T") & ~low        # the cluster says "T"; each cell then decides its own subset
        n_t = n_t_cluster = int(is_t.sum())
        sub = _t_subset(ex, n)
        label[is_t] = sub[is_t]
        t_no_subset = int((sub[is_t] == "T").sum())
        for c in clusters:
            c["n"] = int(n_in_cluster[cl_names.index(c["cluster"])])
            if c["lineage"] == "T":
                in_c = is_t & (c_idx == cl_names.index(c["cluster"]))
                c["t_split"] = {l: int((sub[in_c] == l).sum()) for l in ("CD8T", "CD4T", "Treg", "T")}
    label[low] = "LowCount"
    note = np.asarray(note, dtype=object)
    note[low] = "total counts < min_counts"

    adata.obs["annotation"] = pd.Categorical(label.astype(str), categories=LABELS)
    adata.obs["annotation_ratio"] = np.asarray(out_ratio).astype(np.float32)
    adata.obs["annotation_note"] = pd.Categorical(note.astype(str))
    adata.obsm["lineage_ratio"] = ratio.astype(np.float32)
    adata.obsm["lineage_fired"] = fired.astype(np.uint8)
    clusters.sort(key=lambda c: _nat_key(c["cluster"]))
    adata.uns["annotation_params"] = json.dumps({
        "method": method, "clusters": clusters, "absent": absent,
        "cluster_rules": {"min_scaled": CL_MIN_SCALED, "marker_frac": CL_MARKER_FRAC, "marker_frac_2": CL_MARKER_FRAC_2, "tiebreak": CL_TIEBREAK, "absent_mult": CL_ABSENT_MULT,
                          "ambig": CL_AMBIG, "ref_min_cells": CL_REF_MIN_CELLS, "ref_min_frac": CL_REF_MIN_FRAC,
                          "t_subset": T_SUBSET_RULE, "nk_t_ratio": CL_NK_T_RATIO,
                          "nk_min_specific": CL_NK_MIN_SPECIFIC, "unreliable_gap": CL_UNRELIABLE_GAP},
        "reference": ref_info,
        "lineages": LINEAGES, "markers_used": used, "markers_requested": LINEAGE_MARKERS, "unavailable": unavailable,
        "single_marker": solo, "ambient_ig_excluded": {"Plasma": ig_excluded},
        "threshold_per_1000": dict(zip(LINEAGES, thr.round(4).tolist())),
        "background_per_1000": dict(zip(LINEAGES, bg.round(4).tolist())), "background_source": bg_source,
        "min_counts": min_counts, "cd3_genes": cd3, "t_subset_genes_present": extras,
        "t_default_cd4": t_default, "t_no_subset": t_no_subset, "n_t_cluster": n_t_cluster, "n_t": n_t})
    return adata


def _cell_unreliable(spec, real, avail, counts, min_group):
    out = {}
    for l in avail:
        bad = [o for o in real if counts[o] >= min_group and spec.loc[l, o] > UNRELIABLE_FRAC]
        if len(bad) >= UNRELIABLE_LABELS:
            out[l] = bad
    return out


def _cluster_unreliable(clusters):
    """lineage -> [reason] when the median (best - runner-up scaled value) over the clusters labelled with it is
    < CL_UNRELIABLE_GAP."""
    by = {}
    for c in clusters:
        if c["lineage"]:
            by.setdefault(c["lineage"], []).append(c["best_scaled"] - c["second_scaled"])
    return {l: [f"median best-minus-runner-up {np.median(v):.2f} over {len(v)} cluster(s), below {CL_UNRELIABLE_GAP:g}"]
            for l, v in by.items() if np.median(v) < CL_UNRELIABLE_GAP}


def annotate_report(adata, min_group=MIN_GROUP):
    """QC of the labels: counts per label, markers used, clusters table (cluster method), marker-leakage table,
    unreliable lineages."""
    if "annotation_params" not in adata.uns:
        raise ValueError("run bio.annotate.annotate first")
    P = json.loads(adata.uns["annotation_params"])
    method = P.get("method", "cell")
    clusters = P.get("clusters", [])
    lab = pd.Categorical(adata.obs["annotation"].astype(str), categories=LABELS)
    vc = np.bincount(lab.codes[lab.codes >= 0], minlength=len(LABELS))
    counts = {l: int(c) for l, c in zip(LABELS, vc)}
    fired = np.asarray(adata.obsm["lineage_fired"]).astype(bool)
    real = [l for l in LABELS if l not in ("Ambiguous", "LowCount", "Unclustered") and counts[l] > 0]
    avail = [l for l in P["lineages"] if l not in P["unavailable"]]
    # fraction of cells of each real label (columns) in which each lineage's markers fire (rows): one label-indicator
    # x fired product, no per-label masks
    lab_idx = {l: i for i, l in enumerate(LABELS)}
    ok = lab.codes >= 0
    oh = sparse.csr_matrix((np.ones(int(ok.sum())), (lab.codes[ok], np.flatnonzero(ok))), shape=(len(LABELS), len(lab)))
    frac = (oh @ fired.astype(np.float64)).T / np.maximum(vc, 1)[None, :]         # [lineage, label]
    spec = pd.DataFrame(np.nan, index=avail, columns=real)
    for l in avail:
        li = P["lineages"].index(l)
        for o in real:
            if LABEL_LINEAGE.get(o, o) != l:  # own labels are not a leakage question
                spec.loc[l, o] = frac[li, lab_idx[o]]
    unreliable = {}
    if method in ("cell", "mixed"):
        unreliable.update(_cell_unreliable(spec, real, avail, counts, min_group))
    if method in ("cluster", "mixed"):
        for l, v in _cluster_unreliable(clusters).items():
            unreliable[l] = unreliable.get(l, []) + v
    note = adata.obs["annotation_note"].astype(str)
    cl_df = pd.DataFrame(
        [{**({"section": c["section"]} if "section" in c else {}), "cluster": c["cluster"], "cells": c["n"],
          "label": c["label"], "best": c["best_scaled"],
          "top3": ", ".join(f"{k} {v:.1f}" for k, v in c["top3"]) + ("" if not c["note"] else f"  [{c['note']}]")
                  + ("" if not c.get("t_split") else "  [cells: " + ", ".join(f"{l} {v}" for l, v in c["t_split"].items() if v) + "]")}
         for c in clusters])
    absent = P.get("absent_sections", {l: [1, 1] for l in P.get("absent", [])})
    return {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars), "labels": counts, "method": method,
            "clusters": cl_df, "absent": absent, "cluster_rules": P.get("cluster_rules"),
            "markers_used": P["markers_used"], "markers_requested": P["markers_requested"],
            "unavailable": P["unavailable"], "threshold_per_1000": P["threshold_per_1000"],
            "single_marker": P.get("single_marker", {}), "single_marker_sections": P.get("single_marker_sections"),
            "ambient_ig_excluded": P.get("ambient_ig_excluded", {}),
            "background_source": P["background_source"], "specificity": spec, "unreliable": unreliable,
            "n_no_lineage": int((note == "no lineage fired").sum()),
            "n_tie": int((note == "tie").sum()), "n_no_cluster": int((note == "no cluster id").sum()),
            "t_default_cd4": P["t_default_cd4"], "t_no_subset": P.get("t_no_subset", 0), "n_t_cluster": P.get("n_t_cluster", 0),
            "n_t": P["n_t"], "reference": P.get("reference"), "t_subset_genes_present": P["t_subset_genes_present"],
            "min_group": min_group}


def report_text(rep, name=""):
    """Human-readable version of annotate_report."""
    n = rep["n_cells"]
    method = rep.get("method", "cell")
    if method == "cluster":
        how = ("Labels are cluster-level marker calls: each 10x graph-based cluster is labelled from how specific the "
               "lineage marker sets are in it, and cells inherit the label.")
    elif method == "mixed":
        how = ("Labels are marker calls; some sections were labelled by 10x cluster, others cell by cell (see the "
               "section blocks).")
    else:
        how = "Labels are marker calls (deterministic per-cell scoring, no clustering)."
    L = [f"Annotation report{': ' + name if name else ''}",
         f"{n} cells, {rep['n_genes']} genes. {how} They do not separate malignant from normal epithelium."]
    R = rep.get("cluster_rules") or {}
    if method in ("cluster", "mixed") and R:
        L.append(f"Method: cluster. For each lineage, a cluster's mean per-1000 marker expression is scaled by the lineage's largest "
                 f"cluster mean in the section (1.0 = the cluster where the lineage is highest, however common it is). A cluster is "
                 f"called lineage L when L has the highest scaled value, it is >= {R['min_scaled']:g}, the cluster mean is above the "
                 f"lineage floor, at least {MIN_MARKERS} of L's markers (or the single specific marker) are each >= "
                 f"{R['marker_frac']:g} of their largest cluster mean, and the runner-up scaled value is < {R['ambig']:g} x the best "
                 f"(T/NK and Macrophage/Monocyte/DC do not count against each other); otherwise Ambiguous. A lineage whose largest "
                 f"cluster mean is below {R['absent_mult']:g} x its floor is absent from the section and is not called. "
                 f"The lineage's largest cluster mean (and each marker's) is taken over clusters with at least "
                 f"max({R.get('ref_min_cells', CL_REF_MIN_CELLS)}, {100 * R.get('ref_min_frac', CL_REF_MIN_FRAC):g}% of the "
                 f"section's non-LowCount cells) cells (all clusters if none is that large; a lineage that is below its absent floor in all of them uses all clusters), "
                 f"so a tiny outlier cluster cannot set the scale; smaller clusters are still labelled and can scale above 1.0. Cells under "
                 f"min_counts are left out of the means and labelled LowCount.")
    elif method == "cell":
        L.append("Method: cell (no cluster_10x column, or --clusters none).")
    L += ["", f"Background: {rep['background_source']}", "", "Cells per label:"]
    for l, c in rep["labels"].items():
        L.append(f"  {l:12s} {c:9d}  {100 * c / max(n, 1):5.1f}%")
    if method == "cell":
        L.append(f"  (Ambiguous = {rep['n_no_lineage']} with no lineage firing + {rep['n_tie']} ties between lineages)")
    else:
        L.append(f"  (Ambiguous = {rep['n_no_lineage']} in clusters where no lineage was specific enough or its markers were not "
                 f"enriched + {rep['n_tie']} in clusters tied between lineages)")
        L.append(f"  (Unclustered = {rep['n_no_cluster']} cells ({100 * rep['n_no_cluster'] / max(n, 1):.1f}%) that the clusters file "
                 "does not list; they are not labelled cell by cell, and role 'other' keeps them as cells)")
    cl = rep.get("clusters")
    if method in ("cluster", "mixed") and cl is not None and len(cl):
        L += ["", "Clusters (n cells in cluster, label, best scaled value, top three scaled lineage values; 1.0 = the lineage's strongest cluster):"]
        if len(cl) <= MAX_CLUSTER_ROWS:
            t = cl.copy()
            t["best"] = t["best"].map(lambda v: f"{v:.2f}")
            L.append(t.rename(columns={"top3": "top three scaled lineage values"}).to_string(index=False))
        else:
            L.append(f"  {len(cl)} clusters pooled over sections; see the per-section Clusters tables.")
            L.append("  " + ", ".join(f"{l} {int((cl['label'] == l).sum())}" for l in LABELS if (cl["label"] == l).any())
                     + " clusters")
    ab = rep.get("absent") or {}
    if method in ("cluster", "mixed") and ab:
        L.append("Absent from the section (the lineage's largest cluster mean is below "
                 f"{R.get('absent_mult', CL_ABSENT_MULT):g} x its floor, so it cannot be called there): "
                 + "; ".join(f"{l}" + (f" ({v[0]} of {v[1]} sections)" if v[1] > 1 else "") for l, v in ab.items()) + ".")
    L += ["", "Markers used (in panel / requested; missing ones are listed):"]
    for l, used in rep["markers_used"].items():
        miss = [g for g in rep["markers_requested"][l] if g not in used]
        flag = "  UNAVAILABLE (<2 markers in panel, this lineage cannot be called)" if l in rep["unavailable"] else ""
        if l in rep["single_marker"]:
            flag = f"  SINGLE-MARKER CALL (fires at >={SOLO_MIN_COUNTS} counts of {rep['single_marker'][l]})"
        L.append(f"  {l:11s} {len(used)}/{len(rep['markers_requested'][l])}  {', '.join(used) or '-'}"
                 + (f"   [missing: {', '.join(miss)}]" if miss else "") + flag
                 + f"   threshold {rep['threshold_per_1000'][l]:g}/1000")
    solo = rep["single_marker"]
    if solo:
        ns = rep.get("single_marker_sections")
        L += ["", "Called from a single marker (the panel has no second marker; needs >=2 counts of it, so these calls are "
              "less robust than multi-marker calls): "
              + "; ".join(f"{l} ({g}" + (f", in {ns[l][0]} of {ns[l][1]} sections" if ns and l in ns else "") + ")"
                          for l, g in solo.items())]
    else:
        L += ["", "Lineages called from a single marker: none."]
    for l, gs in rep["ambient_ig_excluded"].items():
        if gs:
            L.append(f"Not used for {l} (ambient immunoglobulin transcripts, present in tumour tissue everywhere): {', '.join(gs)}.")
    L.append("Plasma uses JCHAIN and MZB1 only, and needs both (XBP1 and SDC1 are also expressed by epithelium and other cells).")
    if method == "cell":
        L += ["", "T subsets: CD8T if CD8A+CD8B > CD4; Treg if FOXP3 detected or IL2RA and CTLA4 both detected; else CD4T."]
    else:
        L += ["", "T subsets (" + ("cluster method: the cluster decides T, each cell its subset; cell method in sections without clusters"
                                   if method == "mixed" else "the cluster decides T, each cell inside it decides its subset")
              + "): Treg if FOXP3 count >= 2 and CD8A+CD8B = 0 (panels without FOXP3: IL2RA >= 2 and CTLA4 >= 2, same CD8 exclusion); else CD8T if CD8A+CD8B > CD4; "
                "else CD4T if CD4 > CD8A+CD8B; else T (subset unknown). "
                f"A cluster is NK only if its T value is < {CL_NK_T_RATIO:g}x its NK value and {CL_NK_MIN_SPECIFIC} NK markers "
                "besides GNLY are enriched; otherwise T."]
    have = rep["t_subset_genes_present"]
    L.append(f"  genes in panel: {', '.join(have) or 'none'}")
    n_tc = rep.get("n_t_cluster", 0) if method in ("cluster", "mixed") else 0
    n_tl = rep["n_t"] - n_tc
    if n_tc:
        L.append(f"  {rep['t_no_subset']} of {n_tc} T cells ({100 * rep['t_no_subset'] / n_tc:.0f}%) have no subset evidence "
                 "(no FOXP3, and CD8A+CD8B = CD4, e.g. no counts of either) and are labelled T; treat the CD4T/CD8T/Treg split as soft.")
    if n_tl:
        L.append(f"  {rep['t_default_cd4']} of {n_tl} T cells ({100 * rep['t_default_cd4'] / n_tl:.0f}%) have no CD4, "
                 "CD8A/B, FOXP3 or IL2RA+CTLA4 counts and were labelled CD4T by default; treat the CD4T/CD8T split as soft.")
    spec = rep["specificity"]
    star = f"'*' = above {100 * UNRELIABLE_FRAC:.0f}%; " if method != "cluster" else "no '*' flags here, the cluster method flags lineages by cluster ratio (below); "
    L += ["", f"Marker leakage (per-cell firing of each lineage's markers inside each label): % of cells of each label (columns) "
              f"in which a lineage's markers fire (rows). A lineage's own labels are left blank; {star}labels with fewer than "
              f"{rep['min_group']} cells are shown but not counted."]
    tab = spec.copy().astype(object)
    for l in spec.index:
        for o in spec.columns:
            v = spec.loc[l, o]
            tab.loc[l, o] = "" if np.isnan(v) else f"{100 * v:.1f}" + ("*" if v > UNRELIABLE_FRAC and method != "cluster" else "")
    tab.columns = [f"{c}({rep['labels'][c]})" for c in spec.columns]
    L.append(tab.to_string())
    L.append("")
    if rep["unreliable"]:
        if method == "cluster":
            L.append(f"UNRELIABLE lineages (median of best minus runner-up scaled value over their clusters below {CL_UNRELIABLE_GAP:g}):")
        elif method == "mixed":
            L.append(f"UNRELIABLE lineages (cell method: markers fire in >{100 * UNRELIABLE_FRAC:.0f}% of cells of {UNRELIABLE_LABELS}+ "
                     f"other labels; cluster method: median best-minus-runner-up scaled value below {CL_UNRELIABLE_GAP:g}):")
        else:
            L.append(f"UNRELIABLE lineages (markers fire in >{100 * UNRELIABLE_FRAC:.0f}% of cells of {UNRELIABLE_LABELS}+ other labels):")
        for l, bad in rep["unreliable"].items():
            L.append(f"  {l}: {', '.join(bad)}. Its labels may be mostly background or spill-over; do not trust tests that use it.")
    else:
        L.append("No lineage flagged unreliable.")
    return "\n".join(L) + "\n"


def merge_annotation_params(adatas, names=None):
    """uns['annotation_params'] for a cohort stacked from separately annotated sections, so annotate_report works on
    the pooled object. Markers = union over sections; a lineage is unavailable only if it is in every section;
    thresholds are the mean over sections (each section was scored with its own background). Cluster tables are
    pooled with a `section` field (names, default 0..n-1); method = cluster / cell if all sections agree, else mixed."""
    Ps = [json.loads(a.uns["annotation_params"]) for a in adatas]
    names = [str(x) for x in (names if names is not None else range(len(Ps)))]
    lin = Ps[0]["lineages"]
    methods = {P.get("method", "cell") for P in Ps}
    return json.dumps({
        "method": methods.pop() if len(methods) == 1 else "mixed",
        "clusters": [{**c, "section": nm} for P, nm in zip(Ps, names) for c in P.get("clusters", [])],
        "cluster_rules": next((P["cluster_rules"] for P in Ps if P.get("cluster_rules")), None),
        "absent_sections": {l: [sum(l in P.get("absent", []) for P in Ps), len(Ps)]
                            for l in lin if any(l in P.get("absent", []) for P in Ps)},
        "lineages": lin, "markers_requested": Ps[0]["markers_requested"],
        "markers_used": {l: [g for g in LINEAGE_MARKERS[l] if any(g in P["markers_used"][l] for P in Ps)] for l in lin},
        "unavailable": [l for l in lin if all(l in P["unavailable"] for P in Ps)],
        "single_marker": {l: g for P in Ps for l, g in P.get("single_marker", {}).items()},
        "single_marker_sections": {l: [sum(l in P.get("single_marker", {}) for P in Ps), len(Ps)]
                                   for l in {l for P in Ps for l in P.get("single_marker", {})}},
        "ambient_ig_excluded": {"Plasma": sorted({g for P in Ps for g in P.get("ambient_ig_excluded", {}).get("Plasma", [])})},
        "threshold_per_1000": {l: round(float(np.mean([P["threshold_per_1000"][l] for P in Ps])), 4) for l in lin},
        "background_source": f"per section (own negative-control rate in each of {len(Ps)} sections); "
                             "thresholds below are means over sections",
        "min_counts": Ps[0]["min_counts"], "t_subset_genes_present": sorted({g for P in Ps for g in P["t_subset_genes_present"]}),
        "t_default_cd4": sum(P["t_default_cd4"] for P in Ps), "t_no_subset": sum(P.get("t_no_subset", 0) for P in Ps),
        "n_t_cluster": sum(P.get("n_t_cluster", 0) for P in Ps), "n_t": sum(P["n_t"] for P in Ps),
        "reference": [P.get("reference") for P in Ps if P.get("reference")] or None})
