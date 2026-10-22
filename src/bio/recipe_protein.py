"""Locked recipe, protein-imaging version v3p.1 (specification: docs/recipe-v3.md).

v3p = v2p.1 (frozen 2026-10-20, every v2p.1 column keeps its name and exact value) plus tumor compartments,
composition, T-cell subsets, nest-restricted proximity, distance gradients and macrophage-state tests.
The whole-image v2p.1 loop below is unchanged; the v3p outputs are computed in separate passes afterwards.
v3p.1 (2026-10-21): the nest null pool excludes tumor cells; sparse contact tests also report pooled sums
(_near_sum, _far_n, _far_sum); TmacTouch_* is whole-image, nest_TmacTouch_* is the nest-only version.

v2p -> v2p.1: shuffles use the shared PRNG (bio.prng) so the browser version reproduces them exactly.
Method unchanged; only the random draws differ.

Per field (image): nearest-neighbor distances, shuffle null, contact effects on T-state proteins.
Protein state = within-field percentile rank (0-100), so brightness differences between fields cancel.
Per patient: pool fields -> one row. Across patients: mean, SE, sign count (see meta()).
"""
import math

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

from bio.compartment import COMPARTMENTS, ORDER as COMP_ORDER, compartments_for
from bio.prng import Mulberry32, choose, fnv1a

VERSION = "v3p.1"
NEAR_UM, ZONE_UM, N_SHUFFLE, MIN_PER_PATIENT = 10, 30, 20, 20
SUBS = ["CD8", "CD4", "Treg"]                    # T-cell subsets (df column "sub"; "T" = subset unknown)
IMMUNE = ["T", "neut", "mac", "B"]               # denominator of the immune-share columns
GRAD_EDGES = [0, 5, 10, 20, 30, 50, 100]         # um; bin i = [edge i, edge i+1), last bin = 100+
GRAD_LABELS = ["0-5", "5-10", "10-20", "20-30", "30-50", "50-100", "100+"]
GRAD_SOURCES = ["neut", "mac"]
SUM_KEYS = ["n", "near_sum", "far_n", "far_sum", "diff"]  # v3p.1 column suffixes of the nest and macrophage-state contact tests


def _nearest(xy, src, tgt):
    if tgt.sum() == 0:
        return np.full(src.sum(), np.inf)
    return cKDTree(xy[tgt]).query(xy[src], k=1)[0]


PAIRS = [("T", "neut"), ("B", "neut"), ("T", "mac")]
SOURCES = ["neut", "mac", "stroma"]


def patient_rows(df, types, states, bleed, pairs=PAIRS, sources=SOURCES, mac_states=None):
    """df: one row per cell with columns patient, field, x_um, y_um, type, and protein columns.
    Optional df columns: sub (T subset CD8/CD4/Treg/T, default T) and comp (nest/margin/far/unknown; default: computed
    from geometry by bio.compartment).
    types: dict role -> list of type labels (roles: T, neut, mac, B, stroma).
    states: list of protein columns describing T state. bleed: dict source role -> protein column.
    mac_states: dict macrophage-state name -> df column (protein or gene-score value per cell), default none."""
    df = df.reset_index(drop=True)
    for col in states + list(bleed.values()):
        df[col + "_pct"] = df.groupby(["patient", "field"])[col].rank(pct=True) * 100  # image = patient + image name
    role = {r: df["type"].isin(v).to_numpy() for r, v in types.items()}
    acc = {}
    for (pat, field), g in df.groupby(["patient", "field"]):
        ix = g.index.to_numpy()
        xy = g[["x_um", "y_um"]].to_numpy()
        R = {r: m[ix] for r, m in role.items()}
        a = acc.setdefault(pat, {})
        for tgt, src in pairs:
            if tgt not in R or src not in R or R[tgt].sum() == 0 or R[src].sum() == 0:
                continue
            real = (_nearest(xy, R[tgt], R[src]) < ZONE_UM).sum()
            pool = np.where(~R[tgt])[0]
            rng = Mulberry32(fnv1a(f"{pat}|{field}|{tgt}|{src}"))
            null = []
            for _ in range(N_SHUFFLE):
                m = np.zeros(len(ix), bool); m[choose(rng, pool, int(R[src].sum()))] = True
                null.append((_nearest(xy, R[tgt], m) < ZONE_UM).sum())
            k = f"{tgt}_near_{src}"
            a[k + "_real"] = a.get(k + "_real", 0) + real
            a[k + "_null"] = a.get(k + "_null", 0) + np.mean(null)
        for src in sources:
            if src not in R or R["T"].sum() == 0 or R[src].sum() == 0:
                continue
            d = _nearest(xy, R["T"], R[src])
            t = g[R["T"]]
            near, far = d < NEAR_UM, d >= ZONE_UM
            for col in states + ([bleed[src]] if src in bleed else []):
                v = t[col + "_pct"].to_numpy()
                a.setdefault(f"{src}Touch_{col}_near", []).extend(v[near])
                a.setdefault(f"{src}Touch_{col}_far", []).extend(v[far])
    rows = []
    for pat, a in acc.items():
        r = {"patient": pat, "recipe": VERSION}
        for k in list(a):
            if k.endswith("_real"):
                base = k[:-5]
                r[f"{base}_log2"] = np.log2((a[k] + 1) / (a[base + "_null"] + 1))
            if k.endswith("_near"):
                base = k[:-5]
                n, f = np.array(a[k]), np.array(a[base + "_far"])
                r[f"{base}_n"] = len(n)
                r[f"{base}_diff"] = n.mean() - f.mean() if len(n) >= MIN_PER_PATIENT and len(f) >= MIN_PER_PATIENT else np.nan
        rows.append(r)
    legacy = pd.DataFrame(rows)  # v2p.1 columns, in their original order
    new = _v3p(df, role, states, bleed, pairs, sources, mac_states or {}, list(acc))
    return pd.concat([legacy, new], axis=1)


# ---------------------------------------------------------------------------------------------------------------
# v3p passes (docs/recipe-v3.md). Plain loops over images; each image is patient + image name.
# ---------------------------------------------------------------------------------------------------------------

def grad_bin(d):
    """Bin index 0..6 for distances d (um): [0,5) [5,10) [10,20) [20,30) [30,50) [50,100) [100,inf)."""
    return np.searchsorted(GRAD_EDGES, d, side="right") - 1


def _stratified_mask(n, pool_by_comp, k_by_comp, rng):
    """Boolean mask of n cells with k_by_comp[c] cells drawn from pool_by_comp[c] for each compartment c, drawing the
    compartments in the fixed order nest, margin, far, unknown, all from the one rng."""
    m = np.zeros(n, bool)
    for c in COMP_ORDER:
        k = min(k_by_comp[c], len(pool_by_comp[c]))
        if k > 0:
            m[choose(rng, pool_by_comp[c], k)] = True
    return m


def _v3p(df, role, states, bleed, pairs, sources, mac_states, patients):
    """Returns one row per patient (in the order of `patients`) with the v3p columns:
    composition, nest-restricted proximity, distance gradient, nest contact tests, macrophage state next to T cells.
    df has the _pct columns of states/bleed already; role: dict role -> boolean array over df rows."""
    n_all = len(df)
    for col in mac_states.values():
        df[col + "_pct"] = df.groupby(["patient", "field"])[col].rank(pct=True) * 100
    pct = {c: df[c + "_pct"].to_numpy() for c in states + list(bleed.values()) + list(mac_states.values())}
    sub_all = df["sub"].astype(str).to_numpy() if "sub" in df.columns else np.full(n_all, "T")
    compartments_for(df)  # keeps an explicit df["comp"], otherwise geometry
    comp_all = df["comp"].astype(str).to_numpy()
    have_comp = bool((comp_all != "unknown").any())
    have_nest = bool((comp_all == "nest").any())
    is_T = role["T"] if "T" in role else np.zeros(n_all, bool)
    is_tum = (df["type"] == "tumor").to_numpy()
    subs = [s for s in SUBS if (is_T & (sub_all == s)).any()]
    comp_roles = [r for r in IMMUNE if r in role]
    names = [x for x in ["T"] + subs + ["neut", "mac", "B"] if x in subs or x in comp_roles]
    grad_pairs = [("T", src) for src in GRAD_SOURCES if src in role]
    nest_pairs = list(pairs) + [p for p in grad_pairs if p not in pairs] + [(s, src) for s in subs for src in GRAD_SOURCES if src in role]

    acc = {pat: {} for pat in patients}
    cn = {pat: {} for pat in patients}  # composition counts: (compartment, key) -> count
    for (pat, field), g in df.groupby(["patient", "field"]):
        ix = g.index.to_numpy()
        xy = g[["x_um", "y_um"]].to_numpy()
        R = {r: m[ix] for r, m in role.items()}
        cm = comp_all[ix]
        sub = sub_all[ix]
        a, counts = acc[pat], cn[pat]
        n = len(ix)

        # 1. composition (pooled counts; percent is taken after pooling over the patient's images)
        if have_comp:
            for c in COMPARTMENTS:
                inc = cm == c
                counts[(c, "cells")] = counts.get((c, "cells"), 0) + int(inc.sum())
                for r in comp_roles:
                    counts[(c, r)] = counts.get((c, r), 0) + int((R[r] & inc).sum())
                for s in subs:
                    counts[(c, s)] = counts.get((c, s), 0) + int((R["T"] & inc & (sub == s)).sum())

        nest = cm == "nest"
        if have_comp and nest.any() and "T" in R:
            # 2 + 3. nest-restricted proximity and distance gradient, stratified shuffle null
            tmask = dict(R)
            for s in subs:
                tmask[s] = R["T"] & (sub == s)
            for tgt, src in nest_pairs:
                if tgt not in tmask or src not in R:
                    continue
                Tn = tmask[tgt] & nest  # measured cells: this type inside the nest
                if Tn.sum() == 0 or R[src].sum() == 0:
                    continue
                d = _nearest(xy, Tn, R[src])
                real = (d < ZONE_UM).sum()
                # v3p.1: the pool holds non-target, non-tumor cells (immune and stromal cells share the non-epithelial
                # strands inside a nest; tumor cells are never a source label in these tests)
                pool_by_comp = {c: np.where(~tmask[tgt] & ~is_tum[ix] & (cm == c))[0] for c in COMP_ORDER}
                k_by_comp = {c: int((R[src] & (cm == c)).sum()) for c in COMP_ORDER}
                rng = Mulberry32(fnv1a(f"{pat}|{field}|{tgt}|{src}|nest"))
                null = []
                null_bins = np.zeros(len(GRAD_LABELS))
                for _ in range(N_SHUFFLE):
                    m = _stratified_mask(n, pool_by_comp, k_by_comp, rng)
                    dn = _nearest(xy, Tn, m)
                    null.append((dn < ZONE_UM).sum())
                    null_bins += np.bincount(grad_bin(dn), minlength=len(GRAD_LABELS))
                k = f"nest_{tgt}_near_{src}"
                a[k + "_real"] = a.get(k + "_real", 0) + real
                a[k + "_null"] = a.get(k + "_null", 0) + np.mean(null)
                if (tgt, src) in grad_pairs:  # same shuffles as the line above
                    bi = grad_bin(d)
                    real_bins = np.bincount(bi, minlength=len(GRAD_LABELS))
                    for i, lab in enumerate(GRAD_LABELS):
                        gk = f"grad_{src}_{lab}"
                        a[gk + "_real"] = a.get(gk + "_real", 0) + int(real_bins[i])
                        a[gk + "_null"] = a.get(gk + "_null", 0) + null_bins[i] / N_SHUFFLE
                        for col in states:
                            a.setdefault(f"{gk}_{col}_mean", []).extend(pct[col][ix][Tn][bi == i])

            # 4. contact tests on nest T cells (and per subset)
            Tn = R["T"] & nest
            sub_n = sub[Tn]
            for src in sources:
                if src not in R or Tn.sum() == 0 or R[src].sum() == 0:
                    continue
                d = _nearest(xy, Tn, R[src])
                near, far = d < NEAR_UM, d >= ZONE_UM
                for col in states + ([bleed[src]] if src in bleed else []):
                    v = pct[col][ix][Tn]
                    a.setdefault(f"nest_{src}Touch_{col}_near", []).extend(v[near])
                    a.setdefault(f"nest_{src}Touch_{col}_far", []).extend(v[far])
                for s in subs:
                    sm = sub_n == s
                    for col in states:
                        v = pct[col][ix][Tn]
                        a.setdefault(f"nest_{src}Touch_{s}_{col}_near", []).extend(v[near & sm])
                        a.setdefault(f"nest_{src}Touch_{s}_{col}_far", []).extend(v[far & sm])

        # 5. macrophage state next to T cells: all macrophages of the image (TmacTouch_*), and nest macrophages only
        #    (nest_TmacTouch_*, only when the dataset has a nest)
        if mac_states and "mac" in R and "T" in R and R["T"].sum() > 0 and R["mac"].sum() > 0:
            Mm = R["mac"]
            d = _nearest(xy, Mm, R["T"])
            near, far = d < NEAR_UM, d >= ZONE_UM
            inn = nest[Mm]  # nest flag of each macrophage, same order as d
            for name, col in mac_states.items():
                v = pct[col][ix][Mm]
                a.setdefault(f"TmacTouch_{name}_near", []).extend(v[near])
                a.setdefault(f"TmacTouch_{name}_far", []).extend(v[far])
                if have_nest and inn.any():
                    a.setdefault(f"nest_TmacTouch_{name}_near", []).extend(v[near & inn])
                    a.setdefault(f"nest_TmacTouch_{name}_far", []).extend(v[far & inn])

    rows = []
    for pat in patients:
        a, counts = acc[pat], cn[pat]
        r = {}
        if have_comp:
            for c in COMPARTMENTS:
                cells = counts.get((c, "cells"), 0)
                imm = sum(counts.get((c, x), 0) for x in comp_roles)
                r[f"comp_{c}_cells"] = cells
                for x in names:
                    r[f"comp_{c}_pct_{x}"] = 100 * counts.get((c, x), 0) / cells if cells > 0 else np.nan
                for x in names:
                    r[f"comp_{c}_immshare_{x}"] = 100 * counts.get((c, x), 0) / imm if imm > 0 else np.nan
        for k in list(a):
            if k.endswith("_real"):
                base = k[:-5]
                r[f"{base}_log2"] = np.log2((a[k] + 1) / (a[base + "_null"] + 1))
            if k.endswith("_mean"):
                v = np.array(a[k])
                r[k[:-5]] = v.mean() if len(v) >= MIN_PER_PATIENT else np.nan
            if k.endswith("_near"):
                base = k[:-5]
                nr, f = np.array(a[k]), np.array(a[base + "_far"])
                r[f"{base}_n"] = len(nr)
                # v3p.1: pooled sums for the sparse tests (every key in this pass is a nest or macrophage-state key)
                r[f"{base}_near_sum"] = math.fsum(nr)
                r[f"{base}_far_n"] = len(f)
                r[f"{base}_far_sum"] = math.fsum(f)
                r[f"{base}_diff"] = nr.mean() - f.mean() if len(nr) >= MIN_PER_PATIENT and len(f) >= MIN_PER_PATIENT else np.nan
        rows.append(r)
    out = pd.DataFrame(rows, index=range(len(patients)))
    # fixed column order (independent of which patient has what): composition, nest proximity, gradient, nest contact, macrophage
    canon = []
    if have_comp:
        for c in COMPARTMENTS:
            canon += [f"comp_{c}_cells"] + [f"comp_{c}_pct_{x}" for x in names] + [f"comp_{c}_immshare_{x}" for x in names]
        canon += [f"nest_{t}_near_{s}_log2" for t, s in nest_pairs]
        for src in GRAD_SOURCES:
            for lab in GRAD_LABELS:
                canon += [f"grad_{src}_{lab}_log2"] + [f"grad_{src}_{lab}_{col}" for col in states]
        for src in sources:
            cols = states + ([bleed[src]] if src in bleed else [])
            canon += [f"nest_{src}Touch_{col}_{k}" for col in cols for k in SUM_KEYS]
            canon += [f"nest_{src}Touch_{s}_{col}_{k}" for s in subs for col in states for k in SUM_KEYS]
    canon += [f"TmacTouch_{name}_{k}" for name in mac_states for k in SUM_KEYS]
    canon += [f"nest_TmacTouch_{name}_{k}" for name in mac_states for k in SUM_KEYS]
    ordered = [c for c in canon if c in out.columns]
    return out[ordered + [c for c in out.columns if c not in ordered]]


def meta(rows):
    """Across patients: mean, SE, number of patients with effect < 0, Wilcoxon p (null: median 0)."""
    out = []
    for c in rows.columns:
        if not (c.endswith("_log2") or c.endswith("_diff")):
            continue
        v = rows[c].dropna().to_numpy()
        if len(v) < 5:
            continue
        p = wilcoxon(v).pvalue if np.any(v != 0) else np.nan
        out.append(dict(metric=c, n_patients=len(v), mean=v.mean(), se=v.std(ddof=1) / np.sqrt(len(v)),
                        n_negative=int((v < 0).sum()), wilcoxon_p=p))
    return pd.DataFrame(out)
