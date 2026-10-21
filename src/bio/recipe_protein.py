"""Locked recipe, protein-imaging version (v2p.1, frozen 2026-10-20). Same distances and logic as recipe v2.

v2p -> v2p.1: shuffles use the shared PRNG (bio.prng) so the browser version reproduces them exactly.
Method unchanged; only the random draws differ.

Per field (image): nearest-neighbor distances, shuffle null, contact effects on T-state proteins.
Protein state = within-field percentile rank (0-100), so brightness differences between fields cancel.
Per patient: pool fields -> one row. Across patients: mean, SE, sign count (see meta()).
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

from bio.prng import Mulberry32, choose, fnv1a

VERSION = "v2p.1"
NEAR_UM, ZONE_UM, N_SHUFFLE, MIN_PER_PATIENT = 10, 30, 20, 20


def _nearest(xy, src, tgt):
    if tgt.sum() == 0:
        return np.full(src.sum(), np.inf)
    return cKDTree(xy[tgt]).query(xy[src], k=1)[0]


PAIRS = [("T", "neut"), ("B", "neut"), ("T", "mac")]
SOURCES = ["neut", "mac", "stroma"]


def patient_rows(df, types, states, bleed, pairs=PAIRS, sources=SOURCES):
    """df: one row per cell with columns patient, field, x_um, y_um, type, and protein columns.
    types: dict role -> list of type labels (roles: T, neut, mac, B, stroma).
    states: list of protein columns describing T state. bleed: dict source role -> protein column."""
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
    return pd.DataFrame(rows)


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
