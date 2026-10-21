"""Locked analysis recipe (v2 frozen 2026-10-20). Do not tune per dataset; bump the version instead.

v1 -> v2: T-state scores per 1000 counts were diluted by neighbor spillover (fibroblast-contact control
showed the same "effect" as macrophage contact). v2 scores T state as log ratio to T-identity genes.

Input: AnnData with raw counts in .X, obs['annotation'] (cell type), obs['sample'], obsm['spatial'] in µm,
optional obs['region'] (tissue-region label used as a covariate).
Output: one row of metrics per sample (see docs/recipe.md).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree

VERSION = "v2"
NEAR_UM = 10            # "touching"
ZONE_UM = 30            # "nearby" (preprint: 1-3 cell diameters)
TOP_FRAC = 0.10         # activated neutrophils = top 10% activation ratio
MIN_IDENT = 3           # neutrophil identity counts needed to trust its ratio
MIN_CELLS = 300         # minimum cells of each type to run a test
N_SHUFFLE = 50
T_TYPES = ["CD4T", "CD8T", "Tprolif"]
NEUT_ACT = ["AQP9", "BCL2A1", "NCF2", "OSM", "TREM1"]          # Veglia 2021, neutrophil-enriched
NEUT_IDENT = ["CSF3R", "FPR1", "FCGR3B", "CXCR1", "CXCR2"]
T_IDENT = ["CD3D", "CD3E", "CD3G", "CD2", "CD247"]
STATES = {
    "prolif": ["MKI67", "TOP2A", "CDK1", "PCNA", "MCM2", "TYMS"],
    "exh": ["PDCD1", "LAG3", "TIGIT", "CTLA4", "TOX"],
    "cyto": ["GZMB", "GZMA", "GZMK", "PRF1", "IFNG"],
}


class _Cells:
    def __init__(self, a):
        self.a = a
        self.X = a.X.tocsc()
        self.tot = np.maximum(np.asarray(self.X.sum(1)).ravel(), 1)
        self.lab = a.obs["annotation"].astype(str).to_numpy()
        self.xy = np.asarray(a.obsm["spatial"])
        self.region = a.obs["region"].astype(str).to_numpy() if "region" in a.obs else np.full(a.n_obs, "all")
        self.rng = np.random.default_rng(0)

    def counts(self, genes):
        idx = [self.a.var_names.get_loc(g) for g in genes if g in self.a.var_names]
        return np.asarray(self.X[:, idx].sum(1)).ravel() if idx else np.zeros(self.a.n_obs)

    def per1k(self, genes):
        return 1000 * self.counts(genes) / self.tot

    def nearest(self, src, tgt):
        d, _ = cKDTree(self.xy[tgt]).query(self.xy[src], k=1)
        return d

    def n_within(self, centers, tgt, r):
        tree = cKDTree(self.xy[tgt])
        return np.array([len(x) for x in tree.query_ball_point(self.xy[centers], r)])


def _depletion(c, target, source, pool):
    """Share of target cells within ZONE_UM of a source cell, real vs source labels shuffled within pool.
    Returns log2(real / null mean): <0 = target avoids source."""
    real = np.mean(c.nearest(target, source) < ZONE_UM)
    idx, k = np.where(pool)[0], source.sum()
    null = []
    for _ in range(N_SHUFFLE):
        m = np.zeros(len(c.lab), bool); m[c.rng.choice(idx, k, replace=False)] = True
        null.append(np.mean(c.nearest(target, m) < ZONE_UM))
    return np.log2(real / np.mean(null)), real, np.mean(null)


def _contact_effects(c, target, source, prefix):
    """Adjusted difference in T state for target cells <NEAR_UM from a source cell vs >= ZONE_UM."""
    d = c.nearest(target, source)
    near_t = c.n_within(target, target, 50)
    tid = c.counts(T_IDENT)[target]
    df = pd.DataFrame({k: np.log((c.counts(g)[target] + 0.5) / (tid + 0.5)) for k, g in STATES.items()})
    keep = tid >= MIN_IDENT  # need enough T-identity counts for a stable ratio
    df, d, near_t = df[keep].copy(), d[keep], near_t[keep]
    target = np.where(target)[0][keep]
    df["near"] = (d < NEAR_UM).astype(int)
    df["mid"] = ((d >= NEAR_UM) & (d < ZONE_UM)).astype(int)
    df["region"] = c.region[target]
    df["tdens"] = np.log1p(near_t)
    df["cd8"] = (c.lab[target] == "CD8T").astype(int)
    df["ltot"] = np.log(c.tot[target])
    out = {}
    reg = " + C(region)" if df.region.nunique() > 1 else ""
    for k in STATES:
        m = smf.ols(f"{k} ~ near + mid{reg} + tdens + cd8 + ltot", df).fit(cov_type="HC1")
        # log-ratio scale: exp(beta)-1 = % change in state genes relative to T-identity genes
        out[f"{prefix}_{k}_pct"] = 100 * (np.exp(m.params["near"]) - 1)
        out[f"{prefix}_{k}_se_pct"] = 100 * m.bse["near"]
    return out


def _activated_neutrophils(c, T, N_all):
    act, ident = c.counts(NEUT_ACT), c.counts(NEUT_IDENT)
    N = N_all & (ident >= MIN_IDENT)
    if N.sum() < MIN_CELLS:
        return {"actN_n": int(N.sum())}
    ratio = np.log((act + 0.5) / (ident + 0.5))
    z = (ratio[N] - ratio[N].mean()) / ratio[N].std()
    cov = {"T30": c.n_within(N, T, ZONE_UM), "all30": c.n_within(N, np.ones(len(c.lab), bool), ZONE_UM)}
    for t in ["Epithelial", "Neutrophil", "Macrophage", "Fibroblast"]:
        cov[f"{t}30"] = c.n_within(N, c.lab == t, ZONE_UM)
    df = pd.DataFrame({"z": z, "region": c.region[N], **cov})
    reg = " + C(region)" if df.region.nunique() > 1 else ""
    out = {"actN_n": int(N.sum())}
    for name, f in [("raw", "T30 ~ z"), ("adj", f"T30 ~ z{reg} + all30 + Epithelial30 + Neutrophil30 + Macrophage30 + Fibroblast30")]:
        m = smf.glm(f, df, family=sm.families.NegativeBinomial()).fit()
        out[f"actN_T30_{name}_log2_per_sd"] = m.params["z"] / np.log(2)
        out[f"actN_T30_{name}_se"] = m.bse["z"] / np.log(2)
    hi = ratio[N] > np.quantile(ratio[N], 1 - TOP_FRAC)
    out["actN_top10_T30_mean"] = df.T30[hi].mean()
    out["actN_rest_T30_mean"] = df.T30[~hi].mean()
    return out


def run_sample(a, sample):
    c = _Cells(a)
    lab = c.lab
    T = np.isin(lab, T_TYPES)
    N = lab == "Neutrophil"
    M = lab == "Macrophage"
    B = np.isin(lab, ["B", "Plasma"])
    F = lab == "Fibroblast"
    row = {"sample": sample, "recipe": VERSION, "n_cells": a.n_obs, "n_T": int(T.sum()), "n_neut": int(N.sum()),
           "n_mac": int(M.sum()), "n_B": int(B.sum())}
    if T.sum() >= MIN_CELLS and N.sum() >= MIN_CELLS:
        row["T_near_neut_log2"], row["T_near_neut_real"], row["T_near_neut_null"] = _depletion(c, T, N, ~T)
        row.update(_contact_effects(c, T, N, "neutTouch"))
        row.update(_activated_neutrophils(c, T, N))
        if B.sum() >= MIN_CELLS:  # control: do B/plasma cells avoid neutrophils the same way?
            row["CTRL_B_near_neut_log2"] = _depletion(c, B, N, ~B)[0]
    if T.sum() >= MIN_CELLS and M.sum() >= MIN_CELLS:
        row["T_near_mac_log2"] = _depletion(c, T, M, ~T)[0]
        row.update(_contact_effects(c, T, M, "macTouch"))
    if T.sum() >= MIN_CELLS and F.sum() >= MIN_CELLS:  # control: fibroblast contact
        row.update(_contact_effects(c, T, F, "CTRL_fibTouch"))
    return row
