"""Data-level checks on a cell table. Each returns (name, status, detail); status is pass | warn | fail.
A run is refused if any check fails (override with --force, which is recorded in the report)."""
import numpy as np
from scipy.spatial import cKDTree

from bio.compartment import compartments_for

from . import SUBS

MIN_ROLE = 300
MIN_NEST_FRAC, MIN_NEST_T = 0.20, 300
MIN_SPACING, MAX_SPACING = 4.0, 15.0  # µm; real IMC (CellProfiler segmentation, 1 µm/px) has a median of 4.9


def run_checks(df, m):
    out = []
    # 1. units: neighboring cells should be ~4-15 µm apart once pixel size is applied
    nn = []
    for _, g in df.groupby("field"):
        if len(g) > 20:
            xy = g[["x_um", "y_um"]].to_numpy()
            nn.append(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
    med = float(np.median(nn)) if nn else np.nan
    out.append(("cell spacing (µm)", "pass" if MIN_SPACING <= med <= MAX_SPACING else "fail",
                f"median nearest-cell distance {med:.1f} µm (expected {MIN_SPACING:g}–{MAX_SPACING:g}). Wrong pixel_size_um if outside."))
    # 2. roles present
    counts = df["type"].value_counts().to_dict()
    for r in ["T", "neut", "mac", "B"]:
        n = counts.get(r, 0)
        out.append((f"role {r}", "pass" if n >= MIN_ROLE else ("warn" if n > 0 else "fail" if r in ("T",) else "warn"),
                    f"{n} cells" + ("" if n >= MIN_ROLE else f" (<{MIN_ROLE}: tests involving {r} will be weak or skipped)")))
    unm = sorted(df.loc[df["type"] == "unmapped", "raw_type"].unique())
    out.append(("unmapped labels", "warn" if unm else "pass",
                (f"{len(unm)} labels not mapped (treated as 'other'): {', '.join(unm[:12])}" if unm else "all labels mapped")))
    # 2b. tumor cells and compartments (recipe v3p)
    n_tumor = counts.get("tumor", 0)
    if n_tumor == 0:
        out.append(("tumor cells", "warn", "no label maps to the tumor role: compartments are unavailable, only whole-image results are computed"))
    else:
        explicit = "comp" in df.columns
        d = df if explicit else compartments_for(df[["patient", "field", "x_um", "y_um", "type"]].copy())
        unk = d.groupby(["patient", "field"])["comp"].apply(lambda c: bool((c == "unknown").all()))
        out.append(("tumor cells", "pass",
                    f"{n_tumor} tumor cells; compartments from " +
                    ("the annotation column" if explicit else "geometry (nest = >=50% tumor cells within 50 µm; margin = within 100 µm of a nest cell)") +
                    (f"; {int(unk.sum())} of {len(unk)} images have no usable compartment and are skipped by compartment metrics" if unk.any() else "")))
        known = d[d["comp"] != "unknown"]
        if len(known):
            frac = float((known["comp"] == "nest").mean())
            n_t_nest = int(((known["comp"] == "nest") & (known["type"] == "T")).sum())
            bad = frac < MIN_NEST_FRAC or n_t_nest < MIN_NEST_T
            out.append(("nest size", "warn" if bad else "pass",
                        f"nest holds {100 * frac:.0f}% of cells and {n_t_nest} T cells" +
                        (f" (warn below {100 * MIN_NEST_FRAC:.0f}% of cells or {MIN_NEST_T} T cells: nest tests will be weak)" if bad else "")))
        else:
            out.append(("nest size", "warn", "no image has a compartment: nest tests are skipped"))
    # 2c. T subsets
    if "sub" in df.columns and counts.get("T", 0) > 0:
        tsub = df.loc[df["type"] == "T", "sub"].value_counts().to_dict()
        have = [f"{s} ({tsub[s]})" for s in SUBS if tsub.get(s, 0) > 0]
        out.append(("T subsets", "pass" if have else "warn",
                    ("present: " + ", ".join(have) + f"; T cells without subset: {tsub.get('T', 0)}") if have
                    else "no T subset in the labels: subset tests are skipped"))
    # 3. each role's own marker should be high in that role, compared with T cells
    for col in [c for c in df.columns if c.startswith("bleed_")]:
        r = col[len("bleed_"):]
        if counts.get(r, 0) == 0 or counts.get("T", 0) == 0:
            continue
        pct = df.groupby("field")[col].rank(pct=True) * 100
        gap = pct[df["type"] == r].median() - pct[df["type"] == "T"].median()
        out.append((f"marker sanity: {r}", "pass" if gap > 20 else "fail",
                    f"{col} median percentile in {r} minus in T = {gap:.0f} (expect >20; low = mislabeled role or wrong marker)"))
    # 4. states available
    have = [s for s in m.get("states", {}) if s in df.columns and df[s].notna().any()]
    out.append(("T-state markers", "pass" if have else "warn", f"available: {', '.join(have) or 'none'}"))
    # 5. design
    pats = df["patient"].nunique()
    out.append(("patients", "pass" if pats >= 5 else "warn",
                f"{pats} patients, {df['field'].nunique()} fields (<5 patients: cross-patient statistics unreliable)"))
    return out


def blocking(checks):
    return [c for c in checks if c[1] == "fail"]
