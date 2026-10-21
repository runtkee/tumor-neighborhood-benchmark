"""Data-level checks on a cell table. Each returns (name, status, detail); status is pass | warn | fail.
A run is refused if any check fails (override with --force, which is recorded in the report)."""
import numpy as np
from scipy.spatial import cKDTree

MIN_ROLE = 300


def run_checks(df, m):
    out = []
    # 1. units: neighboring cells should be ~5-15 µm apart once pixel size is applied
    nn = []
    for _, g in df.groupby("field"):
        if len(g) > 20:
            xy = g[["x_um", "y_um"]].to_numpy()
            nn.append(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
    med = float(np.median(nn)) if nn else np.nan
    out.append(("cell spacing (µm)", "pass" if 5 <= med <= 15 else "fail",
                f"median nearest-cell distance {med:.1f} µm (expected 5–15). Wrong pixel_size_um if outside."))
    # 2. roles present
    counts = df["type"].value_counts().to_dict()
    for r in ["T", "neut", "mac", "B"]:
        n = counts.get(r, 0)
        out.append((f"role {r}", "pass" if n >= MIN_ROLE else ("warn" if n > 0 else "fail" if r in ("T",) else "warn"),
                    f"{n} cells" + ("" if n >= MIN_ROLE else f" (<{MIN_ROLE}: tests involving {r} will be weak or skipped)")))
    unm = sorted(df.loc[df["type"] == "unmapped", "raw_type"].unique())
    out.append(("unmapped labels", "warn" if unm else "pass",
                (f"{len(unm)} labels not mapped (treated as 'other'): {', '.join(unm[:12])}" if unm else "all labels mapped")))
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
