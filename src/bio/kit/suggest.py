"""Draft a manifest from a file by pattern matching (no AI, no network). Every guess is marked '# CHECK'.
The human reviews and edits the YAML; `check` then validates it against the data."""
import re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

COLUMN_HINTS = {
    "patient": r"^(patient|patients|patient_id|sampleid|donor|case)",
    "field": r"^(image|file ?name|fov|roi|region|sample_id|spot|core)",
    "x": r"^(x|pos_x|centroid_?x|x:x|cell x position|centerx)",
    "y": r"^(y|pos_y|centroid_?y|y:y|cell y position|centery)",
    "type": r"(cluster ?name|cell_?type|celltypes|phenotype|annotation|class)",
}
ROLE_HINTS = [  # first match wins
    ("drop", r"dirt|artifact|undefined|unassigned|^nan$"),
    ("B", r"\bb[ _-]?cell|^b$|plasma|cd20"),
    ("T", r"t[ _-]?cell|^cd4|^cd8|treg|cd3\+|^t$|tprolif|cd4t|cd8t"),
    ("neut", r"neutro|granulo|pmn|mdsc"),
    ("mac", r"macro|mac|kupffer|cd68|cd163"),
    ("stroma", r"stroma|fibro|mural|caf|smooth"),
    ("tumor", r"tumou?r|epithel|cancer|carcinoma|hepatocyte"),
]
STATE_HINTS = {"prolif": r"ki-?67|mki67", "PD1": r"pd-?1\b|pdcd1", "GZMB": r"granzyme|gzmb|grzb",
               "LAG3": r"lag-?3", "ICOS": r"\bicos", "apoptosis": r"cleaved.?parp|caspase"}
BLEED_HINTS = {"neut": r"mpo|cd15|cd66b|ceacam8", "mac": r"cd68", "stroma": r"\bsma\b|vimentin|pdgfrb"}


def _first(cols, pattern):
    for c in cols:
        if re.search(pattern, str(c).strip().lower()):
            return c
    return None


def draft(path, dataset):
    if str(path).endswith(".h5ad"):
        df = _h5ad_frame(path)
        cols = list(df.columns)
        guess = {k: _first(cols, p) for k, p in COLUMN_HINTS.items()}
    else:
        cols = list(pd.read_csv(path, nrows=0).columns)
        guess = {k: _first(cols, p) for k, p in COLUMN_HINTS.items()}
        df = pd.read_csv(path, usecols=[v for v in guess.values() if v])  # every row of the key columns
    lines = [f"dataset: {dataset}", "modality: protein  # CHECK (protein | rna)", f"path: {path}", "columns:"]
    for k, v in guess.items():
        lines.append(f"  {k}: {v!r}  # CHECK" if v else f"  {k}: null  # FILL IN")
    px = "1.0"
    if guess["x"] and guess["y"] and guess["field"]:
        nn = []
        for _, g in df.groupby(guess["field"]):
            xy = g[[guess["x"], guess["y"]]].to_numpy(float)
            if len(xy) > 20:
                nn.append(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))
        if nn:
            px = f"{8.0 / np.median(nn):.3f}"
            lines.append(f"pixel_size_um: {px}  # CHECK: estimated assuming ~8 µm between cells; use the instrument's value")
    if px == "1.0":
        lines.append("pixel_size_um: 1.0  # FILL IN")
    lines.append("roles:  # raw label -> T | neut | mac | B | stroma | tumor | other | drop")
    if guess["type"]:
        for lab, n in df[guess["type"]].astype(str).value_counts().items():
            role = next((r for r, p in ROLE_HINTS if re.search(p, lab.lower())), "other")
            lines.append(f"  {lab!r}: {role}  # CHECK ({n} cells)")
    lines.append("states:  # canonical T-state -> column")
    for s, p in STATE_HINTS.items():
        c = _first(cols, p)
        lines.append(f"  {s}: {c!r}  # CHECK" if c else f"  # {s}: not found")
    lines.append("bleed:  # role -> that role's own marker column (used to measure signal bleed into T cells)")
    for r, p in BLEED_HINTS.items():
        c = _first(cols, p)
        lines.append(f"  {r}: {c!r}  # CHECK" if c else f"  # {r}: not found")
    lines.append("# extra_roles: {pmn_mdsc: neut}   # private-only roles and their nearest public proxy")
    lines.append("# extra_states: [lmp]               # private-only T states")
    return "\n".join(lines) + "\n"


def _h5ad_frame(path):
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    return a.obs.reset_index(drop=True)
