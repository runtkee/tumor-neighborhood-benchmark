"""Run the locked protein-style recipe (v2p) on a kit cell table -> one row per patient."""
from bio.recipe_protein import PAIRS, SOURCES, VERSION, patient_rows

from . import ROLES


def compute(df, m):
    roles_present = set(df["type"].unique())
    types = {r: [r] for r in roles_present if r in ROLES or r in m.get("extra_roles", {})}
    states = [s for s in m.get("states", {}) if s in df.columns and df[s].notna().any()]
    bleed = {r: f"bleed_{r}" for r in m.get("bleed", {}) if f"bleed_{r}" in df.columns}
    extra = list(m.get("extra_roles", {}))
    pairs = PAIRS + [("T", r) for r in extra]
    sources = SOURCES + extra
    rows = patient_rows(df, types, states, bleed, pairs=pairs, sources=sources)
    rows.insert(0, "dataset", m["dataset"])
    rows["recipe"] = VERSION
    return rows
