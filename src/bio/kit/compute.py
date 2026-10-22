"""Run the locked protein-style recipe (v3p) on a kit cell table -> one row per patient."""
from bio.recipe_protein import PAIRS, SOURCES, VERSION, patient_rows

from . import ROLES


def compute(df, m):
    roles_present = set(df["type"].unique())
    types = {r: [r] for r in roles_present if r in ROLES or r in m.get("extra_roles", {})}
    # a readout is usable only with >= 2 distinct numeric values (same rule as web/recipe.js usable())
    usable = lambda c: c in df.columns and df[c].nunique(dropna=True) >= 2
    states = [s for s in m.get("states", {}) if usable(s)]
    bleed = {r: f"bleed_{r}" for r in m.get("bleed", {}) if usable(f"bleed_{r}")}
    mac_states = {s: f"macstate_{s}" for s in m.get("mac_states", {}) if usable(f"macstate_{s}")}
    extra = list(m.get("extra_roles", {}))
    pairs = PAIRS + [("T", r) for r in extra]
    sources = SOURCES + extra
    rows = patient_rows(df, types, states, bleed, pairs=pairs, sources=sources, mac_states=mac_states)
    rows.insert(0, "dataset", m["dataset"])
    rows["recipe"] = VERSION
    return rows
