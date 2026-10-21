"""Dataset manifest: the human-reviewed mapping from a file's columns and labels onto the kit's roles.

Example (protein, CSV):
    dataset: CODEX_CRC_Schurch2020
    modality: protein                 # protein | rna
    path: data/codex_crc/CRC_clusters_neighborhoods_markers.csv
    columns: {patient: patients, field: File Name, x: "X:X", y: "Y:Y", type: ClusterName}
    pixel_size_um: 0.377
    roles: {"CD8+ T cells": T, granulocytes: neut, dirt: drop, ...}   # raw label -> role ("drop" removes)
    states: {prolif: "Ki67 - proliferation:Cyc_5_ch_4", ...}          # canonical state -> column
    bleed: {neut: "CD15 - granulocytes:Cyc_14_ch_2", ...}             # role -> that role's own marker
    extra_roles: {pmn_mdsc: neut}     # private-only roles -> nearest public proxy role (optional)
For modality rna (h5ad): columns.x / columns.y may be "obsm:spatial:0" / "obsm:spatial:1", states and bleed
take gene lists, and identity_genes gives the T-identity genes that RNA states are expressed relative to.
"""
from pathlib import Path

import yaml

from . import ROLES, STATES

REQUIRED = ["dataset", "modality", "path", "columns", "pixel_size_um", "roles"]


def load(path):
    m = yaml.safe_load(Path(path).read_text())
    m["_manifest_path"] = str(path)
    return m


def problems(m):
    """Structural problems that block a run (data-level checks live in checks.py)."""
    out = [f"missing key: {k}" for k in REQUIRED if k not in m]
    if out:
        return out
    if m["modality"] not in ("protein", "rna"):
        out.append("modality must be protein or rna")
    for k in ["patient", "field", "x", "y", "type"]:
        if k not in m["columns"]:
            out.append(f"columns.{k} missing")
    allowed = set(ROLES) | {"drop"} | set(m.get("extra_roles", {}))
    bad = sorted({r for r in m["roles"].values() if r not in allowed})
    if bad:
        out.append(f"unknown roles: {bad} (allowed: {sorted(allowed)})")
    for s in m.get("states", {}):
        if s not in STATES and s not in m.get("extra_states", []):
            out.append(f"unknown state {s!r}; add it to extra_states to mark it private-only")
    for r, proxy in m.get("extra_roles", {}).items():
        if proxy not in ROLES:
            out.append(f"extra role {r}: proxy {proxy!r} is not a public role")
    return out
