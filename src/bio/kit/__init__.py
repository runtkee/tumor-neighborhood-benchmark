"""Local benchmark kit: run the locked recipe on a private dataset and compare it with the bundled public
reference, without any data leaving the machine. See docs/kit.md."""
ROLES = ["T", "neut", "mac", "B", "stroma", "tumor", "other"]
STATES = ["prolif", "PD1", "GZMB", "LAG3", "ICOS", "apoptosis"]
SUBS = ["CD8", "CD4", "Treg"]  # a role value may carry a T subset: "T:CD8" -> type T, sub CD8
MAC_STATES = ["PDL1", "IDO1", "CD163", "HLADR", "ARG1", "M2"]
DEFAULT_MAC_IDENTITY_GENES = ["CD68", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC"]
COMPARTMENT_LABELS = ["nest", "margin", "far"]
