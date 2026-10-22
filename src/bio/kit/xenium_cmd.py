"""`python -m bio.kit xenium <bundle> [<bundle> ...] --dataset NAME --out DIR [--patient ID | --patient-from parent|dirname]
                              [--areas <coordinates.csv|none> ...] [--clusters <clusters.csv|none> ...] [--jobs N] [--resume] [--no-keep-parts]`
`python -m bio.kit xenium --samplesheet <deposit>/samplesheet.csv --dataset NAME --out DIR [--regions R ...]
                              [--clusters <clusters.csv|none> ...] [--jobs N] [--resume] [--no-keep-parts]`

Reads Xenium outs bundles, labels cells from marker genes, and writes in DIR:
  NAME.h5ad                    cells with obs annotation / patient / sample, obsm spatial (µm)
  NAME_annotation_report.txt   label counts, clusters table (cluster method), markers used, marker-leakage table (per section, plus pooled for a cohort)
  NAME.yaml                    draft manifest for `python -m bio.kit check|run`
One bundle = one patient. Several bundles are annotated one by one (own background and marker-leakage table per
section) and then stacked into one h5ad; the manifest reads patient and field from obs columns.

A slide that carries several patients is cut with `--areas`: one `<region>_coordinates.csv` per bundle (same order;
`none` = take the bundle whole). Each cell goes to the selection polygon that contains its centroid, cells outside every
polygon are dropped (and counted), and each selection becomes one section named after the selection. `--samplesheet`
does this for a whole deposit: it reads samplesheet.csv (slide, region, sample_areas) and expands it into bundles
(`<folder of the samplesheet>/<region>`) and areas.

Labelling is by 10x cluster when the bundle has analysis/clustering/gene_expression_graphclust/clusters.csv (each cluster
is called from how specific the lineage marker sets are in it, cells inherit the label; best for sparse panels) and cell by
cell otherwise. `--clusters` gives an explicit clusters.csv (Barcode,Cluster) per bundle, in bundle order, `none` = no
clusters for that bundle (cell-by-cell labels). On a slide cut by `--areas`, the clusters are evaluated inside each cut section.

Work is split in two stages. Stage A reads each region and cuts it into sections (at most min(jobs, 3) regions at once:
a region can hold 1M cells) and writes one raw per-section file. Stage B annotates the sections in `--jobs` workers,
largest section first, and writes the annotated section and its report to `<out>/parts/<dataset>/<section>.h5ad` / `.txt`
(atomically, so a killed run never leaves a half-written part). One line per finished section goes to stderr. The main
process then assembles the cohort from the parts in samplesheet / section order, so the output is identical whatever the
completion order and whatever `--jobs` is. `--resume` reuses the parts that are complete (same inputs, same code) and
only does the rest; `--no-keep-parts` deletes the parts after a successful run. Default `--jobs` is max(1, CPUs - 1).
A worker needs about 1 GB per 200k cells of its section (a 1M-cell section about 5 GB), plus the stacked cohort in the
main process: on a small-memory machine pass a lower `--jobs`.

obs columns written for every section: patient = the full section name (sections stay separate images), sample,
tissue (`_c` core, `_m` margin, `_n` normal, else unknown) and patient_id (the name without that suffix, so the
core / margin / normal sections of one patient share it).
"""
import hashlib
import json
import multiprocessing as mp
import os
import re
import shutil
import sys
import time
import warnings
import csv
import re as _re
from pathlib import Path

from bio.annotate import annotate, annotate_report, merge_annotation_params, report_text, LABELS
from bio.xenium import GRAPHCLUST, MATRIX, _find, concat_sections, cut_by_areas, read_areas, read_bundle

IDENTITY = ["CD3D", "CD3E", "CD3G", "CD2", "CD247"]
STATES = {"prolif": ["MKI67", "TOP2A", "CDK1", "PCNA", "MCM2", "TYMS"], "PD1": ["PDCD1"], "GZMB": ["GZMB"],
          "LAG3": ["LAG3"], "ICOS": ["ICOS"]}
BLEED = {"neut": ["CSF3R", "CXCR1", "CXCR2", "FCGR3B", "FPR1"],
         "mac": ["CD68", "CD163", "C1QA", "C1QB", "C1QC", "MRC1", "APOE"],
         "stroma": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PDGFRB"]}
MAC_STATES = {"PDL1": ["CD274"], "IDO1": ["IDO1"], "M2": ["CD163", "MRC1", "MSR1"]}
MAC_IDENTITY = ["CD68", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC"]
ROLES = {"CD4T": "T:CD4", "CD8T": "T:CD8", "Treg": "T:Treg", "T": "T", "Neutrophil": "neut", "Macrophage": "mac",
         "B": "B", "Plasma": "B", "Fibroblast": "stroma", "SMC": "stroma", "Epithelial": "tumor",
         "Endothelial": "other", "Monocyte": "other", "DC": "other", "Mast": "other", "NK": "other",
         "Unclustered": "other", "Ambiguous": "other", "LowCount": "other"}  # Unclustered/Ambiguous/LowCount stay as cells (they take up space) but have no role
LABEL_LINEAGE = {"CD4T": "T", "CD8T": "T", "Treg": "T"}   # "T" is its own lineage


def add_args(p):
    p.add_argument("bundle", nargs="*", help="Xenium outs folder(s) (cell_feature_matrix.h5 + cells.csv.gz/parquet) or .h5ad file(s); "
                   "several = one patient each. Not needed with --samplesheet")
    p.add_argument("--areas", nargs="+", metavar="CSV",
                   help="cut slides into patients: one <region>_coordinates.csv (Selection,X,Y polygons in µm) per bundle, "
                        "in the same order as the bundles; 'none' = use that bundle whole. Give it after the bundles")
    p.add_argument("--samplesheet", help="deposit samplesheet.csv (slide, region, sample_areas): reads every region "
                   "as <folder of the samplesheet>/<region> and cuts it with its coordinates file when sample_areas names one")
    p.add_argument("--xenium-dir", help="with --samplesheet: folder holding the region folders and coordinates files "
                   "(default: the folder of the samplesheet)")
    p.add_argument("--regions", nargs="+", metavar="REGION", help="with --samplesheet: only these regions")
    p.add_argument("--clusters", nargs="+", metavar="CSV",
                   help="10x clusters.csv (Barcode,Cluster) per bundle, in bundle order; 'none' = no clusters for that bundle "
                        "(cell-by-cell labels). Default: analysis/clustering/gene_expression_graphclust/clusters.csv inside each "
                        "bundle when present, labelling whole clusters instead of single cells (better on sparse panels)")
    p.add_argument("--jobs", type=int, default=None, metavar="N",
                   help="workers for the annotation stage (default: CPUs - 1; capped at the CPU count). The read-and-cut stage "
                        "uses at most min(N, 3) workers. A worker needs about 1 GB per 200k cells of its section, so on a "
                        "small-memory machine pass a lower N. Output is identical for any N")
    p.add_argument("--resume", action="store_true",
                   help="reuse the complete section files in <out>/parts/<dataset>/ from an earlier run (same inputs and "
                        "code) and only do the sections that are missing")
    p.add_argument("--no-keep-parts", dest="keep_parts", action="store_false",
                   help="delete <out>/parts/<dataset>/ after a successful run (default: keep it, so --resume works and the "
                        "annotated sections are available one by one)")
    p.add_argument("--dataset", required=True, help="name used for the output files, e.g. XENIUM_LUNG_5K")
    p.add_argument("--out", required=True, help="output folder")
    p.add_argument("--patient", help="patient ID for a single bundle (default: the dataset name)")
    p.add_argument("--patient-from", choices=["parent", "dirname"], default="parent",
                   help="several bundles: patient = last folder name that is not 'outs' (parent, default) or the "
                        "final folder name exactly (dirname)")
    p.add_argument("--min-counts", type=int, default=20, help="cells with fewer counts are LowCount (default 20)")


def patient_name(path, mode="parent"):
    """Patient ID from a bundle path. parent: last component that is not 'outs' (.../P01/outs -> P01); dirname: final component."""
    p = Path(path)
    if p.suffix == ".h5ad":
        return p.stem
    parts = [x for x in p.resolve().parts if x != os.sep]
    if mode == "dirname":
        return parts[-1]
    keep = [x for x in parts if x != "outs"]
    return keep[-1] if keep else parts[-1]


_SUFFIX = _re.compile(r"^(?P<id>.+)_(?P<t>[cmn])$")
TISSUE = {"c": "core", "m": "margin", "n": "normal"}


def region_info(name):
    """(patient_id, tissue) from a section name: run1_crc_b_m -> ('run1_crc_b', 'margin'). No _c/_m/_n suffix ->
    (name, 'unknown'). The suffix meaning is our reading of the names; the deposit does not state it."""
    m = _SUFFIX.match(str(name))
    return (m["id"], TISSUE[m["t"]]) if m else (str(name), "unknown")


def read_samplesheet(path, xenium_dir=None, regions=None):
    """samplesheet.csv (slide, region, sample_areas) -> [(bundle_dir, areas_csv or None)] in file order."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"samplesheet {path} not found")
    base = Path(xenium_dir) if xenium_dir else path.parent
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or not {"region", "sample_areas"} <= set(rows[0]):
        raise ValueError(f"{path.name}: expected the columns slide, region, sample_areas (found {list(rows[0]) if rows else 'no rows'})")
    out, missing = [], []
    for r in rows:
        reg = (r["region"] or "").strip()
        if not reg or (regions and reg not in regions):
            continue
        sa = (r["sample_areas"] or "").strip()
        areas = None if sa.lower() in ("", "none", "nan", "na") else base / sa
        if not (base / reg).is_dir():
            missing.append(f"{base / reg} (region folder)")
        if areas is not None and not areas.is_file():
            missing.append(f"{areas} (coordinates file)")
        out.append((base / reg, areas))
    if regions:
        absent = sorted(set(regions) - {Path(b).name for b, _ in out})
        if absent:
            raise ValueError(f"regions not in the samplesheet: {', '.join(absent)}")
    if missing:
        raise ValueError("missing: " + "; ".join(missing))
    if not out:
        raise ValueError(f"{path.name}: no regions selected")
    return out


def _q(s):
    return json.dumps(str(s))  # a quoted string is valid YAML


def _genes(adata, genes):
    return [g for g in genes if g in adata.var_names]


def _flow(genes):
    return "[" + ", ".join(genes) + "]"


def build_manifest(adata, dataset, patient, h5ad_path, rep, per_sample=False, section_unreliable=None, notes=None):
    """Draft manifest text. Gene lists are cut to the panel; states with no gene in the panel are left out.
    per_sample: patient and field are read from the obs columns `patient` / `sample` (cohort), not constants."""
    section_unreliable = section_unreliable or {}
    flagged = set(rep["unreliable"]) | {l for ls in section_unreliable.values() for l in ls}
    L = [f"# Draft written by `python -m bio.kit xenium`. Labels are marker calls; review the roles and the marker-leakage table in",
         f"# {dataset}_annotation_report.txt before running. Epithelial -> tumor does not separate malignant from normal cells."]
    if flagged:
        where = "; ".join(f"{s}: {', '.join(ls)}" for s, ls in section_unreliable.items() if ls)
        L.append(f"# WARNING: unreliable lineages: {', '.join(sorted(flagged))}"
                 + (f" (per section: {where})" if where else "") + ". Their roles are kept below; tests that use them are not trustworthy.")
    L += [f"# {x}" for x in (notes or [])]
    if rep["unavailable"]:
        L.append(f"# NOTE: lineages with <2 markers in this panel (and no specific single marker) cannot be called: {', '.join(rep['unavailable'])}.")
    cols = ("{patient: patient, field: sample, " if per_sample else f"{{patient: {_q(patient)}, field: {_q(dataset)}, ")
    path = os.path.relpath(Path(h5ad_path).resolve(), Path.cwd())  # relative: manifests get committed
    L += [f"dataset: {_q(dataset)}", "modality: rna", f"path: {_q(path)}",
          f"columns: {cols}x: \"obsm:spatial:0\", y: \"obsm:spatial:1\", type: annotation}}"
          + ("   # one patient per section: obs columns patient / sample" if per_sample else "   # one patient, one section"),
          "pixel_size_um: 1.0   # Xenium centroids are already in µm"]
    ident = _genes(adata, IDENTITY)
    L.append(f"identity_genes: {_flow(ident)}" + ("" if ident else "   # FILL IN: no T-identity gene in this panel"))
    L.append("roles:")
    for lab, role in ROLES.items():
        flag = "   # UNRELIABLE labels, see report" if LABEL_LINEAGE.get(lab, lab) in flagged else ""
        L.append(f"  {lab}: {_q(role)}{flag}")
    L.append("states:")
    for s, genes in STATES.items():
        g = _genes(adata, genes)
        L.append(f"  {s}: {_flow(g)}" if g else f"  # {s}: no gene in this panel")
    L.append("bleed:")
    for r, genes in BLEED.items():
        g = _genes(adata, genes)
        L.append(f"  {r}: {_flow(g)}" if g else f"  # {r}: no gene in this panel")
    ms = {k: _genes(adata, v) for k, v in MAC_STATES.items()}
    L.append("mac_states:" if any(ms.values()) else "# mac_states: no gene in this panel")
    for k, g in ms.items():
        L.append(f"  {k}: {_flow(g)}" if g else f"  # {k}: no gene in this panel")
    mi = _genes(adata, MAC_IDENTITY)
    L.append(f"mac_identity_genes: {_flow(mi)}" if mi else "# mac_identity_genes: no gene in this panel")
    return "\n".join(L) + "\n"


def _clusters_arg(a, n):
    if not a.clusters:
        return [None] * n
    if len(a.clusters) != n:
        sys.exit(f"--clusters needs one entry per bundle ({n} bundles, {len(a.clusters)} entries; use 'none' to skip one)")
    return ["none" if x.lower() == "none" else x for x in a.clusters]


def _plan(a):
    """-> (bundles, areas, clusters) from the positional bundles / --areas / --clusters / --samplesheet, or exit with a
    plain message. clusters entries: None = automatic, 'none' = skip, else a path."""
    if a.samplesheet:
        if a.bundle or a.areas or a.patient:
            sys.exit("--samplesheet already names the bundles and areas; do not also give bundles, --areas or --patient")
        try:
            plan = read_samplesheet(a.samplesheet, a.xenium_dir, a.regions)
        except ValueError as e:
            sys.exit(f"samplesheet: {e}")
        return [str(b) for b, _ in plan], [x for _, x in plan], _clusters_arg(a, len(plan))
    if a.regions or a.xenium_dir:
        sys.exit("--regions / --xenium-dir only work with --samplesheet")
    if not a.bundle:
        sys.exit("give at least one Xenium bundle (or --samplesheet)")
    areas = [None] * len(a.bundle)
    if a.areas:
        if len(a.areas) != len(a.bundle):
            sys.exit(f"--areas needs one entry per bundle ({len(a.bundle)} bundles, {len(a.areas)} entries; use 'none' to skip one)")
        areas = [None if x.lower() == "none" else x for x in a.areas]
    return a.bundle, areas, _clusters_arg(a, len(a.bundle))


class RegionError(Exception):
    """A bundle / region could not be read, cut or labelled (message is shown to the user as is)."""


def _safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", str(name)) or "_"


def _stat(p):
    try:
        st = os.stat(p)
        return [str(Path(p).resolve()), st.st_size, st.st_mtime_ns]
    except OSError:
        return [str(p), None, None]


def _code_hash():
    """Parts are only reused by the code that wrote them: hash of the labelling, reading and cutting sources."""
    h = hashlib.sha256()
    for m in ("bio.annotate", "bio.xenium"):
        h.update(Path(__import__(m, fromlist=["x"]).__file__).read_bytes())
    return h.hexdigest()[:16]


def _region_key(task):
    """Fingerprint of everything that decides what a region's sections look like."""
    b = Path(task["bundle"])
    src = b if b.suffix == ".h5ad" else (_find(b, [MATRIX]) or b)
    cl = task["cluster"]
    cl_src = None if (cl is not None and str(cl).lower() == "none") else (cl if cl is not None else (_find(b, [GRAPHCLUST]) if b.is_dir() else None))
    blob = json.dumps([_stat(src), _stat(task["area"]) if task["area"] else None,
                       "none" if cl_src is None else _stat(cl_src), task["min_counts"], task["name"], _code_hash()])
    return hashlib.sha256(blob.encode()).hexdigest()[:32]


def _section_key(rkey, name):
    return hashlib.sha256(f"{rkey}|{name}".encode()).hexdigest()[:32]


def _part_paths(pdir, safe):
    return Path(pdir) / f"{safe}.h5ad", Path(pdir) / f"{safe}.txt"


def _part_ok(pdir, safe, key):
    """A part is complete when both files exist (the .txt is renamed into place last) and carry this key."""
    h5, txt = _part_paths(pdir, safe)
    try:
        return h5.is_file() and txt.is_file() and txt.read_text(encoding="utf-8").split("\n", 1)[0] == f"# part-key {key}"
    except OSError:
        return False


def _atomic_text(path, text):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _cut_region(task):
    """Stage A. Read one bundle, cut it by its areas file if it has one, write one raw (unannotated) file per section that
    still needs annotating and a region record (<pdir>/region_NNN.json). Top-level (picklable) so a spawned worker can run
    it. Returns the record plus dict keys `skipped` (region not read: all its sections are complete parts) and `raw`
    (section name -> raw file). task: index, bundle, area, cluster, name, min_counts, safe {section: file stem}, pdir,
    resume, rkey."""
    t0 = time.time()
    b, ar = task["bundle"], task["area"]
    pdir = Path(task["pdir"])
    rec_path = pdir / f"region_{task['index']:03d}.json"
    if task["resume"] and rec_path.is_file():
        try:
            rec = json.loads(rec_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            rec = None
        if rec and rec.get("rkey") == task["rkey"] and all(_part_ok(pdir, x["safe"], x["key"]) for x in rec["sections"]):
            return {**rec, "skipped": True, "raw": {}, "seconds": time.time() - t0}
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ad_ = read_bundle(b, clusters=task["cluster"] if task["cluster"] is not None else "auto")
        wmsg = [str(x.message) for x in w if issubclass(x.category, UserWarning)]
    except (FileNotFoundError, ValueError, ImportError) as e:
        raise RegionError(f"cannot read {b}: {e}") from None
    if "spatial" not in ad_.obsm:
        raise RegionError(f"cannot use {b}: no obsm['spatial'] coordinates (µm)")
    n_cells = int(ad_.n_obs)
    note, warn = None, bool(wmsg)
    if ar is None:
        sections = [(task["name"], ad_, False)]
    else:
        try:
            parts, info = cut_by_areas(ad_, read_areas(ar))
        except (FileNotFoundError, ValueError) as e:
            raise RegionError(f"cannot use areas {ar} for {b}: {e}") from None
        n_in = info["n_cells"] - info["n_outside"]
        note = (f"{Path(b).name} cut by {Path(ar).name}: {info['n_cells']} cells, {info['n_outside']} outside every "
                f"selection (dropped), {n_in} kept in {len(parts)} selection(s): "
                + ", ".join(f"{n} {c}" for n, c in info["counts"].items())
                + (f"; {info['n_overlap']} cells lie in more than one polygon and went to the first" if info["n_overlap"] else "")
                + (f"; no cells in: {', '.join(info['empty'])} (left out)" if info["empty"] else ""))
        warn = warn or bool(info["n_outside"] or info["empty"] or info["n_overlap"])
        sections = [(n, s_, True) for n, s_ in parts]
        del ad_
    rdir = pdir / "raw"
    rdir.mkdir(parents=True, exist_ok=True)
    recs, raw = [], {}
    for name, sec, cut in sections:
        safe = task["safe"][name]
        key = _section_key(task["rkey"], name)
        recs.append({"name": name, "safe": safe, "key": key, "cells": int(sec.n_obs), "cut": cut})
        if not (task["resume"] and _part_ok(pdir, safe, key)):      # a complete part needs no raw file
            tmp, fin = rdir / f"{safe}.tmp.h5ad", rdir / f"{safe}.h5ad"
            sec.write_h5ad(tmp)
            os.replace(tmp, fin)
            raw[name] = str(fin)
    rec = {"index": task["index"], "rkey": task["rkey"], "bundle": str(b), "n_cells": n_cells, "note": note, "warn": warn,
           "messages": wmsg, "sections": recs}
    _atomic_text(rec_path, json.dumps(rec))
    return {**rec, "skipped": False, "raw": raw, "seconds": time.time() - t0}


def _annotate_section(item):
    """Stage B. Annotate one raw section file and write <pdir>/<safe>.h5ad and .txt atomically; remove the raw file.
    item: name, safe, key, raw, pdir, min_counts, cells. Returns (name, cells, seconds)."""
    import anndata as ad
    t0 = time.time()
    name = item["name"]
    sec = ad.read_h5ad(item["raw"])
    try:
        annotate(sec, min_counts=item["min_counts"])
    except ValueError as e:
        raise RegionError(f"cannot annotate {name}: {e}") from None
    text = f"# part-key {item['key']}\n" + report_text(annotate_report(sec), name)
    h5, txt = _part_paths(item["pdir"], item["safe"])
    tmp = Path(item["pdir"]) / f"{item['safe']}.tmp.h5ad"
    sec.write_h5ad(tmp)
    os.replace(tmp, h5)
    txt.unlink(missing_ok=True)       # the key file goes last: a part is complete only once it is there
    _atomic_text(txt, text)
    Path(item["raw"]).unlink(missing_ok=True)
    return name, int(item["cells"]), time.time() - t0


def _map(fn, items, workers, on_result):
    """fn over items, results to on_result as they finish. workers <= 1 runs in this process; else a spawn pool that
    hands the items out in the order given (chunksize 1), so a sorted list is processed largest first."""
    if workers <= 1 or len(items) <= 1:
        for it in items:
            on_result(fn(it))
    else:
        with mp.get_context("spawn").Pool(min(workers, len(items))) as pool:
            for r in pool.imap_unordered(fn, items, chunksize=1):
                on_result(r)


def _run_stages(a, tasks, jobs, pdir):
    """Stage A (read + cut, at most min(jobs, 3) workers) then stage B (annotate, `jobs` workers, largest section first),
    then load the finished sections in samplesheet / section order. -> (region records in task order, reports,
    [(name, adata, was_cut)], area notes). Progress lines go to stderr."""
    import anndata as ad
    recs = {}

    def done_a(r):
        recs[r["index"]] = r
        nm = Path(r["bundle"]).name
        if r["skipped"]:
            print(f"skip cutting {nm}: its {len(r['sections'])} section(s) are complete parts [{len(recs)}/{len(tasks)}]",
                  file=sys.stderr, flush=True)
        else:
            print(f"read {nm}: {r['n_cells']} cells -> {len(r['sections'])} section(s), {r['seconds']:.0f} s "
                  f"[{len(recs)}/{len(tasks)}]", file=sys.stderr, flush=True)

    _map(_cut_region, tasks, min(jobs, 3), done_a)
    order = [recs[i] for i in range(len(tasks))]
    flat = [(r, x) for r in order for x in r["sections"]]
    total, finished, items = len(flat), 0, []
    for r, x in flat:
        raw = r["raw"].get(x["name"])
        if raw is None:
            finished += 1
            print(f"kept {x['name']} (from parts) [{finished}/{total}]", file=sys.stderr, flush=True)
        else:
            items.append({"name": x["name"], "safe": x["safe"], "key": x["key"], "raw": raw, "pdir": str(pdir),
                          "min_counts": a.min_counts, "cells": x["cells"]})
    items.sort(key=lambda it: -it["cells"])     # stable: ties stay in section order; the pool hands items out in this order

    def done_b(res):
        nonlocal finished
        finished += 1
        print(f"done {res[0]}: {res[1]} cells, {res[2]:.0f} s [{finished}/{total}]", file=sys.stderr, flush=True)

    _map(_annotate_section, items, jobs, done_b)

    sections, reps = [], []
    for r, x in flat:                           # always from the part files, so a resumed run assembles the same bytes
        sec = ad.read_h5ad(_part_paths(pdir, x["safe"])[0])
        sections.append((x["name"], sec, x["cut"]))
        reps.append(annotate_report(sec))
    return order, reps, sections, [r["note"] for r in order if r["note"]]


def run(a):
    if not re.fullmatch(r"[A-Za-z0-9._-]+", a.dataset):
        sys.exit("--dataset may contain only letters, digits, '.', '_' and '-' (it names the output files)")
    if a.jobs is not None and a.jobs < 1:
        sys.exit("--jobs must be 1 or more")
    bundles, areas, clusters = _plan(a)
    n_bundles = len(bundles)
    if n_bundles > 1 and a.patient:
        sys.exit("--patient names a single bundle; with several bundles the patient comes from the folder name (see --patient-from)")

    tasks = []
    for i, (b, ar, cl) in enumerate(zip(bundles, areas, clusters)):
        nm = (patient_name(b, "dirname") if a.samplesheet  # the region folder is the section name, even for one region
              else patient_name(b, a.patient_from) if n_bundles > 1 else (a.patient or a.dataset))
        if ar is not None:  # fail early on a bad areas file, and on section names that would clash
            try:
                names_here = list(dict.fromkeys(n for n, _ in read_areas(ar)))
            except (FileNotFoundError, ValueError) as e:
                sys.exit(f"cannot use areas {ar} for {b}: {e}")
        else:
            names_here = [nm]
        tasks.append({"index": i, "bundle": b, "area": ar, "cluster": cl, "name": nm, "min_counts": a.min_counts,
                      "names_here": names_here})
    planned = [n for t in tasks for n in t["names_here"]]
    if len(set(planned)) < len(planned):
        dup = sorted({n for n in planned if planned.count(n) > 1})
        sys.exit(f"two sections have the same name {dup}; try --patient-from dirname or rename the folders")
    stems, used_stems = {}, set()          # section name -> file stem for the part files (unique, filesystem-safe)
    for n in planned:
        st, k = _safe(n), 1
        while st in used_stems:
            st, k = f"{_safe(n)}-{k}", k + 1
        stems[n] = st
        used_stems.add(st)

    cpus = os.cpu_count() or 1
    if a.jobs is None:
        jobs = max(1, cpus - 1)
    else:
        jobs = min(a.jobs, cpus)
        if a.jobs > cpus:
            print(f"NOTE: --jobs {a.jobs} capped at the CPU count ({cpus})", file=sys.stderr)

    out = Path(a.out)
    pdir = out / "parts" / a.dataset
    pdir.mkdir(parents=True, exist_ok=True)
    for t in tasks:
        t["safe"] = {n: stems[n] for n in t.pop("names_here")}
        t["pdir"], t["resume"], t["rkey"] = str(pdir), bool(a.resume), _region_key(t)
    hint = f"\nFinished sections are kept in {pdir}; rerun the same command with --resume to continue from them."
    try:
        recs, reps, sections, area_notes = _run_stages(a, tasks, jobs, pdir)
    except RegionError as e:
        sys.exit(str(e) + hint)
    except SystemExit as e:
        if e.code not in (0, None):
            print(hint.strip(), file=sys.stderr)
        raise
    except Exception:
        print(hint.strip(), file=sys.stderr)
        raise

    for r in recs:
        for m in r["messages"]:
            print(f"WARNING: {m}", file=sys.stderr)
        if r["note"]:
            print(f"WARNING: {r['note']}" if r["warn"] else r["note"], file=sys.stderr)
    if not sections:
        sys.exit("no cells left: every selection is empty" + hint)
    names = [n for n, _, _ in sections]
    multi = len(names) > 1
    per_sample = multi or any(c for _, _, c in sections)

    for name, ad_, cut in sections:
        pid, tissue = region_info(name)
        ad_.obs["patient"] = name
        ad_.obs["sample"] = name if (multi or cut) else a.dataset
        ad_.obs["patient_id"] = pid
        ad_.obs["tissue"] = tissue
    sections = [s_ for _, s_, _ in sections]

    notes = []
    if multi:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            adata = concat_sections(sections, names)
        for x in w:
            print(f"WARNING: {x.message}", file=sys.stderr)
        diff = json.loads(adata.uns["panel_differences"])
        if diff:
            notes.append("Panels differ between sections; a gene missing from a panel counts as 0 there:")
            notes += [f"  {n}: {', '.join(g)}" for n, g in diff.items()]
        adata.obs["annotation"] = adata.obs["annotation"].astype(str).astype("category").cat.set_categories(LABELS)
        adata.obs["annotation_note"] = adata.obs["annotation_note"].astype(str).astype("category")
        for c in ("patient", "sample", "patient_id", "tissue"):
            adata.obs[c] = adata.obs[c].astype(str)
        if "cluster_10x" in adata.obs:  # 10x cluster ids are per slide / section; "" = no clusters for that section
            adata.obs["cluster_10x"] = adata.obs["cluster_10x"].astype(object).where(adata.obs["cluster_10x"].notna(), "") \
                .astype(str).astype("category")
        adata.uns["annotation_params"] = merge_annotation_params(sections, names)
        rep = annotate_report(adata)
    else:
        adata, rep = sections[0], reps[0]

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    h5 = out / f"{a.dataset}.h5ad"
    adata.write_h5ad(h5)

    # what the manifest and the report say about how sections were made and named
    info_lines = list(area_notes)
    info_lines.append("obs columns: patient = full section name (sections stay separate images); sample; "
                      "patient_id = the name without its _c / _m / _n suffix (links the sections of one patient); "
                      "tissue = core (_c) / margin (_m) / normal (_n) / unknown. The suffix meaning is our reading of the "
                      "names, not stated by the deposit.")
    ti = {n: region_info(n) for n in names}
    info_lines.append("section -> patient_id / tissue: " + "; ".join(f"{n} -> {p_}/{t}" for n, (p_, t) in ti.items()))

    if multi:
        parts = [f"Annotation report: {a.dataset} ({len(names)} sections, {adata.n_obs} cells)", ""]
        parts += info_lines + [""]
        parts += ["Sections: " + "; ".join(
            f"{n} ({r['n_cells']} cells" + (f", unreliable: {', '.join(r['unreliable'])}" if r["unreliable"] else "") + ")"
            for n, r in zip(names, reps)), ""] + (notes + [""] if notes else [])
        for n, r in zip(names, reps):
            parts += ["=" * 78, f"SECTION {n}", "=" * 78, report_text(r, n)]
        parts += ["=" * 78, "COHORT (all sections pooled; the leakage table pools cells across sections)", "=" * 78,
                  report_text(rep, a.dataset)]
        text = "\n".join(parts)
    else:
        text = ("\n".join(info_lines) + "\n\n" if area_notes else "") + report_text(rep, a.dataset)
    (out / f"{a.dataset}_annotation_report.txt").write_text(text)
    sec_unrel = {n: list(r["unreliable"]) for n, r in zip(names, reps)} if multi else {}
    (out / f"{a.dataset}.yaml").write_text(
        build_manifest(adata, a.dataset, names[0], h5, rep, per_sample=per_sample, section_unreliable=sec_unrel,
                       notes=[x.replace("\n", " ") for x in info_lines]))
    print(text)
    shutil.rmtree(pdir / "raw", ignore_errors=True)
    if a.keep_parts:
        kept = f"\n      {pdir}/  (annotated sections; --resume reuses them)"
    else:
        shutil.rmtree(pdir, ignore_errors=True)
        try:
            (out / "parts").rmdir()
        except OSError:
            pass
        kept = ""
    print(f"wrote {h5}\n      {out / (a.dataset + '_annotation_report.txt')}\n      {out / (a.dataset + '.yaml')}{kept}")
    print("Review the roles and the marker-leakage table before running.")
