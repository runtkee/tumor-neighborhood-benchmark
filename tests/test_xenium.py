"""Xenium input path on a synthetic 10x-format bundle: reader, marker annotation, specificity flag, CLI."""
import gzip
import json
import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from bio.annotate import LINEAGE_MARKERS, annotate, annotate_report
from bio.kit import manifest
from bio.kit.__main__ import main
from bio.xenium import read_bundle

PLANTED = ["Epithelial", "Fibroblast", "Endothelial", "SMC", "CD8T", "CD4T", "Treg", "B", "Plasma", "Macrophage",
           "Monocyte", "DC", "Neutrophil", "Mast", "NK"]
PROPORTION = np.array([.25, .15, .06, .04, .07, .07, .03, .05, .03, .08, .02, .02, .06, .03, .04])
T_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD247"]
OTHER = ["CD14", "CSF1R", "MRC1", "CD274", "IDO1", "MSR1", "MKI67", "PDCD1", "GZMB", "LAG3", "ICOS", "APOE", "ACTB", "GAPDH"]


def _panel():
    genes = []
    for l, m in LINEAGE_MARKERS.items():
        genes += m if l == "Neutrophil" else m[:5]
    genes += ["CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4"] + OTHER
    return list(dict.fromkeys(genes))


def _lineage(label):
    return "T" if label in ("CD8T", "CD4T", "Treg") else label


MIXED_ID = "99"


def _make(tmp_path, seed=0, n=2000, parquet=False, drop=None, ambient=None, outlier=False, clusters=False, proportion=None):
    """drop: gene or collection of genes left out of the panel; ambient: {gene: mean counts in every cell};
    clusters: also write analysis/clustering/gene_expression_graphclust/clusters.csv (Barcode,Cluster) whose clusters
    are the planted labels (ids 1..15), except that the first 30 Plasma and 30 Mast cells form one mixed cluster (id 99)
    (clusters="pure": no mixed cluster), and a kmeans_5 file that must not be read.
    proportion: share of each PLANTED label (default PROPORTION)."""
    rng = np.random.default_rng(seed)
    dropped = {drop} if isinstance(drop, str) else set(drop or ())
    genes = [x for x in _panel() if x not in dropped] + [x for x in (ambient or {})]
    g = {x: i for i, x in enumerate(genes)}
    pr = PROPORTION if proportion is None else np.asarray(proportion, float)
    truth = rng.choice(PLANTED, n, p=pr / pr.sum())
    M = rng.poisson(0.03, (n, len(genes))).astype(np.int32)           # noise everywhere
    M[:, [g["ACTB"], g["GAPDH"]]] += rng.poisson(40, (n, 2))            # housekeeping: realistic library size (>20 counts)
    for i, t in enumerate(truth):
        for x in LINEAGE_MARKERS[_lineage(t)]:
            if x in g:
                M[i, g[x]] += rng.poisson(3)
        if t == "CD8T":
            M[i, [g["CD8A"], g["CD8B"]]] += rng.poisson(3, 2)
        if t in ("CD4T", "Treg"):
            M[i, g["CD4"]] += rng.poisson(3)
        if t == "Treg":
            M[i, g["FOXP3"]] += rng.poisson(3)
    for x, rate in (ambient or {}).items():
        M[:, g[x]] += rng.poisson(rate, n)
    neg = rng.poisson(0.02, (n, 20)).astype(np.int32)
    unassigned = rng.poisson(0.5, (n, 2)).astype(np.int32)              # must be ignored
    full = np.hstack([M, neg, unassigned])
    names = genes + [f"NegControlProbe_{i}" for i in range(10)] + [f"NegControlCodeword_{i}" for i in range(10)] \
        + ["UnassignedCodeword_0", "UnassignedCodeword_1"]
    types = ["Gene Expression"] * len(genes) + ["Negative Control Probe"] * 10 + ["Negative Control Codeword"] * 10 \
        + ["Unassigned Codeword"] * 2
    barcodes = [f"cell{i:05d}-1" for i in range(n)]
    xy = rng.uniform(0, 600, (n, 2))                                    # µm; ~6.7 µm between cells
    if outlier:
        xy[0] = (-50.0, -50.0)                                          # cell00000-1 sits outside every polygon

    bundle = tmp_path / "outs"
    bundle.mkdir()
    csc = sparse.csc_matrix(full.T)                                     # features x cells, as 10x writes it
    with h5py.File(bundle / "cell_feature_matrix.h5", "w") as f:
        m = f.create_group("matrix")
        m["data"], m["indices"], m["indptr"] = csc.data, csc.indices.astype(np.int64), csc.indptr.astype(np.int64)
        m["shape"] = np.array(csc.shape, dtype=np.int32)
        m["barcodes"] = np.array(barcodes, dtype="S")
        fe = m.create_group("features")
        fe["name"] = np.array(names, dtype="S")
        fe["id"] = np.array([f"id{i}" for i in range(len(names))], dtype="S")
        fe["feature_type"] = np.array(types, dtype="S")
        fe["genome"] = np.array(["Unknown"] * len(names), dtype="S")
    cells = pd.DataFrame({"cell_id": barcodes[::-1], "x_centroid": xy[::-1, 0], "y_centroid": xy[::-1, 1],  # reversed: must be aligned
                          "transcript_counts": full[::-1].sum(1), "cell_area": 50.0})
    if parquet:
        cells.to_parquet(bundle / "cells.parquet")
    else:
        with gzip.open(bundle / "cells.csv.gz", "wt") as fh:
            cells.to_csv(fh, index=False)
    if clusters:
        cid = np.array([str(PLANTED.index(t) + 1) for t in truth], dtype=object)
        for lab in (("Plasma", "Mast") if clusters != "pure" else ()):
            cid[np.where(truth == lab)[0][:30]] = MIXED_ID
        d = bundle / "analysis" / "clustering" / "gene_expression_graphclust"
        d.mkdir(parents=True)
        pd.DataFrame({"Barcode": barcodes, "Cluster": cid}).to_csv(d / "clusters.csv", index=False)
        k = bundle / "analysis" / "clustering" / "gene_expression_kmeans_5_clusters"
        k.mkdir()
        pd.DataFrame({"Barcode": barcodes, "Cluster": rng.integers(1, 6, n)}).to_csv(k / "clusters.csv", index=False)
    return bundle, truth, xy, genes, neg


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    return _make(tmp_path_factory.mktemp("xen"))


def test_read_bundle(synth):
    bundle, truth, xy, genes, neg = synth
    a = read_bundle(bundle)
    assert a.shape == (2000, len(genes)) and list(a.var_names) == genes     # negatives and unassigned are not in X
    assert sparse.isspmatrix_csr(a.X) and a.X.dtype == np.float32
    np.testing.assert_allclose(a.obsm["spatial"], xy)                       # aligned despite reversed cells table
    assert a.obs_names[0] == "cell00000-1"
    assert (a.obs["neg_control_n_features"] == 20).all()
    np.testing.assert_array_equal(a.obs["neg_control_counts"].to_numpy(), neg.sum(1))
    assert "cell_area" in a.obs and "transcript_counts" in a.obs
    # parent folder works too, and an h5ad is returned as is
    assert read_bundle(bundle.parent).n_obs == 2000
    p = bundle.parent / "x.h5ad"
    a.write_h5ad(p)
    assert read_bundle(p).shape == a.shape


def test_read_bundle_parquet_and_errors(tmp_path):
    bundle, *_ = _make(tmp_path, parquet=True)
    assert read_bundle(bundle).n_obs == 2000
    with pytest.raises(FileNotFoundError, match="neither a directory"):
        read_bundle(tmp_path / "nothing")
    (tmp_path / "empty").mkdir()
    with pytest.raises(FileNotFoundError, match="cell_feature_matrix.h5"):
        read_bundle(tmp_path / "empty")
    (bundle / "cells.parquet").unlink()
    with pytest.raises(FileNotFoundError, match="cells"):
        read_bundle(bundle)


def test_annotate_recovers_planted(synth):
    bundle, truth, *_ = synth
    a = annotate(read_bundle(bundle))
    acc = (a.obs["annotation"].astype(str).to_numpy() == truth).mean()
    assert acc >= 0.90, acc
    for lab in PLANTED:                                                   # every class is recovered, not just the big ones
        m = truth == lab
        assert (a.obs["annotation"].astype(str).to_numpy()[m] == lab).mean() >= 0.80, lab
    rep = annotate_report(a)
    assert rep["unreliable"] == {}
    assert sum(rep["labels"].values()) == 2000
    assert set(rep["markers_used"]["Neutrophil"]) == set(LINEAGE_MARKERS["Neutrophil"])
    assert rep["specificity"].loc["Neutrophil", "Epithelial"] < 0.10
    assert np.isnan(rep["specificity"].loc["T", "CD8T"])                  # own labels are not scored


def test_low_count_and_no_negative_controls(synth):
    bundle, *_ = synth
    a = read_bundle(bundle)
    X = a.X.toarray()
    X[:10, :] = 0                                                         # ten empty cells
    a.X = sparse.csr_matrix(X)
    a = a[:, :].copy()
    a.obs = a.obs.drop(columns=["neg_control_counts", "neg_control_n_features"])
    annotate(a)
    assert (a.obs["annotation"].astype(str).iloc[:10] == "LowCount").all()
    assert "fixed floor" in annotate_report(a)["background_source"]


def test_specificity_flags_planted_neutrophil_bleed(synth):
    bundle, *_ = synth
    a = read_bundle(bundle)
    rng = np.random.default_rng(5)
    cols = [list(a.var_names).index(g) for g in LINEAGE_MARKERS["Neutrophil"]]
    X = a.X.toarray()
    X[:, cols] += rng.random((a.n_obs, len(cols))) < 0.15                  # every neutrophil marker seen in 15% of every cell
    a.X = sparse.csr_matrix(X)
    rep = annotate_report(annotate(a))
    assert "Neutrophil" in rep["unreliable"] and len(rep["unreliable"]["Neutrophil"]) >= 3
    assert set(rep["unreliable"]) == {"Neutrophil"}                        # only the contaminated lineage


def test_normalized_input_is_refused(synth):
    a = read_bundle(synth[0])
    a.X = a.X * 0.37
    with pytest.raises(ValueError, match="raw integer counts"):
        annotate(a)


def test_cli_writes_outputs_and_manifest(synth, tmp_path, capsys):
    bundle, *_ = synth
    out = tmp_path / "prepared"
    main(["xenium", str(bundle), "--dataset", "SYN_XEN", "--out", str(out), "--patient", "P7"])
    txt = capsys.readouterr().out
    assert "Review the roles and the marker-leakage table before running." in txt
    h5, rpt, yml = out / "SYN_XEN.h5ad", out / "SYN_XEN_annotation_report.txt", out / "SYN_XEN.yaml"
    assert h5.exists() and rpt.exists() and yml.exists()
    assert "Marker leakage" in rpt.read_text() and "Cells per label" in rpt.read_text()
    assert "Method: cell" in rpt.read_text()                                  # no clusters.csv in this bundle
    a = ad.read_h5ad(h5)
    assert {"annotation", "patient", "sample"} <= set(a.obs.columns)
    assert (a.obs["patient"] == "P7").all() and (a.obs["sample"] == "SYN_XEN").all()
    assert a.obsm["spatial"].shape == (2000, 2)

    m = manifest.load(yml)
    assert m["modality"] == "rna" and m["dataset"] == "SYN_XEN" and m["pixel_size_um"] == 1.0
    assert m["path"] == os.path.relpath(h5.resolve(), os.getcwd()) and not os.path.isabs(m["path"])
    assert m["columns"] == {"patient": "P7", "field": "SYN_XEN", "x": "obsm:spatial:0", "y": "obsm:spatial:1", "type": "annotation"}
    assert m["roles"]["CD8T"] == "T:CD8" and m["roles"]["CD4T"] == "T:CD4" and m["roles"]["Treg"] == "T:Treg"
    assert m["roles"]["Neutrophil"] == "neut" and m["roles"]["Macrophage"] == "mac" and m["roles"]["Epithelial"] == "tumor"
    assert set(m["roles"]) == set(__import__("bio.annotate", fromlist=["LABELS"]).LABELS)
    assert m["mac_states"] == {"PDL1": ["CD274"], "IDO1": ["IDO1"], "M2": ["CD163", "MRC1", "MSR1"]}
    assert m["mac_identity_genes"] == ["CD68", "CD14", "CSF1R", "C1QA", "C1QB", "C1QC"]
    assert m["identity_genes"] == ["CD3D", "CD3E", "CD3G", "CD2", "CD247"]
    assert m["states"]["PD1"] == ["PDCD1"] and m["bleed"]["neut"][0] == "CSF3R"
    probs = manifest.problems(m)
    # the "T:CD8" role values are accepted only once the role-subset change in manifest.py has landed
    assert all(p.startswith("unknown roles") and "T:" in p for p in probs), probs
    if not probs:
        from bio.kit.importers import load_cells
        df = load_cells(m)
        assert len(df) == 2000 and df["patient"].eq("P7").all()


def test_cli_bad_dataset_name_and_missing_bundle(tmp_path):
    with pytest.raises(SystemExit):
        main(["xenium", str(tmp_path), "--dataset", "bad name", "--out", str(tmp_path / "o")])
    (tmp_path / "empty").mkdir()
    with pytest.raises(SystemExit, match="cell_feature_matrix.h5"):
        main(["xenium", str(tmp_path / "empty"), "--dataset", "OK", "--out", str(tmp_path / "o")])


def test_cli_two_bundles_cohort(tmp_path, capsys):
    for pat, seed, drop in (("P01", 1, None), ("P02", 2, "APOE")):       # second panel lacks one gene
        (tmp_path / pat).mkdir()
        _make(tmp_path / pat, seed=seed, n=800, drop=drop)
    out = tmp_path / "kit"
    main(["xenium", str(tmp_path / "P01" / "outs"), str(tmp_path / "P02" / "outs"),
          "--dataset", "XEN_CRC", "--out", str(out)])
    cap = capsys.readouterr()
    txt = cap.out
    assert "panels differ" in cap.err and "APOE" in cap.err                # warned on stderr
    assert "SECTION P01" in txt and "SECTION P02" in txt and "COHORT" in txt
    assert "Panels differ" in txt and "P02: APOE" in txt
    a = ad.read_h5ad(out / "XEN_CRC.h5ad")
    assert a.n_obs == 1600 and a.obs_names.is_unique
    assert set(a.obs["patient"]) == {"P01", "P02"} and (a.obs["patient"] == a.obs["sample"]).all()
    assert a.obs["patient"].value_counts().to_dict() == {"P01": 800, "P02": 800}
    assert a.n_vars == len(_panel()) and "APOE" in a.var_names           # union of genes
    j = list(a.var_names).index("APOE")
    X = sparse.csr_matrix(a.X)
    p2 = (a.obs["patient"] == "P02").to_numpy()
    assert X[p2][:, j].sum() == 0 and X[~p2][:, j].sum() > 0              # missing from P02 -> 0 there, present in P01
    assert a.obsm["spatial"].shape == (1600, 2)
    truth_ok = (a.obs["annotation"].astype(str) != "Ambiguous").mean()
    assert truth_ok > 0.9
    rpt = (out / "XEN_CRC_annotation_report.txt").read_text()
    assert rpt.count("Cells per label") == 3 and "Marker leakage" in rpt   # two sections + cohort

    m = manifest.load(out / "XEN_CRC.yaml")
    assert m["columns"] == {"patient": "patient", "field": "sample", "x": "obsm:spatial:0", "y": "obsm:spatial:1", "type": "annotation"}
    assert m["path"] == os.path.relpath((out / "XEN_CRC.h5ad").resolve(), os.getcwd())
    assert m["roles"]["CD8T"] == "T:CD8"
    assert all(p.startswith("unknown roles") and "T:" in p for p in manifest.problems(m))
    # single-patient flag is refused for several bundles; duplicate names are refused
    with pytest.raises(SystemExit, match="single bundle"):
        main(["xenium", str(tmp_path / "P01" / "outs"), str(tmp_path / "P02" / "outs"), "--dataset", "X",
              "--out", str(out), "--patient", "A"])


def test_patient_name():
    from bio.kit.xenium_cmd import patient_name
    assert patient_name("data/crc/P01/outs") == "P01" and patient_name("data/crc/P02") == "P02"
    assert patient_name("data/crc/P01/outs", "dirname") == "outs" and patient_name("data/crc/P09_run", "dirname") == "P09_run"
    assert patient_name("x/P3.h5ad") == "P3"


def test_single_marker_lineage_epcam_only(tmp_path):
    bundle, truth, *_ = _make(tmp_path, drop=["KRT8", "KRT18", "KRT19", "CDH1"])
    a = annotate(read_bundle(bundle))
    rep = annotate_report(a)
    assert rep["markers_used"]["Epithelial"] == ["EPCAM"] and "Epithelial" not in rep["unavailable"]
    assert rep["single_marker"] == {"Epithelial": "EPCAM"}
    lab = a.obs["annotation"].astype(str).to_numpy()
    assert (lab[truth == "Epithelial"] == "Epithelial").mean() >= 0.70         # needs >=2 counts: Poisson(3) gives ~80%
    assert (truth[lab == "Epithelial"] == "Epithelial").mean() >= 0.95
    from bio.annotate import report_text
    txt = report_text(rep)
    assert "SINGLE-MARKER CALL" in txt and "Called from a single marker" in txt and "Epithelial (EPCAM)" in txt
    # a single marker that is not in SOLO_OK stays unavailable (Neutrophil with CSF3R only)
    (tmp_path / "x").mkdir()
    b2, *_ = _make(tmp_path / "x", drop=[g for g in LINEAGE_MARKERS["Neutrophil"] if g != "CSF3R"])
    r2 = annotate_report(annotate(read_bundle(b2)))
    assert "Neutrophil" in r2["unavailable"] and "Neutrophil" not in r2["single_marker"]
    # ACTA2 is specific enough for smooth muscle / myofibroblasts (SOLO_OK)
    (tmp_path / "y").mkdir()
    b3, *_ = _make(tmp_path / "y", drop=["MYH11", "TAGLN", "DES"])
    r3 = annotate_report(annotate(read_bundle(b3)))
    assert r3["single_marker"].get("SMC") == "ACTA2" and "SMC" not in r3["unavailable"]


def test_single_marker_needs_two_counts(tmp_path):
    bundle, *_ = _make(tmp_path, drop=["KRT8", "KRT18", "KRT19", "CDH1"], n=600)
    a = read_bundle(bundle)
    j = list(a.var_names).index("EPCAM")
    X = a.X.toarray()
    X[:, j] = 0
    X[:300, j] = 1                                                           # exactly one count of EPCAM
    a.X = sparse.csr_matrix(X)
    a = annotate(a)
    fired = np.asarray(a.obsm["lineage_fired"])[:, list(json_lineages(a)).index("Epithelial")]
    assert fired[:300].sum() == 0


def json_lineages(a):
    import json
    return json.loads(a.uns["annotation_params"])["lineages"]


def test_ambient_immunoglobulin_does_not_make_plasma(tmp_path):
    amb = {"IGKC": 6.0, "IGHG1": 4.0, "XBP1": 3.0, "SDC1": 3.0}
    bundle, truth, *_ = _make(tmp_path, ambient=amb)
    a = annotate(read_bundle(bundle))
    rep = annotate_report(a)
    assert rep["markers_used"]["Plasma"] == ["JCHAIN", "MZB1"] and rep["markers_requested"]["Plasma"] == ["JCHAIN", "MZB1"]
    assert rep["ambient_ig_excluded"] == {"Plasma": ["IGKC", "IGHG1"]}      # XBP1 / SDC1 (ambient here) are not Plasma markers
    lab = a.obs["annotation"].astype(str).to_numpy()
    assert (lab[truth != "Plasma"] == "Plasma").mean() < 0.02               # no over-calling from ambient Ig
    assert (lab[truth == "Plasma"] == "Plasma").mean() >= 0.80             # both JCHAIN and MZB1 must be detected (~0.9^2)
    assert "Plasma" not in rep["unreliable"]
    from bio.annotate import report_text
    assert "ambient immunoglobulin" in report_text(rep) and "IGKC, IGHG1" in report_text(rep)


# ----------------------------------------------------------------------------------------------- sample areas
COORDS = ('#Selection names:,"run1_crc_a_c, run1_crc_b_m"\n'
          '"#Areas (µm^2): 180000.00, 180000.00"\n'
          "Selection,X,Y\n"
          "run1_crc_a_c,-10,-10\nrun1_crc_a_c,300,-10\nrun1_crc_a_c,300,610\nrun1_crc_a_c,-10,610\n"
          "run1_crc_b_m,300,-10\nrun1_crc_b_m,610,-10\nrun1_crc_b_m,610,610\nrun1_crc_b_m,300,610\n")


def test_region_info():
    from bio.kit.xenium_cmd import region_info
    assert region_info("run1_crc_b_m") == ("run1_crc_b", "margin")
    assert region_info("run2_crc_e_c") == ("run2_crc_e", "core")
    assert region_info("run2_crc_f_n") == ("run2_crc_f", "normal")
    assert region_info("slide5_region1") == ("slide5_region1", "unknown")
    assert region_info("P7") == ("P7", "unknown")


def test_read_areas_format(tmp_path):
    from bio.xenium import read_areas
    f = tmp_path / "c.csv"
    f.write_text(COORDS + "run1_crc_a_c,0,0\nrun1_crc_a_c,5,0\nrun1_crc_a_c,5,5\n")   # a second polygon of selection a, not consecutive
    polys = read_areas(f)
    assert [n for n, _ in polys] == ["run1_crc_a_c", "run1_crc_b_m", "run1_crc_a_c"]
    assert polys[0][1].shape == (4, 2) and polys[1][1][2].tolist() == [610.0, 610.0] and polys[2][1].shape == (3, 2)
    f.write_text("Selection,X,Y\nA,1,2\nA,3,4\n")
    with pytest.raises(ValueError, match="fewer than 3 vertices"):
        read_areas(f)
    f.write_text("a,b,c\n1,2,3\n")
    with pytest.raises(ValueError, match="Selection,X,Y"):
        read_areas(f)
    f.write_text("Selection,X,Y\nA,1,x\n")
    with pytest.raises(ValueError, match="not numbers"):
        read_areas(f)


def test_point_in_polygon_backends_agree():
    from bio.xenium import points_in_polygon
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 2 * np.pi, 400))                               # a star-like concave polygon, hundreds of vertices
    r = 100 + 60 * np.sin(7 * t)
    verts = np.c_[200 + r * np.cos(t), 200 + r * np.sin(t)]
    pts = rng.uniform(0, 400, (20000, 2))
    a, b = points_in_polygon(pts, verts, "matplotlib"), points_in_polygon(pts, verts, "numpy")
    assert (a != b).sum() <= 2 and 0.1 < a.mean() < 0.6                         # only points a hair from an edge may differ
    assert not points_in_polygon(np.array([[500.0, 500.0]]), verts, "numpy").any()


def test_cli_areas_cut_one_slide(tmp_path, capsys):
    bundle, truth, xy, genes, neg = _make(tmp_path, outlier=True)
    coords = tmp_path / "slide_coordinates.csv"
    coords.write_text(COORDS)
    out = tmp_path / "kit"
    main(["xenium", str(bundle), "--areas", str(coords), "--dataset", "SYN_CUT", "--out", str(out)])
    cap = capsys.readouterr()
    assert "1 outside every selection (dropped)" in cap.out and "WARNING" in cap.err
    a = ad.read_h5ad(out / "SYN_CUT.h5ad")
    assert a.n_obs == 1999 and "cell00000-1_run1_crc_a_c" not in a.obs_names and "cell00000-1_run1_crc_b_m" not in a.obs_names
    assert set(a.obs["patient"]) == {"run1_crc_a_c", "run1_crc_b_m"} and (a.obs["patient"] == a.obs["sample"]).all()
    left = a.obsm["spatial"][:, 0] < 300
    assert (a.obs["patient"].to_numpy()[left] == "run1_crc_a_c").all() and (a.obs["patient"].to_numpy()[~left] == "run1_crc_b_m").all()
    assert (a.obs.loc[a.obs["patient"] == "run1_crc_a_c", "tissue"] == "core").all()
    assert (a.obs.loc[a.obs["patient"] == "run1_crc_b_m", "tissue"] == "margin").all()
    assert set(a.obs["patient_id"]) == {"run1_crc_a", "run1_crc_b"}
    assert a.obs["patient"].value_counts().sum() == 1999 and a.obs_names.is_unique
    rpt = (out / "SYN_CUT_annotation_report.txt").read_text()
    assert "SECTION run1_crc_a_c" in rpt and "1 outside every selection" in rpt and "patient_id" in rpt
    ytxt = (out / "SYN_CUT.yaml").read_text()
    assert "# obs columns: patient = full section name" in ytxt and "run1_crc_b_m -> run1_crc_b/margin" in ytxt
    assert "outside every selection" in ytxt
    m = manifest.load(out / "SYN_CUT.yaml")
    assert m["columns"]["patient"] == "patient" and m["columns"]["field"] == "sample"


def test_cli_areas_none_and_count_mismatch(tmp_path):
    bundle, *_ = _make(tmp_path, n=500)
    out = tmp_path / "o"
    main(["xenium", str(bundle), "--areas", "none", "--dataset", "SYN_NONE", "--out", str(out)])
    a = ad.read_h5ad(out / "SYN_NONE.h5ad")
    assert a.n_obs == 500 and set(a.obs["patient"]) == {"SYN_NONE"} and set(a.obs["tissue"]) == {"unknown"}
    with pytest.raises(SystemExit, match="one entry per bundle"):
        main(["xenium", str(bundle), "--areas", "none", "none", "--dataset", "X", "--out", str(out)])
    with pytest.raises(SystemExit, match="cannot use areas"):
        main(["xenium", str(bundle), "--areas", str(tmp_path / "missing.csv"), "--dataset", "X", "--out", str(out)])


def _deposit(root):
    root.mkdir()
    for reg, seed in (("slide1_region1", 1), ("run2_crc_e_c", 2)):
        (root / reg).mkdir()
        _make(root / reg, seed=seed, n=700)
    (root / "slide1_region1_coordinates.csv").write_text(COORDS)
    (root / "samplesheet.csv").write_text("slide,region,sample_areas\ncrca_slide1,slide1_region1,slide1_region1_coordinates.csv\n"
                                          "crca_slide3,run2_crc_e_c,none\n")
    return root


def test_cli_samplesheet_whole_deposit(tmp_path, capsys):
    dep = _deposit(tmp_path / "deposit")
    out = tmp_path / "kit"
    main(["xenium", "--samplesheet", str(dep / "samplesheet.csv"), "--dataset", "DEP", "--out", str(out)])
    a = ad.read_h5ad(out / "DEP.h5ad")
    assert a.n_obs == 1400                                                    # every synthetic cell lies inside a polygon
    assert a.obs["patient"].value_counts().index.isin(["run1_crc_a_c", "run1_crc_b_m", "run2_crc_e_c"]).all()
    assert set(a.obs["patient"]) == {"run1_crc_a_c", "run1_crc_b_m", "run2_crc_e_c"}
    assert a.obs.groupby("patient")["tissue"].first().to_dict() == {"run1_crc_a_c": "core", "run1_crc_b_m": "margin", "run2_crc_e_c": "core"}
    assert a.obs.groupby("patient")["patient_id"].first().to_dict() == {"run1_crc_a_c": "run1_crc_a", "run1_crc_b_m": "run1_crc_b", "run2_crc_e_c": "run2_crc_e"}
    # --regions restricts; unknown regions, missing folders and mixing modes are refused
    out2 = tmp_path / "kit2"
    main(["xenium", "--samplesheet", str(dep / "samplesheet.csv"), "--regions", "run2_crc_e_c", "--dataset", "DEP1", "--out", str(out2)])
    assert set(ad.read_h5ad(out2 / "DEP1.h5ad").obs["patient"]) == {"run2_crc_e_c"}
    with pytest.raises(SystemExit, match="not in the samplesheet"):
        main(["xenium", "--samplesheet", str(dep / "samplesheet.csv"), "--regions", "nope", "--dataset", "X", "--out", str(out2)])
    with pytest.raises(SystemExit, match="do not also give"):
        main(["xenium", str(dep / "run2_crc_e_c"), "--samplesheet", str(dep / "samplesheet.csv"), "--dataset", "X", "--out", str(out2)])
    (dep / "slide1_region1_coordinates.csv").unlink()
    with pytest.raises(SystemExit, match="coordinates file"):
        main(["xenium", "--samplesheet", str(dep / "samplesheet.csv"), "--dataset", "X", "--out", str(out2)])


# ----------------------------------------------------------------------------------------------- cluster method
def test_plasma_needs_both_markers(tmp_path):
    bundle, truth, *_ = _make(tmp_path, n=1500)
    a = read_bundle(bundle)
    X = a.X.toarray()
    only_one = np.where(truth == "Plasma")[0]
    X[np.ix_(only_one, [list(a.var_names).index("MZB1")])] = 0               # JCHAIN alone must not make a Plasma cell
    a.X = sparse.csr_matrix(X)
    assert (annotate(a).obs["annotation"].astype(str).to_numpy()[only_one] == "Plasma").sum() == 0


def test_read_bundle_loads_graphclust_only(tmp_path):
    bundle, truth, *_ = _make(tmp_path, clusters=True)
    a = read_bundle(bundle)
    csv = pd.read_csv(bundle / "analysis/clustering/gene_expression_graphclust/clusters.csv", dtype=str).set_index("Barcode")["Cluster"]
    assert a.obs["cluster_10x"].astype(str).tolist() == csv.loc[a.obs_names].tolist()   # aligned, graphclust not kmeans
    assert "clusters_file" in a.uns["xenium_bundle"]
    assert "cluster_10x" not in read_bundle(bundle, clusters="none").obs
    assert "cluster_10x" in read_bundle(bundle.parent).obs                     # parent folder finds it too
    other = tmp_path / "other.csv"
    pd.DataFrame({"Barcode": ["zzz-1"], "Cluster": ["1"]}).to_csv(other, index=False)
    with pytest.raises(ValueError, match="none of its"):
        read_bundle(bundle, clusters=other)
    with pytest.raises(FileNotFoundError, match="clusters file"):
        read_bundle(bundle, clusters=tmp_path / "missing.csv")
    plain, *_ = _make(tmp_path / "x") if (tmp_path / "x").mkdir() is None else (None,)
    assert "cluster_10x" not in read_bundle(plain).obs


def test_cluster_method_recovers_planted_and_mixed_is_ambiguous(tmp_path):
    bundle, truth, *_ = _make(tmp_path, n=3000, clusters=True)
    a = annotate(read_bundle(bundle))                                          # auto -> cluster
    P = json.loads(a.uns["annotation_params"])
    assert P["method"] == "cluster" and len(P["clusters"]) == 16
    lab = a.obs["annotation"].astype(str).to_numpy()
    cid = a.obs["cluster_10x"].astype(str).to_numpy()
    mixed = cid == MIXED_ID
    assert mixed.sum() == 60 and (lab[mixed] == "Ambiguous").all()
    ok = ~mixed & (a.obs["annotation"].astype(str) != "LowCount").to_numpy()
    sub = np.isin(truth, ["CD8T", "CD4T", "Treg"])
    # T subsets are per cell now: a T cell with no subset evidence (e.g. zero CD4 counts) is "T", not a wrong subset
    assert (lab[ok & ~sub] == truth[ok & ~sub]).mean() >= 0.95 and (lab[ok] == truth[ok]).mean() >= 0.90
    for t in PLANTED:                                                           # every planted class, T subsets included
        m = ok & (truth == t)
        if t in ("CD8T", "CD4T", "Treg"):
            assert (lab[m] == t).mean() >= (0.60 if t == "Treg" else 0.80) and np.isin(lab[m], [t, "T"]).mean() >= (0.65 if t == "Treg" else 0.88), t   # Treg needs FOXP3 >= 2 (fixture: Poisson(3))
        else:
            assert (lab[m] == t).mean() >= 0.95, t
    tab = {c["cluster"]: c for c in P["clusters"]}
    assert tab[MIXED_ID]["label"] == "Ambiguous" and tab[MIXED_ID]["note"] in ("tie", "no lineage fired")
    assert tab["1"]["label"] == "Epithelial" and tab["1"]["best_scaled"] == 1.0 and len(tab["1"]["top3"]) == 3
    rep = annotate_report(a)
    assert rep["method"] == "cluster" and rep["unreliable"] == {}
    assert len(rep["clusters"]) == 16
    from bio.annotate import report_text
    txt = report_text(rep)
    assert "Method: cluster" in txt and "Clusters (n cells" in txt and "Marker leakage" in txt and "scaled" in txt and ("[tie]" in txt or "[no lineage fired]" in txt)
    # method switches
    b = annotate(read_bundle(bundle), method="cell")
    assert json.loads(b.uns["annotation_params"])["method"] == "cell"
    with pytest.raises(ValueError, match="cluster_10x"):
        annotate(read_bundle(bundle, clusters="none"), method="cluster")
    assert json.loads(annotate(read_bundle(bundle, clusters="none")).uns["annotation_params"])["method"] == "cell"
    with pytest.raises(ValueError, match="method must be"):
        annotate(read_bundle(bundle), method="nope")


def test_cluster_unreliable_flag_uses_scaled_gap():
    from bio.annotate import _cluster_unreliable
    c = lambda lin, b, s2: {"lineage": lin, "best_scaled": b, "second_scaled": s2}   # noqa: E731
    flagged = _cluster_unreliable([c("Neutrophil", 0.6, 0.5), c("Neutrophil", 0.55, 0.45), c("Neutrophil", 1.0, 0.1),
                                   c("B", 1.0, 0.1), c("B", 0.9, 0.2), c("", 0.9, 0.9)])
    assert set(flagged) == {"Neutrophil"} and "median best-minus-runner-up 0.10" in flagged["Neutrophil"][0]
    assert _cluster_unreliable([c("T", 0.5, 0.37)]) == {"T": ["median best-minus-runner-up 0.13 over 1 cluster(s), below 0.25"]}


def _dominant_proportion(share):
    p = np.full(len(PLANTED), (1 - share) / (len(PLANTED) - 1))
    p[PLANTED.index("Epithelial")] = share
    return p


@pytest.mark.parametrize("share", [0.80, 0.96])
def test_cluster_method_labels_a_dominant_lineage(tmp_path, share):
    n = 5000 if share < 0.9 else 10000
    bundle, truth, *_ = _make(tmp_path, n=n, clusters="pure", proportion=_dominant_proportion(share))
    a = annotate(read_bundle(bundle))
    P = json.loads(a.uns["annotation_params"])
    tab = {c["cluster"]: c for c in P["clusters"]}
    assert tab["1"]["label"] == "Epithelial" and tab["1"]["best_scaled"] == 1.0         # the old ratio rule gave Ambiguous here
    lab = a.obs["annotation"].astype(str).to_numpy()
    ok = lab != "LowCount"
    assert (lab[ok & (truth == "Epithelial")] == "Epithelial").mean() >= 0.99
    called = [t for t in PLANTED if t != "Epithelial" and (truth == t).sum() >= 15]
    assert len(called) >= 8
    for t in called:                                                                    # the rare lineages are still found
        if t in ("CD8T", "CD4T", "Treg"):                                               # subsets are per cell: the T lineage must be found
            assert np.isin(lab[ok & (truth == t)], ["CD8T", "CD4T", "Treg", "T"]).mean() >= 0.90, t
        else:
            assert (lab[ok & (truth == t)] == t).mean() >= 0.90, t
    assert "Epithelial" not in P["absent"] and (share > 0.9 or P["absent"] == [])
    assert annotate_report(a)["unreliable"] == {}


def test_cluster_method_does_not_call_an_absent_lineage(tmp_path):
    pr = PROPORTION.copy()
    pr[PLANTED.index("Neutrophil")] = 0                                                  # no neutrophils in this section
    bundle, truth, *_ = _make(tmp_path, n=3000, clusters="pure", proportion=pr)
    a = annotate(read_bundle(bundle))
    P = json.loads(a.uns["annotation_params"])
    assert "Neutrophil" in P["absent"]
    assert "Neutrophil" not in set(a.obs["annotation"].astype(str))
    assert "Neutrophil" not in {c["label"] for c in P["clusters"]}                      # noise in other clusters is not a call
    from bio.annotate import report_text
    txt = report_text(annotate_report(a))
    assert "Absent from the section" in txt and "Neutrophil" in txt.split("Absent from the section")[1].split("\n")[0]
    lab = a.obs["annotation"].astype(str).to_numpy()
    assert (lab[(truth == "Macrophage")] == "Macrophage").mean() >= 0.90                # others are unaffected


def test_cluster_cells_inherit_lowcount_and_missing_ids(tmp_path):
    bundle, truth, *_ = _make(tmp_path, n=1500, clusters=True)
    a = read_bundle(bundle)
    X = a.X.toarray()
    X[:10, :] = 0
    a.X = sparse.csr_matrix(X)
    cl = a.obs["cluster_10x"].astype(object)
    cl.iloc[10:20] = np.nan                                                    # ten cells the clusters file lacks
    a.obs["cluster_10x"] = pd.Categorical(cl)
    a = annotate(a)
    lab = a.obs["annotation"].astype(str)
    assert (lab.iloc[:10] == "LowCount").all()
    assert (lab.iloc[10:20] == "Unclustered").all() and (a.obs["annotation_note"].astype(str).iloc[10:20] == "no cluster id").all()
    rep = annotate_report(a)
    assert rep["n_no_cluster"] == 10 and rep["labels"]["Unclustered"] == 10 and "Unclustered" not in rep["specificity"].columns
    from bio.annotate import report_text
    assert "Unclustered = 10 cells" in report_text(rep)


def test_clusters_file_that_omits_barcodes_end_to_end(tmp_path, capsys):
    bundle, truth, *_ = _make(tmp_path, n=2000, clusters=True)
    f = bundle / "analysis/clustering/gene_expression_graphclust/clusters.csv"
    full = pd.read_csv(f, dtype=str)
    keep = full.drop(index=range(0, 2000, 20))                                 # 5% of the cells have no cluster
    keep.to_csv(f, index=False)
    a = read_bundle(bundle)
    assert a.obs["cluster_10x"].isna().sum() == 100
    assert json.loads(a.uns["xenium_bundle"])["clusters_cells_without_cluster"] == 100
    out = tmp_path / "o"
    main(["xenium", str(bundle), "--dataset", "OM", "--out", str(out)])
    h = ad.read_h5ad(out / "OM.h5ad")
    lab = h.obs["annotation"].astype(str)
    gone = ~h.obs_names.isin(keep["Barcode"])
    assert gone.sum() == 100 and (lab[gone] == "Unclustered").all()            # no per-cell fallback
    assert (lab[~gone] == "Unclustered").sum() == 0
    assert "Unclustered = 100 cells (5.0%)" in (out / "OM_annotation_report.txt").read_text()
    assert manifest.load(out / "OM.yaml")["roles"]["Unclustered"] == "other"
    # the same through a cohort: cut slide + clusters of the whole slide
    coords = tmp_path / "c.csv"
    coords.write_text(COORDS)
    out2 = tmp_path / "o2"
    main(["xenium", str(bundle), "--areas", str(coords), "--dataset", "OM2", "--out", str(out2)])
    h2 = ad.read_h5ad(out2 / "OM2.h5ad")
    assert (h2.obs["annotation"].astype(str) == "Unclustered").sum() == 100


def test_cli_clusters_default_none_and_explicit_with_areas(tmp_path, capsys):
    bundle, truth, *_ = _make(tmp_path, clusters=True)
    coords = tmp_path / "slide_coordinates.csv"
    coords.write_text(COORDS)
    out = tmp_path / "o1"
    main(["xenium", str(bundle), "--dataset", "CL", "--out", str(out)])           # auto-loaded
    rpt = (out / "CL_annotation_report.txt").read_text()
    assert "Method: cluster" in rpt and "Clusters (n cells" in rpt
    assert ad.read_h5ad(out / "CL.h5ad").obs["cluster_10x"].notna().all()
    out = tmp_path / "o2"
    main(["xenium", str(bundle), "--clusters", "none", "--dataset", "CL", "--out", str(out)])
    rpt = (out / "CL_annotation_report.txt").read_text()
    assert "Method: cell" in rpt and "Clusters (n cells" not in rpt
    assert "cluster_10x" not in ad.read_h5ad(out / "CL.h5ad").obs
    # explicit file elsewhere + areas
    ext = tmp_path / "elsewhere.csv"
    csv = pd.read_csv(bundle / "analysis/clustering/gene_expression_graphclust/clusters.csv")
    csv.to_csv(ext, index=False)
    import shutil
    shutil.rmtree(bundle / "analysis")
    out = tmp_path / "o3"
    main(["xenium", str(bundle), "--areas", str(coords), "--clusters", str(ext), "--dataset", "CLA", "--out", str(out)])
    a = ad.read_h5ad(out / "CLA.h5ad")
    assert set(a.obs["patient"]) == {"run1_crc_a_c", "run1_crc_b_m"}
    rpt = (out / "CLA_annotation_report.txt").read_text()
    assert rpt.count("Method: cluster") >= 3 and "SECTION run1_crc_a_c" in rpt          # two sections + cohort
    lab = a.obs["annotation"].astype(str).to_numpy()
    cid = a.obs["cluster_10x"].astype(str).to_numpy()
    ok = (cid != MIXED_ID) & (lab != "LowCount")
    tr = truth[pd.Index(ad.read_h5ad(out / "CLA.h5ad").obs_names.str.replace(r"_run1_crc_[ab]_[cm]$", "", regex=True))
               .map({f"cell{i:05d}-1": i for i in range(len(truth))})]
    to_lin = lambda x: np.where(np.isin(x, ["CD8T", "CD4T", "Treg"]), "T", x)           # T subsets are per cell (checked elsewhere)
    assert (to_lin(lab[ok]) == to_lin(tr[ok])).mean() >= 0.90
    # one --clusters entry per bundle
    with pytest.raises(SystemExit, match="one entry per bundle"):
        main(["xenium", str(bundle), "--clusters", "none", "none", "--dataset", "X", "--out", str(out)])


def test_cli_jobs_matches_serial(tmp_path, capsys):
    for pat, seed in (("P01", 1), ("P02", 2)):
        (tmp_path / pat).mkdir()
        _make(tmp_path / pat, seed=seed, n=600, clusters=True)
    b = [str(tmp_path / "P01" / "outs"), str(tmp_path / "P02" / "outs")]
    outs = {}
    for jobs in (1, 2):
        outs[jobs] = tmp_path / f"kit{jobs}"
        main(["xenium", *b, "--dataset", "J", "--out", str(outs[jobs]), "--jobs", str(jobs)])
        err = capsys.readouterr().err
        assert err.count("done ") == 2 and "cells" in err and " s [" in err        # one progress line per bundle
    a1, a2 = ad.read_h5ad(outs[1] / "J.h5ad"), ad.read_h5ad(outs[2] / "J.h5ad")
    pd.testing.assert_frame_equal(a1.obs, a2.obs)
    assert (outs[1] / "J_annotation_report.txt").read_text() == (outs[2] / "J_annotation_report.txt").read_text()
    assert (outs[1] / "J.yaml").read_text().replace("kit1", "kit") == (outs[2] / "J.yaml").read_text().replace("kit2", "kit")
    assert (a1.X != a2.X).nnz == 0
    with pytest.raises(SystemExit, match="--jobs must be"):
        main(["xenium", *b, "--dataset", "J", "--out", str(outs[1]), "--jobs", "0"])


def test_jobs_capped_at_cpu_count(tmp_path, capsys, monkeypatch):
    bundle, *_ = _make(tmp_path, n=500)
    monkeypatch.setattr(os, "cpu_count", lambda: 1)
    main(["xenium", str(bundle), "--dataset", "C", "--out", str(tmp_path / "o"), "--jobs", "8"])
    assert "capped at the CPU count (1)" in capsys.readouterr().err


# ----------------------------------------------------------------------------- regression, sections, parts, resume
def test_annotate_and_report_outputs_are_unchanged(tmp_path):
    """Stored digests of the labelling output (made with the per-cell-loop implementation before it was vectorised):
    the sparse-product scoring and the groupby-style leakage table must not change a single value."""
    import hashlib
    from bio.annotate import annotate_report, report_text
    _make(tmp_path, seed=11, n=3000, clusters=True)
    a0 = read_bundle(tmp_path / "outs")
    h = lambda x: hashlib.sha256(x.encode()).hexdigest()[:16]
    # Updated for the robust reference maximum + per-cell T subsets + the new "T" label: the cell-method labels and the
    # fired matrix are unchanged (Treg rule v2: FOXP3 >= 2 and no CD8 changed the cluster labels and the rule text in params); its params (new rule keys) and report (a "T" row in the label counts) differ; the cluster labels differ. Cell method (rule now shared with the cluster method: Treg FOXP3 >= 2 and no CD8, no evidence -> T): its labels, params (t_default_cd4 counts T) and report (rule text, T row) differ; fired is unchanged.
    expect = {"cell": {"fired": "c04ee0e8b5c22d51", "labels": "41aaa5d574432f01", "params": "534fae21d4201747", "report": "7802701856d3f434"},
              "cluster": {"fired": "c04ee0e8b5c22d51", "labels": "51497821bc9803fd", "params": "47545c5adf7bbb9e", "report": "c02077d4baf63f22"}}
    for meth, exp in expect.items():
        a = a0.copy()
        annotate(a, method=meth)
        rep = annotate_report(a)
        got = {"fired": hashlib.sha256(np.ascontiguousarray(a.obsm["lineage_fired"]).tobytes()).hexdigest()[:16],
               "labels": h(",".join(a.obs["annotation"].astype(str))), "params": h(a.uns["annotation_params"]),
               "report": h(report_text(rep, "X"))}
        assert got == exp, meth
    # Mast is 109 (was 49 + 60 Ambiguous): the true Mast cluster has 49 cells, one under the size rule, so the 60-cell
    # Plasma/Mast mixed cluster became the Mast reference and is now called Mast (Plasma scaled 0.53)
    assert rep["labels"]["Epithelial"] == 760 and rep["labels"]["Mast"] == 109 and rep["labels"]["T"] == 14


def test_points_in_polygon_numpy_matches_a_plain_loop():
    from bio.xenium import _contains_numpy
    rng = np.random.default_rng(8)
    t = np.sort(rng.uniform(0, 2 * np.pi, 60))
    r = 50 + 20 * np.sin(5 * t)
    verts = np.c_[100 + r * np.cos(t), 100 + r * np.sin(t)]
    pts = rng.uniform(0, 200, (500, 2))

    def one(p):
        ins, n = False, len(verts)
        for i in range(n):
            (x1, y1), (x2, y2) = verts[i], verts[(i + 1) % n]
            if y1 != y2 and ((y1 > p[1]) != (y2 > p[1])) and p[0] < x1 + (p[1] - y1) * (x2 - x1) / (y2 - y1):
                ins = not ins
        return ins
    import bio.xenium as bx
    old, bx._PIP_CELLS = bx._PIP_CELLS, 700            # force several chunks
    try:
        got = _contains_numpy(pts, verts)
    finally:
        bx._PIP_CELLS = old
    assert got.tolist() == [one(p) for p in pts]


def _deposit_clusters(root):
    """One slide cut in two sections (900 cells) and one whole region (500 cells), both with 10x clusters."""
    root.mkdir()
    for reg, seed, n in (("slide1_region1", 1, 900), ("run2_crc_e_c", 2, 500)):
        (root / reg).mkdir()
        _make(root / reg, seed=seed, n=n, clusters=True)
    (root / "slide1_region1_coordinates.csv").write_text(COORDS)
    (root / "samplesheet.csv").write_text("slide,region,sample_areas\ncrca_slide1,slide1_region1,slide1_region1_coordinates.csv\n"
                                          "crca_slide3,run2_crc_e_c,none\n")
    return root / "samplesheet.csv"


def _same_cohort(o1, o2, name, tag1="k1", tag2="k2"):
    a1, a2 = ad.read_h5ad(o1 / f"{name}.h5ad"), ad.read_h5ad(o2 / f"{name}.h5ad")
    pd.testing.assert_frame_equal(a1.obs, a2.obs)
    assert (a1.X != a2.X).nnz == 0 and list(a1.var_names) == list(a2.var_names)
    np.testing.assert_array_equal(a1.obsm["spatial"], a2.obsm["spatial"])
    assert (o1 / f"{name}_annotation_report.txt").read_text() == (o2 / f"{name}_annotation_report.txt").read_text()
    assert (o1 / f"{name}.yaml").read_text().replace(tag1, "kit") == (o2 / f"{name}.yaml").read_text().replace(tag2, "kit")


def test_cli_jobs_sections_largest_first_and_identical(tmp_path, capsys):
    import re
    sheet = _deposit_clusters(tmp_path / "dep")
    outs, lines = {}, {}
    for jobs in (1, 3):
        outs[jobs] = tmp_path / f"k{jobs}"
        main(["xenium", "--samplesheet", str(sheet), "--dataset", "D", "--out", str(outs[jobs]), "--jobs", str(jobs)])
        lines[jobs] = re.findall(r"^done (\S+): (\d+) cells, \d+ s \[(\d+)/(\d+)\]$", capsys.readouterr().err, re.M)
        assert len(lines[jobs]) == 3 and {int(x[3]) for x in lines[jobs]} == {3}      # one line per section
        assert sorted(int(x[2]) for x in lines[jobs]) == [1, 2, 3]
    cells = [int(x[1]) for x in lines[1]]
    assert cells == sorted(cells, reverse=True)                                      # serial: largest section first
    assert set(x[0] for x in lines[3]) == {"run1_crc_a_c", "run1_crc_b_m", "run2_crc_e_c"}
    _same_cohort(outs[1], outs[3], "D", "k1", "k3")
    a = ad.read_h5ad(outs[1] / "D.h5ad")
    assert list(dict.fromkeys(a.obs["patient"])) == ["run1_crc_a_c", "run1_crc_b_m", "run2_crc_e_c"]   # samplesheet order
    assert not (outs[1] / "parts" / "D" / "raw").exists()                              # scratch files cleaned on success
    assert sorted(f.name for f in (outs[1] / "parts" / "D").glob("*.h5ad")) == \
        ["run1_crc_a_c.h5ad", "run1_crc_b_m.h5ad", "run2_crc_e_c.h5ad"]
    assert (outs[1] / "parts" / "D" / "run2_crc_e_c.txt").read_text().startswith("# part-key ")


def test_default_jobs_is_cpus_minus_one_and_stage_a_is_capped(tmp_path, monkeypatch):
    import bio.kit.xenium_cmd as xc
    sheet = _deposit_clusters(tmp_path / "dep")
    seen = []
    real = xc._map
    monkeypatch.setattr(xc, "_map", lambda fn, items, workers, cb: (seen.append((fn.__name__, workers)), real(fn, items, 1, cb))[1])
    monkeypatch.setattr(os, "cpu_count", lambda: 10)
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "D", "--out", str(tmp_path / "o")])
    assert seen == [("_cut_region", 3), ("_annotate_section", 9)]
    seen.clear()
    monkeypatch.setattr(os, "cpu_count", lambda: 1)
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "D", "--out", str(tmp_path / "o1")])
    assert seen == [("_cut_region", 1), ("_annotate_section", 1)]


def test_resume_is_identical_and_picks_up_after_an_interrupted_run(tmp_path, capsys):
    import re
    import shutil
    sheet = _deposit_clusters(tmp_path / "dep")
    out = tmp_path / "kit"
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "R", "--out", str(out)])
    capsys.readouterr()
    first = tmp_path / "first"
    first.mkdir()
    for f in ("R.h5ad", "R_annotation_report.txt", "R.yaml"):
        shutil.copy(out / f, first / f)
    # all parts complete: delete the cohort outputs, resume -> nothing is recomputed, same cohort
    for f in ("R.h5ad", "R_annotation_report.txt", "R.yaml"):
        (out / f).unlink()
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "R", "--out", str(out), "--resume"])
    err = capsys.readouterr().err
    assert len(re.findall(r"^kept \S+ \(from parts\)", err, re.M)) == 3 and "done " not in err and err.count("skip cutting") == 2
    shutil.copy(out / "R.yaml", first / "R2.yaml")
    assert (first / "R2.yaml").read_text() == (first / "R.yaml").read_text()
    _same_cohort(first, out, "R", "first", "kit")
    # interrupted after the first section: one part exists, the rest are missing (and a half-written one is not a part)
    parts = out / "parts" / "R"
    for f in list(parts.iterdir()):
        if f.stem != "run1_crc_a_c" and f.suffix in (".h5ad", ".txt"):
            f.unlink()
    (parts / "run1_crc_b_m.h5ad").write_bytes(b"half written")           # no .txt: not complete
    (parts / "run2_crc_e_c.tmp.h5ad").write_bytes(b"killed mid-write")
    for f in ("R.h5ad", "R_annotation_report.txt", "R.yaml"):
        (out / f).unlink()
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "R", "--out", str(out), "--resume", "--jobs", "2"])
    err = capsys.readouterr().err
    assert err.count("kept run1_crc_a_c (from parts)") == 1 and err.count("done ") == 2 and "skip cutting" not in err
    _same_cohort(first, out, "R", "first", "kit")
    # changed settings do not reuse the old parts
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "R", "--out", str(out), "--resume", "--min-counts", "30"])
    err = capsys.readouterr().err
    assert not re.findall(r"^kept ", err, re.M) and err.count("done ") == 3
    # --no-keep-parts removes the parts after a good run
    out2 = tmp_path / "kit2"
    main(["xenium", "--samplesheet", str(sheet), "--dataset", "R", "--out", str(out2), "--no-keep-parts"])
    assert (out2 / "R.h5ad").exists() and not (out2 / "parts").exists()


def test_failed_run_keeps_parts_and_says_where(tmp_path, capsys):
    sheet = _deposit_clusters(tmp_path / "dep")
    (tmp_path / "dep" / "run2_crc_e_c" / "outs" / "cells.csv.gz").unlink()
    out = tmp_path / "kit"
    with pytest.raises(SystemExit) as e:
        main(["xenium", "--samplesheet", str(sheet), "--dataset", "F", "--out", str(out)])
    assert "cannot read" in str(e.value) and str(out / "parts" / "F") in str(e.value) and "--resume" in str(e.value)


# ------------------------------------------------------------------ cluster-level T subsets and NK vs T (synthetic sections)
def _t_section(spec, seed=1, with_foxp3=True, extra_genes=(), n_epi=3000):
    """spec: {cluster id: (n cells, {gene: (probability detected, mean count when detected)})}. A big epithelial cluster
    dilutes the all-cell means as in real sections. Returns the annotated AnnData; clusters are the spec keys."""
    rng = np.random.default_rng(seed)
    genes = ["EPCAM", "KRT8", "KRT18", "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4",
             "NCAM1", "KLRD1", "KLRF1", "GNLY", "NCR1", "ACTB", "GAPDH"]
    if not with_foxp3:
        genes.remove("FOXP3")
    genes += list(extra_genes)
    spec = {**({"epi": (n_epi, {"EPCAM": (.9, 4), "KRT8": (.9, 4), "KRT18": (.9, 4)})} if n_epi else {}), **spec}
    blocks, cl = [], []
    for name, (n, expr) in spec.items():
        M = rng.poisson(0.004, (n, len(genes))).astype(np.float32)
        M[:, genes.index("ACTB")] += rng.poisson(12, n)
        M[:, genes.index("GAPDH")] += rng.poisson(12, n)
        for g, (p, mu) in expr.items():
            if g in genes:
                on = rng.random(n) < p
                M[:, genes.index(g)] += np.where(on, 1 + rng.poisson(mu, n), 0)
        blocks.append(M)
        cl += [name] * n
    X = np.vstack(blocks)
    a = ad.AnnData(sparse.csr_matrix(X), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cluster_10x"] = pd.Categorical(cl)
    return annotate(a, method="cluster")


def _lab(a):
    return {c["cluster"]: c["label"] for c in json.loads(a.uns["annotation_params"])["clusters"]}


CD3 = {"CD3D": (.8, 3), "CD3E": (.8, 3), "CD3G": (.7, 2), "CD2": (.6, 2)}


def _split(a):
    """{cluster: {label: share of the cluster's cells}} from the cell labels."""
    df = pd.DataFrame({"c": a.obs["cluster_10x"].astype(str), "l": a.obs["annotation"].astype(str)})
    return {c: g["l"].value_counts(normalize=True).to_dict() for c, g in df.groupby("c")}


def test_cluster_t_subsets_are_decided_per_cell_inside_t_clusters():
    spec = {
        "cd4_sparse": (300, {**CD3, "CD4": (.7, 3), "FOXP3": (.05, 2)}),             # 5% of cells carry FOXP3
        "treg": (300, {**CD3, "CD4": (.7, 3), "FOXP3": (.6, 3), "IL2RA": (.5, 2)}),  # 60%
        "cd8": (300, {**CD3, "CD8A": (.8, 4), "CD8B": (.7, 3)}),
        "cd4_b": (300, {**CD3, "CD4": (.7, 3)}),
    }
    a = _t_section(spec)
    assert set(_lab(a).values()) >= {"T", "Epithelial"} and {_lab(a)[k] for k in spec} == {"T"}   # the cluster only says T
    sp = _split(a)
    assert sp["cd4_sparse"]["CD4T"] > 0.6 and sp["cd4_sparse"].get("Treg", 0) < 0.15
    assert 0.5 < sp["treg"]["Treg"] < 0.7
    assert sp["cd8"]["CD8T"] > 0.85 and sp["cd8"].get("Treg", 0) < 0.03 and sp["cd8"].get("CD4T", 0) < 0.03
    assert sp["cd4_b"]["CD4T"] > 0.6 and sp["cd4_b"].get("Treg", 0) < 0.03
    assert sp["epi"]["Epithelial"] > 0.95
    P = json.loads(a.uns["annotation_params"])
    assert P["t_no_subset"] == sum(P["clusters"][i]["t_split"]["T"] for i in range(len(P["clusters"])) if "t_split" in P["clusters"][i])
    assert P["n_t_cluster"] == P["n_t"] > 0 and "FOXP3 count >= 2" in P["cluster_rules"]["t_subset"]


def test_one_mixed_t_cluster_yields_all_four_labels_in_the_planted_proportions():
    """One T cluster holding planted CD8, CD4, Treg and marker-less T cells: each cell gets its own label."""
    rng = np.random.default_rng(5)
    genes = ["EPCAM", "KRT8", "KRT18", "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4", "ACTB", "GAPDH"]
    g = {x: i for i, x in enumerate(genes)}
    n_epi, plan = 3000, {"CD8T": 400, "CD4T": 300, "Treg": 100, "T": 200}
    rows, truth = [], []
    M = rng.poisson(0.004, (n_epi, len(genes))).astype(np.float32)
    M[:, [g["ACTB"], g["GAPDH"]]] += rng.poisson(12, (n_epi, 2))
    M[:, [g["EPCAM"], g["KRT8"], g["KRT18"]]] += rng.poisson(4, (n_epi, 3))
    rows.append(M)
    truth += ["epi"] * n_epi
    for lab, n in plan.items():
        M = rng.poisson(0.004, (n, len(genes))).astype(np.float32)
        M[:, [g["ACTB"], g["GAPDH"]]] += rng.poisson(12, (n, 2))
        M[:, [g["CD3D"], g["CD3E"], g["CD3G"], g["CD2"]]] += 2 + rng.poisson(2, (n, 4))   # every T cell is clearly CD3+
        if lab == "CD8T":
            M[:, g["CD8A"]] += 2 + rng.poisson(2, n)
        if lab in ("CD4T", "Treg"):
            M[:, g["CD4"]] += 2 + rng.poisson(2, n)
        if lab == "Treg":
            M[:, g["FOXP3"]] += 2 + rng.poisson(2, n)
        rows.append(M)
        truth += [lab] * n
    a = ad.AnnData(sparse.csr_matrix(np.vstack(rows)), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cluster_10x"] = pd.Categorical(["epi"] * n_epi + ["tcl"] * sum(plan.values()))        # ONE T cluster
    a = annotate(a, method="cluster")
    assert _lab(a)["tcl"] == "T"
    truth = np.array(truth)
    lab = a.obs["annotation"].astype(str).to_numpy()
    t = truth != "epi"
    for k, n in plan.items():
        got = (lab[t] == k).sum()
        assert abs(got - n) <= 0.10 * n, (k, got, n)
        assert (lab[truth == k] == k).mean() >= 0.9, k
    assert (lab[~t] == "Epithelial").mean() > 0.98
    P = json.loads(a.uns["annotation_params"])
    assert P["n_t"] == 1000 and 190 <= P["t_no_subset"] <= 215
    from bio.annotate import report_text
    txt = report_text(annotate_report(a))
    assert "have no subset evidence" in txt and "labelled T" in txt


def test_cluster_t_subsets_without_foxp3_use_il2ra_and_ctla4():
    spec = {"treg": (300, {**CD3, "CD4": (.7, 3), "IL2RA": (.7, 4), "CTLA4": (.7, 4)}),
            "il2ra_only": (300, {**CD3, "CD4": (.7, 3), "IL2RA": (.7, 4)}),
            "cd4": (300, {**CD3, "CD4": (.7, 3), "IL2RA": (.05, 1), "CTLA4": (.05, 1)})}
    sp = _split(_t_section(spec, with_foxp3=False))
    assert sp["treg"]["Treg"] > 0.4 and sp["il2ra_only"].get("Treg", 0) < 0.05 and sp["il2ra_only"]["CD4T"] > 0.6
    assert sp["cd4"].get("Treg", 0) < 0.03 and sp["cd4"]["CD4T"] > 0.6


def test_cluster_nk_vs_cytotoxic_t():
    nk_markers = {"KLRD1": (.8, 4), "KLRF1": (.7, 3), "NCAM1": (.6, 3), "GNLY": (.9, 6)}
    spec = {
        "cd8_cyto": (300, {**CD3, "CD8A": (.8, 4), "CD8B": (.7, 3), "GNLY": (.9, 6), "KLRD1": (.2, 1)}),   # (d)
        "nk": (300, nk_markers),                                                                          # (e) no CD3
        "cd4": (300, {**CD3, "CD4": (.7, 3)}),
    }
    a = _t_section(spec)
    lab, sp = _lab(a), _split(a)
    assert lab["cd8_cyto"] == "T" and sp["cd8_cyto"]["CD8T"] > 0.85 and lab["nk"] == "NK" and lab["cd4"] == "T" and sp["cd4"]["CD4T"] > 0.6
    # GNLY alone is not enough for NK
    lab = _lab(_t_section({**spec, "gnly_only": (300, {"GNLY": (.9, 6)})}))
    assert lab["gnly_only"] != "NK"


def test_tiny_outlier_cluster_does_not_become_the_reference():
    """A 3-cell cluster with 10x the EPCAM of the real epithelial clusters used to set the 1.0 reference, so every real
    epithelial cluster scaled below 0.4 and went Ambiguous."""
    epi = {"EPCAM": (.9, 4), "KRT8": (.9, 4), "KRT18": (.9, 4)}
    spec = {"epi_b": (800, epi), "epi_c": (500, epi), "t": (300, {**CD3, "CD4": (.7, 3)}),
            "outlier": (3, {"EPCAM": (1.0, 40), "KRT8": (1.0, 40), "KRT18": (1.0, 40)})}
    a = _t_section(spec)
    lab = _lab(a)
    assert lab["epi"] == lab["epi_b"] == lab["epi_c"] == "Epithelial", lab
    assert _row(a, "epi")["best_scaled"] >= 0.9 and _row(a, "outlier")["label"] == "Epithelial"
    assert _row(a, "outlier")["best_scaled"] > 1.0                                          # not capped at 1.0
    ref = json.loads(a.uns["annotation_params"])["reference"]
    assert ref["min_n"] == 50.0 and ref["n_reference_clusters"] == 4 and ref["n_clusters"] == 5 and not ref["fallback_all"]
    # the same section with the old reference (all clusters) loses the epithelium
    import bio.annotate as an
    old, an.CL_REF_MIN_CELLS, an.CL_REF_MIN_FRAC = (an.CL_REF_MIN_CELLS, an.CL_REF_MIN_FRAC), 0, 0
    try:
        assert _lab(_t_section(spec))["epi"] != "Epithelial"
    finally:
        an.CL_REF_MIN_CELLS, an.CL_REF_MIN_FRAC = old
    # no cluster meets the size rule: fall back to all clusters
    rng_spec = {"epi": (40, epi), "t": (30, {**CD3, "CD4": (.7, 3)})}
    a = _t_section(rng_spec, n_epi=0)
    assert _lab(a)["epi"] == "Epithelial" and _lab(a)["t"] == "T"
    ref = json.loads(a.uns["annotation_params"])["reference"]
    assert ref["fallback_all"] and ref["n_reference_clusters"] == 2
    # a lineage that only forms small clusters is still called (reference falls back to all clusters for it)
    a = _t_section({"epi_b": (800, epi), "neut_small": (30, {"CSF3R": (.9, 6), "FCGR3B": (.9, 6), "CXCR2": (.9, 6)})},
                   extra_genes=["CSF3R", "FCGR3B", "CXCR2"])
    assert _lab(a)["neut_small"] == "Neutrophil" and "Neutrophil" in json.loads(a.uns["annotation_params"])["reference"]["fallback_lineages"]


def test_marteau_style_manifest_maps_t_to_role_t(tmp_path):
    from bio.annotate import LABELS
    from bio.kit.xenium_cmd import ROLES
    assert "T" in LABELS and ROLES["T"] == "T"                                       # plain T: subset unknown
    m = manifest.load("kit_manifests/MARTEAU_CRC_IMC_core.yaml")
    assert m["roles"]["Tcells_undefined"] == "T"                                     # the Marteau IMC manifest uses the same plain role
    assert manifest.problems(m) == []
    assert manifest.problems({**m, "roles": {**m["roles"], "T": "T"}}) == []        # the xenium kit's new "T" label as a raw type


def _row(a, cl):
    return {c["cluster"]: c for c in json.loads(a.uns["annotation_params"])["clusters"]}[cl]


def test_cluster_two_marker_lineage_called_at_low_fraction_of_max():
    # Plasma has exactly 2 panel markers (JCHAIN, MZB1): both need >= 0.3 (not 0.5) of their largest cluster mean, and
    # the cluster only needs scaled >= 0.4. This one sits at ~0.4-0.5 of the strongest plasma cluster.
    spec = {"plasma_hi": (300, {"JCHAIN": (.9, 9), "MZB1": (.9, 9)}),
            "plasma_lo": (300, {"JCHAIN": (.9, 2), "MZB1": (.9, 2.5)})}
    a = _t_section(spec, extra_genes=["JCHAIN", "MZB1"])
    assert _lab(a)["plasma_hi"] == "Plasma"
    r = _row(a, "plasma_lo")
    assert 0.4 <= r["best_scaled"] < 0.5 and r["label"] == "Plasma", r
    # below 0.3 x of the max for one of the two markers: not called
    spec["plasma_lo"] = (300, {"JCHAIN": (.9, 9), "MZB1": (.05, 1)})
    r = _row(_t_section(spec, extra_genes=["JCHAIN", "MZB1"]), "plasma_lo")
    assert r["label"] == "Ambiguous" and r["note"] == "no lineage fired"


def test_cluster_epithelial_at_low_scaled_value_is_called():
    spec = {"epi_lo": (400, {"EPCAM": (.9, 1.6), "KRT8": (.9, 1.6), "KRT18": (.05, 1)})}
    a = _t_section(spec)
    r = _row(a, "epi_lo")
    assert 0.4 <= r["best_scaled"] < 0.5 and r["label"] == "Epithelial", r


def test_cluster_tie_is_broken_by_firing_share_or_kept():
    ex = ["MS4A1", "CD79A", "CCR7", "CD1C"]
    # B cells that also carry CCR7: the DC scaled value ties with B, but DC (needs 2 markers) rarely fires per cell
    a = _t_section({"b": (400, {"MS4A1": (.9, 5), "CD79A": (.9, 5), "CCR7": (.9, 5)})}, extra_genes=ex)
    r = _row(a, "b")
    assert r["top3"][0][1] >= 0.99 and r["top3"][1][0] == "DC" and r["top3"][1][1] >= 0.75 * r["top3"][0][1]   # a tie by scaled value
    assert r["label"] == "B" and r["note"] == "tiebreak: fired"
    # cells that really fire both lineages in the same cells (firing shares within 1.5x): the tie stays
    a = _t_section({"bdc": (400, {"MS4A1": (.9, 5), "CD79A": (.9, 5), "CCR7": (.9, 5), "CD1C": (.9, 5)})}, extra_genes=ex)
    r = _row(a, "bdc")
    assert r["label"] == "Ambiguous" and r["note"] == "tie", r


def test_cluster_acta2_only_panel_calls_smc():
    a = _t_section({"smc": (300, {"ACTA2": (.9, 8)}), "t": (300, {**CD3, "CD4": (.7, 3)})}, extra_genes=["ACTA2"])
    P = json.loads(a.uns["annotation_params"])
    assert P["single_marker"].get("SMC") == "ACTA2"
    assert _lab(a)["smc"] == "SMC"


def test_t_subset_treg_needs_two_foxp3_and_no_cd8():
    from bio.annotate import _t_subset
    f = lambda **k: _t_subset({g: np.array([v], np.float32) for g, v in k.items()}, 1)[0]
    assert f(CD4=1, FOXP3=1) == "CD4T"                       # one FOXP3 count is background
    assert f(CD4=1, CD8A=2, FOXP3=2) == "CD8T"               # FOXP3 with CD8A is not a Treg
    assert f(CD4=1, CD8B=2, FOXP3=3) == "CD8T"
    assert f(CD4=1, FOXP3=2) == "Treg" and f(FOXP3=2) == "Treg"
    assert f(IL2RA=1, CTLA4=1, CD4=1) == "CD4T"             # no-FOXP3 panels: same >= 2 and CD8 exclusion
    assert f(IL2RA=2, CTLA4=2, CD4=1, CD8A=1) == "T"   # CD8 = CD4 tie, not Treg
    assert f(IL2RA=2, CTLA4=2, CD4=1) == "Treg"


def test_cell_method_t_subsets_use_the_cluster_rule():
    """Cell method: no subset evidence -> plain T (not CD4T); one FOXP3 count is not a Treg; FOXP3 >= 2 with CD8 is not a Treg."""
    rng = np.random.default_rng(9)
    genes = ["EPCAM", "KRT8", "KRT18", "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B", "FOXP3", "IL2RA", "CTLA4", "ACTB", "GAPDH"]
    g = {x: i for i, x in enumerate(genes)}
    plan = {"none": 200, "foxp3_1": 200, "treg": 200, "foxp3_cd8": 200, "cd8": 200, "cd4": 200}
    rows, kinds = [], []
    n_epi = 1500
    M = rng.poisson(0.004, (n_epi, len(genes))).astype(np.float32)
    M[:, [g["ACTB"], g["GAPDH"]]] += rng.poisson(12, (n_epi, 2))
    M[:, [g["EPCAM"], g["KRT8"], g["KRT18"]]] += rng.poisson(4, (n_epi, 3))
    rows.append(M)
    kinds += ["epi"] * n_epi
    for k, n in plan.items():
        M = np.zeros((n, len(genes)), np.float32)
        M[:, [g["ACTB"], g["GAPDH"]]] += rng.poisson(12, (n, 2))
        M[:, [g["CD3D"], g["CD3E"], g["CD3G"], g["CD2"]]] += 3
        if k == "foxp3_1":
            M[:, g["FOXP3"]] = 1
        if k == "treg":
            M[:, g["FOXP3"]] = 3
        if k == "foxp3_cd8":
            M[:, g["FOXP3"]] = 3
            M[:, g["CD8A"]] = 3
        if k == "cd8":
            M[:, g["CD8A"]] = 3
        if k == "cd4":
            M[:, g["CD4"]] = 3
        rows.append(M)
        kinds += [k] * n
    a = ad.AnnData(sparse.csr_matrix(np.vstack(rows)), var=pd.DataFrame(index=genes))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a = annotate(a, method="cell")
    lab, kinds = a.obs["annotation"].astype(str).to_numpy(), np.array(kinds)
    want = {"none": "T", "foxp3_1": "T", "treg": "Treg", "foxp3_cd8": "CD8T", "cd8": "CD8T", "cd4": "CD4T"}
    for k, w in want.items():
        assert (lab[kinds == k] == w).mean() > 0.95, (k, w)
    P = json.loads(a.uns["annotation_params"])
    assert P["method"] == "cell" and P["t_default_cd4"] == (lab == "T").sum()
    from bio.annotate import report_text
    txt = report_text(annotate_report(a))
    assert "labelled CD4T by default" not in txt and "have no subset evidence" in txt
