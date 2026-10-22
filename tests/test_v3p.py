"""Recipe v3p: tumor compartments, composition, T subsets, nest-restricted proximity, gradients, macrophage state.
(1) the v2p.1 columns must be bit-identical to the frozen snapshot; (2) synthetic data with a planted tumor nest."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from bio.compartment import ORDER, assign, compartments_for
from bio.kit import checks, manifest
from bio.kit.compute import compute
from bio.kit.importers import load_cells
from bio.prng import Mulberry32, choose, fnv1a
from bio.recipe_protein import GRAD_LABELS, N_SHUFFLE, VERSION, ZONE_UM, _nearest, grad_bin
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "codex_subset_3pt.csv"
SNAP = ROOT / "tests" / "fixtures" / "v2p1_snapshot_codex3pt.json"


# ---------------------------------------------------------------- v2p.1 columns unchanged ---------------------------
def test_v2p1_columns_bit_identical_and_recipe_is_v3p1():
    m = manifest.load(ROOT / "kit_manifests" / "CODEX_CRC_Schurch2020.yaml")
    m.update(path=str(FIX), dataset="T")
    assert manifest.problems(m) == []
    rows = compute(load_cells(m), m)
    snap = json.loads(SNAP.read_text())
    assert VERSION == "v3p.1" and set(rows["recipe"]) == {"v3p.1"}
    legacy_cols = [c for c in snap[0] if c != "recipe"]
    assert list(rows.columns[:len(snap[0])]) == list(snap[0])  # same columns, same order, legacy first
    got = rows.set_index("patient")
    assert sorted(got.index) == sorted(r["patient"] for r in snap)
    for rec in snap:
        for c in legacy_cols:
            if c == "patient":
                continue
            want, have = rec[c], got.loc[rec["patient"], c]
            if want is None:
                assert pd.isna(have), (rec["patient"], c)
            elif c == "dataset":
                assert have == want
            else:
                assert abs(have - want) <= 1e-12, (rec["patient"], c, have, want)


# ---------------------------------------------------------------- synthetic nest --------------------------------------
def _write(tmp_path, px=0.5, n_pat=6, tumor_free=False, annotate=False):
    """Each image is 1000 x 600 um. Left half (x < 500) is a tumor nest (72% tumor cells), right half is immune-rich stroma.
    The nest has three strips: A (x < 170) macrophages but no T cells, B (170-330) T cells but no macrophages,
    C (330-500) both. That leaves nest T cells far from macrophages and nest macrophages far from T cells, so both
    contact tests have enough cells. Planted: 30% of nest T cells sit right next to a macrophage (shows in the 0-5 um gradient bin); nest T cells < 10 um
    from a macrophage have higher Ki67; macrophages < 10 um from a T cell have higher PD-L1."""
    rng = np.random.default_rng(11)
    labs_l = ["CD8", "CD4", "Treg", "Tgen", "Neu", "Mac", "Bcell", "Tum"]
    p_l = np.array([.06, .04, .02, .02, .03, .08, .03, .72])
    labs_r = labs_l + ["Str"]
    p_r = np.array([.17, .12, .03, .03, .12, .12, .21, .10, .10])
    t_labs = ["CD8", "CD4", "Treg", "Tgen"]
    frames = []
    for p in range(n_pat + int(tumor_free)):
        free = p == n_pat
        left = np.column_stack([rng.uniform(0, 500, 2300), rng.uniform(0, 600, 2300)])
        right = np.column_stack([rng.uniform(500, 1000, 800), rng.uniform(0, 600, 800)])
        xy = np.vstack([left, right])
        lab = np.concatenate([rng.choice(labs_l, 2300, p=p_l), rng.choice(labs_r, 800, p=p_r)])
        lab = np.where(np.isin(lab, t_labs) & (xy[:, 0] < 170), "Bcell", lab)  # strip A: no T cells
        lab = np.where((lab == "Mac") & (xy[:, 0] >= 170) & (xy[:, 0] < 330), "Neu", lab)  # strip B: no macrophages
        if free:
            lab = np.where(lab == "Tum", "Str", lab)
        is_t = np.isin(lab, t_labs)
        mac = np.where(lab == "Mac")[0]
        t_nest = np.where(is_t & (xy[:, 0] < 500))[0]
        mac_c = mac[(xy[mac, 0] >= 330) & (xy[mac, 0] < 500)]
        move = rng.choice(t_nest, int(.3 * len(t_nest)), replace=False)
        xy[move] = xy[rng.choice(mac_c, len(move))] + rng.normal(0, 3, (len(move), 2))
        dm = np.min(np.linalg.norm(xy[:, None, :] - xy[mac][None, :, :], axis=2), axis=1)
        dt = np.min(np.linalg.norm(xy[:, None, :] - xy[is_t][None, :, :], axis=2), axis=1)
        ki = rng.normal(0, 1, len(lab))
        ki[is_t & (dm < 10) & (xy[:, 0] < 500)] += 2
        pdl1 = rng.normal(0, 1, len(lab))
        pdl1[(lab == "Mac") & (dt < 10)] += 2
        cd68 = rng.normal(0, 1, len(lab)); cd68[lab == "Mac"] += 3
        mpo = rng.normal(0, 1, len(lab)); mpo[lab == "Neu"] += 3
        frames.append(pd.DataFrame({
            "pt": f"P{p}", "img": f"P{p}_1", "cx": xy[:, 0] / px, "cy": xy[:, 1] / px, "lab": lab,
            "Ki67": ki, "PDL1": pdl1, "CD68": cd68, "MPO": mpo, "region": np.where(xy[:, 0] < 500, "Tumor", "Stroma")}))
    pd.concat(frames).to_csv(tmp_path / "cells.csv", index=False)
    m = {"dataset": "SYN", "modality": "protein", "path": str(tmp_path / "cells.csv"),
         "columns": {"patient": "pt", "field": "img", "x": "cx", "y": "cy", "type": "lab"},
         "pixel_size_um": px,
         "roles": {"CD8": "T:CD8", "CD4": "T:CD4", "Treg": "T:Treg", "Tgen": "T", "Neu": "neut", "Mac": "mac", "Bcell": "B",
                   "Tum": "tumor", "Str": "stroma"},
         "states": {"prolif": "Ki67"}, "bleed": {"mac": "CD68", "neut": "MPO"}, "mac_states": {"PDL1": "PDL1"}}
    if annotate:
        m["columns"]["compartment"] = "region"
        m["compartments"] = {"Tumor": "nest", "Stroma": "far"}
    (tmp_path / "m.yaml").write_text(yaml.safe_dump(m))
    return manifest.load(tmp_path / "m.yaml")


@pytest.fixture(scope="module")
def syn(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("syn")
    m = _write(tmp, tumor_free=True)
    df = load_cells(m)
    return m, df, compute(df, m)


def test_manifest_and_checks_accept_subsets(syn):
    m, df, _ = syn
    assert manifest.problems(m) == []
    res = {n: (s, d) for n, s, d in checks.run_checks(df, m)}
    assert not checks.blocking(list(map(lambda kv: (kv[0], kv[1][0], kv[1][1]), res.items())))
    assert res["tumor cells"][0] == "pass" and "1 of 7 images" in res["tumor cells"][1]
    assert res["nest size"][0] == "pass"
    assert res["T subsets"][0] == "pass" and all(s in res["T subsets"][1] for s in ("CD8", "CD4", "Treg"))
    assert set(df.loc[df["type"] == "T", "sub"]) == {"CD8", "CD4", "Treg", "T"}


def test_compartment_counts_add_up(syn):
    _, df, rows = syn
    got = rows.set_index("patient")
    n_cells = df.groupby("patient").size()
    for pat in n_cells.index:
        total = sum(got.loc[pat, f"comp_{c}_cells"] for c in ("nest", "margin", "far"))
        if pat == "P6":  # image with no tumor cells: compartments unknown, no compartment metrics
            assert total == 0
            assert pd.isna(got.loc[pat, "nest_T_near_mac_log2"]) and pd.isna(got.loc[pat, "comp_nest_pct_T"])
        else:
            assert total == n_cells[pat]
            assert got.loc[pat, "comp_nest_cells"] > 0.8 * 2300  # most of the left half is nest
            assert got.loc[pat, "comp_margin_cells"] > 0 and got.loc[pat, "comp_far_cells"] > 0


def test_composition_percentages(syn):
    _, df, rows = syn
    got = rows[rows["patient"] != "P6"].set_index("patient")
    for c in ("nest", "margin", "far"):
        t = got[f"comp_{c}_pct_T"]
        subs = sum(got[f"comp_{c}_pct_{s}"] for s in ("CD8", "CD4", "Treg"))
        imm = sum(got[f"comp_{c}_pct_{r}"] for r in ("T", "neut", "mac", "B"))
        share = sum(got[f"comp_{c}_immshare_{r}"] for r in ("T", "neut", "mac", "B"))
        assert (subs <= t + 1e-9).all() and (imm <= 100).all()
        assert np.allclose(share, 100)  # immune shares add to 100 by construction
    assert (got["comp_nest_pct_T"] < got["comp_far_pct_T"]).all()  # the nest is tumor-dominated, the far zone is immune-rich
    assert (got["comp_nest_pct_T"] > 0).all()


def test_subset_labels_land_in_subset_columns(syn):
    _, df, rows = syn
    d = compartments_for(df.copy())
    got = rows.set_index("patient")
    for pat, g in d.groupby("patient"):
        if pat == "P6":
            continue
        nest = g[g["comp"] == "nest"]
        for s in ("CD8", "CD4", "Treg"):
            want = 100 * ((nest["type"] == "T") & (nest["sub"] == s)).sum() / len(nest)
            assert want > 0 and abs(got.loc[pat, f"comp_nest_pct_{s}"] - want) < 1e-9
    for c in ("nest_CD8_near_mac_log2", "nest_CD4_near_neut_log2", "nest_Treg_near_mac_log2"):
        assert c in rows.columns


def test_planted_effect_in_nest(syn):
    _, _, rows = syn
    r = rows[rows["patient"] != "P6"]
    assert (r["nest_macTouch_prolif_diff"] > 10).all()          # planted: T cells touching macrophages in the nest
    assert r["nest_T_near_mac_log2"].notna().all()              # (its sign is not tested: the nest has strips with and without macrophages)
    assert (r["TmacTouch_PDL1_diff"] > 10).all()                # planted: macrophages touching T cells carry more PD-L1
    assert (r["TmacTouch_PDL1_n"] >= 20).all()
    assert (r["nest_TmacTouch_PDL1_diff"] > 10).all() and (r["nest_TmacTouch_PDL1_n"] >= 20).all()


def test_gradient_bins(syn):
    _, _, rows = syn
    r = rows[rows["patient"] != "P6"]
    assert (r["grad_mac_0-5_log2"] > 0).all()  # planted clustering shows in the closest bin
    assert (r["grad_mac_0-5_log2"] > r["grad_mac_20-30_log2"] + 1).all()
    for src in ("neut", "mac"):
        for lab in GRAD_LABELS:
            assert f"grad_{src}_{lab}_log2" in rows.columns and f"grad_{src}_{lab}_prolif" in rows.columns
    # planted Ki67 effect: T cells within 10 um of a macrophage have a higher mean Ki67 percentile than those 20-30 um away
    assert (r["grad_mac_0-5_prolif"] > r["grad_mac_20-30_prolif"] + 15).all()


def test_grad_bin_edges_and_counts():
    d = np.array([0, 4.99, 5, 9.99, 10, 19.9, 20, 29.9, 30, 49.9, 50, 99.9, 100, 500.0])
    assert grad_bin(d).tolist() == [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6]
    rng = np.random.default_rng(3)
    d = rng.exponential(40, 1000)
    counts = np.bincount(grad_bin(d), minlength=len(GRAD_LABELS))
    assert counts.sum() == len(d) and counts[:4].sum() == (d < 30).sum()  # bins 0-5..20-30 are exactly "within 30 um"


def test_assign_geometry():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.uniform(0, 200, 800), rng.uniform(400, 600, 300)])
    y = rng.uniform(0, 200, len(x))
    tumor = np.concatenate([np.ones(800, bool), np.zeros(300, bool)])
    comp = assign(x, y, tumor)
    assert set(comp) == {"nest", "far"}
    assert (comp[:800] == "nest").mean() > 0.9
    assert (comp[800:] == "far").all()  # 200+ um from the nest
    # a stripe just outside the nest is margin
    x2, y2, t2 = np.append(x, [260.0, 295.0]), np.append(y, [100.0, 100.0]), np.append(tumor, [False, False])
    assert assign(x2, y2, t2)[-2:].tolist() == ["margin", "margin"]


def test_explicit_compartment_column_wins(tmp_path):
    m = _write(tmp_path, n_pat=2, annotate=True)
    assert manifest.problems(m) == []
    df = load_cells(m)
    assert set(df["comp"]) == {"nest", "far"}
    rows = compute(df, m)
    for pat, g in df.groupby("patient"):
        r = rows.set_index("patient").loc[pat]
        assert r["comp_nest_cells"] == (g["comp"] == "nest").sum() == (g["x_um"] < 500).sum()
        assert r["comp_margin_cells"] == 0 and r["comp_far_cells"] == (g["x_um"] >= 500).sum()


def test_no_tumor_label_still_runs_whole_image(tmp_path):
    m = _write(tmp_path, n_pat=2)
    m["roles"]["Tum"] = "other"
    df = load_cells(m)
    res = {n: s for n, s, _ in checks.run_checks(df, m)}
    assert res["tumor cells"] == "warn"
    rows = compute(df, m)
    assert not [c for c in rows.columns if c.startswith(("comp_", "nest_", "grad_"))]
    assert "macTouch_prolif_diff" in rows.columns and "TmacTouch_PDL1_diff" in rows.columns


def test_deterministic(syn):
    m, df, rows = syn
    pd.testing.assert_frame_equal(rows, compute(df, m))


# ---------------------------------------------------------------- v3p.1: tumor-free null pool -------------------------
def _write_segregated(tmp_path, n_pat=6, n_imm=1200, n_tum=3600, p_mac=0.03):
    """Planted case for the nest null pool. One image per patient, 600 x 600 um (pixel size 1). Immune and stromal cells sit only
    in six vertical strands (x mod 100 < 40) and tumor cells only in the gaps between them (x mod 100 >= 40). Tumor density is twice
    the immune density, so about two thirds of the cells within 50 um of a strand cell are tumor cells and the strands count as
    tumor nest. Inside the strands the macrophage labels are random (3% of cells), so macrophages are NOT enriched around T cells.
    Old pool (non-T cells of the compartment, tumor included): most shuffled macrophages land in the tumor gaps, far from the
    strand T cells, so the null is too low and the real share looks like enrichment (log2 about +1). New pool (tumor cells
    excluded): the shuffled macrophages land in the strands like the real ones, so log2 is about 0."""
    rng = np.random.default_rng(5)
    frames = []
    for p in range(n_pat):
        xi = 100 * rng.integers(0, 6, n_imm) + rng.uniform(0, 40, n_imm)
        xt = 100 * rng.integers(0, 6, n_tum) + rng.uniform(40, 100, n_tum)
        lab = rng.choice(["T", "mac", "neut", "B", "stroma"], n_imm, p=[.3, p_mac, .1, .3 - p_mac, .3])
        frames.append(pd.DataFrame({
            "pt": f"P{p}", "img": f"P{p}_1", "cx": np.r_[xi, xt], "cy": np.r_[rng.uniform(0, 600, n_imm), rng.uniform(0, 600, n_tum)],
            "lab": np.r_[lab, np.full(n_tum, "tumor")], "Ki67": rng.normal(0, 1, n_imm + n_tum)}))
    pd.concat(frames).to_csv(tmp_path / "cells.csv", index=False)
    m = {"dataset": "SEG", "modality": "protein", "path": str(tmp_path / "cells.csv"),
         "columns": {"patient": "pt", "field": "img", "x": "cx", "y": "cy", "type": "lab"}, "pixel_size_um": 1,
         "roles": {r: r for r in ("T", "mac", "neut", "B", "stroma", "tumor")}, "states": {"prolif": "Ki67"}}
    (tmp_path / "m.yaml").write_text(yaml.safe_dump(m))
    return manifest.load(tmp_path / "m.yaml")


def _nest_t_mac_log2(df, exclude_tumor):
    """Independent re-implementation of nest_T_near_mac_log2 for the one-image-per-patient planted case; exclude_tumor=False is the
    pre-v3p.1 pool (all non-T cells of the compartment), True the v3p.1 pool (non-T, non-tumor). Same seed, same draws."""
    df = compartments_for(df.copy())
    out = {}
    for (pat, field), g in df.groupby(["patient", "field"]):
        xy, cm, ty = g[["x_um", "y_um"]].to_numpy(), g["comp"].to_numpy(), g["type"].to_numpy()
        T, M, tum = ty == "T", ty == "mac", ty == "tumor"
        Tn = T & (cm == "nest")
        real = (_nearest(xy, Tn, M) < ZONE_UM).sum()
        rng = Mulberry32(fnv1a(f"{pat}|{field}|T|mac|nest"))
        pools = {c: np.where(~T & ((~tum) if exclude_tumor else True) & (cm == c))[0] for c in ORDER}
        k = {c: int((M & (cm == c)).sum()) for c in ORDER}
        null = []
        for _ in range(N_SHUFFLE):
            mask = np.zeros(len(g), bool)
            for c in ORDER:
                kk = min(k[c], len(pools[c]))
                if kk > 0:
                    mask[choose(rng, pools[c], kk)] = True
            null.append((_nearest(xy, Tn, mask) < ZONE_UM).sum())
        out[pat] = float(np.log2((real + 1) / (np.mean(null) + 1)))
    return out


def test_nest_null_pool_excludes_tumor_cells(tmp_path):
    m = _write_segregated(tmp_path)
    df = load_cells(m)
    rows = compute(df, m).set_index("patient")
    assert (rows["comp_nest_pct_T"] > 5).all() and (rows["comp_nest_cells"] > 0.9 * df.groupby("patient").size() - 500).all()  # strands are nest
    old, new = _nest_t_mac_log2(df, False), _nest_t_mac_log2(df, True)
    for pat in rows.index:
        assert abs(rows.loc[pat, "nest_T_near_mac_log2"] - new[pat]) < 1e-12          # the recipe uses the new pool
        assert old[pat] > 0.8, (pat, old[pat])                                          # old pool: spurious enrichment
        assert abs(new[pat]) < 0.3, (pat, new[pat])                                     # new pool: none
    assert abs(np.median(list(new.values()))) < 0.1 and np.median(list(old.values())) > 0.9


# ---------------------------------------------------------------- v3p.1: pooled sums ----------------------------------
def test_pooled_sums_consistent_with_diff(syn):
    _, _, rows = syn
    bases = sorted(c[:-9] for c in rows.columns if c.endswith("_near_sum"))
    assert bases and all(b.startswith(("nest_", "TmacTouch_")) for b in bases)
    assert not [c for c in rows.columns if c.endswith("_near_sum") and c.startswith(("neutTouch", "macTouch", "stromaTouch"))]  # legacy untouched
    checked = 0
    for b in bases:
        cols = [f"{b}_{k}" for k in ("n", "near_sum", "far_n", "far_sum", "diff")]
        i = list(rows.columns).index(cols[0])
        assert list(rows.columns[i:i + 5]) == cols                                      # column order
        for _, r in rows.iterrows():
            if r[cols[0]] >= 20 and r[cols[2]] >= 20:
                assert abs(r[cols[4]] - (r[cols[1]] / r[cols[0]] - r[cols[3]] / r[cols[2]])) < 1e-9
                checked += 1
            else:
                assert pd.isna(r[cols[4]])                                              # 20-cell rule for diff
    assert checked > 20
    # independent recomputation of one test for the first patient: nest T cells touching a macrophage, Ki67 percentile
    m, df, _ = syn
    d = compartments_for(df.copy())
    g = d[(d["patient"] == "P0")]
    pct = g["prolif"].rank(pct=True).to_numpy() * 100
    xy, ty, cm = g[["x_um", "y_um"]].to_numpy(), g["type"].to_numpy(), g["comp"].to_numpy()
    Tn = (ty == "T") & (cm == "nest")
    dist = cKDTree(xy[ty == "mac"]).query(xy[Tn], k=1)[0]
    v = pct[Tn]
    r0 = rows.set_index("patient").loc["P0"]
    assert r0["nest_macTouch_prolif_n"] == (dist < 10).sum() and r0["nest_macTouch_prolif_far_n"] == (dist >= 30).sum()
    assert abs(r0["nest_macTouch_prolif_near_sum"] - math.fsum(v[dist < 10])) < 1e-9
    assert abs(r0["nest_macTouch_prolif_far_sum"] - math.fsum(v[dist >= 30])) < 1e-9


def test_macrophage_state_whole_image_and_nest(syn):
    _, df, rows = syn
    r = rows[rows["patient"] != "P6"]
    assert (r["TmacTouch_PDL1_n"] >= r["nest_TmacTouch_PDL1_n"]).all()                # nest macrophages are a subset of all macrophages
    assert (r["TmacTouch_PDL1_far_n"] >= r["nest_TmacTouch_PDL1_far_n"]).all()
    assert (r["TmacTouch_PDL1_n"] > r["nest_TmacTouch_PDL1_n"]).any()                 # (the stroma half also has macrophages)
    # whole-image count checked independently on P1: macrophages < 10 um from any T cell
    g = df[df["patient"] == "P1"]
    xy, ty = g[["x_um", "y_um"]].to_numpy(), g["type"].to_numpy()
    d = cKDTree(xy[ty == "T"]).query(xy[ty == "mac"], k=1)[0]
    got = rows.set_index("patient").loc["P1"]
    assert got["TmacTouch_PDL1_n"] == (d < 10).sum() and got["TmacTouch_PDL1_far_n"] == (d >= 30).sum()
    p6 = rows.set_index("patient").loc["P6"]                                          # image without tumor cells: whole image only
    assert p6["TmacTouch_PDL1_n"] > 0 and pd.isna(p6["nest_TmacTouch_PDL1_n"])


def test_tmactouch_without_tumor_label_has_no_nest_version(tmp_path):
    m = _write(tmp_path, n_pat=2)
    m["roles"]["Tum"] = "other"
    rows = compute(load_cells(m), m)
    for k in ("n", "near_sum", "far_n", "far_sum", "diff"):
        assert f"TmacTouch_PDL1_{k}" in rows.columns
    assert (rows["TmacTouch_PDL1_n"] > 0).all()
    assert not [c for c in rows.columns if c.startswith("nest_")]
