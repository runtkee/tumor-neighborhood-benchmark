"""Benchmark kit command line. Runs fully offline.

  python -m bio.kit suggest <file> --dataset NAME > manifest.yaml   draft a manifest to review (# CHECK marks guesses)
  python -m bio.kit check <manifest.yaml>                            validate manifest + data, print checks
  python -m bio.kit run <manifest.yaml> --out DIR [--force] [--exclude-ref NAME ...]
                                                                      compute, compare to the reference, write report.html
  python -m bio.kit xenium <bundle_dir_or_h5ad> [...] --dataset NAME --out DIR [--patient ID | --patient-from parent|dirname]
                                      [--areas <coordinates.csv|none> ...] | --samplesheet <samplesheet.csv> [--regions R ...]
                                                                      read a 10x Xenium bundle, label cells from markers, write h5ad + report + draft manifest
  python -m bio.kit build-reference <manifest.yaml> ...              (maintainers) add public datasets to the bundled reference
"""
import argparse
import secrets
import sys
from pathlib import Path

import pandas as pd

from . import checks as ck
from . import manifest as mf
from .compute import compute
from .importers import load_cells
from .report import build
from .suggest import draft

REF_DIR = Path(__file__).parent / "reference"


def _load_checked(path, force=False):
    m = mf.load(path)
    probs = mf.problems(m)
    if probs:
        sys.exit("manifest problems:\n  " + "\n  ".join(probs))
    df = load_cells(m)
    res = ck.run_checks(df, m)
    for n, s, d in res:
        print(f"[{s:4s}] {n}: {d}")
    bad = ck.blocking(res)
    if bad and not force:
        sys.exit(f"\n{len(bad)} check(s) failed; fix the manifest (or pass --force, which the report will flag).")
    return m, df, res, bool(bad)


def _reference(exclude):
    files = sorted(REF_DIR.glob("*.csv"))
    ref = pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else pd.DataFrame(columns=["dataset", "patient"])
    ref = ref[~ref["dataset"].isin(exclude)]
    src = (REF_DIR / "SOURCES.html").read_text() if (REF_DIR / "SOURCES.html").exists() else ""
    return ref, src


def main(argv=None):
    p = argparse.ArgumentParser(prog="bio.kit")
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("suggest"); s.add_argument("file"); s.add_argument("--dataset", default="MY_DATASET")
    c = sub.add_parser("check"); c.add_argument("manifest")
    r = sub.add_parser("run"); r.add_argument("manifest"); r.add_argument("--out", required=True)
    r.add_argument("--force", action="store_true"); r.add_argument("--exclude-ref", nargs="*", default=[])
    from .xenium_cmd import add_args as _xen_args
    _xen_args(sub.add_parser("xenium"))
    b = sub.add_parser("build-reference"); b.add_argument("manifests", nargs="+")
    b.add_argument("--force", action="store_true", help="build even if checks fail (maintainers; say why in the manifest)")
    a = p.parse_args(argv)

    if a.cmd == "suggest":
        print(draft(a.file, a.dataset), end="")
    elif a.cmd == "check":
        _load_checked(a.manifest, force=True)
    elif a.cmd == "run":
        m, df, res, forced = _load_checked(a.manifest, a.force)
        rows = compute(df, m)
        ref, src = _reference(a.exclude_ref + [m["dataset"]])
        out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
        rows.to_csv(out / "patients_full.csv", index=False)  # stays local: real patient IDs
        # random codes, not hashes: a hash of a short ID ("1".."35") can be reversed by guessing
        codes = {p: "P-" + secrets.token_hex(4) for p in rows["patient"].unique()}
        share = rows.copy()
        share["patient"] = share["patient"].map(codes)
        share.to_csv(out / "share_summary.csv", index=False)  # per-patient effect sizes only, no cell data
        pd.DataFrame({"code": list(codes.values()), "your_patient_id": list(codes.keys())}).to_csv(out / "patient_key_KEEP_PRIVATE.csv", index=False)
        path = build(out, m, res, forced, rows, ref, src)
        print(f"\nreport: {path}\nshareable summary (no cell data, random patient codes): {out / 'share_summary.csv'}"
              f"\nkey linking codes to patients (keep private): {out / 'patient_key_KEEP_PRIVATE.csv'}")
    elif a.cmd == "xenium":
        from .xenium_cmd import run as _xen_run
        _xen_run(a)
    elif a.cmd == "build-reference":
        REF_DIR.mkdir(exist_ok=True)
        for mp in a.manifests:
            m, df, _, forced = _load_checked(mp, a.force)
            if forced:
                print(f"WARNING: {m['dataset']} built with --force despite failed checks")
            rows = compute(df, m)
            rows.to_csv(REF_DIR / f"{m['dataset']}.csv", index=False)
            print(f"reference: {m['dataset']} ({rows['patient'].nunique()} patients)")


if __name__ == "__main__":
    main()
