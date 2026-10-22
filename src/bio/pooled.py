"""Pooled estimate for the sparse contact tests (recipe v3p.1). Ported 1:1 to `pooled()` in web/findings.js.

Per-patient contact differences need 20 cells in each arm, so inside tumor nests only a few patients qualify.
The recipe therefore also reports, per patient and per test, the sum of the within-image percentiles in each arm
(`<key>_near_sum`, `<key>_far_sum`) and the arm sizes (`<key>_n`, `<key>_far_n`). Pooled over patients:

    diff = sum(near_sum) / sum(n)  -  sum(far_sum) / sum(far_n)

Each cell counts once, so a patient with many cells weighs more than one with few. That is a pooled estimate, not a
per-patient one, and the app says so. The 95% interval is a patient bootstrap: 200 resamples of patients with
replacement, the draws seeded by fnv1a("pooled|" + key) with the shared Mulberry32, 2.5 and 97.5 percentiles
(linear interpolation) of the resampled estimates.
"""
import numpy as np

from bio.prng import Mulberry32, fnv1a

N_BOOT = 200


def _num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


def pooled(rows, key):
    """rows: DataFrame or list of dicts, one per patient; key: the test without suffix, e.g. "nest_neutTouch_PD1".
    Uses the patients with n > 0 or far_n > 0. Returns None when no patient has the sum columns, else a dict:
    diff, lo, hi (95% bootstrap interval; NaN when an arm is empty), n_near, n_far (cells), patients (patients with
    near cells), used (patients contributing to either arm)."""
    recs = rows.to_dict("records") if hasattr(rows, "to_dict") else list(rows)
    P = []
    for r in recs:
        v = [_num(r.get(f"{key}{s}")) for s in ("_n", "_near_sum", "_far_n", "_far_sum")]
        if None in v or not (v[0] > 0 or v[2] > 0):
            continue
        P.append(v)
    if not P:
        return None

    def est(sn, ss, fn, fs):
        return ss / sn - fs / fn if sn > 0 and fn > 0 else float("nan")

    sn = ss = fn = fs = 0.0
    for a, b, c, d in P:
        sn += a; ss += b; fn += c; fs += d
    diff = est(sn, ss, fn, fs)
    rng = Mulberry32(fnv1a("pooled|" + key))
    m, vals = len(P), []
    for _ in range(N_BOOT):
        bn = bs = bf = bg = 0.0
        for _ in range(m):
            a, b, c, d = P[int(rng.random() * m)]
            bn += a; bs += b; bf += c; bg += d
        e = est(bn, bs, bf, bg)
        if np.isfinite(e):
            vals.append(e)
    lo, hi = (float(x) for x in np.percentile(vals, [2.5, 97.5])) if len(vals) >= 2 else (float("nan"), float("nan"))
    return dict(diff=diff, lo=lo, hi=hi, n_near=int(sn), n_far=int(fn), patients=sum(1 for p in P if p[0] > 0), used=m)
