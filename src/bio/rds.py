"""Minimal reader for R SingleCellExperiment .rds files without R, via the `rdata` low-level parser.

Navigates S4 slots and pulls out: assays (dense numeric matrices), colData (DataFrame), rownames.
Only handles the plain (non-DelayedArray) case; fails loudly otherwise.
"""
import numpy as np
import pandas as pd
import rdata

T = rdata.parser.RObjectType


def deref(o):
    """Follow R serialization back-references (REFSXP)."""
    while o is not None and o.info.type == T.REF:
        o = o.referenced_object
    return o


def _sym(tag):
    tag = deref(tag)
    return tag.value.value.decode() if tag is not None else None


def pairlist(pl):
    """R pairlist -> list of (name, RObject)."""
    out = []
    pl = deref(pl)
    while pl is not None and pl.info.type == T.LIST:
        out.append((_sym(pl.tag), deref(pl.value[0])))
        pl = deref(pl.value[1]) if len(pl.value) > 1 else None
    return out


def attrs(o):
    o = deref(o)
    return dict(pairlist(o.attributes)) if o.attributes is not None else {}


def strings(o):
    """STRSXP -> list of str."""
    return [None if deref(c).value is None else deref(c).value.decode() for c in deref(o).value]


def vector(o):
    """Atomic vector -> numpy / list; factors -> categorical strings."""
    a = attrs(o)
    if o.info.type == T.STR:
        return np.array(strings(o), dtype=object)
    vals = np.asarray(o.value)
    if "levels" in a:  # factor
        lev = np.array(strings(a["levels"]), dtype=object)
        return pd.Categorical.from_codes(np.where(vals > 0, vals - 1, -1), lev)
    return vals


def named_list(o):
    """VECSXP with names -> dict name -> RObject."""
    names = strings(attrs(o)["names"]) if "names" in attrs(o) else [str(i) for i in range(len(o.value))]
    return dict(zip(names, [deref(v) for v in o.value]))


def slots(o):
    return attrs(o)


def dataframe(o):
    """S4 DataFrame (S4Vectors) -> pandas DataFrame."""
    s = slots(o)
    cols = named_list(s["listData"])
    df = pd.DataFrame({k: vector(v) for k, v in cols.items()})
    if "rownames" in s and s["rownames"].info.type == T.STR:
        df.index = strings(s["rownames"])
    return df


def matrix(o):
    """Dense numeric matrix (with dim attribute) -> numpy array (rows x cols)."""
    a = attrs(o)
    if "dim" not in a:
        raise ValueError(f"not a dense matrix (type {o.info.type}); DelayedArray/sparse not supported")
    nr, nc = np.asarray(a["dim"].value)
    return np.asarray(o.value).reshape((nc, nr)).T  # R is column-major
