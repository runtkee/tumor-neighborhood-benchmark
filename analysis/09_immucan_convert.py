"""Convert IMMUcan panel-1 sce_raw.rds (SingleCellExperiment) to h5ad without R."""
import sys
import numpy as np
import pandas as pd
import anndata as ad
import rdata
from bio import rds

src, dst = sys.argv[1], sys.argv[2]
o = rdata.parser.parse_file(src).object
s = rds.slots(o)
col = rds.dataframe(s["colData"])
assays = rds.named_list(rds.slots(rds.slots(s["assays"])["data"])["listData"])
print("assays:", list(assays), "| colData:", list(col.columns))
row = rds.dataframe(s["elementMetadata"])
names = rds.slots(s["NAMES"]) if s["NAMES"].info.type != rds.T.SYM else None
mats = {k: rds.matrix(v) for k, v in assays.items()}
first = next(iter(mats.values()))
print("matrix shape (markers x cells):", first.shape)
markers = None
dn = rds.attrs(assays[list(assays)[0]]).get("dimnames")
if dn is not None:
    markers = rds.strings(rds.deref(dn).value[0])
for k in ["name", "channel", "clean_target", "Target", "marker"]:
    if markers is None and k in row.columns:
        markers = list(row[k].astype(str)); break
if markers is None and "dimnames" in rds.attrs(assays[list(assays)[0]]):
    dn = rds.attrs(assays[list(assays)[0]])["dimnames"].value
    markers = rds.strings(dn[0])
print("rowData columns:", list(row.columns))
a = ad.AnnData(X=mats.get("counts", first).T.astype(np.float32), obs=col.reset_index(drop=True),
               var=pd.DataFrame(index=markers if markers else [f"m{i}" for i in range(first.shape[0])]))
for k, m in mats.items():
    a.layers[k] = m.T.astype(np.float32)
a.var = a.var.join(row.reset_index(drop=True).set_index(a.var.index), how="left") if len(row) == a.n_vars else a.var
a.obs.columns = a.obs.columns.astype(str)
for c in a.obs.columns:
    if a.obs[c].dtype == object:
        a.obs[c] = a.obs[c].astype(str)
a.write_h5ad(dst)
print(a)
