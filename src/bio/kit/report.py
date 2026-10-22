"""Self-contained HTML report: private patients plotted against the public reference, claim by claim.
No JavaScript, no external files, no network: images are embedded PNGs."""
import base64
import html
import io
from datetime import date
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bio.pooled import pooled

from . import MAC_STATES, STATES

STATE_LABEL = {"prolif": "dividing (Ki-67)", "PD1": "PD-1", "GZMB": "granzyme B", "LAG3": "LAG-3", "ICOS": "ICOS",
               "apoptosis": "dying (cleaved PARP)"}


def claims(extra_roles, extra_states):
    st = STATES + list(extra_states)
    c = [
        ("T cells avoid neutrophils", "Share of T cells within 30 µm of a neutrophil vs shuffled labels (log2; below 0 = avoidance).",
         ["T_near_neut_log2"], ["B_near_neut_log2"],
         "Control: B cells. If they avoid neutrophils as much as T cells do, the effect is about lymphocyte zones, not T cells."),
        ("T cells gather near macrophages", "Share of T cells within 30 µm of a macrophage vs shuffled labels (log2).",
         ["T_near_mac_log2"], [], ""),
        ("State of T cells touching neutrophils", "Percentile-point difference, T cells <10 µm vs ≥30 µm from a neutrophil.",
         [f"neutTouch_{s}_diff" for s in st], ["neutTouch_bleed_neut_diff"],
         "Bleed: the neutrophil's own marker measured inside touching T cells. Large bleed means contact comparisons are unreliable."),
        ("State of T cells touching macrophages", "Percentile-point difference, T cells <10 µm vs ≥30 µm from a macrophage.",
         [f"macTouch_{s}_diff" for s in st], ["macTouch_bleed_mac_diff"], "Bleed as above, for the macrophage marker."),
        ("Control: T cells touching stroma", "Same contact test with stromal cells. Effects here flag technical contact artifacts.",
         [f"stromaTouch_{s}_diff" for s in st], ["stromaTouch_bleed_stroma_diff"], ""),
    ]
    c += [
        ("Tumor nest: T cells vs neutrophils", "Same share test restricted to T cells inside the tumor nest; the shuffle keeps the number of neutrophils per compartment.",
         ["nest_T_near_neut_log2", "nest_CD8_near_neut_log2", "nest_CD4_near_neut_log2", "nest_Treg_near_neut_log2"], ["nest_B_near_neut_log2"],
         "Control: B cells in the nest."),
        ("Tumor nest: T cells vs macrophages", "Share of nest T cells within 30 µm of a macrophage vs stratified shuffle (log2).",
         ["nest_T_near_mac_log2", "nest_CD8_near_mac_log2", "nest_CD4_near_mac_log2", "nest_Treg_near_mac_log2"], [], ""),
        ("Tumor nest: state of T cells touching neutrophils", "Percentile-point difference, nest T cells <10 µm vs ≥30 µm from a neutrophil.",
         [f"nest_neutTouch_{s}_diff" for s in st], ["nest_neutTouch_bleed_neut_diff"], ""),
        ("Tumor nest: state of T cells touching macrophages", "Percentile-point difference, nest T cells <10 µm vs ≥30 µm from a macrophage.",
         [f"nest_macTouch_{s}_diff" for s in st], ["nest_macTouch_bleed_mac_diff"], ""),
        ("Macrophage state next to T cells", "Percentile-point difference of the macrophage marker, all macrophages <10 µm from a T cell vs ≥30 µm.",
         [f"TmacTouch_{s}_diff" for s in MAC_STATES], [], ""),
        ("Tumor nest: macrophage state next to T cells", "Same comparison for macrophages inside the tumor nest only.",
         [f"nest_TmacTouch_{s}_diff" for s in MAC_STATES], [], ""),
    ]
    for r, proxy in extra_roles.items():
        c.append((f"Private-only: T cells near {r} (public proxy: {proxy})",
                  f"Public datasets have no '{r}' label; the closest public role '{proxy}' is shown for context.",
                  [f"T_near_{r}_log2"] + [f"{r}Touch_{s}_diff" for s in st], [f"T_near_{proxy}_log2"], ""))
    return c


def _composition(priv):
    """Median [IQR] percent of cells per compartment (one value per patient), or '' when no compartments were computed."""
    roles = [r for r in ["T", "CD8", "CD4", "Treg", "neut", "mac", "B"] if f"comp_nest_pct_{r}" in priv.columns]
    if not roles or "comp_nest_cells" not in priv.columns:
        return ""
    f = lambda v: "—" if not len(v) else f"{v.median():.1f} [{v.quantile(.25):.1f}–{v.quantile(.75):.1f}]"
    out = ["<h2>Composition by compartment</h2><p class='mut'>Percent of all cells in the compartment; median [IQR] over your patients. "
           "Nest = tumor-dominated neighbourhood, margin = within 100 µm of it, far = the rest.</p>",
           "<div class='wrap'><table><tr><th>compartment</th><th>cells (median)</th>" + "".join(f"<th>% {r}</th>" for r in roles) + "</tr>"]
    for c in ["nest", "margin", "far"]:
        cells = priv[f"comp_{c}_cells"].dropna()
        out.append(f"<tr><td>{c}</td><td>{cells.median():.0f}</td>" +
                   "".join(f"<td>{f(priv[f'comp_{c}_pct_{r}'].dropna())}</td>" for r in roles) + "</tr>")
    out.append("</table></div>")
    return "".join(out)


def _pooled_table(metrics, ref, priv, priv_name):
    """Pooled estimate (recipe v3p.1) for contact tests that carry the sums: one row per metric and dataset. '' when none do."""
    f = lambda x: "—" if x is None or pd.isna(x) else f"{x:+.1f}"
    body = []
    for m in metrics:
        if not m.endswith("_diff"):
            continue
        base = m[:-5]
        sets = [(d, g) for d, g in ref.groupby("dataset")] if f"{base}_near_sum" in ref else []
        sets.append((f"YOURS {priv_name}", priv))
        for d, g in sets:
            p = pooled(g, base) if f"{base}_near_sum" in g else None
            if p is None:
                continue
            body.append(f"<tr><td>{html.escape(base)}</td><td>{html.escape(str(d))}</td><td>{p['n_near']:,}</td><td>{p['patients']}</td>"
                        f"<td><b>{f(p['diff'])}</b></td><td>{f(p['lo'])} to {f(p['hi'])}</td></tr>")
    if not body:
        return ""
    return ("<p class='mut'>Pooled across patients (every cell counts once; 95% interval from a patient bootstrap). "
            "Pooled estimate, not per patient.</p><div class='wrap'><table><tr><th>test</th><th>dataset</th><th>cells touching</th>"
            "<th>patients</th><th>pooled difference (pct-pts)</th><th>95% interval</th></tr>" + "".join(body) + "</table></div>")


def _png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _strip(metric, ref, priv, priv_name):
    groups = [(d, g[metric].dropna().to_numpy()) for d, g in ref.groupby("dataset")] if metric in ref else []
    groups = [(d, v) for d, v in groups if len(v)]
    pv = priv[metric].dropna().to_numpy() if metric in priv else np.array([])
    if not groups and not len(pv):
        return None
    fig, ax = plt.subplots(figsize=(1.3 * (len(groups) + 1) + 1.5, 2.6))
    rng = np.random.default_rng(0)
    for i, (d, v) in enumerate(groups):
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=9, color="#8a8f98", alpha=0.6, linewidths=0)
        ax.hlines(np.median(v), i - 0.3, i + 0.3, color="#33363b", lw=2)
    j = len(groups)
    if len(pv):
        ax.scatter(j + rng.uniform(-0.18, 0.18, len(pv)), pv, s=16, color="#c2410c", alpha=0.85, linewidths=0)
        ax.hlines(np.median(pv), j - 0.3, j + 0.3, color="#7c2d12", lw=2.5)
    ax.axhline(0, color="#bbb", lw=0.8, zorder=0)
    ax.set_xticks(range(j + 1), [d.replace("_", "\n", 1) for d, _ in groups] + [f"YOURS\n{priv_name}"], fontsize=7)
    ax.set_title(metric, fontsize=8, loc="left")
    ax.tick_params(axis="y", labelsize=7)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    return _png(fig)


def _row(metric, ref, priv):
    rv = ref[metric].dropna() if metric in ref else pd.Series(dtype=float)
    pv = priv[metric].dropna() if metric in priv else pd.Series(dtype=float)
    pm = pv.median() if len(pv) else np.nan
    pct = (rv < pm).mean() * 100 if len(rv) and len(pv) else np.nan
    agree = ("—" if not len(rv) or not len(pv) or rv.median() == 0 else
             ("same direction" if np.sign(pm) == np.sign(rv.median()) else "opposite"))
    f = lambda x: "—" if pd.isna(x) else f"{x:+.2f}"
    return (f"<tr><td>{html.escape(metric)}</td><td>{len(rv)}</td><td>{f(rv.median() if len(rv) else np.nan)}</td>"
            f"<td>{f(rv.quantile(.25) if len(rv) else np.nan)} to {f(rv.quantile(.75) if len(rv) else np.nan)}</td>"
            f"<td>{len(pv)}</td><td><b>{f(pm)}</b></td><td>{'—' if pd.isna(pct) else f'{pct:.0f}th'}</td><td>{agree}</td></tr>")


CSS = """
:root{--bg:#fbfaf8;--fg:#1d1f22;--mut:#5d636b;--line:#e2e0db;--card:#fff;--accent:#c2410c}
@media (prefers-color-scheme:dark){:root{--bg:#16181b;--fg:#e8e6e3;--mut:#9aa0a6;--line:#2c2f33;--card:#1e2125;--accent:#fb923c}}
body{background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,Segoe UI,sans-serif;margin:0;padding:24px 16px}
main{max-width:1000px;margin:0 auto}h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:32px 0 6px}
.mut{color:var(--mut)}.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0}td,th{border-bottom:1px solid var(--line);padding:4px 6px;text-align:left}
th{color:var(--mut);font-weight:600}.pass{color:#15803d}.warn{color:#b45309}.fail{color:#b91c1c;font-weight:700}
img{max-width:100%;height:auto;background:#fff;border-radius:4px}.plots{display:flex;flex-wrap:wrap;gap:8px}
.wrap{overflow-x:auto}code{font-size:12px}
"""


def build(out_dir, m, checks, forced, priv, ref, ref_sources):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    name = m["dataset"]
    parts = [f"<!doctype html><html><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>Benchmark report</title><style>{CSS}</style></head><body><main>",
             f"<h1>{html.escape(name)} vs public reference</h1>",
             f"<p class='mut'>Generated {date.today()} · recipe {html.escape(str(priv['recipe'].iloc[0]) if len(priv) else '?')} · "
             f"{priv['patient'].nunique()} private patients · reference: {', '.join(sorted(ref['dataset'].unique()))} "
             f"({ref['patient'].nunique()} patients). Everything was computed locally; no data left this machine.</p>",
             "<div class='card'><b>How to read this.</b> Grey dots are public patients (one per patient, line = median); orange dots are "
             "your patients. Each claim has its control next to it. A private effect counts as consistent with public data when it "
             "falls inside the public spread <i>and</i> the control doesn't show the same effect.</div>",
             "<h2>Data checks</h2><div class='wrap'><table><tr><th>check</th><th>status</th><th>detail</th></tr>"]
    for n, s, d in checks:
        parts.append(f"<tr><td>{html.escape(n)}</td><td class='{s}'>{s}</td><td>{html.escape(d)}</td></tr>")
    parts.append("</table></div>")
    if forced:
        parts.append("<p class='fail'>Run was forced despite failed checks. Treat results as unreliable.</p>")
    parts.append(_composition(priv))
    for title, desc, metrics, controls, note in claims(m.get("extra_roles", {}), m.get("extra_states", [])):
        mets = [x for x in metrics + controls if x in priv.columns or x in ref.columns]
        if not mets:
            continue
        parts.append(f"<h2>{html.escape(title)}</h2><p class='mut'>{html.escape(desc)}</p>")
        if note:
            parts.append(f"<p class='mut'>{html.escape(note)}</p>")
        parts.append("<div class='wrap'><table><tr><th>metric</th><th>public n</th><th>public median</th><th>public IQR</th>"
                     "<th>your n</th><th>your median</th><th>your median's percentile among public</th><th>direction</th></tr>")
        parts += [_row(x, ref, priv) for x in mets]
        parts.append("</table></div>")
        parts.append(_pooled_table(mets, ref, priv, name))
        parts.append("<div class='plots'>")
        for x in mets:
            img = _strip(x, ref, priv, name)
            if img:
                parts.append(f"<img alt='{html.escape(x)}' src='data:image/png;base64,{img}'>")
        parts.append("</div>")
    parts.append("<h2>Reference sources</h2><div class='card'>" + ref_sources + "</div>")
    parts.append("<h2>Manifest used</h2><div class='card wrap'><pre><code>" +
                 html.escape(Path(m["_manifest_path"]).read_text()) + "</code></pre></div></main></body></html>")
    (out / "report.html").write_text("\n".join(parts))
    return out / "report.html"
