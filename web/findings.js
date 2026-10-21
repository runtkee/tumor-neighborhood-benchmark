// Top-level findings per dataset, from its per-patient rows. Same rules for public and loaded datasets.
// A direction counts as consistent when at least 70% of patients (and at least 5) agree.
(function (root) {
  "use strict";
  const MIN_N = 5, AGREE = 0.7, BLEED_HIGH = 15;
  const STATE = { prolif: "dividing (Ki-67)", PD1: "PD-1", GZMB: "granzyme B", LAG3: "LAG-3", ICOS: "ICOS", apoptosis: "dying (cleaved PARP)" };
  const vals = (rows, k) => rows.map((r) => r[k]).filter((v) => typeof v === "number" && Number.isFinite(v));
  const median = (a) => { const s = [...a].sort((x, y) => x - y), n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : NaN; };
  const sgn = (x, d = 2) => (x > 0 ? "+" : "") + x.toFixed(d);

  // direction of one metric across patients: "down" | "up" | "mixed" | "few"
  function direction(rows, k) {
    const v = vals(rows, k);
    if (v.length < MIN_N) return { dir: "few", n: v.length };
    const down = v.filter((x) => x < 0).length, up = v.filter((x) => x > 0).length, med = median(v);
    const dir = down / v.length >= AGREE && med < 0 ? "down" : up / v.length >= AGREE && med > 0 ? "up" : "mixed";
    return { dir, n: v.length, down, up, med };
  }

  function proximity(rows, metric, what, controlMetric) {
    const d = direction(rows, metric);
    if (d.dir === "few") return { key: metric, status: "na", text: d.n === 0 ? `${what[0].toUpperCase() + what.slice(1)}: not measured in this dataset.` : `${what[0].toUpperCase() + what.slice(1)}: too few patients to say (${d.n}).` };
    const nn = `${d.dir === "down" ? d.down : d.dir === "up" ? d.up : Math.max(d.down, d.up)} of ${d.n} patients`;
    let text, status;
    if (d.dir === "down") { text = `Fewer ${what} than chance in ${nn} (median ${sgn(d.med)} log2).`; status = "down"; }
    else if (d.dir === "up") { text = `More ${what} than chance in ${nn} (median ${sgn(d.med)} log2).`; status = "up"; }
    else { text = `${what[0].toUpperCase() + what.slice(1)}: no consistent pattern (${d.down} of ${d.n} fewer, median ${sgn(d.med)}).`; status = "mixed"; }
    if (controlMetric && (status === "down" || status === "up")) {
      const c = direction(rows, controlMetric);
      if (c.dir === status) text += " B cells show the same pattern, so it isn't specific to T cells.";
      else if (c.dir !== "few") text += " B cells don't, so it looks specific to T cells.";
    }
    return { key: metric, status, text };
  }

  function contact(rows, prefix, partner, single) {
    const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))].filter((k) => k.startsWith(prefix + "Touch_") && k.endsWith("_diff") && !k.includes("bleed"));
    if (!keys.length) return null;
    const found = [], mixed = [];
    for (const k of keys) {
      const st = k.slice(prefix.length + 6, -5), d = direction(rows, k);
      if (d.dir === "down" || d.dir === "up") found.push(`${STATE[st] || st} ${d.dir === "down" ? "lower" : "higher"} (${d.dir === "down" ? d.down : d.up} of ${d.n} patients, median ${sgn(d.med, 1)} pct-pts)`);
      else if (d.dir === "mixed") mixed.push(st);
    }
    const bleed = vals(rows, `${prefix}Touch_bleed_${prefix}_diff`), b = bleed.length ? median(bleed) : NaN;
    const caution = Number.isFinite(b) && b > BLEED_HIGH ? ` Caution: ${single} signal bleeds into touching T cells (+${b.toFixed(0)} pct-pts), so these contact results are unreliable.` : "";
    if (!found.length && !mixed.length) return { key: prefix + "Touch", status: "na", text: `T cells touching ${partner}: too few patients to say.` };
    return { key: prefix + "Touch", status: found.length ? "change" : "none",
      text: (found.length ? `T cells touching ${partner}: ${found.join("; ")}.` : `T cells touching ${partner}: no consistent change in ${mixed.length} readouts.`) + caution };
  }

  function summarize(rows) {
    if (!rows || !rows.length) return [];
    const out = [
      proximity(rows, "T_near_neut_log2", "T cells near neutrophils", "B_near_neut_log2"),
      proximity(rows, "T_near_mac_log2", "T cells near macrophages"),
      contact(rows, "neut", "neutrophils", "neutrophil"),
      contact(rows, "mac", "macrophages", "macrophage"),
    ];
    if (rows.some((r) => "T_near_neut_subset_log2" in r)) out.push(proximity(rows, "T_near_neut_subset_log2", "T cells near the neutrophil subset"));
    return out.filter(Boolean);
  }

  // which findings change when the pixel size is 20% smaller or larger (only used when it was estimated)
  function sensitivity(base, alternatives) {
    const changed = new Set();
    for (const alt of alternatives) {
      const a = summarize(alt), m = new Map(a.map((f) => [f.key, f.status]));
      for (const f of base) if (m.has(f.key) && m.get(f.key) !== f.status) changed.add(f.key);
    }
    return changed;
  }

  const api = { summarize, sensitivity, direction };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Findings = api;
})(typeof self !== "undefined" ? self : this);
