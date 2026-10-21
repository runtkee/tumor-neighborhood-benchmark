// Browser/Node port of the locked protein recipe (src/bio/recipe_protein.py, v2p.1) and the kit's import + checks.
// Must stay numerically identical to the Python version: tests/test_parity.py runs both on the same data.
(function (root) {
  "use strict";
  const VERSION = "v2p.1";
  const NEAR_UM = 10, ZONE_UM = 30, N_SHUFFLE = 20, MIN_PER_PATIENT = 20, MIN_ROLE = 300;
  const ROLES = ["T", "neut", "mac", "B", "stroma", "tumor", "other"];
  const STATES = ["prolif", "PD1", "GZMB", "LAG3", "ICOS", "apoptosis"];
  const PAIRS = [["T", "neut"], ["B", "neut"], ["T", "mac"]];
  const SOURCES = ["neut", "mac", "stroma"];

  // ---- shared PRNG (bio/prng.py) ----
  function fnv1a(s) {
    let h = 0x811c9dc5;
    const bytes = new TextEncoder().encode(s);
    for (const b of bytes) h = Math.imul(h ^ b, 0x01000193) >>> 0;
    return h >>> 0;
  }
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a) >>> 0;
      t = ((t + (Math.imul(t ^ (t >>> 7), 61 | t) >>> 0)) >>> 0) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function choose(rng, pool, k) {
    const p = Int32Array.from(pool), n = p.length;
    for (let i = 0; i < k; i++) {
      const j = i + Math.floor(rng() * (n - i));
      const tmp = p[i]; p[i] = p[j]; p[j] = tmp;
    }
    return p.subarray(0, k);
  }

  // ---- CSV / TSV (streaming: onRow(fields[]) per line; quotes, CRLF, BOM handled; delimiter auto-detected) ----
  function prepText(t) { return t.charCodeAt(0) === 0xfeff ? t.slice(1) : t; }
  function detectDelim(text) {
    const nl = text.search(/\r?\n/), h = nl > 0 ? text.slice(0, nl) : text;
    const c = { ",": 0, "\t": 0, ";": 0 };
    for (const ch of h) if (ch in c) c[ch]++;
    const [best, n] = Object.entries(c).sort((a, b) => b[1] - a[1])[0];
    return n > 0 ? best : ",";
  }
  function eachCSVRow(text, onRow, delim) {
    const D = (delim || ",").charCodeAt(0);
    let row = [], field = "", q = false, i = 0, n = text.length, start = 0;
    while (i < n) {
      const c = text.charCodeAt(i);
      if (q) {
        if (c === 34) { if (text.charCodeAt(i + 1) === 34) { field += '"'; i += 2; continue; } q = false; i++; continue; }
        field += text[i]; i++; continue;
      }
      if (c === 34) { q = true; i++; continue; }
      if (c === D) { row.push(field); field = ""; i++; continue; }
      if (c === 10 || c === 13) {
        row.push(field); field = ""; if (!(row.length === 1 && row[0] === "")) onRow(row); row = [];
        if (c === 13 && text.charCodeAt(i + 1) === 10) i++;
        i++; continue;
      }
      start = i;  // fast path: copy an unquoted run in one slice
      while (i < n) { const d = text.charCodeAt(i); if (d === D || d === 10 || d === 13 || d === 34) break; i++; }
      field += text.slice(start, i);
    }
    if (field.length || row.length) { row.push(field); onRow(row); }
  }
  function csvHeader(text) {
    text = prepText(text);
    let h = null;
    const nl = text.search(/\r?\n/);
    eachCSVRow(nl > 0 ? text.slice(0, nl) : text, (r) => { if (!h) h = r; }, detectDelim(text));
    return h || [];
  }
  const NA = /^(nan|na|n\/a|#n\/a|null|none|-)$/i;
  const num = (s, dc) => (s === "" || s === undefined || NA.test(s) ? NaN : Number(dc ? s.replace(",", ".") : s));
  const TRUE = /^(1|1\.0+|true|yes|y|pos|positive|\+)$/i;
  const truthy = (s) => s !== undefined && (TRUE.test(String(s).trim()) || (Number(s) > 0));
  const FILE_COL = "__file__";  // virtual column: the file each row came from

  // patient ID from a value: whole, text before the first separator, or the first match of a pattern
  function derive(value, rule) {
    if (!rule || rule.mode === "whole" || !rule.mode) return value;
    if (rule.mode === "before") { const k = value.indexOf(rule.sep); return k > 0 ? value.slice(0, k) : value; }
    if (rule.mode === "regex") { try { const mm = value.match(new RegExp(rule.regex)); return mm ? (mm[1] !== undefined ? mm[1] : mm[0]) : value; } catch (e) { return value; } }
    return value;
  }

  // ---- import: files [{name, text}] (or one CSV string) + config (kit manifest shape) -> cell table ----
  // config extras: columns.x / columns.y may be {min, max} (box centre); phenotype_rules (ordered) replace a label
  // column; patient_rule / field_rule derive patient / image IDs from a column (incl. "__file__"); min_cells_per_field
  // drops tiny images. Note: shuffle seeds use the image ID, so renaming images changes the random draws slightly.
  function importCells(files, m) {
    if (typeof files === "string") files = [{ name: "data.csv", text: files }];
    const c = m.columns, px = Number(m.pixel_size_um);
    const usable = (m.phenotype_rules || []).filter((ru) => (ru.positive || []).length);  // a rule needs a positive marker
    const rules = m.phenotype_rules ? usable : null;
    const cells = { patient: [], field: [], x: [], y: [], raw_type: [], type: [], cols: {} };
    const stateCols = Object.entries(m.states || {}), bleedCols = Object.entries(m.bleed || {}).map(([k, v]) => ["bleed_" + k, v]);
    for (const [k] of stateCols.concat(bleedCols)) cells.cols[k] = [];
    for (const f of files) {
      const text = prepText(f.text), delim = detectDelim(text);
      const dc = delim === ";" && m.decimal_comma !== false;
      let H = null, len = 0, get, coord, extra, ruleIdx, it;
      eachCSVRow(text, (r) => {
        if (!H) {
          H = {}; r.forEach((h, i) => (H[h.trim()] = i)); len = r.length;
          const col = (name) => { if (name === FILE_COL) return -1; if (!(name in H)) throw new Error(`column "${name}" not found in ${f.name}`); return H[name]; };
          get = (idx, row) => (idx === -1 ? f.name : row[idx]);
          const ax = (spec) => (spec && typeof spec === "object" ? [col(spec.min), col(spec.max)] : [col(spec)]);
          coord = { x: ax(c.x), y: ax(c.y) };
          extra = stateCols.concat(bleedCols).map(([k, name]) => [k, col(name)]);
          ruleIdx = rules ? rules.map((ru) => ({ role: ru.role, pos: (ru.positive || []).map(col), neg: (ru.negative || []).map(col) })) : null;
          it = rules ? null : (c.type ? col(c.type) : null);
          H.__ip = col(c.patient); H.__if = col(c.field);
          return;
        }
        if (r.length < len) return;
        let raw, role;
        if (ruleIdx) {
          const k = ruleIdx.findIndex((ru) => ru.pos.every((i) => truthy(get(i, r))) && ru.neg.every((i) => !truthy(get(i, r))));
          raw = k >= 0 ? `rule ${k + 1}: ${ruleIdx[k].role}` : "no rule matched";
          role = k >= 0 ? ruleIdx[k].role : "other";
        } else {
          raw = it === null ? "(no cell type column)" : get(it, r);
          role = Object.prototype.hasOwnProperty.call(m.roles || {}, raw) ? m.roles[raw] : "unmapped";
        }
        if (role === "drop") return;
        const v = (ix) => (ix.length === 2 ? (num(get(ix[0], r), dc) + num(get(ix[1], r), dc)) / 2 : num(get(ix[0], r), dc));
        cells.patient.push(derive(get(H.__ip, r), m.patient_rule)); cells.field.push(derive(get(H.__if, r), m.field_rule));
        cells.x.push(v(coord.x) * px); cells.y.push(v(coord.y) * px);
        cells.raw_type.push(raw); cells.type.push(role);
        for (const [k, i] of extra) cells.cols[k].push(num(get(i, r), dc));
      }, delim);
    }
    const minF = Number(m.min_cells_per_field || 0);
    if (minF > 0) {
      const n = {}; for (const f of cells.field) n[f] = (n[f] || 0) + 1;
      const keep = cells.field.map((f) => n[f] >= minF);
      const filt = (a) => a.filter((_, i) => keep[i]);
      for (const k of ["patient", "field", "x", "y", "raw_type", "type"]) cells[k] = filt(cells[k]);
      for (const k of Object.keys(cells.cols)) cells.cols[k] = filt(cells.cols[k]);
      cells.dropped_fields = Object.values(n).filter((x) => x < minF).length;
    }
    cells.n = cells.type.length;
    return cells;
  }

  // ---- spatial grid: nearest distance from each source point to a target set ----
  function makeGrid(xs, ys, idx, size) {
    const g = new Map();
    for (const i of idx) {
      const k = Math.floor(xs[i] / size) + "," + Math.floor(ys[i] / size);
      let b = g.get(k); if (!b) g.set(k, (b = [])); b.push(i);
    }
    return { g, size };
  }
  // exact nearest distance if it is < maxR, else Infinity (all recipe comparisons use thresholds <= maxR)
  function nearestWithin(grid, xs, ys, i, maxR, skipSelf) {
    const { g, size } = grid, cx = Math.floor(xs[i] / size), cy = Math.floor(ys[i] / size);
    const rr = Math.ceil(maxR / size);
    let best = Infinity;
    for (let dx = -rr; dx <= rr; dx++) for (let dy = -rr; dy <= rr; dy++) {
      const b = g.get(cx + dx + "," + (cy + dy)); if (!b) continue;
      for (const j of b) {
        if (skipSelf && j === i) continue;
        const ddx = xs[i] - xs[j], ddy = ys[i] - ys[j];
        const d = Math.sqrt(ddx * ddx + ddy * ddy);
        if (d < best) best = d;
      }
    }
    return best < maxR ? best : Infinity;
  }

  // ---- percentile rank within field, ties averaged, NaN excluded (pandas rank(pct=True)) ----
  function pctRankByField(values, fieldIdx) {
    const out = new Float64Array(values.length).fill(NaN);
    for (const idx of fieldIdx.values()) {
      const valid = idx.filter((i) => !Number.isNaN(values[i]));
      valid.sort((a, b) => values[a] - values[b]);
      const n = valid.length;
      let s = 0;
      while (s < n) {
        let e = s; while (e + 1 < n && values[valid[e + 1]] === values[valid[s]]) e++;
        const avg = (s + 1 + e + 1) / 2;
        for (let k = s; k <= e; k++) out[valid[k]] = (avg / n) * 100;
        s = e + 1;
      }
    }
    return out;
  }

  function groupIndex(keys) {
    const m = new Map();
    keys.forEach((k, i) => { let b = m.get(k); if (!b) m.set(k, (b = [])); b.push(i); });
    return m;
  }
  const mean = (a) => { let s = 0; for (const v of a) s += v; return s / a.length; };

  // ---- recipe: one row per patient ----
  function patientRows(cells, m) {
    const extra = Object.keys(m.extra_roles || {});
    const present = new Set(cells.type);
    const roles = [...present].filter((r) => ROLES.includes(r) || extra.includes(r));
    const states = Object.keys(m.states || {}).filter((s) => s in cells.cols && cells.cols[s].some((v) => !Number.isNaN(v)));
    const bleed = {}; for (const r of Object.keys(m.bleed || {})) if (("bleed_" + r) in cells.cols) bleed[r] = "bleed_" + r;
    const pairs = PAIRS.concat(extra.map((r) => ["T", r]));
    const sources = SOURCES.concat(extra);
    const byField = groupIndex(cells.field);
    const pct = {};
    for (const col of states.concat(Object.values(bleed))) pct[col] = pctRankByField(cells.cols[col], byField);
    const xs = cells.x, ys = cells.y;
    const byPF = groupIndex(cells.patient.map((p, i) => p + "\u0000" + cells.field[i]));
    const acc = new Map();
    for (const [key, ix] of byPF) {
      const [pat, field] = key.split("\u0000");
      if (!acc.has(pat)) acc.set(pat, {});
      const a = acc.get(pat);
      const R = {}; for (const r of roles) R[r] = ix.filter((i) => cells.type[i] === r);
      for (const [tgt, src] of pairs) {
        if (!R[tgt] || !R[src] || !R[tgt].length || !R[src].length) continue;
        const countNear = (srcIdx) => {
          const grid = makeGrid(xs, ys, srcIdx, ZONE_UM);
          let c = 0; for (const i of R[tgt]) if (nearestWithin(grid, xs, ys, i, ZONE_UM, false) < ZONE_UM) c++;
          return c;
        };
        const real = countNear(R[src]);
        const tgtSet = new Set(R[tgt]);
        const pool = []; ix.forEach((i, pos) => { if (!tgtSet.has(i)) pool.push(pos); });
        const rng = mulberry32(fnv1a(`${pat}|${field}|${tgt}|${src}`));
        let sum = 0;
        for (let s = 0; s < N_SHUFFLE; s++) sum += countNear(Array.from(choose(rng, pool, R[src].length), (pos) => ix[pos]));
        const k = `${tgt}_near_${src}`;
        a[k + "_real"] = (a[k + "_real"] || 0) + real;
        a[k + "_null"] = (a[k + "_null"] || 0) + sum / N_SHUFFLE;
      }
      for (const src of sources) {
        if (!R[src] || !R.T || !R.T.length || !R[src].length) continue;
        const grid = makeGrid(xs, ys, R[src], ZONE_UM);
        const cols = states.concat(src in bleed ? [bleed[src]] : []);
        for (const i of R.T) {
          const d = nearestWithin(grid, xs, ys, i, ZONE_UM, false);
          const near = d < NEAR_UM, far = d >= ZONE_UM;
          if (!near && !far) continue;
          for (const col of cols) {
            const k = `${src}Touch_${col}`;
            (a[k + "_near"] ||= []); (a[k + "_far"] ||= []);
            (near ? a[k + "_near"] : a[k + "_far"]).push(pct[col][i]);
          }
        }
        // make sure keys exist even when no T cell was near/far in this field (Python extends with empty lists)
        for (const col of cols) { const k = `${src}Touch_${col}`; (a[k + "_near"] ||= []); (a[k + "_far"] ||= []); }
      }
    }
    const rows = [];
    for (const [pat, a] of acc) {
      const r = { dataset: m.dataset, patient: pat, recipe: VERSION };
      for (const k of Object.keys(a)) {
        if (k.endsWith("_real")) { const b = k.slice(0, -5); r[b + "_log2"] = Math.log2((a[k] + 1) / (a[b + "_null"] + 1)); }
        if (k.endsWith("_near")) {
          const b = k.slice(0, -5), n = a[k], f = a[b + "_far"];
          r[b + "_n"] = n.length;
          r[b + "_diff"] = n.length >= MIN_PER_PATIENT && f.length >= MIN_PER_PATIENT ? mean(n) - mean(f) : NaN;
        }
      }
      rows.push(r);
    }
    return rows;
  }

  // ---- checks (src/bio/kit/checks.py) ----
  function median(a) { const s = Float64Array.from(a).sort(); const n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : NaN; }
  // median over fields of the median nearest-cell distance (same units as cells.x / cells.y)
  function medianSpacing(cells) {
    const byField = groupIndex(cells.field);
    const nn = [];
    for (const idx of byField.values()) {
      if (idx.length <= 20) continue;
      const finite = idx.filter((i) => Number.isFinite(cells.x[i]) && Number.isFinite(cells.y[i]));
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (const i of finite) { x0 = Math.min(x0, cells.x[i]); x1 = Math.max(x1, cells.x[i]); y0 = Math.min(y0, cells.y[i]); y1 = Math.max(y1, cells.y[i]); }
      const span = Math.max(x1 - x0, y1 - y0);
      const size = Math.max(1e-6, Math.sqrt((span * span) / Math.max(finite.length, 1)) * 2);  // ~2 cells per grid bin, any unit
      const grid = makeGrid(cells.x, cells.y, finite, size);
      const d = finite.map((i) => { let r = size, v = Infinity; while (v === Infinity && r < size * 4096) { v = nearestWithin(grid, cells.x, cells.y, i, r, true); r *= 2; } return v; });
      nn.push(median(d));
    }
    return nn.length ? median(nn) : NaN;
  }
  function runChecks(cells, m) {
    const out = [];
    const byField = groupIndex(cells.field);
    const med = medianSpacing(cells);
    out.push(["cell spacing (µm)", med >= 5 && med <= 15 ? "pass" : "fail", `median nearest-cell distance ${med.toFixed(1)} µm (expected 5–15). Wrong pixel size if outside.`]);
    const counts = {}; for (const t of cells.type) counts[t] = (counts[t] || 0) + 1;
    for (const r of ["T", "neut", "mac", "B"]) {
      const n = counts[r] || 0;
      out.push([`role ${r}`, n >= MIN_ROLE ? "pass" : (n === 0 && r === "T" ? "fail" : "warn"),
        `${n} cells` + (n >= MIN_ROLE ? "" : ` (<${MIN_ROLE}: tests involving ${r} will be weak or skipped)`)]);
    }
    const unm = [...new Set(cells.raw_type.filter((_, i) => cells.type[i] === "unmapped"))].sort();
    out.push(["unmapped labels", unm.length ? "warn" : "pass", unm.length ? `${unm.length} labels not mapped (treated as 'other'): ${unm.slice(0, 12).join(", ")}` : "all labels mapped"]);
    for (const col of Object.keys(cells.cols).filter((k) => k.startsWith("bleed_"))) {
      const r = col.slice(6);
      if (!counts[r] || !counts.T) continue;
      const p = pctRankByField(cells.cols[col], byField);
      const inR = [], inT = [];
      for (let i = 0; i < cells.n; i++) { if (Number.isNaN(p[i])) continue; if (cells.type[i] === r) inR.push(p[i]); else if (cells.type[i] === "T") inT.push(p[i]); }
      const gap = median(inR) - median(inT);
      out.push([`marker sanity: ${r}`, gap > 20 ? "pass" : "fail", `${r}'s own marker: median percentile in ${r} minus in T = ${gap.toFixed(0)} (expect >20; low = mislabeled role or wrong marker)`]);
    }
    const have = Object.keys(m.states || {}).filter((s) => s in cells.cols && cells.cols[s].some((v) => !Number.isNaN(v)));
    out.push(["T-state markers", have.length ? "pass" : "warn", `available: ${have.join(", ") || "none"}`]);
    const pats = new Set(cells.patient).size;
    out.push(["patients", pats >= 5 ? "pass" : "warn", `${pats} patients, ${byField.size} fields (<5 patients: cross-patient statistics unreliable)`]);
    return out;
  }

  const api = { VERSION, ROLES, STATES, fnv1a, mulberry32, choose, eachCSVRow, csvHeader, detectDelim, prepText, derive, truthy, FILE_COL, importCells, patientRows, runChecks, medianSpacing, pctRankByField };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Recipe = api;
})(typeof self !== "undefined" ? self : this);
