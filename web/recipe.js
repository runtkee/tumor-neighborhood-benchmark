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

  // ---- CSV / TSV parser: stateful, so it can be fed a file in chunks (quotes, CRLF and "" escapes may span chunks) ----
  function prepText(t) { return t.charCodeAt(0) === 0xfeff ? t.slice(1) : t; }
  function detectDelim(text) {
    const nl = text.search(/\r?\n/), h = nl > 0 ? text.slice(0, nl) : text;
    const c = { ",": 0, "\t": 0, ";": 0 };
    for (const ch of h) if (ch in c) c[ch]++;
    const [best, n] = Object.entries(c).sort((a, b) => b[1] - a[1])[0];
    return n > 0 ? best : ",";
  }
  function createParser(delim, onRow) {
    const D = (delim || ",").charCodeAt(0);
    let row = [], field = "", q = false, quotePending = false, skipLF = false;
    const emit = () => { row.push(field); field = ""; if (!(row.length === 1 && row[0] === "")) onRow(row); row = []; };
    function feed(text) {
      let i = 0; const n = text.length;
      if (skipLF && n && text.charCodeAt(0) === 10) i = 1;
      skipLF = false;
      if (quotePending) {  // previous chunk ended on a quote inside a quoted field
        quotePending = false;
        if (i < n && text.charCodeAt(i) === 34) { field += '"'; i++; } else q = false;
      }
      while (i < n) {
        const c = text.charCodeAt(i);
        if (q) {
          if (c === 34) {
            if (i + 1 === n) { quotePending = true; i++; continue; }
            if (text.charCodeAt(i + 1) === 34) { field += '"'; i += 2; continue; }
            q = false; i++; continue;
          }
          const st = i; while (i < n && text.charCodeAt(i) !== 34) i++;
          field += text.slice(st, i); continue;
        }
        if (c === 34) { q = true; i++; continue; }
        if (c === D) { row.push(field); field = ""; i++; continue; }
        if (c === 10 || c === 13) {
          emit();
          if (c === 13) { if (i + 1 < n) { if (text.charCodeAt(i + 1) === 10) i++; } else skipLF = true; }
          i++; continue;
        }
        const st = i;  // fast path: copy an unquoted run in one slice
        while (i < n) { const d = text.charCodeAt(i); if (d === D || d === 10 || d === 13 || d === 34) break; i++; }
        field += text.slice(st, i);
      }
    }
    function end() { if (quotePending) { quotePending = false; q = false; } if (field.length || row.length) emit(); }
    return { feed, end };
  }
  function eachCSVRow(text, onRow, delim) { const p = createParser(delim, onRow); p.feed(text); p.end(); }
  function csvHeader(text) {
    text = prepText(text);
    let h = null;
    const nl = text.search(/\r?\n/);
    eachCSVRow(nl > 0 ? text.slice(0, nl) : text, (r) => { if (!h) h = r; }, detectDelim(text));
    return (h || []).map((x) => x.trim());
  }
  const NA = /^(nan|na|n\/a|#n\/a|null|none|-)$/i;
  const num = (s, dc) => (s === "" || s === undefined || NA.test(s) ? NaN : Number(dc ? s.replace(",", ".") : s));
  const TRUE = /^(1|1\.0+|true|yes|y|pos|positive|\+)$/i;
  const truthy = (s) => s !== undefined && (TRUE.test(String(s).trim()) || (Number(s) > 0));
  const FILE_COL = "__file__";  // virtual column: the file each row came from

  // What kind of file is this, from its first bytes? Returns null if it looks like text.
  function sniff(bytes, name) {
    const b = bytes, s = (k) => String.fromCharCode(...b.slice(0, k));
    if (b[0] === 0x50 && b[1] === 0x4b) return { title: `${name} is an Excel workbook or a zip archive`, fix: "Open it in Excel and use File > Save As > CSV, or unzip it, then load the CSV." };
    if (b[0] === 0xd0 && b[1] === 0xcf) return { title: `${name} is an old-format Excel (.xls) file`, fix: "Save it as CSV from Excel, then load the CSV." };
    if (s(4) === "\x89HDF") return { title: `${name} is an HDF5 file (for example .h5ad)`, fix: "This page reads cell tables (CSV or tab-separated). For .h5ad files use the Python kit (python -m bio.kit), or export the cell table to CSV." };
    if (s(4) === "%PDF") return { title: `${name} is a PDF`, fix: "Load the cell table (CSV) exported from your image analysis software." };
    if (b[0] === 0x1f && b[1] === 0x8b) return { title: `${name} is gzip-compressed`, fix: "Decompress it first (double-click, or gunzip), then load the CSV." };
    if (b.slice(0, 4096).includes(0)) return { title: `${name} is not a text file`, fix: "Export the cell table as CSV or tab-separated text." };
    return null;
  }

  // patient / image ID from a value: whole, text before the first separator, or the first match of a pattern
  function derive(value, rule) {
    if (!rule || rule.mode === "whole" || !rule.mode) return value;
    if (rule.mode === "before") { const k = value.indexOf(rule.sep); return k > 0 ? value.slice(0, k) : value; }
    if (rule.mode === "regex") { try { const mm = value.match(new RegExp(rule.regex)); return mm ? (mm[1] !== undefined ? mm[1] : mm[0]) : value; } catch (e) { return value; } }
    return value;
  }

  // growable Float64Array: millions of cells without the overhead of JS number arrays
  function F64() {
    let a = new Float64Array(4096), n = 0;
    return { push(v) { if (n === a.length) { const b = new Float64Array(a.length * 2); b.set(a); a = b; } a[n++] = v; }, done() { return a.slice(0, n); } };
  }
  function interner() { const m = new Map(); return (s) => { let v = m.get(s); if (v === undefined) { m.set(s, s); v = s; } return v; }; }
  const DUP_SAMPLE = 300000;

  // ---- import: one importer, fed files row by row (from a whole string or a stream) ----
  // config (kit manifest shape) extras: columns.x / columns.y may be {min, max} (box centre); phenotype_rules (ordered)
  // replace a label column; patient_rule / field_rule derive IDs from a column (incl. "__file__"); min_cells_per_field
  // drops tiny images. Cells with missing coordinates are dropped (and counted). Shuffle seeds use the image ID.
  function createImporter(m) {
    const c = m.columns, px = Number(m.pixel_size_um);
    const rules = m.phenotype_rules ? m.phenotype_rules.filter((ru) => (ru.positive || []).length) : null;  // a rule needs a positive marker
    const stateCols = Object.entries(m.states || {}), bleedCols = Object.entries(m.bleed || {}).map(([k, v]) => ["bleed_" + k, v]);
    const valCols = stateCols.concat(bleedCols);
    const P = [], Fd = [], RT = [], TY = [], X = F64(), Y = F64(), V = {};
    for (const [k] of valCols) V[k] = F64();
    const I = interner();
    const diag = { files: [], rows: 0, short: 0, coordBad: 0, idEmpty: { patient: 0, field: 0 }, dropped: 0, dup: 0, dupChecked: 0,
                   colStats: {}, ruleCols: {}, ruleHits: rules ? rules.map(() => 0) : null, noRule: 0, coordEx: [] };
    for (const [k, name] of valCols) diag.colStats[k] = { col: name, n: 0, nan: 0, textEx: [], distinct: new Set(), min: Infinity, max: -Infinity };
    const seen = new Set();
    let cur = null;
    const stat = (st, raw, v) => {
      if (Number.isNaN(v)) { st.nan++; if (raw !== "" && raw !== undefined && !NA.test(raw) && st.textEx.length < 3 && !st.textEx.includes(raw)) st.textEx.push(raw); return; }
      st.n++; if (v < st.min) st.min = v; if (v > st.max) st.max = v; if (st.distinct.size < 12) st.distinct.add(v);
    };
    function startFile(name, delim) {
      cur = { name, delim, dc: delim === ";" && m.decimal_comma !== false, H: null, info: { name, rows: 0, short: 0, missing: [] } };
      diag.files.push(cur.info);
    }
    function header(r) {
      const H = {}; r.forEach((h, i) => { const k = h.replace(/^﻿/, "").trim(); if (!(k in H)) H[k] = i; });
      const missing = [];
      const col = (name) => { if (name === FILE_COL) return -1; if (!(name in H)) { missing.push(name); return -2; } return H[name]; };
      const ax = (spec) => (spec && typeof spec === "object" ? [col(spec.min), col(spec.max)] : [col(spec)]);
      cur.len = r.length;
      cur.ip = col(c.patient); cur.if = col(c.field); cur.cx = ax(c.x); cur.cy = ax(c.y);
      cur.extra = valCols.map(([k, name]) => [k, col(name)]);
      cur.rules = rules ? rules.map((ru) => ({ role: ru.role, pos: (ru.positive || []).map((n) => [n, col(n)]), neg: (ru.negative || []).map((n) => [n, col(n)]) })) : null;
      cur.it = rules ? null : (c.type ? col(c.type) : null);
      cur.info.missing = [...new Set(missing)];
      if (cur.info.missing.length) {
        const e = new Error(`${cur.name} is missing column${cur.info.missing.length > 1 ? "s" : ""}: ${cur.info.missing.map((x) => `"${x}"`).join(", ")}`);
        e.fix = "All files must have the mapped columns. Check that every file is the same export type, or re-map the columns using a file that has them.";
        throw e;
      }
      cur.H = H;
    }
    function row(r) {
      if (!cur.H) return header(r);
      diag.rows++; cur.info.rows++;
      if (r.length < cur.len) { diag.short++; cur.info.short++; return; }
      const g = (i) => (i === -1 ? cur.name : r[i]);
      let raw, role;
      if (cur.rules) {
        let k = -1;
        for (let j = 0; j < cur.rules.length && k < 0; j++) {
          const ru = cur.rules[j];
          if (ru.pos.every(([, i]) => truthy(g(i))) && ru.neg.every(([, i]) => !truthy(g(i)))) k = j;
        }
        if (diag.rows <= 200000) for (const ru of cur.rules) for (const [name, i] of ru.pos.concat(ru.neg)) {
          const st = (diag.ruleCols[name] ||= { n: 0, pos: 0, distinct: new Set() }), v = g(i);
          if (v !== "" && v !== undefined) { st.n++; if (truthy(v)) st.pos++; if (st.distinct.size < 6) st.distinct.add(String(v).trim()); }
        }
        if (k >= 0) { diag.ruleHits[k]++; raw = `rule ${k + 1}: ${cur.rules[k].role}`; role = cur.rules[k].role; } else { diag.noRule++; raw = "no rule matched"; role = "other"; }
      } else {
        raw = cur.it === null ? "(no cell type column)" : g(cur.it);
        role = Object.prototype.hasOwnProperty.call(m.roles || {}, raw) ? m.roles[raw] : "unmapped";
      }
      if (role === "drop") { diag.dropped++; return; }
      const dc = cur.dc;
      const v = (ix) => (ix.length === 2 ? (num(g(ix[0]), dc) + num(g(ix[1]), dc)) / 2 : num(g(ix[0]), dc));
      const x = v(cur.cx) * px, y = v(cur.cy) * px;
      if (!Number.isFinite(x) || !Number.isFinite(y)) {
        diag.coordBad++;
        if (diag.coordEx.length < 3) diag.coordEx.push(cur.cx.concat(cur.cy).map((i) => g(i)).join(" / "));
        return;
      }
      const pat = derive(g(cur.ip), m.patient_rule), fld = derive(g(cur.if), m.field_rule);
      if (pat === "" || pat === undefined) diag.idEmpty.patient++;
      if (fld === "" || fld === undefined) diag.idEmpty.field++;
      if (diag.dupChecked < DUP_SAMPLE) { diag.dupChecked++; const key = pat + "\u0000" + fld + "\u0000" + x + "\u0000" + y; if (seen.has(key)) diag.dup++; else seen.add(key); }
      P.push(I(pat)); Fd.push(I(fld)); RT.push(I(raw)); TY.push(I(role)); X.push(x); Y.push(y);
      for (const [k, i] of cur.extra) { const rv = g(i), val = num(rv, dc); V[k].push(val); stat(diag.colStats[k], rv, val); }
    }
    function finish() {
      seen.clear();
      let cells = { patient: P, field: Fd, x: X.done(), y: Y.done(), raw_type: RT, type: TY, cols: {}, diag };
      for (const [k] of valCols) cells.cols[k] = V[k].done();
      const minF = Number(m.min_cells_per_field || 0);
      diag.droppedFields = 0;
      if (minF > 0) {
        const n = new Map(); for (let i = 0; i < Fd.length; i++) { const key = P[i] + "\u0000" + Fd[i]; n.set(key, (n.get(key) || 0) + 1); }
        const keep = new Uint8Array(Fd.length); let kept = 0;
        for (let i = 0; i < Fd.length; i++) if (n.get(P[i] + "\u0000" + Fd[i]) >= minF) { keep[i] = 1; kept++; }
        diag.droppedFields = [...n.values()].filter((x) => x < minF).length;
        diag.droppedSmall = Fd.length - kept;
        if (kept < Fd.length) {
          const fa = (a) => a.filter((_, i) => keep[i]);
          cells = { ...cells, patient: fa(P), field: fa(Fd), raw_type: fa(RT), type: fa(TY), x: cells.x.filter((_, i) => keep[i]), y: cells.y.filter((_, i) => keep[i]) };
          for (const [k] of valCols) cells.cols[k] = cells.cols[k].filter((_, i) => keep[i]);
        }
      }
      cells.n = cells.type.length;
      return cells;
    }
    return { startFile, row, finish, diag };
  }
  // whole-string input (Node, tests, small files): files [{name, text}] or one CSV string
  function importCells(files, m) {
    if (typeof files === "string") files = [{ name: "data.csv", text: files }];
    const imp = createImporter(m);
    for (const f of files) { const text = prepText(f.text), d = detectDelim(text); imp.startFile(f.name, d); eachCSVRow(text, imp.row, d); }
    return imp.finish();
  }
  // streaming input (browser File / Blob): reads in chunks, never holds a whole file as one string
  async function importFiles(files, m, onProgress) {
    const imp = createImporter(m);
    const total = files.reduce((a, f) => a + f.size, 0); let doneBytes = 0;
    for (const f of files) {
      const head = prepText(await f.slice(0, 1 << 16).text()), d = detectDelim(head);
      imp.startFile(f.name, d);
      const parser = createParser(d, imp.row);
      const reader = f.stream().getReader(), dec = new TextDecoder("utf-8");
      for (;;) {
        const { value, done } = await reader.read();
        if (done) break;
        doneBytes += value.length;
        parser.feed(dec.decode(value, { stream: true }));
        if (onProgress) onProgress(doneBytes / total);
      }
      parser.feed(dec.decode()); parser.end();
    }
    return imp.finish();
  }

  // ---- spatial grid: nearest distance from each source point to a target set ----
  // numeric bin keys (cx * 2^22 + cy) instead of strings; bins only choose candidates, distances are exact
  const GK = 4194304;
  function makeGrid(xs, ys, idx, size) {
    const g = new Map();
    for (const i of idx) {
      const k = Math.floor(xs[i] / size) * GK + Math.floor(ys[i] / size);
      let b = g.get(k); if (!b) g.set(k, (b = [])); b.push(i);
    }
    return { g, size };
  }
  // exact nearest distance if it is < maxR, else Infinity (all recipe comparisons use thresholds <= maxR)
  function nearestWithin(grid, xs, ys, i, maxR, skipSelf) {
    const { g, size } = grid, cx = Math.floor(xs[i] / size), cy = Math.floor(ys[i] / size);
    const rr = Math.ceil(maxR / size), xi = xs[i], yi = ys[i];
    let best = Infinity;
    for (let dx = -rr; dx <= rr; dx++) {
      const base = (cx + dx) * GK + cy;
      for (let dy = -rr; dy <= rr; dy++) {
        const b = g.get(base + dy); if (!b) continue;
        for (let q = 0; q < b.length; q++) {
          const j = b[q];
          if (skipSelf && j === i) continue;
          const ddx = xi - xs[j], ddy = yi - ys[j];
          const d = Math.sqrt(ddx * ddx + ddy * ddy);
          if (d < best) best = d;
        }
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

  // an image is identified by patient + image name, so "Core 1" in two patients stays two images
  function imageKeys(cells) { return cells.patient.map((p, i) => p + "\u0000" + cells.field[i]); }
  // a readout is usable if it has at least two different numeric values (constant columns carry no information)
  function usable(cells, k) {
    const a = cells.cols[k]; if (!a) return false;
    let first = NaN;
    for (const v of a) { if (Number.isNaN(v)) continue; if (Number.isNaN(first)) first = v; else if (v !== first) return true; }
    return false;
  }
  function groupIndex(keys) {
    const m = new Map();
    keys.forEach((k, i) => { let b = m.get(k); if (!b) m.set(k, (b = [])); b.push(i); });
    return m;
  }
  const mean = (a) => { let s = 0; for (const v of a) s += v; return s / a.length; };

  // ---- recipe: one row per patient ----
  function patientRows(cells, m, onProgress) {
    const extra = Object.keys(m.extra_roles || {});
    const present = new Set(cells.type);
    const roles = [...present].filter((r) => ROLES.includes(r) || extra.includes(r));
    const states = Object.keys(m.states || {}).filter((s) => usable(cells, s));
    const bleed = {}; for (const r of Object.keys(m.bleed || {})) if (usable(cells, "bleed_" + r)) bleed[r] = "bleed_" + r;
    const pairs = PAIRS.concat(extra.map((r) => ["T", r]));
    const sources = SOURCES.concat(extra);
    const byField = groupIndex(imageKeys(cells));
    const pct = {};
    for (const col of states.concat(Object.values(bleed))) pct[col] = pctRankByField(cells.cols[col], byField);
    const xs = cells.x, ys = cells.y;
    const byPF = groupIndex(imageKeys(cells));
    const acc = new Map();
    let done = 0;
    for (const [key, ix] of byPF) {
      if (onProgress) onProgress(done++ / byPF.size);
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
    const byField = groupIndex(imageKeys(cells));
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
  // Each check: [name, status (pass | warn | fail), what we found, how to fix it]. A fail blocks the comparison.
  function runChecks(cells, m) {
    const out = [], d = cells.diag || {}, add = (n, s, detail, fix = "") => out.push([n, s, detail, fix]);
    const pct = (a, b) => (b ? (100 * a) / b : 0), f1 = (x) => x.toFixed(1);
    // ---- rows ----
    if (d.rows !== undefined) {
      add("rows read", d.rows > 0 ? "pass" : "fail", `${d.rows.toLocaleString()} data rows in ${d.files.length} file${d.files.length > 1 ? "s" : ""}`,
        d.rows > 0 ? "" : "The file has a header but no data rows. Check that you exported the cell (object) table, not a summary.");
      if (d.short) {
        const worst = d.files.filter((f) => f.short).map((f) => `${f.name}: ${f.short.toLocaleString()}`).slice(0, 3).join("; ");
        add("incomplete rows", pct(d.short, d.rows) > 5 ? "fail" : "warn", `${d.short.toLocaleString()} rows (${f1(pct(d.short, d.rows))}%) have fewer columns than the header and were skipped (${worst})`,
          "Usually a truncated export, or a delimiter inside unquoted text. Re-export the file; if it was edited in Excel, save as CSV UTF-8.");
      }
      if (d.coordBad) add("missing coordinates", pct(d.coordBad, d.rows) > 5 ? "fail" : "warn",
        `${d.coordBad.toLocaleString()} cells (${f1(pct(d.coordBad, d.rows))}%) have no usable x/y and were left out. Example values: ${d.coordEx.join(" | ")}`,
        "Check that x and y point to numeric position columns (for HALO, XMin/XMax and YMin/YMax).");
      if (d.idEmpty.patient || d.idEmpty.field) add("empty IDs", pct(d.idEmpty.patient + d.idEmpty.field, cells.n) > 5 ? "fail" : "warn",
        `${d.idEmpty.patient.toLocaleString()} cells without a patient ID, ${d.idEmpty.field.toLocaleString()} without an image ID`,
        "Choose different patient / image columns, or take them from the file name.");
    }
    if (!cells.n) { add("cells after mapping", "fail", "No cells are left to analyse.", "Check the cell-type mapping (are all labels set to Exclude?) and the minimum cells per image."); return out; }
    if (d.dup && d.dupChecked) add("duplicate cells", pct(d.dup, d.dupChecked) > 1 ? "fail" : "warn",
      `${d.dup.toLocaleString()} of ${d.dupChecked.toLocaleString()} checked cells repeat another cell exactly (same patient, image and position)`,
      "The same file may be loaded twice, or two exports overlap. Load each slide once.");
    if (d.droppedFields) add("small images skipped", "warn", `${d.droppedFields} images with fewer than ${m.min_cells_per_field} cells were left out (${(d.droppedSmall || 0).toLocaleString()} cells)`,
      "Lower the minimum if these images are real.");
    // ---- geometry ----
    const byImg = groupIndex(imageKeys(cells));
    const spans = [];
    for (const idx of byImg.values()) {
      let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
      for (const i of idx) { const x = cells.x[i], y = cells.y[i]; if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y; }
      spans.push(Math.max(x1 - x0, y1 - y0));
    }
    const span = median(spans), med = medianSpacing(cells);
    add("cell spacing", med >= 5 && med <= 15 ? "pass" : "fail", `median distance to the nearest cell: ${Number.isFinite(med) ? f1(med) : "?"} µm (expected 5–15)`,
      med === 0 ? "Many cells sit exactly on top of each other: check for duplicated files or rows." :
      med < 5 ? `Distances look too small. If your coordinates are in pixels, the pixel size may be too low; if they are already in µm, set pixel size to 1. Implied pixel size: ${f1(8 / (med / Number(m.pixel_size_um)))} µm.` :
      med > 15 ? `Distances look too large. Check the pixel size (implied: ${(8 / (med / Number(m.pixel_size_um))).toFixed(3)} µm), and that x and y are positions, not areas or IDs.` : "");
    add("image size", span >= 100 && span <= 100000 ? "pass" : "fail", `median image width ${span >= 1000 ? f1(span / 1000) + " mm" : f1(span) + " µm"} across ${byImg.size} images`,
      span < 100 ? "Images look tiny. The image column may split one image into many (each cell its own image?), or the pixel size is wrong." :
      span > 100000 ? "Images look larger than 10 cm. The pixel size or coordinate columns are probably wrong, or several slides share one image ID." : "");
    const imgNames = new Map(); for (const k of byImg.keys()) { const f = k.split("\u0000")[1]; imgNames.set(f, (imgNames.get(f) || 0) + 1); }
    const shared = [...imgNames.entries()].filter(([, n]) => n > 1);
    if (shared.length) add("image names repeat across patients", "pass", `${shared.length} image names occur in more than one patient (e.g. "${shared[0][0]}"). They are kept as separate images.`);
    // ---- roles ----
    const counts = {}; for (const t of cells.type) counts[t] = (counts[t] || 0) + 1;
    for (const r of ["T", "neut", "mac", "B"]) {
      const n = counts[r] || 0;
      add(`cells: ${({ T: "T cells", neut: "neutrophils", mac: "macrophages", B: "B cells" })[r]}`, n >= MIN_ROLE ? "pass" : (r === "T" && n === 0 ? "fail" : "warn"),
        `${n.toLocaleString()} cells` + (n >= MIN_ROLE ? "" : ` (fewer than ${MIN_ROLE}: tests involving them will be weak or skipped)`),
        n >= MIN_ROLE ? "" : (r === "T" && n === 0 ? "Every test needs T cells. Map at least one label (or rule) to T cell." : "Map the matching label(s) or add a rule, if your panel has this cell type."));
    }
    const otherShare = pct((counts.other || 0) + (counts.unmapped || 0), cells.n);
    if (otherShare > 80) add("most cells unassigned", "warn", `${f1(otherShare)}% of cells are Other or unmapped`, "Check the role mapping or the marker rules.");
    const unm = [...new Set(cells.raw_type.filter((_, i) => cells.type[i] === "unmapped"))].sort();
    add("labels", unm.length ? "warn" : "pass", unm.length ? `${unm.length} labels not mapped (treated as Other): ${unm.slice(0, 10).join(", ")}` : "all labels mapped",
      unm.length ? "Assign them a role if they matter." : "");
    if (d.ruleHits) {
      const rules = (m.phenotype_rules || []).filter((ru) => (ru.positive || []).length);
      rules.forEach((ru, k) => {
        const n = d.ruleHits[k];
        if (n === 0) add(`rule ${k + 1} (${ru.role})`, "warn", "matched no cells", "Check the chosen marker columns and whether a rule above it already takes these cells.");
        else if (pct(n, d.rows) > 90) add(`rule ${k + 1} (${ru.role})`, "fail", `matched ${f1(pct(n, d.rows))}% of all cells`, "Almost every cell matches: the marker column is probably an intensity rather than a yes/no flag.");
      });
      for (const [name, st] of Object.entries(d.ruleCols)) {
        const share = pct(st.pos, st.n);
        if (st.distinct.size > 3 && share > 50) add(`rule marker "${name}"`, "fail", `${f1(share)}% of cells count as positive, and it has many different values (e.g. ${[...st.distinct].slice(0, 3).join(", ")})`,
          "Rules need yes/no columns (0/1, Positive/Negative). For HALO use the '… Positive' columns, not '… Intensity'.");
      }
    }
    // ---- readouts ----
    for (const [k, st] of Object.entries(d.colStats || {})) {
      const label = k.startsWith("bleed_") ? `bleed marker (${k.slice(6)})` : `readout ${k}`;
      const tot = st.n + st.nan;
      if (st.n === 0) add(label, "warn", `"${st.col}" has no numbers` + (st.textEx.length ? ` (values look like: ${st.textEx.join(", ")})` : ""),
        "Left out of the analysis. Choose a numeric intensity column for it, or set it to not available.");
      else if (st.distinct.size < 2) add(label, "warn", `"${st.col}" has the same value (${st.min}) in every cell`, "Left out of the analysis: a constant readout cannot show differences.");
      else {
        if (pct(st.nan, tot) > 20) add(label, "warn", `"${st.col}" is missing or non-numeric in ${f1(pct(st.nan, tot))}% of cells` + (st.textEx.length ? ` (e.g. ${st.textEx.join(", ")})` : ""), "Those cells are ignored for this readout.");
        if (st.distinct.size === 2) add(label, "warn", `"${st.col}" has only two values (${[...st.distinct].join(" / ")}), so it looks like a yes/no flag`, "It still works, but percentile comparisons are coarse. Prefer an intensity column if you have one.");
      }
    }
    for (const col of Object.keys(cells.cols).filter((k) => k.startsWith("bleed_"))) {
      const r = col.slice(6);
      if (!counts[r] || !counts.T || !usable(cells, col)) continue;
      const p = pctRankByField(cells.cols[col], byImg);
      const inR = [], inT = [];
      for (let i = 0; i < cells.n; i++) { if (Number.isNaN(p[i])) continue; if (cells.type[i] === r) inR.push(p[i]); else if (cells.type[i] === "T") inT.push(p[i]); }
      const gap = median(inR) - median(inT);
      add(`marker sanity: ${r}`, gap > 20 ? "pass" : "fail", `${r}'s own marker sits at the ${median(inR).toFixed(0)}th percentile in ${r} vs ${median(inT).toFixed(0)}th in T cells (gap ${gap.toFixed(0)}, expect >20)`,
        gap > 20 ? "" : `Either the ${r} label/rule catches the wrong cells, or "${(m.bleed || {})[r]}" is not ${r}'s marker. Fix the mapping, or set the bleed marker to not available.`);
    }
    const have = Object.keys(m.states || {}).filter((s) => usable(cells, s));
    add("T-cell readouts", have.length ? "pass" : "warn", `usable: ${have.join(", ") || "none"}`, have.length ? "" : "Without readouts only the proximity tests run.");
    // ---- patients ----
    const byPat = groupIndex(cells.patient), weak = [];
    for (const [p, idx] of byPat) { let t = 0; for (const i of idx) if (cells.type[i] === "T") t++; if (t < 20) weak.push(p); }
    add("patients", byPat.size >= 5 ? "pass" : "warn", `${byPat.size} patients, ${byImg.size} images` + (weak.length ? `; ${weak.length} patients have fewer than 20 T cells (${weak.slice(0, 5).join(", ")})` : ""),
      byPat.size >= 5 ? (weak.length ? "Those patients will show blanks for most tests." : "") : "With fewer than 5 patients, cross-patient comparisons are unreliable. Check the patient rule if you expected more.");
    if (byPat.size > 0.5 * byImg.size && byPat.size > 50) add("patient rule", "warn", `${byPat.size} patients for ${byImg.size} images`, "Almost one patient per image: if images come from fewer patients, change the patient rule.");
    return out;
  }

  const api = { VERSION, ROLES, STATES, fnv1a, mulberry32, choose, createParser, eachCSVRow, csvHeader, detectDelim, prepText, sniff, derive, truthy, FILE_COL,
    createImporter, importCells, importFiles, patientRows, runChecks, medianSpacing, pctRankByField, usable, imageKeys };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Recipe = api;
})(typeof self !== "undefined" ? self : this);
