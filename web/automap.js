// Automatic mapping: decide how to read a cell table from its column names AND sample values, with no user input.
// autoMap(samples) -> { config, report: [{what, decision, evidence, confidence}], fails: [{title, why, fix}] }
// samples: [{ name, header: [...], rows: [[...], ...], size }] (first few MB of each file).
// mapLabels(counts) then maps every distinct cell-type label (from a full pass) to a role.
// The mapping is shown to the user for transparency; it is not editable. If a required decision can't be made, it fails.
(function (root) {
  "use strict";
  const R = root.Recipe || (typeof require !== "undefined" ? require("./recipe.js") : null);
  const clean = (s) => String(s ?? "").replace(/^﻿/, "").trim();
  const lc = (s) => clean(s).toLowerCase();

  // ---- marker vocabulary (header text or label text) ----
  const M = {
    CD3: /\bcd3(?![0-9])/i, CD8: /\bcd8(?![0-9])/i, CD4: /\bcd4(?![0-9])/i, CD2: /\bcd2(?![0-9])/i, CD5: /\bcd5(?![0-9])/i, CD7: /\bcd7(?![0-9])/i,
    FOXP3: /foxp-?3/i, CD20: /\bcd20(?![0-9])/i, CD19: /\bcd19(?![0-9])/i, CD79A: /\bcd79a?\b/i, CD138: /\bcd138\b|syndecan/i,
    CD66B: /\bcd66b\b|ceacam8/i, CD15: /\bcd15(?![0-9])/i, MPO: /\bmpo\b|myeloperoxidase/i, LOX1: /\blox-?1\b|\bolr1\b/i,
    CD68: /\bcd68\b/i, CD163: /\bcd163\b/i, CD206: /\bcd206\b|\bmrc1\b/i,
    PANCK: /pan-?ck|pan.?cytokeratin|cytokeratin|\bck(?![a-z0-9])|\bkrt\d*/i, EPCAM: /epcam/i, ECAD: /e-?cad(herin)?\b|\bcdh1\b/i, AMACR: /amacr/i,
    SOX10: /sox-?10/i, MELANA: /melan-?a|mart-?1/i, CDX2: /cdx-?2/i, TTF1: /ttf-?1/i,
    SMA: /(?:^|[^a-z])(?:a|α|alpha)?-?sma\b|smooth muscle actin|acta2/i, PDGFRB: /pdgfr-?b|pdgfrβ/i, FAP: /\bfap\b/i,
    KI67: /ki-?67|mki67/i, PD1: /\bpd-?1\b|pdcd1\b|\bcd279\b/i, GZMB: /granzyme.?b|gzmb|grzb/i, LAG3: /lag-?3|\bcd223\b/i, ICOS: /\bicos\b|\bcd278\b/i,
    CPARP: /cleaved.?parp|\bc-?parp\b|cleaved.?caspase|caspase.?3/i,
  };
  const LINEAGE = { T: ["CD3", "CD8", "CD4", "CD2", "CD5", "CD7", "FOXP3"], B: ["CD20", "CD19", "CD79A", "CD138"], neut: ["CD66B", "MPO", "CD15"],
    mac: ["CD68", "CD163", "CD206"], tumor: ["PANCK", "EPCAM", "ECAD", "AMACR", "SOX10", "MELANA", "CDX2", "TTF1"], stroma: ["SMA", "PDGFRB", "FAP"] };
  const STATE_MARKER = { prolif: "KI67", PD1: "PD1", GZMB: "GZMB", LAG3: "LAG3", ICOS: "ICOS", apoptosis: "CPARP" };
  const BLEED_MARKERS = { neut: ["CD66B", "MPO", "CD15"], mac: ["CD68", "CD163"], stroma: ["SMA", "PDGFRB"] };
  const markersIn = (text) => Object.keys(M).filter((k) => M[k].test(text));

  // ---- label -> role ----
  const WORDS = [
    ["drop", /\b(dirt|artifact|artefact|debris|exclude[d]?|background)\b|^(undefined|unassigned|nan|none|na|-)$/],
    ["neut_subset", /pmn-?mdsc|\bg-?mdsc|lox-?1\s*\+/],
    ["neut", /neutro|granulocyt|\bpmns?\b/],
    ["mac", /macrophag|\btams?\b|kupffer|histiocyt|\bm[12]\b|\bmacs?\b|^mac|mono.?mac/],
    ["B", /\bb[ _-]?cells?\b|\bbcells?\b|b[ -]?lymph|plasma|^b$/],
    ["T", /\bt[ _-]?cells?\b|\btcells?\b|t[- ]?lymph|cytotoxic|\bt ?helper|\bth[0-9]*\b|\btregs?\b|regulatory t|\bctls?\b|^cd[348]\s*t|^cd[48]t$|tprolif|^t$/],
    ["tumor", /tumou?r|cancer|carcinoma|malignan|epitheli|melanoma|neoplas|hepatocyte/],
    ["stroma", /stroma|fibroblast|\bcafs?\b|smooth muscle|myofib|mural|pericyte|\bsmcs?\b|mesench/],
    ["other", /endotheli|vessel|vascul|\bnk\b|natural killer|dendritic|\bdcs?\b|\bpdc\b|mast|monocyte|lymphatic|nerve|adipo|\bimmune\b|\bother\b|unknown/],
  ];
  function labelRole(label) {
    const t = lc(label);
    if (/ \/ | and |&|\bbnt\b/.test(t)) return ["other", "mixed label"];
    for (const [role, re] of WORDS) if (re.test(t)) return [role, "word"];
    // marker phenotype: "CD3+CD8+", "CD68+CD163+", or a bare marker name like "CD8"
    const pos = (t.match(/[a-z0-9α-]+\+/g) || []).map((x) => x.slice(0, -1));
    const toks = pos.length ? pos : (Object.keys(M).some((k) => M[k].test(t)) ? [t] : []);
    const lin = new Set();
    for (const tok of toks) for (const [role, ms] of Object.entries(LINEAGE)) if (ms.some((k) => M[k].test(tok))) lin.add(role);
    if (lin.has("neut") && toks.some((x) => M.LOX1.test(x))) return ["neut_subset", "markers"];
    if (lin.size === 1) return [[...lin][0], "markers"];
    if (lin.size > 1) return ["other", "markers of several cell types"];
    return ["other", "not recognised"];
  }
  function mapLabels(counts) {
    const roles = {}, detail = [];
    for (const [lab, n] of counts) { const [r, how] = labelRole(lab); roles[lab] = r; detail.push([lab, n, r, how]); }
    return { roles, detail };
  }

  // ---- column profiling from sample values ----
  function profile(sample) {
    const H = sample.header.map(clean), rows = sample.rows, n = rows.length;
    return H.map((name, j) => {
      let nonEmpty = 0, numeric = 0, ints = 0, bin = 0, mono = 0, prev = null, mn = Infinity, mx = -Infinity;
      const distinct = new Map();
      for (const r of rows) {
        const v = r[j]; if (v === undefined || v === "") continue;
        nonEmpty++;
        if (distinct.size < 600) distinct.set(v, (distinct.get(v) || 0) + 1);
        const x = Number(String(v).replace(/^\s+|\s+$/g, ""));
        if (Number.isFinite(x)) { numeric++; if (Number.isInteger(x)) ints++; if (x < mn) mn = x; if (x > mx) mx = x; if (prev !== null && x === prev + 1) mono++; prev = x; }
        if (/^(0|1|0\.0|1\.0|true|false|yes|no|pos|neg|positive|negative|\+|-)$/i.test(String(v).trim())) bin++;
      }
      return { j, name, low: name.toLowerCase(), n, nonEmpty, numeric: nonEmpty ? numeric / nonEmpty : 0, integer: numeric ? ints / numeric : 0,
               binary: nonEmpty ? bin / nonEmpty : 0, distinct: distinct.size, values: distinct, min: mn, max: mx, sequential: numeric > 10 ? mono / numeric : 0 };
    });
  }
  const top = (cols) => cols.sort((a, b) => b.score - a.score)[0];

  // known instrument pixel sizes (µm) for snapping an estimate
  const INSTRUMENTS = [[0.12, "CosMx"], [0.18, "CosMx (2021)"], [0.2125, "Xenium"], [0.25, "40x scanner"], [0.325, "Zeiss / 20x"], [0.377, "CODEX (Keyence 20x)"],
    [0.4962, "Vectra Polaris 20x"], [0.5, "20x scanner"], [0.65, "10x scanner"], [1.0, "imaging mass cytometry / already in µm"]];

  function autoMap(samples) {
    const report = [], fails = [];
    const say = (what, decision, evidence, confidence = "high") => report.push({ what, decision, evidence, confidence });
    const fail = (title, why, fix) => fails.push({ title, why, fix });
    // profile all files together (same header), so a column constant within each file but varying across files is seen
    const s0 = { name: samples[0].name, header: samples[0].header, rows: [], files: [] };
    const per = Math.max(500, Math.floor(20000 / samples.length));
    for (const s of samples) {
      if (s.header.map(clean).join("\u0001") !== samples[0].header.map(clean).join("\u0001")) continue;
      for (const r of s.rows.slice(0, per)) { s0.rows.push(r); s0.files.push(s.name); }
    }
    const P = profile(s0);
    const cfg = { dataset: "", modality: "protein", columns: {}, pixel_size_um: 1, roles: {}, states: {}, bleed: {}, extra_roles: {}, min_cells_per_field: 0 };
    const multi = samples.length > 1;

    // 1. positions: numeric, spread out, not a running index; named like coordinates; boxes if min/max pairs exist
    const isPos = (c) => c.numeric > 0.99 && c.distinct > Math.min(50, c.nonEmpty * 0.2) && c.max - c.min > (c.nonEmpty > 100 ? 20 : 0) && c.sequential < 0.5;
    const axisScore = (c, a) => {
      const t = c.low; let s = 0;
      if (new RegExp(`^${a}$|^${a}:${a}$|^pos_?${a}$|^${a}_?(um|µm|px|pixel|position|coord)`).test(t)) s += 4;
      if (new RegExp(`(centroid|center|centre|position|pos|coord)\\b.{0,12}\\b${a}\\b|\\b${a}\\b.{0,3}(centroid|position|coord)|cell ${a} position|center${a}|centroid_?${a}`).test(t)) s += 4;
      if (/\b(min|max)\b|min$|max$|^min|^max/.test(t)) s -= 6;
      return s;
    };
    const boxCol = (a, mm) => P.find((c) => isPos(c) && new RegExp(`^${a}[ _]?${mm}$|^${mm}[ _]?${a}$|bbox.*${a}.*${mm}|${a}${mm}`, "i").test(c.name));
    const pick = (a) => { const cands = P.filter(isPos).map((c) => ({ ...c, score: axisScore(c, a) })).filter((c) => c.score > 0); return cands.length ? top(cands) : null; };
    let cx = pick("x"), cy = pick("y");
    const bx = [boxCol("x", "min"), boxCol("x", "max")], by = [boxCol("y", "min"), boxCol("y", "max")];
    if (cx && cy) {
      cfg.columns.x = cx.name; cfg.columns.y = cy.name;
      say("Cell position", `x = "${cx.name}", y = "${cy.name}"`, `named as coordinates; numeric with ${cx.distinct.toLocaleString()}+ distinct values`);
    } else if (bx.every(Boolean) && by.every(Boolean)) {
      cfg.columns.x = { min: bx[0].name, max: bx[1].name }; cfg.columns.y = { min: by[0].name, max: by[1].name };
      say("Cell position", `centre of the box ${bx[0].name}–${bx[1].name}, ${by[0].name}–${by[1].name}`, "the file stores each cell as a bounding box (HALO style)");
    } else {
      fail("No cell positions found", "No pair of numeric columns is named like x/y coordinates (e.g. X, Centroid X µm, Pos_X, XMin/XMax).",
        "Export the cell (object) table with centroid or bounding-box coordinates.");
    }

    // 2. units / pixel size
    const posNames = [cfg.columns.x, cfg.columns.y].flatMap((v) => (v && typeof v === "object" ? [v.min, v.max] : [v])).filter(Boolean);
    const unitUm = posNames.some((nm) => /µm|μm|\bum\b|micron/i.test(nm));
    const unitPx = posNames.some((nm) => /\bpx\b|pixel/i.test(nm));
    const allNames = P.map((c) => c.low).join(" | ");
    cfg._unitKnown = unitUm;
    if (unitUm) { cfg.pixel_size_um = 1; say("Units", "positions are in µm", `column name says so ("${posNames[0]}")`); }
    else if (/cyc_\d+_ch_\d+/.test(allNames)) { cfg._unitKnown = true; cfg.pixel_size_um = 0.377; say("Units", "0.377 µm per pixel (CODEX)", "marker columns use CODEX processor naming (Cyc_N_ch_N); 0.377 µm is the CODEX Keyence 20x pixel size", "medium"); }
    else if (/cell x position/.test(allNames) && /\(pixels\)/.test(allNames)) { cfg._unitKnown = true; cfg.pixel_size_um = 0.4962; say("Units", "0.4962 µm per pixel (inForm on Vectra Polaris, 20x)", "inForm export that reports positions in pixels", "medium"); }
    else if (/cell x position/.test(allNames) && /\(microns\)/.test(allNames)) { cfg._unitKnown = true; cfg.pixel_size_um = 1; say("Units", "positions are in µm", "inForm export that reports microns"); }

    // 3. patient column: categorical, named like a patient / donor / case
    const patRe = /patient|donor|\bcase\b|subject|participant|individual|specimen|^sample ?id$|sampleid|^pt$|person/;
    const patC = P.filter((c) => patRe.test(c.low) && c.nonEmpty === c.n && c.distinct >= 1 && !(c.numeric > 0.99 && c.distinct > 2000))
      .map((c) => ({ ...c, score: /patient|donor|subject|case/.test(c.low) ? 2 : 1 }));
    const patCol = patC.length ? top(patC) : null;
    // 4. image column: categorical, named like an image / region / core, and nested inside one patient
    const imgRe = /^(image|image ?location|image ?name|file ?name|filename|fov|roi|region|core|spot|slide|tma|tile|section|sample ?name|path|acquisition|scene|field)$|image|fov|\broi\b|core|tma/;
    const nestsIn = (img, pat) => { const m = new Map(); for (const r of s0.rows) { const f = r[img.j]; if (!m.has(f)) m.set(f, new Set()); m.get(f).add(r[pat.j]); } return [...m.values()].every((s) => s.size === 1); };
    const imgC = P.filter((c) => imgRe.test(c.low) && c.nonEmpty === c.n && c.distinct >= 1 && !(c.numeric > 0.99 && c.distinct > 500) && (!patCol || c.j === patCol.j || nestsIn(c, patCol)))
      .map((c) => ({ ...c, score: (/^(image|fov|roi|core|spot)$|image ?(location|name)?$|file ?name|sample ?name/.test(c.low) ? 3 : 1) + (c.distinct > 1 ? 1 : 0) - (c.distinct === 1 && multi ? 2 : 0) }));
    const imgCol = imgC.length ? top(imgC) : null;
    if (imgCol && (imgCol.distinct > 1 || !multi)) { cfg.columns.field = imgCol.name; say("Image", `"${imgCol.name}"`, `${imgCol.distinct} different values in the first part of the file`, imgCol.distinct > 1 ? "high" : "medium"); }
    else if (multi) { cfg.columns.field = R.FILE_COL; say("Image", "each file is one image", "no image / region / core column found", "medium"); }
    else if (patCol) { cfg.columns.field = patCol.name; say("Image", `one image per patient ("${patCol.name}")`, "no image / region / core column; cells of different patients are kept apart", "medium"); }
    else { cfg.columns.field = R.FILE_COL; say("Image", "the whole file is one image", "no image / region / core or patient column found", "low"); }
    if (patCol) { cfg.columns.patient = patCol.name; say("Patient", `"${patCol.name}"`, `${patCol.distinct} patients in the first part of the file`); }
    else {
      const fj = cfg.columns.field === R.FILE_COL ? -1 : P.find((c) => c.name === cfg.columns.field).j;
      const fieldVals = (row, i) => (fj < 0 ? s0.files[i] : row[fj]);
      // several files: file names usually carry the patient / slide, image names the core; so try file names first
      const sources = (multi ? [[R.FILE_COL, samples.map((s) => s.name), "file"]] : [])
        .concat(cfg.columns.field !== R.FILE_COL ? [[cfg.columns.field, [...new Set(s0.rows.map((r, i) => fieldVals(r, i)))], "image"]] : []);
      let done = false;
      for (const [col, names, kind] of sources) {
        for (const sep of ["_", "-", " ", "."]) {
          const g = new Set(names.map((v) => { const k = v.indexOf(sep); return k > 0 ? v.slice(0, k) : v; }));
          if (g.size >= 2 && g.size < names.length) {
            cfg.columns.patient = col; cfg.patient_rule = { mode: "before", sep };
            say("Patient", `text before the first "${sep}" of the ${kind} name`, `gives ${g.size} patients for ${names.length} ${kind}s (e.g. ${[...g].slice(0, 3).join(", ")})`, "medium");
            done = true; break;
          }
        }
        if (done) break;
      }
      if (!done && multi) { cfg.columns.patient = R.FILE_COL; say("Patient", "each file is one patient", "no patient column, and file names have no shared prefix", "medium"); }
      else if (!done) { cfg.columns.patient = "__one__"; say("Patient", "all cells are treated as one patient", "no patient column found and image names carry no patient prefix", "low"); }
    }

    // 5. cell types: a label column (recognised values) beats marker flags; dye-numbered flags can't be used
    const REGION_NAME = /neighbo|niche|community|\bcn\b|cn_|zone|compartment|tissue|category|domain|region|classifier label|analysis/;
    const REGION_VAL = /enriched|\bbulk\b|boundary|border|niche|zone|region|community|interface|front/;
    const labelCands = P.filter((c) => c.numeric < 0.5 && c.distinct >= 2 && c.distinct <= 500).map((c) => {
      const vals = [...c.values.keys()];
      const known = vals.filter((v) => labelRole(v)[0] !== "other").length / vals.length;
      const named = /cell ?type|phenotype|class(ification)?$|^class|annotation|cluster ?name|lineage|population|^label$|celltype|cell_label|^type$/.test(c.low) ? 1 : 0;
      const region = (REGION_NAME.test(c.low) ? 1 : 0) + (vals.filter((v) => REGION_VAL.test(lc(v))).length / vals.length > 0.3 ? 1 : 0);
      return { ...c, known, score: named * 2 + known * 3 - region * 3 };
    }).filter((c) => (c.known >= 0.25 || (c.score >= 2 && c.known > 0)) && (c.distinct >= 3 || c.score >= 2) && [...c.values.keys()].some((v) => String(v).length > 1));
    // value check: do cells given a role actually carry that role's marker? (percentile of the marker, role vs the rest)
    const markerCols = {};
    for (const [role, keys] of Object.entries({ T: ["CD3", "CD8", "CD4"], mac: ["CD68", "CD163"], neut: ["CD66B", "MPO", "CD15"], B: ["CD20", "CD19"], tumor: ["PANCK", "EPCAM", "ECAD"], stroma: ["SMA"] }))
      for (const k of keys) { const c = P.find((c) => M[k].test(c.name) && c.numeric > 0.95 && c.distinct > 10); if (c) { markerCols[role] = c; break; } }
    const agreement = (cand) => {
      const gaps = [];
      for (const [role, mc] of Object.entries(markerCols)) {
        const v = s0.rows.map((r) => Number(r[mc.j])), order = v.map((x, i) => [x, i]).filter((a) => Number.isFinite(a[0])).sort((a, b) => a[0] - b[0]);
        const pct = new Float64Array(v.length).fill(NaN); order.forEach(([, i], k) => (pct[i] = (100 * (k + 1)) / order.length));
        const inR = [], out = [];
        s0.rows.forEach((r, i) => { if (Number.isNaN(pct[i])) return; (labelRole(r[cand.j])[0] === role ? inR : out).push(pct[i]); });
        if (inR.length >= 20 && out.length >= 20) gaps.push(inR.reduce((a, b) => a + b, 0) / inR.length - out.reduce((a, b) => a + b, 0) / out.length);
      }
      return gaps.length ? gaps.reduce((a, b) => a + b, 0) / gaps.length : null;
    };
    if (labelCands.length > 1) for (const c of labelCands) { c.agree = agreement(c); if (c.agree !== null) c.score += c.agree / 10; }
    const flagCands = P.filter((c) => c.binary > 0.98 && c.distinct <= 3);
    const dyeFlags = flagCands.filter((c) => /dye ?\d+|opal ?\d+|channel ?\d+|ch\d+/i.test(c.name) && !markersIn(c.name).length);
    if (labelCands.length) {
      const L = top(labelCands);
      cfg.columns.type = L.name;
      const others = labelCands.filter((c) => c !== L).map((c) => `"${c.name}"`);
      say("Cell types", `labels in "${L.name}"`, `${Math.round(L.known * 100)}% of its ${L.distinct} values are recognised cell types` +
        (L.agree != null ? `; cells labelled with a type carry that type's marker (average +${L.agree.toFixed(0)} percentile points)` : "") +
        (others.length ? `; chosen over ${others.join(", ")}` : ""), L.known > 0.5 ? "high" : "medium");
    } else {
      // rules from yes/no columns named after lineage markers (e.g. "CD8 Positive", "CD3+ (Opal 520) Positive Classification")
      const byRole = {};
      for (const c of flagCands) {
        const ms = markersIn(c.name);
        for (const [role, list] of Object.entries(LINEAGE)) for (const k of list) if (ms.includes(k)) (byRole[role] ||= []).push([k, c.name]);
        if (ms.includes("LOX1")) (byRole.neut_subset ||= []).push(["LOX1", c.name]);
      }
      const order = ["neut_subset", "T", "B", "mac", "neut", "tumor", "stroma"], rules = [];
      for (const role of order) for (const [k, col] of byRole[role] || []) {
        if (role === "neut_subset") { const n = (byRole.neut || [])[0]; if (n) rules.push({ role, positive: [col, n[1]], negative: [], why: `LOX-1+ and ${n[0]}+` }); }
        else rules.push({ role, positive: [col], negative: [], why: `${k}+` });
      }
      if (rules.length && rules.some((r) => r.role === "T")) {
        cfg.columns.type = null; cfg.phenotype_rules = rules.map(({ why, ...r }) => r);
        say("Cell types", `from yes/no marker columns: ${rules.map((r) => `${r.why} → ${r.role}`).join("; ")}`, "no cell-type label column; first matching rule wins", "medium");
      } else if (dyeFlags.length) {
        const inputs = P.find((c) => c.low === "analysis inputs"), dyes = inputs ? [...String([...inputs.values.keys()][0] || "").matchAll(/(Dye \d+):([^;]+);/g)].filter((x) => x[2].trim()).map((x) => `${x[1]} = ${x[2].trim()}`) : [];
        fail("Cell types can't be identified", `The file marks cells positive per dye (${dyeFlags.slice(0, 3).map((c) => `"${c.name}"`).join(", ")}${dyeFlags.length > 3 ? ", …" : ""}) but never says which antibody each dye is` + (dyes.length ? ` (it only records fluorophores: ${dyes.join(", ")})` : "") + ".",
          "Re-export with marker names in the column headers (e.g. 'CD8 Positive'), or add a phenotype / cell-type column in HALO before exporting.");
      } else {
        fail("Cell types can't be identified", "No column holds cell-type labels this app recognises, and no yes/no columns are named after lineage markers (CD3, CD8, CD68, CD66b, …).",
          "Export a phenotype or cell-type column from your analysis software. This app does not assign cell types from raw intensities.");
      }
    }

    // 6. readouts and each cell type's own marker: numeric, varied, named after the marker; prefer mean intensity
    const statScore = (c, key) => {
      let s = 0; const t = c.low;
      if (/mean|average|avg/.test(t)) s += 2; if (/median/.test(t)) s += 1;
      if (/\b(min|max|std|stdev|sd|total|sum|integrated|positive|classification|area|%|percent|count|threshold)\b/.test(t)) s -= 3;
      if (key === "KI67" && /nucle/.test(t)) s += 1;
      if (["PD1", "LAG3", "ICOS", "CD66B", "CD15", "CD68", "CD163"].includes(key) && /membran|cell\b|whole/.test(t)) s += 1;
      if (["GZMB", "MPO"].includes(key) && /cyto/.test(t)) s += 1;
      return s;
    };
    const markerCol = (key) => {
      const re = M[key];
      const cs = P.filter((c) => re.test(c.name) && c.numeric > 0.95 && c.distinct > 2 && c.max > c.min && !(key === "PD1" && /pd-?l1|cd274/i.test(c.name)))
        .map((c) => ({ ...c, score: statScore(c, key) }));
      if (cs.length) return [top(cs).name, "intensity"];
      const fl = flagCands.filter((c) => re.test(c.name));
      return fl.length ? [fl[0].name, "yes/no flag"] : null;
    };
    const found = [];
    for (const [st, key] of Object.entries(STATE_MARKER)) { const c = markerCol(key); if (c) { cfg.states[st] = c[0]; found.push(`${st} = "${c[0]}"${c[1] === "yes/no flag" ? " (yes/no flag, coarse)" : ""}`); } }
    say("T-cell readouts", found.length ? found.join("; ") : "none found", found.length ? "columns named after the marker; mean intensity preferred" : "no Ki-67, PD-1, granzyme B, LAG-3, ICOS or cleaved-PARP columns", found.length ? "high" : "medium");
    const bl = [];
    for (const [role, keys] of Object.entries(BLEED_MARKERS)) for (const k of keys) { const c = markerCol(k); if (c && c[1] === "intensity") { cfg.bleed[role] = c[0]; bl.push(`${role} = "${c[0]}"`); break; } }
    say("Each cell type's own marker (bleed check)", bl.length ? bl.join("; ") : "none found", "used to measure signal bleeding into touching T cells", bl.length ? "high" : "medium");
    return { config: cfg, report, fails, profile: P };
  }

  // after labels are known (full pass), finish the config and estimate pixel size if units are unknown
  function finishConfig(cfg, labelCounts, spacingAt1) {
    const report = [], fails = [];
    const say = (what, decision, evidence, confidence = "high") => report.push({ what, decision, evidence, confidence });
    if (cfg.columns.type) {
      const { roles, detail } = mapLabels(labelCounts);
      cfg.roles = roles;
      const total = detail.reduce((a, d) => a + d[1], 0), by = {};
      for (const [, n, r] of detail) by[r] = (by[r] || 0) + n;
      say("Cell-type labels", `${detail.length} labels → ${Object.entries(by).sort((a, b) => b[1] - a[1]).map(([r, n]) => `${r} ${Math.round((100 * n) / total)}%`).join(", ")}`, "matched by cell-type words and marker names (e.g. CD3+CD8+ → T)");
      cfg._labels = detail;
      if (roles && Object.values(roles).includes("neut_subset")) cfg.extra_roles = { neut_subset: "neut" };
    } else if (cfg.phenotype_rules && cfg.phenotype_rules.some((r) => r.role === "neut_subset")) cfg.extra_roles = { neut_subset: "neut" };
    if (cfg.extra_roles.neut_subset && cfg.bleed.neut) cfg.bleed.neut_subset = cfg.bleed.neut;
    if (!cfg._unitKnown) {
      if (!(spacingAt1 > 0 && Number.isFinite(spacingAt1))) {
        fails.push({ title: "Cell spacing can't be measured", why: "Positions don't vary enough to find neighbouring cells.", fix: "Check that the x / y columns hold cell positions." });
      } else {
        const est = 8 / spacingAt1;
        const near = INSTRUMENTS.filter(([v]) => Math.abs(v - est) / v <= 0.2);
        cfg._pxEstimated = true;
        if (near.length === 1) { cfg.pixel_size_um = near[0][0]; say("Units", `${near[0][0]} µm per pixel (${near[0][1]}, estimated)`, `the file doesn't state units; cell spacing implies about ${est.toFixed(3)} µm per pixel, and this is the only common value within 20%`, "medium"); }
        else { cfg.pixel_size_um = +est.toFixed(4); say("Units", `about ${est.toFixed(3)} µm per pixel (estimated)`, `the file doesn't state units; estimated from cell spacing (cells about 8 µm apart)` + (near.length ? `. Could be ${near.map(([v, n]) => `${v} (${n})`).join(" or ")}` : "") + ". Results are also recomputed at ±20% to show whether conclusions depend on it.", "low"); }
      }
    }
    return { report, fails };
  }

  // Sample rows from across a whole file (files are often sorted by cell type or image, so the start is biased):
  // 16 evenly spaced chunks; each chunk after the first starts at the next line break.
  const CHUNKS = 16, CHUNK = 256 << 10;
  function sampleRows(readChunk, size, maxRows = 20000) {
    return (async () => {
      const R0 = R, rows = []; let header = null, delim = ",";
      if (size <= CHUNKS * CHUNK) {
        const text = R0.prepText(await readChunk(0, size)); delim = R0.detectDelim(text);
        R0.eachCSVRow(text, (r) => { if (!header) header = r; else if (rows.length < maxRows) rows.push(r); }, delim);
        return { header: header || [], rows, delim, whole: true };
      }
      for (let k = 0; k < CHUNKS; k++) {
        const start = Math.floor((k * (size - CHUNK)) / (CHUNKS - 1));
        let text = await readChunk(start, Math.min(size, start + CHUNK));
        if (k === 0) { text = R0.prepText(text); delim = R0.detectDelim(text); }
        else { const nl = text.indexOf("\n"); text = nl >= 0 ? text.slice(nl + 1) : ""; }
        text = text.slice(0, text.lastIndexOf("\n") + 1);
        let first = k === 0;
        R0.eachCSVRow(text, (r) => { if (first) { header = r; first = false; return; } if (header && r.length === header.length && rows.length < maxRows) rows.push(r); }, delim);
      }
      return { header: header || [], rows, delim, whole: false };
    })();
  }

  const api = { autoMap, finishConfig, mapLabels, labelRole, profile, sampleRows, MARKERS: M };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.AutoMap = api;
})(typeof self !== "undefined" ? self : this);
