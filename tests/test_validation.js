// Validation tests for web/recipe.js: parser chunking, streaming == whole string, and bad-data diagnostics.
// Run: node tests/test_validation.js   (exit 1 on failure)
const R = require("../web/recipe.js");
const assert = require("assert");
let failed = 0;
const test = (name, fn) => Promise.resolve().then(fn).then(() => console.log("ok  ", name), (e) => { failed++; console.log("FAIL", name, "\n   ", e.message); });

// synthetic tidy table: 6 patients x 2 images, ~8 µm spacing
function table({ patients = 6, perImage = 900, extraCols = {}, rowHook } = {}) {
  const head = ["patient", "image", "x", "y", "type", "ki67", "cd68", "mpo"].concat(Object.keys(extraCols));
  const types = ["T", "T", "neut", "mac", "B", "tumor", "tumor", "stroma"];
  const lines = [head.join(",")];
  let seed = 7; const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let p = 0; p < patients; p++) for (let im = 0; im < 2; im++) for (let i = 0; i < perImage; i++) {
    const t = types[i % types.length];
    let row = ["P" + p, "img" + im, (i % 30) * 8 + rnd() * 2, Math.floor(i / 30) * 8 + rnd() * 2, t, (rnd() * 100).toFixed(3),
      (t === "mac" ? 80 + rnd() * 20 : rnd() * 40).toFixed(2), (t === "neut" ? 80 + rnd() * 20 : rnd() * 40).toFixed(2)]
      .concat(Object.values(extraCols).map((f) => f(i, t)));
    if (rowHook) row = rowHook(row, p, im, i);
    if (row) lines.push(row.join(","));
  }
  return lines.join("\n") + "\n";
}
const CFG = { dataset: "T", pixel_size_um: 1, columns: { patient: "patient", field: "image", x: "x", y: "y", type: "type" },
  roles: { T: "T", neut: "neut", mac: "mac", B: "B", tumor: "tumor", stroma: "stroma" },
  states: { prolif: "ki67" }, bleed: { mac: "cd68", neut: "mpo" } };
const status = (checks, name) => { const c = checks.find((x) => x[0].startsWith(name)); return c ? c[1] : "absent"; };
const run = (text, cfg = CFG) => { const cells = R.importCells(text, cfg); return { cells, checks: R.runChecks(cells, cfg) }; };

(async () => {
  await test("clean data passes every check", () => {
    const { checks } = run(table());
    const bad = checks.filter((c) => c[1] === "fail"); assert.deepStrictEqual(bad, []);
  });
  await test("parser: random chunk sizes give identical rows (quotes, \"\" escapes, CRLF)", () => {
    const text = 'a,b,c\r\n"x,1","he said ""hi""",3\r\n4,"multi\nline",6\r\n7,8,"9"\r\n';
    const whole = []; R.eachCSVRow(text, (r) => whole.push(r.slice()), ",");
    for (let trial = 0; trial < 300; trial++) {
      const got = []; const p = R.createParser(",", (r) => got.push(r.slice()));
      let i = 0; while (i < text.length) { const k = 1 + Math.floor(Math.random() * 5); p.feed(text.slice(i, i + k)); i += k; }
      p.end(); assert.deepStrictEqual(got, whole);
    }
  });
  await test("streaming import == whole-string import", async () => {
    const text = "﻿" + table().replace(/\n/g, "\r\n");
    const a = R.importCells(text, CFG);
    const blob = new Blob([text]);
    const b = await R.importFiles([{ name: "data.csv", size: blob.size, slice: (...x) => blob.slice(...x), stream: () => blob.stream() }], CFG);
    assert.strictEqual(a.n, b.n); assert.deepStrictEqual(Array.from(a.x), Array.from(b.x));
    assert.deepStrictEqual(R.patientRows(a, CFG), R.patientRows(b, CFG));
  });
  await test("Excel / HDF5 / gzip / PDF files are recognised", () => {
    assert.match(R.sniff(new Uint8Array([0x50, 0x4b, 3, 4]), "a.xlsx").title, /Excel/);
    assert.match(R.sniff(new Uint8Array([0x89, 72, 68, 70]), "a.h5ad").title, /HDF5/);
    assert.match(R.sniff(new Uint8Array([0x1f, 0x8b, 8, 0]), "a.gz").title, /gzip/);
    assert.match(R.sniff(new Uint8Array([37, 80, 68, 70]), "a.pdf").title, /PDF/);
    assert.strictEqual(R.sniff(new TextEncoder().encode("a,b\n1,2\n"), "a.csv"), null);
  });
  await test("wrong pixel size fails the spacing check and suggests a value", () => {
    const { checks } = run(table(), { ...CFG, pixel_size_um: 0.25 });
    assert.strictEqual(status(checks, "cell spacing"), "fail");
    assert.match(checks.find((c) => c[0] === "cell spacing")[3], /pixel size/);
  });
  await test("cell spacing: 4.5 µm passes, 3.5 µm fails (lower bound 4 µm)", () => {
    const text = table(); const med1 = R.medianSpacing(R.importCells(text, { ...CFG, pixel_size_um: 1 }));
    for (const [target, want] of [[4.5, "pass"], [3.5, "fail"]]) {
      const cfg = { ...CFG, pixel_size_um: target / med1 };
      const c = R.runChecks(R.importCells(text, cfg), cfg).find((x) => x[0] === "cell spacing");
      assert.strictEqual(c[1], want, `${target} µm: ${c[2]}`); assert.match(c[2], /4–15/);
    }
  });
  await test("truncated rows: skipped, counted, and fail above 5%", () => {
    const { checks } = run(table({ rowHook: (r, p, im, i) => (i % 10 === 0 ? r.slice(0, 4) : r) }));
    assert.strictEqual(status(checks, "incomplete rows"), "fail");
  });
  await test("text in coordinates: cells dropped, explained", () => {
    const { cells, checks } = run(table({ rowHook: (r, p, im, i) => (i % 50 === 0 ? (r[2] = "n/a", r) : r) }));
    assert.strictEqual(status(checks, "missing coordinates"), "warn"); assert.ok(cells.n < 6 * 2 * 900);
  });
  await test("same file loaded twice is flagged as duplicates", () => {
    const t = table(); const cells = R.importCells([{ name: "a.csv", text: t }, { name: "b.csv", text: t }], CFG);
    assert.strictEqual(status(R.runChecks(cells, CFG), "duplicate cells"), "fail");
  });
  await test("constant readout is reported and left out of the analysis", () => {
    const cfg = { ...CFG, states: { prolif: "ki67", PD1: "flat" } };
    const { cells, checks } = run(table({ extraCols: { flat: () => "5" } }), cfg);
    assert.strictEqual(status(checks, "readout PD1"), "warn");
    assert.ok(!Object.keys(R.patientRows(cells, cfg)[0]).some((k) => k.includes("PD1")));
  });
  await test("text readout ('Positive') is reported with examples", () => {
    const cfg = { ...CFG, states: { prolif: "ki67", PD1: "flag" } };
    const { checks } = run(table({ extraCols: { flag: (i) => (i % 3 ? "Negative" : "Positive") } }), cfg);
    const c = checks.find((x) => x[0] === "readout PD1"); assert.match(c[2], /Positive/);
  });
  await test("intensity column used in a marker rule fails", () => {
    const cfg = { ...CFG, columns: { ...CFG.columns, type: null }, phenotype_rules: [{ role: "T", positive: ["ki67"] }] };
    const { checks } = run(table(), cfg);
    assert.ok(checks.some((c) => c[0].startsWith("rule") && c[1] === "fail"));
  });
  await test("mislabelled role fails the marker sanity check", () => {
    const cfg = { ...CFG, roles: { ...CFG.roles, mac: "B", B: "mac" } };  // swapped
    assert.strictEqual(status(run(table(), cfg).checks, "marker sanity: mac"), "fail");
  });
  await test("missing column in a second file stops with the file name and a fix", () => {
    const t2 = table().replace("ki67", "Ki-67 %");
    assert.throws(() => R.importCells([{ name: "a.csv", text: table() }, { name: "slide2.csv", text: t2 }], CFG),
      (e) => /slide2\.csv/.test(e.message) && /ki67/.test(e.message) && !!e.fix);
  });
  await test("no T cells fails; everything excluded fails with explanation", () => {
    assert.strictEqual(status(run(table(), { ...CFG, roles: { ...CFG.roles, T: "other" } }).checks, "cells: T cells"), "fail");
    const all = Object.fromEntries(Object.keys(CFG.roles).map((k) => [k, "drop"]));
    assert.strictEqual(status(run(table(), { ...CFG, roles: all }).checks, "cells after mapping"), "fail");
  });
  await test("same image name in two patients stays two images", () => {
    const { cells, checks } = run(table());  // every patient has img0 and img1
    assert.strictEqual(new Set(R.imageKeys(cells)).size, 12);
    assert.strictEqual(status(checks, "image names repeat"), "pass");
  });
  await test("too few patients warns", () => {
    assert.strictEqual(status(run(table({ patients: 3 })).checks, "patients"), "warn");
  });
  await test("v3p checks: tumor cells, nest size, T subsets", () => {
    const { checks } = run(table());
    assert.strictEqual(status(checks, "tumor cells"), "pass"); assert.ok(["pass", "warn"].includes(status(checks, "nest size")));
    assert.strictEqual(status(checks, "T subsets"), "warn");  // plain T labels only
    const sub = run(table(), { ...CFG, roles: { ...CFG.roles, T: "T:CD8" } }).checks;
    assert.strictEqual(status(sub, "T subsets"), "pass");
    assert.strictEqual(status(run(table(), { ...CFG, roles: { ...CFG.roles, tumor: "other" } }).checks, "tumor cells"), "warn");
  });
  await test("v3p columns exist with tumor cells and vanish without", () => {
    const cfg = { ...CFG, roles: { ...CFG.roles, T: "T:CD8" } };
    const keys = (c) => Object.keys(R.patientRows(R.importCells(table(), c), c)[0]);
    assert.ok(keys(cfg).includes("comp_nest_pct_CD8"));
    assert.ok(!keys({ ...cfg, roles: { ...cfg.roles, tumor: "other" } }).some((k) => /^(comp_|nest_|grad_)/.test(k)));
  });
  await test("v3p.1 columns: sums on nest and macrophage-state tests only, nest macrophages only with a nest", () => {
    const cfg = { ...CFG, mac_states: { PDL1: "pdl1" } };
    // left half of each image is tumor (80%) with every fifth cell kept as it was, so there is a nest with immune cells in it
    const text = table({ extraCols: { pdl1: (i, t) => (t === "mac" ? 50 + (i % 7) : 10 + (i % 5)) },
      rowHook: (row, p, im, i) => { if (i % 30 < 15 && i % 5 !== 0) row[4] = "tumor"; return row; } });
    const rows = R.patientRows(R.importCells(text, cfg), cfg), r = rows[0], keys = Object.keys(r);
    assert.strictEqual(R.VERSION, "v3p.1"); assert.strictEqual(r.recipe, "v3p.1");
    for (const base of ["nest_neutTouch_prolif", "nest_macTouch_bleed_mac", "TmacTouch_PDL1", "nest_TmacTouch_PDL1"]) {
      const i = keys.indexOf(base + "_n");
      assert.ok(i >= 0, base); assert.deepStrictEqual(keys.slice(i, i + 5), ["_n", "_near_sum", "_far_n", "_far_sum", "_diff"].map((x) => base + x));
    }
    assert.ok(!keys.some((k) => /^(neut|mac|stroma)Touch_.*_(near_sum|far_n|far_sum)$/.test(k)));   // legacy whole-image tests have no sums
    assert.ok(rows.every((x) => x.TmacTouch_PDL1_n >= x.nest_TmacTouch_PDL1_n));
    const noTumor = { ...cfg, roles: { ...cfg.roles, tumor: "other" } };
    const k2 = Object.keys(R.patientRows(R.importCells(text, noTumor), noTumor)[0]);
    assert.ok(k2.includes("TmacTouch_PDL1_near_sum") && !k2.some((k) => k.startsWith("nest_")));
  });
  await test("pooled(): closed-form estimate, bootstrap from the shared generator, skips patients without data", () => {
    const F = require("../web/findings.js"), key = "nest_neutTouch_PD1", rng0 = R.mulberry32(R.fnv1a("pooled|" + key));
    const rows = []; let seed = 11; const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
    for (let i = 0; i < 40; i++) { const n = i % 4 ? Math.floor(rnd() * 8) : 0, f = Math.floor(rnd() * 300);
      rows.push({ patient: "P" + i, [key + "_n"]: n, [key + "_near_sum"]: n * (40 + 10 * rnd()), [key + "_far_n"]: f, [key + "_far_sum"]: f * (45 + 10 * rnd()) }); }
    rows.push({ patient: "bad", [key + "_n"]: 3, [key + "_near_sum"]: NaN, [key + "_far_n"]: 5, [key + "_far_sum"]: 1 });
    const use = rows.filter((r) => Number.isFinite(r[key + "_near_sum"]) && (r[key + "_n"] > 0 || r[key + "_far_n"] > 0));
    const S = (xs, k) => xs.reduce((a, r) => a + r[key + k], 0);
    const p = F.pooled(rows, key);
    assert.ok(Math.abs(p.diff - (S(use, "_near_sum") / S(use, "_n") - S(use, "_far_sum") / S(use, "_far_n"))) < 1e-12);
    assert.strictEqual(p.n_near, S(use, "_n")); assert.strictEqual(p.used, use.length); assert.strictEqual(p.patients, use.filter((r) => r[key + "_n"] > 0).length);
    const est = []; for (let b = 0; b < 200; b++) { let a = 0, bs = 0, f = 0, g = 0;
      for (let j = 0; j < use.length; j++) { const r = use[Math.floor(rng0() * use.length)]; a += r[key + "_n"]; bs += r[key + "_near_sum"]; f += r[key + "_far_n"]; g += r[key + "_far_sum"]; }
      if (a > 0 && f > 0) est.push(bs / a - g / f); }
    est.sort((x, y) => x - y); const q = (f) => { const pos = (est.length - 1) * f, lo = Math.floor(pos); return est[lo] + (est[Math.min(lo + 1, est.length - 1)] - est[lo]) * (pos - lo); };
    assert.ok(Math.abs(p.lo - q(0.025)) < 1e-12 && Math.abs(p.hi - q(0.975)) < 1e-12);
    assert.strictEqual(F.pooled(rows, "missing_test"), null);
  });
  await test("findings: pooled sentence when too few patients qualify, 'too rare' under 50 cells, nothing without sums", () => {
    const F = require("../web/findings.js"), key = "nest_neutTouch_PD1";
    const mk = (nNear, withSums = true) => Array.from({ length: 30 }, (_, i) => {
      const n = i < 20 ? nNear : 0, r = { patient: "P" + i, T_near_neut_log2: -0.2, [key + "_n"]: n, [key + "_diff"]: NaN };
      if (withSums) Object.assign(r, { [key + "_near_sum"]: n * 40, [key + "_far_n"]: 100, [key + "_far_sum"]: 100 * 50 });
      return r; });
    const text = (rows) => F.summarize(rows).map((f) => f.text).join("\n");
    const t = text(mk(5));      // 100 near cells in 20 patients, near arm 10 pct-pts below the far arm
    assert.match(t, /Inside tumor nests, T cells touching neutrophils are rare \(100 cells in 20 patients\); pooled across patients, PD-1 10\.0 pct-pts lower \(95% CI \u221210\.0 to \u221210\.0\)\. Pooled estimate, not per patient\./);
    assert.match(text(mk(2)), /T cells touching neutrophils are too rare to compare \(40 cells\)\./);
    assert.ok(!/pooled|too rare/.test(text(mk(5, false))));
  });
  // pooled sentences name only the markers whose CI excludes 0 (here each patient has the same values, so the CI is the estimate itself)
  const MARKERS = ["prolif", "PD1", "GZMB", "LAG3", "ICOS", "apoptosis"];
  const poolRows = (spec) => Array.from({ length: 20 }, (_, i) => {
    const r = { patient: "P" + i, T_near_neut_log2: -0.2 };
    for (const [grp, eff] of Object.entries(spec)) for (const m of MARKERS) {
      const base = "nest_neutTouch_" + (grp === "all" ? "" : grp + "_") + m, d = eff[m] || 0;
      Object.assign(r, { [base + "_n"]: 5, [base + "_near_sum"]: 5 * (50 + d), [base + "_far_n"]: 100, [base + "_far_sum"]: 100 * 50, [base + "_diff"]: NaN });
    }
    return r; });
  const poolTexts = (spec) => require("../web/findings.js").summarize(poolRows(spec)).filter((f) => /TouchPooled/.test(f.key));
  await test("findings: pooled sentence names only markers with CI excluding 0 and counts the rest", () => {
    const [t] = poolTexts({ all: { GZMB: 7, PD1: -3 } });
    assert.strictEqual(t.status, "change");
    assert.match(t.text, /pooled across patients, PD-1 3\.0 pct-pts lower \(95% CI −3\.0 to −3\.0\); granzyme B 7\.0 pct-pts higher \(95% CI 7\.0 to 7\.0\); 4 other markers include 0\. Pooled estimate, not per patient\.$/);
    assert.ok(!/LAG-3|ICOS|dividing/.test(t.text));
    assert.match(poolTexts({ all: { GZMB: 7 } })[0].text, /; 5 other markers include 0\./);
    assert.match(poolTexts({ all: { GZMB: 7, PD1: 2, LAG3: 2, ICOS: 2, prolif: 2, apoptosis: 2 } })[0].text, /apoptosis|dying/);   // all six named, no "other" clause
    assert.ok(!/other marker/.test(poolTexts({ all: { GZMB: 7, PD1: 2, LAG3: 2, ICOS: 2, prolif: 2, apoptosis: 2 } })[0].text));
    const none = poolTexts({ all: {} })[0];
    assert.strictEqual(none.status, "none"); assert.match(none.text, /pooled across patients, no marker differs from 0 \(6 tested\)\. Pooled estimate, not per patient\./);
    assert.ok(!/includes 0|CI/.test(none.text));
  });
  await test("findings: subset sentences only when they differ from the all-T sentence; otherwise one folded line, in order", () => {
    const all = { GZMB: 7, PD1: -3 };
    const f = poolTexts({ all, CD8: all, CD4: { GZMB: -5 }, Treg: {} });
    assert.deepStrictEqual(f.map((x) => x.key), ["nest_neutTouchPooled", "nest_neutTouchPooled_CD4", "nest_neutTouchPooled_subsets"]);
    assert.match(f[1].text, /^Inside tumor nests, CD4 T cells touching neutrophils are rare .*granzyme B 5\.0 pct-pts lower/);
    assert.match(f[2].text, /^Inside tumor nests, CD8 and Treg subsets touching neutrophils show the same picture\.$/);
    const same = poolTexts({ all, CD8: all, CD4: all, Treg: all });   // nothing differs: one folded line, no subset sentence
    assert.deepStrictEqual(same.map((x) => x.key), ["nest_neutTouchPooled", "nest_neutTouchPooled_subsets"]);
    assert.match(same[1].text, /CD8, CD4 and Treg subsets touching neutrophils show the same picture\./);
    assert.deepStrictEqual(poolTexts({ all: {}, CD8: {}, CD4: {}, Treg: {} }).map((x) => x.key), ["nest_neutTouchPooled", "nest_neutTouchPooled_subsets"]);
    // all-T names nothing but CD8 does: CD8 gets its sentence; a subset marker in the same direction as all-T is folded
    const g = poolTexts({ all: {}, CD8: { PD1: 4 }, CD4: {} });
    assert.deepStrictEqual(g.map((x) => x.key), ["nest_neutTouchPooled", "nest_neutTouchPooled_CD8", "nest_neutTouchPooled_subsets"]);
    assert.match(g[2].text, /^Inside tumor nests, CD4 subset touching neutrophils shows the same picture\.$/);
    // too rare / too few patients: subsets say nothing, not even a fold
    const rare = poolRows({ all: { GZMB: 7 }, CD8: { GZMB: 7 } }).map((r) => { for (const k of Object.keys(r)) if (/CD8_.*_n$/.test(k) && !/far_n$/.test(k)) r[k] = 1; return r; });
    assert.ok(!require("../web/findings.js").summarize(rare).some((x) => /Pooled_(CD8|subsets)/.test(x.key)));
  });
  await test("findings: under 5 patients each question is answered with the measured values, marked descriptive", () => {
    const F = require("../web/findings.js"), base = "nest_neutTouch_GZMB";
    const r = { patient: "P1", T_near_neut_log2: -0.8, B_near_neut_log2: -0.9, T_near_mac_log2: 0.1, comp_nest_cells: 5000, comp_nest_pct_T: 2,
      [base + "_n"]: 40, [base + "_near_sum"]: 40 * 58, [base + "_far_n"]: 200, [base + "_far_sum"]: 200 * 50, [base + "_diff"]: 8 };
    const f = F.summarize([r]), t = f.map((x) => x.text).join("\n");
    assert.match(t, /^T cells near neutrophils: -0\.80 log2, fewer than chance \(1 patient, descriptive\)\. B cells too \(-0\.90\), so not specific to T cells\.$/m);
    assert.match(t, /^T cells near macrophages: \+0\.10 log2, about chance \(1 patient, descriptive\)\.$/m);
    assert.match(t, /^Inside tumor nests, T cells touching neutrophils: granzyme B higher \(\+8\.0 pct-pts\) \(40 cells, 1 patient, descriptive\)\.$/m);
    assert.match(t, /no cell-death readout/);
    assert.deepStrictEqual(f.slice(0, 2).map((x) => x.status), ["down", "none"]);
    r[base + "_n"] = 3;
    assert.match(F.summarize([r]).map((x) => x.text).join("\n"), /touching neutrophils: only 3 cells, too few to compare\./);
  });
  await test("findings: B-cell control wording has three tiers (same, weaker, different)", () => {
    const F = require("../web/findings.js");
    const text = (tneg, bneg, bSign = -1, tSign = -1) => F.summarize(Array.from({ length: 30 }, (_, i) =>
      ({ patient: "P" + i, T_near_neut_log2: tSign * 0.3, B_near_neut_log2: (i < bneg ? bSign : -bSign) * 0.2 })))[0].text;
    assert.match(text(30, 30), /^Fewer T cells near neutrophils than chance in 30 of 30 patients.* B cells show the same pattern, so it isn't specific to T cells\.$/);
    assert.match(text(30, 21), /same pattern, so it isn't specific/);                       // 70% is consistent
    assert.match(text(30, 20), / B cells show a weaker version of it \(20 of 30\), so specificity is unclear\.$/);
    assert.match(text(30, 17), / B cells show a weaker version of it \(17 of 30\), so specificity is unclear\.$/);   // 17/30 = 57%
    assert.match(text(30, 16), / B cells don't, so it looks specific to T cells\.$/);          // 16/30 = 53%
    assert.match(text(30, 5), / B cells don't, so it looks specific to T cells\.$/);
    assert.match(text(30, 18, 1, 1), /^More T cells near neutrophils than chance.* B cells show a weaker version of it \(18 of 30\)/);   // up direction
    assert.match(text(30, 5, 1, 1), /B cells don't/);
  });
  // ---- kit results tables (share_summary.csv, patients_full.csv, kit/reference/*.csv) ----
  const A = require("../web/automap.js"), fs = require("fs"), path = require("path");
  const KIT_HEAD = "dataset,patient,recipe,T_near_neut_log2,B_near_neut_log2,T_near_mac_log2,neutTouch_PD1_n,neutTouch_PD1_diff,comp_nest_cells,grad_neut_0-5_log2,note";
  const kitRows = (ds, n, recipe = "v3p.1") => Array.from({ length: n }, (_, i) => `${ds},P${i},${recipe},${(-0.3 + i / 50).toFixed(3)},0.1,0.2,${10 + i},${i % 2 ? "" : "2.5"},${100 + i},-0.1,hello`);
  await test("kit results: header detection needs dataset, patient, recipe and a _log2/_diff column", () => {
    const h = (s) => A.isKitHeader(s.split(","));
    assert.ok(h(KIT_HEAD)); assert.ok(h("dataset,patient,recipe,x_diff")); assert.ok(h("Dataset,Patient,Recipe,T_near_mac_log2"));   // header names are matched without case
    assert.ok(!h("dataset,patient,T_near_neut_log2"));              // no recipe
    assert.ok(!h("patient,recipe,T_near_neut_log2"));               // no dataset
    assert.ok(!h("dataset,patient,recipe,comp_nest_cells"));        // no _log2/_diff
    assert.ok(!h("patient,image,x,y,type,ki67"));                   // cell table
    assert.ok(!A.isKitHeader(KIT_HEAD.split(","), A.KIT_MAX_BYTES + 1));   // too big to be a results table
  });
  await test("kit results: a cell table that has a patient column is not misdetected", () => {
    const head = table().split("\n")[0].split(",");
    assert.ok(head.includes("patient") && !A.isKitHeader(head));
    assert.ok(!A.isKitHeader(["dataset", "patient", "image", "x", "y", "type"]));       // dataset + patient, no recipe
    assert.ok(!A.isKitHeader(["dataset", "patient", "recipe", "x", "y", "type"]));      // recipe but no result column
    assert.strictEqual(A.readKitResults(table(), "cells.csv", R.VERSION), null);
    // a lookalike header with more than KIT_MAX_ROWS rows is read as a cell table
    const big = [KIT_HEAD].concat(kitRows("D", A.KIT_MAX_ROWS + 1)).join("\n");
    assert.strictEqual(A.readKitResults(big, "big.csv", R.VERSION), null);
    assert.ok(A.readKitResults([KIT_HEAD].concat(kitRows("D", A.KIT_MAX_ROWS)).join("\n"), "ok.csv", R.VERSION));
  });
  await test("kit results: one entry per dataset value, numbers as the importer would read them, NaN for empty", () => {
    const text = "﻿" + [KIT_HEAD].concat(kitRows("DS_A", 6), kitRows("DS_B", 7), ["DS_B,,v3p.1,1,1,1,1,1,1,1,x"]).join("\r\n") + "\r\n";
    const out = A.readKitResults(text, "share_summary.csv", R.VERSION);
    assert.deepStrictEqual(out.datasets.map((d) => [d.id, d.rows.length]), [["DS_A", 6], ["DS_B", 7]]);   // blank patient skipped
    assert.strictEqual(out.blank, 1);
    const r = out.datasets[0].rows;
    assert.strictEqual(r[0].dataset, "DS_A"); assert.strictEqual(r[0].patient, "P0");
    assert.strictEqual(r[0].T_near_neut_log2, -0.3); assert.strictEqual(r[2].neutTouch_PD1_n, 12);
    assert.ok(Number.isNaN(r[1].neutTouch_PD1_diff) && r[0].neutTouch_PD1_diff === 2.5);   // empty = NaN
    assert.ok("comp_nest_cells" in r[0] && "grad_neut_0-5_log2" in r[0]);
    assert.ok(!("note" in r[0]) && !("recipe" in r[0]));                                     // other columns ignored
    for (const s of ["NA", "nan", "-", "n/a", "  ", "abc"]) {
      const o = A.readKitResults([KIT_HEAD, `D,P1,v3p.1,${s},0,0,0,0,0,0,x`].join("\n"), "x.csv", R.VERSION);
      assert.ok(Number.isNaN(o.datasets[0].rows[0].T_near_neut_log2), s);
    }
  });
  await test("kit results: checks (recipe version, patients, measures, contents) and the read-out report", () => {
    const get = (text) => A.readKitResults(text, "f.csv", R.VERSION).datasets[0];
    const chk = (d, n) => d.checks.find((c) => c[0] === n);
    const same = get([KIT_HEAD].concat(kitRows("D", 6)).join("\n"));
    assert.strictEqual(chk(same, "recipe version")[1], "pass"); assert.strictEqual(chk(same, "patients")[1], "pass");
    assert.match(chk(same, "measures present")[2], /whole image.*tumor nest|whole image/); assert.match(chk(same, "measures present")[2], /Not in this file: .*macrophage markers/);
    assert.match(chk(same, "contents")[2], /per-patient summary numbers only/);
    const other = get([KIT_HEAD].concat(kitRows("D", 3, "v2p.1")).join("\n"));
    assert.strictEqual(chk(other, "recipe version")[1], "warn"); assert.match(chk(other, "recipe version")[2], /v2p\.1.*v3p\.1/); assert.match(chk(other, "recipe version")[3], /may not be comparable/);
    assert.strictEqual(chk(other, "patients")[1], "warn");
    const dup = get([KIT_HEAD].concat(kitRows("D", 6), kitRows("D", 1)).join("\n"));
    assert.match(chk(dup, "patients")[2], /repeats a patient/);
    assert.ok(same.report.some((r) => r.what === "Computation" && /None/.test(r.decision)));
    assert.ok(same.report.some((r) => r.what === "Columns used" && /dataset, patient, recipe/.test(r.decision)));
  });
  await test("kit results: the shipped reference tables are detected and give findings", () => {
    const F = require("../web/findings.js"), dir = path.join(__dirname, "..", "src", "bio", "kit", "reference");
    for (const [f, n] of [["CODEX_CRC_Schurch2020.csv", 35], ["IMMUcan_IMC_panel1.csv", 179]]) {
      const out = A.readKitResults(fs.readFileSync(path.join(dir, f), "utf8"), f, R.VERSION);
      assert.strictEqual(out.datasets.length, 1); assert.strictEqual(out.datasets[0].rows.length, n);
      assert.ok(F.summarize(out.datasets[0].rows, { loaded: true, results: true }).some((x) => x.status !== "na" && x.key !== "compartments"));
    }
  });
  if (failed) { console.log(`\n${failed} failed`); process.exit(1); } else console.log("\nall validation tests passed");
})();
