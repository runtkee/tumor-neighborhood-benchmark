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
  if (failed) { console.log(`\n${failed} failed`); process.exit(1); } else console.log("\nall validation tests passed");
})();
