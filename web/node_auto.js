// Full automatic pipeline in Node, same steps as the browser: node web/node_auto.js <out.json> <file> [file ...]
const fs = require("fs"), path = require("path");
const R = require("./recipe.js"), A = require("./automap.js");
const [outPath, ...paths] = process.argv.slice(2);
async function analyse(files) {
  const samples = [];
  for (const f of files) { const s = await A.sampleRows(async (a, b) => f.text.slice(a, b), f.text.length); samples.push({ name: f.name, ...s, size: f.text.length }); }
  const bad = files.map((f) => R.sniff(new TextEncoder().encode(f.text.slice(0, 4096)), f.name)).filter(Boolean);
  if (bad.length) return { fails: bad.map((b) => ({ title: b.title, why: "", fix: b.fix })) };
  const am = A.autoMap(samples);
  if (am.fails.length) return { report: am.report, fails: am.fails };
  const cfg = am.config; cfg.dataset = path.basename(files[0].name).replace(/\.[^.]+$/, "");
  const counts = new Map();
  if (cfg.columns.type) for (const f of files) {
    const t = R.prepText(f.text); let H = null;
    R.eachCSVRow(t, (r) => { if (H === null) { H = r.map((x) => x.trim()).indexOf(cfg.columns.type); return; } counts.set(r[H], (counts.get(r[H]) || 0) + 1); }, R.detectDelim(t));
  }
  const sampleCells = R.importCells(samples.map((s) => ({ name: s.name, text: [s.header].concat(s.rows).map((r) => r.map((v) => (/[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v)).join(",")).join("\n") })), { ...cfg, pixel_size_um: 1, roles: {} });
  const fin = A.finishConfig(cfg, [...counts.entries()].sort((a, b) => b[1] - a[1]), cfg._unitKnown ? null : R.medianSpacing(sampleCells));
  const report = am.report.concat(fin.report);
  if (fin.fails.length) return { report, fails: fin.fails, config: cfg };
  const cells = R.importCells(files, cfg);
  const checks = R.runChecks(cells, cfg);
  const failed = checks.filter((c) => c[1] === "fail");
  if (failed.length) return { report, config: cfg, checks, fails: failed.map((c) => ({ title: c[0], why: c[2], fix: c[3] })) };
  return { report, config: cfg, checks, rows: R.patientRows(cells, cfg), fails: [] };
}
const files = paths.map((p) => ({ name: path.basename(p), text: fs.readFileSync(p, "latin1").length && fs.readFileSync(p, "utf8") }));
(async () => { const t0 = Date.now(); const out = await analyse(files); out.seconds = (Date.now() - t0) / 1000;
fs.writeFileSync(outPath, JSON.stringify(out, (k, v) => (v instanceof Map || v instanceof Set ? undefined : v)));
for (const r of out.report || []) console.log(`  [${r.confidence}] ${r.what}: ${r.decision}  — ${r.evidence}`);
for (const f of out.fails) console.log(`  FAIL ${f.title}: ${f.why} → ${f.fix}`);
if (out.rows) console.log(`  OK ${out.rows.length} patients, ${out.seconds}s`); })();
