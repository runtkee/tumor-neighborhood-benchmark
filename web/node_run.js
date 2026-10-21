// Run the browser recipe in Node: node web/node_run.js [--stream] <config.json> <out.json> <file> [file ...]
// --stream reads files through the same chunked path the browser uses (Blob.stream), instead of whole strings.
const fs = require("fs");
const path = require("path");
const R = require("./recipe.js");
const args = process.argv.slice(2);
const stream = args[0] === "--stream" ? args.shift() : null;
const [cfgPath, outPath, ...paths] = args;
const m = JSON.parse(fs.readFileSync(cfgPath, "utf8"));
(async () => {
  const t0 = Date.now();
  let cells;
  try {
    if (stream) {
      const blobs = await Promise.all(paths.map(async (p) => { const b = await fs.openAsBlob(p); b.name = path.basename(p); return Object.assign(b, { name: path.basename(p) }); }));
      const named = blobs.map((b, i) => ({ name: path.basename(paths[i]), size: b.size, slice: (...a) => b.slice(...a), stream: () => b.stream() }));
      cells = await R.importFiles(named, m);
    } else {
      cells = R.importCells(paths.map((p) => ({ name: path.basename(p), text: fs.readFileSync(p, "utf8") })), m);
    }
  } catch (e) {
    fs.writeFileSync(outPath, JSON.stringify({ error: e.message, fix: e.fix || "" }));
    console.log("error:", e.message); return;
  }
  const checks = R.runChecks(cells, m);
  const blocked = checks.some((c) => c[1] === "fail");
  const rows = R.patientRows(cells, m);
  const mem = process.memoryUsage();
  fs.writeFileSync(outPath, JSON.stringify({ checks, rows, n_cells: cells.n, blocked, seconds: (Date.now() - t0) / 1000, heap_mb: mem.heapUsed / 1e6, rss_mb: mem.rss / 1e6 }));
  console.log(`cells ${cells.n}, patients ${rows.length}, ${(Date.now() - t0) / 1000}s, rss ${(mem.rss / 1e6).toFixed(0)} MB`);
})();
