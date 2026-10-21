// Run the browser recipe in Node: node web/node_run.js <config.json> <out.json> <file> [file ...]
const fs = require("fs");
const path = require("path");
const R = require("./recipe.js");
const [cfgPath, outPath, ...paths] = process.argv.slice(2);
const m = JSON.parse(fs.readFileSync(cfgPath, "utf8"));
const t0 = Date.now();
const files = paths.map((p) => ({ name: path.basename(p), text: fs.readFileSync(p, "utf8") }));
const cells = R.importCells(files, m);
const checks = R.runChecks(cells, m);
const rows = R.patientRows(cells, m);
fs.writeFileSync(outPath, JSON.stringify({ checks, rows, n_cells: cells.n, seconds: (Date.now() - t0) / 1000 }));
console.log(`cells ${cells.n}, patients ${rows.length}, ${(Date.now() - t0) / 1000}s`);
