// Prints Findings.pooled(rows, key) as JSON (NaN -> null). Used by tests/test_pooled.py.
// Usage: node tests/pooled_parity.js rows.json key
const fs = require("fs");
const F = require("../web/findings.js");
const rows = JSON.parse(fs.readFileSync(process.argv[2], "utf8"));
console.log(JSON.stringify(F.pooled(rows, process.argv[3]), (k, v) => (typeof v === "number" && !Number.isFinite(v) ? null : v)));
