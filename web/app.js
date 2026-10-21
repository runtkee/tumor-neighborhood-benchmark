// Tumor Neighborhood Benchmark: UI. All computation runs in a Web Worker built from recipe.js; nothing is uploaded.
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const REF = JSON.parse(document.getElementById("ref-data").textContent);
  const RECIPE_SRC = document.getElementById("recipe-src").textContent;
  const EXAMPLE = document.getElementById("example-data") ? document.getElementById("example-data").textContent : "";
  const EXAMPLE_CFG = JSON.parse(document.getElementById("example-cfg").textContent);
  const FILE_COL = "__file__", FILE_LABEL = "(file name)";

  const ROLE_OPTS = [["T", "T cell"], ["neut", "Neutrophil"], ["neut_subset", "Private-only neutrophil subset"], ["mac", "Macrophage"],
    ["B", "B / plasma cell"], ["stroma", "Stroma / fibroblast / smooth muscle"], ["tumor", "Tumor / epithelium"], ["other", "Other"], ["drop", "Exclude (artifact, dirt)"]];
  const STATE_SLOTS = [["prolif", "Dividing (Ki-67)"], ["PD1", "PD-1"], ["GZMB", "Granzyme B"], ["LAG3", "LAG-3"], ["ICOS", "ICOS"], ["apoptosis", "Dying (cleaved PARP)"]];
  const BLEED_SLOTS = [["neut", "Neutrophil marker (MPO, CD15, CD66b)"], ["mac", "Macrophage marker (CD68)"], ["stroma", "Stroma marker (SMA)"]];
  const RULE_MODES = [["whole", "the whole value"], ["before_", 'text before the first "_"'], ["before-", 'text before the first "-"'],
    ["before ", "text before the first space"], ["regex", "a pattern (regular expression)"]];

  // same patterns as src/bio/kit/suggest.py, plus HALO / QuPath / inForm column names
  const COL_HINTS = { patient: /^(patient|patients|patient_id|sampleid|sample name|donor|case)/, field: /^(image location|image|file ?name|fov|roi|region|sample_id|spot|core|path)/,
    x: /^(x|pos_x|centroid_?x|x:x|cell x position|centroid x µm|centroid x um|centerx)$/, y: /^(y|pos_y|centroid_?y|y:y|cell y position|centroid y µm|centroid y um|centery)$/,
    type: /(cluster ?name|cell_?type|celltypes|phenotype|annotation|^class$)/ };
  const ROLE_HINTS = [["drop", /dirt|artifact|undefined|unassigned|^nan$/], ["B", /\bb[ _-]?cell|^b$|plasma|cd20/],
    ["T", /t[ _-]?cell|^cd4|^cd8|treg|cd3\+|^t$|tprolif|cd4t|cd8t/], ["neut_subset", /mdsc|lox-?1/], ["neut", /neutro|granulo|pmn/],
    ["mac", /macro|mac|kupffer|cd68|cd163/], ["stroma", /stroma|fibro|mural|caf|smooth/], ["tumor", /tumou?r|epithel|cancer|carcinoma|hepatocyte|amacr/]];
  const STATE_HINTS = { prolif: /ki-?67|mki67/, PD1: /pd-?1\b|pdcd1/, GZMB: /granzyme|gzmb|grzb/, LAG3: /lag-?3/, ICOS: /\bicos/, apoptosis: /cleaved.?parp|caspase/ };
  const BLEED_HINTS = { neut: /mpo|cd15|cd66b|ceacam8/, mac: /cd68/, stroma: /\bsma\b|asma|pdgfrb/ };
  const first = (cols, re) => cols.find((c) => re.test(String(c).trim().toLowerCase())) || "";

  const S = { header: [], labels: [], checks: null, rows: null, example: false, refOn: new Set(REF.datasets.map((d) => d.id)), forced: false, files: [] };

  // ---- worker ----
  const workerSrc = RECIPE_SRC + "\n;(" + function () {
    let files = [], cells = null, cfgKey = "";
    self.onmessage = async (e) => {
      const { id, cmd } = e.data;
      try {
        let out;
        if (cmd === "load") {
          files = e.data.files ? await Promise.all(e.data.files.map(async (f) => ({ name: f.name, text: await f.text() }))) : [{ name: e.data.name, text: e.data.text }];
          cells = null;
          const headers = files.map((f) => Recipe.csvHeader(f.text));
          const sample = []; const t0 = Recipe.prepText(files[0].text);
          let k = 0; Recipe.eachCSVRow(t0.slice(0, 3e6), (r) => { if (k++ > 0 && sample.length < 4000) sample.push(r); }, Recipe.detectDelim(t0));
          let rows = 0;
          for (const f of files) { for (let i = f.text.indexOf("\n"); i !== -1; i = f.text.indexOf("\n", i + 1)) rows++; rows--; if (!f.text.endsWith("\n")) rows++; }
          const mismatch = headers.slice(1).filter((h) => h.join("\u0001") !== headers[0].join("\u0001")).length;
          out = { header: headers[0], sample, rows, files: files.map((f) => f.name), delim: Recipe.detectDelim(t0), mismatch };
        } else if (cmd === "labels") {
          const counts = new Map();
          for (const f of files) {
            const t = Recipe.prepText(f.text); let H = null;
            Recipe.eachCSVRow(t, (r) => { if (H === null) { H = r.map((x) => x.trim()).indexOf(e.data.col); return; } if (H < 0) return; counts.set(r[H], (counts.get(r[H]) || 0) + 1); }, Recipe.detectDelim(t));
          }
          out = [...counts.entries()].sort((a, b) => b[1] - a[1]);
        } else if (cmd === "preview") {  // patient / image IDs and per-role counts under the current config
          const c = Recipe.importCells(files, e.data.cfg);
          const roleCounts = {}; for (const t of c.type) roleCounts[t] = (roleCounts[t] || 0) + 1;
          const pats = [...new Set(c.patient)], fields = [...new Set(c.field)];
          out = { patients: pats.length, patientEx: pats.slice(0, 4), fields: fields.length, fieldEx: fields.slice(0, 3), roleCounts, dropped: c.dropped_fields || 0,
                  spacing: e.data.spacing ? Recipe.medianSpacing(c) : null };
        } else if (cmd === "check") {
          cells = Recipe.importCells(files, e.data.cfg); cfgKey = JSON.stringify(e.data.cfg);
          out = Recipe.runChecks(cells, e.data.cfg);
          if (cells.dropped_fields) out.push(["small images skipped", "warn", `${cells.dropped_fields} images had fewer cells than the minimum and were left out`]);
        } else if (cmd === "run") {
          if (!cells || cfgKey !== JSON.stringify(e.data.cfg)) cells = Recipe.importCells(files, e.data.cfg);
          out = Recipe.patientRows(cells, e.data.cfg);
        }
        self.postMessage({ id, ok: true, out });
      } catch (err) { self.postMessage({ id, ok: false, error: String(err && err.message || err) }); }
    };
  }.toString() + ")();";
  const worker = new Worker(URL.createObjectURL(new Blob([workerSrc], { type: "text/javascript" })));
  let seq = 0; const pending = new Map();
  worker.onmessage = (e) => { const p = pending.get(e.data.id); pending.delete(e.data.id); e.data.ok ? p.res(e.data.out) : p.rej(new Error(e.data.error)); };
  const call = (cmd, data) => new Promise((res, rej) => { const id = ++seq; pending.set(id, { res, rej }); worker.postMessage({ id, cmd, ...data }); });

  // ---- helpers ----
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(v) ? "—" : (v > 0 ? "+" : "") + v.toFixed(d));
  function select(id, options, value, none, attrs = "") {
    const opts = (none ? [["", none === true ? "— not available —" : none]] : []).concat(options.map((o) => (Array.isArray(o) ? o : [o, o === FILE_COL ? FILE_LABEL : o])));
    return `<select id="${id}" ${attrs}>` + opts.map(([v, l]) => `<option value="${esc(v)}"${v === value ? " selected" : ""}>${esc(l)}</option>`).join("") + "</select>";
  }
  function multi(id, options, values) {
    return `<select id="${id}" multiple size="3">` + options.map((o) => `<option value="${esc(o)}"${values.includes(o) ? " selected" : ""}>${esc(o)}</option>`).join("") + "</select>";
  }
  const status = (el, html) => ($(el).innerHTML = html);
  const setStepEnabled = (id, on) => $(id).setAttribute("aria-disabled", on ? "false" : "true");
  const typeMode = () => ($("#mode-rules").checked ? "rules" : "label");

  // ---- step 1: load ----
  async function load(input, isExample) {
    const names = input.files ? input.files.map((f) => f.name) : [input.name];
    status("#load-status", `<span class="spinner"></span><span>Reading ${esc(names.length > 1 ? names.length + " files" : names[0])}…</span>`);
    S.example = isExample; S.checks = S.rows = null; renderChecks(); renderResults();
    const info = await call("load", input);
    S.header = info.header; S.sample = info.sample; S.files = info.files;
    const delim = { ",": "comma-separated", "\t": "tab-separated", ";": "semicolon-separated" }[info.delim];
    status("#load-status", `<span class="chip pass">loaded</span><span class="mono">${esc(names.length > 1 ? names.length + " files" : names[0])}</span>` +
      `<span class="muted">${info.rows.toLocaleString()} rows · ${info.header.length} columns · ${delim}</span>` +
      (info.mismatch ? `<span class="chip warn">${info.mismatch} file(s) have different columns from the first; matching columns are used by name</span>` : "") +
      (isExample ? `<span class="example-note">Example: 6 patients from the public CODEX colon cancer dataset (Schürch et al. 2020), standing in for private data. CODEX is left out of the reference while this example is loaded.</span>` : ""));
    if (isExample) S.refOn.delete("CODEX_CRC_Schurch2020"); else S.refOn = new Set(REF.datasets.map((d) => d.id));
    await buildMapping(isExample ? EXAMPLE_CFG : null);
  }
  $("#pick").addEventListener("click", () => $("#file").click());
  $("#file").addEventListener("change", (e) => { const f = [...e.target.files]; if (f.length) load({ files: f }, false).catch(showError); });
  const drop = $("#drop");
  ["dragenter", "dragover"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add("over"); }));
  ["dragleave", "drop"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
  drop.addEventListener("drop", (e) => { const f = [...e.dataTransfer.files]; if (f.length) load({ files: f }, false).catch(showError); });
  function showError(err) { status("#run-status", `<span class="chip fail">error</span><span>${esc(err.message)}</span>`); }

  // ---- step 2: map ----
  function idBlock(key, label, col, rule) {
    const mode = !rule || rule.mode === "whole" ? "whole" : rule.mode === "before" ? "before" + rule.sep : "regex";
    return `<div class="field"><label for="col-${key}">${label} comes from</label>${select("col-" + key, [FILE_COL].concat(S.header), col, "— choose —")}` +
      `<label for="rule-${key}">using</label>${select("rule-" + key, RULE_MODES, mode, false)}` +
      `<input type="text" id="re-${key}" placeholder="e.g. ^P(\\d+)_  (first group is used)" value="${esc(rule && rule.regex || "")}"${mode === "regex" ? "" : " hidden"}>` +
      `<span class="hint" id="pv-${key}"></span></div>`;
  }
  function xyBlock(axis, spec) {
    const box = spec && typeof spec === "object";
    return `<div class="field"><label for="${axis}-mode">Cell ${axis}</label>${select(axis + "-mode", [["one", "one column (centre)"], ["box", "box: min and max columns"]], box ? "box" : "one", false)}` +
      `<span id="${axis}-one"${box ? " hidden" : ""}>${select("col-" + axis, S.header, box ? "" : spec || "", "— choose —")}</span>` +
      `<span id="${axis}-box" class="row"${box ? "" : " hidden"}>${select(axis + "-min", S.header, box ? spec.min : "", "min")}${select(axis + "-max", S.header, box ? spec.max : "", "max")}</span></div>`;
  }
  function guessConfig() {
    const H = S.header, low = (re) => first(H, re);
    const hasBox = H.includes("XMin") && H.includes("XMax");
    const posCols = H.filter((h) => / positive( classification)?$/i.test(h));
    const typeCol = low(COL_HINTS.type);
    const pat = low(COL_HINTS.patient), fld = low(COL_HINTS.field);
    return {
      columns: { patient: pat || FILE_COL, field: fld || FILE_COL, type: typeCol,
        x: hasBox ? { min: "XMin", max: "XMax" } : low(COL_HINTS.x), y: hasBox ? { min: "YMin", max: "YMax" } : low(COL_HINTS.y) },
      patient_rule: pat ? { mode: "whole" } : { mode: "before", sep: "_" }, field_rule: { mode: "whole" },
      mode: !typeCol && posCols.length ? "rules" : "label",
      phenotype_rules: !typeCol && posCols.length ? [{ role: "T", positive: [], negative: [] }] : [],
      states: Object.fromEntries(STATE_SLOTS.map(([k]) => [k, low(STATE_HINTS[k])])),
      bleed: Object.fromEntries(BLEED_SLOTS.map(([k]) => [k, low(BLEED_HINTS[k])])),
    };
  }
  async function buildMapping(preset) {
    const g = preset ? { ...preset, mode: preset.phenotype_rules ? "rules" : "label" } : guessConfig();
    $("#map-ids").innerHTML = idBlock("patient", "Patient", g.columns.patient, g.patient_rule) + idBlock("field", "Image", g.columns.field, g.field_rule);
    $("#map-xy").innerHTML = xyBlock("x", g.columns.x) + xyBlock("y", g.columns.y);
    $("#col-type-slot").innerHTML = select("col-type", S.header, g.columns.type || "", "— choose —");
    (g.mode === "rules" ? $("#mode-rules") : $("#mode-label")).checked = true;
    S.rules = (g.phenotype_rules || []).map((r) => ({ role: r.role, positive: r.positive || [], negative: r.negative || [] }));
    renderRules(); showTypeMode();
    $("#map-states").innerHTML = STATE_SLOTS.map(([k, label]) =>
      `<div class="field"><label for="st-${k}">${label}</label>${select("st-" + k, S.header, (g.states || {})[k] || "", true)}</div>`).join("") +
      `<div class="field"><label for="st-extra">Private-only T-cell readout (e.g. LMP)</label>${select("st-extra", S.header, g.extra || "", true)}</div>`;
    $("#map-bleed").innerHTML = BLEED_SLOTS.map(([k, label]) =>
      `<div class="field"><label for="bl-${k}">${label}</label>${select("bl-" + k, S.header, (g.bleed || {})[k] || "", true)}</div>`).join("");
    const dyes = dyeNames(); $("#dye-names").textContent = dyes ? "From the file's Analysis Inputs: " + dyes : "";
    if (preset) { $("#px").value = preset.pixel_size_um; $("#dataset-name").value = preset.dataset; $("#min-cells").value = preset.min_cells_per_field ?? 50; }
    setStepEnabled("#step-map", true); setStepEnabled("#step-check", true);
    if (typeMode() === "label") await refreshLabels(preset ? preset.roles : null);
    await refreshPreview(!preset);
    if (preset) $("#px-hint").textContent = "Instrument value for CODEX (0.377 µm).";
    if (preset) await runChecks().then(() => S.checks && !S.checks.some((c) => c[1] === "fail") && runCompare());
  }
  function dyeNames() {  // HALO stores "Dye 2:Opal 620;..." in Analysis Inputs
    const i = S.header.indexOf("Analysis Inputs"); if (i < 0 || !S.sample.length) return "";
    const m = [...String(S.sample[0][i]).matchAll(/(Dye \d+):([^;]+);/g)].filter((x) => x[2].trim()).map((x) => `${x[1]} = ${x[2].trim()}`);
    return m.join(", ");
  }
  function showTypeMode() { $("#type-label").hidden = typeMode() !== "label"; $("#type-rules").hidden = typeMode() !== "rules"; }
  function renderRules() {
    const pos = S.header.filter((h) => /positive|class|pos$/i.test(h));
    const cols = pos.length ? pos : S.header;
    $("#rules tbody").innerHTML = S.rules.map((r, i) =>
      `<tr><td>${select("rr-" + i, ROLE_OPTS.filter(([v]) => v !== "other"), r.role, false)}</td><td>${multi("rp-" + i, cols, r.positive)}</td>` +
      `<td>${multi("rn-" + i, cols, r.negative)}</td><td><button type="button" data-del="${i}" aria-label="Remove rule ${i + 1}">Remove</button></td></tr>`).join("") ||
      `<tr><td colspan="4" class="muted">No rules yet. Add one per cell type.</td></tr>`;
  }
  function readRules() {
    return S.rules.map((_, i) => ({ role: $("#rr-" + i).value,
      positive: [...$("#rp-" + i).selectedOptions].map((o) => o.value), negative: [...$("#rn-" + i).selectedOptions].map((o) => o.value) }));
  }
  $("#btn-add-rule").addEventListener("click", () => { S.rules = readRules(); S.rules.push({ role: "T", positive: [], negative: [] }); renderRules(); });
  $("#rules").addEventListener("click", (e) => { const i = e.target.dataset.del; if (i !== undefined) { S.rules = readRules(); S.rules.splice(+i, 1); renderRules(); } });
  $("#btn-preview-rules").addEventListener("click", () => refreshPreview(false).catch(showError));
  document.querySelectorAll('input[name="typemode"]').forEach((el) => el.addEventListener("change", async () => {
    showTypeMode(); if (typeMode() === "label") await refreshLabels(null); else if (!S.rules.length) { S.rules = [{ role: "T", positive: [], negative: [] }]; renderRules(); }
  }));
  $("#step-map").addEventListener("change", async (e) => {
    const id = e.target.id || "";
    if (/^(x|y)-mode$/.test(id)) { const a = id[0], box = e.target.value === "box"; $("#" + a + "-one").hidden = box; $("#" + a + "-box").hidden = !box; }
    if (/^rule-(patient|field)$/.test(id)) $("#re-" + id.slice(5)).hidden = e.target.value !== "regex";
    if (id === "col-type") await refreshLabels(null);
    if (/^(col|rule|re)-(patient|field)$/.test(id) || /^(x|y)-(mode|min|max)$/.test(id) || /^col-(x|y)$/.test(id) || id === "min-cells")
      await refreshPreview(/^(x|y)-|^col-(x|y|field)$/.test(id)).catch(() => {});
  });
  async function refreshLabels(presetRoles) {
    const col = $("#col-type").value;
    if (!col) { $("#roles tbody").innerHTML = `<tr><td colspan="3" class="muted">Choose the cell type column above.</td></tr>`; S.labels = []; return; }
    S.labels = await call("labels", { col });
    $("#roles tbody").innerHTML = S.labels.map(([lab, n], i) => {
      const g = presetRoles ? presetRoles[lab] || "other" : (ROLE_HINTS.find(([, re]) => re.test(String(lab).toLowerCase())) || ["other"])[0];
      return `<tr><td class="mono">${esc(lab)}</td><td class="num">${n.toLocaleString()}</td><td>${select("role-" + i, ROLE_OPTS, g, false)}</td></tr>`;
    }).join("");
  }
  async function refreshPreview(estimate) {
    const cfg = readConfig(true);
    if (!cfg) return;
    const pv = await call("preview", { cfg: { ...cfg, pixel_size_um: estimate ? 1 : cfg.pixel_size_um }, spacing: estimate });
    $("#pv-patient").textContent = `${pv.patients} patients, e.g. ${pv.patientEx.join(", ")}`;
    $("#pv-field").textContent = `${pv.fields} images, e.g. ${pv.fieldEx.join(", ")}` + (pv.dropped ? ` (${pv.dropped} skipped as too small)` : "");
    if (typeMode() === "rules") $("#rule-counts").textContent = Object.entries(pv.roleCounts).map(([r, n]) => `${r}: ${n.toLocaleString()}`).join(" · ");
    if (estimate && pv.spacing > 0 && Number.isFinite(pv.spacing)) {
      const est = 8 / pv.spacing; $("#px").value = est.toFixed(3);
      $("#px-hint").textContent = `Estimated ${est.toFixed(3)} µm/px, assuming cells sit about 8 µm apart. Replace it with your instrument's value (use 1 if coordinates are already in µm).`;
    }
  }
  function readRule(key) {
    const v = $("#rule-" + key).value;
    if (v === "whole") return { mode: "whole" };
    if (v.startsWith("before")) return { mode: "before", sep: v.slice(6) };
    return { mode: "regex", regex: $("#re-" + key).value };
  }
  function readAxis(a) { return $("#" + a + "-mode").value === "box" ? { min: $("#" + a + "-min").value, max: $("#" + a + "-max").value } : $("#col-" + a).value; }
  function readConfig(partial) {
    const v = (id) => ($(id) ? $(id).value : "");
    const columns = { patient: v("#col-patient"), field: v("#col-field"), x: readAxis("x"), y: readAxis("y"), type: typeMode() === "label" ? v("#col-type") : null };
    const okAxis = (s) => (typeof s === "object" ? s.min && s.max : s);
    if (!columns.patient || !columns.field || !okAxis(columns.x) || !okAxis(columns.y)) return partial ? null : { error: "Choose patient, image, x and y columns." };
    const roles = {}; S.labels.forEach(([lab], i) => { if ($("#role-" + i)) roles[lab] = v("#role-" + i); });
    const phenotype_rules = typeMode() === "rules" ? readRules() : null;
    if (!partial && typeMode() === "label" && !columns.type) return { error: "Choose the cell type column, or switch to marker rules." };
    if (!partial && phenotype_rules && !phenotype_rules.some((r) => r.positive.length)) return { error: "Each rule needs at least one positive marker." };
    const states = {}; STATE_SLOTS.forEach(([k]) => { if (v("#st-" + k)) states[k] = v("#st-" + k); });
    const extra_states = []; if (v("#st-extra")) { states.private_readout = v("#st-extra"); extra_states.push("private_readout"); }
    const bleed = {}; BLEED_SLOTS.forEach(([k]) => { if (v("#bl-" + k)) bleed[k] = v("#bl-" + k); });
    const usesSubset = phenotype_rules ? phenotype_rules.some((r) => r.role === "neut_subset") : Object.values(roles).includes("neut_subset");
    const extra_roles = usesSubset ? { neut_subset: "neut" } : {};
    if (usesSubset && bleed.neut) bleed.neut_subset = bleed.neut;
    return { dataset: v("#dataset-name") || "MY_COHORT", modality: "protein", columns, pixel_size_um: Number(v("#px")),
      patient_rule: readRule("patient"), field_rule: readRule("field"), roles, phenotype_rules, states, extra_states, bleed, extra_roles,
      min_cells_per_field: Number(v("#min-cells")) || 0 };
  }

  // ---- step 3: check ----
  async function runChecks() {
    const cfg = readConfig(false);
    if (cfg.error) { S.checks = [["mapping", "fail", cfg.error]]; renderChecks(); return; }
    status("#run-status", `<span class="spinner"></span><span>Checking…</span>`);
    S.checks = await call("check", { cfg }); status("#run-status", ""); renderChecks();
  }
  function renderChecks() {
    const tb = $("#checks tbody");
    if (!S.checks) { tb.innerHTML = `<tr><td class="muted">Checks appear here: cell spacing (catches a wrong pixel size), cells per role, whether each cell type's own marker is actually high in it, available readouts, number of patients.</td></tr>`; $("#btn-run").disabled = true; return; }
    tb.innerHTML = S.checks.map(([n, s, d]) => `<tr><td><span class="chip ${s}">${s}</span></td><td>${esc(n)}</td><td class="muted">${esc(d)}</td></tr>`).join("");
    $("#btn-run").disabled = S.checks.some((c) => c[1] === "fail") && !$("#force").checked;
  }
  $("#btn-check").addEventListener("click", () => runChecks().catch(showError));
  $("#force").addEventListener("change", renderChecks);

  // ---- step 4: compare ----
  async function runCompare() {
    const cfg = readConfig(false);
    S.forced = S.checks.some((c) => c[1] === "fail");
    status("#run-status", `<span class="spinner"></span><span>Computing distances and shuffles…</span>`);
    const t0 = performance.now();
    S.rows = await call("run", { cfg }); S.cfg = cfg;
    status("#run-status", `<span class="chip pass">done</span><span class="muted">${S.rows.length} patients in ${((performance.now() - t0) / 1000).toFixed(1)} s</span>`);
    setStepEnabled("#step-compare", true); renderResults(); renderSummary();
  }
  $("#btn-run").addEventListener("click", () => runCompare().catch(showError));

  const CLAIMS = [
    { title: "Do T cells keep away from neutrophils?", desc: "Share of T cells within 30 µm of a neutrophil, compared with shuffled labels. Below 0 means fewer T cells nearby than chance.",
      main: ["T_near_neut_log2"], control: ["B_near_neut_log2"], ctrlText: "Control: B cells. If they keep away just as much, the pattern is about lymphocyte zones, not T cells specifically.", unit: "log2" },
    { title: "Do T cells gather near macrophages?", desc: "Share of T cells within 30 µm of a macrophage, compared with shuffled labels. Above 0 means more T cells nearby than chance.",
      main: ["T_near_mac_log2"], control: [], unit: "log2" },
    { title: "Are T cells touching neutrophils different?", desc: "Percentile-point difference in each readout: T cells under 10 µm from a neutrophil minus T cells 30 µm or more away.",
      main: "neutTouch", control: ["neutTouch_bleed_neut_diff"], ctrlText: "Bleed: the neutrophil marker measured inside touching T cells. A large value means contact comparisons are unreliable.", unit: "pct-pts" },
    { title: "Are T cells touching macrophages different?", desc: "Same contact comparison for macrophages.", main: "macTouch", control: ["macTouch_bleed_mac_diff"],
      ctrlText: "Bleed: the macrophage marker inside touching T cells.", unit: "pct-pts" },
    { title: "Control: T cells touching stroma", desc: "Same contact comparison for stromal cells. Effects here flag technical contact artifacts rather than biology.", main: "stromaTouch",
      control: ["stromaTouch_bleed_stroma_diff"], unit: "pct-pts", isControl: true },
    { title: "Private-only: T cells near your neutrophil subset", desc: "Your subset has no public counterpart. All public neutrophils are shown as the closest reference.",
      main: ["T_near_neut_subset_log2"], also: "neut_subsetTouch", control: ["T_near_neut_log2"], ctrlText: "Reference shown: all neutrophils (public), the closest public role.", unit: "log2", privateOnly: true },
  ];
  const STATE_NAME = { prolif: "dividing (Ki-67)", PD1: "PD-1", GZMB: "granzyme B", LAG3: "LAG-3", ICOS: "ICOS", apoptosis: "dying (cleaved PARP)", private_readout: "your private readout", bleed_neut: "bleed", bleed_mac: "bleed", bleed_stroma: "bleed", bleed_neut_subset: "bleed" };
  function metricLabel(m) {
    const t = m.match(/^(\w+?)Touch_(.+)_diff$/); if (t) return STATE_NAME[t[2]] || t[2];
    return { B_near_neut_log2: "B cells near neutrophils", T_near_neut_log2: "T cells near neutrophils", T_near_mac_log2: "T cells near macrophages", T_near_neut_subset_log2: "T cells near your subset" }[m] || m;
  }
  const quant = (a, q) => { const s = [...a].sort((x, y) => x - y); const p = (s.length - 1) * q, lo = Math.floor(p); return s.length ? s[lo] + (s[Math.min(lo + 1, s.length - 1)] - s[lo]) * (p - lo) : NaN; };
  const vals = (rows, m) => rows.map((r) => r[m]).filter((v) => typeof v === "number" && !Number.isNaN(v));
  function verdict(m) {
    const mine = vals(S.rows, m), pub = vals(REF.rows.filter((r) => S.refOn.has(r.dataset)), m);
    if (!mine.length) return ["info", "no data in your file"];
    if (mine.length < 5) return ["warn", `only ${mine.length} patient${mine.length > 1 ? "s" : ""}`];
    if (pub.length < 5) return ["info", "no public counterpart"];
    const med = quant(mine, 0.5), lo = quant(pub, 0.1), hi = quant(pub, 0.9);
    if (med < lo) return ["fail", "below public range"]; if (med > hi) return ["fail", "above public range"];
    return ["pass", "within public range"];
  }
  function niceTicks(lo, hi, n) {
    const step0 = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(step0))), f = step0 / mag;
    const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
    const out = []; for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9; t += step) out.push(+t.toFixed(10));
    return out;
  }
  function strip(m) {
    const refRows = REF.rows.filter((r) => S.refOn.has(r.dataset));
    const groups = REF.datasets.filter((d) => S.refOn.has(d.id)).map((d) => [d.short, vals(refRows.filter((r) => r.dataset === d.id), m)]).filter((g) => g[1].length);
    const mine = vals(S.rows, m);
    if (!groups.length && !mine.length) return "";
    const all = groups.flatMap((g) => g[1]).concat(mine, [0]);
    let lo = quant(all, 0.01), hi = quant(all, 0.99); if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    const cols = groups.concat([["Yours", mine]]);
    const W = 340, H = 170, L = 40, R = 8, T = 10, B = 40, cw = (W - L - R) / cols.length;
    const y = (v) => T + (H - T - B) * (1 - (Math.min(Math.max(v, lo), hi) - lo) / (hi - lo));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(metricLabel(m))}">`;
    for (const t of niceTicks(lo, hi, 4)) svg += `<line x1="${L}" x2="${W - R}" y1="${y(t)}" y2="${y(t)}" stroke="var(--line)" stroke-width="${t === 0 ? 1.4 : 0.8}"/>` +
      `<text x="${L - 5}" y="${y(t) + 3.5}" text-anchor="end" font-size="9.5" fill="var(--muted)" font-family="var(--f-data)">${t}</text>`;
    cols.forEach(([name, v], i) => {
      const cx = L + cw * (i + 0.5), yours = i === cols.length - 1, col = yours ? "var(--accent)" : "var(--ref)";
      v.forEach((val, j) => { const jit = ((j * 0.6180339887) % 1 - 0.5) * cw * 0.5;
        svg += `<circle cx="${(cx + jit).toFixed(1)}" cy="${y(val).toFixed(1)}" r="${yours ? 3 : 2.2}" fill="${col}" fill-opacity="${yours ? 0.9 : 0.55}"/>`; });
      if (v.length) { const md = quant(v, 0.5); svg += `<line x1="${cx - cw * 0.32}" x2="${cx + cw * 0.32}" y1="${y(md)}" y2="${y(md)}" stroke="${yours ? "var(--accent)" : "var(--ref-strong)"}" stroke-width="2.2"/>`; }
      const words = String(name).split(" ");
      svg += `<text x="${cx}" y="${H - B + 14}" text-anchor="middle" font-size="9.5" fill="${yours ? "var(--accent)" : "var(--muted)"}" font-family="var(--f-body)">${esc(words[0])}</text>`;
      if (words.length > 1) svg += `<text x="${cx}" y="${H - B + 26}" text-anchor="middle" font-size="9.5" fill="var(--muted)" font-family="var(--f-body)">${esc(words.slice(1).join(" "))}</text>`;
    });
    return svg + "</svg>";
  }
  function plotBlock(m, unit) {
    const [cls, text] = verdict(m), mine = vals(S.rows, m);
    return `<div class="plot"><div class="row"><b>${esc(metricLabel(m))}</b><span class="chip ${cls}">${esc(text)}</span>` +
      `<span class="muted num">yours: ${fmt(quant(mine, 0.5))} ${unit}</span></div>${strip(m)}</div>`;
  }
  function renderResults() {
    const box = $("#claims");
    $("#refs").innerHTML = "<span class='eyebrow'>Reference</span> " + REF.datasets.map((d) =>
      `<label for="ref-${d.id}"><input type="checkbox" id="ref-${d.id}" data-ref="${esc(d.id)}"${S.refOn.has(d.id) ? " checked" : ""}>${esc(d.label)} <span class="muted num">(${d.n})</span></label>`).join("");
    $("#refs").querySelectorAll("input").forEach((el) => el.addEventListener("change", () => { el.checked ? S.refOn.add(el.dataset.ref) : S.refOn.delete(el.dataset.ref); renderResults(); }));
    $("#forced-note").innerHTML = S.rows && S.forced ? `<p class="chip fail">Run forced despite a failed check. Treat these results as unreliable.</p>` : "";
    if (!S.rows) { box.innerHTML = `<p class="muted">Results appear here after step 3.</p>`; return; }
    const keys = new Set(S.rows.flatMap((r) => Object.keys(r)));
    const pick = (prefix) => [...keys].filter((k) => k.startsWith(prefix + "_") && k.endsWith("_diff") && !k.includes("bleed")).sort();
    box.innerHTML = CLAIMS.map((c) => {
      let mains = Array.isArray(c.main) ? c.main : pick(c.main);
      if (c.also) mains = mains.concat(pick(c.also));
      mains = mains.filter((k) => keys.has(k));
      if (!mains.length) return "";
      const ctrls = (c.control || []).filter((k) => keys.has(k) || REF.rows.some((r) => k in r));
      return `<article class="claim${c.isControl ? " control" : ""}"><div class="claim-head"><h3>${esc(c.title)}</h3>${c.privateOnly ? '<span class="chip yours">private-only</span>' : ""}${c.isControl ? '<span class="chip info">control</span>' : ""}</div>` +
        `<p class="muted">${esc(c.desc)}</p><div class="plots">${mains.map((m) => plotBlock(m, c.unit)).join("")}</div>` +
        (ctrls.length ? `<p class="hint">${esc(c.ctrlText || "")}</p><div class="plots">${ctrls.map((m) => plotBlock(m, c.unit)).join("")}</div>` : "") + "</article>";
    }).join("");
  }

  // ---- summary export: random patient codes (not hashes: hashes of short IDs can be reversed by guessing) ----
  function randomCode() { const b = new Uint8Array(4); crypto.getRandomValues(b); return "P-" + [...b].map((x) => x.toString(16).padStart(2, "0")).join(""); }
  function renderSummary() {
    const cols = [...new Set(S.rows.flatMap((r) => Object.keys(r)))].filter((k) => k !== "patient");
    const used = new Set(), key = [["code", "your_patient_id"].join(",")];
    const lines = [["patient"].concat(cols).join(",")];
    const q = (s) => `"${String(s ?? "").replace(/"/g, '""')}"`;
    for (const r of S.rows) {
      let code; do code = randomCode(); while (used.has(code)); used.add(code);
      key.push([code, q(r.patient)].join(","));
      lines.push([code].concat(cols.map((k) => (typeof r[k] === "number" ? (Number.isNaN(r[k]) ? "" : r[k]) : q(r[k])))).join(","));
    }
    $("#summary").value = lines.join("\n"); $("#key").value = key.join("\n");
  }
  $("#btn-copy").addEventListener("click", () => {
    const t = $("#summary").value;
    navigator.clipboard.writeText(t).then(() => ($("#copy-status").textContent = "Copied."), () => { $("#summary").select(); $("#copy-status").textContent = "Text selected; press Ctrl/Cmd+C to copy."; });
  });
  $("#btn-download").addEventListener("click", () => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([$("#summary").value], { type: "text/csv" }));
    a.download = "summary.csv"; document.body.appendChild(a); a.click(); a.remove();
  });

  // ---- boot ----
  if (window.__HOSTED) $("#btn-download").hidden = true;  // hosted viewers block downloads; Copy still works
  $("#recipe-version").textContent = REF.recipe;
  $("#ref-count").textContent = REF.datasets.filter((d) => d.id !== "CODEX_CRC_Schurch2020").reduce((a, d) => a + d.n, 0) + "+";
  setStepEnabled("#step-map", false); setStepEnabled("#step-check", false); setStepEnabled("#step-compare", false);
  renderChecks(); renderResults();
  if (EXAMPLE) load({ name: "example_codex_6_patients.csv", text: EXAMPLE }, true).catch(showError);
})();
