// Tumor Neighborhood Benchmark: UI. Loaded files are read and analysed in a Web Worker; nothing is uploaded.
(function () {
  "use strict";
  const $ = (s, el = document) => el.querySelector(s);
  const REF = JSON.parse(document.getElementById("ref-data").textContent);
  const WORKER_SRC = ["recipe-src", "automap-src"].map((id) => document.getElementById(id).textContent).join("\n;\n");
  const LOADED_COLORS = ["var(--accent)", "var(--c2)", "var(--c3)", "var(--c4)"];
  // One hue per study, one shade per dataset; a public dataset not listed here falls back to the reference green.
  const DS_COLOR = {
    IMMUcan_IMC_panel1: "--ds-immucan", CODEX_CRC_Schurch2020: "--ds-codex",
    SPATCH_COAD_xenium5k: "--ds-spatch-1", SPATCH_COAD_cosmx6k: "--ds-spatch-2", SPATCH_OV_xenium5k: "--ds-spatch-3", SPATCH_HCC_xenium5k: "--ds-spatch-4",
    MARTEAU_CRC_IMC_core: "--ds-marteau-1", MARTEAU_CRC_IMC_margin: "--ds-marteau-2",
    CERVILLA_COLORECTAL_xenium5k: "--ds-cervilla-1", CERVILLA_LUNG_xenium5k: "--ds-cervilla-2", CERVILLA_BREAST_xenium5k: "--ds-cervilla-3", CERVILLA_OVARIAN_xenium5k: "--ds-cervilla-4",
    XENIUM10X_LUNG_prime5k: "--ds-10x-lung",
    MARTEAU_XENIUM_core: "--ds-mxen-1", MARTEAU_XENIUM_margin: "--ds-mxen-2", MARTEAU_XENIUM_normal: "--ds-mxen-3",
  };
  const HIGHLIGHT = new Set(["down", "up", "change"]);  // findings shown on the card; the rest sit under "more"
  const openMore = new Set();
  const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const fmt = (v, d = 2) => (v === null || v === undefined || !Number.isFinite(v) ? "—" : (v > 0 ? "+" : "") + v.toFixed(d));

  // ---- datasets ----
  const DS = [];  // {id, label, short, kind, status, rows, findings, sens, report, checks, fails, meta, color, n, progress}
  const selected = new Set();
  for (const d of REF.datasets) {
    const rows = d.data.map((v) => { const r = { dataset: d.id }; d.cols.forEach((c, i) => { r[c] = v[i]; }); return r; });  // columnar in the page, row objects here
    DS.push({ id: d.id, label: d.label, short: d.short, kind: "public", status: "ready", rows, n: d.n, findings: Findings.summarize(rows), color: `var(${DS_COLOR[d.id] || "--ref"})` });
  }
  ["IMMUcan_IMC_panel1", "CODEX_CRC_Schurch2020"].forEach((id) => selected.add(id));

  // ---- worker: the whole automatic pipeline ----
  const workerBody = function () {
    const progress = (id, stage, f) => self.postMessage({ id, progress: { stage, f } });
    const asBlob = (f) => ({ name: f.name, size: f.size, slice: (...a) => f.slice(...a), stream: () => f.stream() });
    const toCSV = (header, rows) => [header].concat(rows).map((r) => r.map((v) => (/[",\n\r]/.test(v) ? '"' + String(v).replace(/"/g, '""') + '"' : v)).join(",")).join("\n");
    const held = new Map();  // datasets waiting for a pixel size
    let key = null;
    // A per-patient results table written by the Python kit: shown as it is, nothing is computed (no pixel size question either).
    async function tryKitResults(id, files) {
      const isKit = [];
      for (const f of files) isKit.push(AutoMap.isKitHeader(Recipe.csvHeader(await f.slice(0, 1 << 18).text()), f.size));
      if (!isKit.some(Boolean)) return null;
      progress(id, "Reading results table", 0);
      const parts = [];
      for (let i = 0; i < files.length; i++) {
        const out = isKit[i] ? AutoMap.readKitResults(await files[i].slice(0, files[i].size).text(), files[i].name, Recipe.VERSION) : null;
        if (!out) {
          if (files.length === 1 && isKit[i]) return null;  // too many rows for a results table: read it as a cell table
          return { fails: [{ title: files.length === 1 ? `${files[i].name} is not a results table` : "These files mix results tables and cell tables", why: "A kit results table has one row per patient. This file has too many rows.", fix: "Add kit results tables and cell tables as separate datasets." }] };
        }
        parts.push(out);
      }
      const merged = new Map();  // dataset id -> {rows, report, checks}; several files with the same dataset id are stacked
      for (const p of parts) for (const d of p.datasets) {
        const m = merged.get(d.id);
        if (!m) merged.set(d.id, { ...d }); else { m.rows = m.rows.concat(d.rows); m.checks = d.checks; m.report = d.report; }
      }
      if (parts.length > 1) for (const m of merged.values()) {  // recompute the patient check on the stacked rows
        const dup = m.rows.length - new Set(m.rows.map((r) => r.patient)).size;
        if (dup) m.checks = m.checks.map((c) => (c[0] === "patients" ? [c[0], "warn", `${m.rows.length - dup} patients; ${dup} rows repeat a patient already listed (files stacked)`, "Check that the files hold different runs."] : c));
      }
      const sets = [...merged.values()].filter((d) => d.rows.length);
      if (!sets.length) return { fails: [{ title: "The results table has no patient rows", why: "Every row had an empty patient ID.", fix: "Check that the file is the kit's patients_full.csv or share_summary.csv." }] };
      return { kitResults: sets };
    }
    async function analyse(id, fileList) {
      const files = fileList.map(asBlob), fails = [];
      for (const f of files) {
        if (f.size === 0) { fails.push({ title: `${f.name} is empty`, why: "The file has no content.", fix: "Export the cell table again." }); continue; }
        const bad = Recipe.sniff(new Uint8Array(await f.slice(0, 4096).arrayBuffer()), f.name);
        if (bad) fails.push({ title: bad.title, why: "", fix: bad.fix });
      }
      if (fails.length) return { fails };
      const kit = await tryKitResults(id, files);
      if (kit) return kit;
      progress(id, "Sampling", 0);
      const samples = [];
      for (const f of files) {
        const s = await AutoMap.sampleRows((a, b) => f.slice(a, b).text(), f.size);
        if (s.header.length < 4) return { fails: [{ title: `${f.name}: found ${s.header.length} column${s.header.length === 1 ? "" : "s"}`, why: "No comma, tab or semicolon separated columns.", fix: "Export the cell table as CSV or tab-separated text." }] };
        samples.push({ name: f.name, size: f.size, ...s });
      }
      const am = AutoMap.autoMap(samples);
      if (am.fails.length) return { report: am.report, fails: am.fails };
      const cfg = am.config;
      const counts = new Map();
      if (cfg.columns.type) {
        let done = 0; const total = files.reduce((a, f) => a + f.size, 0);
        for (const f of files) {
          const d = Recipe.detectDelim(Recipe.prepText(await f.slice(0, 1 << 16).text()));
          let H = null;
          const p = Recipe.createParser(d, (r) => { if (H === null) { H = r.map((x) => x.replace(/^﻿/, "").trim()).indexOf(cfg.columns.type); return; } if (H >= 0) counts.set(r[H], (counts.get(r[H]) || 0) + 1); });
          const rd = f.stream().getReader(), dec = new TextDecoder();
          for (;;) { const { value, done: dn } = await rd.read(); if (dn) break; done += value.length; p.feed(dec.decode(value, { stream: true })); progress(id, "Reading cell types", done / total); }
          p.feed(dec.decode()); p.end();
        }
      }
      let spacing = null;
      if (!cfg._unitKnown) spacing = Recipe.medianSpacing(Recipe.importCells(samples.map((s) => ({ name: s.name, text: toCSV(s.header, s.rows) })), { ...cfg, pixel_size_um: 1, roles: {} }));
      const fin = AutoMap.finishConfig(cfg, [...counts.entries()].sort((a, b) => b[1] - a[1]), spacing);
      const report = am.report.concat(fin.report);
      if (fin.fails.length) return { report, fails: fin.fails };
      if (cfg._pxEstimated) {  // the file doesn't state units: ask the user before measuring distances
        held.set(key, { files, cfg, report });
        return { report, needPx: { estimate: cfg._pxEstimate, candidates: cfg._pxCandidates } };
      }
      return measure(id, files, cfg, report);
    }
    async function measure(id, files, cfg, report) {
      const cells = await Recipe.importFiles(files, cfg, (f) => progress(id, "Reading cells", f));
      progress(id, "Checking", 1);
      const checks = Recipe.runChecks(cells, cfg);
      const failed = checks.filter((c) => c[1] === "fail");
      const meta = { labels: cfg._labels || null, rules: cfg.phenotype_rules || null, cells: cells.n };
      if (failed.length) return { report, checks, meta, fails: failed.map((c) => ({ title: c[0], why: c[2], fix: c[3] })) };
      const rows = Recipe.patientRows(cells, cfg, (f) => progress(id, "Measuring distances", f));
      let alt = null;
      if (cfg._pxEstimated && cells.n <= 3e6) {
        alt = [];
        for (const fct of [0.8, 1.25]) {
          const c2 = { ...cells, x: cells.x.map((v) => v * fct), y: cells.y.map((v) => v * fct) };
          alt.push(Recipe.patientRows(c2, cfg, (f) => progress(id, `Re-measuring at ${fct < 1 ? "−20%" : "+25%"} pixel size`, f)));
        }
      }
      return { report, checks, rows, alt, meta, pxEstimated: !!cfg._pxEstimated, big: cells.n > 3e6 };
    }
    async function resume(id, k, px) {
      const h = held.get(k); held.delete(k);
      if (!h) throw new Error("This dataset is no longer waiting for a pixel size. Add the file again.");
      if (px) {
        h.cfg.pixel_size_um = px; h.cfg._pxEstimated = false;
        h.report = h.report.filter((r) => r.what !== "Units").concat([{ what: "Units", decision: `${px} µm per pixel`, evidence: "entered for this dataset (the file doesn't state units)", confidence: "high" }]);
      }
      return measure(id, h.files, h.cfg, h.report);
    }
    self.onmessage = async (e) => {
      const { id } = e.data; key = e.data.key;
      try { self.postMessage({ id, ok: true, out: e.data.cmd === "resume" ? await resume(id, key, e.data.px) : await analyse(id, e.data.files) }); }
      catch (err) { self.postMessage({ id, ok: false, error: String(err && err.message || err), fix: err && err.fix || "" }); }
    };
  };
  const worker = new Worker(URL.createObjectURL(new Blob([WORKER_SRC + "\n;(" + workerBody.toString() + ")();"], { type: "text/javascript" })));
  let seq = 0; const pending = new Map();
  worker.onmessage = (e) => {
    const p = pending.get(e.data.id); if (!p) return;
    if (e.data.progress) { p.onProgress(e.data.progress); return; }
    pending.delete(e.data.id);
    if (e.data.ok) p.res(e.data.out); else { const err = new Error(e.data.error); err.fix = e.data.fix; p.rej(err); }
  };
  worker.onerror = (e) => {
    for (const p of pending.values()) { const err = new Error("The analysis stopped" + (e.message ? `: ${e.message}` : "") + "."); err.fix = "The browser may have run out of memory. Split the cohort into groups of patients and add each group as its own dataset."; p.rej(err); }
    pending.clear();
  };
  const send = (msg, onProgress) => new Promise((res, rej) => { const id = ++seq; pending.set(id, { res, rej, onProgress }); worker.postMessage({ id, ...msg }); });
  const analyse = (key, files, onProgress) => send({ cmd: "analyse", key, files }, onProgress);
  const resume = (key, px, onProgress) => send({ cmd: "resume", key, px }, onProgress);

  // ---- add data ----
  function datasetName(files) {
    const base = files.map((f) => f.name.replace(/\.[^.]+$/, ""));
    if (base.length === 1) return base[0];
    let p = base[0]; for (const b of base) while (!b.startsWith(p)) p = p.slice(0, -1);
    p = p.replace(/[_\-. ]+$/, "");
    return (p.length >= 3 ? p : "dataset") + ` (${files.length} files)`;
  }
  let queue = Promise.resolve();
  function addFiles(files) {
    if (!files.length) return;
    const nLoaded = DS.filter((d) => d.kind !== "public").length;
    const d = { id: "loaded-" + Date.now() + "-" + nLoaded, label: datasetName(files), kind: "loaded", status: "analysing", color: LOADED_COLORS[nLoaded % LOADED_COLORS.length],
      progress: { stage: "Waiting", f: 0 } };
    d.short = d.label.length > 16 ? d.label.slice(0, 15) + "…" : d.label;
    DS.push(d); render();
    queue = queue.then(() => settle(d, analyse(d.id, files, (p) => { d.progress = p; renderPanel(); })));  // one dataset at a time
  }
  // A kit results table replaces its placeholder with one dataset per `dataset` value in the file.
  function adoptResults(placeholder, sets) {
    const at = DS.indexOf(placeholder), others = DS.filter((x) => x !== placeholder && x.kind !== "public").length;
    const entries = sets.map((s, i) => {
      let id = "results-" + s.id, k = 2; while (DS.some((x) => x.id === id)) id = `results-${s.id}-${k++}`;
      const rows = s.rows;
      const e = { id, label: s.id, short: s.id.length > 16 ? s.id.slice(0, 15) + "…" : s.id, kind: "results", status: "ready", rows, n: new Set(rows.map((r) => r.patient)).size,
        findings: Findings.summarize(rows, { loaded: true, results: true }), sens: new Set(), report: s.report, checks: s.checks, fails: [], meta: null,
        color: LOADED_COLORS[(others + i) % LOADED_COLORS.length] };
      selected.add(id);
      return e;
    });
    DS.splice(at < 0 ? DS.length : at, 1, ...entries);
  }
  async function settle(d, job) {
    try {
      const out = await job;
      if (out.kitResults) { adoptResults(d, out.kitResults); render(); return; }
      Object.assign(d, out, { needPx: out.needPx || null });
      if (out.fails && out.fails.length) d.status = "failed";
      else if (out.needPx) d.status = "needs-px";
      else {
        d.status = "ready"; d.n = d.rows.length; d.findings = Findings.summarize(d.rows, { loaded: true });
        d.sens = d.alt ? Findings.sensitivity(d.findings, d.alt) : new Set();
        selected.add(d.id); renderSummary(d);
      }
    } catch (err) { d.status = "failed"; d.fails = [{ title: err.message, why: "", fix: err.fix || "" }]; }
    render();
  }
  function answerPx(d, px) {
    d.status = "analysing"; d.progress = { stage: "Measuring", f: 0 }; render();
    queue = queue.then(() => settle(d, resume(d.id, px, (p) => { d.progress = p; renderPanel(); })));
  }
  $("#pick").addEventListener("click", () => $("#file").click());
  $("#file").addEventListener("change", (e) => { addFiles([...e.target.files]); e.target.value = ""; });
  const drop = $("#drop");
  ["dragenter", "dragover"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add("over"); }));
  ["dragleave", "drop"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
  drop.addEventListener("drop", (e) => addFiles([...e.dataTransfer.files]));

  // ---- left panel ----
  const STATUS_CHIP = { ready: ["pass", "ready"], analysing: ["info", "analysing"], failed: ["fail", "couldn't analyse"], "needs-px": ["warn", "needs pixel size"] };
  function pxForm(d) {
    const np = d.needPx, near = np.candidates.filter((c) => c.near);
    return `<form class="pxform" data-px="${esc(d.id)}"><p><b>How big is one pixel?</b> This file gives positions in pixels but doesn't say how large a pixel is. ` +
      `It sets the scale for the 10 µm and 30 µm distances. Your scanner or microscope settings show it (often in µm per pixel).</p>` +
      (near.length ? `<div class="row">${near.map((c) => `<button type="button" data-pick="${c.um}">${c.um} µm · ${esc(c.name)}</button>`).join("")}</div>` : "") +
      `<div class="row"><label for="px-${esc(d.id)}">µm per pixel</label><input type="number" id="px-${esc(d.id)}" step="0.0001" min="0.05" max="5" value="${np.estimate}" required>` +
      `<button type="submit" class="primary">Use this value</button></div>` +
      `<p class="hint">Cell spacing suggests about ${np.estimate} µm. Not sure? <button type="button" class="linkish" data-estimate="1">Use the estimate</button>: results are then also recomputed at ±20% and fragile findings are flagged.</p></form>`;
  }
  // Sidebar blurb, generated from the embedded reference (fields source/detail/unit/part come from SOURCE in web/build.py).
  // Studies where every dataset has fewer patients than the per-patient rule needs (MIN_N = 5) are merged into one "shown for comparison" item.
  function refNote(datasets) {
    const WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"], MIN_N = 5;
    const groups = new Map();
    for (const d of datasets) {
      if (!groups.has(d.source)) groups.set(d.source, { source: d.source, detail: d.detail || "", unit: d.unit || "patient", ds: [] });
      groups.get(d.source).ds.push(d);
    }
    const plural = (n, unit) => `${n} ${unit}${n === 1 ? "" : "s"}`;
    const big = [...groups.values()].filter((g) => g.ds.some((d) => d.n >= MIN_N)), small = [...groups.values()].filter((g) => !big.includes(g));
    const items = big.map((g, i) => {
      const parts = g.ds.filter((d) => d.part), n = g.ds.reduce((a, d) => a + d.n, 0);
      const count = parts.length > 1 ? `${parts.map((d) => `${d.n} ${esc(d.part)}`).join(" and ")} ${g.unit}s` : i === 0 ? plural(n, g.unit) : String(n);
      return `${esc(g.source)} (${count}${g.detail ? ", " + esc(g.detail) : ""})`;
    });
    if (small.length) {
      const all = small.flatMap((g) => g.ds), same = all.every((d) => d.n === all[0].n), unit = small[0].unit;
      items.push(`${small.map((g) => esc(g.source)).join(" and ")} sections (${same ? plural(all[0].n, unit) + " each" : plural(all.reduce((a, d) => a + d.n, 0), unit) + " in total"}, shown for comparison)`);
    }
    const total = datasets.reduce((a, d) => a + d.n, 0), k = groups.size;
    return `<li class="group-note"><b>Public reference</b>: ${total} patients from ${WORDS[k] || k} published studies, all measured with the same recipe: ${items.join("; ")}. ` +
      `Only per-patient summary numbers are included, not the cell data.</li>`;
  }
  const REF_NOTE = refNote(REF.datasets);
  function renderPanel() {
    // Public datasets too small for summary findings (under 5 patients) go to the bottom, below loaded data.
    const tooSmall = (d) => d.kind === "public" && d.n < 5;
    $("#dsets").innerHTML = REF_NOTE + DS.filter((d) => !tooSmall(d)).concat(DS.filter(tooSmall)).map((d) => {
      const [cls, txt] = STATUS_CHIP[d.status];
      const box = `<input type="checkbox" id="sel-${esc(d.id)}" data-id="${esc(d.id)}"${selected.has(d.id) ? " checked" : ""}${d.status !== "ready" ? " disabled" : ""}>`;
      let body = "";
      if (d.status === "needs-px") body = pxForm(d);
      else if (d.status === "analysing") body = `<div class="bar"><i style="width:${Math.round((d.progress.f || 0) * 100)}%"></i></div><span class="hint">${esc(d.progress.stage)}…</span>`;
      else if (d.status === "failed") body = `<ul class="fails">${(d.fails || []).slice(0, 3).map((f) => `<li><b>${esc(f.title)}</b>${f.fix ? `<br><span class="hint">${esc(f.fix)}</span>` : ""}</li>`).join("")}</ul>`;
      else if (d.findings.every((f) => f.status === "na" || f.key === "compartments")) body = `<p class="hint">${d.n < 5 ? `Only ${d.n} patient${d.n === 1 ? "" : "s"}: shown in the charts, too few for summary findings.` : "No summary findings for this dataset."}</p>`;
      else {
        const li = (f) => `<li class="f-${f.status}">${esc(f.text)}${d.sens && d.sens.has(f.key) ? ' <span class="chip warn">depends on estimated pixel size</span>' : ""}</li>`;
        const shown = d.findings.filter((f) => f.status !== "na");
        let top = shown.filter((f) => HIGHLIGHT.has(f.status)).slice(0, 3); if (!top.length) top = shown.slice(0, 2);
        const rest = shown.filter((f) => !top.includes(f));
        body = `<ul class="finds">${top.map(li).join("")}</ul>` + (rest.length ? `<details class="more-finds" data-more="${esc(d.id)}"${openMore.has(d.id) ? " open" : ""}><summary>${rest.length} more detail${rest.length === 1 ? "" : "s"}</summary><ul class="finds">${rest.map(li).join("")}</ul></details>` : "");
      }
      return `<li class="dset${selected.has(d.id) ? " on" : ""}"><div class="dset-head"><label for="sel-${esc(d.id)}">${box}<i class="swatch" style="background:${d.color}"></i><span class="dname">${esc(d.label)}</span></label>` +
        `<span class="chip ${cls}">${txt}</span></div>` +
        `<div class="dmeta">${d.kind === "public" ? "public reference" : d.kind === "results" ? "kit results" : "loaded"}${d.n ? ` · ${d.n} patient${d.n === 1 ? "" : "s"}` : ""}` +
        `${d.kind !== "public" && d.status !== "analysing" ? ` · <button type="button" class="linkish" data-remove="${esc(d.id)}">remove</button>` : ""}</div>${body}</li>`;
    }).join("");
  }
  $("#dsets").addEventListener("toggle", (e) => { const id = e.target.dataset && e.target.dataset.more; if (id) e.target.open ? openMore.add(id) : openMore.delete(id); }, true);  // keep open across re-renders
  $("#dsets").addEventListener("change", (e) => { const id = e.target.dataset.id; if (!id) return; e.target.checked ? selected.add(id) : selected.delete(id); render(); });
  $("#dsets").addEventListener("click", (e) => {
    const form = e.target.closest("form.pxform");
    if (form && e.target.dataset.pick) { form.querySelector("input").value = e.target.dataset.pick; form.querySelector("input").focus(); return; }
    if (form && e.target.dataset.estimate) { answerPx(DS.find((d) => d.id === form.dataset.px), null); return; }
    const id = e.target.dataset.remove; if (!id) return; const i = DS.findIndex((d) => d.id === id); if (i >= 0) DS.splice(i, 1); selected.delete(id); render();
  });
  $("#dsets").addEventListener("submit", (e) => {
    const form = e.target.closest("form.pxform"); if (!form) return;
    e.preventDefault();
    const v = Number(form.querySelector("input").value);
    if (!(v >= 0.05 && v <= 5)) { form.querySelector("input").setCustomValidity("Enter a value between 0.05 and 5 µm per pixel."); form.reportValidity(); return; }
    answerPx(DS.find((d) => d.id === form.dataset.px), v);
  });
  $("#dsets").addEventListener("input", (e) => { if (e.target.matches(".pxform input")) e.target.setCustomValidity(""); });

  // ---- comparison plots ----
  // scope: "whole" = every cell in the image (recipe v2p.1 measures); "nest" = inside the tumor nest (recipe v3p).
  const NEST_DESC = "Same test, restricted to cells inside the tumor nest. The chance level comes from shuffling labels within each compartment.";
  const CLAIMS = [
    { title: "Do T cells keep away from neutrophils?", scope: "whole", desc: "Share of T cells within 30 µm of a neutrophil, compared with shuffled labels. Below 0 means fewer T cells nearby than chance.",
      main: ["T_near_neut_log2"], control: ["B_near_neut_log2"], ctrlText: "Control: B cells. If they keep away just as much, the pattern is about lymphocyte zones, not T cells specifically.", unit: "log2" },
    { title: "Do T cells keep away from neutrophils?", scope: "nest", desc: "Share of T cells inside the tumor nest that have a neutrophil within 30 µm, compared with labels shuffled within each compartment. The same for CD8, CD4 and regulatory T cells.",
      main: ["nest_T_near_neut_log2", "nest_CD8_near_neut_log2", "nest_CD4_near_neut_log2", "nest_Treg_near_neut_log2"], control: ["nest_B_near_neut_log2"], ctrlText: "Control: B cells in the nest.", unit: "log2" },
    { title: "Do T cells gather near macrophages?", scope: "whole", desc: "Share of T cells within 30 µm of a macrophage, compared with shuffled labels. Above 0 means more T cells nearby than chance.", main: ["T_near_mac_log2"], unit: "log2" },
    { title: "Do T cells gather near macrophages?", scope: "nest", desc: "Same test inside the tumor nest, for all T cells and for CD8, CD4 and regulatory T cells.",
      main: ["nest_T_near_mac_log2", "nest_CD8_near_mac_log2", "nest_CD4_near_mac_log2", "nest_Treg_near_mac_log2"], unit: "log2" },
    { gradient: true },
    { title: "Are T cells touching neutrophils different?", scope: "whole", desc: "Percentile-point difference in each readout: T cells under 10 µm from a neutrophil minus T cells 30 µm or more away.",
      main: "neutTouch", control: ["neutTouch_bleed_neut_diff"], ctrlText: "Bleed: the neutrophil's own marker measured inside touching T cells. Large values mean contact comparisons are unreliable.", unit: "pct-pts" },
    { title: "Are T cells touching neutrophils different?", scope: "nest", desc: "Same contact comparison for T cells inside the tumor nest.", main: "nest_neutTouch", bySubset: true, control: ["nest_neutTouch_bleed_neut_diff"], ctrlText: "Bleed: the neutrophil marker inside touching T cells in the nest.", unit: "pct-pts" },
    { title: "Are T cells touching macrophages different?", scope: "whole", desc: "Same contact comparison for macrophages.", main: "macTouch", control: ["macTouch_bleed_mac_diff"], ctrlText: "Bleed: the macrophage marker inside touching T cells.", unit: "pct-pts" },
    { title: "Are T cells touching macrophages different?", scope: "nest", desc: "Same contact comparison for T cells inside the tumor nest.", main: "nest_macTouch", bySubset: true, control: ["nest_macTouch_bleed_mac_diff"], ctrlText: "Bleed: the macrophage marker inside touching T cells in the nest.", unit: "pct-pts" },
    { title: "Do macrophages touching T cells look suppressive?", scope: "whole", desc: "Percentile-point difference in each macrophage marker: macrophages under 10 µm from a T cell minus macrophages 30 µm or more from any T cell. All macrophages in the image.",
      main: "TmacTouch", unit: "pct-pts" },
    { title: "Do macrophages touching T cells look suppressive?", scope: "nest", desc: "Same comparison for macrophages inside the tumor nest only.",
      main: "nest_TmacTouch", unit: "pct-pts" },
    { title: "Control: T cells touching stroma", scope: "whole", desc: "Same contact comparison for stromal cells. Effects here flag technical contact artifacts rather than biology.", main: "stromaTouch", control: ["stromaTouch_bleed_stroma_diff"], unit: "pct-pts", isControl: true },
    { title: "Control: T cells touching stroma", scope: "nest", desc: "Same stromal control inside the tumor nest.", main: "nest_stromaTouch", bySubset: true, control: ["nest_stromaTouch_bleed_stroma_diff"], unit: "pct-pts", isControl: true },
    { title: "T cells near the neutrophil subset (loaded data only)", scope: "whole", desc: "A neutrophil subset found in loaded data (for example LOX-1+). Public datasets have no such label; all neutrophils are shown as the closest reference.",
      main: ["T_near_neut_subset_log2"], also: "neut_subsetTouch", control: ["T_near_neut_log2"], ctrlText: "Reference: all neutrophils.", unit: "log2" },
  ];
  const SCOPE_CHIP = { whole: '<span class="chip info">Whole image</span>', nest: '<span class="chip nest">Tumor nest</span>' };
  const STATE_NAME = { prolif: "dividing (Ki-67)", PD1: "PD-1", GZMB: "granzyme B", LAG3: "LAG-3", ICOS: "ICOS", apoptosis: "dying (cleaved PARP)" };
  const MAC_NAME = { PDL1: "PD-L1", IDO1: "IDO1", CD163: "CD163", HLADR: "HLA-DR", ARG1: "arginase 1", M2: "M2-like score" };
  const SUB_NAME = { CD8: "CD8 T cells", CD4: "CD4 T cells", Treg: "Tregs" };
  const GRAD = ["0-5", "5-10", "10-20", "20-30", "30-50", "50-100", "100+"];
  function metricLabel(m) {
    const mac = m.match(/^(?:nest_)?TmacTouch_(.+)_diff$/); if (mac) return `${MAC_NAME[mac[1]] || mac[1]} on macrophages touching T cells`;
    const t = m.match(/^(?:nest_)?(\w+?)Touch_(?:(CD8|CD4|Treg)_)?(.+)_diff$/);
    if (t) { const st = t[3].startsWith("bleed") ? "bleed (own marker in T cells)" : STATE_NAME[t[3]] || t[3]; return t[2] ? `${SUB_NAME[t[2]]}: ${st}` : st; }
    const p = m.match(/^nest_(T|B|CD8|CD4|Treg)_near_(neut|mac)_log2$/);
    if (p) return `${{ T: "T cells", B: "B cells" }[p[1]] || SUB_NAME[p[1]]} near ${p[2] === "neut" ? "neutrophils" : "macrophages"}`;
    return { B_near_neut_log2: "B cells near neutrophils", T_near_neut_log2: "T cells near neutrophils", T_near_mac_log2: "T cells near macrophages", T_near_neut_subset_log2: "T cells near the subset" }[m] || m;
  }
  const quant = (a, q) => { const s = [...a].sort((x, y) => x - y); const p = (s.length - 1) * q, lo = Math.floor(p); return s.length ? s[lo] + (s[Math.min(lo + 1, s.length - 1)] - s[lo]) * (p - lo) : NaN; };
  const vals = (rows, m) => (rows || []).map((r) => r[m]).filter((v) => typeof v === "number" && Number.isFinite(v));
  function niceTicks(lo, hi, n) {
    const step0 = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(step0))), f = step0 / mag;
    const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
    const out = []; for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9; t += step) out.push(+t.toFixed(10));
    return out;
  }
  function strip(m, sets) {
    const cols = sets.map((d) => [d, vals(d.rows, m)]).filter(([, v]) => v.length);
    if (!cols.length) return "";
    const all = cols.flatMap(([, v]) => v).concat([0]);
    let lo = quant(all, 0.01), hi = quant(all, 0.99); if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.08; lo -= pad; hi += pad;
    const W = 425, H = 212, L = 50, R = 10, T = 12, B = 50, cw = (W - L - R) / cols.length;  // 25% larger than before
    const y = (v) => T + (H - T - B) * (1 - (Math.min(Math.max(v, lo), hi) - lo) / (hi - lo));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(metricLabel(m))}">`;
    for (const t of niceTicks(lo, hi, 4)) svg += `<line x1="${L}" x2="${W - R}" y1="${y(t)}" y2="${y(t)}" stroke="var(--line)" stroke-width="${t === 0 ? 1.6 : 0.9}"/>` +
      `<text x="${L - 6}" y="${y(t) + 4}" text-anchor="end" font-size="12" fill="var(--muted)" font-family="var(--f-data)">${t}</text>`;
    cols.forEach(([d, v], i) => {  // median line first, so the patient dots sit on top of it
      const cx = L + cw * (i + 0.5), loaded = d.kind !== "public";
      const md = quant(v, 0.5);
      svg += `<line x1="${cx - cw * 0.32}" x2="${cx + cw * 0.32}" y1="${y(md)}" y2="${y(md)}" stroke="${d.color}" stroke-width="2.6" stroke-opacity="0.75"/>`;
      v.forEach((val, j) => { const jit = ((j * 0.6180339887) % 1 - 0.5) * cw * 0.5;
        svg += `<circle cx="${(cx + jit).toFixed(1)}" cy="${y(val).toFixed(1)}" r="${loaded ? 4 : 3.4}" fill="${d.color}" fill-opacity="0.85" stroke="var(--surface)" stroke-width="0.8"/>`; });
      const words = String(d.short).split(" ");
      svg += `<text x="${cx}" y="${H - B + 18}" text-anchor="middle" font-size="12" fill="${loaded ? "var(--ink)" : "var(--muted)"}" font-family="var(--f-body)">${esc(words[0])}</text>`;
      if (words.length > 1) svg += `<text x="${cx}" y="${H - B + 33}" text-anchor="middle" font-size="12" fill="var(--muted)" font-family="var(--f-body)">${esc(words.slice(1).join(" "))}</text>`;
    });
    return svg + "</svg>";
  }
  function verdicts(m, sets) {  // each selected loaded dataset vs the selected public datasets pooled
    const pub = sets.filter((d) => d.kind === "public").flatMap((d) => vals(d.rows, m));
    const anyPublic = sets.some((d) => d.kind === "public");
    return sets.filter((d) => d.kind !== "public").map((d) => {
      const v = vals(d.rows, m);
      let cls = "info", txt = "no data";
      if (v.length && v.length < 5) { cls = "warn"; txt = `only ${v.length} patient${v.length > 1 ? "s" : ""}`; }
      else if (v.length && pub.length < 5) txt = anyPublic ? "public reference not yet computed for this measure" : "no public counterpart selected";
      else if (v.length) { const md = quant(v, 0.5), lo = quant(pub, 0.1), hi = quant(pub, 0.9); [cls, txt] = md < lo ? ["fail", "below public range"] : md > hi ? ["fail", "above public range"] : ["pass", "within public range"]; }
      return `<span class="vd"><i class="swatch" style="background:${d.color}"></i>${esc(d.short)}: <span class="chip ${cls}">${esc(txt)}</span> <span class="num muted">median ${fmt(quant(v, 0.5))}</span></span>`;
    }).join("");
  }
  // Pooled estimate under the sparse nest contact charts (recipe v3p.1): shown only for datasets that carry the pooled sums.
  const POOL_MIN_CELLS = 50;   // same rule as the pooled sentences in findings.js
  const POOLED = new WeakMap();   // rows array -> {test: pooled()}
  function pooledFor(d, base) {
    let c = POOLED.get(d.rows); if (!c) POOLED.set(d.rows, (c = {}));
    return base in c ? c[base] : (c[base] = Findings.pooled(d.rows, base));
  }
  function pooledLine(m, sets) {
    if (!m.endsWith("_diff") || !/^(nest_\w+Touch|TmacTouch|nest_TmacTouch)_/.test(m)) return "";
    const base = m.slice(0, -5), parts = []; let est = false;
    for (const d of sets) {
      const p = pooledFor(d, base);
      if (!p || p.n_near < 1) continue;
      const cnt = `${p.n_near.toLocaleString("en-US")} cells, ${p.patients} patient${p.patients === 1 ? "" : "s"}`, dot = `<i class="swatch" style="background:${d.color}"></i>${esc(d.short)}`;
      if (p.n_near < POOL_MIN_CELLS || p.patients < 5 || !Number.isFinite(p.diff)) { parts.push(`<span class="num">${dot}</span> <span class="muted">too few to pool (${cnt})</span>`); continue; }
      const ci = Number.isFinite(p.lo) && Number.isFinite(p.hi) ? `, 95% CI ${fmt(p.lo, 1)} to ${fmt(p.hi, 1)}` : "";
      est = true;
      parts.push(`<span class="num">${dot} ${fmt(p.diff, 1)}</span>${ci} <span class="muted">(${cnt})</span>`);
    }
    return parts.length ? `<div class="hint pooled">Pooled across patients${est ? " (pct-pts)" : ""}: ${parts.join(" · ")}.${est ? " Pooled estimate, not per patient." : ""}</div>` : "";
  }
  const plotBlock = (m, sets) => `<div class="plot"><b>${esc(metricLabel(m))}</b>${strip(m, sets)}<div class="vds">${verdicts(m, sets)}</div>${pooledLine(m, sets)}</div>`;

  // ---- tumor compartments: composition table and distance-gradient chart (recipe v3p) ----
  const COMP_NOTE = "Tumor nest = tumor cells make up at least half of the cells within 50 µm. Chance level is computed by shuffling labels within each compartment. " +
    "On whole sections the nest can include normal epithelium when the labels do not separate malignant cells. A nest proximity value can read as avoidance when the structure lies inside the nest; " +
    "the distance gradient and contact tests are the more direct readouts there.";
  const hasComp = (d) => d.rows.some((r) => "comp_nest_cells" in r);
  const fmtPct = (v) => (v >= 10 ? String(Math.round(v)) : v.toFixed(1));
  function compositionTable(d) {
    const COLS = [["T", "T cells"], ["CD8", "CD8"], ["CD4", "CD4"], ["Treg", "Treg"], ["neut", "neutrophils"], ["mac", "macrophages"], ["B", "B cells"]];
    const cols = COLS.filter(([k]) => d.rows.some((r) => `comp_nest_pct_${k}` in r || `comp_margin_pct_${k}` in r || `comp_far_pct_${k}` in r));
    const head = `<tr><th>compartment</th><th class="num">cells</th>${cols.map(([, n]) => `<th class="num">${esc(n)}</th>`).join("")}</tr>`;
    const body = [["nest", "Tumor nest"], ["margin", "Margin (within 100 µm of the nest)"], ["far", "Far tissue"]].map(([c, name]) => {
      const cells = vals(d.rows, `comp_${c}_cells`);
      const tds = cols.map(([k]) => { const v = vals(d.rows, `comp_${c}_pct_${k}`); return v.length ? `<td class="num"><b>${fmtPct(quant(v, 0.5))}</b> <span class="muted">${fmtPct(quant(v, 0.25))}–${fmtPct(quant(v, 0.75))}</span></td>` : '<td class="num muted">—</td>'; });
      return `<tr><td>${name}</td><td class="num">${cells.length ? Math.round(quant(cells, 0.5)).toLocaleString() : "—"}</td>${tds.join("")}</tr>`;
    }).join("");
    return `<div class="comptable"><div class="comp-cap"><i class="swatch" style="background:${d.color}"></i><b>${esc(d.label)}</b> <span class="muted num">(${d.n} patients)</span></div>` +
      `<div class="tablewrap"><table>${head}${body}</table></div></div>`;
  }
  function compartmentBlock(sets, top) {
    const with_ = sets.filter(hasComp), without = sets.filter((d) => !hasComp(d));
    if (top && !with_.length) return "";
    if (!top && with_.length) return "";
    return `<article class="claim"><div class="claim-head"><h3>Who is in the tumor? Immune composition by compartment</h3>${SCOPE_CHIP.nest}</div>` +
      (with_.length ? `<p class="muted">Percent of the cells in each compartment (median across patients, then the middle half of patients). Each patient counts once.</p>${with_.map(compositionTable).join("")}` : "") +
      (without.length ? `<ul class="nocomp">${without.map((d) => `<li><i class="swatch" style="background:${d.color}"></i><b>${esc(d.label)}</b>: ${d.kind === "loaded" ? "no tumor cells were found, so compartment results are not available." : d.kind === "results" ? "this results table has no compartment results." : "compartment results not yet computed for this dataset."}</li>`).join("")}</ul>` : "") +
      `<p class="hint">${esc(COMP_NOTE)}</p></article>`;
  }
  const PUB_DASH = ["", "7 4", "2 4", "10 3 2 3"];
  function gradSvg(src, sets) {
    const pubSeen = [];
    const series = sets.map((d) => {
      const pts = GRAD.map((lab) => { const v = vals(d.rows, `grad_${src}_${lab}_log2`); return v.length >= 5 ? { med: quant(v, 0.5), lo: quant(v, 0.25), hi: quant(v, 0.75) } : null; });
      if (d.kind === "public") pubSeen.push(d);
      return { d, pts, dash: d.kind === "public" ? PUB_DASH[(pubSeen.length - 1) % PUB_DASH.length] : "" };
    }).filter((s) => s.pts.some(Boolean));
    if (!series.length) return "";
    const all = series.flatMap((s) => s.pts.filter(Boolean).flatMap((p) => [p.lo, p.hi, p.med])).concat([0]);
    let lo = Math.min(...all), hi = Math.max(...all); if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
    const pad = (hi - lo) * 0.1; lo -= pad; hi += pad;
    const W = 500, H = 320, L = 52, R = 78, T = 14, B = 70, x = (i) => L + ((W - L - R) * (i + 0.5)) / GRAD.length, y = (v) => T + (H - T - B) * (1 - (v - lo) / (hi - lo));
    let svg = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Nest T cells by distance to the nearest ${src === "neut" ? "neutrophil" : "macrophage"}">`;
    for (const t of niceTicks(lo, hi, 5)) svg += `<line x1="${L}" x2="${W - R}" y1="${y(t)}" y2="${y(t)}" stroke="var(--line)" stroke-width="${t === 0 ? 1.8 : 0.9}"/>` +
      `<text x="${L - 6}" y="${y(t) + 4}" text-anchor="end" font-size="12" fill="var(--muted)" font-family="var(--f-data)">${t}</text>`;
    GRAD.forEach((lab, i) => { svg += `<text x="${x(i)}" y="${H - B + 20}" text-anchor="middle" font-size="12" fill="var(--muted)" font-family="var(--f-data)">${lab === "100+" ? "100+" : lab.replace("-", "–")}</text>`; });
    svg += `<text x="${(L + W - R) / 2}" y="${H - 14}" text-anchor="middle" font-size="12" fill="var(--muted)" font-family="var(--f-body)">distance to the nearest ${src === "neut" ? "neutrophil" : "macrophage"} (µm)</text>` +
      `<text transform="translate(14 ${(T + H - B) / 2}) rotate(-90)" text-anchor="middle" font-size="12" fill="var(--muted)" font-family="var(--f-body)">log2 vs chance</text>`;
    for (const { d, pts, dash } of series) {
      const loaded = d.kind !== "public", col = d.color, runs = []; let cur = [];
      pts.forEach((p, i) => { if (p) cur.push([i, p]); else if (cur.length) { runs.push(cur); cur = []; } }); if (cur.length) runs.push(cur);
      for (const r of runs) {
        if (r.length > 1) svg += `<polygon points="${r.map(([i, p]) => `${x(i).toFixed(1)},${y(p.hi).toFixed(1)}`).concat(r.slice().reverse().map(([i, p]) => `${x(i).toFixed(1)},${y(p.lo).toFixed(1)}`)).join(" ")}" fill="${col}" fill-opacity="0.13"/>`;
        svg += `<polyline points="${r.map(([i, p]) => `${x(i).toFixed(1)},${y(p.med).toFixed(1)}`).join(" ")}" fill="none" stroke="${col}" stroke-width="3" ${dash ? `stroke-dasharray="${dash}"` : ""} stroke-linejoin="round"/>`;
        for (const [i, p] of r) svg += `<circle cx="${x(i).toFixed(1)}" cy="${y(p.med).toFixed(1)}" r="${loaded ? 5 : 4.5}" fill="${col}" stroke="var(--surface)" stroke-width="1.2"/>`;
      }
      const lastRun = runs[runs.length - 1], [li, lp] = lastRun[lastRun.length - 1], name = String(d.short);
      svg += `<text x="${(x(li) + 9).toFixed(1)}" y="${(y(lp.med) + 4).toFixed(1)}" font-size="12" fill="var(--ink)" font-family="var(--f-body)">${esc(name.length > 12 ? name.slice(0, 11) + "…" : name)}</text>`;
    }
    return svg + "</svg>";
  }
  function gradientBlock(sets) {
    const panels = [["neut", "neutrophil"], ["mac", "macrophage"]].map(([src, nm]) => ({ src, nm, svg: gradSvg(src, sets) })).filter((p) => p.svg);
    if (!panels.length) return "";
    return `<article class="claim"><div class="claim-head"><h3>How does the T-cell pattern change with distance?</h3>${SCOPE_CHIP.nest}</div>` +
      `<p class="muted">For T cells inside the tumor nest: how many have their nearest neutrophil (or macrophage) in each distance band, compared with labels shuffled within each compartment. ` +
      `Above 0 means more T cells at that distance than chance. The first four bands together are "within 30 µm". Lines are the median across patients; the band is the middle half of patients.</p>` +
      `<div class="plots wide">${panels.map((p) => `<div class="plot"><b>Nest T cells by distance to the nearest ${p.nm}</b>${p.svg}</div>`).join("")}</div>` +
      `<p class="hint">A band needs at least 5 patients. Association in tissue snapshots, not mechanism.</p></article>`;
  }
  function renderClaims() {
    const sets = DS.filter((d) => selected.has(d.id) && d.status === "ready");
    $("#legend").innerHTML = sets.map((d) => `<span><i class="swatch" style="background:${d.color}"></i>${esc(d.label)} <span class="muted num">(${d.n})</span></span>`).join("");
    if (!sets.length) { $("#claims").innerHTML = `<p class="muted">Tick one or more datasets on the left to see them here.</p>`; return; }
    const keys = new Set(sets.flatMap((d) => d.rows.flatMap((r) => Object.keys(r))));
    const pick = (prefix) => [...keys].filter((k) => k.startsWith(prefix + "_") && k.endsWith("_diff") && !k.includes("bleed")).sort();
    const subOf = (k, prefix) => { const m = k.slice(prefix.length + 1).match(/^(CD8|CD4|Treg)_/); return m ? m[1] : null; };
    const html = [compartmentBlock(sets, true)];
    for (const c of CLAIMS) {
      if (c.gradient) { html.push(gradientBlock(sets)); continue; }
      let mains = Array.isArray(c.main) ? c.main : pick(c.main);
      if (c.also) mains = mains.concat(pick(c.also));
      mains = mains.filter((k) => keys.has(k));
      let subs = [];
      if (c.bySubset) { subs = ["CD8", "CD4", "Treg"].map((s) => [s, mains.filter((k) => subOf(k, c.main) === s)]).filter(([, ks]) => ks.length); mains = mains.filter((k) => !subOf(k, c.main)); }
      if (!mains.length && !subs.length) continue;
      const ctrls = (c.control || []).filter((k) => keys.has(k));
      html.push(`<article class="claim${c.isControl ? " control" : ""}"><div class="claim-head"><h3>${esc(c.title)}</h3><span class="row">${SCOPE_CHIP[c.scope]}${c.isControl ? '<span class="chip info">control</span>' : ""}</span></div>` +
        `<p class="muted">${esc(c.desc)}</p>${mains.length ? `<div class="plots">${mains.map((m) => plotBlock(m, sets)).join("")}</div>` : ""}` +
        subs.map(([s, ks]) => `<details><summary>${esc(SUB_NAME[s])} only (${ks.length} readouts)</summary><div class="plots">${ks.map((m) => plotBlock(m, sets)).join("")}</div></details>`).join("") +
        (ctrls.length ? `<p class="hint">${esc(c.ctrlText || "")}</p><div class="plots">${ctrls.map((m) => plotBlock(m, sets)).join("")}</div>` : "") + "</article>");
    }
    html.push(compartmentBlock(sets, false));
    $("#claims").innerHTML = html.join("");
  }

  // ---- details for loaded datasets: how it was read (read-only), checks, shareable summary ----
  const CONF = { high: "pass", medium: "info", low: "warn" };
  function renderDetails() {
    const loaded = DS.filter((d) => d.kind !== "public" && d.status !== "analysing" && (selected.has(d.id) || d.status === "failed"));
    $("#details").innerHTML = loaded.map((d) => {
      const report = (d.report || []).map((r) => `<tr><td>${esc(r.what)}</td><td>${esc(r.decision)}</td><td class="muted">${esc(r.evidence)}</td><td><span class="chip ${CONF[r.confidence] || "info"}">${esc(r.confidence)}</span></td></tr>`).join("");
      const labels = d.meta && d.meta.labels ? `<details><summary>Cell-type labels (${d.meta.labels.length})</summary><div class="tablewrap"><table><tr><th>label</th><th class="num">cells</th><th>read as</th></tr>` +
        d.meta.labels.map(([l, n, r]) => `<tr><td class="mono">${esc(l)}</td><td class="num">${n.toLocaleString()}</td><td>${esc(r.startsWith("T:") ? `T cell (${r.slice(2)})` : r)}</td></tr>`).join("") + "</table></div></details>" : "";
      const checks = (d.checks || []).map(([n, s, det, fix]) => `<tr><td><span class="chip ${s}">${s === "fail" ? "failed" : s === "warn" ? "note" : "ok"}</span></td><td>${esc(n)}</td><td>${esc(det)}</td><td class="muted">${esc(fix || "")}</td></tr>`).join("");
      const fails = d.status === "failed" ? `<div class="failbox"><h4>Why this dataset couldn't be analysed</h4><ul>${(d.fails || []).map((f) => `<li><b>${esc(f.title)}</b>${f.why ? `<br>${esc(f.why)}` : ""}${f.fix ? `<br><span class="muted">What to do: ${esc(f.fix)}</span>` : ""}</li>`).join("")}</ul></div>` : "";
      const share = d.status === "ready" && d.kind === "loaded" && window.__HOSTED ? `<p class="hint">Summary export is available in the downloadable page (docs/index.html) and the Python kit.</p>` : d.status === "ready" && d.kind === "loaded" ? `<details><summary>Share a summary (optional)</summary><p class="hint">Per-patient effect sizes only, with each patient replaced by a random code. No cell-level data. Nothing is sent anywhere.</p>` +
        `<div class="row"><button type="button" data-copy="${esc(d.id)}">Copy summary</button><button type="button" data-dl="${esc(d.id)}">Download summary.csv</button><span class="muted" id="cs-${esc(d.id)}"></span></div>` +
        `<textarea class="summary" readonly aria-label="Summary CSV" id="sum-${esc(d.id)}">${esc(d.summary || "")}</textarea>` +
        `<details><summary>Key: code → patient ID in your file (keep private)</summary><textarea class="summary" readonly aria-label="Patient key">${esc(d.key || "")}</textarea></details></details>` : "";
      return `<section class="detail"><h3><i class="swatch" style="background:${d.color}"></i>${esc(d.label)}</h3>${fails}` +
        (report ? `<details${d.status === "failed" ? " open" : ""}><summary>How this file was read</summary><div class="tablewrap"><table><tr><th>what</th><th>decision</th><th>evidence</th><th>confidence</th></tr>${report}</table></div>${labels}</details>` : "") +
        (checks ? `<details${d.status === "failed" ? " open" : ""}><summary>Checks (${(d.checks || []).length})</summary><div class="tablewrap"><table>${checks}</table></div></details>` : "") +
        (d.pxEstimated ? `<p class="hint">The file doesn't state units, so the pixel size was estimated. ${d.big ? "This dataset is too large to recompute at other pixel sizes." : "Results were recomputed at −20% and +25%; findings that change are marked in the dataset list."}</p>` : "") +
        share + "</section>";
    }).join("");
  }
  function randomCode() { const b = new Uint8Array(4); crypto.getRandomValues(b); return "P-" + [...b].map((x) => x.toString(16).padStart(2, "0")).join(""); }
  function renderSummary(d) {
    const cols = [...new Set(d.rows.flatMap((r) => Object.keys(r)))].filter((k) => k !== "patient" && k !== "dataset");
    const q = (s) => `"${String(s ?? "").replace(/"/g, '""')}"`;
    const used = new Set(), key = ["code,patient_id_in_your_file"], lines = [["dataset", "patient"].concat(cols).join(",")];
    for (const r of d.rows) {
      let code; do code = randomCode(); while (used.has(code)); used.add(code);
      key.push([code, q(r.patient)].join(","));
      lines.push([q(d.label), code].concat(cols.map((k) => (typeof r[k] === "number" ? (Number.isFinite(r[k]) ? r[k] : "") : q(r[k])))).join(","));
    }
    d.summary = lines.join("\n"); d.key = key.join("\n");
  }
  $("#details").addEventListener("click", (e) => {
    const c = e.target.dataset.copy, dl = e.target.dataset.dl;
    if (c) { const ta = document.getElementById("sum-" + c); navigator.clipboard.writeText(ta.value).then(() => (document.getElementById("cs-" + c).textContent = "Copied."), () => { ta.select(); document.getElementById("cs-" + c).textContent = "Text selected: press Ctrl/Cmd+C."; }); }
    if (dl) { const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([document.getElementById("sum-" + dl).value], { type: "text/csv" })); a.download = "summary.csv"; document.body.appendChild(a); a.click(); a.remove(); }
  });

  function render() { renderPanel(); renderClaims(); renderDetails(); }
  $("#recipe-version").textContent = REF.recipe;  // the public reference is measured with the same recipe as the page
  render();
})();
