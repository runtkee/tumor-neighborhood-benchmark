// Top-level findings per dataset, from its per-patient rows. Same rules for public and loaded datasets.
// A direction counts as consistent when at least 70% of patients (and at least 5) agree (MIN_N, AGREE).
// A B-cell control "shows a weaker version" when its median has the same sign and 55-70% of patients agree (WEAK_AGREE).
(function (root) {
  "use strict";
  const MIN_N = 5, AGREE = 0.7, WEAK_AGREE = 0.55, BLEED_HIGH = 15, POOL_MIN_CELLS = 50, N_BOOT = 200;
  const STATE = { prolif: "dividing (Ki-67)", PD1: "PD-1", GZMB: "granzyme B", LAG3: "LAG-3", ICOS: "ICOS", apoptosis: "dying (cleaved PARP)" };
  const vals = (rows, k) => rows.map((r) => r[k]).filter((v) => typeof v === "number" && Number.isFinite(v));
  const median = (a) => { const s = [...a].sort((x, y) => x - y), n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2) : NaN; };
  const sgn = (x, d = 2) => (x > 0 ? "+" : "") + x.toFixed(d);
  const SUBS = ["CD8", "CD4", "Treg"], SUB_NAME = { CD8: "CD8 T cells", CD4: "CD4 T cells", Treg: "Tregs" };
  const MAC_STATE = { PDL1: "PD-L1", IDO1: "IDO1", CD163: "CD163", HLADR: "HLA-DR", ARG1: "arginase 1", M2: "M2-like score" };
  const GRAD = ["0-5", "5-10", "10-20", "20-30", "30-50", "50-100", "100+"];
  const binText = (lab) => (lab === "100+" ? "beyond 100 µm" : `${lab} µm`);
  const hasKey = (rows, k) => rows.some((r) => k in r);
  const pctText = (v) => (v >= 10 ? String(Math.round(v)) : v.toFixed(1)) + "%";

  // Per-patient cell counts for one compartment c ("nest" | "margin" | "far") and group k ("T", "neut", ...):
  // percent x cells / 100, for the patients that have both numbers. Counts are what the table and the nest sentence show by default.
  function perPatientCounts(rows, c, k) {
    const out = [];
    for (const r of rows) { const p = finite(r[`comp_${c}_pct_${k}`]), n = finite(r[`comp_${c}_cells`]); if (p !== null && n !== null) out.push((p * n) / 100); }
    return out;
  }
  const countText = (v) => Math.round(v).toLocaleString("en-US");   // cells are whole numbers

  // ---- pooled estimate for the sparse contact tests (recipe v3p.1; src/bio/pooled.py is the Python twin) ----
  // Same generators as web/recipe.js and src/bio/prng.py (tests/test_validation.js checks they stay identical).
  function fnv1a(s) {
    let h = 0x811c9dc5;
    for (const b of new TextEncoder().encode(s)) h = Math.imul(h ^ b, 0x01000193) >>> 0;
    return h >>> 0;
  }
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  const finite = (v) => (typeof v === "number" && Number.isFinite(v) ? v : null);
  // key = the test without suffix, e.g. "nest_neutTouch_PD1". Pooled diff = sum(near_sum)/sum(n) - sum(far_sum)/sum(far_n) over the
  // patients with n > 0 or far_n > 0; 95% interval = patient bootstrap (200 resamples, rng seeded by fnv1a("pooled|" + key)).
  // Returns null when no patient has the sum columns.
  function pooled(rows, key) {
    const P = [];
    for (const r of rows) {
      const v = ["_n", "_near_sum", "_far_n", "_far_sum"].map((x) => finite(r[key + x]));
      if (v.includes(null) || !(v[0] > 0 || v[2] > 0)) continue;
      P.push(v);
    }
    if (!P.length) return null;
    const est = (sn, ss, fn, fs) => (sn > 0 && fn > 0 ? ss / sn - fs / fn : NaN);
    let sn = 0, ss = 0, fn = 0, fs = 0;
    for (const [a, b, c, d] of P) { sn += a; ss += b; fn += c; fs += d; }
    const diff = est(sn, ss, fn, fs), rng = mulberry32(fnv1a("pooled|" + key)), m = P.length, vals = [];
    for (let i = 0; i < N_BOOT; i++) {
      let bn = 0, bs = 0, bf = 0, bg = 0;
      for (let j = 0; j < m; j++) { const [a, b, c, d] = P[Math.floor(rng() * m)]; bn += a; bs += b; bf += c; bg += d; }
      const e = est(bn, bs, bf, bg);
      if (Number.isFinite(e)) vals.push(e);
    }
    vals.sort((x, y) => x - y);
    const q = (f) => { const p = (vals.length - 1) * f, lo = Math.floor(p), hi = Math.min(lo + 1, vals.length - 1); return vals[lo] + (vals[hi] - vals[lo]) * (p - lo); };
    const ok = vals.length >= 2;
    return { diff, lo: ok ? q(0.025) : NaN, hi: ok ? q(0.975) : NaN, n_near: sn, n_far: fn, patients: P.filter((p) => p[0] > 0).length, used: m };
  }
  const minus = (v, d = 1) => (v < 0 ? "\u2212" : "") + Math.abs(v).toFixed(d);
  const ciText = (p) => (Number.isFinite(p.lo) && Number.isFinite(p.hi) ? ` (95% CI ${minus(p.lo)} to ${minus(p.hi)}${p.lo <= 0 && p.hi >= 0 ? ", includes 0" : ""})` : "");
  const excludesZero = (p) => Number.isFinite(p.lo) && Number.isFinite(p.hi) && (p.lo > 0 || p.hi < 0);
  // One sentence for tests that fail the per-patient rule (under 5 patients with 20 cells per arm). items: [label, pooled()] pairs;
  // subject: "T cells touching neutrophils". Names only the markers whose 95% CI excludes 0; the rest are counted.
  // Returns {text, status, sig} or null (nothing to say for this group); sig = Map(label -> +1 | -1) of the named markers.
  function pooledSentence(lead, subject, items, always, caution = "") {
    items = items.filter(([, p]) => p);
    if (!items.length) return null;
    const cells = Math.max(...items.map(([, p]) => p.n_near)), pats = Math.max(...items.map(([, p]) => p.patients));
    const L = lead ? `${lead}${subject}` : subject[0].toUpperCase() + subject.slice(1);
    if (cells < POOL_MIN_CELLS) return always ? { status: "info", text: `${L} are too rare to compare (${cells} cells).` } : null;
    if (pats < MIN_N) return always ? { status: "info", text: `${L}: too few patients to pool (${cells} cells in ${pats} patients).` } : null;
    const tested = items.filter(([, p]) => Number.isFinite(p.diff));
    if (!tested.length) return null;
    const hit = tested.filter(([, p]) => excludesZero(p)), rest = tested.length - hit.length;
    const named = hit.map(([lab, p]) => `${lab} ${minus(Math.abs(p.diff))} pct-pts ${p.diff < 0 ? "lower" : "higher"}${ciText(p)}`);
    const tail = !hit.length ? `no marker differs from 0 (${tested.length} tested)` : rest ? `${rest} other marker${rest === 1 ? "" : "s"} include${rest === 1 ? "s" : ""} 0` : "";
    return { status: hit.length ? "change" : "none", sig: new Map(hit.map(([lab, p]) => [lab, Math.sign(p.diff)])),
      text: `${L} are rare (${cells.toLocaleString("en-US")} cells in ${pats} patients); pooled across patients, ${hit.length ? [...named, ...(tail ? [tail] : [])].join("; ") : tail}. Pooled estimate, not per patient.${caution}` };
  }

  // direction of one metric across patients: "down" | "up" | "mixed" | "few"
  function direction(rows, k) {
    const v = vals(rows, k);
    if (v.length < MIN_N) return { dir: "few", n: v.length };
    const down = v.filter((x) => x < 0).length, up = v.filter((x) => x > 0).length, med = median(v);
    const dir = down / v.length >= AGREE && med < 0 ? "down" : up / v.length >= AGREE && med > 0 ? "up" : "mixed";
    return { dir, n: v.length, down, up, med };
  }

  // lead (optional): "Inside tumor nests: " puts the sentence in the nest scope; without it the sentences are the whole-image ones
  function proximity(rows, metric, what, controlMetric, lead = "") {
    const d = direction(rows, metric), W = lead ? what : what[0].toUpperCase() + what.slice(1);
    if (d.dir === "few") return { key: metric, status: "na", text: d.n === 0 ? `${lead}${W}: not measured in this dataset.` : `${lead}${W}: too few patients to say (${d.n}).` };
    const nn = `${d.dir === "down" ? d.down : d.dir === "up" ? d.up : Math.max(d.down, d.up)} of ${d.n} patients`;
    let text, status;
    if (d.dir === "down") { text = `${lead}${lead ? "fewer" : "Fewer"} ${what} than chance in ${nn} (median ${sgn(d.med)} log2).`; status = "down"; }
    else if (d.dir === "up") { text = `${lead}${lead ? "more" : "More"} ${what} than chance in ${nn} (median ${sgn(d.med)} log2).`; status = "up"; }
    else { text = `${lead}${W}: no consistent pattern (${d.down} of ${d.n} fewer, median ${sgn(d.med)}).`; status = "mixed"; }
    if (controlMetric && (status === "down" || status === "up")) {
      const c = direction(rows, controlMetric);
      if (c.dir === status) text += " B cells show the same pattern, so it isn't specific to T cells.";
      else if (c.dir !== "few") {
        // weaker version: the control's median has the same sign as the main result and 55-70% of patients agree
        const sign = status === "down" ? -1 : 1, agree = status === "down" ? c.down : c.up;
        if (Math.sign(c.med) === sign && agree / c.n >= WEAK_AGREE) text += ` B cells show a weaker version of it (${agree} of ${c.n}), so specificity is unclear.`;
        else text += " B cells don't, so it looks specific to T cells.";
      }
    }
    return { key: metric, status, text };
  }

  function contact(rows, prefix, partner, single) {
    const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))].filter((k) => k.startsWith(prefix + "Touch_") && k.endsWith("_diff") && !k.includes("bleed"));
    if (!keys.length) return null;
    const found = [], mixed = [];
    for (const k of keys) {
      const st = k.slice(prefix.length + 6, -5), d = direction(rows, k);
      if (d.dir === "down" || d.dir === "up") found.push(`${STATE[st] || st} ${d.dir === "down" ? "lower" : "higher"} (${d.dir === "down" ? d.down : d.up} of ${d.n} patients, median ${sgn(d.med, 1)} pct-pts)`);
      else if (d.dir === "mixed") mixed.push(st);
    }
    const bleed = vals(rows, `${prefix}Touch_bleed_${prefix}_diff`), b = bleed.length ? median(bleed) : NaN;
    const caution = Number.isFinite(b) && b > BLEED_HIGH ? ` Caution: ${single} signal bleeds into touching T cells (+${b.toFixed(0)} pct-pts), so these contact results are unreliable.` : "";
    if (!found.length && !mixed.length) return { key: prefix + "Touch", status: "na", text: `T cells touching ${partner}: too few patients to say.` };
    return { key: prefix + "Touch", status: found.length ? "change" : "none",
      text: (found.length ? `T cells touching ${partner}: ${found.join("; ")}.` : `T cells touching ${partner}: no consistent change in ${mixed.length} readouts.`) + caution };
  }

  // ---- recipe v3p: tumor nest, composition, macrophage state, distance gradient ----
  function nestContact(rows, src, partner, single) {
    const pre = `nest_${src}Touch_`;
    const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))].filter((k) => k.startsWith(pre) && k.endsWith("_diff") && !k.includes("bleed"));
    const groups = { all: { found: [], mixed: [], few: [] } }; for (const s of SUBS) groups[s] = { found: [], mixed: [], few: [] };
    for (const k of keys) {
      const rest = k.slice(pre.length, -5), m = rest.match(/^(CD8|CD4|Treg)_(.+)$/), g = groups[m ? m[1] : "all"], st = m ? m[2] : rest, d = direction(rows, k);
      if (d.dir === "down" || d.dir === "up") g.found.push(`${STATE[st] || st} ${d.dir === "down" ? "lower" : "higher"} (${d.dir === "down" ? d.down : d.up} of ${d.n} patients, median ${sgn(d.med, 1)} pct-pts)`);
      else if (d.dir === "mixed") g.mixed.push(st);
      else if (hasKey(rows, k.slice(0, -5) + "_near_sum")) g.few.push([STATE[st] || st, k.slice(0, -5)]);   // v3p.1: too few patients per patient, pooled instead
    }
    const out = [], a = groups.all;
    if (a.found.length || a.mixed.length) {
      const bleed = vals(rows, `${pre}bleed_${src}_diff`), b = bleed.length ? median(bleed) : NaN;
      const caution = Number.isFinite(b) && b > BLEED_HIGH ? ` Caution: ${single} signal bleeds into touching T cells (+${b.toFixed(0)} pct-pts), so these contact results are unreliable.` : "";
      out.push({ key: `nest_${src}Touch`, status: a.found.length ? "change" : "none",
        text: (a.found.length ? `Inside tumor nests, T cells touching ${partner}: ${a.found.join("; ")}.` : `Inside tumor nests, T cells touching ${partner}: no consistent change in ${a.mixed.length} readouts.`) + caution });
    }
    for (const s of SUBS) if (groups[s].found.length) out.push({ key: `nest_${src}Touch_${s}`, status: "change", text: `Inside tumor nests, ${SUB_NAME[s]} touching ${partner}: ${groups[s].found.join("; ")}.` });
    // pooled fallback for the readouts that fail the per-patient rule (the nest has few T cells touching these cells)
    const pb = pooled(rows, `${pre}bleed_${src}`), pcaution = pb && pb.diff > BLEED_HIGH ? ` Caution: ${single} signal bleeds into touching T cells (+${pb.diff.toFixed(0)} pct-pts), so these contact results are unreliable.` : "";
    // The all-T sentence first. A subset gets its own sentence only when it names a marker that the all-T sentence does not
    // (not named there, or the other direction); otherwise the subsets are folded into one line.
    const sent = {};
    for (const g of ["all"].concat(SUBS)) {
      if (!groups[g].few.length) continue;
      const subject = g === "all" ? `T cells touching ${partner}` : `${SUB_NAME[g]} touching ${partner}`;
      sent[g] = pooledSentence("Inside tumor nests, ", subject, groups[g].few.map(([lab, base]) => [lab, pooled(rows, base)]), g === "all" && rows.length >= MIN_N, pcaution);
    }
    const strip = ({ sig, ...f }) => f, allSig = (sent.all && sent.all.sig) || new Map(), folded = [];
    if (sent.all) out.push({ key: `nest_${src}TouchPooled`, ...strip(sent.all) });
    for (const g of SUBS) {
      const sg = sent[g];
      if (!sg || !sg.sig) continue;   // not enough cells or patients: nothing to say
      if ([...sg.sig].some(([lab, sign]) => allSig.get(lab) !== sign)) out.push({ key: `nest_${src}TouchPooled_${g}`, ...strip(sg) });
      else folded.push(g);
    }
    if (folded.length) out.push({ key: `nest_${src}TouchPooled_subsets`, status: "info",
      text: `Inside tumor nests, ${folded.length === 1 ? `${folded[0]} subset touching ${partner} shows` : `${folded.slice(0, -1).join(", ")} and ${folded[folded.length - 1]} subsets touching ${partner} show`} the same picture.` });
    return out;
  }
  function nestSentences(rows) {
    const L = "Inside tumor nests: ", out = [];
    for (const [src, partner, control] of [["neut", "neutrophils", "nest_B_near_neut_log2"], ["mac", "macrophages", null]]) {
      const key = `nest_T_near_${src}_log2`;
      if (hasKey(rows, key)) out.push(proximity(rows, key, `T cells near ${partner}`, control && hasKey(rows, control) ? control : null, L));
      const quiet = [];
      for (const s of SUBS) {
        const k = `nest_${s}_near_${src}_log2`;
        if (!hasKey(rows, k)) continue;
        const d = direction(rows, k);
        if (d.dir === "down" || d.dir === "up") out.push(proximity(rows, k, `${SUB_NAME[s]} near ${partner}`, null, L)); else quiet.push(SUB_NAME[s]);
      }
      if (quiet.length) out.push({ key: `nest_sub_${src}`, status: "mixed", text: `${L}no consistent pattern (or too few patients) for ${quiet.join(", ")} near ${partner}.` });
    }
    const med = (k) => median(vals(rows, k));
    if (vals(rows, "comp_nest_pct_T").length >= MIN_N) {
      const subs = SUBS.filter((s) => vals(rows, `comp_nest_pct_${s}`).length >= MIN_N).map((s) => `${s} ${pctText(med(`comp_nest_pct_${s}`))}`);
      const parts = [`${pctText(med("comp_nest_pct_T"))} T cells${subs.length ? ` (${subs.join(", ")})` : ""}`];
      if (vals(rows, "comp_nest_pct_neut").length >= MIN_N) parts.push(`${pctText(med("comp_nest_pct_neut"))} neutrophils`);
      if (vals(rows, "comp_nest_pct_mac").length >= MIN_N) parts.push(`${pctText(med("comp_nest_pct_mac"))} macrophages`);
      const cmed = (k) => median(perPatientCounts(rows, "nest", k)), cn = (k) => perPatientCounts(rows, "nest", k).length >= MIN_N;
      let textCount;
      if (cn("T")) {
        const csubs = SUBS.filter((s) => cn(s)).map((s) => `${s} ${countText(cmed(s))}`);
        const cparts = [`${countText(cmed("T"))} T cells${csubs.length ? ` (${csubs.join(", ")})` : ""}`];
        if (cn("neut")) cparts.push(`${countText(cmed("neut"))} neutrophils`);
        if (cn("mac")) cparts.push(`${countText(cmed("mac"))} macrophages`);
        textCount = `Tumor nest: median ${cparts.join(", ")} per patient.`;
      }
      const text = `Tumor nest: median ${parts.join(", ")}.`;
      out.push({ key: "compNest", status: "info", text, textPct: text, ...(textCount ? { textCount } : {}) });
    }
    out.push(...nestContact(rows, "neut", "neutrophils", "neutrophil"), ...nestContact(rows, "mac", "macrophages", "macrophage"));
    for (const [src, partner] of [["neut", "neutrophil"], ["mac", "macrophage"]]) {
      const bins = GRAD.map((lab) => [lab, direction(rows, `grad_${src}_${lab}_log2`)]).filter(([, d]) => d.dir !== "few");
      if (!bins.length) continue;
      const [lab, d] = bins.reduce((a, b) => (Math.abs(b[1].med) > Math.abs(a[1].med) ? b : a));
      const nn = `${d.med < 0 ? d.down : d.up} of ${d.n} patients`;
      out.push({ key: `grad_${src}`, status: d.dir === "mixed" ? "mixed" : d.dir,
        text: `T cells in tumor nests by distance to the nearest ${partner}: the largest shift from chance is at ${binText(lab)}, ` +
          (d.dir === "mixed" ? `but patients disagree (${nn} ${d.med < 0 ? "below" : "above"} chance, median ${sgn(d.med)} log2).` : `${d.med < 0 ? "fewer" : "more"} T cells than chance in ${nn} (median ${sgn(d.med)} log2).`) });
    }
    return out;
  }
  // Macrophage markers next to T cells. prefix "TmacTouch" = all macrophages of the image, "nest_TmacTouch" = macrophages inside the nest.
  function macState(rows, prefix = "TmacTouch") {
    const nest = prefix.startsWith("nest_"), lead = nest ? "Inside tumor nests, macrophages touching T cells" : "Macrophages touching T cells";
    const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))].filter((k) => k.startsWith(prefix + "_") && k.endsWith("_diff"));
    if (!keys.length) return [];
    const found = [], mixed = [], few = [];
    for (const k of keys) {
      const name = k.slice(prefix.length + 1, -5), d = direction(rows, k);
      if (d.dir === "down" || d.dir === "up") found.push(`${MAC_STATE[name] || name} ${d.dir === "down" ? "lower" : "higher"} in ${d.dir === "down" ? d.down : d.up} of ${d.n} patients (median ${sgn(d.med, 1)} pct-pts)`);
      else if (d.dir === "mixed") mixed.push(name);
      else if (hasKey(rows, k.slice(0, -5) + "_near_sum")) few.push([MAC_STATE[name] || name, k.slice(0, -5)]);
    }
    const out = [];
    if (!found.length && !mixed.length) { if (!nest) out.push({ key: prefix, status: "na", text: `${lead}: too few patients to say.` }); }
    else out.push({ key: prefix, status: found.length ? "change" : "none",
      text: found.length ? `${lead}: ${found.join("; ")}.` : `${lead}: no consistent change in ${mixed.length} markers.` });
    if (few.length) {
      const sent = pooledSentence("", lead[0].toLowerCase() + lead.slice(1), few.map(([lab, base]) => [lab, pooled(rows, base)]), !found.length && !mixed.length && rows.length >= MIN_N);
      if (sent) { const { sig, ...f } = sent; out.push({ key: prefix + "Pooled", ...f }); }
    }
    return out;
  }

  // ---- under 5 patients: no per-patient rule is possible, so each original question is answered with the measured values ----
  // A log2 value counts as "fewer/more than chance" beyond ±0.3 and a contact difference is named beyond ±5 pct-pts (FEW_LOG2, FEW_PCT).
  // Status up/down/change only marks the sentence as a highlight; every sentence says it is descriptive.
  const FEW_LOG2 = 0.3, FEW_PCT = 5, FEW_MIN_CELLS = 20;
  function fewSummary(rows) {
    const n = rows.length, tag = n === 1 ? "1 patient, descriptive" : `${n} patients, median, descriptive`;
    const one = (k) => { const v = vals(rows, k); return v.length ? median(v) : null; };
    const word = (v) => (v <= -FEW_LOG2 ? "fewer than chance" : v >= FEW_LOG2 ? "more than chance" : "about chance");
    const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))], out = [];
    const prox = (key, what, ctrl, lead = "") => {
      const v = one(key); if (v === null) return;
      const W = lead ? what : what[0].toUpperCase() + what.slice(1);
      let text = `${lead}${W}: ${sgn(v)} log2, ${word(v)} (${tag}).`;
      const c = ctrl ? one(ctrl) : null;
      if (c !== null && Math.abs(v) >= FEW_LOG2) text += word(c) === word(v) ? ` B cells too (${sgn(c)}), so not specific to T cells.` : ` B cells ${sgn(c)} (${word(c)}).`;
      out.push({ key, status: v <= -FEW_LOG2 ? "down" : v >= FEW_LOG2 ? "up" : "none", text });
    };
    const subs = (src, partner, lead) => {
      const parts = SUBS.map((s) => [s, one(`nest_${s}_near_${src}_log2`)]).filter(([, v]) => v !== null);
      if (parts.length) out.push({ key: `nest_sub_${src}`, status: "info", text: `${lead}subsets near ${partner}: ${parts.map(([s, v]) => `${s} ${sgn(v)}`).join(", ")} log2 (${tag}).` });
    };
    // markers: [label, value, cells or null]; names the ones beyond ±FEW_PCT
    const markers = (key, subject, items, bleedV, single) => {
      items = items.filter(([, v]) => v !== null && Number.isFinite(v));
      if (!items.length) return;
      const cells = items.map(([, , c]) => c).filter((c) => c !== null), nc = cells.length ? Math.max(...cells) : null;
      const S = subject[0].toUpperCase() + subject.slice(1), cnt = nc !== null ? `${nc.toLocaleString("en-US")} cells, ` : "";
      if (nc !== null && nc < FEW_MIN_CELLS) { out.push({ key, status: "info", text: `${S}: only ${nc} cell${nc === 1 ? "" : "s"}, too few to compare.` }); return; }
      const big = items.filter(([, v]) => Math.abs(v) >= FEW_PCT), rest = items.length - big.length;
      const named = big.map(([lab, v]) => `${lab} ${v < 0 ? "lower" : "higher"} (${sgn(v, 1)} pct-pts)`);
      const tail = !big.length ? `all ${items.length} readouts within ±${FEW_PCT} pct-pts` : rest ? `${rest} other${rest === 1 ? "" : "s"} within ±${FEW_PCT}` : "";
      const caution = bleedV !== null && bleedV > BLEED_HIGH ? ` Caution: ${single} signal bleeds into the neighbouring cell (+${bleedV.toFixed(0)} pct-pts), so treat as unreliable.` : "";
      out.push({ key, status: big.length ? "change" : "none", text: `${S}: ${[...named, ...(tail ? [tail] : [])].join("; ")} (${cnt}${tag}).${caution}` });
    };
    const viaPool = (base) => { const p = pooled(rows, base); return p ? [p.diff, p.n_near] : [one(base + "_diff"), null]; };
    const listFor = (prefix, names) => keys.filter((k) => k.startsWith(prefix) && k.endsWith("_diff") && !k.includes("bleed") && !/_(CD8|CD4|Treg)_/.test(k.slice(prefix.length - 1)))
      .map((k) => { const nm = k.slice(prefix.length, -5), [v, c] = viaPool(k.slice(0, -5)); return [names[nm] || nm, v, c]; });
    const bleedOf = (base) => { const p = pooled(rows, base); return p ? p.diff : one(base + "_diff"); };

    // Project 1, question 1: do T cells keep away from neutrophils? (B cells as control)
    prox("T_near_neut_log2", "T cells near neutrophils", "B_near_neut_log2");
    prox("nest_T_near_neut_log2", "T cells near neutrophils", "nest_B_near_neut_log2", "Inside tumor nests, ");
    subs("neut", "neutrophils", "Inside tumor nests, T-cell ");
    // question 2: do T cells gather near macrophages?
    prox("T_near_mac_log2", "T cells near macrophages");
    prox("nest_T_near_mac_log2", "T cells near macrophages", null, "Inside tumor nests, ");
    subs("mac", "macrophages", "Inside tumor nests, T-cell ");
    // question 3: are T cells touching neutrophils / macrophages different (dying, exhausted, less active)?
    for (const [src, partner, single] of [["neut", "neutrophils", "neutrophil"], ["mac", "macrophages", "macrophage"]]) {
      markers(`${src}Touch`, `T cells touching ${partner}`, listFor(`${src}Touch_`, STATE), bleedOf(`${src}Touch_bleed_${src}`), single);
      markers(`nest_${src}Touch`, `inside tumor nests, T cells touching ${partner}`, listFor(`nest_${src}Touch_`, STATE), bleedOf(`nest_${src}Touch_bleed_${src}`), single);
    }
    if (!keys.some((k) => /Touch_apoptosis_diff$/.test(k))) out.push({ key: "death", status: "info", text: "T-cell death near neutrophils: no cell-death readout on this panel." });
    // question 4: do macrophages touching T cells look suppressive?
    markers("TmacTouch", "macrophages touching T cells", listFor("TmacTouch_", MAC_STATE), null, "");
    markers("nest_TmacTouch", "inside tumor nests, macrophages touching T cells", listFor("nest_TmacTouch_", MAC_STATE), null, "");
    // where the cells are
    const t = one("comp_nest_pct_T");
    if (t !== null) {
      const parts = [`${pctText(t)} T cells`], nu = one("comp_nest_pct_neut"), ma = one("comp_nest_pct_mac"), mT = one("comp_margin_pct_T");
      if (nu !== null) parts.push(`${pctText(nu)} neutrophils`); if (ma !== null) parts.push(`${pctText(ma)} macrophages`);
      const nn = one("comp_nest_cells"), size = nn !== null ? ` (${Math.round(nn).toLocaleString("en-US")} cells)` : "";
      const cnt = (k) => { const p = one(`comp_nest_pct_${k}`); return nn !== null && p !== null ? countText((p * nn) / 100) : null; };
      const cparts = [`${cnt("T")} T cells`]; if (cnt("neut")) cparts.push(`${cnt("neut")} neutrophils`); if (cnt("mac")) cparts.push(`${cnt("mac")} macrophages`);
      const ftext = nn !== null && nn < 200 ? `Tumor nest is tiny${size}, so nest results are unreliable (${tag}).` : `Tumor nest${size}: ${parts.join(", ")}${mT !== null ? `; T cells ${pctText(mT)} at the margin` : ""} (${tag}).`;
      const small = nn !== null && nn < 200;
      out.push({ key: "compNest", status: "info", text: ftext, ...(!small && nn !== null && cnt("T") ? { textPct: ftext, textCount: `Tumor nest${size}: ${cparts.join(", ")} (${tag}).` } : {}) });
    } else if (!keys.some((k) => /^(comp_|nest_)/.test(k))) out.push({ key: "compartments", status: "info", text: "No tumor nest found in this section, so nest results are not available." });
    return out;
  }

  // opts.loaded: the dataset came from the user's own file (explains why compartment results can be missing)
  // "111 of 145 patients" has a natural percent form, "77% of 145 patients". The count form stays in `text`; the percent form goes in `textPct`.
  const OF_PATIENTS = /(\d[\d,]*) of (\d[\d,]*) patients/g;
  function withPercentForms(list) {
    return list.map((f) => {
      if (f.textPct || !OF_PATIENTS.test(f.text)) return f;
      OF_PATIENTS.lastIndex = 0;
      return { ...f, textPct: f.text.replace(OF_PATIENTS, (_, a, b) => `${Math.round((100 * a.replace(/,/g, "")) / b.replace(/,/g, ""))}% of ${b} patients`) };
    });
  }
  function summarize(rows, opts = {}) {
    if (!rows || !rows.length) return [];
    if (rows.length < MIN_N) return fewSummary(rows);
    const out = [
      proximity(rows, "T_near_neut_log2", "T cells near neutrophils", "B_near_neut_log2"),
      proximity(rows, "T_near_mac_log2", "T cells near macrophages"),
      contact(rows, "neut", "neutrophils", "neutrophil"),
      contact(rows, "mac", "macrophages", "macrophage"),
    ];
    if (rows.some((r) => "T_near_neut_subset_log2" in r)) out.push(proximity(rows, "T_near_neut_subset_log2", "T cells near the neutrophil subset"));
    const keys = new Set(rows.flatMap((r) => Object.keys(r)));
    if ([...keys].some((k) => /^(comp_|nest_|grad_)/.test(k))) out.push(...nestSentences(rows));
    else out.push({ key: "compartments", status: "info", text: opts.results ? "This results table has no compartment results." : opts.loaded ? "No tumor cells were found in this dataset, so compartment results are not available." : "Compartment results not yet computed for this dataset." });
    out.push(...macState(rows, "TmacTouch"), ...macState(rows, "nest_TmacTouch"));
    return withPercentForms(out.filter(Boolean));
  }

  // which findings change when the pixel size is 20% smaller or larger (only used when it was estimated)
  function sensitivity(base, alternatives) {
    const changed = new Set();
    for (const alt of alternatives) {
      const a = summarize(alt), m = new Map(a.map((f) => [f.key, f.status]));
      for (const f of base) if (m.has(f.key) && m.get(f.key) !== f.status) changed.add(f.key);
    }
    return changed;
  }

  const api = { summarize, sensitivity, direction, pooled, perPatientCounts, countText };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.Findings = api;
})(typeof self !== "undefined" ? self : this);
