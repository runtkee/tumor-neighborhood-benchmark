// End-to-end check of the built page in headless Chrome (Chrome DevTools Protocol, no extra packages).
// Usage: node tests/browser_smoke.js <path/to/page.html> [file-to-load ...]
// Loads the page; it runs its built-in example. Reports console errors, step statuses and rendered results.
const { spawn } = require("child_process");
const path = require("path"), os = require("os"), fs = require("fs");
const CHROME = process.env.CHROME || "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const page = path.resolve(process.argv[2]);
const uploads = process.argv.slice(3).map((p) => path.resolve(p));
const dir = fs.mkdtempSync(path.join(os.tmpdir(), "tnb-"));
const chrome = spawn(CHROME, ["--headless=new", "--disable-gpu", "--no-first-run", ...(process.getuid && process.getuid() === 0 ? ["--no-sandbox"] : []), `--user-data-dir=${dir}`, "--remote-debugging-port=0", "about:blank"], { stdio: ["ignore", "ignore", "pipe"] });
let buf = "";
chrome.stderr.on("data", async (d) => {
  buf += d; const m = buf.match(/DevTools listening on (ws:\/\/\S+)/); if (!m || chrome.started) return; chrome.started = true;
  const list = await (await fetch(m[1].replace("ws://", "http://").replace(/\/devtools\/browser\/.*/, "/json/list"))).json();
  const ws = new WebSocket(list.find((t) => t.type === "page").webSocketDebuggerUrl);
  let id = 0; const wait = new Map(); const errors = [];
  const send = (method, params = {}) => new Promise((res) => { const i = ++id; wait.set(i, res); ws.send(JSON.stringify({ id: i, method, params })); });
  ws.onmessage = (e) => { const msg = JSON.parse(e.data); if (msg.id && wait.has(msg.id)) { wait.get(msg.id)(msg.result || msg); wait.delete(msg.id); }
    if (msg.method === "Runtime.exceptionThrown") errors.push(msg.params.exceptionDetails.exception?.description || msg.params.exceptionDetails.text);
    if (msg.method === "Runtime.consoleAPICalled" && msg.params.type === "error") errors.push(msg.params.args.map((a) => a.value ?? a.description).join(" "));
    if (msg.method === "Log.entryAdded" && msg.params.entry.level === "error") errors.push(msg.params.entry.text); };
  ws.onopen = async () => {
    await send("Runtime.enable"); await send("Log.enable"); await send("Page.enable"); await send("DOM.enable");
    await send("Page.navigate", { url: "file://" + page });
    const evalJS = async (expr) => (await send("Runtime.evaluate", { expression: expr, returnByValue: true, awaitPromise: true })).result?.value;
    const until = async (expr, ms) => { const t0 = Date.now(); while (Date.now() - t0 < ms) { if (await evalJS(expr)) return true; await new Promise((r) => setTimeout(r, 250)); } return false; };
    const report = async (label) => {
      const r = await evalJS(`({datasets: [...document.querySelectorAll('#dsets .dset')].map(li=>li.innerText.replace(/\\s+/g,' ').slice(0,420)),
        claims: document.querySelectorAll('#claims .claim').length, plots: document.querySelectorAll('#claims svg').length,
        details: [...document.querySelectorAll('#details .detail')].map(d=>d.innerText.replace(/\\s+/g,' ').slice(0,300))})`);
      console.log(`\n== ${label}\n` + JSON.stringify(r, null, 1));
      return r;
    };
    const okStart = await until("document.querySelectorAll('#claims .claim').length > 0", 60000);
    await report("start (public reference) " + (okStart ? "rendered" : "DID NOT RENDER"));
    for (const f of uploads) {
      const before = await evalJS("document.querySelectorAll('#dsets .dset').length");
      const doc = await send("DOM.getDocument"); const node = await send("DOM.querySelector", { nodeId: doc.root.nodeId, selector: "#file" });
      await send("DOM.setFileInputFiles", { files: [f], nodeId: node.nodeId });
      await until(`document.querySelectorAll('#dsets .dset').length > ${before} && ![...document.querySelectorAll('#dsets .chip')].some(c=>/analysing/i.test(c.innerText))`, 180000);
      if (await evalJS("[...document.querySelectorAll('#dsets .chip')].some(c=>/needs pixel size/i.test(c.innerText))")) {
        await report("added " + path.basename(f) + " (asks for pixel size)");
        await evalJS("(() => { const f = [...document.querySelectorAll('form.pxform')].pop(); f.querySelector('input').value = '0.5'; f.requestSubmit(); })()");
        await until("![...document.querySelectorAll('#dsets .chip')].some(c=>/analysing|needs pixel size/i.test(c.innerText))", 180000);
        await report("after entering 0.5 µm per pixel");
      } else await report("added " + path.basename(f));
    }
    // machine-readable end state: one line per dataset in the list, and the rendered details
    const fin = await evalJS(`[...document.querySelectorAll('#dsets .dset')].map(li=>({name:(li.querySelector('.dname')||{}).innerText||'',meta:(li.querySelector('.dmeta')||{}).innerText||'',chips:[...li.querySelectorAll('.chip')].map(c=>c.innerText).join('|'),findings:li.querySelectorAll('.finds li').length,px:!!li.querySelector('form.pxform')})).concat([{details:[...document.querySelectorAll('#details .detail')].map(d=>d.innerText.replace(/\\s+/g,' ')),privacy:document.querySelector('.privacy').innerText,plots:document.querySelectorAll('#claims svg').length}])`);
    console.log("FINAL_STATE " + JSON.stringify(fin));
    const okExample = okStart;
    console.log("\nconsole errors:", errors.length ? errors : "none");
    if (process.env.SHOT) { await send("Emulation.setDeviceMetricsOverride", { width: 1440, height: Number(process.env.SHOT_H) || 1800, deviceScaleFactor: 1, mobile: false }); await new Promise((r) => setTimeout(r, 800)); const s = await send("Page.captureScreenshot", { format: "png" }); fs.writeFileSync(process.env.SHOT, Buffer.from(s.data, "base64")); }
    ws.close(); chrome.kill(); setTimeout(() => { try { fs.rmSync(dir, { recursive: true, force: true }); } catch (e) {} }, 500);
    process.exit(okExample && !errors.length ? 0 : 1);
  };
});
