// End-to-end check of the built page in headless Chrome (Chrome DevTools Protocol, no extra packages).
// Usage: node tests/browser_smoke.js <path/to/page.html> [file-to-load ...]
// Loads the page; it runs its built-in example. Reports console errors, step statuses and rendered results.
const { spawn } = require("child_process");
const path = require("path"), os = require("os"), fs = require("fs");
const CHROME = process.env.CHROME || "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome";
const page = path.resolve(process.argv[2]);
const uploads = process.argv.slice(3).map((p) => path.resolve(p));
const dir = fs.mkdtempSync(path.join(os.tmpdir(), "tnb-"));
const chrome = spawn(CHROME, ["--headless=new", "--disable-gpu", "--no-first-run", `--user-data-dir=${dir}`, "--remote-debugging-port=0", "about:blank"], { stdio: ["ignore", "ignore", "pipe"] });
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
      const r = await evalJS(`({load: document.querySelector('#load-status')?.innerText, run: document.querySelector('#run-status')?.innerText,
        checks: [...document.querySelectorAll('#checks tr')].map(t=>t.innerText.replace(/\\s+/g,' ')).slice(0,40),
        claims: document.querySelectorAll('#claims .claim').length, plots: document.querySelectorAll('#claims svg').length,
        chips: [...document.querySelectorAll('#claims .chip')].map(c=>c.innerText), summaryLines: (document.querySelector('#summary')?.value||'').split('\\n').length})`);
      console.log(`\n== ${label}\n` + JSON.stringify(r, null, 1));
      return r;
    };
    const okExample = await until("document.querySelectorAll('#claims .claim').length > 0", 90000);
    await report("example (built-in) " + (okExample ? "rendered" : "DID NOT RENDER"));
    for (const f of uploads) {  // load a user file through the real file input
      const doc = await send("DOM.getDocument"); const node = await send("DOM.querySelector", { nodeId: doc.root.nodeId, selector: "#file" });
      await send("DOM.setFileInputFiles", { files: [f], nodeId: node.nodeId });
      await until("document.querySelector('#load-status').innerText.match(/loaded|can.t use|stopped/)", 60000);
      await new Promise((r) => setTimeout(r, 1500));
      if (await evalJS("/can.t use/i.test(document.querySelector('#load-status').innerText)")) { await report("upload " + path.basename(f) + " (rejected at load)"); continue; }
      await evalJS("document.querySelector('#btn-check').click()");
      await until("document.querySelectorAll('#checks tr').length > 2 || document.querySelector('#run-status').innerText.match(/stopped/)", 120000);
      await report("upload " + path.basename(f));
    }
    console.log("\nconsole errors:", errors.length ? errors : "none");
    ws.close(); chrome.kill(); fs.rmSync(dir, { recursive: true, force: true });
    process.exit(okExample && !errors.length ? 0 : 1);
  };
});
