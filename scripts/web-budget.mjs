// Checks the in-browser wall clock for a 1000-trial run.
//
// Phase 9 requires the run to finish in Chrome inside 20 s and to score what
// the native CLI scores. A headless browser will not hold its page lifecycle
// open while work proceeds on workers, so with `?report=1` the page posts its
// result back to this origin and the harness waits for that.
//
//   cd web && npm run build
//   node scripts/web-budget.mjs [workers...]

import { createServer } from "node:http";
import { spawn } from "node:child_process";
import { readFile } from "node:fs/promises";
import { extname, join } from "node:path";

const ROOT = new URL("../web/dist/", import.meta.url);
const PORT = 5190;
const CHROME =
  process.env.CHROME ?? "C:/Program Files/Google/Chrome/Application/chrome.exe";
const BUDGET_SECONDS = 20;

const TYPES = {
  ".html": "text/html",
  ".js": "text/javascript",
  ".css": "text/css",
  ".wasm": "application/wasm",
  ".json": "application/json",
  ".svg": "image/svg+xml",
};

const asked = process.argv.slice(2).map(Number).filter((n) => n > 0);
const plan = asked.length > 0 ? asked : [1, 4, 8];

let finish;
const server = createServer(async (request, response) => {
  const url = new URL(request.url ?? "/", `http://localhost:${PORT}`);
  if (url.pathname === "/benchmark-done") {
    response.writeHead(204).end();
    finish?.(Object.fromEntries(url.searchParams));
    return;
  }
  const relative = url.pathname === "/" ? "index.html" : url.pathname.slice(1);
  try {
    const file = new URL(relative.replace(/^\/+/, ""), ROOT);
    const body = await readFile(file);
    response.writeHead(200, {
      "content-type": TYPES[extname(file.pathname)] ?? "application/octet-stream",
    });
    response.end(body);
  } catch {
    response.writeHead(404).end("not found");
  }
});

await new Promise((resolve) => server.listen(PORT, resolve));
console.log(`serving web/dist on :${PORT}, 1000 trials of nominal at seed 42`);
console.log("");
console.log(
  `  ${"workers".padStart(8)}${"seconds".padStart(10)}${"score".padStart(9)}` +
    `${"wrong".padStart(7)}${"median \u2033".padStart(11)}${"solve ms".padStart(10)}`,
);

let failures = 0;
const fastest = Math.max(...plan);

for (const workers of plan) {
  const query = new URLSearchParams({
    autorun: "1",
    report: "1",
    trials: "1000",
    seed: "42",
    preset: "nominal",
    workers: String(workers),
  });
  const chrome = spawn(
    CHROME,
    [
      "--headless",
      "--disable-gpu",
      "--no-sandbox",
      `--user-data-dir=${join(process.env.TEMP ?? ".", `st-budget-${workers}`)}`,
      `http://localhost:${PORT}/?${query.toString()}`,
    ],
    { stdio: "ignore" },
  );

  const result = await new Promise((resolve) => {
    finish = resolve;
    setTimeout(() => resolve(undefined), 300000);
  });
  chrome.kill();

  if (!result) {
    console.log(`  ${String(workers).padStart(8)}   timed out`);
    failures += 1;
    continue;
  }

  const seconds = Number(result.seconds);
  const over = workers === fastest && seconds > BUDGET_SECONDS;
  console.log(
    `  ${String(workers).padStart(8)}${seconds.toFixed(2).padStart(10)}` +
      `${(Number(result.score).toFixed(2) + "%").padStart(9)}` +
      `${String(result.wrong).padStart(7)}` +
      `${Number(result.median).toFixed(2).padStart(11)}` +
      `${Number(result.solve_ms).toFixed(3).padStart(10)}` +
      (over ? `   <-- over the ${BUDGET_SECONDS} s budget` : ""),
  );
  // The single-worker row is there to show what the parallelism buys; only the
  // widest configuration has to meet the budget.
  if (over) failures += 1;
  if (Number(result.wrong) !== 0) {
    console.log("           WRONG_CONFIDENT must be zero");
    failures += 1;
  }
}

server.close();
console.log("");
console.log(failures === 0 ? "web budget: PASS" : `web budget: FAIL (${failures})`);
process.exitCode = failures === 0 ? 0 : 1;
