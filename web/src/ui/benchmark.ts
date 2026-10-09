// The Benchmark tab: run many trials across a pool of workers, score them,
// chart them, and let a failure be opened for inspection.

import { drawAll } from "../render/charts.js";
import { compiledModule, type FromWorker, type SimConfig, type TrialRow } from "../wasm.js";

/** Trials per message, so progress arrives while a run is in flight. */
const CHUNK = 25;

/** Rows per page of the trial table. */
const PAGE = 25;

/** Columns of the trial table: key, heading, and how to show a value. */
const COLUMNS: {
  key: keyof TrialRow | "solve_ms";
  head: string;
  show: (row: TrialRow) => string;
}[] = [
  { key: "index", head: "#", show: (r) => String(r.index) },
  { key: "outcome", head: "outcome", show: (r) => r.outcome.replace("_", " ") },
  { key: "total_arcsec", head: "total ″", show: (r) => fixed(r.total_arcsec, 3) },
  { key: "cross_arcsec", head: "cross ″", show: (r) => fixed(r.cross_arcsec, 3) },
  { key: "roll_arcsec", head: "roll ″", show: (r) => fixed(r.roll_arcsec, 3) },
  { key: "detected", head: "detected", show: (r) => String(r.detected) },
  { key: "claimed", head: "claimed", show: (r) => String(r.claimed) },
  { key: "correct", head: "correct", show: (r) => String(r.correct) },
  { key: "matched", head: "matched", show: (r) => String(r.matched) },
  {
    key: "residual_rms_px",
    head: "resid px",
    show: (r) => fixed(r.residual_rms_px, 4),
  },
  { key: "focal_scale", head: "focal", show: (r) => fixed(r.focal_scale, 6) },
  { key: "tries", head: "tries", show: (r) => String(r.tries) },
  { key: "solve_ms", head: "solve ms", show: (r) => fixed(solveMs(r), 3) },
];

function fixed(value: number, places: number): string {
  return Number.isFinite(value) ? value.toFixed(places) : "—";
}

function solveMs(row: TrialRow): number {
  return row.centroid_ms + row.identify_ms + row.attitude_ms + row.verify_ms;
}

function need<T extends HTMLElement>(id: string): T {
  const element = document.getElementById(id);
  if (!element) throw new Error(`missing element #${id}`);
  return element as T;
}

/** Median, 95th percentile and mean of a sample. */
function spread(values: number[]): { median: number; p95: number; mean: number } {
  if (values.length === 0) return { median: NaN, p95: NaN, mean: NaN };
  const sorted = [...values].sort((a, b) => a - b);
  const at = (q: number) =>
    sorted[Math.min(sorted.length - 1, Math.max(0, Math.ceil(q * sorted.length) - 1))]!;
  return {
    median: sorted[Math.floor(sorted.length / 2)]!,
    p95: at(0.95),
    mean: sorted.reduce((a, b) => a + b, 0) / sorted.length,
  };
}

/** One pool run, so a cancel can be told from a finish. */
interface Run {
  workers: Worker[];
  cancelled: boolean;
}

export class BenchmarkTab {
  private rows: TrialRow[] = [];
  private config: SimConfig;
  private preset = "nominal";
  private baseSeed = "42";
  private sortKey: string = "index";
  private sortDescending = false;
  private page = 0;
  private run: Run | undefined;
  /** Opens a trial in the other tab. */
  private readonly openTrial: (seed: string, index: number, preset: string) => void;
  /** Fetches a preset's configuration from the Rust side. */
  private readonly configFor: (preset: string) => Promise<SimConfig>;

  constructor(
    defaults: SimConfig,
    configFor: (preset: string) => Promise<SimConfig>,
    openTrial: (seed: string, index: number, preset: string) => void,
  ) {
    this.config = defaults;
    this.configFor = configFor;
    this.openTrial = openTrial;

    need("b-run").addEventListener("click", () => void this.start());
    need("b-cancel").addEventListener("click", () => this.cancel());
    need("b-compare").addEventListener("click", () => void this.comparePresets());
    need("b-csv").addEventListener("click", () => this.exportCsv());
    need("b-json").addEventListener("click", () => this.exportJson());
    need("b-prev").addEventListener("click", () => this.turnPage(-1));
    need("b-next").addEventListener("click", () => this.turnPage(1));

    // A sensible default for the worker count, which the user can override.
    const cores = navigator.hardwareConcurrency || 4;
    need<HTMLInputElement>("b-workers").value = String(Math.min(8, Math.max(1, cores - 1)));

    this.buildHead();
    this.render();
  }

  /**
   * Starts a run from URL parameters, for reproducing one from a link and for
   * checking the wall-clock budget without a hand on the mouse:
   * `?trials=1000&seed=42&preset=nominal&workers=4&autorun=1`.
   */
  async autorun(params: URLSearchParams): Promise<void> {
    const set = (id: string, value: string | null) => {
      if (value !== null) (need(id) as HTMLInputElement | HTMLSelectElement).value = value;
    };
    set("b-trials", params.get("trials"));
    set("b-seed", params.get("seed"));
    set("b-preset", params.get("preset"));
    set("b-workers", params.get("workers"));
    const started = performance.now();
    await this.start();

    // Report back to the origin when asked. A headless browser will not hold
    // its page lifecycle open for work done on workers, and the wall clock here
    // is the number that has to come in under 20 s, so a budget check needs the
    // page to tell it rather than the other way round.
    if (params.get("report")) {
      const rows = this.rows;
      const correct = rows.filter((r) => r.outcome === "CORRECT").length;
      const errors = rows
        .filter((r) => Number.isFinite(r.total_arcsec))
        .map((r) => r.total_arcsec);
      const query = new URLSearchParams({
        seconds: ((performance.now() - started) / 1000).toFixed(3),
        trials: String(rows.length),
        score: ((100 * correct) / Math.max(1, rows.length)).toFixed(4),
        wrong: String(rows.filter((r) => r.outcome === "WRONG_CONFIDENT").length),
        median: String(spread(errors).median),
        solve_ms: String(spread(rows.map(solveMs)).mean),
        workers: need<HTMLInputElement>("b-workers").value,
      });
      await fetch("/benchmark-done?" + query.toString()).catch(() => {});
    }
  }

  /** Runs `trials` trials across the pool, resolving when they are all in. */
  private async execute(
    trials: number,
    baseSeed: string,
    config: SimConfig,
    workerCount: number,
    onProgress: (done: number) => void,
  ): Promise<TrialRow[]> {
    const module = await compiledModule();
    const collected: TrialRow[] = [];
    const run: Run = { workers: [], cancelled: false };
    this.run = run;

    // Contiguous ranges, so each worker's indices are disjoint. Trials are
    // independent and seeded from their index, so the union is exactly what one
    // long run would produce.
    const count = Math.max(1, Math.min(workerCount, trials));
    const share = Math.ceil(trials / count);
    const ranges: [number, number][] = [];
    for (let start = 0; start < trials; start += share) {
      ranges.push([start, Math.min(share, trials - start)]);
    }

    await Promise.all(
      ranges.map(
        ([start, length]) =>
          new Promise<void>((resolve, reject) => {
            const worker = new Worker(new URL("../worker.ts", import.meta.url), {
              type: "module",
            });
            run.workers.push(worker);
            worker.onmessage = (event: MessageEvent<FromWorker>) => {
              const message = event.data;
              if (message.kind === "ready") {
                worker.postMessage({
                  kind: "run",
                  baseSeed,
                  start,
                  count: length,
                  chunk: CHUNK,
                  config,
                });
              } else if (message.kind === "rows") {
                collected.push(...message.rows);
                onProgress(collected.length);
              } else if (message.kind === "done") {
                worker.terminate();
                resolve();
              } else if (message.kind === "error") {
                worker.terminate();
                reject(new Error(message.message));
              }
            };
            worker.onerror = (event) => {
              worker.terminate();
              reject(new Error(event.message || "worker failed"));
            };
            worker.postMessage({ kind: "init", module });
          }),
      ),
    );

    this.run = undefined;
    collected.sort((a, b) => a.index - b.index);
    return collected;
  }

  async start(): Promise<void> {
    const trials = Math.max(1, Number(need<HTMLInputElement>("b-trials").value) || 1000);
    this.baseSeed = need<HTMLInputElement>("b-seed").value.trim() || "42";
    this.preset = need<HTMLSelectElement>("b-preset").value;
    const workers = Math.max(1, Number(need<HTMLInputElement>("b-workers").value) || 4);
    this.config = await this.configFor(this.preset);

    const progress = need<HTMLProgressElement>("b-progress");
    const status = need("b-status");
    need<HTMLButtonElement>("b-run").disabled = true;
    need<HTMLButtonElement>("b-cancel").disabled = false;
    need("b-compare-out").hidden = true;
    progress.max = trials;
    progress.value = 0;
    this.rows = [];
    this.page = 0;

    const started = performance.now();
    status.textContent = `Running ${trials} trials on ${workers} worker${workers === 1 ? "" : "s"}…`;

    try {
      const rows = await this.execute(trials, this.baseSeed, this.config, workers, (done) => {
        progress.value = done;
        const elapsed = (performance.now() - started) / 1000;
        status.textContent =
          `${done} / ${trials} trials · ${elapsed.toFixed(1)} s · ` +
          `${(done / Math.max(0.001, elapsed)).toFixed(0)} trials/s`;
      });
      this.rows = rows;
      const elapsed = (performance.now() - started) / 1000;
      status.textContent =
        `${rows.length} trials in ${elapsed.toFixed(2)} s on ${workers} worker${workers === 1 ? "" : "s"} ` +
        `(${(rows.length / Math.max(0.001, elapsed)).toFixed(0)} trials/s).`;
    } catch (error) {
      status.textContent = `Failed: ${error instanceof Error ? error.message : String(error)}`;
    } finally {
      need<HTMLButtonElement>("b-run").disabled = false;
      need<HTMLButtonElement>("b-cancel").disabled = true;
      progress.value = progress.max;
      this.render();
    }
  }

  private cancel(): void {
    const run = this.run;
    if (!run) return;
    run.cancelled = true;
    for (const worker of run.workers) worker.terminate();
    this.run = undefined;
    need("b-status").textContent = "Cancelled.";
    need<HTMLButtonElement>("b-run").disabled = false;
    need<HTMLButtonElement>("b-cancel").disabled = true;
  }

  /** Runs every preset at the current trial count and tabulates the scores. */
  private async comparePresets(): Promise<void> {
    const trials = Math.max(1, Number(need<HTMLInputElement>("b-trials").value) || 1000);
    const workers = Math.max(1, Number(need<HTMLInputElement>("b-workers").value) || 4);
    const seed = need<HTMLInputElement>("b-seed").value.trim() || "42";
    const status = need("b-status");
    need<HTMLButtonElement>("b-run").disabled = true;
    need<HTMLButtonElement>("b-compare").disabled = true;

    const results: { preset: string; rows: TrialRow[]; seconds: number }[] = [];
    try {
      for (const preset of ["easy", "nominal", "hard", "brutal"]) {
        status.textContent = `Comparing presets: ${preset}…`;
        const config = await this.configFor(preset);
        const started = performance.now();
        const rows = await this.execute(trials, seed, config, workers, () => {});
        results.push({ preset, rows, seconds: (performance.now() - started) / 1000 });
      }
      status.textContent = `Compared four presets at ${trials} trials each.`;
    } catch (error) {
      status.textContent = `Failed: ${error instanceof Error ? error.message : String(error)}`;
    } finally {
      need<HTMLButtonElement>("b-run").disabled = false;
      need<HTMLButtonElement>("b-compare").disabled = false;
    }

    const host = need("b-compare-out");
    host.hidden = results.length === 0;
    host.textContent = "";
    if (results.length === 0) return;
    const table = document.createElement("table");
    const head = document.createElement("tr");
    for (const text of ["preset", "score", "wrong", "rejected", "no sol.", "median ″", "seconds"]) {
      const th = document.createElement("th");
      th.textContent = text;
      head.append(th);
    }
    table.append(head);
    for (const { preset, rows, seconds } of results) {
      const correct = rows.filter((r) => r.outcome === "CORRECT").length;
      const errors = rows.filter((r) => Number.isFinite(r.total_arcsec)).map((r) => r.total_arcsec);
      const cells = [
        preset,
        `${((100 * correct) / Math.max(1, rows.length)).toFixed(1)}%`,
        String(rows.filter((r) => r.outcome === "WRONG_CONFIDENT").length),
        String(rows.filter((r) => r.outcome === "REJECTED").length),
        String(rows.filter((r) => r.outcome === "NO_SOLUTION").length),
        fixed(spread(errors).median, 2),
        seconds.toFixed(1),
      ];
      const tr = document.createElement("tr");
      for (const text of cells) {
        const td = document.createElement("td");
        td.textContent = text;
        tr.append(td);
      }
      table.append(tr);
    }
    host.append(table);
  }

  // --- presentation ---

  private buildHead(): void {
    const head = need("b-table").querySelector("thead");
    if (!head) return;
    head.textContent = "";
    const row = document.createElement("tr");
    for (const column of COLUMNS) {
      const th = document.createElement("th");
      th.textContent = column.head;
      th.addEventListener("click", () => {
        if (this.sortKey === column.key) this.sortDescending = !this.sortDescending;
        else {
          this.sortKey = column.key;
          this.sortDescending = false;
        }
        this.page = 0;
        this.render();
      });
      row.append(th);
    }
    head.append(row);
  }

  private sorted(): TrialRow[] {
    const key = this.sortKey;
    const value = (row: TrialRow): number | string =>
      key === "solve_ms"
        ? solveMs(row)
        : (row[key as keyof TrialRow] as number | string);
    const rows = [...this.rows].sort((a, b) => {
      const left = value(a);
      const right = value(b);
      if (typeof left === "string" || typeof right === "string") {
        return String(left).localeCompare(String(right));
      }
      // Not-a-number sorts last, so unsolved trials do not head the table.
      const safe = (v: number) => (Number.isFinite(v) ? v : Number.POSITIVE_INFINITY);
      return safe(left) - safe(right);
    });
    return this.sortDescending ? rows.reverse() : rows;
  }

  private turnPage(by: number): void {
    const pages = Math.max(1, Math.ceil(this.rows.length / PAGE));
    this.page = Math.min(pages - 1, Math.max(0, this.page + by));
    this.render();
  }

  private render(): void {
    const rows = this.rows;
    const correct = rows.filter((r) => r.outcome === "CORRECT").length;
    const wrong = rows.filter((r) => r.outcome === "WRONG_CONFIDENT").length;
    const errors = rows.filter((r) => Number.isFinite(r.total_arcsec)).map((r) => r.total_arcsec);
    const errorSpread = spread(errors);
    const claimed = rows.reduce((t, r) => t + r.claimed, 0);
    const right = rows.reduce((t, r) => t + r.correct, 0);
    const available = rows.reduce((t, r) => t + r.available, 0);

    const set = (id: string, text: string) => {
      need(id).textContent = text;
    };
    set("b-score", rows.length ? `${((100 * correct) / rows.length).toFixed(1)}%` : "—");
    set("b-wrong", rows.length ? String(wrong) : "—");
    set("b-rejected", rows.length ? String(rows.filter((r) => r.outcome === "REJECTED").length) : "—");
    set("b-nosol", rows.length ? String(rows.filter((r) => r.outcome === "NO_SOLUTION").length) : "—");
    set("b-precision", claimed ? (right / claimed).toFixed(4) : "—");
    set("b-recall", available ? (right / available).toFixed(4) : "—");
    set("b-median", rows.length ? fixed(errorSpread.median, 2) : "—");
    set("b-p95", rows.length ? fixed(errorSpread.p95, 2) : "—");
    set("b-solve", rows.length ? fixed(spread(rows.map(solveMs)).mean, 3) : "—");
    // The one number that must be zero gets a red tile when it is not.
    need("tile-wrong").classList.toggle("alarm", wrong > 0);

    drawAll(rows);

    const body = need("b-table").querySelector("tbody");
    if (body) {
      body.textContent = "";
      const sorted = this.sorted();
      const from = this.page * PAGE;
      for (const row of sorted.slice(from, from + PAGE)) {
        const tr = document.createElement("tr");
        tr.addEventListener("click", () => {
          this.openTrial(this.baseSeed, row.index, this.preset);
        });
        for (const column of COLUMNS) {
          const td = document.createElement("td");
          td.textContent = column.show(row);
          if (column.key === "outcome") td.className = `o-${row.outcome.toLowerCase()}`;
          tr.append(td);
        }
        body.append(tr);
      }
      const pages = Math.max(1, Math.ceil(sorted.length / PAGE));
      set("b-page", rows.length ? `page ${this.page + 1} of ${pages}` : "—");
    }
  }

  // --- export ---

  private download(name: string, text: string, mime: string): void {
    const url = URL.createObjectURL(new Blob([text], { type: mime }));
    const anchor = document.createElement("a");
    anchor.href = url;
    anchor.download = name;
    anchor.click();
    // Revoked on the next tick, once the click has been handled.
    window.setTimeout(() => URL.revokeObjectURL(url), 0);
  }

  private exportCsv(): void {
    if (this.rows.length === 0) return;
    const keys: (keyof TrialRow)[] = [
      "index", "seed", "outcome", "total_arcsec", "cross_arcsec", "roll_arcsec",
      "claimed", "correct", "available", "detected", "used", "matched",
      "residual_rms_px", "focal_scale", "tries",
      "simulate_ms", "centroid_ms", "identify_ms", "attitude_ms", "verify_ms",
    ];
    const lines = [keys.join(",")];
    for (const row of this.rows) {
      lines.push(
        keys
          .map((key) => {
            const value = row[key];
            return typeof value === "number" && !Number.isFinite(value) ? "" : String(value);
          })
          .join(","),
      );
    }
    this.download(
      `bench_${this.preset}_${this.baseSeed}.csv`,
      lines.join("\n"),
      "text/csv",
    );
  }

  private exportJson(): void {
    if (this.rows.length === 0) return;
    const rows = this.rows;
    const count = (outcome: TrialRow["outcome"]) =>
      rows.filter((r) => r.outcome === outcome).length;
    const errors = rows.filter((r) => Number.isFinite(r.total_arcsec)).map((r) => r.total_arcsec);
    const claimed = rows.reduce((t, r) => t + r.claimed, 0);
    const right = rows.reduce((t, r) => t + r.correct, 0);
    const available = rows.reduce((t, r) => t + r.available, 0);
    const summary = {
      preset: this.preset,
      base_seed: this.baseSeed,
      trials: rows.length,
      score_percent: (100 * count("CORRECT")) / rows.length,
      wrong_confident: count("WRONG_CONFIDENT"),
      outcomes: {
        CORRECT: count("CORRECT"),
        WRONG_CONFIDENT: count("WRONG_CONFIDENT"),
        REJECTED: count("REJECTED"),
        NO_SOLUTION: count("NO_SOLUTION"),
      },
      id_precision: claimed ? right / claimed : 0,
      id_recall: available ? right / available : 0,
      error_arcsec: { total: spread(errors) },
      solve_ms: spread(rows.map(solveMs)),
    };
    this.download(
      `bench_${this.preset}_${this.baseSeed}.json`,
      JSON.stringify(summary, null, 2),
      "application/json",
    );
  }
}
