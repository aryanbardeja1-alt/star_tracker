// The Tracking tab: run a slew, and compare the predicted-window search
// against the full lost-in-space search on the very same frames.

import type { Engine } from "../engine.js";
import type { SimConfig, TrackRun, TrackStep } from "../wasm.js";

/** Reads an element by id, or throws: a missing node is a programming error. */
function need<T extends HTMLElement>(id: string): T {
  const element = document.getElementById(id);
  if (!element) throw new Error(`missing element #${id}`);
  return element as T;
}

const OUTCOME_COLOUR: Record<string, string> = {
  CORRECT: "#4ade80",
  WRONG_CONFIDENT: "#f87171",
  REJECTED: "#fbbf24",
  NO_SOLUTION: "#94a3b8",
};

/** Axis and label colours, matching the other charts. */
const AXIS = "#334155";
const TEXT = "#94a3b8";

interface Series {
  values: number[];
  colour: string;
  label: string;
}

/**
 * Draws one or more series against frame index.
 *
 * `log` puts the y axis on a log scale, which is the only way to show a
 * tracked time and a lost-in-space time on one pair of axes.
 */
function lineChart(
  canvas: HTMLCanvasElement,
  series: Series[],
  unit: string,
  log: boolean,
): void {
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const { width, height } = canvas;
  ctx.clearRect(0, 0, width, height);

  const pad = { left: 56, right: 12, top: 14, bottom: 28 };
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;

  const finite = series.flatMap((s) => s.values.filter((v) => Number.isFinite(v) && (!log || v > 0)));
  if (finite.length === 0) return;
  const lowest = Math.min(...finite);
  const highest = Math.max(...finite);
  // A flat series would otherwise collapse to a zero-height axis.
  const lo = log ? Math.max(lowest / 2, 1e-6) : Math.min(0, lowest);
  const hi = highest > lo ? highest * (log ? 2 : 1.1) : lo + 1;

  const frames = Math.max(...series.map((s) => s.values.length));
  const x = (i: number) => pad.left + (frames <= 1 ? 0 : (i / (frames - 1)) * plotW);
  const y = (v: number) => {
    const t = log
      ? (Math.log10(v) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo))
      : (v - lo) / (hi - lo);
    return pad.top + plotH - t * plotH;
  };

  ctx.strokeStyle = AXIS;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, pad.top + plotH);
  ctx.lineTo(pad.left + plotW, pad.top + plotH);
  ctx.stroke();

  ctx.fillStyle = TEXT;
  ctx.font = "11px ui-monospace, monospace";
  ctx.textAlign = "right";
  for (const fraction of [0, 0.5, 1]) {
    const value = log
      ? 10 ** (Math.log10(lo) + fraction * (Math.log10(hi) - Math.log10(lo)))
      : lo + fraction * (hi - lo);
    const at = y(value);
    ctx.fillText(value < 10 ? value.toFixed(2) : value.toFixed(0), pad.left - 6, at + 4);
    ctx.strokeStyle = AXIS;
    ctx.globalAlpha = 0.4;
    ctx.beginPath();
    ctx.moveTo(pad.left, at);
    ctx.lineTo(pad.left + plotW, at);
    ctx.stroke();
    ctx.globalAlpha = 1;
  }
  ctx.textAlign = "center";
  ctx.fillText("frame", pad.left + plotW / 2, height - 8);
  ctx.save();
  ctx.translate(14, pad.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(unit, 0, 0);
  ctx.restore();

  for (const line of series) {
    ctx.strokeStyle = line.colour;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    let open = false;
    line.values.forEach((value, index) => {
      if (!Number.isFinite(value) || (log && value <= 0)) {
        open = false;
        return;
      }
      if (open) ctx.lineTo(x(index), y(value));
      else ctx.moveTo(x(index), y(value));
      open = true;
    });
    ctx.stroke();
  }

  // Legend, top right.
  ctx.textAlign = "left";
  series.forEach((line, at) => {
    const top = pad.top + 4 + at * 14;
    ctx.fillStyle = line.colour;
    ctx.fillRect(pad.left + plotW - 96, top, 10, 3);
    ctx.fillStyle = TEXT;
    ctx.fillText(line.label, pad.left + plotW - 82, top + 5);
  });
}

export class TrackingTab {
  private readonly engine: Engine;
  private readonly configFor: (preset: string) => Promise<SimConfig>;
  private running = false;

  constructor(engine: Engine, configFor: (preset: string) => Promise<SimConfig>) {
    this.engine = engine;
    this.configFor = configFor;
    need("t-run").addEventListener("click", () => void this.run());
  }

  /** Runs a slew and shows it. */
  async run(): Promise<void> {
    if (this.running) return;
    this.running = true;
    const status = need("t-status");
    const button = need<HTMLButtonElement>("t-run");
    button.disabled = true;
    status.textContent = "slewing…";

    try {
      const preset = need<HTMLSelectElement>("t-preset").value;
      const seed = need<HTMLInputElement>("t-seed").value.trim() || "42";
      const frames = Number(need<HTMLInputElement>("t-frames").value) || 60;
      const rate = Number(need<HTMLInputElement>("t-rate").value);

      const config = await this.configFor(preset);
      const run = await this.engine.tracking(seed, config, {
        frames,
        rate_deg_s: rate,
        interval_s: 0.1,
        axis: [1, 0, 0],
        search_margin_arcsec: 60,
        field_slack_deg: 1,
      });
      this.show(run);
      status.textContent = "";
    } catch (error) {
      status.textContent = error instanceof Error ? error.message : String(error);
    } finally {
      button.disabled = false;
      this.running = false;
    }
  }

  private show(run: TrackRun): void {
    need("t-held").textContent = `${run.correct} / ${run.frames}`;
    need("t-wrong").textContent = String(run.wrong_confident);
    need("t-tile-wrong").classList.toggle("alarm", run.wrong_confident > 0);
    need("t-reacq").textContent = String(run.reacquisitions);
    need("t-error").textContent = `${run.median_arcsec.toFixed(2)}″`;
    need("t-tracked").textContent = `${run.track_ms.toFixed(3)} ms`;
    need("t-lost").textContent = `${run.lost_ms.toFixed(3)} ms`;
    need("t-speedup").textContent = `${run.speedup.toFixed(1)}×`;

    need("t-note").textContent =
      `The slew turns ${run.step_deg.toFixed(3)}° between frames, and the tracker ` +
      `searches ${run.window_arcsec.toFixed(0)}″ around its prediction — the window ` +
      `does not grow with the rate, because the prediction absorbs it. The first ` +
      `two frames are searched in full: a rate needs two attitudes before it can ` +
      `be predicted, so tracking begins at the third. In the browser ` +
      `performance.now() is clamped to 0.1 ms, so a tracked frame often lands ` +
      `on that floor and the real speedup is larger than the one shown: the ` +
      `native CLI measures the same slew at 0.034 ms a frame.`;

    lineChart(
      need<HTMLCanvasElement>("t-timing"),
      [
        { values: run.steps.map((s) => s.lost_ms), colour: "#f59e0b", label: "whole sky" },
        { values: run.steps.map((s) => s.track_ms), colour: "#38bdf8", label: "tracked" },
      ],
      "ms",
      true,
    );
    lineChart(
      need<HTMLCanvasElement>("t-error-chart"),
      [{ values: run.steps.map((s) => s.total_arcsec), colour: "#4ade80", label: "total error" }],
      "arcsec",
      false,
    );

    this.fillTable(run.steps);
  }

  private fillTable(steps: TrackStep[]): void {
    const table = need<HTMLTableElement>("t-table");
    table.textContent = "";
    const head = table.insertRow();
    for (const label of [
      "frame",
      "outcome",
      "predicted",
      "detected",
      "claimed",
      "correct",
      "error ″",
      "tracked ms",
      "sky ms",
    ]) {
      const th = document.createElement("th");
      th.textContent = label;
      head.appendChild(th);
    }
    for (const step of steps) {
      const row = table.insertRow();
      const cells = [
        String(step.index),
        step.outcome.replace("_", " "),
        step.reacquired ? "reacquired" : step.aided ? "yes" : "full search",
        String(step.detected),
        String(step.claimed),
        String(step.correct),
        Number.isFinite(step.total_arcsec) ? step.total_arcsec.toFixed(2) : "—",
        step.aided ? step.track_ms.toFixed(3) : "—",
        step.lost_ms.toFixed(3),
      ];
      cells.forEach((text, at) => {
        const cell = row.insertCell();
        cell.textContent = text;
        if (at === 1) cell.style.color = OUTCOME_COLOUR[step.outcome] ?? TEXT;
      });
    }
  }
}
