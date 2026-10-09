// The benchmark charts, on plain canvas.
//
// Five small plots, one visual language: the same axis treatment, the same
// muted gridlines, the same palette, so they read as one instrument panel
// rather than five unrelated pictures. Outcome colour is consistent with the
// rest of the page — green correct, red wrong-confident, amber rejected, grey
// no-solution — because that mapping is load-bearing, not decorative.

import type { TrialRow } from "../wasm.js";

const INK = "#c9d2e3";
const DIM = "#76829a";
const GRID = "#1b2230";
const SERIES = "#7fd1ff";

export const OUTCOME_COLOUR: Record<TrialRow["outcome"], string> = {
  CORRECT: "#52d98b",
  WRONG_CONFIDENT: "#ff6b6b",
  REJECTED: "#ffc457",
  NO_SOLUTION: "#76829a",
};

interface Frame {
  context: CanvasRenderingContext2D;
  left: number;
  right: number;
  top: number;
  bottom: number;
}

/** Clears a canvas and returns the plotting box. */
function begin(canvas: HTMLCanvasElement, leftPad = 46): Frame {
  const context = canvas.getContext("2d");
  if (!context) throw new Error("2D canvas unavailable");
  context.fillStyle = "#05070b";
  context.fillRect(0, 0, canvas.width, canvas.height);
  context.font = "10px ui-monospace, monospace";
  context.textBaseline = "middle";
  return {
    context,
    left: leftPad,
    right: canvas.width - 12,
    top: 14,
    bottom: canvas.height - 26,
  };
}

function axes(frame: Frame, xLabel: string, yLabel: string): void {
  const { context, left, right, top, bottom } = frame;
  context.strokeStyle = GRID;
  context.lineWidth = 1;
  context.beginPath();
  context.moveTo(left, top);
  context.lineTo(left, bottom);
  context.lineTo(right, bottom);
  context.stroke();
  context.fillStyle = DIM;
  context.textAlign = "right";
  context.fillText(xLabel, right, bottom + 14);
  context.textAlign = "left";
  context.fillText(yLabel, left - 40, top - 4);
}

function noData(canvas: HTMLCanvasElement): void {
  const frame = begin(canvas);
  frame.context.fillStyle = DIM;
  frame.context.textAlign = "center";
  frame.context.fillText(
    "no data yet",
    (frame.left + frame.right) / 2,
    (frame.top + frame.bottom) / 2,
  );
}

/** Histogram of total attitude error, on a log x-axis. */
export function errorHistogram(canvas: HTMLCanvasElement, rows: TrialRow[]): void {
  const values = rows
    .filter((r) => Number.isFinite(r.total_arcsec) && r.total_arcsec > 0)
    .map((r) => r.total_arcsec);
  if (values.length === 0) return noData(canvas);

  // Decade-aligned bins: errors here span arcseconds to degrees.
  const lowest = Math.max(1e-3, Math.min(...values));
  const highest = Math.max(...values);
  const from = Math.floor(Math.log10(lowest));
  const to = Math.ceil(Math.log10(highest));
  const perDecade = 4;
  const bins = Math.max(1, (to - from) * perDecade);
  const counts = new Array<number>(bins).fill(0);
  for (const value of values) {
    const at = Math.min(bins - 1, Math.max(0, Math.floor((Math.log10(value) - from) * perDecade)));
    counts[at] = (counts[at] ?? 0) + 1;
  }

  const frame = begin(canvas);
  axes(frame, "arcsec", "trials");
  const { context, left, right, top, bottom } = frame;
  const peak = Math.max(...counts, 1);
  const width = (right - left) / bins;

  for (const [at, count] of counts.entries()) {
    if (count === 0) continue;
    const height = ((bottom - top) * count) / peak;
    context.fillStyle = SERIES;
    context.globalAlpha = 0.85;
    context.fillRect(left + at * width + 0.5, bottom - height, Math.max(1, width - 1), height);
  }
  context.globalAlpha = 1;

  // One tick per decade, which is as dense as this width can carry.
  context.fillStyle = DIM;
  context.textAlign = "center";
  for (let decade = from; decade <= to; decade += 1) {
    const x = left + (decade - from) * perDecade * width;
    context.fillText(`1e${decade}`, x, bottom + 14);
  }
  context.textAlign = "right";
  context.fillStyle = INK;
  context.fillText(String(peak), left - 4, top + 4);
}

/** Cross-boresight against roll, which shows how anisotropic the error is. */
export function crossRollScatter(canvas: HTMLCanvasElement, rows: TrialRow[]): void {
  const points = rows.filter(
    (r) => Number.isFinite(r.cross_arcsec) && Number.isFinite(r.roll_arcsec),
  );
  if (points.length === 0) return noData(canvas);

  const frame = begin(canvas);
  axes(frame, "roll ″", "cross ″");
  const { context, left, right, top, bottom } = frame;

  // Log-log: the spread covers several decades on both axes.
  const floor = 1e-2;
  const xs = points.map((p) => Math.max(floor, p.roll_arcsec));
  const ys = points.map((p) => Math.max(floor, p.cross_arcsec));
  const x0 = Math.floor(Math.log10(Math.min(...xs)));
  const x1 = Math.ceil(Math.log10(Math.max(...xs)));
  const y0 = Math.floor(Math.log10(Math.min(...ys)));
  const y1 = Math.ceil(Math.log10(Math.max(...ys)));
  const toX = (v: number) =>
    left + ((Math.log10(v) - x0) / Math.max(1e-9, x1 - x0)) * (right - left);
  const toY = (v: number) =>
    bottom - ((Math.log10(v) - y0) / Math.max(1e-9, y1 - y0)) * (bottom - top);

  for (const [at, point] of points.entries()) {
    context.fillStyle = OUTCOME_COLOUR[point.outcome];
    context.globalAlpha = 0.5;
    context.beginPath();
    context.arc(toX(xs[at]!), toY(ys[at]!), 1.9, 0, Math.PI * 2);
    context.fill();
  }
  context.globalAlpha = 1;

  context.fillStyle = DIM;
  context.textAlign = "center";
  for (let decade = x0; decade <= x1; decade += 1) {
    context.fillText(`1e${decade}`, toX(10 ** decade), bottom + 14);
  }
  context.textAlign = "right";
  for (let decade = y0; decade <= y1; decade += 1) {
    context.fillText(`1e${decade}`, left - 4, toY(10 ** decade));
  }
}

/** Outcome counts. */
export function outcomeBars(canvas: HTMLCanvasElement, rows: TrialRow[]): void {
  if (rows.length === 0) return noData(canvas);
  const order: TrialRow["outcome"][] = [
    "CORRECT",
    "WRONG_CONFIDENT",
    "REJECTED",
    "NO_SOLUTION",
  ];
  const counts = order.map((o) => rows.filter((r) => r.outcome === o).length);

  const frame = begin(canvas, 46);
  axes(frame, "", "trials");
  const { context, left, right, top, bottom } = frame;
  const peak = Math.max(...counts, 1);
  const slot = (right - left) / order.length;

  for (const [at, count] of counts.entries()) {
    const height = ((bottom - top) * count) / peak;
    const x = left + at * slot + slot * 0.18;
    const width = slot * 0.64;
    context.fillStyle = OUTCOME_COLOUR[order[at]!];
    context.fillRect(x, bottom - height, width, height);
    context.fillStyle = INK;
    context.textAlign = "center";
    if (count > 0) context.fillText(String(count), x + width / 2, bottom - height - 7);
    context.fillStyle = DIM;
    // Short labels: the full names do not fit at this width.
    const short = ["correct", "wrong", "rejected", "no sol."][at]!;
    context.fillText(short, x + width / 2, bottom + 14);
  }
}

/** Solve time per frame, rendering excluded. */
export function timingHistogram(canvas: HTMLCanvasElement, rows: TrialRow[]): void {
  const values = rows.map(
    (r) => r.centroid_ms + r.identify_ms + r.attitude_ms + r.verify_ms,
  );
  if (values.length === 0) return noData(canvas);

  const highest = Math.max(...values, 0.001);
  const bins = 32;
  const counts = new Array<number>(bins).fill(0);
  for (const value of values) {
    const at = Math.min(bins - 1, Math.floor((value / highest) * bins));
    counts[at] = (counts[at] ?? 0) + 1;
  }

  const frame = begin(canvas);
  axes(frame, "ms", "trials");
  const { context, left, right, top, bottom } = frame;
  const peak = Math.max(...counts, 1);
  const width = (right - left) / bins;
  context.fillStyle = SERIES;
  context.globalAlpha = 0.85;
  for (const [at, count] of counts.entries()) {
    if (count === 0) continue;
    const height = ((bottom - top) * count) / peak;
    context.fillRect(left + at * width + 0.5, bottom - height, Math.max(1, width - 1), height);
  }
  context.globalAlpha = 1;

  context.fillStyle = DIM;
  context.textAlign = "center";
  for (const fraction of [0, 0.5, 1]) {
    context.fillText(
      (highest * fraction).toFixed(2),
      left + fraction * (right - left),
      bottom + 14,
    );
  }
  context.textAlign = "right";
  context.fillStyle = INK;
  context.fillText(String(peak), left - 4, top + 4);
}

/** Success rate against how many stars the centroider found. */
export function successByDetected(canvas: HTMLCanvasElement, rows: TrialRow[]): void {
  if (rows.length === 0) return noData(canvas);

  // Bucket by detections, which run to a couple of hundred.
  const buckets = new Map<number, { total: number; correct: number }>();
  const step = 10;
  for (const row of rows) {
    const key = Math.floor(row.detected / step) * step;
    const bucket = buckets.get(key) ?? { total: 0, correct: 0 };
    bucket.total += 1;
    if (row.outcome === "CORRECT") bucket.correct += 1;
    buckets.set(key, bucket);
  }
  const keys = [...buckets.keys()].sort((a, b) => a - b);

  const frame = begin(canvas);
  axes(frame, "stars detected", "% correct");
  const { context, left, right, top, bottom } = frame;
  const lowest = keys[0] ?? 0;
  const highest = (keys[keys.length - 1] ?? 0) + step;
  const toX = (v: number) =>
    left + ((v - lowest) / Math.max(1, highest - lowest)) * (right - left);
  const toY = (fraction: number) => bottom - fraction * (bottom - top);

  context.strokeStyle = GRID;
  for (const fraction of [0.25, 0.5, 0.75, 1]) {
    context.beginPath();
    context.moveTo(left, toY(fraction));
    context.lineTo(right, toY(fraction));
    context.stroke();
    context.fillStyle = DIM;
    context.textAlign = "right";
    context.fillText(String(fraction * 100), left - 4, toY(fraction));
  }

  // Bar height is the rate; opacity carries how many trials back it up, so a
  // bucket holding three trials does not read as confidently as one holding
  // three hundred.
  const busiest = Math.max(...keys.map((k) => buckets.get(k)!.total), 1);
  const width = (right - left) / Math.max(1, keys.length);
  for (const [at, key] of keys.entries()) {
    const bucket = buckets.get(key)!;
    const rate = bucket.correct / bucket.total;
    context.globalAlpha = 0.25 + 0.75 * Math.sqrt(bucket.total / busiest);
    context.fillStyle = OUTCOME_COLOUR.CORRECT;
    const height = rate * (bottom - top);
    context.fillRect(left + at * width + 1, bottom - height, Math.max(1, width - 2), height);
    context.globalAlpha = 1;
  }

  context.fillStyle = DIM;
  context.textAlign = "center";
  for (const key of [keys[0], keys[Math.floor(keys.length / 2)], keys[keys.length - 1]]) {
    if (key === undefined) continue;
    context.fillText(String(key), toX(key + step / 2), bottom + 14);
  }
}

/** Draws every chart from the rows gathered so far. */
export function drawAll(rows: TrialRow[]): void {
  const at = (id: string) => document.getElementById(id) as HTMLCanvasElement | null;
  const hist = at("c-hist");
  const scatter = at("c-scatter");
  const outcomes = at("c-outcomes");
  const timing = at("c-timing");
  const detected = at("c-detected");
  if (hist) errorHistogram(hist, rows);
  if (scatter) crossRollScatter(scatter, rows);
  if (outcomes) outcomeBars(outcomes, rows);
  if (timing) timingHistogram(timing, rows);
  if (detected) successByDetected(detected, rows);
}
