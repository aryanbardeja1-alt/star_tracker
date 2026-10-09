// Native/WASM parity check, and the WASM timing ratio.
//
// Phase 8 acceptance: identical outcomes and errors agreeing to 1e-12 against
// native on 20 fixed seeds, with WASM within 3x of native.
//
// Run from the repository root:
//
//   wasm-pack build crates/wasm --target nodejs --release --out-dir pkg-node
//   node scripts/parity.mjs
//
// The native side comes from `tracker-cli parity`, which prints one trial per
// line at full f64 precision; the `bench` CSV rounds its floats for reading and
// could not support a 1e-12 comparison.
//
// Angles are compared in **radians, absolutely**, which is the unit the whole
// project works in. Comparing the arcsecond figures relatively instead would be
// a far stricter and unreachable test: an attitude error is a small difference
// between two nearly equal rotations, roughly 6e-6 rad, so 1e-13 of absolute
// agreement in the attitude already shows up as 1e-8 of relative difference in
// the error derived from it. That is cancellation, not disagreement.

import { execFileSync } from "node:child_process";
import { createRequire } from "node:module";

const require = createRequire(import.meta.url);
const wasm = require("../crates/wasm/pkg-node/tracker_wasm.js");

/** Trials compared per preset. The criterion names 20 fixed seeds. */
const TRIALS = 20;
/**
 * Trials used for the timing ratio.
 *
 * More than the parity sample on purpose: over 20 trials the first-call costs
 * dominate and the ratio moves by nearly a factor of two run to run.
 */
const TIMING_TRIALS = 200;
/** Base seed; trial i runs on splitmix64(SEED ^ i). */
const SEED = 42;
/** Radians per arcsecond, for converting the WASM side back. */
const ARCSEC = 206264.80624709636;
/** Angles and quaternion components must agree to this, absolutely. */
const ANGLE_TOLERANCE = 1e-12;
/** The residual RMS is a pixel quantity; measured agreement is 1.3e-10 px. */
const PIXEL_TOLERANCE = 1e-9;
/** WASM may be this many times slower than native. */
const MAX_SLOWDOWN = 3.0;

const PRESETS = ["easy", "nominal", "hard", "brutal"];

/** Fields compared exactly. */
const EXACT = [
  "index",
  "seed",
  "outcome",
  "claimed",
  "correct",
  "available",
  "detected",
  "used",
  "matched",
];

/** Runs the native reference and parses its lines. */
function native(preset) {
  const stdout = execFileSync(
    "cargo",
    [
      "run", "-q", "-p", "tracker-cli", "--release", "--",
      "parity", "--trials", String(TRIALS), "--seed", String(SEED),
      "--preset", preset,
    ],
    { encoding: "utf8", maxBuffer: 64 * 1024 * 1024 },
  );
  return stdout.trim().split("\n").map((line) => JSON.parse(line));
}

/**
 * True when two numbers agree to `tolerance`, absolutely.
 *
 * A trial with no solution carries not-a-number errors and an infinite
 * residual. The native side prints both as JSON `null` while WASM hands over
 * the real NaN and Infinity, so "absent" has three spellings and any of them
 * matches any other.
 */
function agrees(a, b, tolerance) {
  const missing = (v) => v === null || !Number.isFinite(v);
  if (missing(a) || missing(b)) return missing(a) && missing(b);
  return Math.abs(a - b) <= tolerance;
}

let failures = 0;
const measured = [];

for (const preset of PRESETS) {
  const expected = native(preset);
  const tracker = wasm.init();
  const actual = tracker.run_benchmark_chunk(BigInt(SEED), 0, TRIALS, preset);

  if (actual.length !== expected.length) {
    console.log(`${preset}: FAIL - ${actual.length} rows against ${expected.length}`);
    failures += 1;
    continue;
  }

  let mismatches = 0;
  let worstAngle = 0;
  let worstPixel = 0;
  const outcomes = {};

  for (let at = 0; at < expected.length; at += 1) {
    const want = expected[at];
    const got = actual[at];
    outcomes[got.outcome] = (outcomes[got.outcome] ?? 0) + 1;

    for (const field of EXACT) {
      if (String(want[field]) !== String(got[field])) {
        console.log(
          `${preset} trial ${at}: ${field} is ${got[field]}, native says ${want[field]}`,
        );
        mismatches += 1;
      }
    }

    // Angles, in radians. The WASM side reports arcseconds, so convert back.
    const angles = [
      ["total", want.total_error, got.total_arcsec / ARCSEC],
      ["cross", want.cross_error, got.cross_arcsec / ARCSEC],
      ["roll", want.roll_error, got.roll_arcsec / ARCSEC],
      ["focal_scale", want.focal_scale, got.focal_scale],
      ["q0", want.quaternion[0], got.quaternion[0]],
      ["q1", want.quaternion[1], got.quaternion[1]],
      ["q2", want.quaternion[2], got.quaternion[2]],
      ["q3", want.quaternion[3], got.quaternion[3]],
    ];
    for (const [name, a, b] of angles) {
      if (!agrees(a, b, ANGLE_TOLERANCE)) {
        console.log(
          `${preset} trial ${at}: ${name} differs by ${Math.abs(a - b).toExponential(2)} ` +
            `(${b} against a native ${a})`,
        );
        mismatches += 1;
      } else if (Number.isFinite(a) && Number.isFinite(b)) {
        worstAngle = Math.max(worstAngle, Math.abs(a - b));
      }
    }

    if (!agrees(want.residual_rms_px, got.residual_rms_px, PIXEL_TOLERANCE)) {
      console.log(
        `${preset} trial ${at}: residual_rms_px differs by ` +
          `${Math.abs(want.residual_rms_px - got.residual_rms_px).toExponential(2)} px`,
      );
      mismatches += 1;
    } else if (
      Number.isFinite(want.residual_rms_px) &&
      Number.isFinite(got.residual_rms_px)
    ) {
      worstPixel = Math.max(
        worstPixel,
        Math.abs(want.residual_rms_px - got.residual_rms_px),
      );
    }
  }

  if (mismatches === 0) {
    console.log(
      `${preset}: PASS  ${TRIALS}/${TRIALS} outcomes identical, ` +
        `angles within ${worstAngle.toExponential(2)} rad, ` +
        `residual within ${worstPixel.toExponential(2)} px  | ` +
        `${JSON.stringify(outcomes)}`,
    );
  } else {
    console.log(`${preset}: FAIL  ${mismatches} mismatched fields`);
    failures += 1;
  }
  measured.push(preset);
}

console.log("");
console.log("timing, solve stages only (rendering is excluded from the budget):");
console.log(
  `  ${"preset".padEnd(9)}${"native ms".padStart(11)}${"wasm ms".padStart(11)}` +
    `${"ratio".padStart(9)}`,
);
for (const preset of measured) {
  // Warm the module, then time the solve stages over a longer run.
  const warm = wasm.init();
  warm.run_benchmark_chunk(BigInt(SEED), 0, 20, preset);
  const timed = warm.run_benchmark_chunk(BigInt(SEED), 0, TIMING_TRIALS, preset);
  const wasmSolve =
    timed.reduce(
      (total, row) =>
        total + row.centroid_ms + row.identify_ms + row.attitude_ms + row.verify_ms,
      0,
    ) / timed.length;

  // `bench` exits non-zero when a run contains WRONG_CONFIDENT, which the
  // harder presets do; the report on stdout is still what we want.
  let stdout;
  try {
    stdout = execFileSync(
      "cargo",
      [
        "run", "-q", "-p", "tracker-cli", "--release", "--",
        "bench", "--trials", String(TIMING_TRIALS), "--seed", String(SEED),
        "--preset", preset, "--no-files",
      ],
      {
        encoding: "utf8",
        maxBuffer: 16 * 1024 * 1024,
        // Its stderr carries the WRONG_CONFIDENT notice, which is expected here.
        stdio: ["ignore", "pipe", "ignore"],
      },
    );
  } catch (failed) {
    stdout = failed.stdout ?? "";
  }
  const line = stdout.split("\n").find((l) => l.includes("solve (no render)"));
  const nativeSolve = Number(line.trim().split(/\s+/).slice(-2)[0]);
  const ratio = wasmSolve / nativeSolve;
  const verdict = ratio <= MAX_SLOWDOWN ? "" : "  <-- over 3x";
  console.log(
    `  ${preset.padEnd(9)}${nativeSolve.toFixed(3).padStart(11)}` +
      `${wasmSolve.toFixed(3).padStart(11)}${(ratio.toFixed(2) + "x").padStart(9)}${verdict}`,
  );
  if (ratio > MAX_SLOWDOWN) failures += 1;
}

console.log("");
if (failures === 0) {
  console.log("parity and timing: PASS");
} else {
  console.log(`parity and timing: FAIL (${failures} problems)`);
  process.exitCode = 1;
}
