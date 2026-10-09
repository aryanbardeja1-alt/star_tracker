// Loading the WASM module, once, and sharing it with the workers.
//
// `WebAssembly.compileStreaming` compiles the binary a single time; the
// resulting `Module` is structured-cloneable, so each worker gets the compiled
// code by `postMessage` and only has to instantiate it. Each instance still
// builds its own pair database, which is unavoidable: the 10 MB table lives in
// that instance's linear memory.

import { type Tracker, initSync, version } from "./pkg/tracker_wasm.js";
import wasmUrl from "./pkg/tracker_wasm_bg.wasm?url";

/** The compiled module, shared with every worker. */
let compiled: WebAssembly.Module | undefined;

/** Compiles the module once and returns it. */
export async function compiledModule(): Promise<WebAssembly.Module> {
  if (!compiled) {
    compiled = await WebAssembly.compileStreaming(fetch(wasmUrl));
  }
  return compiled;
}

/** Instantiates synchronously from an already-compiled module, inside a worker. */
export function trackerFromModule(module: WebAssembly.Module): Tracker {
  initSync({ module });
  // The glue's `Tracker` is available once `initSync` has run.
  return new TrackerCtor();
}

// Imported eagerly so the worker does not need a dynamic import after
// `initSync`, which would be a second module evaluation.
import { Tracker as TrackerCtor } from "./pkg/tracker_wasm.js";

export { version };

// --- the shapes the Rust side sends back ---------------------------------

/** An attitude, in the forms the UI shows. */
export interface Attitude {
  quaternion: [number, number, number, number];
  ra_deg: number;
  dec_deg: number;
  roll_deg: number;
  focal_scale: number;
}

/** One detected blob. */
export interface CentroidDto {
  x: number;
  y: number;
  flux: number;
  pixels: number;
}

/** One identification. */
export interface MatchDto {
  observed: number;
  catalog_index: number;
  id: number;
  /** IAU proper name, where the star has one. */
  name: string | null;
  from_pyramid: boolean;
}

/** The truth-blind self-check. */
export interface SelfCheckDto {
  matched: number;
  residual_rms_px: number;
  pyramid_stars: number;
  claims_explained: number;
  claims: number;
  confident: boolean;
}

/** What `solve_frame` returns. */
export interface SolveReport {
  detected: number;
  centroids: CentroidDto[];
  matches: MatchDto[];
  attitude: Attitude | null;
  self_check: SelfCheckDto | null;
  identify_tries: number;
  identify_ms: number;
}

/** One truth star. */
export interface TruthStar {
  id: number;
  x: number;
  y: number;
  mag: number;
  dropped: boolean;
}

/** Ground truth for a frame. */
export interface FrameTruth {
  attitude: Attitude;
  stars: TruthStar[];
  false_stars: [number, number][];
}

/** One scored trial. */
export interface TrialRow {
  index: number;
  seed: string;
  outcome: "CORRECT" | "WRONG_CONFIDENT" | "REJECTED" | "NO_SOLUTION";
  total_arcsec: number;
  cross_arcsec: number;
  roll_arcsec: number;
  claimed: number;
  correct: number;
  available: number;
  detected: number;
  used: number;
  matched: number;
  residual_rms_px: number;
  focal_scale: number;
  quaternion: [number, number, number, number];
  tries: number;
  simulate_ms: number;
  centroid_ms: number;
  identify_ms: number;
  attitude_ms: number;
  verify_ms: number;
}

/** A trial plus what a drawing of it needs. */
export interface TrialDetail {
  trial: TrialRow;
  truth: TruthStar[];
  false_stars: [number, number][];
  true_attitude: Attitude;
}

/** The solver's configuration, as the Rust `SimConfig`. */
export interface SimConfig {
  width: number;
  height: number;
  fov_deg: number;
  psf_sigma_px: number;
  bit_depth: number;
  f0_electrons_per_s: number;
  exposure_s: number;
  background_e: number;
  read_noise_e: number;
  gain_e_per_adu: number;
  sim_mag_limit: number;
  db_mag_cut: number;
  mag_jitter: number;
  focal_error_frac: number;
  k1: number;
  false_star_range: [number, number];
  dropout_prob: number;
  hot_pixels: number;
  bright_blob: boolean;
  max_centroids: number;
  background_stride: number;
  centroid_k_sigma: number;
  centroid_min_pixels: number;
  centroid_max_pixels: number;
  centroid_sigma_px: number;
  id_k_sigma: number;
  calib_margin_arcsec: number;
  id_focal_tolerance_frac: number;
  id_max_tries: number;
  verify_match_px: number;
  verify_min_matches: number;
  verify_max_rms_px: number;
  err_threshold_arcsec: number;
  truth_match_px: number;
}

/** How a slew is simulated and how wide the tracker searches. */
export interface TrackConfig {
  frames: number;
  rate_deg_s: number;
  interval_s: number;
  axis: [number, number, number];
  search_margin_arcsec: number;
  field_slack_deg: number;
}

/** One frame of a tracking sequence. */
export interface TrackStep {
  index: number;
  outcome: string;
  aided: boolean;
  reacquired: boolean;
  detected: number;
  claimed: number;
  correct: number;
  matched: number;
  total_arcsec: number;
  track_ms: number;
  lost_ms: number;
  ra_deg: number;
  dec_deg: number;
  true_ra_deg: number;
  true_dec_deg: number;
}

/** A whole tracking run. */
export interface TrackRun {
  steps: TrackStep[];
  frames: number;
  correct: number;
  wrong_confident: number;
  reacquisitions: number;
  track_ms: number;
  lost_ms: number;
  speedup: number;
  median_arcsec: number;
  step_deg: number;
  window_arcsec: number;
}

/** Messages the main thread sends a worker. */
export type ToWorker =
  | { kind: "init"; module: WebAssembly.Module }
  | { kind: "config"; id: number; preset: string }
  | {
      kind: "track";
      id: number;
      baseSeed: string;
      config: SimConfig | string;
      track: TrackConfig | null;
    }
  | {
      kind: "frame";
      id: number;
      baseSeed: string;
      index: number;
      config: SimConfig | string;
    }
  | {
      kind: "run";
      baseSeed: string;
      start: number;
      count: number;
      chunk: number;
      config: SimConfig | string;
    };

/** A simulated and solved frame, with the pixels the main thread draws. */
export interface FrameResult {
  detail: TrialDetail;
  report: SolveReport;
  image: Uint16Array;
  width: number;
  height: number;
}

/** Messages a worker sends back. */
export type FromWorker =
  | {
      kind: "ready";
      version: string;
      catalogStars: number;
      databasePairs: number;
    }
  | { kind: "config"; id: number; config: SimConfig }
  | { kind: "track"; id: number; run: TrackRun }
  | ({ kind: "frame"; id: number } & FrameResult)
  | { kind: "rows"; rows: TrialRow[] }
  | { kind: "done" }
  | { kind: "error"; id?: number; message: string };
