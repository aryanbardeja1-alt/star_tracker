// The Single Pattern tab: simulate one frame, solve it, reveal the stages.

import { SensorView, type Stage } from "../render/sensor.js";
import { SkyView } from "../render/sky.js";
import { starLabel } from "./names.js";
import type {
  SelfCheckDto,
  SimConfig,
  TrialDetail,
  TrialRow,
  TruthStar,
} from "../wasm.js";
import type { Engine } from "../engine.js";

/** Milliseconds each reveal stage holds before the next. */
const STAGE_HOLD = 420;

const STAGES: Stage[] = ["image", "centroids", "pyramid", "matched", "verified"];

/** The fields the Advanced panel exposes, and how to show them. */
const ADVANCED: { key: keyof SimConfig; label: string; step: number }[] = [
  { key: "fov_deg", label: "FOV deg", step: 0.5 },
  { key: "psf_sigma_px", label: "PSF sigma px", step: 0.1 },
  { key: "exposure_s", label: "Exposure s", step: 0.01 },
  { key: "background_e", label: "Background e-", step: 10 },
  { key: "read_noise_e", label: "Read noise e-", step: 1 },
  { key: "sim_mag_limit", label: "Sim mag limit", step: 0.1 },
  { key: "db_mag_cut", label: "DB mag cut", step: 0.1 },
  { key: "mag_jitter", label: "Mag jitter", step: 0.01 },
  { key: "focal_error_frac", label: "Focal error", step: 0.001 },
  { key: "k1", label: "Distortion k1", step: 0.005 },
  { key: "dropout_prob", label: "Dropout", step: 0.01 },
  { key: "hot_pixels", label: "Hot pixels", step: 1 },
  { key: "max_centroids", label: "Max centroids", step: 1 },
  { key: "centroid_k_sigma", label: "Detect k sigma", step: 0.5 },
  { key: "centroid_sigma_px", label: "Assumed sigma px", step: 0.01 },
  { key: "id_k_sigma", label: "ID k sigma", step: 0.5 },
  { key: "calib_margin_arcsec", label: "Calib margin ″", step: 5 },
  { key: "verify_min_matches", label: "Min matches", step: 1 },
  { key: "verify_max_rms_px", label: "Max RMS px", step: 0.1 },
  { key: "err_threshold_arcsec", label: "Err threshold ″", step: 10 },
];

/** Reads an element by id, or throws: a missing node is a programming error. */
function need<T extends HTMLElement>(id: string): T {
  const element = document.getElementById(id);
  if (!element) throw new Error(`missing element #${id}`);
  return element as T;
}

export class PatternTab {
  private readonly engine: Engine;
  private readonly sensor: SensorView;
  private readonly sky: SkyView;
  private config: SimConfig;
  private timer: number | undefined;
  private labels = new Map<number, string>();
  private claimRight: boolean[] = [];

  constructor(engine: Engine, defaults: SimConfig) {
    this.engine = engine;
    this.config = defaults;
    this.sensor = new SensorView(need<HTMLCanvasElement>("sensor"));
    this.sky = new SkyView(need<HTMLCanvasElement>("sky"));

    this.buildAdvanced();
    need("generate").addEventListener("click", () => {
      // A fresh base seed, so "random" really does move off the last field.
      need<HTMLInputElement>("seed").value = String(
        Math.floor(Math.random() * 2 ** 48),
      );
      void this.run(true);
    });
    need("identify").addEventListener("click", () => void this.run(true));
    need<HTMLSelectElement>("preset").addEventListener("change", () => {
      void this.presetChanged();
    });
    for (const id of ["seed", "trial-index"]) {
      need(id).addEventListener("change", () => void this.run(false));
    }
  }

  /** The base seed and trial index currently shown. */
  current(): { seed: string; index: number } {
    return {
      seed: need<HTMLInputElement>("seed").value.trim() || "42",
      index: Number(need<HTMLInputElement>("trial-index").value) || 0,
    };
  }

  /** Opens a particular trial, as clicking a benchmark row does. */
  async open(seed: string, index: number, preset: string): Promise<void> {
    need<HTMLInputElement>("seed").value = seed;
    need<HTMLInputElement>("trial-index").value = String(index);
    const select = need<HTMLSelectElement>("preset");
    if (select.value !== preset) {
      select.value = preset;
      await this.presetChanged(false);
    }
    await this.run(true);
  }

  private async presetChanged(rerun = true): Promise<void> {
    const preset = need<HTMLSelectElement>("preset").value;
    this.config = await this.engine.config(preset);
    this.fillAdvanced();
    if (rerun) await this.run(true);
  }

  private buildAdvanced(): void {
    const host = need("advanced-fields");
    host.textContent = "";
    for (const field of ADVANCED) {
      const label = document.createElement("label");
      label.textContent = field.label;
      const input = document.createElement("input");
      input.type = "number";
      input.step = String(field.step);
      input.dataset.key = field.key;
      input.addEventListener("change", () => {
        const value = Number(input.value);
        if (Number.isFinite(value)) {
          (this.config as unknown as Record<string, number>)[field.key] = value;
          void this.run(true);
        }
      });
      label.append(input);
      host.append(label);
    }
    this.fillAdvanced();
  }

  private fillAdvanced(): void {
    for (const input of document.querySelectorAll<HTMLInputElement>(
      "#advanced-fields input",
    )) {
      const key = input.dataset.key as keyof SimConfig | undefined;
      if (!key) continue;
      const value = (this.config as unknown as Record<string, unknown>)[key];
      if (typeof value === "number") input.value = String(value);
    }
  }

  /** Simulates and solves the current seed, then reveals the stages. */
  async run(animate: boolean): Promise<void> {
    if (this.timer !== undefined) {
      window.clearTimeout(this.timer);
      this.timer = undefined;
    }
    const { seed, index } = this.current();

    // One round-trip to the worker: it scores the trial exactly as the
    // benchmark would and returns the frame it drew, so what is shown is what
    // was scored. The pixels are transferred rather than copied.
    const { detail, report, image, width, height } = await this.engine.frame(
      seed,
      index,
      this.config,
    );

    // Labels: an IAU name where the star has one, otherwise the HIP number.
    this.labels = new Map(
      report.matches.map((m) => [m.observed, starLabel(m)]),
    );

    // Was each claim right? Judged the same way the Rust ground-truth check
    // does: the star named must actually sit on that centroid.
    const rendered = detail.truth.filter((s: TruthStar) => !s.dropped);
    this.claimRight = [];
    for (const match of report.matches) {
      const centroid = report.centroids[match.observed];
      const star = rendered.find((s) => s.id === match.id);
      this.claimRight[match.observed] = Boolean(
        centroid &&
          star &&
          Math.hypot(star.x - centroid.x, star.y - centroid.y) <=
            this.config.truth_match_px,
      );
    }

    this.sensor.show({
      width,
      height,
      image,
      centroids: report.centroids,
      matches: report.matches,
      truth: detail.truth,
      falseStars: detail.false_stars,
      claimRight: this.claimRight,
      labels: this.labels,
    });

    this.sky.draw({
      truth: { raDeg: detail.true_attitude.ra_deg, decDeg: detail.true_attitude.dec_deg },
      estimate:
        detail.trial.outcome === "NO_SOLUTION"
          ? undefined
          : { raDeg: estimateRa(detail), decDeg: estimateDec(detail) },
      fieldRadiusDeg: halfDiagonalDeg(this.config),
    });

    this.showResult(detail.trial, detail, report.self_check);
    this.reveal(animate ? 0 : STAGES.length - 1, animate);
  }

  /** Walks the reveal stages, highlighting each chip as it lands. */
  private reveal(from: number, animate: boolean): void {
    const chips = [...document.querySelectorAll<HTMLElement>("#stages .chip")];
    const step = (at: number) => {
      this.sensor.setStage(STAGES[at]!);
      for (const [index, chip] of chips.entries()) {
        chip.classList.toggle("on", index <= at);
      }
      if (animate && at + 1 < STAGES.length) {
        this.timer = window.setTimeout(() => step(at + 1), STAGE_HOLD);
      }
    };
    step(from);
  }

  private showResult(
    row: TrialRow,
    detail: TrialDetail,
    check: SelfCheckDto | null,
  ): void {
    const badge = need("outcome-badge");
    badge.textContent = "";
    const span = document.createElement("span");
    span.className = `badge ${row.outcome.toLowerCase()}`;
    span.textContent = row.outcome.replace("_", " ");
    badge.append(span);

    const q = (v: [number, number, number, number]) =>
      v.map((x) => x.toFixed(6)).join("  ");
    const solved = row.outcome !== "NO_SOLUTION";
    const rows: [string, string][] = [
      ["Seed", row.seed],
      ["True quaternion", q(detail.true_attitude.quaternion)],
      ["Estimated quaternion", solved ? q(row.quaternion) : "—"],
      [
        "True RA / Dec / roll",
        `${detail.true_attitude.ra_deg.toFixed(4)}  ${detail.true_attitude.dec_deg.toFixed(4)}  ${detail.true_attitude.roll_deg.toFixed(4)}`,
      ],
      [
        "Estimated RA / Dec / roll",
        solved
          ? `${estimateRa(detail).toFixed(4)}  ${estimateDec(detail).toFixed(4)}  ${estimateRoll(detail).toFixed(4)}`
          : "—",
      ],
      ["Total error ″", solved ? row.total_arcsec.toFixed(3) : "—"],
      ["Cross-boresight ″", solved ? row.cross_arcsec.toFixed(3) : "—"],
      ["Roll ″", solved ? row.roll_arcsec.toFixed(3) : "—"],
      ["Stars detected", String(row.detected)],
      ["Fed to identifier", String(row.used)],
      ["Identified", String(row.claimed)],
      ["Correct", String(row.correct)],
      ["In database, visible", String(row.available)],
      [
        "Self-check",
        check
          ? `${check.confident ? "confident" : "refused"}, ${check.matched} matched, ${check.residual_rms_px.toFixed(3)} px RMS`
          : "—",
      ],
      ["Fitted focal scale", solved ? row.focal_scale.toFixed(6) : "—"],
      ["Pyramid tries", String(row.tries)],
    ];
    fillTable(need("result-table"), rows);

    fillTable(need("timing-table"), [
      ["Simulate", `${row.simulate_ms.toFixed(2)} ms`],
      ["Centroid", `${row.centroid_ms.toFixed(3)} ms`],
      ["Identify", `${row.identify_ms.toFixed(3)} ms`],
      ["Attitude", `${row.attitude_ms.toFixed(3)} ms`],
      ["Verify", `${row.verify_ms.toFixed(3)} ms`],
      [
        "Solve, render excluded",
        `${(row.centroid_ms + row.identify_ms + row.attitude_ms + row.verify_ms).toFixed(3)} ms`,
      ],
    ]);
  }
}

function fillTable(table: HTMLElement, rows: [string, string][]): void {
  table.textContent = "";
  const body = document.createElement("tbody");
  for (const [key, value] of rows) {
    const tr = document.createElement("tr");
    const th = document.createElement("th");
    th.textContent = key;
    const td = document.createElement("td");
    td.textContent = value;
    tr.append(th, td);
    body.append(tr);
  }
  table.append(body);
}

/** Half the diagonal field, degrees, from the configuration. */
function halfDiagonalDeg(config: SimConfig): number {
  const f = config.width / 2 / Math.tan((config.fov_deg * Math.PI) / 360);
  const corner = Math.hypot((config.width - 1) / 2, (config.height - 1) / 2);
  return (Math.atan(corner / f) * 180) / Math.PI;
}

// The estimated boresight comes from the trial's quaternion. Converting here
// rather than in Rust keeps `TrialRow` to the numbers the table needs.
function estimateMatrix(detail: TrialDetail): number[] {
  const [w, x, y, z] = detail.trial.quaternion;
  // Row-major rotation with b = A r.
  return [
    1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
    2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
    2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y),
  ];
}

function estimateRa(detail: TrialDetail): number {
  const a = estimateMatrix(detail);
  // The boresight's inertial direction is the third row.
  const ra = Math.atan2(a[7]!, a[6]!);
  return ((ra * 180) / Math.PI + 360) % 360;
}

function estimateDec(detail: TrialDetail): number {
  const a = estimateMatrix(detail);
  return (Math.atan2(a[8]!, Math.hypot(a[6]!, a[7]!)) * 180) / Math.PI;
}

function estimateRoll(detail: TrialDetail): number {
  const a = estimateMatrix(detail);
  return (Math.atan2(a[5]!, a[2]!) * 180) / Math.PI;
}
