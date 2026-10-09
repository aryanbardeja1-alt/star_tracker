// The sensor view: the frame, plus whatever the current reveal stage shows.
//
// The frame is painted once into an offscreen canvas at full resolution; zoom
// and pan then only redraw the visible sub-rectangle. That keeps panning cheap
// no matter how far in the view is zoomed, and means the million-pixel stretch
// runs once per frame rather than once per interaction.

import { estimateStretch, toRgba } from "./stretch.js";
import type { CentroidDto, MatchDto, TruthStar } from "../wasm.js";

/** How far through the reveal the view is. */
export type Stage = "image" | "centroids" | "pyramid" | "matched" | "verified";

/** Everything the view can draw. */
export interface SensorScene {
  width: number;
  height: number;
  image: Uint16Array;
  centroids: CentroidDto[];
  matches: MatchDto[];
  /** Truth, revealed only at the verified stage. */
  truth: TruthStar[];
  falseStars: [number, number][];
  /** Which claims were right. Empty until verification. */
  claimRight: boolean[];
  /** Label per centroid index: a name or a HIP number. */
  labels: Map<number, string>;
}

const COLOURS = {
  centroid: "#7fd1ff",
  pyramid: "#c792ea",
  matched: "#52d98b",
  wrong: "#ff6b6b",
  missed: "#ffc457",
  falseStar: "#ff8ad8",
};

export class SensorView {
  private readonly canvas: HTMLCanvasElement;
  private readonly context: CanvasRenderingContext2D;
  private readonly offscreen: HTMLCanvasElement;
  private readonly offscreenContext: CanvasRenderingContext2D;
  private scene: SensorScene | undefined;
  private stage: Stage = "image";
  /** Visible region in image pixels. */
  private view = { x: 0, y: 0, w: 1, h: 1 };
  private dragging: { x: number; y: number } | undefined;

  constructor(canvas: HTMLCanvasElement) {
    this.canvas = canvas;
    const context = canvas.getContext("2d", { alpha: false });
    const offscreen = document.createElement("canvas");
    const offscreenContext = offscreen.getContext("2d", { alpha: false });
    if (!context || !offscreenContext) throw new Error("2D canvas unavailable");
    this.context = context;
    this.offscreen = offscreen;
    this.offscreenContext = offscreenContext;
    this.attach();
  }

  /** Paints a new frame and resets the view to fit it. */
  show(scene: SensorScene): void {
    this.scene = scene;
    this.canvas.width = scene.width;
    this.canvas.height = scene.height;
    this.offscreen.width = scene.width;
    this.offscreen.height = scene.height;

    const rgba = new Uint8ClampedArray(scene.image.length * 4);
    toRgba(scene.image, estimateStretch(scene.image), rgba);
    this.offscreenContext.putImageData(
      new ImageData(rgba, scene.width, scene.height),
      0,
      0,
    );

    this.view = { x: 0, y: 0, w: scene.width, h: scene.height };
    this.draw();
  }

  /** Sets the reveal stage and repaints. */
  setStage(stage: Stage): void {
    this.stage = stage;
    this.draw();
  }

  /** Centres the view on a point at the given zoom, clamped to the frame. */
  focus(x: number, y: number, zoom: number): void {
    const scene = this.scene;
    if (!scene) return;
    const w = Math.max(16, scene.width / zoom);
    const h = Math.max(16, scene.height / zoom);
    this.view = {
      x: clamp(x - w / 2, 0, scene.width - w),
      y: clamp(y - h / 2, 0, scene.height - h),
      w,
      h,
    };
    this.draw();
  }

  private attach(): void {
    this.canvas.addEventListener(
      "wheel",
      (event) => {
        if (!this.scene) return;
        event.preventDefault();
        const scene = this.scene;
        const point = this.toImage(event);
        const factor = Math.exp(event.deltaY * 0.0015);
        const w = clamp(this.view.w * factor, 24, scene.width);
        const h = clamp(this.view.h * factor, 24, scene.height);
        // Keep the cursor over the same pixel while zooming.
        const fx = (point.x - this.view.x) / this.view.w;
        const fy = (point.y - this.view.y) / this.view.h;
        this.view = {
          x: clamp(point.x - fx * w, 0, scene.width - w),
          y: clamp(point.y - fy * h, 0, scene.height - h),
          w,
          h,
        };
        this.draw();
      },
      { passive: false },
    );

    this.canvas.addEventListener("pointerdown", (event) => {
      this.dragging = this.toImage(event);
      this.canvas.classList.add("dragging");
      this.canvas.setPointerCapture(event.pointerId);
    });
    this.canvas.addEventListener("pointermove", (event) => {
      if (!this.dragging || !this.scene) return;
      const scene = this.scene;
      const point = this.toImage(event);
      this.view.x = clamp(this.view.x - (point.x - this.dragging.x), 0, scene.width - this.view.w);
      this.view.y = clamp(this.view.y - (point.y - this.dragging.y), 0, scene.height - this.view.h);
      this.draw();
    });
    const release = (event: PointerEvent) => {
      this.dragging = undefined;
      this.canvas.classList.remove("dragging");
      if (this.canvas.hasPointerCapture(event.pointerId)) {
        this.canvas.releasePointerCapture(event.pointerId);
      }
    };
    this.canvas.addEventListener("pointerup", release);
    this.canvas.addEventListener("pointercancel", release);
    this.canvas.addEventListener("dblclick", () => {
      if (!this.scene) return;
      this.view = { x: 0, y: 0, w: this.scene.width, h: this.scene.height };
      this.draw();
    });
  }

  /** Pointer position in image pixels. */
  private toImage(event: { clientX: number; clientY: number }): { x: number; y: number } {
    const box = this.canvas.getBoundingClientRect();
    return {
      x: this.view.x + ((event.clientX - box.left) / box.width) * this.view.w,
      y: this.view.y + ((event.clientY - box.top) / box.height) * this.view.h,
    };
  }

  private draw(): void {
    const scene = this.scene;
    const { context, canvas } = this;
    context.fillStyle = "#05070b";
    context.fillRect(0, 0, canvas.width, canvas.height);
    if (!scene) return;

    const scale = canvas.width / this.view.w;
    context.imageSmoothingEnabled = this.view.w > canvas.width;
    context.drawImage(
      this.offscreen,
      this.view.x,
      this.view.y,
      this.view.w,
      this.view.h,
      0,
      0,
      canvas.width,
      canvas.height,
    );

    // Overlay geometry is in image pixels, so push the same transform rather
    // than converting every coordinate by hand.
    context.save();
    context.setTransform(scale, 0, 0, scale, -this.view.x * scale, -this.view.y * scale);
    // Keep strokes and text a constant size on screen.
    const px = 1 / scale;
    context.lineWidth = 1.4 * px;
    const font = `${Math.round(11 * px * 1000) / 1000}px ui-monospace, monospace`;

    const order: Stage[] = ["image", "centroids", "pyramid", "matched", "verified"];
    const reached = (stage: Stage) => order.indexOf(this.stage) >= order.indexOf(stage);

    if (reached("centroids")) {
      context.strokeStyle = COLOURS.centroid;
      for (const centroid of scene.centroids) {
        ring(context, centroid.x, centroid.y, 7 * px);
      }
    }

    if (reached("pyramid")) {
      context.strokeStyle = COLOURS.pyramid;
      context.lineWidth = 2.2 * px;
      const pyramid = scene.matches.filter((m) => m.from_pyramid);
      for (const m of pyramid) {
        const c = scene.centroids[m.observed];
        if (c) ring(context, c.x, c.y, 11 * px);
      }
      // Join them, so the four-star figure the confirmation rests on is visible.
      if (pyramid.length >= 2) {
        context.beginPath();
        for (const [at, m] of pyramid.entries()) {
          const c = scene.centroids[m.observed];
          if (!c) continue;
          if (at === 0) context.moveTo(c.x, c.y);
          else context.lineTo(c.x, c.y);
        }
        context.closePath();
        context.globalAlpha = 0.5;
        context.stroke();
        context.globalAlpha = 1;
      }
      context.lineWidth = 1.4 * px;
    }

    if (reached("matched")) {
      context.font = font;
      context.textBaseline = "middle";
      for (const m of scene.matches) {
        const c = scene.centroids[m.observed];
        if (!c) continue;
        const verified = reached("verified") && scene.claimRight.length > 0;
        const right = scene.claimRight[m.observed] !== false;
        context.strokeStyle = verified ? (right ? COLOURS.matched : COLOURS.wrong) : COLOURS.matched;
        ring(context, c.x, c.y, 9 * px);
        const mark = verified ? (right ? "✓ " : "✗ ") : "";
        context.fillStyle = context.strokeStyle;
        context.fillText(
          mark + (scene.labels.get(m.observed) ?? String(m.id)),
          c.x + 13 * px,
          c.y,
        );
      }
    }

    if (reached("verified")) {
      // Truth stars nothing claimed: dashed, so a miss is visibly a miss.
      context.strokeStyle = COLOURS.missed;
      context.setLineDash([4 * px, 3 * px]);
      const claimedNear = scene.matches
        .map((m) => scene.centroids[m.observed])
        .filter((c): c is CentroidDto => Boolean(c));
      for (const star of scene.truth) {
        if (star.dropped) continue;
        const near = claimedNear.some(
          (c) => Math.hypot(c.x - star.x, c.y - star.y) <= 2.0,
        );
        if (!near) ring(context, star.x, star.y, 10 * px);
      }
      context.setLineDash([]);

      // And the injected false stars, now that truth may be shown.
      context.strokeStyle = COLOURS.falseStar;
      for (const [x, y] of scene.falseStars) {
        cross(context, x, y, 8 * px);
      }
    }

    context.restore();
  }
}

function ring(context: CanvasRenderingContext2D, x: number, y: number, r: number): void {
  context.beginPath();
  context.arc(x, y, r, 0, Math.PI * 2);
  context.stroke();
}

function cross(context: CanvasRenderingContext2D, x: number, y: number, r: number): void {
  context.beginPath();
  context.moveTo(x - r, y - r);
  context.lineTo(x + r, y + r);
  context.moveTo(x + r, y - r);
  context.lineTo(x - r, y + r);
  context.stroke();
}

function clamp(value: number, low: number, high: number): number {
  return value < low ? low : value > high ? high : value;
}
