// The all-sky context map, in Mollweide.
//
// Mollweide is equal-area, which is the property that matters here: it shows
// honestly how little of the sky a 20 degree field covers, which a conformal
// projection would exaggerate near the poles.

/** A point on the sky, in degrees. */
export interface SkyPoint {
  raDeg: number;
  decDeg: number;
}

/** What to draw. */
export interface SkyScene {
  /** The true boresight, if truth is known. */
  truth?: SkyPoint;
  /** The estimated boresight, if there is a solution. */
  estimate?: SkyPoint;
  /** Half the diagonal field, degrees, for the footprint. */
  fieldRadiusDeg?: number;
}

/**
 * Mollweide, solving `2t + sin 2t = pi sin(dec)` for the auxiliary angle.
 *
 * Newton converges in a handful of steps everywhere except exactly at the
 * poles, where the derivative vanishes and the answer is just +/- pi/2.
 */
function project(raDeg: number, decDeg: number): { x: number; y: number } {
  // Centre the map on RA 180 so the 0/360 seam sits at the edges.
  let lambda = ((raDeg - 180) * Math.PI) / 180;
  const phi = (decDeg * Math.PI) / 180;
  if (lambda > Math.PI) lambda -= 2 * Math.PI;
  if (lambda < -Math.PI) lambda += 2 * Math.PI;

  let theta = phi;
  if (Math.abs(Math.abs(phi) - Math.PI / 2) < 1e-9) {
    theta = phi > 0 ? Math.PI / 2 : -Math.PI / 2;
  } else {
    const target = Math.PI * Math.sin(phi);
    for (let step = 0; step < 20; step += 1) {
      const f = 2 * theta + Math.sin(2 * theta) - target;
      const df = 2 + 2 * Math.cos(2 * theta);
      if (Math.abs(df) < 1e-12) break;
      const next = theta - f / df;
      if (Math.abs(next - theta) < 1e-12) {
        theta = next;
        break;
      }
      theta = next;
    }
  }
  // Unit-radius form: x spans [-2, 2] and y spans [-1, 1].
  return { x: (2 / Math.PI) * lambda * Math.cos(theta), y: Math.sin(theta) };
}

export class SkyView {
  private readonly canvas: HTMLCanvasElement;
  private readonly context: CanvasRenderingContext2D;

  constructor(canvas: HTMLCanvasElement) {
    const context = canvas.getContext("2d", { alpha: false });
    if (!context) throw new Error("2D canvas unavailable");
    this.canvas = canvas;
    this.context = context;
  }

  draw(scene: SkyScene): void {
    const { context, canvas } = this;
    const pad = 14;
    const halfWidth = (canvas.width - 2 * pad) / 4;
    const halfHeight = (canvas.height - 2 * pad) / 2;
    const scale = Math.min(halfWidth, halfHeight);
    const cx = canvas.width / 2;
    const cy = canvas.height / 2;
    const place = (p: SkyPoint) => {
      const { x, y } = project(p.raDeg, p.decDeg);
      return { x: cx + x * scale, y: cy - y * scale };
    };

    context.fillStyle = "#05070b";
    context.fillRect(0, 0, canvas.width, canvas.height);

    // Graticule: meridians every 30 degrees, parallels every 30.
    context.strokeStyle = "#1b2230";
    context.lineWidth = 1;
    for (let ra = 0; ra <= 360; ra += 30) {
      context.beginPath();
      for (let dec = -90; dec <= 90; dec += 2) {
        const p = place({ raDeg: ra, decDeg: dec });
        if (dec === -90) context.moveTo(p.x, p.y);
        else context.lineTo(p.x, p.y);
      }
      context.stroke();
    }
    for (let dec = -60; dec <= 60; dec += 30) {
      context.beginPath();
      for (let ra = 0; ra <= 360; ra += 2) {
        const p = place({ raDeg: ra, decDeg: dec });
        if (ra === 0) context.moveTo(p.x, p.y);
        else context.lineTo(p.x, p.y);
      }
      context.stroke();
    }
    // The outline itself.
    context.strokeStyle = "#2a3344";
    context.beginPath();
    for (let dec = -90; dec <= 90; dec += 1) {
      const p = place({ raDeg: 0.001, decDeg: dec });
      if (dec === -90) context.moveTo(p.x, p.y);
      else context.lineTo(p.x, p.y);
    }
    for (let dec = 90; dec >= -90; dec -= 1) {
      const p = place({ raDeg: 359.999, decDeg: dec });
      context.lineTo(p.x, p.y);
    }
    context.closePath();
    context.stroke();

    // The field footprint, as a small circle of the given radius about truth.
    const centre = scene.truth ?? scene.estimate;
    if (centre && scene.fieldRadiusDeg) {
      context.strokeStyle = "rgba(127, 209, 255, 0.55)";
      context.lineWidth = 1.2;
      context.beginPath();
      const radius = (scene.fieldRadiusDeg * Math.PI) / 180;
      const dec0 = (centre.decDeg * Math.PI) / 180;
      const ra0 = (centre.raDeg * Math.PI) / 180;
      let started = false;
      let previous = { x: 0, y: 0 };
      for (let step = 0; step <= 72; step += 1) {
        // Walk a small circle on the sphere, then project each point.
        const bearing = (step / 72) * 2 * Math.PI;
        const dec = Math.asin(
          Math.sin(dec0) * Math.cos(radius) +
            Math.cos(dec0) * Math.sin(radius) * Math.cos(bearing),
        );
        const ra =
          ra0 +
          Math.atan2(
            Math.sin(bearing) * Math.sin(radius) * Math.cos(dec0),
            Math.cos(radius) - Math.sin(dec0) * Math.sin(dec),
          );
        const p = place({
          raDeg: ((ra * 180) / Math.PI + 360) % 360,
          decDeg: (dec * 180) / Math.PI,
        });
        // A footprint straddling the seam would otherwise draw a line across
        // the whole map.
        if (started && Math.abs(p.x - previous.x) > scale) {
          context.stroke();
          context.beginPath();
          started = false;
        }
        if (!started) {
          context.moveTo(p.x, p.y);
          started = true;
        } else {
          context.lineTo(p.x, p.y);
        }
        previous = p;
      }
      context.stroke();
    }

    const dot = (p: SkyPoint, colour: string, radius: number) => {
      const at = place(p);
      context.fillStyle = colour;
      context.beginPath();
      context.arc(at.x, at.y, radius, 0, Math.PI * 2);
      context.fill();
    };
    if (scene.truth) dot(scene.truth, "#7fd1ff", 4);
    if (scene.estimate) dot(scene.estimate, "#ffa657", 2.6);

    context.fillStyle = "#76829a";
    context.font = "10px ui-monospace, monospace";
    context.fillText("RA 0", pad - 6, cy - 4);
    context.fillText("RA 360", canvas.width - pad - 34, cy - 4);
    context.fillText("+90", cx - 9, pad + 2);
    context.fillText("-90", cx - 9, canvas.height - pad + 2);
  }
}
