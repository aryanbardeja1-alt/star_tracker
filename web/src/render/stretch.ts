// Turning a 16-bit sensor frame into something an eye can read.
//
// A linear ramp is useless here: the background sits around 20 ADU, a faint
// star peaks near 80, and a bright one saturates at 4095, so a linear map
// leaves everything but the brightest few stars indistinguishable from black.
// An asinh stretch is the usual astronomical answer — logarithmic in the
// highlights, linear through zero, so it keeps faint stars visible without
// clipping the bright ones.

/** How the frame is mapped to grey. */
export interface Stretch {
  /** Background level, in ADU. */
  black: number;
  /** ADU above the background that maps to mid grey. */
  softening: number;
  /** ADU above the background that maps to white. */
  white: number;
}

/**
 * Estimates a stretch from the frame itself.
 *
 * The black point is the median, which the stars cannot move because they
 * cover well under a percent of the frame. The white point is a high
 * percentile rather than the maximum, so one saturated star or a hot pixel
 * does not crush everything else into the floor.
 */
export function estimateStretch(image: Uint16Array): Stretch {
  // Sub-sample: a percentile does not need a million points.
  const step = Math.max(1, Math.floor(image.length / 40000));
  const sample: number[] = [];
  for (let at = 0; at < image.length; at += step) sample.push(image[at]!);
  sample.sort((a, b) => a - b);

  const at = (q: number) => sample[Math.min(sample.length - 1, Math.floor(q * sample.length))]!;
  const black = at(0.5);
  const noise = Math.max(1, at(0.84) - black);
  return {
    black,
    // A few noise widths: enough that noise stays dark and a faint star lifts.
    softening: noise * 6,
    white: Math.max(noise * 40, at(0.9995) - black),
  };
}

/**
 * Writes the frame into `out` as RGBA, applying the stretch.
 *
 * `out` must hold `4 * image.length` bytes. Grey is used rather than a colour
 * map: these are single-band intensities, and a false-colour ramp would imply
 * information that is not there.
 */
export function toRgba(image: Uint16Array, stretch: Stretch, out: Uint8ClampedArray): void {
  const { black, softening, white } = stretch;
  // asinh(x / s) normalised so the white point lands on 1.
  const top = Math.asinh(Math.max(white, 1) / softening);
  const gain = 255 / (top > 0 ? top : 1);

  for (let at = 0, rgba = 0; at < image.length; at += 1, rgba += 4) {
    const above = image[at]! - black;
    const grey = above <= 0 ? 0 : Math.asinh(above / softening) * gain;
    const value = grey > 255 ? 255 : grey;
    out[rgba] = value;
    out[rgba + 1] = value;
    out[rgba + 2] = value;
    out[rgba + 3] = 255;
  }
}
