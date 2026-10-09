# Star Tracker

A lost-in-space star tracker: it takes a picture of the night sky and works out
which way the camera was pointing, with no prior knowledge of its orientation.
The core is Rust, compiled both natively and to WebAssembly, so the same code
that runs the benchmark runs in the browser.

**Pipeline:** image → centroids → unit vectors → star identification → attitude
→ verification.

## Results

A thousand simulated frames at the nominal preset, seed 42:

| | |
|---|---|
| Correct | **99.90%** |
| Confidently wrong | **0** |
| Median attitude error | **1.13″** (p95 4.01″) |
| Identification precision / recall | 1.0000 / 0.998 |
| Solve time, rendering excluded | 1.8 ms per frame |

Across the four difficulty presets, 500 trials each:

| Preset | False stars | Dropout | Focal error | Distortion | Correct | Confidently wrong |
|---|---|---|---|---|---|---|
| easy | 0 | 0% | 0 | 0 | 100% | 0 |
| nominal | 0–1 | 5% | 0 | 0 | 99.9% | 0 |
| hard | 0–3 | 10% | 0.2% | k₁ = 0.01 | 99.6% | 0 |
| brutal | 2–5 | 20% | 0.5% | k₁ = 0.03 | 8.0% | 0 |

`brutal` is deliberately past the design point — it exists to show the failure
mode is *refusal*, not a confident lie. The number that matters most in this
project is the confidently-wrong column, and its target is zero everywhere.

**Tracking mode.** Once an attitude is known, the next frame is a much smaller
problem. Identification runs 2.2–5.0× faster on `nominal` and **13.8× faster on
`hard`**, holding every frame of a 60-frame slew at rates up to 60 °/s.

## The double check

Every solution is judged twice, and the two layers are kept strictly apart.

The **self-check** is truth-blind — it is what a real tracker can do in flight.
It reprojects the catalogue through the estimated attitude, counts how many
observed centroids that explains, measures the residual RMS, and requires that
every star the solver named actually lands on the centroid it was named for.
Only then is a solution called CONFIDENT.

The **ground-truth check** uses the simulator's truth, which the solver never
sees, and compares the claimed catalogue IDs and the attitude against it.

Each trial ends as one of:

| Outcome | Meaning |
|---|---|
| `CORRECT` | confident, every ID right, attitude within threshold |
| `WRONG_CONFIDENT` | confident but wrong — the worst case, target rate zero |
| `REJECTED` | a solution was found and the self-check refused it |
| `NO_SOLUTION` | identification found nothing |

## How it works

**Catalogue.** Hipparcos (VizieR I/239), cut at V ≤ 6.5, proper motion
propagated from epoch J1991.25 to J2000.0 in vector form. Stars closer together
than the detector can separate are merged, because the catalogue must not claim
two stars where the sensor can only ever produce one centroid. 8,763 entries,
committed as a 284 KB binary.

**Simulation.** Stars are rendered with a Gaussian PSF integrated over each
pixel — sampled at pixel centres it would bias the centroids — plus Poisson shot
noise, read noise, background and saturation. The simulator is deliberately not
a mirror of the solver's assumptions: it renders fainter than the solver's
database, jitters magnitudes, and on the harder presets perturbs the focal
length and adds radial distortion the solver does not model.

**Identification.** The Pyramid algorithm (Mortari et al. 2004), over a
k-vector-indexed database of 849,829 star pairs. Three pair lookups propose a
triangle, a fourth star confirms it, and a confirmed pyramid must then go on to
explain the rest of the frame — a coincidence can match four stars, but it
cannot predict ten.

**Attitude.** Wahba's problem by SVD, `A = U diag(1, 1, det U · det V) Vᵀ`,
with Davenport's q-method as an independent cross-check. The focal length is
fitted from the data, because a camera whose focal length is 0.2% off puts a
star pair out by about 8 arcseconds per degree of separation.

## Running it

```bash
cargo test --workspace --release                 # 149 tests
cargo run -p tracker-cli --release -- bench --trials 1000 --preset nominal --seed 42
cargo run -p tracker-cli --release -- track --frames 60 --rate 5 --preset nominal
```

The website:

```bash
wasm-pack build crates/wasm --target web --release --out-dir ../../web/src/pkg
cd web && npm install && npm run dev
```

It has three tabs: a single pattern stepped through stage by stage, a
thousand-trial benchmark run across a pool of web workers, and tracking.

## Native and WebAssembly agree

The browser figures are only worth anything if they are the same figures. All
randomness comes from a seeded ChaCha8 generator — trial *i* uses
`splitmix64(base_seed ^ i)` and nothing else — so the two targets must produce
identical outcomes, and a parity check enforces it: 80 trials across all four
presets, every outcome identical, attitudes agreeing to within 1e-12 radians.

That same seeding is why the browser benchmark can be split across workers: a
trial depends on its index alone, so a pool over disjoint ranges returns exactly
what one long run would. A thousand trials complete in about 5 seconds.

## Layout

```
crates/core/     the library: no I/O, no panics, no unsafe
crates/cli/      native benchmark and tracking runner
crates/wasm/     wasm-bindgen wrapper
tools/           catalogue builder
web/             Vite + TypeScript site, canvas rendering, no framework
docs/SPEC.md         conventions, frames, units and the maths reference
docs/BUILD_PLAN.md   phase-by-phase record of what was built and measured
```

`docs/BUILD_PLAN.md` is worth a look if you want the reasoning rather than the
result — every acceptance criterion, the measurements behind it, and the dead
ends, including the ones that were wrong the first time.

## Data

Star positions and magnitudes are real measurements from the **Hipparcos**
catalogue (ESA, 1997), retrieved from VizieR at CDS Strasbourg. Proper names
come from the **IAU Working Group on Star Names** catalogue. The images are
entirely synthetic.

## Licence

MIT — see [LICENSE](LICENSE).
