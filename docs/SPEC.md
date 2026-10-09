# Star Tracker — Specification

Lost-in-space star tracker: image → centroids → unit vectors → star identification → attitude (Wahba) → verification.
Core in Rust, compiled natively (CLI benchmark) and to WebAssembly (browser demo + 1000-trial scoring).
Build plan and phase status live in `docs/BUILD_PLAN.md`. Read it at the start of every session and update its checkboxes when a phase is done.

## Repository layout

```
star-tracker/
├── docs/SPEC.md
├── docs/BUILD_PLAN.md          # phases, acceptance criteria, progress log
├── data/
│   ├── raw/                    # downloaded catalogue (gitignored)
│   └── catalog.bin             # preprocessed binary catalogue (committed, small)
├── tools/build_catalog/        # Rust bin: raw catalogue → catalog.bin
├── crates/
│   ├── core/                   # pure library, no I/O, no wasm deps
│   │   └── src/{catalog,camera,simulate,centroid,pairdb,identify,attitude,verify,bench,math}.rs
│   ├── cli/                    # native benchmark runner
│   └── wasm/                   # wasm-bindgen wrapper around core
└── web/                        # Vite + TypeScript site
    └── src/{main.ts,worker.ts,render/*.ts,ui/*.ts,pkg/ (wasm-pack output, gitignored)}
```

## Commands

```bash
cargo test --workspace                        # all tests
cargo test -p tracker-core --release          # slow statistical tests run in release
cargo bench -p tracker-core                   # criterion benchmarks
cargo run -p tracker-cli --release -- bench --trials 1000 --preset nominal --seed 42
cargo run -p tracker-cli --release -- track --frames 60 --rate 5 --preset nominal --seed 42
cargo run -p build-catalog --release          # regenerate data/catalog.bin
wasm-pack build crates/wasm --target web --release --out-dir ../../web/src/pkg
cd web && npm install && npm run dev          # local site
cd web && npm run build                       # static build → web/dist
```

Run `cargo fmt` and `cargo clippy --workspace -- -D warnings` before declaring any phase done.

## Conventions (do not change without updating this file)

**Units.** All internal angles are in radians and use f64. Convert to degrees or arcseconds only at the UI/CLI boundary.

**Inertial frame.** ICRS/J2000. A catalogue star's unit vector is
`r = (cos δ cos α, cos δ sin α, sin δ)`.

**Camera frame.** This is the OpenCV convention and is right-handed.
- +z is the boresight.
- +x points along increasing column (right in the image).
- +y points along increasing row (down in the image).

**Pixels.** Pixel centres sit at integer coordinates (col, row). The principal point is `((W-1)/2, (H-1)/2)`.

**Projection** (inverse projection is the opposite direction):
- Camera vector → pixel: `x = cx + f·bx/bz`, `y = cy + f·by/bz`. Only valid when `bz > 0`.
- Pixel → camera vector: `b = normalize(x − cx, y − cy, f)`.

**Attitude.** A is a rotation matrix with `b = A r`, mapping inertial vectors into the camera frame.
- Quaternions are Hamilton convention, scalar-first `[w, x, y, z]`.
- Always normalise quaternions and keep `w ≥ 0`.

**Randomness.** All randomness comes from a seeded `rand_chacha::ChaCha8Rng`. The seed for trial i is `splitmix64(base_seed ^ i)`. Native and WASM must give bit-identical trial outcomes for the same seed. Never use thread_rng or OS entropy.

## Maths reference

**Angular separation.** Compute as `atan2(|a×b|, a·b)`, never as `acos(a·b)`, which loses precision for small angles.

**Centroid.** `x̄ = Σ Iᵢxᵢ / Σ Iᵢ` over background-subtracted pixels in an 8-connected blob.

**Pair database.**
1. Take every catalogue pair (i < j) with separation ≤ the FOV diagonal plus a margin.
2. Store it as `(angle: f32, i: u16, j: u16)`, sorted by angle.
3. Build a k-vector index (Mortari) for O(1) range lookup.

Search tolerance is `ε(θ) = k_sigma · σ_angle + calib_margin + focal_frac · θ`, for an observed
separation θ. Here σ_angle ≈ `√2 · σ_centroid_px / f`.

The term in θ is there because a focal length that is off by a fraction `d` stretches every
observed angle by about `d`, so the error in a pair's angle grows with the angle and no constant
tolerance can cover it. Measured on the `hard` preset: a 0.2% focal error puts a pair out by
**8 arcseconds per degree** of separation, which at 20° is 160″ against a 50″ constant.

Search the tight tolerance (`focal_frac = 0`) first and widen only if that finds nothing. A wider
window admits wrong pairs, and with them triangles a fourth star can confirm by coincidence, so
frames that never needed it must not pay for it.

**Identification.** Use the Pyramid algorithm (Mortari et al. 2004).
1. Find triangle candidates from three pair lookups whose catalogue indices close the loop.
2. Confirm each candidate with a fourth star, checking all three of its angles against the triangle.
3. Iterate star combinations in the Mortari order, so that one false star doesn't block every attempt.
4. After a confirmed pyramid, compute a provisional attitude.
5. Using that attitude, identify the remaining observed stars by nearest catalogue neighbour within tolerance.
6. A confirmed pyramid must then explain the frame: if that sweep adds no star at all, the pyramid
   was a coincidence, so discard it and resume the triple walk rather than returning it. A right
   pyramid sweeps up ten or so stars; every wrong one measured swept up exactly zero. Only skip
   this when there is no other observed star left to predict.

**Wahba / attitude.** Build `B = Σ wᵢ bᵢ rᵢᵀ`, with weights `wᵢ = 1/σᵢ²`.
- Take the SVD `B = U S Vᵀ`.
- Then `A = U · diag(1, 1, det U · det V) · Vᵀ`.

The diagonal correction is mandatory; without it you can get a reflection instead of a rotation. QUEST is optional, as a cross-check only.

**Uniform random attitude.** Draw 4 independent N(0,1) samples and normalise them into a quaternion.

**Attitude error.**
- Compute `ΔA = A_est · A_trueᵀ`, convert it to a rotation vector φ in the camera frame, and take the total error as `|φ|`.
- Cross-boresight error is `√(φx² + φy²)`.
- Roll error is `|φz|`.
- Report all of these in arcseconds.

**PSF rendering.** Use a Gaussian integrated over each pixel with erf, not sampled at pixel centres; sampling biases the centroids.
- Star flux is `F0 · 10^(−0.4 m) · t_exp`.
- Add Poisson shot noise, Gaussian read noise and a background level.
- Clip at full well.

## The "double check" (verification) — two independent layers

**Self-check (truth-blind).** This is what a real tracker would do.
1. Reproject catalogue stars with `A_est`.
2. Count the observed centroids matched within tolerance.
3. Compute the residual RMS.

A solution is CONFIDENT only if matched ≥ `min_matches`, residual RMS ≤ `max_rms`, and at least 4 stars came from the pyramid.

**Ground-truth check.** This uses the simulator's truth. The verifier must never call identification internals.
- Compare each claimed catalogue ID against the simulator's truth ID for that centroid.
- Compute the attitude error against `A_true`.

**Trial outcomes.**

| Outcome | Meaning |
|---|---|
| CORRECT | CONFIDENT, all claimed IDs are right, and attitude error is below `err_threshold` |
| WRONG_CONFIDENT | CONFIDENT but with any wrong ID or a large error. **This is the worst case and the target rate is 0.** |
| REJECTED | A solution was found but the self-check refused it |
| NO_SOLUTION | Identification found nothing |

**Score.** The primary score is the percentage of trials that are CORRECT. Always show WRONG_CONFIDENT separately and prominently.

## Avoiding an "inverse crime"

The simulator must not be a perfect mirror of the solver's assumptions:
- Render stars down to the simulation limiting magnitude (default 6.5). This is fainter than the database cut (default 6.0), so the image contains stars that aren't in the database.
- Jitter rendered magnitudes (σ = 0.1 mag).
- Optionally perturb the true focal length and add radial distortion that the solver doesn't know about (controlled by the preset).
- Inject false stars and random dropouts according to the preset.

## Default camera and presets

**Camera.** 1024×1024 pixels, 20° FOV (f ≈ 2904 px), PSF σ = 1.0 px, 12-bit output, database magnitude cut 6.0, keep the brightest 15 centroids.

The 15 is the *identification* budget, because that search is combinatorial in the number of
stars. The attitude fit, the sweep that follows it and both verification layers see every
detection: those stages are linear, and roll -- the weakly constrained axis in a narrow field --
improves with both the number of stars and their spread. The brightest are a prefix of the
flux-sorted detections, so an index into them is already an index into the whole list.

| Preset | False stars | Dropout | Read noise | Focal-length error | Distortion |
|---|---|---|---|---|---|
| easy | 0 | 0% | low | 0 | 0 |
| nominal | 0–1 | 5% | normal | 0 | 0 |
| hard | 0–3 | 10% | high | 0.2% | small k1 |
| brutal | 2–5 | 20% | very high | 0.5% | k1 + bright-object blob |

All parameters live in one serialisable `SimConfig` struct, shared by the CLI and the web.

## Performance targets

- Native, nominal preset: ≤ 5 ms per frame end-to-end (render excluded) on a laptop.
- WASM: within 3× of native.
- The 1000-trial web benchmark must finish in under 20 s, including rendering.
- Pair database build at start-up: ≤ 300 ms in WASM. Otherwise precompute it into `catalog.bin`.

## Coding rules

- `tracker-core` is a pure library: no I/O, no printing, no `unwrap()` or `expect()` outside tests. Return `Result` with a `thiserror` enum.
- No allocation in per-star hot loops. Reuse buffers through a `Workspace` struct passed by `&mut`.
- No `unsafe`.
- Approved dependencies: `nalgebra`, `rand`, `rand_chacha`, `rand_distr`, `thiserror`, `serde`, `bincode`, `wasm-bindgen`, `serde-wasm-bindgen`, `criterion`, `proptest`, `libm` (erf). Ask the user before adding anything else.
- Every public function gets a doc comment stating its frames and units.
- The web front end is vanilla TypeScript with a canvas. No UI framework unless the user asks.

## Testing rules

**Required tests:**
- Projection round-trip.
- Wahba recovers the exact attitude from noiseless vectors, and stays within the expected error under noise (proptest).
- `det(A) = +1` always.
- Pyramid achieves 100% identification on noiseless vector-level frames.
- Centroid bias is below 0.02 px on rendered single stars at random sub-pixel offsets.
- Native/WASM parity on 20 fixed seeds.
- A golden benchmark: seed 42, nominal preset, 1000 trials. Store its numbers in `docs/BUILD_PLAN.md`; regressions must be explained.

**Discipline:**
- Never weaken a test or a threshold to make it pass. If a target looks wrong, stop and explain why to the user.
- When a bug is found, first add a failing test that reproduces it, then fix it.

## Common pitfalls (check these first when something breaks)

1. Sign or frame mix-ups between image rows and +y, or using A versus Aᵀ. Re-read the Conventions section.
2. Missing the determinant correction in the SVD.
3. Using `acos` for small angles.
4. Forgetting that RA wraps around at 2π.
5. Pair tolerance too tight, so nothing is found, or too loose, so the search explodes combinatorially.
6. Indexing the pair database by visible-star index instead of catalogue index.
7. Native/WASM drift from f32/f64 mixing or nondeterministic iteration order (e.g. `HashMap` iteration).
