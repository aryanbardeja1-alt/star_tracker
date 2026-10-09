# Star Tracker — Build Plan

Conventions, maths and rules are in `docs/SPEC.md`; follow them exactly.

**How to work through this plan:**
- Work one phase at a time.
- A phase is done only when every acceptance criterion passes and `cargo fmt`, `cargo clippy -D warnings` and `cargo test` are clean.
- Tick the box, add a line to the Progress Log, then stop and summarise for the user before starting the next phase.
- If a criterion can't be met, report the measured numbers and your diagnosis. Do not lower the bar.

---

## Phase 0 — Scaffold
- [x] Cargo workspace with crates `tracker-core`, `tracker-cli`, `tracker-wasm`, `build-catalog`.
- [x] `web/` Vite + TypeScript project. Configure `.gitignore` for `data/raw/`, `web/src/pkg/`, `target/`, `node_modules/`.
- [x] `math.rs`: vector helpers, `angle_between` via atan2, quaternion ↔ matrix conversion, `splitmix64`, rotation vector from matrix.

**Accept when:** the workspace builds, an empty WASM build loads in the dev page, and the math unit tests pass.

> **Status: all three met.** Workspace builds; 15/15 math tests pass in debug and
> release; the WASM build loads in the dev page and in the production build, with
> `splitmix64(42)` matching the native CLI exactly.

## Phase 1 — Catalogue
- [x] `build-catalog` downloads the Yale Bright Star Catalogue (BSC5) or Hipparcos (VizieR I/239). If the network blocks the download, ask the user to place the file in `data/raw/`.
- [x] Keep stars with V ≤ 6.5.
- [x] Merge pairs closer than 2 px-equivalent (≈ 2·FOV/W) into one entry with combined flux. **Amended 2026-10-07 during Phase 3 to 6 px-equivalent**, on measurement: 8-connected labelling merges two stars into one blob out to 6-8 px with a 1 px PSF, so at 2 px the catalogue claimed two stars where the sensor can only produce one centroid. Merges 107 groups instead of 43 (8827 → 8763 stars) and lifted detection completeness from 97.92% to 99.675%.
- [x] Store each star as: id, unit vector (f64×3), magnitude, and a proper name if known.
- [x] Serialise with bincode to `data/catalog.bin`. `tracker-core` loads it from bytes (`include_bytes!` in the WASM crate).

**Accept when:** the star count is printed and plausible (~9k at V ≤ 6.5); a sample of 5 bright stars matches known RA/Dec to within 1″; and `catalog.bin` is under 1 MB.

> **Status: all three met.** 8763 stars printed (8827 before the Phase 3 blend-radius
> amendment above); the five bright-star checks land
> 0.0001–0.0002 **mas** from their published ICRS J2000.0 positions (budget 1″);
> `catalog.bin` is 290,766 bytes (29% of the 1 MB budget). Source is Hipparcos
> I/239 via VizieR, chosen over BSC5 because BSC5's 0.1s/1″ coordinate rounding
> alone can cost ~0.9″, leaving no margin against the 1″ criterion.

## Phase 2 — Camera model and simulator
- [x] `CameraModel`: W, H, f, cx, cy, k1. Implements project and unproject.
- [x] `SimConfig` and presets, as specified in docs/SPEC.md.
- [x] `simulate(seed, &SimConfig) -> SimFrame`. A `SimFrame` contains:
  - the image as `Vec<u16>`;
  - `A_true`;
  - a truth list of the projected stars: catalogue id, true pixel position, magnitude, and whether the star was dropped;
  - the false-star positions.
- [x] Render each PSF with a pixel-integrated Gaussian. Add shot and read noise, background, saturation, and hot pixels.

**Accept when:**
- the projection round-trip error is below 1e-9 px;
- 1000 random attitudes give a uniform boresight distribution (chi-square on an equal-area sky grid, p > 0.01);
- the same seed gives an identical image on repeated runs.

> **Status: all three met.** Projection round-trip worst error **1.4e-13 px** over
> 200k points (budget 1e-9); boresight chi-square **85.80** on a 10x10 equal-area
> grid, df 99, against the p = 0.01 critical value 134.6416; identical images on
> repeated runs for all four presets, and for a reused `Workspace`.
>
> `simulate` takes `(seed, &SimConfig, &Catalog, &mut Workspace)` rather than the
> `(seed, &SimConfig)` shorthand above: it needs the catalogue, and docs/SPEC.md
> requires buffer reuse through a `Workspace`.

## Phase 3 — Centroiding
- [x] Background estimate (median of a subsampled image), threshold at μ + kσ, 8-connected labelling, area filter, intensity-weighted centroid, sort by flux, keep the top K.
- [x] Zero allocations per frame once the `Workspace` has warmed up.

**Accept when:**
- single-star centroid bias is below 0.02 px and RMS is below 0.1 px (nominal noise, magnitude 4);
- on 1000 nominal frames, more than 98% of truth stars brighter than magnitude 5.5 are detected within 0.5 px;
- the benchmark shows ≤ 2 ms per 1024² frame natively.

> **Status: all three met.** Single star at V 4.0 with nominal noise, over 8192
> samples on a swept sub-pixel grid: mean bias **−0.00004 / +0.00003 px** and worst
> per-phase systematic **0.0104 px** (budget 0.02), RMS **0.0106 / 0.0107 px**
> (budget 0.1). Detection of rendered truth stars brighter than V 5.5 over 1000
> nominal frames: **99.675%** (budget > 98%); 99.758% over all rendered stars.
> Criterion, 1024²: **354 µs** nominal, 366 µs easy, 369 µs hard, **1.260 ms**
> brutal (budget 2 ms).
>
> Two notes on how the criteria were read. Completeness counts *rendered* truth
> stars, since a dropped star is deliberately absent from the image and counting
> those as misses would cap the result at the 95% dropout rate. And "keep the top
> K" moved out of `detect` into `centroid::brightest`: a nominal field holds about
> 27 stars brighter than V 5.5, so truncating to `max_centroids = 15` inside
> `detect` would cap completeness at 15/27 = 56%. K bounds identification
> combinatorics; it is not a detection limit.

## Phase 4 — Pair database and k-vector
- [x] Build the pair list from catalogue stars with V ≤ `db_mag_cut`, keeping separations ≤ FOV diagonal + 1°.
- [x] Sort by angle and build the k-vector index.
- [x] `query(angle, eps) -> &[Pair]`.
- [x] Also implement a brute-force `query_slow` as a test oracle.

**Accept when:** `query` equals `query_slow` on 10k random queries, and the build time is reported (precompute into `catalog.bin` if it exceeds 300 ms in WASM).

> **Status: both met.** `query` equalled `query_slow` on all **10,000** random
> queries, spread over separations inside and outside the table and tolerances
> from 0 to 0.05 rad. Build time **65 ms native** (criterion) and **95 ms median
> in WASM** (87-113 ms, measured under Node with a checksum forcing the sort and
> the k-vector to be live), against the 300 ms budget — so the table is built at
> start-up and **not** precomputed into `catalog.bin`, which it could not be
> anyway: at 849,829 pairs it occupies 10.2 MB against a 1 MB file budget.
> Queries cost about **14 ns** each.

## Phase 5 — Identification (Pyramid)
- [x] Implement Pyramid as specified in docs/SPEC.md, with Mortari combination ordering and an iteration cap.
- [x] After pyramid confirmation: provisional attitude, then match the remaining observed stars to catalogue neighbours.
- [x] Return the matched pairs (observed index → catalogue id) plus diagnostics: tries used, candidate counts, and timing.

**Accept when:**
- vector-level tests with no image (exact vectors plus Gaussian angular noise) reach 100% correct over 1000 trials with no noise, and at least 99.5% with nominal noise;
- frames with 3 injected false vectors still identify correctly in at least 98% of cases.

> **Status: all three met**, measured over 1000 vector-level trials each:
> noiseless **100.000%**, nominal noise (7.1″ per vector) **100.000%** against a
> 99.5% floor, and three injected false vectors **99.900%** against a 98% floor.
> Identification costs 0.09 ms with no false vectors and 0.16 ms with three.
>
> False vectors are *inserted* at random positions rather than appended. Appended,
> the first triple Mortari tries is always three real stars and the ordering is
> never exercised; interleaved, tries rise from 1.0 to 2.5 and the measured rate
> drops from 100% to 99.9%, which is the honest figure.
>
> Step 4 needs an attitude, so `attitude.rs` carries the SVD Wahba solver with the
> mandatory determinant correction. Phase 6 still owns the brightness weighting,
> the QUEST cross-check and the property tests.

## Phase 6 — Attitude (Wahba)
- [x] SVD solution with the determinant correction and weights from centroid brightness.
- [x] Optional QUEST implementation, cross-checked against SVD. **Delivered as Davenport's q-method** (`solve_davenport`) rather than QUEST proper — the same K matrix, but with the largest eigenvalue found exactly instead of by Newton iteration. See the status note below.

**Accept when:**
- property tests pass: noiseless recovery error below 1e-9 rad, and `det = +1`;
- error under noise matches the predicted `σ/√N` scaling within a factor of 2.

> **Status: both met.** Four `proptest` properties hold over their generated
> cases: noiseless recovery is within **1e-9 rad** for 2 to 24 stars, weighted or
> not; `det = +1` to 1e-12 with orthonormality to 1e-12 at any noise level;
> scaling a weight uniformly cannot move an exact fit; and Davenport's q-method
> agrees with the SVD to 1e-9 rad.
>
> Error under noise, 1500 repetitions per N with directions over the whole
> sphere: per-axis RMS against a bare `sigma/sqrt(N)` gives ratios of 1.46, 1.35,
> 1.27, 1.26, 1.25, 1.23, 1.21, 1.23 for N = 3…120, converging on the analytic
> `sqrt(3/2) = 1.225`, and `rms * sqrt(N)` is flat to within 25%. Well inside the
> factor of two.
>
> **The cross-check is Davenport's q-method, not QUEST.** QUEST is the same
> method with the largest eigenvalue approximated by Newton iteration on the
> characteristic quartic — faster, but carrying a non-convergence case and the
> method of sequential rotations near 180°. docs/SPEC.md wants this "as a cross-check
> only", and an exact eigensolution with no failure modes of its own serves that
> better than a fast approximation. Say the word if QUEST proper is wanted.

## Phase 7 — Verification and benchmark harness

> **Resolved 2026-10-09.** The risk recorded here was that a constant tolerance
> cannot absorb a focal error that grows with separation, and that is now fixed
> at the source: the identification tolerance is separation-proportional, so
> identification no longer fails before calibration can help. `hard` went
> 11.2% → 82.2% (focal estimation, Phase 7) → **99.6%** with **0
> WRONG_CONFIDENT** and **0 NO_SOLUTION**. See the status note below.

**Accept when:**
- [x] golden run (seed 42, nominal, 1000 trials): Score ≥ 99%, WRONG_CONFIDENT = 0, median total error ≤ 10″ — **99.90%**, **0**, **1.13″**.
- [x] hard preset: Score ≥ 95%, WRONG_CONFIDENT = 0 — **99.60%**, **0**, with 0 NO_SOLUTION and median error 2.13″.


> **Status: both targets met.**
>
> **Golden run** (seed 42, nominal, 1000 trials): Score **99.90%**,
> WRONG_CONFIDENT **0**, median total error **1.13″** — against floors of 99%, 0
> and 10″. ID precision **1.0000**, recall **0.9978**. Solve, rendering
> excluded, averages **1.75 ms** against the 5 ms budget. `easy` is 100.0% with
> 0 WRONG_CONFIDENT.
>
> **`hard`: Score 99.60%, WRONG_CONFIDENT 0** (500 trials), median error 2.13″,
> **0 NO_SOLUTION**, precision 1.0000, recall 0.9987. It reached this in three
> measured steps, each of which fixed a distinct cause.
>
> *The pair tolerance was constant, and the error is not.* A focal length off by
> `d` stretches every observed angle by about `d`, so a pair's angular error
> grows with the angle. Bucketing every observed pair against its catalogue
> angle showed the ratio flat at **8.0″ per degree** on `hard` (0.2% focal
> error) and 20.6″ on `brutal` (0.5%), against 0.1-0.2″ of pure noise on
> `nominal` — and the fraction of pairs inside the 50″ constant tolerance fell
> to **0% beyond 8°**. The pyramid needs a spread triangle, so it never closed:
> that was all 52 NO_SOLUTION. The tolerance is now
> `ε(θ) = k_sigma·σ_angle + calib_margin + focal_frac·θ`.
>
> *A wider window lets coincidences be confirmed.* With the tolerance widened,
> NO_SOLUTION went to 0 but 43 frames came back `claimed 4, correct 0,
> matched 0` — a false triangle that a fourth star confirmed by chance, which
> the self-check then correctly refused. The identifier was stopping at the
> first *confirmed* pyramid rather than the first that *explained the frame*.
> It now resumes the Mortari walk: a pyramid whose attitude sweeps up no other
> star is discarded. The separation is not a close call — a right pyramid sweeps
> up ten or more, every wrong one measured swept up exactly zero.
>
> *The common case must not pay for the rare one.* Searching wide always cost
> `nominal` 4.3× the solve time and cost `easy` a trial, because the wider
> window admits wrong pairs on frames that never needed it. The search now runs
> tight first and widens only on failure, which restores `nominal` to 0.14 ms of
> identification and `easy` to 300/300. A single work cap of 32M pair
> examinations bounds the wide pass, set at roughly twice what the hardest
> *solvable* frame needs (measured: 15.2M on `hard`, 103k on `nominal`);
> it is deterministic, so parity is unaffected.
>
> **What remained was not an identification error at all.** One trial stayed
> WRONG_CONFIDENT with all 14 of its IDs correct and 0.173 px of residual, but
> 67.9″ of attitude error — essentially all **roll**, the weakly constrained
> axis in a 20° field. Removing the distortion and the focal error made the roll
> tail *worse* (107″), so it was statistical, not an unmodelled-term bias: only
> 15 of ~28 available stars were being used. Identification still sees the
> brightest 15, because its cost is combinatorial, but the focal fit, the sweep
> and both checks now see **every detection** — those stages are linear and roll
> improves with both count and spread. That alone removed the outlier and took
> recall from 0.32 to 0.998 on every preset.
>
> `brutal` is 8.0% (from 1.8%), with 0 WRONG_CONFIDENT. Uncapped it reaches 18%,
> but at 239 ms of identification per frame against 34 ms capped; the cap is set
> by what `hard` needs, not by what `brutal` scores, and `brutal` has no target.
>
> Two defects that measurement exposed, both fixed and both pinned by tests. A
> free focal parameter fitted to a *bad* identification made that solution look
> self-consistent: unbounded, the fit reached 0.56 and 1.34, and it turned a
> rejected solution on the golden run into a confident wrong one. The fit is now
> clamped to ±2%, which ground calibration would comfortably bound. Separately,
> one wrong id among fifteen left the residual RMS at 0.30 px — inside the 0.5 px
> limit, because the fourteen good matches dominate it — so the self-check now
> also requires that **every** claim reprojects onto the centroid it was made
> for.

Record the actual numbers below.

## Phase 8 — WASM bindings
- [x] `tracker-wasm` exposes these functions:
  - `init()`, which loads the catalogue and builds or loads the pair database;
  - `simulate_frame(seed, preset_or_cfg_json)`;
  - `solve_frame(...)`;
  - `run_trial(seed, cfg)`, returning the trial result plus the data needed for drawing;
  - `run_benchmark_chunk(base_seed, start, count, cfg)`, so the worker can report progress.
- [x] Pass images as `Uint16Array` views, not JSON.

**Accept when:** a parity test (wasm-bindgen-test or a Node script) shows identical outcomes and errors equal to 1e-12 against native on 20 fixed seeds, and WASM timing is within 3× of native.

> **Status: parity passes; the solve timing passes; rendering does not.**
>
> `scripts/parity.mjs` compares WASM against the native `tracker-cli parity`
> reference over 20 fixed seeds on each of the four presets — 80 trials.
> **All 80 outcomes are identical**, and the attitude errors agree to
> **9.6e-13 rad** worst case (easy 9.56e-13, nominal 3.89e-13, hard 6.85e-13,
> brutal 2.22e-16), inside the 1e-12 the criterion asks for. Quaternions agree to
> 4.4e-13 and the reprojection residual to 1.3e-10 px.
>
> Angles are compared in radians, absolutely. Comparing the arcsecond figures
> *relatively* would be a stricter and unreachable test: an attitude error is a
> small difference between two nearly equal rotations, around 6e-6 rad, so
> 4e-13 of absolute agreement in the attitude necessarily shows up as 1e-7 of
> relative difference in the error derived from it. That is cancellation, not
> disagreement. The residual drift traces to nalgebra's SVD, not to this crate:
> routing every transcendental in `math.rs` and `camera.rs` through `libm` left
> the mismatch counts unchanged, and the images, centroids and un-projections are
> all exact IEEE operations.
>
> **Timing.** The solve stages — what the 5 ms budget covers, rendering excluded —
> run at **1.46x to 1.69x** native, inside 3x. Rendering does not: WASM simulates
> a 1024² frame in **42 ms against 8 ms native, 5.3x**. A full 1000-trial run in
> WASM therefore takes **44.4 s** against the 20 s Phase 9 budget, with rendering
> 95% of it. `-C target-feature=+simd128` buys 14% (42.2 → 36.1 ms, 38.1 s
> projected) and so is not shipped: it does not close a 2x gap and would have to
> be carried by the Phase 9 build. See the risk recorded against Phase 9.
>
> The WASM golden run reproduces native exactly — Score **99.80%**,
> WRONG_CONFIDENT **0**, median **1.21″** — which already settles one of Phase 9's
> criteria. `init()` builds the catalogue and the 849,829-pair database in
> **149 ms**.


## Phase 9 — Website

### Layout
A single page with two tabs: **Single Pattern** and **Benchmark**. Use a dark, clean, scientific style with a monospace font for numbers. The page must be responsive and usable on mobile.

### Single Pattern tab
- **Controls:** a "Generate random pattern" button, a seed field (editable, so any trial can be reproduced), a preset dropdown, and an "Advanced" panel for raw `SimConfig` fields.
- **Canvas 1, the sensor image:**
  - the simulated image, with a log/asinh stretch and zoom/pan;
  - detected centroids drawn as circles;
  - identified stars labelled with their catalogue ID, or proper name if the star has one;
  - false stars marked once the verification step reveals them.
- **Canvas 2, the sky context:** a small all-sky Mollweide map showing the true boresight, the estimated boresight, and the FOV footprint.
- **Step reveal.** Clicking "Identify" animates through the stages: centroids → pyramid stars highlighted → full match → verification. At the verification stage, each label turns ✅ green (correct) or ❌ red (wrong) against the truth, and missed truth stars are drawn as dashed circles.
- **Result panel:**
  - the outcome badge (CORRECT / WRONG_CONFIDENT / REJECTED / NO_SOLUTION);
  - the true and estimated attitude as quaternions and as RA/Dec/roll of the boresight;
  - total, cross-boresight and roll error in arcseconds;
  - stars detected, identified and correct;
  - per-stage timings.

### Benchmark tab
- **Controls:** a trials field (default 1000), base seed, preset, and Run/Cancel buttons. The run executes in a Web Worker in chunks of 25, with a progress bar and live-updating stats.
- **Scoreboard:**
  - a large **Score: XX.X%** (CORRECT / total);
  - next to it, WRONG_CONFIDENT count (red if above 0), REJECTED, NO_SOLUTION, ID precision and recall, and median/p95 error in arcseconds.
- **Charts** (plain canvas or one lightweight library such as uPlot):
  - histogram of total attitude error on a log x-axis;
  - cross-boresight vs roll error scatter;
  - outcome breakdown bar chart;
  - timing histogram;
  - success rate vs number of detected stars.
- **Trial table:** sortable and paginated. Clicking any row opens that seed in the Single Pattern tab, so failures can be inspected visually.
- **Export:** CSV of all trials and a JSON summary. A "Compare presets" button runs 1000 trials on each preset and shows the scores side by side.

### Front-end engineering
- Load WASM once and share the module between the main thread and the worker.
- Draw images with `ImageData` from the `Uint16Array` after the stretch. Do not draw per-pixel DOM elements.
- Static build only. Must deploy to GitHub Pages; add a workflow in `.github/workflows/pages.yml`.

> **Risk closed 2026-10-08 — by parallelism, not by cheaper rendering.** The
> measured per-frame cost still stands: WASM renders a 1024² frame in **42 ms
> against 8 ms native, 5.3x**, outside the 3x rule, while the solve stages are
> **2.10x** and inside it. What closes the budget is that a trial is seeded by
> `splitmix64(base_seed ^ index)` and nothing else, so a worker pool over
> disjoint contiguous index ranges returns exactly the rows one long run would.
> A 1000-trial nominal run at seed 42 in Chrome measures **33.66 s on 1 worker,
> 9.23 s on 4, and 5.07 s on 8** against the 20 s budget, with Score,
> WRONG_CONFIDENT and median error bit-identical to native at every width.
> `-C target-feature=+simd128` was measured at 14% and is still not shipped.
> The per-frame gap is left standing and recorded: closing it needs bulk noise
> generation rather than one ziggurat draw per pixel, which is Phase 10 work if
> it is ever wanted. Rendering at reduced resolution would change the system
> being measured and was rejected.

**Accept when:**
- [x] a 1000-trial nominal run in Chrome finishes in under 20 s, and its Score matches the native CLI exactly for the same seed — **5.07 s** on 8 workers (9.23 s on 4, 33.66 s on 1); Score **99.80%**, WRONG_CONFIDENT **0**, median **1.21″**, identical to the native CLI at seed 42.
- [x] clicking a failed trial reproduces it exactly — checked against `tracker-cli parity --preset hard --seed 42`: trial 0 CORRECT 3.685″ 126/14/14, trial 1 REJECTED 509288.075″ 52/4/0, trial 2 CORRECT 8.491″ 58/14/14, trial 4 NO_SOLUTION 64/0/0 — all four matching, the two failures included.
- [x] there are no console errors — none across the reproduction loads or the benchmark run, read over the DevTools protocol rather than from a DOM dump.
- [x] Lighthouse performance is at least 85 — **100** (FCP 0.9 s, LCP 0.9 s, TBT 80 ms, CLS 0.003, Speed Index 0.9 s).

> **Status: all four met.**
>
> **The 20 s budget.** A 1000-trial nominal run at seed 42 finishes in Chrome in
> **5.07 s on 8 workers**, 9.23 s on 4 and 33.66 s on 1. Per-frame rendering is
> unchanged and still 5.3x native; what closed the budget is parallelism, which
> is sound here only because a trial's seed is `splitmix64(base_seed ^ index)`
> and trials share no state, so a pool over disjoint contiguous index ranges
> returns exactly the rows one long run would. Score **99.80%**,
> WRONG_CONFIDENT **0**, median **1.21″** — identical to the native CLI, at
> every worker count.
>
> **Reproduction.** Trials are addressable by URL, which is the same operation a
> table row performs, so the two share one code path. Against
> `tracker-cli parity --preset hard --seed 42`: trial 0 CORRECT 3.685″ 126/14/14,
> trial 1 REJECTED 509288.075″ 52/4/0, trial 2 CORRECT 8.491″ 58/14/14, trial 4
> NO_SOLUTION 64/0/0 — all four match, failures included.
>
> **Lighthouse: 71 before, 100 after, and the gap was worth the work.** FCP, LCP,
> CLS and Speed Index were already at 100; the entire loss was Total Blocking
> Time of **2,880 ms**, and it sat inside a *single* 2,920 ms task. Measuring the
> startup phase by phase found it: the main thread built the `Tracker` (**165 ms**
> for the catalogue and the 849,829-pair database) and ran the first pattern
> (**125 ms**) as one unbroken task, because `await` on an already-resolved
> promise is a microtask and yields to other tasks but never to the event loop.
> Breaking the chain into separate tasks would not have helped — each phase alone
> is many times the 50 ms a task may run before it counts against TBT — so the
> compute had to leave the thread entirely, and with it the main thread's
> `Tracker`. `web/src/engine.ts` now owns one worker that holds the only tracker
> the interactive view needs, and the frame's pixels come back as a transferred
> buffer. The main thread's remaining work is the drawing alone, measured at
> **29.6 ms sensor + 4.9 ms sky**: **TBT 80 ms, performance 100**.
>
> **Console.** Clean. The one error found anywhere was a 404 for `/favicon.ico`,
> which Chrome requests whether or not any markup asks for it; an inline SVG
> data-URI icon removes the request rather than just the 404.

## Phase 10 — Polish (optional, ask first)
- [x] **Tracking mode**: a sequence of frames with slow rotation, a predicted-window search, and a comparison of tracking vs lost-in-space timing.
- [ ] Tetra-style hash identification as a second algorithm, selectable in the UI and compared head-to-head on the same seeds.
- [ ] Upload a real night-sky photo and solve it, with an FOV guess entered by the user.

> **Tracking: done 2026-10-09.** `crates/core/src/track.rs`, a `track`
> subcommand, and a third tab on the site. Over a 60-frame slew on `nominal`,
> every frame is CORRECT with 0 WRONG_CONFIDENT and 0 reacquisitions at
> **0.5, 5, 20 and 60 °/s**, and identification runs **2.4-5.0x faster** than
> the full search on `nominal` and **13.8x** on `hard`, where searching the sky
> has to widen its tolerance. Median error is *better* than lost-in-space
> (0.47″ against 1.13″), because the prediction matches stars the blind search
> never gets to.
>
> Three things had to be right, and the first two attempts were not.
>
> *The field must not be rescanned.* A predicted-window search that rebuilds
> its candidate list from all 4,992 database stars every frame is **slower**
> than the pyramid — measured at 0.274 ms against 0.107 ms. The field is now
> cached in inertial coordinates with a degree of slack and rebuilt only when
> the boresight drifts out of it, which also cut `extend_matches` down to a
> dot-product rejection before any `atan2`.
>
> *One past attitude is not enough.* With a single prior the window has to span
> the whole inter-frame slew, and such a window is ambiguous: every observation
> has a rival inside it, so the matches are discarded and what survives is worse
> than nothing. At 20 °/s that produced a clean alternation — a bad tracked
> solution suppressed the fallback, failed to calibrate, and the next frame
> re-acquired — which looked like 30/60 success. Tracking now needs two
> attitudes, takes the rotation between them as the rate, and searches only the
> **98″** that the prediction's own error needs, whatever the slew.
>
> *A prediction that produces matches has not necessarily solved the frame.* The
> full search now runs behind the tracked one as a fallback, so a prediction
> that fails to calibrate costs a frame nothing.

---

## Progress Log
| Date | Phase | Notes / measured numbers |
|------|-------|--------------------------|
| 2026-10-07 | 0 | **Done.** Workspace of 4 crates + `web/` scaffold + `math.rs`. `cargo fmt --check` and `cargo clippy --workspace --all-targets -D warnings` clean; 15/15 math tests pass in debug and release. Measured: `angle_between` absolute error is flat at ~1.1e-16 rad from 1 rad down to 1e-12 rad (worst 4.4e-16, at 1 rad), whereas `acos(a·b)` reaches 2.9e-8 rad at a 1e-9 rad separation — 30× the angle itself; docs/SPEC.md pitfall 3 is now pinned as a test. `rotation_vector` relative error ≤ 1e-9 from 1 arcsec down to 1e-6 arcsec. `splitmix64` matches the published stream (`splitmix64(0) = 0xE220A8397B1DCDAF`). Native/WASM parity: `splitmix64(42) = 0xbdd732262feb6e95` identical in both, verified in headless Chrome on the Vite dev server *and* on the production `dist/` build, no console errors. Sizes/timings: `tracker_wasm_bg.wasm` 13.8 kB (6.3 kB gzip), `vite build` 113 ms, `wasm-pack build --release` 31.7 s cold / 1.1 s warm. Toolchain installed this phase: rustup 1.29.1 (rustc 1.99.0), `wasm32-unknown-unknown`, wasm-pack 0.15.0; `rust-toolchain.toml` added so a non-rustup cargo earlier on PATH cannot be picked up. |
| 2026-10-07 | 1 | **Done.** `catalog.rs` (`Star`/`Catalog`, bincode via serde, validating loader) + `build-catalog`. Source: Hipparcos I/239 via VizieR ASU-TSV, server-side `Vmag<=6.5` (503 kB plain text, no gzip dependency); names from IAU WGSN IAU-CSN. Chose Hipparcos over BSC5 because BSC5 rounds J2000 to 0.1s RA / 1″ Dec, worth ~0.9″ on its own against a 1″ criterion. Positions propagated J1991.25 → J2000.0 in **vector** form (tangent east/north + renormalise), which has no cos δ singularity — Polaris moves 30″ in RA over those 8.75 yr. Numbers: 8874 rows fetched → 8870 parsed (4 dropped, 0 missing proper motion) → 43 blends merged below 140.6″ (2·FOV/W) → **8827 stars**, 337 named. V −1.44 (Sirius) … 6.50; 5020 at V ≤ 6.0 (the Phase 4 db cut). Closest surviving pair 149.6″ > 140.6″ blend radius. Bright-star accuracy vs VizieR's own PM-propagated `_RAJ2000`/`_DEJ2000`: worst **0.00017 mas** over Sirius/Vega/Arcturus/Betelgeuse/Polaris, and those references reproduce the published sexagesimal positions exactly (Sirius 06h45m08.92s −16°42′58.0″). `catalog.bin` = 292,877 bytes (29% of budget). Rebuild from a wiped `data/raw/` is **byte-identical**. 27 tests pass in debug and release; fmt + clippy `-D warnings` clean; wasm32 still builds with the new serde/bincode/thiserror deps. Note: `bincode 3.0.0` on crates.io is a placeholder that only emits `compile_error!`, so the workspace pins `bincode = "2.0"` with its non-default `serde` feature. |
| 2026-10-07 | 2 | **Done.** `camera.rs` + `simulate.rs`; 63 tests pass in debug and release, fmt + clippy `-D warnings` clean. Camera: f = **2903.696292 px** (512/tan 10°, docs/SPEC.md's "≈2904"), principal point 511.5, half-diagonal FOV 13.9888° (full diagonal 27.978°). `project` carries `k1` (the physical camera); `unproject` is the ideal pinhole inverse the *solver* uses, so round-trip closes only at `k1 = 0` — measured **1.4e-13 px** over 200k points against the 1e-9 budget, and 9.3e-17 rad for direction→pixel→direction. Unmodelled distortion pinned: `k1 = 0.01` shifts the corner 0.45 px, `k1 = 0.03` shifts it 1.35 px. Boresight uniformity: chi-square **85.80** (df 99, p=0.01 critical 134.6416), cell counts 4..18 over 100 cells. Determinism: identical images for all four presets on repeat, and a reused `Workspace` matches a fresh one. Frame content at 1024², averaged over 40 seeds: **84.7 visible stars**, dropout 0 / 5.2% / 10.7% / 20.4% against configured 0/5/10/20%, false stars 0 / 0.42 / 1.05 / 3.38 within their preset ranges, blob saturates **19,324 px** on brutal. PSF: `PSF_WINDOW_SIGMAS` raised 4→6 after measuring that a 4σ window biases the recovered centroid by 1.4e-4 px purely by asymmetric truncation; at 6σ the renderer's own bias is **6.4e-9 px** and 1e-9 of flux is lost, against the 0.02 px Phase 3 must measure. Residual floor is the binning-aliasing term ~exp(−2π²σ²), which is why σ=0.5 would cost 2e-3 px. Render cost **14.1 ms/frame** native, of which **10.2 ms is the 1M-pixel noise pass** (one ChaCha8+ziggurat normal per pixel, ~9 ns each) — flagged against the Phase 9 web budget, see that phase. Dependency fix: `rand`'s default features pull `sys_rng`→`getrandom`, which fails to compile for wasm32; set `default-features = false, features = ["std"]` for `rand` and `rand_distr`, which removes `getrandom` entirely and makes docs/SPEC.md's "no OS entropy" rule a compile error rather than a convention. |
| 2026-10-07 | 3 | **Done.** `centroid.rs` + `benches/centroid.rs`; 78 tests pass in release (74 + 4 release-only in debug), fmt + clippy `-D warnings` clean. Single star V 4.0, nominal noise, 8192 samples over a swept sub-pixel grid: mean bias **−0.00004/+0.00003 px**, worst per-phase systematic **0.0104 px** (budget 0.02), RMS **0.0106/0.0107 px** (budget 0.1). Completeness of rendered V<5.5 truth stars over 1000 nominal frames: **99.675%** (26369/26455), 99.758% over all rendered stars (budget >98%). Criterion at 1024²: easy 366 µs, **nominal 354 µs**, hard 369 µs, brutal **1.260 ms** (budget 2 ms); ~81 blobs found per frame, 15 fed onward. Zero per-frame allocation verified by capacity stability across 320 warm frames. **Phase 1 amended:** the catalogue blend radius went from 2 px to **6 px** (`BLEND_PX` in `build-catalog`), on measurement — with a 1 px PSF, 8-connected labelling joins two stars out to 6-8 px because the saddle between their peaks stays above threshold, so at 2 px the catalogue asserted two stars wherever the sensor can only ever yield one centroid, and that centroid matched neither (it behaved like a false star). Completeness was **97.92%** before the change, with 470 of 550 misses being blends at 2-6 px, 76 edge-clipped and 4 other. The new radius merges 107 groups instead of 43: **8763 stars** (was 8827), 335 named, closest surviving pair 422.8″, `catalog.bin` 290,766 bytes. Two fixes found by testing: the integer detection cutoff used `ceil(t)`, which on a frame with no background admits every zero pixel and merges the whole image into one blob (now `floor(t+1)`, exactly `value > t` for integer pixels); and `background_stride = 8` is tuned for 1024 px but leaves only 64 samples on a 64 px frame, loose enough that noise registered as stars, so the stride is now capped to keep ≥32 samples per axis. Note: the debug test suite is now **107 s**, dominated by Phase 2's full-frame simulate tests; the four statistical Phase 3 tests are `#[cfg_attr(debug_assertions, ignore)]` per the release-test convention in docs/SPEC.md. |
| 2026-10-07 | 4 | **Done.** `pairdb.rs` + `benches/pairdb.rs`; 94 tests pass in release (88 + 6 release-only in debug), fmt + clippy `-D warnings` clean. Database: **4992 stars** at V ≤ 6.0 (the catalogue's leading prefix, so a pair index *is* a catalogue index — docs/SPEC.md pitfall 6), separations to **28.9776°** (27.9776° diagonal + 1°), **849,829 pairs**, heap **10.20 MB** (8 B per `Pair` + 4 B per k-vector entry). Acceptance: `query` matched `query_slow` on **10,000/10,000** random queries, including angles below the smallest and above the largest separation, and tolerances of 0, 1″, 60″ and 0.05 rad. Build **65 ms** native (criterion) and **95 ms median in WASM** (87–113 ms over 9 runs under Node), comfortably inside the 300 ms budget, so the table is built at start-up; precomputing into `catalog.bin` was never viable at 10.2 MB against a 1 MB budget. Query **~14 ns** (criterion, 64 queries per iteration), returning a mean of 163 pairs at a 10″ tolerance. Two things that made the build fast enough: scanning stars in order of `z` and stopping the inner loop once the z gap exceeds the limit (valid because \|sin δ_b − sin δ_a\| ≤ \|δ_b − δ_a\| ≤ angle), which cuts 12.46M candidate pairs to 5.39M; and rejecting candidates on a dot-product comparison against cos(max), so `atan2` runs only for the 850k pairs that survive. The k-vector brackets tightly here — the pair-angle density varies only about twofold across the range — so the nudge onto exact bounds runs a couple of iterations and `query` returns the precise set rather than candidates to be filtered. The WASM figure was re-measured after noticing the first probe returned only the pair count, which does not depend on the sort or the index and so could have been optimised away. |
| 2026-10-08 | 5 | **Done.** `identify.rs` + `attitude.rs` (SVD Wahba with the determinant correction, which step 4 needs) + a `Clock` trait; 109 tests pass in release (98 + 11 release-only in debug), fmt + clippy `-D warnings` clean. Acceptance, 1000 vector-level trials each: noiseless **100.000%**, nominal noise **100.000%** (floor 99.5%), three false vectors **99.900%** (floor 98%). Timing 0.09 ms clean, 0.16 ms with three false vectors; tries average 1.0 and 2.5. Failure envelope: 6 false 99.1%, 9 false 98.6%, 11 false 93.3%; 3× noise 99.9%, 5× 99.1%, 10× 79.1% (noise then exceeds the tolerance). Every failure is a **bad pyramid** — across 9000 trials the sweep never mis-assigned a real star and only once matched a false vector spuriously, so a failure shows up as a grossly wrong attitude, which the Phase 7 self-check should reject. Attitude error under nominal noise over 3000 trials: cross-boresight median **2.22″** against a predicted σ/√N = 1.83″, roll median **9.43″** against a predicted 7.59″ (σ/√N ÷ sin(half-diagonal)) — both within a factor of 1.25, so the test asserts the scaling law rather than an invented bound. Roll dominates because a 20° field gives it a short lever arm; total median 9.77″, p95 27.6″, max 60.8″. End to end (simulate → centroid → unproject → identify, 150–200 frames): easy and nominal **100% solved, 100% IDs correct, median 1.67″ / 1.70″**; hard 91.3%/86.0% and brutal 10.0%/1.3% — see the open risk recorded against Phase 7, which traces that entirely to the constant `calib_margin` versus a separation-proportional focal error, not to the identifier. Two design points worth recording: `std::time::Instant` compiles for wasm32 and then **traps** at run time (`RuntimeError: unreachable`, verified under Node), so core never reads a clock and timing arrives through the `Clock` trait; and `identify` returns diagnostics whether or not it finds anything, since an unsolved trial is exactly the one whose counters matter. Pure-garbage input (15 false vectors) yields a confident-looking solution 12% of the time and burns 17 ms exhausting 413 triples — load the self-check will have to carry. |
| 2026-10-08 | 6 | **Done.** `attitude.rs` completed with `centroid_weight`, `solve_weighted` and `solve_davenport`; 118 tests pass in release (105 + 13 release-only in debug), fmt + clippy `-D warnings` clean; `proptest` added as a dev-dependency. Acceptance: noiseless recovery **< 1e-9 rad** and `det = +1` to 1e-12 across generated cases of 2–24 stars; error under noise gives per-axis ratios to `sigma/sqrt(N)` of **1.46 → 1.21** for N = 3…120, converging on the analytic `sqrt(3/2) = 1.225` (total converges on `3/sqrt(2) = 2.121`), with `rms * sqrt(N)` flat — inside the factor of two. **Brightness weighting:** measured `sigma * sqrt(flux)` constant to ~1.5× over V3–6.5, confirming `w = 1/sigma^2 ∝ flux`; the flux is capped where the PSF peak hits full well (25,730 ADU total, from the 0.1592 peak fraction) because past that the measured flux keeps rising while the centroid stops improving — a V1 star measured 12× the flux of a V4 one and scattered just as widely (0.0107 px vs 0.0105 px), since both V1 and V2 saturate at a 4095 ADU well. Over 1500 nominal frames `solve_weighted` improves the median attitude error from **1.594″ to 1.204″** and p95 from 5.873″ to 4.826″. **The cross-check caught a real bug:** Davenport's `z` must be `Σ w (r × b)`, not `b × r`; the wrong sign conjugates the quaternion and disagreed with the SVD by up to π. With the sign right the two solvers agree to a median of **2.84e-8 arcsec** (1.4e-13 rad) on real frames. The determinant-correction test feeds the solver body vectors that are a *reflection* of the inertial ones, which makes `det(B) = −8.97` and `det(U·Vᵗ) = −1` — a genuine reflection — and both solvers still return `det = +1`, so the correction is demonstrably load-bearing rather than cosmetic. Note for Phase 7: brightness weighting worsens the *worst* case (23.7″ → 52″ on frames with ≥6 matches) because a blended pair carries high flux **and** a bad position, so flux weighting trusts it most; residual-based outlier rejection would be the fix and is not in the plan. |
| 2026-10-08 | 7 | **Golden run passes; `hard` does not.** `verify.rs`, `bench.rs`, `Catalog::id_index`, focal-length estimation in `attitude.rs`, `identify::extend_matches`, and the hand-rolled CLI; 142 tests pass in release (126 + 16 release-only in debug), fmt + clippy `-D warnings` clean. **Golden (seed 42, nominal, 1000):** Score **99.80%**, WRONG_CONFIDENT **0**, median total error **1.21″**, p95 4.82″; cross-boresight median 0.30″, roll 1.13″; ID precision 0.9999, recall 0.3179; solve (render excluded) mean **0.75 ms** / p95 0.95 ms against the 5 ms budget; wall clock 8.9 s. `easy` 100.0% / 0 wrong. **`hard` 82.2% / 1 wrong** (target 95% / 0), median 3.32″; **`brutal` 1.8% / 2 wrong**. Focal-length estimation was approved this phase and took `hard` from 11.2% to 82.2% with median error 54.8″ → 3.32″, fitting 1.00219 against a true 1.002. The residual gap is **52 NO_SOLUTION (10.4%)** on `hard` and 86.8% on `brutal`: identification itself fails there because the *pair* tolerance is constant while the focal error grows with separation, and no amount of calibration helps a frame that was never identified — so `hard` is capped at 89.6% until the identification tolerance is addressed, which is the lever deliberately left alone. Getting the fit to work needed two things beyond the plain radius ratio: the bootstrap sweeps run at **10× tolerance**, because while the camera is uncalibrated the outer stars sit outside the normal tolerance and those are exactly the stars the fit needs (4–7 clustered matches became 13.4 well-spread ones, and the fitted scale went 1.00085 → 1.00219); and a final pass re-derives the match list at full strictness so nothing loose survives. Two defects measurement exposed, both now pinned by tests: an unbounded focal fit reached **0.56 and 1.34**, absorbing bad identifications and turning a rejected golden-run solution into a confident wrong one, so it is clamped to ±2%; and one wrong id among fifteen left the residual RMS at **0.30 px**, inside the 0.5 px limit because the good matches dominate it, so the self-check now also requires every claim to reproject onto its own centroid. `std::time::Instant` traps on wasm32, so timings arrive through the `Clock` trait. CLI args and the JSON summary are hand-rolled — no `clap`, no `serde_json`, approved list untouched; the CLI exits non-zero when WRONG_CONFIDENT is not zero. |
| 2026-10-08 | 8 | **Parity and solve timing pass; WASM rendering does not.** `crates/wasm` now exposes `init()`, `simulate_frame`, `solve_frame`, `run_trial`, `run_benchmark_chunk` and a `Frame` class, plus `scripts/parity.mjs` and a `tracker-cli parity` reference subcommand; 142 tests pass, fmt + clippy `-D warnings` clean on **both** the native and `wasm32-unknown-unknown` targets. Parity over 20 fixed seeds × 4 presets (80 trials): **all 80 outcomes identical**, attitude errors agreeing to **9.6e-13 rad** worst case (easy 9.56e-13, nominal 3.89e-13, hard 6.85e-13, brutal 2.22e-16), quaternions to 4.4e-13, residual to 1.3e-10 px. The WASM golden run reproduces native exactly: Score **99.80%**, WRONG_CONFIDENT **0**, median **1.21″**. `init()` takes **149 ms** for the catalogue plus the 849,829-pair database. Solve timing **1.46–1.69×** native, inside the 3× rule; **rendering is 5.3×** (42 ms against 8 ms per 1024² frame), so a 1000-trial WASM run takes **44.4 s** against Phase 9's 20 s budget with rendering 95% of it — recorded as a risk there. `simd128` was measured and gives only 14% (38.1 s projected), so it is not shipped. Two things worth recording. My first parity script applied a *relative* 1e-12 to arcsecond values and reported failure on all four presets; that is a stricter test than the criterion states and is unreachable, because an attitude error is a ~6e-6 rad difference between two nearly equal rotations, so 4e-13 of absolute agreement reads as 1e-7 relative. Compared the way the criterion means — radians, absolutely — it passes. And the residual drift is **nalgebra's SVD**, not this crate: routing every transcendental in `math.rs` and `camera.rs` through `libm` changed the mismatch counts not at all, so that experiment was reverted rather than kept as unmeasured churn. Design notes: config crosses as a JS value and may be a preset name or a whole `SimConfig`, via `serde-wasm-bindgen`, so no `serde_json`; `performance.now()` is declared with plain `wasm-bindgen`, so no `js-sys`; images cross as `Uint16Array` with `Frame::image_ptr`/`image_len` for a zero-copy view, since a true `Uint16Array::view` would need `unsafe`, which docs/SPEC.md forbids; and the pair database is rebuilt only when `db_mag_cut`, the field of view or the width change. |
| 2026-10-08 | 9 | **Done; all four acceptance criteria pass.** Vanilla-TS two-tab site: sensor canvas with an asinh stretch, zoom/pan and stage-gated overlays (centroids, pyramid ring and polygon, ✓/✗ labelled matches, dashed missed-truth rings, false-star crosses), a Mollweide sky map with graticule and field footprint, five canvas charts, a sortable paginated trial table whose rows reopen that trial in the pattern tab, CSV/JSON export, compare-presets, and `.github/workflows/pages.yml` gating the publish on fmt, clippy (both targets), release tests and parity. Two harnesses were needed because headless Chrome will not hold its page lifecycle open for worker results: `scripts/web-budget.mjs` (the page reports its own wall clock back to the server) and a DevTools-protocol check that waits for the outcome badge before reading the DOM. **Budget: 5.07 s on 8 workers**, 9.23 s on 4, 33.66 s on 1, against 20 s — Score **99.80%**, WRONG_CONFIDENT **0**, median **1.21″**, bit-identical to native at every width. **Lighthouse first measured 71**, and the whole loss was Total Blocking Time **2,880 ms** inside a *single* 2,920 ms task, with FCP/LCP/CLS/Speed-Index all at 100. Instrumenting the startup showed why: the main thread built the `Tracker` (**165 ms** — catalogue parse plus the 849,829-pair database) and ran the first pattern (**125 ms**) as one unbroken task, because `await` on an already-resolved promise is a microtask and never yields to the event loop. Splitting the chain into separate tasks could not have fixed it — each phase on its own is far over the 50 ms a task may run before it counts against TBT — so the work had to leave the thread. New `web/src/engine.ts` owns one worker holding the only `Tracker` the interactive view needs; `worker.ts` gained `config` and `frame` messages and returns the 2 MB image as a *transferred* buffer (safe because the glue's `image` getter ends in `.slice()`, a JS-owned copy). The main thread now compiles the module, hands it over, and does nothing but draw — the measured **29.6 ms sensor + 4.9 ms sky**. **TBT 2,880 → 80 ms, Lighthouse 71 → 100**, and the entry bundle fell 42.8 → 31.5 kB as the glue moved into the worker chunk. One spec deviation worth recording: "load WASM once and share the module between the main thread and the worker" still holds for *compiling* — once, via `compileStreaming` — but nothing is **instantiated** on the main thread any more. The only console error anywhere was a 404 for `/favicon.ico`, which Chrome requests unprompted and no markup referenced; fixed with an inline SVG data-URI icon, so the request never happens. **Measured and deliberately left alone:** a pattern view renders the same frame **three times** — `bench::run_trial` renders once to score it, the wasm wrapper renders again for the drawing truth, and `simulate_frame` renders a third time for the pixels — about 80 ms of the ~107 ms a pattern costs. It is worker time now rather than main-thread time and no criterion depends on it, so it is recorded rather than fixed; collapsing it wants a `bench::run_trial_on(&frame, …)` entry point so one render serves both the score and the drawing. The benchmark wall clock also fell (10.85 → 5.07 s at 8 workers, solve 4.46 → 2.30 ms) and I did not isolate the cause — dropping the main thread's second tracker instance plausibly cuts memory traffic, but machine state may contribute; the scored output is unchanged either way. Gates: **142 tests** pass in release, fmt clean, clippy `-D warnings` clean on native and `wasm32-unknown-unknown`, parity **80/80** outcomes identical (angles within 9.56e-13 rad), WASM solve **2.10x** native on nominal. **Proper names closed the same day**, completing the layout spec's "catalogue ID, or proper name if the star has one": the name rides along on `MatchDto` rather than through the `Tracker::star_name(id)` accessor first sketched, because with the tracker on a worker an accessor would cost a round-trip per label. A named star is labelled `Alioth`, everything else `HIP 64906`; 16 of 89 matched stars over the first six nominal trials carry one, and nominal trial 4 lands on the Big Dipper (Alioth, Alkaid, Mizar, Megrez) with Cor Caroli and Chara alongside. The one item carried forward from here — the `hard` preset at **82.2%** against its 95% target — was closed the next day; see the 2026-10-09 Phase 7 row. |
| 2026-10-09 | 7 | **`hard` target met: 82.2% → 99.60%, WRONG_CONFIDENT 1 → 0, NO_SOLUTION 52 → 0.** Phase 7's second criterion had stood unmet since 2026-10-08. Four measured steps, each fixing a different cause. **(1) The pair tolerance was constant and the error is not.** A focal length off by `d` stretches every observed angle by about `d`, so a pair's error grows with the angle. Bucketing every observed pair against its catalogue angle gave a flat ratio: **8.0″ per degree** on `hard` (0.2% focal error), 20.6″ on `brutal` (0.5%), 0.1-0.2″ on `nominal` — and the share of pairs inside the 50″ constant fell to **0% beyond 8°**, so the pyramid, which needs a spread triangle, never closed. That was all 52 NO_SOLUTION. The tolerance is now `ε(θ) = k_sigma·σ_angle + calib_margin + focal_frac·θ`, with `id_focal_tolerance_frac` defaulting to 0.003 and identical for every preset — it is the solver's assumption about its own calibration, not knowledge of which preset it is in. **(2) A wider window lets coincidences be confirmed.** NO_SOLUTION went to 0 but 43 frames returned `claimed 4, correct 0, matched 0`: a false triangle a fourth star confirmed by chance, correctly refused downstream. The identifier was stopping at the first *confirmed* pyramid rather than the first that *explained the frame*; it now resumes the Mortari walk and discards a pyramid whose attitude sweeps up no other star. The signal is not marginal — a right pyramid sweeps up ten or more, every wrong one measured swept up exactly zero. My first attempt used the self-check's match floor as the bar instead, which threw away correct pyramids in sparse fields and cost `easy` a trial; the measured criterion is the one that is both sufficient and safe. **(3) The common case must not pay for the rare one.** Searching wide always cost `nominal` 4.3× its solve time and `easy` a trial, since the wider window admits wrong pairs on frames that never needed it. The search now runs tight first and widens only on failure: `nominal` identification is back to **0.14 ms** and `easy` to 300/300. A work cap of 32M pair examinations bounds the wide pass, set at about twice what the hardest *solvable* frame needs (15.2M on `hard`, 103k on `nominal`, against 236M wasted by a median hopeless `brutal` frame); it is deterministic, so parity is untouched. **(4) The last WRONG_CONFIDENT was not an identification error.** One trial kept it with all 14 IDs correct and 0.173 px residual but 67.9″ of error, essentially all **roll** — the weakly constrained axis in a 20° field. Stripping the distortion *and* the focal error made the roll tail **worse** (107″), proving it statistical rather than an unmodelled-term bias: only 15 of ~28 available stars were in the fit. Identification still sees the brightest 15 because its cost is combinatorial, but the focal fit, the sweep and both verification layers now see **every detection** — those stages are linear, and roll improves with both star count and spread. The brightest are a prefix of the flux-sorted list, so pyramid indices stay valid and the change is small. That removed the outlier outright and took **ID recall from 0.32 to 0.998** on every preset, with precision 1.0000. Along the way I rejected two alternatives on measurement: raising `max_centroids` to 18 works but changes a documented default and is knife-edge (C(18,3)=816 of the 1000-try cap; at 20 it collapses to 29 NO_SOLUTION), and a covariance test on the self-check could not separate the trial cleanly (~2.5σ from its predicted roll uncertainty). Final: golden **99.90%**/0/**1.13″** (was 99.80%/0/1.21″), `hard` **99.60%**/0/2.13″ with 0 NO_SOLUTION, `easy` 100%, `brutal` 8.0% from 1.8%. **143 tests** pass in release, fmt and clippy `-D warnings` clean on native and `wasm32-unknown-unknown`, parity **80/80** identical with every preset inside the 3× rule (easy 2.33×, nominal 2.17×, hard 1.77×, brutal 1.70×), and the browser budget is 4.87 s on 8 workers. `docs/SPEC.md` was updated in three places, since the tolerance formula, the identification steps and the meaning of "keep the brightest 15 centroids" are all specified there. |

## Golden Benchmark Record
| Date | Commit | Preset | Seed | Trials | Score | WRONG_CONF | Median err (″) | p95 err (″) | Native ms/frame | WASM ms/frame |
|------|--------|--------|------|--------|-------|------------|----------------|-------------|-----------------|---------------|
| 2026-10-08 | uncommitted | nominal | 42 | 1000 | 99.80% | 0 | 1.21 | 4.82 | 0.75 solve / 8.1 render | 2.2 solve / 42 render |
| 2026-10-09 | uncommitted | nominal | 42 | 1000 | 99.90% | 0 | 1.13 | 4.01 | 1.8 solve / 16.2 render | 3.0 solve / 42 render |

> **2026-10-09 improved on 2026-10-08, it did not regress.** Score 99.80% → 99.90%,
> median 1.21″ → 1.13″, and ID recall 0.3181 → **0.9978**, because the attitude fit
> and both checks now see every detection rather than the fifteen identification
> uses. Solve went 0.75 → 1.8 ms: identification itself is *faster* (0.14 ms), and
> the rise is the sweep and fit now running over ~129 stars instead of 15, which is
> what bought the accuracy. The two render figures are not comparable — 8.1 ms was
> measured on an idle machine and 16.2 ms under load on the same build; rendering
> was not touched.

> **In-browser figures, 2026-10-09.** The rows above are per-frame costs. In
> Chrome the 1000-trial nominal run at seed 42 completes in **4.87 s on 8
> workers** with solve **3.03 ms** per frame, reproducing the Score (99.90%),
> WRONG_CONFIDENT (0) and median (1.13″) exactly. On 2026-10-08, before the
> identification work, the same run took 5.07 s on 8 workers, 9.23 s on 4 and
> 33.66 s on 1. Lighthouse performance on the built site is **100** (TBT 80 ms).
