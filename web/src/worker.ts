// The page's compute worker.
//
// Two jobs run here, and neither belongs on the main thread. A `Tracker` costs
// ~165 ms to build (catalogue parse plus the pair database) and a single
// pattern ~90 ms more, so on the main thread the page would spend a third of a
// second unable to answer a click, right when it first loads.
//
// *Benchmarking*: each worker instantiates the already-compiled module and
// takes a disjoint range of trial indices. That partitioning is safe because a
// trial's seed is `splitmix64(base_seed ^ index)` and nothing else: trials
// share no state, so splitting them across workers yields exactly the rows a
// single long run would, in whatever order they come back. The main thread
// sorts by index.
//
// Rendering is 95% of a trial's cost and the only way past it on this target is
// more cores, so this is where the 20 s budget is actually met.
//
// *One pattern*: the Single Pattern tab asks for one scored trial and the frame
// it drew. The pixels come back as a transferred buffer, so the main thread's
// only remaining work is the drawing itself.

import {
  trackerFromModule,
  version,
  type FromWorker,
  type SimConfig,
  type SolveReport,
  type ToWorker,
  type TrackRun,
  type TrialDetail,
  type TrialRow,
} from "./wasm.js";
import type { Tracker } from "./pkg/tracker_wasm.js";

/**
 * The worker's global scope, with only what this file uses.
 *
 * `DedicatedWorkerGlobalScope` lives in TypeScript's WebWorker lib, which
 * cannot be combined with DOM in one project without conflicts, so the members
 * in play are declared here instead.
 */
interface WorkerScope {
  postMessage(message: FromWorker, transfer?: Transferable[]): void;
  onmessage: ((event: MessageEvent<ToWorker>) => void) | null;
}

const scope = self as unknown as WorkerScope;

let tracker: Tracker | undefined;

scope.onmessage = (event: MessageEvent<ToWorker>) => {
  const message = event.data;
  try {
    if (message.kind === "init") {
      tracker = trackerFromModule(message.module);
      // The counts travel with the handshake because the main thread shows
      // them in the build banner and has no tracker of its own to ask.
      scope.postMessage({
        kind: "ready",
        version: version(),
        catalogStars: tracker.catalog_stars,
        databasePairs: tracker.database_pairs,
      });
      return;
    }

    if (!tracker) throw new Error("worker was asked to run before it was ready");

    if (message.kind === "config") {
      scope.postMessage({
        kind: "config",
        id: message.id,
        config: tracker.default_config(message.preset) as SimConfig,
      });
      return;
    }

    if (message.kind === "track") {
      scope.postMessage({
        kind: "track",
        id: message.id,
        run: tracker.run_tracking(
          BigInt(message.baseSeed),
          message.config,
          message.track,
        ) as TrackRun,
      });
      return;
    }

    if (message.kind === "frame") {
      // `run_trial` scores the trial exactly as the benchmark would, and
      // `simulate_frame` gives the image to draw. Both are driven from the same
      // derived seed, so what is drawn is what was scored.
      const detail = tracker.run_trial(
        BigInt(message.baseSeed),
        message.index,
        message.config,
      ) as TrialDetail;
      const frame = tracker.simulate_frame(
        BigInt(detail.trial.seed),
        message.config,
      );
      // Read once: the getter copies the image out of WASM memory each time.
      const image: Uint16Array = frame.image;
      const report = tracker.solve_frame(image, message.config) as SolveReport;
      // Transferred, not cloned: 2 MB of pixels per pattern would otherwise be
      // copied again on the way out. Solving happens first, because the
      // transfer detaches the buffer.
      scope.postMessage(
        {
          kind: "frame",
          id: message.id,
          detail,
          report,
          image,
          width: frame.width,
          height: frame.height,
        },
        [image.buffer],
      );
      return;
    }

    const { baseSeed, start, count, chunk, config } = message;
    const seed = BigInt(baseSeed);

    // Chunked so progress arrives while the run is in flight rather than only
    // at the end.
    for (let done = 0; done < count; done += chunk) {
      const size = Math.min(chunk, count - done);
      const rows = tracker.run_benchmark_chunk(
        seed,
        start + done,
        size,
        config,
      ) as TrialRow[];
      scope.postMessage({ kind: "rows", rows });
    }
    scope.postMessage({ kind: "done" });
  } catch (error) {
    scope.postMessage({
      kind: "error",
      ...("id" in message ? { id: message.id } : {}),
      message: error instanceof Error ? error.message : String(error),
    });
  }
};
