// The interactive view's compute, kept off the main thread.
//
// One worker owns the only `Tracker` the Single Pattern tab needs. The main
// thread compiles the module, hands it over, and from then on does nothing but
// draw what comes back: building a tracker and simulating a frame together cost
// about 300 ms, which on the main thread is 300 ms of a page that cannot answer
// a click.
//
// The benchmark tab spawns its own pool from the same worker script; this is a
// separate instance, so a pattern can be opened while a run is in flight.

import {
  compiledModule,
  type FrameResult,
  type FromWorker,
  type SimConfig,
  type TrackConfig,
  type TrackRun,
} from "./wasm.js";

/** What the worker reports about the build once it is ready. */
export interface EngineInfo {
  version: string;
  catalogStars: number;
  databasePairs: number;
}

interface Waiting {
  resolve: (value: never) => void;
  reject: (error: Error) => void;
}

/** The handshake's reply is awaited under this id; requests start at 1. */
const HANDSHAKE = 0;

export class Engine {
  private readonly worker: Worker;
  private readonly waiting = new Map<number, Waiting>();
  private next = HANDSHAKE + 1;

  /** Resolves once the worker's tracker exists. */
  readonly info: Promise<EngineInfo>;

  constructor() {
    this.worker = new Worker(new URL("./worker.ts", import.meta.url), {
      type: "module",
    });
    this.worker.onmessage = (event: MessageEvent<FromWorker>) => {
      this.settle(event.data);
    };
    this.info = new Promise<EngineInfo>((resolve, reject) => {
      this.waiting.set(HANDSHAKE, { resolve: resolve as (value: never) => void, reject });
    });
    this.worker.onerror = () => {
      const failed = new Error("the compute worker failed to start");
      for (const [id, waiting] of this.waiting) {
        this.waiting.delete(id);
        waiting.reject(failed);
      }
    };
    void compiledModule().then((module) => {
      this.worker.postMessage({ kind: "init", module });
    });
  }

  /** The default configuration for a preset, as the UI's editable object. */
  async config(preset: string): Promise<SimConfig> {
    await this.info;
    const reply = await this.request<{ config: SimConfig }>((id) => ({
      kind: "config" as const,
      id,
      preset,
    }));
    return reply.config;
  }

  /**
   * Scores trial `index` of `baseSeed` and returns it with the frame it drew.
   *
   * The pixels arrive as a transferred buffer; `image` is row-major
   * `width * height` samples.
   */
  async frame(
    baseSeed: string,
    index: number,
    config: SimConfig,
  ): Promise<FrameResult> {
    await this.info;
    return this.request<FrameResult>((id) => ({
      kind: "frame" as const,
      id,
      baseSeed,
      index,
      config,
    }));
  }

  /**
   * Runs a slew and returns every frame with the figures that summarise it.
   *
   * `track` may be null, which takes the Rust-side defaults.
   */
  async tracking(
    baseSeed: string,
    config: SimConfig,
    track: TrackConfig | null,
  ): Promise<TrackRun> {
    await this.info;
    const reply = await this.request<{ run: TrackRun }>((id) => ({
      kind: "track" as const,
      id,
      baseSeed,
      config,
      track,
    }));
    return reply.run;
  }

  private request<T>(build: (id: number) => object): Promise<T> {
    const id = this.next++;
    return new Promise<T>((resolve, reject) => {
      this.waiting.set(id, { resolve: resolve as (value: never) => void, reject });
      this.worker.postMessage(build(id));
    });
  }

  private settle(message: FromWorker): void {
    const id = "id" in message && message.id !== undefined ? message.id : HANDSHAKE;
    const waiting = this.waiting.get(id);
    if (!waiting) return;
    this.waiting.delete(id);
    if (message.kind === "error") waiting.reject(new Error(message.message));
    else waiting.resolve(message as never);
  }
}
