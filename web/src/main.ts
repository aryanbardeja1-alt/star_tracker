// Page bootstrap: compile the WASM module once, wire the tabs, and start both
// views. The module is shared with the benchmark workers, so it is compiled a
// single time no matter how many of them run.

import { BenchmarkTab } from "./ui/benchmark.js";
import { PatternTab } from "./ui/pattern.js";
import { TrackingTab } from "./ui/tracking.js";
import { Engine } from "./engine.js";
import type { SimConfig } from "./wasm.js";

function need<T extends HTMLElement>(id: string): T {
  const element = document.getElementById(id);
  if (!element) throw new Error(`missing element #${id}`);
  return element as T;
}

/** The panels the nav switches between. */
type Tab = "pattern" | "bench" | "track";

/** Switches the panels, and returns a way to select one by name. */
function wireTabs(): (which: Tab) => void {
  const tabs = {
    pattern: [need("tab-pattern"), need("view-pattern")] as const,
    bench: [need("tab-bench"), need("view-bench")] as const,
    track: [need("tab-track"), need("view-track")] as const,
  };
  const select = (which: Tab) => {
    for (const [name, [button, panel]] of Object.entries(tabs)) {
      const chosen = name === which;
      button.setAttribute("aria-selected", String(chosen));
      panel.hidden = !chosen;
    }
  };
  for (const [name, [button]] of Object.entries(tabs)) {
    button.addEventListener("click", () => select(name as Tab));
  }
  return select;
}

async function start(): Promise<void> {
  const selectTab = wireTabs();
  const build = need("build");
  build.textContent = "loading WASM…";

  // Nothing is instantiated on this thread: the engine owns the tracker on a
  // worker, so the page stays responsive while the catalogue, the pair
  // database and the first pattern are built.
  const engine = new Engine();
  const info = await engine.info;
  build.textContent =
    `v${info.version} · ${info.catalogStars.toLocaleString()} stars · ` +
    `${info.databasePairs.toLocaleString()} pairs`;

  const configFor = (preset: string): Promise<SimConfig> => engine.config(preset);
  const defaults = await configFor("nominal");

  const pattern = new PatternTab(engine, defaults);
  new TrackingTab(engine, configFor);
  const benchmark = new BenchmarkTab(defaults, configFor, (seed, index, preset) => {
    // Clicking a row has to land on the pattern view, or the reproduction it
    // just set up would be invisible.
    selectTab("pattern");
    void pattern.open(seed, index, preset);
  });

  // A particular pattern can be addressed by URL, which is the same operation
  // clicking a benchmark row performs: the trial is identified entirely by its
  // base seed, its index and the preset.
  const params = new URLSearchParams(window.location.search);
  const seed = params.get("seed");
  const trial = params.get("trial");
  const preset = params.get("preset");
  if (seed !== null || trial !== null || preset !== null) {
    await pattern.open(seed ?? "42", Number(trial ?? 0), preset ?? "nominal");
  } else {
    await pattern.run(true);
  }

  // And a benchmark run can be started straight from the URL, which is how the
  // wall-clock budget is checked.
  if (params.get("autorun")) {
    selectTab("bench");
    await benchmark.autorun(params);
  }
}

start().catch((error: unknown) => {
  const build = document.getElementById("build");
  const message = error instanceof Error ? error.message : String(error);
  if (build) build.textContent = `failed: ${message}`;
  console.error(error);
});
