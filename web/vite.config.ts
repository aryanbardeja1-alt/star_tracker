import { defineConfig } from "vite";

export default defineConfig({
  // GitHub Pages serves the site from a subpath; a relative base keeps the
  // built asset URLs correct wherever it is mounted.
  base: "./",
  build: { target: "es2022" },
});
