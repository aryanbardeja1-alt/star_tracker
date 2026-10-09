//! Per-frame timing for the centroider, against the 2 ms budget in the build
//! plan. Run with `cargo bench -p tracker-core`.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use tracker_core::catalog::Catalog;
use tracker_core::centroid::{self, Workspace as CentroidWorkspace};
use tracker_core::simulate::{self, Preset, SimConfig, Workspace as SimWorkspace};

/// The committed catalogue, embedded so the benchmark does no I/O.
const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

fn centroiding(c: &mut Criterion) {
    let catalog = Catalog::from_bytes(CATALOG_BIN).expect("catalog.bin must decode");
    let mut group = c.benchmark_group("centroid/1024x1024");

    for preset in Preset::ALL {
        let cfg = SimConfig::preset(preset);
        // One representative frame per preset, rendered up front: this measures
        // detection only, not simulation.
        let frame = simulate::simulate(42, &cfg, &catalog, &mut SimWorkspace::new());
        let mut ws = CentroidWorkspace::new();
        centroid::detect(&frame.image, frame.width, frame.height, &cfg, &mut ws);

        group.bench_function(preset.as_str(), |b| {
            b.iter(|| {
                let found = centroid::detect(
                    black_box(&frame.image),
                    frame.width,
                    frame.height,
                    &cfg,
                    &mut ws,
                );
                black_box(found.len())
            });
        });
    }
    group.finish();
}

criterion_group!(benches, centroiding);
criterion_main!(benches);
