//! Pair-database build and query timing. The build has a 300 ms budget in
//! WASM; `cargo bench -p tracker-core` reports the native figure, and the
//! Phase 4 log records the measured WASM one.

use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use tracker_core::catalog::Catalog;
use tracker_core::pairdb::PairDb;
use tracker_core::simulate::{Preset, SimConfig};

/// The committed catalogue, embedded so the benchmark does no I/O.
const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

/// One arcsecond in radians.
const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);

fn pair_database(c: &mut Criterion) {
    let catalog = Catalog::from_bytes(CATALOG_BIN).expect("catalog.bin must decode");
    let cfg = SimConfig::preset(Preset::Nominal);

    c.bench_function("pairdb/build", |b| {
        b.iter(|| {
            let db = PairDb::build(black_box(&catalog), &cfg).expect("build");
            black_box(db.pairs().len())
        });
    });

    let db = PairDb::build(&catalog, &cfg).expect("build");

    // A spread of separations across the table, so the measurement is not
    // dominated by one corner of the angle distribution.
    let angles: Vec<f64> = (0..64)
        .map(|step| db.max_angle() * (step as f64 + 0.5) / 64.0)
        .collect();

    let mut group = c.benchmark_group("pairdb/query");
    for eps_arcsec in [5.0, 20.0, 60.0] {
        let eps = eps_arcsec * ARCSEC;
        group.bench_function(format!("{eps_arcsec:.0}arcsec"), |b| {
            b.iter(|| {
                let mut total = 0usize;
                for &angle in &angles {
                    total += db.query(black_box(angle), eps).len();
                }
                black_box(total)
            });
        });
    }
    group.finish();
}

criterion_group!(benches, pair_database);
criterion_main!(benches);
