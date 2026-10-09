//! Star tracker core: image -> centroids -> unit vectors -> identification ->
//! attitude -> verification.
//!
//! A pure library: no I/O and no printing, so the same code runs behind the
//! native CLI and behind wasm-bindgen. Frames, units and conventions are
//! defined in docs/SPEC.md and restated in each module.

#![forbid(unsafe_code)]

pub mod attitude;
pub mod bench;
pub mod camera;
pub mod catalog;
pub mod centroid;
pub mod identify;
pub mod math;
pub mod pairdb;
pub mod simulate;
pub mod track;
pub mod verify;

/// Anything that can go wrong inside `tracker-core`.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// The catalogue bytes could not be decoded as a `Catalog`.
    #[error("catalogue decode failed: {0}")]
    CatalogDecode(String),
    /// The catalogue could not be encoded.
    #[error("catalogue encode failed: {0}")]
    CatalogEncode(String),
    /// The decoded catalogue holds no stars.
    #[error("catalogue contains no stars")]
    CatalogEmpty,
    /// A catalogue entry's direction is not a finite unit vector, which means
    /// the file is stale or corrupt.
    #[error("catalogue star {id} has a non-unit direction (norm {norm})")]
    CatalogNotUnit {
        /// Hipparcos number of the offending entry.
        id: u32,
        /// The norm actually found.
        norm: f64,
    },
    /// A preset name from the CLI or the web UI is not one of the four.
    #[error("unknown preset: {0}")]
    UnknownPreset(String),
    /// A benchmark was asked for zero trials.
    #[error("a benchmark needs at least one trial")]
    EmptyBenchmark,
    /// The magnitude cut admitted more stars than a u16 pair index can address.
    #[error("database would hold {count} stars, more than a u16 index allows")]
    TooManyDatabaseStars {
        /// How many stars the cut admitted.
        count: usize,
    },
}

/// A monotonic clock supplied by the host.
///
/// `std::time::Instant::now()` compiles for `wasm32-unknown-unknown` and then
/// traps at run time, so `tracker-core` never reads a clock of its own. The
/// native side passes [`StdClock`]; the WASM wrapper passes one backed by
/// `performance.now()`.
pub trait Clock {
    /// Nanoseconds since some fixed, arbitrary origin.
    fn now_ns(&self) -> u64;
}

/// A clock that always reads zero, for callers that want no timings.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoClock;

impl Clock for NoClock {
    fn now_ns(&self) -> u64 {
        0
    }
}

/// A clock backed by `std::time::Instant`. Native targets only.
#[cfg(not(target_arch = "wasm32"))]
#[derive(Clone, Debug)]
pub struct StdClock {
    origin: std::time::Instant,
}

#[cfg(not(target_arch = "wasm32"))]
impl StdClock {
    /// Starts a clock whose origin is now.
    pub fn new() -> Self {
        Self {
            origin: std::time::Instant::now(),
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Default for StdClock {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Clock for StdClock {
    fn now_ns(&self) -> u64 {
        self.origin.elapsed().as_nanos() as u64
    }
}

/// `Result` specialised to [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
