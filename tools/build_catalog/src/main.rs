//! Turns the raw Hipparcos catalogue into `data/catalog.bin`.
//!
//! Source is Hipparcos (VizieR I/239) rather than the Yale BSC5: BSC5 stores
//! J2000 positions rounded to 0.1s in RA and 1" in Dec, which alone can put a
//! star ~0.9" from its true place, leaving no margin against the 1" accuracy
//! the build plan requires. Hipparcos gives milliarcsecond positions, but at
//! epoch J1991.25, so proper motion is propagated to J2000.0 here.
//!
//! Proper names come from the IAU WGSN Catalog of Star Names, joined on HIP.
//!
//! Angles are radians and vectors are ICRS/J2000 throughout, per docs/SPEC.md.

use std::error::Error;
use std::f64::consts::PI;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use tracker_core::catalog::{Catalog, Star};
use tracker_core::math::{Vec3, angle_between};

/// Faintest V magnitude kept, matching the simulator's limiting magnitude.
const MAG_LIMIT: f32 = 6.5;

/// Default camera geometry from docs/SPEC.md, used only to size the blend radius.
const FOV_RAD: f64 = 20.0 * PI / 180.0;
const IMAGE_WIDTH_PX: f64 = 1024.0;

/// Separation, in pixels, below which two stars are merged into one entry.
///
/// The build plan's starting value was 2 px, but Phase 3 measured what the
/// detector actually does: with a 1 px PSF, 8-connected labelling joins two
/// stars into a single blob out to 6-8 px, because the saddle between their
/// peaks stays above the detection threshold. At 2 px the catalogue therefore
/// claimed two stars wherever the sensor can only ever produce one centroid,
/// and that centroid matched neither entry -- it behaved like a false star.
/// Six pixels is about 2.5 PSF full-widths and merges only 66 pairs out of
/// 8827 stars, while lifting detection completeness from 97.9% to 99.7%.
const BLEND_PX: f64 = 6.0;

/// That separation as an angle: `BLEND_PX * FOV / W`.
const BLEND_RADIUS: f64 = BLEND_PX * FOV_RAD / IMAGE_WIDTH_PX;

/// Hipparcos positions are given at this epoch; the catalogue stores J2000.0.
const HIP_EPOCH: f64 = 1991.25;
const TARGET_EPOCH: f64 = 2000.0;

/// One arcsecond in radians, for reporting only.
const ARCSEC: f64 = PI / (180.0 * 3600.0);

/// Milliarcseconds to radians.
const MAS: f64 = ARCSEC / 1000.0;

const HIP_URL: &str = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?\
-source=I/239/hip_main&-out=HIP,RAICRS,DEICRS,pmRA,pmDE,Vmag&-out.max=99999&Vmag=%3C=6.5";
const NAMES_URL: &str = "https://www.pas.rochester.edu/~emamajek/WGSN/IAU-CSN.txt";

/// Paths relative to the workspace root, used both to locate and to report
/// the files; the absolute form is built from `CARGO_MANIFEST_DIR` and reads
/// badly in logs.
const HIP_FILE: &str = "data/raw/hip_main_v65.tsv";
const NAMES_FILE: &str = "data/raw/IAU-CSN.txt";
const OUT_FILE: &str = "data/catalog.bin";

/// A parsed Hipparcos row, already propagated to J2000.0.
struct Entry {
    id: u32,
    unit: Vec3,
    mag: f32,
}

fn main() -> Result<(), Box<dyn Error>> {
    let root = PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."));
    let hip_path = root.join(HIP_FILE);
    let names_path = root.join(NAMES_FILE);
    if let Some(parent) = hip_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    ensure_downloaded(&hip_path, HIP_URL, HIP_FILE)?;
    ensure_downloaded(&names_path, NAMES_URL, NAMES_FILE)?;

    let hip_text = std::fs::read_to_string(&hip_path)?;
    let names_text = std::fs::read_to_string(&names_path)?;

    let (entries, stats) = parse_hipparcos(&hip_text)?;
    println!(
        "parsed    {:>6} rows, dropped {} without position/magnitude, \
         {} without proper motion",
        entries.len(),
        stats.dropped,
        stats.no_proper_motion
    );
    println!(
        "kept      {:>6} with V <= {MAG_LIMIT:.1} (propagated J{HIP_EPOCH} -> J{TARGET_EPOCH})",
        entries.len()
    );

    let names = parse_names(&names_text);
    println!(
        "names     {:>6} IAU proper names with a HIP number",
        names.len()
    );

    let stars = merge_blends(entries, &names);
    println!(
        "merged    {:>6} entries into {} after blending below {:.1}\"",
        stats.kept,
        stars.len(),
        BLEND_RADIUS / ARCSEC
    );

    let catalog = Catalog { stars };
    let bytes = catalog.to_bytes()?;
    std::fs::write(root.join(OUT_FILE), &bytes)?;

    report(&catalog, &bytes);
    Ok(())
}

/// Fetches `url` to `path` with curl unless the file is already present.
///
/// curl is used rather than an HTTP crate so the workspace needs no TLS stack
/// for a tool that runs once per catalogue refresh. If curl is missing, the
/// user is told to place the file by hand, as the build plan allows. `label`
/// is the repository-relative path, for readable messages.
fn ensure_downloaded(path: &Path, url: &str, label: &str) -> Result<(), Box<dyn Error>> {
    if path.metadata().is_ok_and(|m| m.len() > 0) {
        println!("cached    {label}");
        return Ok(());
    }

    println!("download  {label}");
    let status = Command::new("curl")
        .args(["-fsSL", "--max-time", "300", "-o"])
        .arg(path)
        .arg(url)
        .status();

    match status {
        Ok(status) if status.success() => Ok(()),
        Ok(status) => {
            // A partial file would otherwise be treated as a valid cache.
            let _ = std::fs::remove_file(path);
            Err(
                format!("curl exited with {status}\n  url: {url}\n  save it as {label} by hand")
                    .into(),
            )
        }
        Err(e) => Err(format!(
            "could not run curl ({e})\n  url: {url}\n  save it as {label} by hand"
        )
        .into()),
    }
}

#[derive(Default)]
struct ParseStats {
    dropped: usize,
    no_proper_motion: usize,
    kept: usize,
}

/// Parses a VizieR ASU-TSV export of I/239/hip_main.
///
/// The layout is a `#`-commented preamble, a column-name line, a units line, a
/// line of dashes, then tab-separated rows. Columns are located by name rather
/// than position so a change in the query's column order cannot go unnoticed.
fn parse_hipparcos(text: &str) -> Result<(Vec<Entry>, ParseStats), Box<dyn Error>> {
    let lines: Vec<&str> = text.lines().collect();
    let separator = lines
        .iter()
        .position(|l| l.starts_with("---"))
        .ok_or("no dashed separator line found; is this a VizieR TSV export?")?;
    let header = lines
        .get(separator.saturating_sub(2))
        .ok_or("TSV header line missing")?;

    let columns: Vec<&str> = header.split('\t').map(str::trim).collect();
    let index = |name: &str| -> Result<usize, String> {
        columns
            .iter()
            .position(|c| *c == name)
            .ok_or_else(|| format!("column {name} missing; found {columns:?}"))
    };
    let (c_hip, c_ra, c_de) = (index("HIP")?, index("RAICRS")?, index("DEICRS")?);
    let (c_pmra, c_pmde, c_vmag) = (index("pmRA")?, index("pmDE")?, index("Vmag")?);

    let dt = TARGET_EPOCH - HIP_EPOCH;
    let mut stats = ParseStats::default();
    let mut entries = Vec::new();

    for line in &lines[separator + 1..] {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').map(str::trim).collect();

        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let (Ok(id), Ok(ra_deg), Ok(dec_deg), Ok(mag)) = (
            get(c_hip).parse::<u32>(),
            get(c_ra).parse::<f64>(),
            get(c_de).parse::<f64>(),
            get(c_vmag).parse::<f32>(),
        ) else {
            stats.dropped += 1;
            continue;
        };
        if mag > MAG_LIMIT {
            stats.dropped += 1;
            continue;
        }

        // Missing proper motion is rare; treating it as zero is better than
        // discarding an otherwise good star.
        let pm_ra = get(c_pmra).parse::<f64>();
        let pm_de = get(c_pmde).parse::<f64>();
        if pm_ra.is_err() || pm_de.is_err() {
            stats.no_proper_motion += 1;
        }
        let pm_ra = pm_ra.unwrap_or(0.0);
        let pm_de = pm_de.unwrap_or(0.0);

        let unit = propagate(
            ra_deg.to_radians(),
            dec_deg.to_radians(),
            pm_ra * MAS,
            pm_de * MAS,
            dt,
        );
        entries.push(Entry { id, unit, mag });
        stats.kept += 1;
    }

    if entries.is_empty() {
        return Err("no usable rows parsed".into());
    }
    Ok((entries, stats))
}

/// Propagates a position by proper motion, in vector form.
///
/// `ra`/`dec` are radians at the catalogue epoch; `pm_ra` is the *projected*
/// motion mu_alpha.cos(delta) and `pm_de` is mu_delta, both radians per year;
/// `dt` is years. Working on the unit sphere rather than in spherical
/// coordinates keeps the result well behaved near the poles, where dividing
/// the RA rate by cos(delta) blows up -- Polaris moves 30" in RA over these
/// 8.75 years.
fn propagate(ra: f64, dec: f64, pm_ra: f64, pm_de: f64, dt: f64) -> Vec3 {
    let (sin_ra, cos_ra) = ra.sin_cos();
    let (sin_dec, cos_dec) = dec.sin_cos();

    let position = Vec3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec);
    // Unit vectors pointing east and north at this position.
    let east = Vec3::new(-sin_ra, cos_ra, 0.0);
    let north = Vec3::new(-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec);

    (position + (east * pm_ra + north * pm_de) * dt).normalize()
}

/// Parses the IAU-CSN name list into `(HIP, name)` pairs, sorted by HIP.
///
/// The file is fixed-width and holds UTF-8 in its Bayer columns, so byte
/// offsets shift from row to row; the ASCII name is taken from the first 18
/// *characters* and HIP from the first whitespace token after character 89,
/// beyond which every column is plain ASCII.
fn parse_names(text: &str) -> Vec<(u32, String)> {
    let mut names: Vec<(u32, String)> = Vec::new();

    for line in text.lines() {
        if line.starts_with('#') || line.starts_with('$') {
            continue;
        }
        let chars: Vec<char> = line.chars().collect();
        if chars.len() < 110 {
            continue;
        }
        let name: String = chars[..18].iter().collect();
        let name = name.trim();
        let tail: String = chars[89..].iter().collect();
        let Some(Ok(hip)) = tail.split_whitespace().next().map(str::parse::<u32>) else {
            continue;
        };
        if !name.is_empty() {
            names.push((hip, name.to_string()));
        }
    }

    // Keep the first name for each HIP, in HIP order, so the join is
    // deterministic regardless of the file's ordering.
    names.sort_by_key(|(hip, _)| *hip);
    names.dedup_by_key(|(hip, _)| *hip);
    names
}

/// Disjoint-set forest over entry indices, for grouping blended stars.
struct DisjointSet {
    parent: Vec<usize>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            // Path halving keeps this near-constant without recursion.
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// Merges entries closer than [`BLEND_RADIUS`] and returns the catalogue,
/// sorted brightest first.
///
/// Merging is transitive: a chain of close stars becomes one entry, which is
/// what the detector would actually see.
fn merge_blends(mut entries: Vec<Entry>, names: &[(u32, String)]) -> Vec<Star> {
    // Sorting by z lets the pair scan stop early: for unit vectors
    // |sin d_b - sin d_a| <= |d_b - d_a| <= angle(a, b), so once the z gap
    // exceeds the radius, no later entry can be within it either.
    entries.sort_by(|a, b| a.unit.z.total_cmp(&b.unit.z));

    let mut sets = DisjointSet::new(entries.len());
    for (i, a) in entries.iter().enumerate() {
        for (offset, b) in entries[i + 1..].iter().enumerate() {
            if b.unit.z - a.unit.z > BLEND_RADIUS {
                break;
            }
            if angle_between(&a.unit, &b.unit) <= BLEND_RADIUS {
                sets.union(i, i + 1 + offset);
            }
        }
    }

    // Group members by root, in index order, so the result is deterministic.
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); entries.len()];
    for i in 0..entries.len() {
        let root = sets.find(i);
        groups[root].push(i);
    }

    let mut stars: Vec<Star> = groups
        .iter()
        .filter(|members| !members.is_empty())
        .map(|members| combine(members, &entries, names))
        .collect();

    // Brightest first, ties by id, so a magnitude-limited subset is a prefix.
    stars.sort_by(|a, b| a.mag.total_cmp(&b.mag).then(a.id.cmp(&b.id)));
    stars
}

/// Combines one blend group into a single catalogue entry.
///
/// Flux adds, the direction is the flux-weighted mean, and the brightest
/// component supplies the id and name -- that is the star a human would call
/// the object by, and the one whose light dominates the blob.
fn combine(members: &[usize], entries: &[Entry], names: &[(u32, String)]) -> Star {
    let flux = |mag: f32| 10f64.powf(-0.4 * f64::from(mag));

    let mut total_flux = 0.0;
    let mut direction = Vec3::zeros();
    for &i in members {
        let f = flux(entries[i].mag);
        total_flux += f;
        direction += entries[i].unit * f;
    }

    let brightest = members
        .iter()
        .copied()
        .reduce(|a, b| {
            if (entries[b].mag, entries[b].id) < (entries[a].mag, entries[a].id) {
                b
            } else {
                a
            }
        })
        .unwrap_or(0);
    let id = entries[brightest].id;

    let unit = direction.normalize();
    Star {
        id,
        unit: [unit.x, unit.y, unit.z],
        mag: (-2.5 * total_flux.log10()) as f32,
        name: names
            .binary_search_by_key(&id, |(hip, _)| *hip)
            .ok()
            .map(|i| names[i].1.clone()),
    }
}

/// Prints the summary the build plan asks for, plus the size budget.
fn report(catalog: &Catalog, bytes: &[u8]) {
    let stars = &catalog.stars;
    let named = stars.iter().filter(|s| s.name.is_some()).count();

    // Closest surviving pair, by the same z-pruned scan as the merge.
    let mut by_z: Vec<&Star> = stars.iter().collect();
    by_z.sort_by(|a, b| a.unit[2].total_cmp(&b.unit[2]));
    let mut closest = f64::INFINITY;
    for (i, a) in by_z.iter().enumerate() {
        for b in &by_z[i + 1..] {
            if b.unit[2] - a.unit[2] > BLEND_RADIUS {
                break;
            }
            closest = closest.min(angle_between(&a.direction(), &b.direction()));
        }
    }

    let brightest = stars.first();
    let faintest = stars.last();
    let mut line = String::new();
    if let (Some(b), Some(f)) = (brightest, faintest) {
        let _ = write!(
            line,
            "V {:.2} ({}) .. {:.2}",
            b.mag,
            b.name.as_deref().unwrap_or("unnamed"),
            f.mag
        );
    }

    println!();
    println!("catalogue {:>6} stars, {named} named", stars.len());
    println!("          magnitude range {line}");
    println!("          closest pair    {:.1}\"", closest / ARCSEC);
    println!(
        "written   {OUT_FILE} ({} bytes, {:.0}% of the 1 MB budget)",
        bytes.len(),
        100.0 * bytes.len() as f64 / 1_000_000.0
    );
}
