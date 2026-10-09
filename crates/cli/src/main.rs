//! Native benchmark runner.
//!
//! ```text
//! cargo run -p tracker-cli --release -- bench --trials 1000 --preset nominal --seed 42
//! ```
//!
//! Prints the report, then writes `bench_<preset>_<seed>.csv` with one row per
//! trial and `bench_<preset>_<seed>.json` with the summary.
//!
//! Arguments are parsed by hand and the JSON is written by hand: neither `clap`
//! nor `serde_json` is on docs/SPEC.md's approved dependency list, and three flags
//! and a flat struct of numbers do not justify asking for them.

use std::error::Error;
use std::fmt::Write as _;
use std::str::FromStr;

use tracker_core::bench::{self, BenchReport, Solver, Spread, TrialResult, Workspaces};
use tracker_core::catalog::Catalog;
use tracker_core::pairdb::PairDb;
use tracker_core::simulate::{Preset, SimConfig};
use tracker_core::track::{self, TrackConfig};
use tracker_core::verify::Outcome;
use tracker_core::{Clock, NoClock, StdClock};

/// The committed catalogue, embedded so the binary needs no data files beside
/// it and cannot disagree with the one the tests ran against.
const CATALOG_BIN: &[u8] = include_bytes!("../../../data/catalog.bin");

const USAGE: &str = "\
usage: tracker-cli <bench|track|parity> [options]

commands:
  bench             run a scored benchmark and report it
  track             run a slew and compare tracking against the full search
  parity            print full-precision trial rows, for the WASM parity check

options:
  --trials N        trials to run (default 1000)
  --frames N        frames in a tracking sequence (track only, default 60)
  --rate DEG        body rate in degrees per second (track only, default 0.5)
  --preset NAME     easy | nominal | hard | brutal (default nominal)
  --seed N          base seed; trial i uses splitmix64(seed ^ i) (default 42)
  --out DIR         where to write the CSV and JSON (default the current directory)
  --no-files        print the report only
  -h, --help        this message
";

fn main() -> Result<(), Box<dyn Error>> {
    match run() {
        Ok(()) => Ok(()),
        Err(message) => {
            eprintln!("error: {message}");
            eprintln!();
            eprint!("{USAGE}");
            std::process::exit(2);
        }
    }
}

/// Parsed command line.
struct Args {
    trials: usize,
    preset: Preset,
    seed: u64,
    out: String,
    write_files: bool,
    /// Frames in a tracking sequence; `track` only.
    frames: usize,
    /// Body rate for a tracking sequence, degrees per second.
    rate_deg_s: f64,
}

fn run() -> Result<(), String> {
    let mut raw = std::env::args().skip(1);
    let command = raw.next().unwrap_or_default();
    if command == "-h" || command == "--help" || command.is_empty() {
        print!("{USAGE}");
        return Ok(());
    }
    if command != "bench" && command != "parity" && command != "track" {
        return Err(format!("unknown command {command:?}"));
    }

    let mut args = Args {
        trials: 1000,
        preset: Preset::Nominal,
        seed: 42,
        out: ".".to_string(),
        write_files: true,
        frames: TrackConfig::default().frames,
        rate_deg_s: TrackConfig::default().rate_deg_s,
    };
    let rest: Vec<String> = raw.collect();
    let mut at = 0;
    while at < rest.len() {
        let flag = rest[at].as_str();
        // Flags that take no value.
        if flag == "--no-files" {
            args.write_files = false;
            at += 1;
            continue;
        }
        if flag == "-h" || flag == "--help" {
            print!("{USAGE}");
            return Ok(());
        }
        let value = rest
            .get(at + 1)
            .ok_or_else(|| format!("{flag} needs a value"))?;
        match flag {
            "--trials" => {
                args.trials = value
                    .parse()
                    .map_err(|_| format!("--trials wants a number, got {value:?}"))?;
            }
            "--preset" => {
                args.preset = Preset::from_str(value).map_err(|e| e.to_string())?;
            }
            "--seed" => {
                args.seed = value
                    .parse()
                    .map_err(|_| format!("--seed wants a number, got {value:?}"))?;
            }
            "--out" => args.out = value.clone(),
            "--frames" => {
                args.frames = value
                    .parse()
                    .map_err(|_| format!("--frames wants a number, got {value:?}"))?;
            }
            "--rate" => {
                args.rate_deg_s = value
                    .parse()
                    .map_err(|_| format!("--rate wants a number, got {value:?}"))?;
            }
            other => return Err(format!("unknown option {other:?}")),
        }
        at += 2;
    }
    if args.trials == 0 {
        return Err("--trials must be at least 1".to_string());
    }

    if command == "track" {
        if args.frames == 0 {
            return Err("--frames must be at least 1".to_string());
        }
        return track_command(&args).map_err(|e| e.to_string());
    }
    if command == "parity" {
        return parity_command(&args).map_err(|e| e.to_string());
    }
    bench_command(&args).map_err(|e| e.to_string())
}

/// Prints one trial per line as JSON, at full `f64` precision.
///
/// The WASM parity check compares against this. The `bench` CSV rounds its
/// floats for reading, which cannot support a 1e-12 comparison, and timings are
/// left out because they are expected to differ between the two targets.
fn parity_command(args: &Args) -> Result<(), Box<dyn Error>> {
    let cfg = SimConfig::preset(args.preset);
    let catalog = Catalog::from_bytes(CATALOG_BIN)?;
    let db = PairDb::build(&catalog, &cfg)?;
    let ids = catalog.id_index();
    let solver = Solver::new(&catalog, &db, &ids);
    let (_, rows) = bench::run_benchmark(args.seed, args.trials, &cfg, &solver, &NoClock)?;

    for row in &rows {
        let [qw, qx, qy, qz] = row.quaternion;
        println!(
            "{{\"index\":{},\"seed\":\"{}\",\"outcome\":\"{}\",\"total_error\":{},\"cross_error\":{},\"roll_error\":{},\"claimed\":{},\"correct\":{},\"available\":{},\"detected\":{},\"used\":{},\"matched\":{},\"residual_rms_px\":{},\"focal_scale\":{},\"quaternion\":[{},{},{},{}]}}",
            row.index,
            row.seed,
            row.outcome,
            number(row.total_error),
            number(row.cross_error),
            number(row.roll_error),
            row.claimed,
            row.correct,
            row.available,
            row.detected,
            row.used,
            row.matched,
            number(row.residual_rms_px),
            number(row.focal_scale),
            number(qw),
            number(qx),
            number(qy),
            number(qz),
        );
    }
    Ok(())
}

fn bench_command(args: &Args) -> Result<(), Box<dyn Error>> {
    let cfg = SimConfig::preset(args.preset);
    let clock = StdClock::new();

    let started = clock.now_ns();
    let catalog = Catalog::from_bytes(CATALOG_BIN)?;
    let db = PairDb::build(&catalog, &cfg)?;
    let ids = catalog.id_index();
    let setup_ms = (clock.now_ns() - started) as f64 / 1e6;

    println!(
        "catalogue {} stars, database {} stars and {} pairs, built in {:.0} ms",
        catalog.stars.len(),
        db.star_count(),
        db.pairs().len(),
        setup_ms
    );
    println!("preset {}", args.preset);
    println!();

    let solver = Solver::new(&catalog, &db, &ids);
    let started = clock.now_ns();
    let (report, rows) = bench::run_benchmark(args.seed, args.trials, &cfg, &solver, &clock)?;
    let wall_s = (clock.now_ns() - started) as f64 / 1e9;

    print!("{}", bench::format_report(&report, &cfg));
    println!();
    println!("  wall clock {wall_s:.2} s for {} trials", args.trials);

    if args.write_files {
        let stem = format!("bench_{}_{}", args.preset, args.seed);
        let csv_path = format!("{}/{stem}.csv", args.out);
        let json_path = format!("{}/{stem}.json", args.out);
        std::fs::write(&csv_path, csv(&rows))?;
        std::fs::write(&json_path, json(&report, &cfg, args.trials))?;
        println!();
        println!("wrote {csv_path}");
        println!("wrote {json_path}");
    }

    // A non-zero exit on the outcome that must never happen, so a CI run or a
    // shell loop notices without parsing the report.
    if report.wrong_confident() > 0 {
        eprintln!();
        eprintln!(
            "error: {} trials were WRONG_CONFIDENT, which must be zero",
            report.wrong_confident()
        );
        std::process::exit(1);
    }
    Ok(())
}

/// An `f64` as JSON.
///
/// An unsolved trial carries NaN errors and an infinite residual, and Rust
/// prints those as `NaN` and `inf`, neither of which a JSON parser accepts.
/// They become `null`, which the parity check reads back as not-a-number.
fn number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:?}")
    } else {
        "null".to_string()
    }
}

/// Every trial as CSV, header first.
fn csv(rows: &[TrialResult]) -> String {
    let mut out = String::with_capacity(rows.len() * 160);
    out.push_str(bench::csv_header());
    out.push('\n');
    for row in rows {
        out.push_str(&bench::csv_row(row));
        out.push('\n');
    }
    out
}

/// The summary as JSON.
///
/// Written by hand because `BenchReport` is a flat struct of numbers and
/// strings, so there is nothing here a serialiser would do better.
fn json(report: &BenchReport, cfg: &SimConfig, trials: usize) -> String {
    let spread = |label: &str, spread: Spread| {
        format!(
            "    \"{label}\": {{ \"median\": {:.6}, \"p95\": {:.6}, \"max\": {:.6}, \
\"mean\": {:.6} }}",
            spread.median, spread.p95, spread.max, spread.mean
        )
    };
    let ms = |s: Spread| s.scaled(1e6);

    let mut out = String::new();
    let _ = writeln!(out, "{{");
    let _ = writeln!(out, "  \"preset\": \"{}\",", preset_name(cfg));
    let _ = writeln!(out, "  \"base_seed\": {},", report.base_seed);
    let _ = writeln!(out, "  \"trials\": {trials},");
    let _ = writeln!(out, "  \"score_percent\": {:.4},", report.score);
    let _ = writeln!(out, "  \"wrong_confident\": {},", report.wrong_confident());
    let _ = writeln!(out, "  \"outcomes\": {{");
    let last = Outcome::ALL.len() - 1;
    for (at, outcome) in Outcome::ALL.iter().enumerate() {
        let _ = writeln!(
            out,
            "    \"{}\": {}{}",
            outcome.as_str(),
            report.count(*outcome),
            if at == last { "" } else { "," }
        );
    }
    let _ = writeln!(out, "  }},");
    let _ = writeln!(out, "  \"id_precision\": {:.6},", report.id_precision);
    let _ = writeln!(out, "  \"id_recall\": {:.6},", report.id_recall);
    let _ = writeln!(out, "  \"error_arcsec\": {{");
    let _ = writeln!(out, "{},", spread("total", report.total_error_arcsec));
    let _ = writeln!(
        out,
        "{},",
        spread("cross_boresight", report.cross_error_arcsec)
    );
    let _ = writeln!(out, "{}", spread("roll", report.roll_error_arcsec));
    let _ = writeln!(out, "  }},");
    let _ = writeln!(out, "  \"stage_ms\": {{");
    let _ = writeln!(out, "{},", spread("simulate", ms(report.timings.simulate)));
    let _ = writeln!(out, "{},", spread("centroid", ms(report.timings.centroid)));
    let _ = writeln!(out, "{},", spread("identify", ms(report.timings.identify)));
    let _ = writeln!(out, "{},", spread("attitude", ms(report.timings.attitude)));
    let _ = writeln!(out, "{},", spread("verify", ms(report.timings.verify)));
    let _ = writeln!(out, "{}", spread("solve", ms(report.timings.solve)));
    let _ = writeln!(out, "  }}");
    let _ = write!(out, "}}");
    out
}

/// The preset a configuration came from, by matching it against the four.
fn preset_name(cfg: &SimConfig) -> &'static str {
    Preset::ALL
        .into_iter()
        .find(|preset| SimConfig::preset(*preset) == *cfg)
        .map_or("custom", Preset::as_str)
}

/// Runs a slew and reports tracking against the full search.
///
/// Both run on every frame from the same directions, so the comparison is of
/// two identifications of one image rather than of two separate runs.
fn track_command(args: &Args) -> Result<(), Box<dyn Error>> {
    let catalog = Catalog::from_bytes(CATALOG_BIN)?;
    let cfg = SimConfig::preset(args.preset);
    let db = PairDb::build(&catalog, &cfg)?;
    let ids = catalog.id_index();
    let solver = Solver::new(&catalog, &db, &ids);
    let track_cfg = TrackConfig {
        frames: args.frames,
        rate_deg_s: args.rate_deg_s,
        ..TrackConfig::default()
    };
    let clock = StdClock::default();
    let mut ws = Workspaces::new();

    println!(
        "catalogue {} stars, database {} stars and {} pairs",
        catalog.stars.len(),
        db.star_count(),
        db.pairs().len(),
    );
    println!("preset {}", args.preset);
    println!();
    println!(
        "{} frames  base seed {}  {} deg/s  {} s apart  step {:.3} deg",
        track_cfg.frames,
        args.seed,
        track_cfg.rate_deg_s,
        track_cfg.interval_s,
        track_cfg.step_rad().to_degrees(),
    );

    let steps = track::run_sequence(args.seed, &cfg, &track_cfg, &solver, &clock, &mut ws);
    let report = track::summarise(&steps);

    const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);
    println!();
    println!(
        "  CORRECT {} of {}   WRONG_CONFIDENT {}   reacquisitions {}",
        report.correct, report.frames, report.wrong_confident, report.reacquisitions,
    );
    println!(
        "  median total error {:.2} arcsec",
        report.median_error / ARCSEC
    );
    println!();
    println!("  identify (median)      tracked      lost-in-space     speedup");
    println!(
        "                      {:>8.3} ms      {:>8.3} ms      {:>6.1}x",
        report.track_ns as f64 / 1.0e6,
        report.lost_ns as f64 / 1.0e6,
        report.speedup(),
    );
    println!();
    println!("  The first two frames are searched in full: a rate needs two attitudes");
    println!("  before it can be predicted, so tracking begins at the third.");

    if args.write_files {
        let path = format!("{}/track_{}_{}.csv", args.out, args.preset, args.seed);
        let mut csv = String::from(
            "index,outcome,aided,reacquired,detected,claimed,correct,matched,total_arcsec,track_ns,lost_ns
",
        );
        for step in &steps {
            csv.push_str(&format!(
                "{},{:?},{},{},{},{},{},{},{:.6},{},{}
",
                step.index,
                step.outcome,
                step.aided,
                step.reacquired,
                step.detected,
                step.claimed,
                step.correct,
                step.matched,
                step.total_error / ARCSEC,
                step.track_ns,
                step.lost_ns,
            ));
        }
        std::fs::write(&path, csv)?;
        println!();
        println!("wrote {path}");
    }
    Ok(())
}
