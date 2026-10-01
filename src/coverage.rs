//! Coverage track construction.
//!
//! Two input paths are supported:
//!   1. A BAM file, read directly with `noodles` (pure Rust — no htslib/C
//!      dependency, no `samtools` binary required on the host).
//!   2. A plain depth/bedGraph-style TSV (`contig<TAB>pos<TAB>depth`, as
//!      produced by `samtools depth`), for pipelines that already have one.
//!
//! Both paths produce the same output: mean depth per fixed-size window per
//! contig, ready to feed straight into the coverage ring of the plot.

use anyhow::{Context, Result};
use noodles::bam;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

/// Per-contig binned mean depth: window_index -> mean depth.
pub type DepthTrack = HashMap<String, Vec<f64>>;

/// Bin BAM alignments into `window`-sized bins per reference sequence,
/// approximating per-base coverage from each alignment's CIGAR-derived
/// reference span (so soft-clips/introns/deletions are handled correctly,
/// unlike a naive "read length" estimate).
pub fn depth_from_bam(
    path: &Path,
    window: u64,
    contig_lengths: &HashMap<String, u64>,
) -> Result<DepthTrack> {
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(path)
        .with_context(|| format!("opening BAM {:?}", path))?;
    let header = reader.read_header().context("reading BAM header")?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect();

    let window = window.max(1);
    let mut bins: HashMap<String, Vec<f64>> = contig_lengths
        .iter()
        .map(|(id, len)| {
            (
                id.clone(),
                vec![0.0f64; ((*len + window - 1) / window) as usize],
            )
        })
        .collect();

    for result in reader.records() {
        let record = result.context("reading BAM record")?;

        let Some(Ok(ref_id)) = record.reference_sequence_id() else {
            continue; // unmapped
        };
        let Some(name) = ref_names.get(ref_id) else {
            continue;
        };
        let Some(bin_vec) = bins.get_mut(name) else {
            continue; // reference not in the FASTA we were given
        };
        let Some(Ok(start)) = record.alignment_start() else {
            continue;
        };
        let start0 = usize::from(start) - 1; // noodles Position is 1-based

        // Walk the CIGAR, accumulating +1 depth per reference base covered
        // by M/D/N/=/X ops (i.e. reference-consuming ops), distributed into
        // the appropriate window bins.
        let mut ref_pos = start0;
        for op in record.cigar().iter() {
            let op = op.context("reading CIGAR op")?;
            if op.kind().consumes_reference() {
                let len = op.len();
                for p in ref_pos..ref_pos + len {
                    let bin_idx = p / window as usize;
                    if let Some(slot) = bin_vec.get_mut(bin_idx) {
                        *slot += 1.0;
                    }
                }
                ref_pos += len;
            }
        }
    }

    // Convert summed per-base depth into mean depth per window.
    for (id, bin_vec) in bins.iter_mut() {
        let len = *contig_lengths.get(id).unwrap_or(&0);
        for (i, v) in bin_vec.iter_mut().enumerate() {
            let bin_start = i as u64 * window;
            let bin_end = (bin_start + window).min(len);
            let bin_len = (bin_end - bin_start).max(1) as f64;
            *v /= bin_len;
        }
    }

    Ok(bins)
}

/// Parse a `samtools depth`-style TSV (`contig\tpos(1-based)\tdepth`) into
/// the same binned representation as [`depth_from_bam`].
pub fn depth_from_tsv(
    path: &Path,
    window: u64,
    contig_lengths: &HashMap<String, u64>,
) -> Result<DepthTrack> {
    let file = File::open(path).with_context(|| format!("opening depth TSV {:?}", path))?;
    let reader = BufReader::new(file);
    let window = window.max(1);

    let mut sums: HashMap<String, Vec<f64>> = contig_lengths
        .iter()
        .map(|(id, len)| {
            (
                id.clone(),
                vec![0.0f64; ((*len + window - 1) / window) as usize],
            )
        })
        .collect();
    let mut counts: HashMap<String, Vec<u64>> = contig_lengths
        .iter()
        .map(|(id, len)| {
            (
                id.clone(),
                vec![0u64; ((*len + window - 1) / window) as usize],
            )
        })
        .collect();

    for line in reader.lines() {
        let line = line?;
        let mut cols = line.split('\t');
        let (Some(contig), Some(pos), Some(depth)) = (cols.next(), cols.next(), cols.next()) else {
            continue;
        };
        let (Ok(pos), Ok(depth)) = (pos.parse::<u64>(), depth.parse::<f64>()) else {
            continue;
        };
        let Some(sum_vec) = sums.get_mut(contig) else {
            continue;
        };
        let bin_idx = ((pos.saturating_sub(1)) / window) as usize;
        if let Some(slot) = sum_vec.get_mut(bin_idx) {
            *slot += depth;
            counts.get_mut(contig).unwrap()[bin_idx] += 1;
        }
    }

    for (id, sum_vec) in sums.iter_mut() {
        let count_vec = &counts[id];
        for (v, c) in sum_vec.iter_mut().zip(count_vec.iter()) {
            if *c > 0 {
                *v /= *c as f64;
            }
        }
    }

    Ok(sums)
}
