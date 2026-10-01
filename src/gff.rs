//! Minimal GFF3 reader.
//!
//! We only need enough of GFF3 to drive a circular plot: seqid, feature
//! type, start/end (1-based, inclusive, per spec), and strand. Attributes
//! are kept as a raw string in case a caller wants `Name=`/`gene=` later.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone)]
pub struct Feature {
    pub seqid: String,
    pub feature_type: String,
    /// 0-based half-open, converted from GFF3's 1-based inclusive coords so
    /// it lines up directly with the FASTA/BAM coordinate space used
    /// elsewhere in this crate.
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub attributes: String,
}

pub fn read_gff3(path: &Path) -> Result<Vec<Feature>> {
    let file = File::open(path).with_context(|| format!("opening GFF3 {:?}", path))?;
    let reader = BufReader::new(file);

    let mut features = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            continue; // malformed / not a feature line (e.g. a stray FASTA block)
        }
        if cols[0] == "##FASTA" {
            break;
        }

        let seqid = cols[0].to_string();
        let feature_type = cols[2].to_string();
        let start_1based: u64 = cols[3]
            .parse()
            .with_context(|| format!("GFF3 line {}: bad start", lineno + 1))?;
        let end_1based: u64 = cols[4]
            .parse()
            .with_context(|| format!("GFF3 line {}: bad end", lineno + 1))?;
        let strand = match cols[6] {
            "+" => Strand::Forward,
            "-" => Strand::Reverse,
            _ => Strand::Unknown,
        };
        let attributes = cols.get(8).unwrap_or(&"").to_string();

        features.push(Feature {
            seqid,
            feature_type,
            start: start_1based.saturating_sub(1),
            end: end_1based,
            strand,
            attributes,
        });
    }

    Ok(features)
}
