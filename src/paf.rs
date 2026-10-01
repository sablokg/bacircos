//! Minimal PAF (Pairwise mApping Format) reader.
//!
//! PAF is minimap2's plain-text alignment output -- one line per alignment,
//! tab-separated, 12 mandatory columns followed by optional `tag:type:value`
//! fields we don't need here. It's the natural input for multi-genome
//! synteny links: `minimap2 -x asm5 genomeA.fasta genomeB.fasta > aln.paf`.
//!
//! Spec: <https://github.com/lh3/miniasm/blob/master/PAF.md>

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone)]
pub struct PafRecord {
    pub qname: String,
    pub qstart: u64,
    pub qend: u64,
    pub strand_reverse: bool,
    pub tname: String,
    pub tstart: u64,
    pub tend: u64,
    /// Number of matching bases in the alignment.
    pub nmatch: u64,
    /// Total alignment block length (matches + mismatches + gaps).
    pub alnlen: u64,
    pub mapq: u8,
}

impl PafRecord {
    /// Fraction of the alignment block that's an exact match, in [0, 1].
    pub fn identity(&self) -> f64 {
        if self.alnlen == 0 {
            0.0
        } else {
            (self.nmatch as f64 / self.alnlen as f64).clamp(0.0, 1.0)
        }
    }
}

pub fn read_paf(path: &Path) -> Result<Vec<PafRecord>> {
    let file = File::open(path).with_context(|| format!("opening PAF {:?}", path))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 12 {
            anyhow::bail!(
                "{:?} line {}: expected >=12 tab-separated columns, got {}",
                path,
                lineno + 1,
                cols.len()
            );
        }
        records.push(PafRecord {
            qname: cols[0].to_string(),
            qstart: cols[2]
                .parse()
                .with_context(|| format!("line {}: qstart", lineno + 1))?,
            qend: cols[3]
                .parse()
                .with_context(|| format!("line {}: qend", lineno + 1))?,
            strand_reverse: cols[4] == "-",
            tname: cols[5].to_string(),
            tstart: cols[7]
                .parse()
                .with_context(|| format!("line {}: tstart", lineno + 1))?,
            tend: cols[8]
                .parse()
                .with_context(|| format!("line {}: tend", lineno + 1))?,
            nmatch: cols[9]
                .parse()
                .with_context(|| format!("line {}: nmatch", lineno + 1))?,
            alnlen: cols[10]
                .parse()
                .with_context(|| format!("line {}: alnlen", lineno + 1))?,
            mapq: cols[11].parse().unwrap_or(0),
        });
    }

    Ok(records)
}
