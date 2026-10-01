//! Minimal FASTA reader for genome/metagenome assemblies.
//!
//! Hand-rolled rather than pulled from a crate: FASTA is a trivial format
//! and this keeps the dependency tree (and MSRV) small.

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/*
Gaurav Sablok
gsablok@proton.me
*/

/// A single contig/chromosome/plasmid sequence from the assembly.
pub struct Contig {
    pub id: String,
    pub seq: Vec<u8>,
}

impl Contig {
    pub fn len(&self) -> u64 {
        self.seq.len() as u64
    }
}

/// The full set of contigs loaded from a FASTA file, in file order (which
/// we preserve as the default ideogram draw order).
pub struct Genome {
    pub contigs: Vec<Contig>,
}

impl Genome {
    pub fn from_fasta(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("opening FASTA {:?}", path))?;
        let reader = BufReader::new(file);

        let mut contigs = Vec::new();
        let mut current_id: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for line in reader.lines() {
            let line = line?;
            if let Some(header) = line.strip_prefix('>') {
                if let Some(id) = current_id.take() {
                    contigs.push(Contig {
                        id,
                        seq: std::mem::take(&mut current_seq),
                    });
                }
                // Contig id is the first whitespace-delimited token, matching
                // samtools/BAM @SQ SN: convention so FASTA and BAM/GFF line up.
                let id = header
                    .split_whitespace()
                    .next()
                    .unwrap_or(header)
                    .to_string();
                current_id = Some(id);
            } else {
                current_seq.extend(line.trim_end().bytes().filter(|b| !b.is_ascii_whitespace()));
            }
        }
        if let Some(id) = current_id.take() {
            contigs.push(Contig {
                id,
                seq: current_seq,
            });
        }

        if contigs.is_empty() {
            anyhow::bail!("no sequences found in FASTA {:?}", path);
        }

        Ok(Genome { contigs })
    }

    #[allow(dead_code)]
    pub fn total_len(&self) -> u64 {
        self.contigs.iter().map(Contig::len).sum()
    }

    pub fn lengths(&self) -> HashMap<String, u64> {
        self.contigs
            .iter()
            .map(|c| (c.id.clone(), c.len()))
            .collect()
    }

    /// GC fraction (0.0-1.0) in fixed-size, non-overlapping windows across a
    /// contig. The final partial window (if any) is included at its actual
    /// (shorter) size.
    pub fn gc_windows(&self, contig_id: &str, window: u64) -> Vec<(u64, u64, f64)> {
        let Some(contig) = self.contigs.iter().find(|c| c.id == contig_id) else {
            return Vec::new();
        };
        let window = window.max(1) as usize;
        let mut out = Vec::new();
        let mut start = 0usize;
        while start < contig.seq.len() {
            let end = (start + window).min(contig.seq.len());
            let slice = &contig.seq[start..end];
            let gc = slice
                .iter()
                .filter(|b| matches!(b.to_ascii_uppercase(), b'G' | b'C'))
                .count();
            let frac = if slice.is_empty() {
                0.0
            } else {
                gc as f64 / slice.len() as f64
            };
            out.push((start as u64, end as u64, frac));
            start = end;
        }
        out
    }
}
