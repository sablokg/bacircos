//! circos-rs: Circos-style circular genome plots for bacterial genomes and
//! metagenome assemblies, from FASTA + GFF3 + (optionally) BAM, and
//! multi-genome synteny/alignment plots from FASTA + PAF.

mod coverage;
mod genome;
mod gff;
mod paf;
mod plot;

/*
Gaurav Sablok
gsablok@proton.me
*/

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};

use genome::Genome;
use plot::{Canvas, GenomeInput, Layout, Ring};

const PALETTE: &[&str] = &[
    "#4C78A8", "#F58518", "#54A24B", "#B279A2", "#E45756", "#72B7B2", "#EECA3B", "#9D755D",
];

#[derive(Clone, Copy, Debug, ValueEnum)]
enum LinkColorBy {
    /// Color each ribbon by its query genome (matches the ideogram colors).
    Genome,
    /// Color each ribbon on a gray->blue->red scale by alignment identity.
    Identity,
}

/// Generate a Circos-style circular plot.
///
/// Single-genome mode (one `--fasta`): contig ideogram + GC-content ring +
/// gene-strand ring (from `--gff`) + read-depth ring (from `--bam` or
/// `--depth-tsv`).
///
/// Multi-genome mode (`--fasta` given more than once): the same, plus a
/// genome-level band and synteny/alignment ribbons between genomes, read
/// from one or more minimap2 PAF files via `--paf`
/// (e.g. `minimap2 -x asm5 genomeA.fasta genomeB.fasta > aln.paf`).
#[derive(Parser, Debug)]
#[command(name = "circos-rs", version, about)]
struct Args {
    /// Genome/metagenome assembly FASTA. Repeat for multiple genomes
    /// (`--fasta a.fasta --fasta b.fasta --fasta c.fasta`) to enable
    /// multi-genome / synteny mode.
    #[arg(long, required = true)]
    fasta: Vec<PathBuf>,

    /// Display label for each `--fasta`, in the same order. Defaults to
    /// the file stem if omitted or shorter than `--fasta`.
    #[arg(long)]
    label: Vec<String>,

    /// GFF3 annotation file(s). One draws a gene ring using that file's
    /// features on whichever plotted contigs they match; repeat to supply
    /// annotations for more than one genome.
    #[arg(long)]
    gff: Vec<PathBuf>,

    /// BAM file(s) for a read-depth ring (pure-Rust reader, no samtools
    /// required). Repeatable; depth from all given BAMs is merged into one
    /// ring. Mutually exclusive with --depth-tsv.
    #[arg(long)]
    bam: Vec<PathBuf>,

    /// `samtools depth`-style TSV file(s) (contig\tpos\tdepth), as an
    /// alternative to --bam. Mutually exclusive with --bam.
    #[arg(long)]
    depth_tsv: Vec<PathBuf>,

    /// minimap2 PAF alignment file(s) for synteny/link ribbons between
    /// genomes. Only meaningful with more than one --fasta. Generate with
    /// e.g. `minimap2 -x asm5 genomeA.fasta genomeB.fasta > aln.paf`.
    #[arg(long)]
    paf: Vec<PathBuf>,

    /// Drop alignments shorter than this many bp (filters spurious short
    /// matches out of noisy whole-genome PAFs).
    #[arg(long, default_value_t = 1000)]
    min_align_len: u64,

    /// Drop alignments below this identity fraction (0.0-1.0).
    #[arg(long, default_value_t = 0.0)]
    min_identity: f64,

    /// Only draw the N longest alignments (by alignment block length), for
    /// readability/performance on very dense whole-genome PAFs.
    #[arg(long)]
    max_links: Option<usize>,

    /// How to color synteny ribbons.
    #[arg(long, value_enum, default_value_t = LinkColorBy::Genome)]
    color_links_by: LinkColorBy,

    /// Output SVG path.
    #[arg(long, default_value = "circos.svg")]
    out: PathBuf,

    /// Plot title.
    #[arg(long, default_value = "")]
    title: String,

    /// Canvas size in SVG user units (square).
    #[arg(long, default_value_t = 900.0)]
    size: f64,

    /// Window size (bp) for the GC-content and coverage tracks.
    #[arg(long, default_value_t = 1000)]
    window: u64,

    /// Angular gap between contigs within one genome, in degrees.
    #[arg(long, default_value_t = 1.5)]
    gap_degrees: f64,

    /// Angular gap between genomes, in degrees (multi-genome mode only).
    #[arg(long, default_value_t = 4.0)]
    genome_gap_degrees: f64,

    /// Only draw the N longest contigs *per genome*.
    #[arg(long)]
    top_contigs: Option<usize>,

    /// Minimum contig length (bp) to include.
    #[arg(long, default_value_t = 0)]
    min_contig_len: u64,
}

struct LoadedGenome {
    label: String,
    genome: Genome,
    contigs: Vec<(String, u64)>, // filtered/sorted, as actually plotted
}

fn main() -> Result<()> {
    let args = Args::parse();

    if !args.bam.is_empty() && !args.depth_tsv.is_empty() {
        anyhow::bail!("pass --bam or --depth-tsv, not both");
    }
    if !args.paf.is_empty() && args.fasta.len() < 2 {
        anyhow::bail!("--paf requires at least two --fasta genomes to link between");
    }

    // --- Load every genome -------------------------------------------------
    let mut loaded = Vec::new();
    for (i, fasta_path) in args.fasta.iter().enumerate() {
        eprintln!("[circos-rs] reading FASTA {:?}", fasta_path);
        let genome = Genome::from_fasta(fasta_path)?;
        let label = args
            .label
            .get(i)
            .cloned()
            .unwrap_or_else(|| default_label(fasta_path));

        let mut contigs: Vec<(String, u64)> = genome
            .contigs
            .iter()
            .map(|c| (c.id.clone(), c.len()))
            .filter(|(_, len)| *len >= args.min_contig_len)
            .collect();
        contigs.sort_by(|a, b| b.1.cmp(&a.1));
        if let Some(n) = args.top_contigs {
            contigs.truncate(n);
        }
        if contigs.is_empty() {
            anyhow::bail!(
                "genome {:?} ({}): no contigs left after filtering",
                fasta_path,
                label
            );
        }
        eprintln!(
            "[circos-rs]   {}: {} contig(s), {} bp",
            label,
            contigs.len(),
            contigs.iter().map(|(_, l)| l).sum::<u64>()
        );
        loaded.push(LoadedGenome {
            label,
            genome,
            contigs,
        });
    }

    // Warn (don't fail) on contig id collisions across genomes -- the
    // layout/link lookups below assume ids are unique across all plotted
    // genomes, which holds for essentially all real assemblies (accessions
    // or assembler-generated contig names don't collide across samples).
    {
        let mut seen: HashMap<&str, &str> = HashMap::new();
        for g in &loaded {
            for (id, _) in &g.contigs {
                if let Some(prev) = seen.insert(id.as_str(), g.label.as_str()) {
                    eprintln!(
                        "[circos-rs] WARNING: contig id {:?} appears in both {:?} and {:?}; links/tracks for it may be ambiguous",
                        id, prev, g.label
                    );
                }
            }
        }
    }

    let is_multi = loaded.len() > 1;

    // --- Build the shared angle layout -------------------------------------
    let genome_inputs: Vec<GenomeInput> = loaded
        .iter()
        .map(|g| GenomeInput {
            label: if is_multi {
                g.label.clone()
            } else {
                String::new()
            },
            contigs: g.contigs.clone(),
        })
        .collect();
    let layout = Layout::new_multi(
        &genome_inputs,
        args.gap_degrees,
        if is_multi {
            args.genome_gap_degrees
        } else {
            0.0
        },
    );

    let plotted_ids: std::collections::HashSet<&str> = loaded
        .iter()
        .flat_map(|g| g.contigs.iter().map(|(id, _)| id.as_str()))
        .collect();
    let mut contig_lengths: HashMap<String, u64> = HashMap::new();
    for g in &loaded {
        contig_lengths.extend(g.genome.lengths());
    }

    // --- Gene annotations ----------------------------------------------------
    let mut features = Vec::new();
    for gff_path in &args.gff {
        eprintln!("[circos-rs] reading GFF3 {:?}", gff_path);
        let all = gff::read_gff3(gff_path)?;
        let kept = all
            .into_iter()
            .filter(|f| f.feature_type == "gene" || f.feature_type == "CDS")
            .filter(|f| plotted_ids.contains(f.seqid.as_str()));
        features.extend(kept);
    }
    if !features.is_empty() {
        eprintln!(
            "[circos-rs] {} gene/CDS features on plotted contigs",
            features.len()
        );
    }

    // --- Depth/coverage (BAM and/or TSV, merged) ------------------------------
    let mut depth_track: coverage::DepthTrack = HashMap::new();
    for bam_path in &args.bam {
        eprintln!(
            "[circos-rs] computing depth from BAM {:?} (pure-Rust noodles reader)",
            bam_path
        );
        let d = coverage::depth_from_bam(bam_path, args.window, &contig_lengths)
            .context("computing depth from BAM")?;
        depth_track.extend(d);
    }
    for tsv_path in &args.depth_tsv {
        eprintln!("[circos-rs] computing depth from TSV {:?}", tsv_path);
        let d = coverage::depth_from_tsv(tsv_path, args.window, &contig_lengths)
            .context("computing depth from depth TSV")?;
        depth_track.extend(d);
    }
    let depth_track = if depth_track.is_empty() {
        None
    } else {
        Some(depth_track)
    };

    // --- Alignments for synteny ribbons ---------------------------------------
    let mut links = Vec::new();
    for paf_path in &args.paf {
        eprintln!("[circos-rs] reading PAF {:?}", paf_path);
        let recs = paf::read_paf(paf_path)?;
        let before = recs.len();
        let kept: Vec<_> = recs
            .into_iter()
            .filter(|r| r.alnlen >= args.min_align_len)
            .filter(|r| r.identity() >= args.min_identity)
            .filter(|r| {
                plotted_ids.contains(r.qname.as_str()) && plotted_ids.contains(r.tname.as_str())
            })
            .collect();
        eprintln!(
            "[circos-rs]   {} / {} alignments kept after filtering",
            kept.len(),
            before
        );
        links.extend(kept);
    }
    if let Some(max) = args.max_links {
        links.sort_by(|a, b| b.alnlen.cmp(&a.alnlen));
        links.truncate(max);
        eprintln!(
            "[circos-rs] capped to top {} link(s) by alignment length (--max-links)",
            links.len()
        );
    }

    render(
        &args,
        &loaded,
        &layout,
        &features,
        depth_track.as_ref(),
        &links,
    )?;

    eprintln!("[circos-rs] wrote {:?}", args.out);
    Ok(())
}

fn render(
    args: &Args,
    loaded: &[LoadedGenome],
    layout: &Layout,
    features: &[gff::Feature],
    depth_track: Option<&coverage::DepthTrack>,
    links: &[paf::PafRecord],
) -> Result<()> {
    let mut canvas = Canvas::new(args.size, &args.title);
    let is_multi = loaded.len() > 1;

    let margin = args.size * 0.1;
    let mut r = args.size / 2.0 - margin;

    // Genome-level band (multi-genome mode only -- no-op / zero-height
    // otherwise since genome labels are empty strings).
    let genome_band = if is_multi {
        let ring = Ring {
            outer_r: r,
            inner_r: r - 10.0,
        };
        r = ring.inner_r - 4.0;
        Some(ring)
    } else {
        None
    };

    let ideogram = Ring {
        outer_r: r,
        inner_r: r - 18.0,
    };
    r = ideogram.inner_r - 6.0;

    let gc_ring = Ring {
        outer_r: r,
        inner_r: r - 40.0,
    };
    r = gc_ring.inner_r - 6.0;

    let gene_ring = if features.is_empty() {
        None
    } else {
        let ring = Ring {
            outer_r: r,
            inner_r: r - 22.0,
        };
        r = ring.inner_r - 6.0;
        Some(ring)
    };

    let depth_ring = depth_track.map(|_| {
        let ring = Ring {
            outer_r: r,
            inner_r: r - 45.0,
        };
        r = ring.inner_r - 6.0;
        ring
    });

    // Links are drawn innermost, anchored at whatever radius is left.
    let link_r = r;

    // 1. Genome band
    if let Some(ring) = &genome_band {
        canvas.genome_band(ring, layout, PALETTE);
    }

    // 2. Ideogram (contig arcs) -- colored per-genome in multi-genome mode
    // (all of one genome's contigs share a color) so genomes read as
    // distinct blocks; colored per-contig in single-genome mode as before.
    for span in &layout.spans {
        let color = PALETTE[span.genome_idx % PALETTE.len()];
        let label = if is_multi {
            short_id(&span.id).to_string()
        } else {
            format!("{} ({} kb)", short_id(&span.id), span.len / 1000)
        };
        canvas.ideogram_arc(&ideogram, span, color, &label);
    }

    // 3. GC-content ring
    canvas.ring_baseline(&gc_ring, "#cccccc");
    for g in loaded {
        for (id, _) in &g.contigs {
            let windows = g.genome.gc_windows(id, args.window);
            let (min_gc, max_gc) = windows
                .iter()
                .fold((1.0f64, 0.0f64), |(mn, mx), (_, _, gc)| {
                    (mn.min(*gc), mx.max(*gc))
                });
            let span = (max_gc - min_gc).max(0.01);
            let scaled: Vec<(u64, u64, f64)> = windows
                .iter()
                .map(|(a, b, gc)| (*a, *b, (gc - min_gc) / span))
                .collect();
            canvas.value_track(&gc_ring, layout, id, &scaled, "#59a14f");
        }
    }

    // 4. Gene track
    if let Some(ring) = &gene_ring {
        canvas.gene_track(ring, layout, features, "#1f77b4", "#d62728", 0.15);
    }

    // 5. Depth/coverage ring
    if let (Some(ring), Some(track)) = (&depth_ring, depth_track) {
        canvas.ring_baseline(ring, "#cccccc");
        let global_max = track
            .values()
            .flat_map(|v| v.iter())
            .cloned()
            .fold(0.0f64, f64::max)
            .max(1.0);
        for g in loaded {
            for (id, _) in &g.contigs {
                if let Some(bins) = track.get(id) {
                    let windows: Vec<(u64, u64, f64)> = bins
                        .iter()
                        .enumerate()
                        .map(|(i, v)| {
                            let start = i as u64 * args.window;
                            (start, start + args.window, v / global_max)
                        })
                        .collect();
                    canvas.value_track(ring, layout, id, &windows, "#e45756");
                }
            }
        }
    }

    // 6. Synteny/alignment ribbons
    for rec in links {
        let (Some(a0), Some(a1)) = (
            layout.angle_for(&rec.qname, rec.qstart),
            layout.angle_for(&rec.qname, rec.qend),
        ) else {
            continue;
        };
        let (Some(b0), Some(b1)) = (
            layout.angle_for(&rec.tname, rec.tstart),
            layout.angle_for(&rec.tname, rec.tend),
        ) else {
            continue;
        };
        let identity = rec.identity();
        let (color, opacity) = match args.color_links_by {
            LinkColorBy::Identity => (plot::identity_color(identity), 0.45),
            LinkColorBy::Genome => {
                let gi = layout.genome_idx_of(&rec.qname).unwrap_or(0);
                (
                    PALETTE[gi % PALETTE.len()].to_string(),
                    0.12 + 0.55 * identity,
                )
            }
        };
        canvas.link_ribbon(link_r, a0, a1, b0, b1, &color, opacity);
    }

    // Legend
    let mut legend_entries: Vec<(&str, &str)> = vec![("#59a14f", "GC content")];
    if gene_ring.is_some() {
        legend_entries.push(("#1f77b4", "Gene (+ strand)"));
        legend_entries.push(("#d62728", "Gene (- strand)"));
    }
    if depth_ring.is_some() {
        legend_entries.push(("#e45756", "Read depth"));
    }
    if !links.is_empty() {
        match args.color_links_by {
            LinkColorBy::Genome => {
                legend_entries.push(("#4C78A8", "Synteny link (colored by query genome)"))
            }
            LinkColorBy::Identity => {
                legend_entries.push(("#c8c8c8", "Link: low identity"));
                legend_entries.push(("#d62728", "Link: high identity"));
            }
        }
    }
    canvas.legend(
        14.0,
        args.size - 14.0 - legend_entries.len() as f64 * 16.0,
        &legend_entries,
    );

    let svg = canvas.finish();
    std::fs::write(&args.out, svg).with_context(|| format!("writing {:?}", args.out))?;
    Ok(())
}

fn default_label(path: &std::path::Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "genome".to_string())
}

fn short_id(id: &str) -> &str {
    if id.len() > 24 {
        &id[..24]
    } else {
        id
    }
}
