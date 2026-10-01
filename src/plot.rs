//! Circular ("Circos-style") SVG renderer.
//!
//! Layout model: contigs are laid end-to-end around the circle in a fixed
//! order, separated by a small angular gap (so a metagenome assembly with
//! many contigs doesn't collapse into a solid ring). Everything downstream
//! (ideogram, GC track, gene arcs, coverage histogram) maps a genomic
//! coordinate `(contig, pos)` to an angle via one shared [`Layout`].

use std::f64::consts::PI;
use std::fmt::Write as _;

use crate::gff::{Feature, Strand};

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Clone)]
pub struct ContigSpan {
    pub id: String,
    pub len: u64,
    pub start_angle: f64,
    pub end_angle: f64,
    /// Index into `Layout::genomes` -- which input genome this contig
    /// belongs to. `0` for the single-genome case.
    pub genome_idx: usize,
}

/// One input genome's angular sector: a contiguous range of the circle
/// containing all of that genome's contigs, used for genome-level coloring,
/// labeling, and (when there's more than one genome) spacing genomes apart
/// more than contigs within a genome are spaced apart.
pub struct GenomeGroup {
    pub label: String,
    pub start_angle: f64,
    pub end_angle: f64,
}

/// One genome's worth of input to [`Layout::new_multi`]: a display label
/// plus its (contig_id, length) pairs, in the order they should be drawn.
pub struct GenomeInput {
    pub label: String,
    pub contigs: Vec<(String, u64)>,
}

/// Maps genomic coordinates to angles (radians, 0 = top / 12 o'clock,
/// increasing clockwise) for one or more genomes' contigs arranged around a
/// full circle. With multiple genomes, each genome gets its own contiguous
/// arc (sized proportionally to that genome's total length relative to the
/// sum of all genomes), separated by `genome_gap_degrees`; contigs within a
/// genome are separated by the smaller `gap_degrees`.
pub struct Layout {
    pub spans: Vec<ContigSpan>,
    pub genomes: Vec<GenomeGroup>,
}

impl Layout {
    /// Single-genome convenience constructor (equivalent to `new_multi`
    /// with one genome and no extra genome-level gap). Kept as public API
    /// for library users even though the CLI always goes through
    /// `new_multi` now.
    #[allow(dead_code)]
    pub fn new(contigs: &[(String, u64)], gap_degrees: f64) -> Self {
        Self::new_multi(
            &[GenomeInput {
                label: String::new(),
                contigs: contigs.to_vec(),
            }],
            gap_degrees,
            0.0,
        )
    }

    pub fn new_multi(genomes: &[GenomeInput], gap_degrees: f64, genome_gap_degrees: f64) -> Self {
        let gap_radians = gap_degrees.to_radians();
        let genome_gap_radians = genome_gap_degrees.to_radians();
        let n_genomes = genomes.len().max(1);

        let genome_totals: Vec<u64> = genomes
            .iter()
            .map(|g| g.contigs.iter().map(|(_, l)| *l).sum::<u64>().max(1))
            .collect();
        let grand_total: u64 = genome_totals.iter().sum::<u64>().max(1);

        let usable_radians = 2.0 * PI - genome_gap_radians * n_genomes as f64;

        let mut spans = Vec::new();
        let mut genome_groups = Vec::with_capacity(genomes.len());
        let mut angle = -PI / 2.0; // start at 12 o'clock

        for (genome_idx, genome) in genomes.iter().enumerate() {
            let genome_start = angle;
            let genome_radians =
                usable_radians * (genome_totals[genome_idx] as f64 / grand_total as f64);
            // Within this genome's arc, lay out its own contigs the same
            // way the single-genome layout does, but scaled to fit
            // `genome_radians` instead of the full circle.
            let n_contigs = genome.contigs.len().max(1);
            let contig_gap_total = gap_radians * n_contigs as f64;
            let usable_contig_radians = (genome_radians - contig_gap_total).max(0.0);
            let genome_total_len = genome_totals[genome_idx];

            let mut inner_angle = angle;
            for (id, len) in &genome.contigs {
                let span_radians = usable_contig_radians * (*len as f64 / genome_total_len as f64);
                spans.push(ContigSpan {
                    id: id.clone(),
                    len: *len,
                    start_angle: inner_angle,
                    end_angle: inner_angle + span_radians,
                    genome_idx,
                });
                inner_angle += span_radians + gap_radians;
            }

            angle = genome_start + genome_radians + genome_gap_radians;
            genome_groups.push(GenomeGroup {
                label: genome.label.clone(),
                start_angle: genome_start,
                end_angle: genome_start + genome_radians,
            });
        }

        Layout {
            spans,
            genomes: genome_groups,
        }
    }

    /// Angle in radians for position `pos` (0-based) within `contig_id`.
    pub fn angle_for(&self, contig_id: &str, pos: u64) -> Option<f64> {
        let span = self.spans.iter().find(|s| s.id == contig_id)?;
        let frac = if span.len == 0 {
            0.0
        } else {
            (pos.min(span.len) as f64) / (span.len as f64)
        };
        Some(span.start_angle + frac * (span.end_angle - span.start_angle))
    }

    pub fn genome_idx_of(&self, contig_id: &str) -> Option<usize> {
        self.spans
            .iter()
            .find(|s| s.id == contig_id)
            .map(|s| s.genome_idx)
    }
}

fn polar(cx: f64, cy: f64, r: f64, angle: f64) -> (f64, f64) {
    (cx + r * angle.cos(), cy + r * angle.sin())
}

/// One ring (track) of the plot. Rings are drawn outer-to-inner in the
/// order they're pushed.
pub struct Ring {
    pub outer_r: f64,
    pub inner_r: f64,
}

pub struct Canvas {
    pub size: f64,
    pub cx: f64,
    pub cy: f64,
    svg: String,
}

impl Canvas {
    pub fn new(size: f64, title: &str) -> Self {
        let mut svg = String::new();
        let _ = write!(
            svg,
            r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" font-family="Helvetica, Arial, sans-serif">
<rect width="{size}" height="{size}" fill="#ffffff"/>
"##,
            size = size
        );
        if !title.is_empty() {
            let _ = write!(
                svg,
                r##"<text x="{cx}" y="24" text-anchor="middle" font-size="18" font-weight="600" fill="#1a1a1a">{title}</text>
"##,
                cx = size / 2.0,
                title = escape_xml(title)
            );
        }
        Canvas {
            size,
            cx: size / 2.0,
            cy: size / 2.0,
            svg,
        }
    }

    pub fn finish(mut self) -> String {
        self.svg.push_str("</svg>\n");
        self.svg
    }

    /// Ideogram arc for one contig: a thick band from `start_angle` to
    /// `end_angle` at radius `[inner_r, outer_r]`.
    pub fn ideogram_arc(&mut self, ring: &Ring, span: &ContigSpan, color: &str, label: &str) {
        let path = arc_band_path(
            self.cx,
            self.cy,
            ring.inner_r,
            ring.outer_r,
            span.start_angle,
            span.end_angle,
        );
        let _ = writeln!(
            self.svg,
            r##"<path d="{path}" fill="{color}" stroke="#333333" stroke-width="0.5"/>"##,
        );

        // Label along the outside of the arc, oriented so it's never
        // upside-down (flip on the bottom half of the circle).
        let mid_angle = (span.start_angle + span.end_angle) / 2.0;
        let label_r = ring.outer_r + 14.0;
        let (lx, ly) = polar(self.cx, self.cy, label_r, mid_angle);
        let mut deg = mid_angle.to_degrees() + 90.0;
        let mut anchor = "middle";
        if !(-90.0..=90.0).contains(&((deg + 360.0) % 360.0 - 180.0)) {
            // keep purely for readability edge case; default anchor is fine
            anchor = "middle";
        }
        if mid_angle.sin() > 0.15 {
            deg += 180.0; // flip text on the lower half so it reads left-to-right
        }
        let _ = writeln!(
            self.svg,
            r##"<text x="{lx:.2}" y="{ly:.2}" transform="rotate({deg:.2} {lx:.2} {ly:.2})" text-anchor="{anchor}" font-size="9" fill="#333333">{label}</text>"##,
            label = escape_xml(label)
        );
    }

    /// Draw a filled histogram-style track from per-window values (already
    /// scaled 0.0-1.0) laid out via `layout`.
    pub fn value_track(
        &mut self,
        ring: &Ring,
        layout: &Layout,
        contig_id: &str,
        windows: &[(u64, u64, f64)],
        color: &str,
    ) {
        for (w_start, w_end, value) in windows {
            let Some(a0) = layout.angle_for(contig_id, *w_start) else {
                continue;
            };
            let Some(a1) = layout.angle_for(contig_id, *w_end) else {
                continue;
            };
            let r = ring.inner_r + (ring.outer_r - ring.inner_r) * value.clamp(0.0, 1.0);
            let path = arc_band_path(self.cx, self.cy, ring.inner_r, r, a0, a1);
            let _ = writeln!(
                self.svg,
                r##"<path d="{path}" fill="{color}" stroke="none"/>"##
            );
        }
    }

    /// Baseline circle for a ring (drawn first, under value tracks).
    pub fn ring_baseline(&mut self, ring: &Ring, color: &str) {
        let _ = writeln!(
            self.svg,
            r##"<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="{color}" stroke-width="0.75" stroke-dasharray="1,2"/>"##,
            cx = self.cx,
            cy = self.cy,
            r = ring.inner_r
        );
    }

    /// One tick/arc per gene feature, colored by strand.
    pub fn gene_track(
        &mut self,
        ring: &Ring,
        layout: &Layout,
        features: &[Feature],
        fwd_color: &str,
        rev_color: &str,
        min_arc_degrees: f64,
    ) {
        for f in features {
            let Some(mut a0) = layout.angle_for(&f.seqid, f.start) else {
                continue;
            };
            let Some(mut a1) = layout.angle_for(&f.seqid, f.end) else {
                continue;
            };
            if a1 < a0 {
                std::mem::swap(&mut a0, &mut a1);
            }
            // Ensure very short genes are still visible.
            let min_rad = min_arc_degrees.to_radians();
            if a1 - a0 < min_rad {
                let mid = (a0 + a1) / 2.0;
                a0 = mid - min_rad / 2.0;
                a1 = mid + min_rad / 2.0;
            }
            let color = match f.strand {
                Strand::Forward => fwd_color,
                Strand::Reverse => rev_color,
                Strand::Unknown => "#888888",
            };
            let path = arc_band_path(self.cx, self.cy, ring.inner_r, ring.outer_r, a0, a1);
            let _ = writeln!(
                self.svg,
                r##"<path d="{path}" fill="{color}" stroke="none"/>"##
            );
        }
    }

    /// Genome-level band drawn one step outside the contig ideogram: one
    /// arc per input genome (spanning all of its contigs), with the
    /// genome's display label centered on it. This is what makes a
    /// multi-genome plot read as "genome A | genome B | genome C" at a
    /// glance, on top of the finer per-contig ideogram ring just inside it.
    pub fn genome_band(&mut self, ring: &Ring, layout: &Layout, colors: &[&str]) {
        for (i, g) in layout.genomes.iter().enumerate() {
            if g.label.is_empty() {
                continue; // single-genome mode: nothing to label
            }
            let color = colors[i % colors.len()];
            let path = arc_band_path(
                self.cx,
                self.cy,
                ring.inner_r,
                ring.outer_r,
                g.start_angle,
                g.end_angle,
            );
            let _ = writeln!(
                self.svg,
                r##"<path d="{path}" fill="{color}" fill-opacity="0.28" stroke="{color}" stroke-width="1"/>"##
            );

            let mid_angle = (g.start_angle + g.end_angle) / 2.0;
            let label_r = ring.outer_r + 16.0;
            let (lx, ly) = polar(self.cx, self.cy, label_r, mid_angle);
            let mut deg = mid_angle.to_degrees() + 90.0;
            if mid_angle.sin() > 0.15 {
                deg += 180.0;
            }
            let _ = writeln!(
                self.svg,
                r##"<text x="{lx:.2}" y="{ly:.2}" transform="rotate({deg:.2} {lx:.2} {ly:.2})" text-anchor="middle" font-size="12" font-weight="600" fill="{color}">{label}</text>"##,
                label = escape_xml(&g.label)
            );
        }
    }

    /// Draw a synteny/alignment ribbon connecting angular range `[a0, a1]`
    /// (the query/source block) to `[b0, b1]` (the target block), both at
    /// radius `r`. The two connecting edges are quadratic Beziers through
    /// the circle's center, which is the standard Circos "link" look: it
    /// reads as a bowtie/ribbon rather than two disconnected arcs.
    pub fn link_ribbon(
        &mut self,
        r: f64,
        a0: f64,
        a1: f64,
        b0: f64,
        b1: f64,
        color: &str,
        opacity: f64,
    ) {
        let large_a = if (a1 - a0).abs() > PI { 1 } else { 0 };
        let large_b = if (b1 - b0).abs() > PI { 1 } else { 0 };
        let (pa0x, pa0y) = polar(self.cx, self.cy, r, a0);
        let (pa1x, pa1y) = polar(self.cx, self.cy, r, a1);
        let (pb0x, pb0y) = polar(self.cx, self.cy, r, b0);
        let (pb1x, pb1y) = polar(self.cx, self.cy, r, b1);
        let path = format!(
            "M {pa0x:.3},{pa0y:.3} \
             A {r:.3},{r:.3} 0 {large_a} 1 {pa1x:.3},{pa1y:.3} \
             Q {cx:.3},{cy:.3} {pb0x:.3},{pb0y:.3} \
             A {r:.3},{r:.3} 0 {large_b} 1 {pb1x:.3},{pb1y:.3} \
             Q {cx:.3},{cy:.3} {pa0x:.3},{pa0y:.3} Z",
            cx = self.cx,
            cy = self.cy,
        );
        let _ = writeln!(
            self.svg,
            r##"<path d="{path}" fill="{color}" fill-opacity="{opacity:.3}" stroke="none"/>"##,
        );
    }

    pub fn legend(&mut self, x: f64, y: f64, entries: &[(&str, &str)]) {
        let mut yy = y;
        for (color, label) in entries {
            let _ = writeln!(
                self.svg,
                r##"<rect x="{x}" y="{yy}" width="10" height="10" fill="{color}"/><text x="{tx}" y="{ty}" font-size="11" fill="#333333">{label}</text>"##,
                x = x,
                yy = yy,
                tx = x + 15.0,
                ty = yy + 9.5,
                label = escape_xml(label)
            );
            yy += 16.0;
        }
    }
}

/// SVG path for an annular sector (a "band") between two radii and two
/// angles — the fundamental primitive for every ring in the plot.
fn arc_band_path(cx: f64, cy: f64, inner_r: f64, outer_r: f64, a0: f64, a1: f64) -> String {
    let large_arc = if (a1 - a0).abs() > PI { 1 } else { 0 };
    let (ox0, oy0) = polar(cx, cy, outer_r, a0);
    let (ox1, oy1) = polar(cx, cy, outer_r, a1);
    let (ix1, iy1) = polar(cx, cy, inner_r, a1);
    let (ix0, iy0) = polar(cx, cy, inner_r, a0);
    format!(
        "M {ox0:.3},{oy0:.3} A {outer_r:.3},{outer_r:.3} 0 {large_arc} 1 {ox1:.3},{oy1:.3} \
         L {ix1:.3},{iy1:.3} A {inner_r:.3},{inner_r:.3} 0 {large_arc} 0 {ix0:.3},{iy0:.3} Z"
    )
}

/// Map an identity fraction (0.0-1.0) to a hex color on a light-gray ->
/// blue -> red gradient, for coloring links by alignment quality when
/// `--color-links-by identity` is requested.
pub fn identity_color(t: f64) -> String {
    let t = t.clamp(0.0, 1.0);
    // Two-segment gradient: gray (0.0) -> blue (0.7) -> red (1.0), so most
    // real alignments (typically >90% identity) land in the blue-to-red
    // range where differences are easiest to see.
    let (r, g, b) = if t < 0.7 {
        let u = t / 0.7;
        lerp_rgb((200, 200, 200), (60, 100, 220), u)
    } else {
        let u = (t - 0.7) / 0.3;
        lerp_rgb((60, 100, 220), (214, 39, 40), u)
    };
    format!("#{r:02x}{g:02x}{b:02x}")
}

fn lerp_rgb(a: (u8, u8, u8), b: (u8, u8, u8), t: f64) -> (u8, u8, u8) {
    let l = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * t).round() as u8;
    (l(a.0, b.0), l(a.1, b.1), l(a.2, b.2))
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
