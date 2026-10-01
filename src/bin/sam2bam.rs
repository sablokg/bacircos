//! Tiny dev-only helper: `sam2bam in.sam out.bam`.
//!
//! Not part of the public `circos-rs` tool -- this exists purely so the
//! test/example fixtures in this repo can be produced without requiring
//! `samtools` on the build machine (everything here is pure-Rust `noodles`).

use std::env;
use std::fs::File;
use std::io::BufReader;

/*
Gaurav Sablok
gsablok@proton.me
 */

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = env::args().collect();
    let in_path = &args[1];
    let out_path = &args[2];

    let mut reader = sam::io::Reader::new(BufReader::new(File::open(in_path)?));
    let header = reader.read_header()?;

    let out = File::create(out_path)?;
    let mut writer = bam::io::Writer::new(out);
    writer.write_header(&header)?;

    for result in reader.record_bufs(&header) {
        let record = result?;
        writer.write_alignment_record(&header, &record)?;
    }

    Ok(())
}
