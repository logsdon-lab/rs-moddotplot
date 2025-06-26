use std::io::{stdout, BufWriter, Write};

use rs_moddotplot::{compute_group_seq_self_identity, compute_self_identity, LocalRow};

fn main() {
    let mut args = std::env::args();
    let infile = args.nth(1).expect("Missing input fasta.");

    let bed = compute_self_identity(infile, None);
    let group_bed = compute_group_seq_self_identity(&bed);
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", LocalRow::header()).unwrap();
    for row in group_bed {
        writeln!(&mut fh, "{}", row.tsv()).unwrap();
    }
}
