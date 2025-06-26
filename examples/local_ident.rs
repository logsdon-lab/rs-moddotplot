use std::io::{stdout, BufWriter, Write};

use rs_moddotplot::{compute_local_seq_self_identity, compute_self_identity, LocalRow};

fn main() {
    let mut args = std::env::args();
    let infile = args.nth(1).expect("Missing input fasta.");

    let bed = compute_self_identity(infile, None);
    let local_bed = compute_local_seq_self_identity(&bed, None);
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", LocalRow::header()).unwrap();
    for row in local_bed {
        writeln!(&mut fh, "{}", row.tsv()).unwrap();
    }
}
