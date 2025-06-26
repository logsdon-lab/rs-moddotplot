use std::io::{stdout, BufWriter, Write};

use rs_moddotplot::{compute_self_identity, Row};

fn main() {
    let mut args = std::env::args();
    let infile = args.nth(1).expect("Missing input fasta.");

    let bed = compute_self_identity(infile, None);
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", Row::header()).unwrap();
    for row in bed {
        writeln!(&mut fh, "{}", row.tsv()).unwrap();
    }
}
