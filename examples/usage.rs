use std::io::{stdout, BufWriter, Write};

use rs_moddotplot::{compute_self_identity, Row};

fn main() {
    let bed = compute_self_identity("data/chm13_chr1.fa", None, 4);
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", Row::header()).unwrap();
    for row in bed {
        writeln!(&mut fh, "{}", row.tsv()).unwrap();
    }
}
