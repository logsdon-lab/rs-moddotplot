use std::io::{stdout, BufWriter, Write};

use rs_moddotplot::{
    compute_self_identity, find_simple_inversion_duplications, EventInterval, SelfIdentConfig,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let infile = args.get(1).expect("Missing input fasta.");
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", EventInterval::header()).unwrap();

    let rows = compute_self_identity(
        infile,
        Some(SelfIdentConfig {
            id_threshold: 99.0,
            ..Default::default()
        }),
    );
    let groups = find_simple_inversion_duplications(&rows);
    // for row in groups {
    //     println!("{}\t{}\t{}\t{:?}", row.chrom, row.start, row.end, row.event)
    // }
}
