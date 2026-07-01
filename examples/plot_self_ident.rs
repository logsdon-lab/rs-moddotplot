#[cfg(feature = "plot")]
use rs_moddotplot::plot_self_ident_tri;
use rs_moddotplot::{compute_group_seq_self_identity, compute_self_identity};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let infile = args.get(1).expect("Missing input fasta.");
    let outfile = args.get(2).expect("Missing output plot.");

    let path = std::path::Path::new(&infile);
    let fname = path
        .file_stem()
        .expect("No filename?")
        .to_os_string()
        .into_string()
        .unwrap();

    let rows = compute_self_identity(infile, None);
    let groups = compute_group_seq_self_identity(&rows);

    #[cfg(feature = "plot")]
    plot_self_ident_tri(&rows, Some(&groups), None, &fname, &outfile, false);
}
