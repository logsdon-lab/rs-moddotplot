use rs_moddotplot::{compute_self_identity, plot_self_ident_tri};

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

    plot_self_ident_tri(&rows, &fname, &outfile, false);
}
