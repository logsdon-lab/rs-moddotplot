use core::str;
use std::path::Path;

use crate::ani::{convert_matrix_to_bed, create_self_matrix};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::io::generate_kmers_from_fasta;
use crate::{io, Row, SelfIdentConfig};

/// Compute self-identity between sequences in a given fasta file.
///
/// # Args
/// * fasta
///     * Fasta input file.
/// * config
///     * Configuration for ANI. Similar to ModDotPlot.
/// * threads
///     * Number of threads.
pub fn compute_self_identity(
    fasta: impl AsRef<Path>,
    config: Option<SelfIdentConfig>,
    threads: usize,
) -> Vec<Row> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let cfg = config.unwrap_or_default();
    let window_size = cfg.window_size;
    let delta = cfg.delta;
    let k = cfg.k;
    let id_threshold = cfg.id_threshold;
    let modimizer = cfg.modimizer;
    let kmers = io::read_kmers(fasta.as_ref(), k);

    kmers
        .into_par_iter()
        .flat_map(|(seq, kmers)| {
            let mtx =
                create_self_matrix(kmers, window_size, delta, k, id_threshold, false, modimizer);
            convert_matrix_to_bed(mtx, window_size, id_threshold, &seq, &seq, true)
        })
        .collect()
}

/// Compute self-identity for a single sequence.
///
/// # Args
/// * seq
///     * Input sequence.
/// * name
///     * Input sequence name.
/// * config
///     * Configuration for ANI. Similar to ModDotPlot.
pub fn compute_seq_self_identity(seq: &str, name: &str, config: Option<SelfIdentConfig>) -> Vec<Row> {
    let cfg = config.unwrap_or_default();
    let window_size = cfg.window_size;
    let delta = cfg.delta;
    let k = cfg.k;
    let id_threshold = cfg.id_threshold;
    let modimizer = cfg.modimizer;

    let kmers = generate_kmers_from_fasta(str::from_utf8(seq.as_bytes()).unwrap(), k);
    let mtx = create_self_matrix(kmers, window_size, delta, k, id_threshold, false, modimizer);
    convert_matrix_to_bed(mtx, window_size, id_threshold, name, name, true)
}

// #[cfg(test)]
// mod test {
//     use crate::compute_self_identity;

//     #[test]
//     fn test_self_ident() {
//         let bed = compute_self_identity("test/chm13_chr1.fa", None, 4);
//     }
// }
