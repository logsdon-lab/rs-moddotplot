/// Rust implementation of ModDotPlot
/// https://github.com/marbl/ModDotPlot/commit/50ecda4eff91acd00584090afd380d4a355be7aa
use std::path::Path;

mod ani;
mod cfg;
mod io;

use ani::{convert_matrix_to_bed, create_self_matrix};
pub use cfg::SelfIdentConfig;
pub use io::Row;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

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
    let kmers = io::read_kmers(fasta.as_ref(), k).unwrap();

    kmers
        .into_par_iter()
        .flat_map(|(seq, kmers)| {
            let mtx =
                create_self_matrix(kmers, window_size, delta, k, id_threshold, false, modimizer)
                    .unwrap();
            convert_matrix_to_bed(mtx, window_size, id_threshold, &seq, &seq, true)
        })
        .collect()
}

#[cfg(test)]
mod test {
    use crate::compute_self_identity;

    #[test]
    fn test_self_ident() {
        let bed = compute_self_identity("test/chm13_chr1.fa", None, 4);
    }
}
