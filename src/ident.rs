use core::str;
use std::collections::HashMap;
use std::path::Path;

use crate::ani::{convert_matrix_to_bed, create_self_matrix};
use crate::cfg::LocalSelfIdentConfig;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::io::{generate_kmers_from_fasta, LocalRow};
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
/// 
/// # Returns
/// * Self-identity BED file matrix as a list of rows.
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

    let kmers = generate_kmers_from_fasta(seq, k);
    let mtx = create_self_matrix(kmers, window_size, delta, k, id_threshold, false, modimizer);
    convert_matrix_to_bed(mtx, window_size, id_threshold, name, name, true)
}


/// Compute the local sequence identity from a set of sequence self-identity matrix [`Row`]s.
/// 
/// # Args
/// * rows
///     * Sequence self-identity matrix rows.
/// * config
///     * [`LocalSelfIdentConfig`] configuration.
/// 
/// # Returns
/// * Local self-identity BED file matrix as a list of rows.
pub fn compute_local_seq_self_identity(
    rows: &[Row],
    config: Option<LocalSelfIdentConfig>
) -> Vec<LocalRow> {
    let cfg = config.unwrap_or_default();
    let window = cfg.window_size;
    let n_bins = cfg.n_bins;
    let ignore_bins = cfg.ignore_bins;
    let Some(chrom) = rows.get(0).map(|row| &row.reference_name) else {
        return vec![];
    };

    let mut aln_mtx: HashMap<usize, HashMap<usize, f32>> = HashMap::new();
    for line in rows {
        let x = line.query_start / window;
        let y = line.reference_start / window;
        let ident = line.perc_id_by_events;
        // Convert position to indices.
        aln_mtx.entry(x)
            .and_modify(|rec| { rec.insert(y, ident); })
            .or_insert_with(|| HashMap::from_iter([(y, ident)]));
    }
    let mut binned_ident = vec![];
    for st_idx in aln_mtx.keys() {
        let start = st_idx * window + 1;
        let end = start + window - 1;
        let band_end_idx = st_idx + n_bins;

        // Within the alignment matrix with a n_bins of 5 and ignore_bands of 2:
        // - '*' is the calculated aln band
        // - '+' is self aln.
        // 4 * * *   +
        // 3 * *   +
        // 2 *   +
        // 1   +
        // 0 +
        //   0 1 2 3 4
        let mut idents = vec![];
        for x in *st_idx..band_end_idx {
            for y in x + ignore_bins..band_end_idx {
                let ident = aln_mtx.get(&x).and_then(|col| col.get(&y)).unwrap_or(&0.0);
                idents.push(ident);
            }
        }
        let n_pos =  idents.len() as f32;
        binned_ident.push(LocalRow {
            chrom: chrom.to_owned(),
            start,
            end,
            avg_perc_id_by_events: idents.into_iter().sum::<f32>() / n_pos,
        });
    }
    binned_ident.sort_by(|r1, r2| r1.start.cmp(&r2.end));
    binned_ident
}

// #[cfg(test)]
// mod test {
//     use crate::compute_self_identity;

//     #[test]
//     fn test_self_ident() {
//         let bed = compute_self_identity("test/chm13_chr1.fa", None, 4);
//     }
// }
