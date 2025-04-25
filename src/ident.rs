use core::str;
use std::path::Path;
use std::usize;

use crate::ani::{convert_matrix_to_bed, create_self_matrix};
use crate::cfg::LocalSelfIdentConfig;
use crate::common::AIndexMap;

use ahash::AHashSet;
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
pub fn compute_seq_self_identity(
    seq: &str,
    name: &str,
    config: Option<SelfIdentConfig>,
) -> Vec<Row> {
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
    config: Option<LocalSelfIdentConfig>,
) -> Vec<LocalRow> {
    let cfg = config.unwrap_or_default();
    let window = cfg.window_size;
    let n_bins = cfg.n_bins;
    let ignore_bins = cfg.ignore_bins;
    let Some(chrom) = rows.first().map(|row| &row.reference_name) else {
        return vec![];
    };

    let mut aln_mtx: AIndexMap<usize, AIndexMap<usize, f32>> = AIndexMap::default();
    for line in rows {
        let x = line.query_start / window;
        let y = line.reference_start / window;
        let ident = line.perc_id_by_events;
        // Convert position to indices.
        aln_mtx
            .entry(x)
            .and_modify(|rec| {
                rec.insert(y, ident);
            })
            .or_insert_with(|| AIndexMap::from_iter([(y, ident)]));
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
        let n_pos = idents.len() as f32;
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


/// cargo test --package rs-moddotplot --lib -- ident::test::test_self_ident --exact --show-output
pub fn compute_nonzero_seq_self_identity(
    rows: &[Row],
) -> Vec<LocalRow> {
    let window = 5000;
    let ignore_band = 20;
    let Some(chrom) = rows.first().map(|row| &row.reference_name) else {
        return vec![];
    };

    let mut binned_ident = vec![];
    let mut aln_mtx: AIndexMap<usize, AIndexMap<usize, f32>> = AIndexMap::default();
    for line in rows {
        let x = line.query_start / window;
        let y = line.reference_start / window;
        let ident = line.perc_id_by_events;
        // Convert position to indices.
        aln_mtx
            .entry(x)
            .and_modify(|rec| {
                rec.insert(y, ident);
            })
            .or_insert_with(|| AIndexMap::from_iter([(y, ident)]));
    }

    // DFS search.
    let mut traveled = AHashSet::new();
    for x in aln_mtx.keys() {
        let y = *x + ignore_band;
        // 7       * * *   +
        // 6       * *   + 
        // 5       *   +    
        // 4 * * *   +
        // 3 * *   +
        // 2 *   +
        // 1   +
        // 0 +
        //   0 1 2 3 4 5 6 7
        if traveled.contains(&(*x,y)) {
            continue;
        }
        let mut positions: Vec<(usize, usize) >= vec![(*x, y)];
        let mut idents: Vec<f32> = vec![];
        let mut max_x = *x;
        while let Some(position) = positions.pop() {
            let (x, y) = position;

            if traveled.contains(&(x,y)) {
                continue;
            }
            // Store position traveled.
            traveled.insert((x, y));

            // Stop if zero.
            println!("{:?}", (x,y));
            let Some(col) = aln_mtx.get(&x) else {
                // Update x since we've gone into region
                max_x = x - 1;
                continue;
            };
            let Some(ident) = col.get(&y) else {
                continue;
            };
            // Add next positions to queue.
            positions.push((x, y + 1));
            positions.push((x + 1, y));
            idents.push(*ident);
        }
        if idents.is_empty() {
            continue;
        }

        let start = x * window + 1;
        let end = max_x * window + 1;
        let n_pos = idents.len() as f32;
        // Calculate average identity within spanned region and min coordinates.
        binned_ident.push(LocalRow {
            chrom: chrom.to_owned(),
            start,
            end,
            avg_perc_id_by_events: idents.into_iter().sum::<f32>() / n_pos,
        });
    }
    binned_ident
}


#[cfg(test)]
mod test {
    use std::{fs::File, io::{BufRead, BufReader, BufWriter, Write}, path::Path};

    use crate::{compute_self_identity, Row};

    use super::compute_nonzero_seq_self_identity;
    
    #[test]
    fn test_self_ident() {
        let path_outfile = Path::new("rows.tsv");
        let rows = if let Ok(mut new_file) = File::create_new(path_outfile).map(BufWriter::new) {
            let rows = compute_self_identity("data/HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509.fa", None, 1);
            for r in rows.iter() {
                writeln!(new_file, "{}", r.tsv()).unwrap();
            }
            rows
        } else {
            let reader = BufReader::new(File::open(path_outfile).unwrap());
            let mut rows = vec![];
            for line in reader.lines() {
                let line = line.unwrap();
                let [qname, qst, qend, rname, rst, rend, ident] = line.trim().split('\t').collect::<Vec<&str>>()[..] else {
                    panic!()
                };
                rows.push(Row {
                    query_name: qname.to_owned(),
                    query_start: qst.parse().unwrap(),
                    query_end: qend.parse().unwrap(),
                    reference_name: rname.to_owned(),
                    reference_start: rst.parse().unwrap(),
                    reference_end: rend.parse().unwrap(),
                    perc_id_by_events: ident.parse().unwrap(),
                });
            }
            rows
        };
        for r in compute_nonzero_seq_self_identity(&rows) {
            println!("{}", r.tsv())
        };       
    }
}
