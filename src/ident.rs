use core::str;
use std::collections::VecDeque;
use std::path::Path;

use crate::ani::{convert_matrix_to_bed, create_self_matrix};
use crate::cfg::LocalSelfIdentConfig;
use crate::common::AIndexMap;

use ahash::AHashSet;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::io::{generate_kmers_from_fasta, read_kmers, LocalRow};
use crate::{Row, SelfIdentConfig};

/// Compute self-identity between sequences in a given fasta file.
///
/// # Args
/// * fasta
///     * Fasta input file.
/// * config
///     * Configuration for ANI. Similar to ModDotPlot.
///
/// # Returns
/// * Self-identity BED file matrix as a list of rows.
pub fn compute_self_identity(fasta: impl AsRef<Path>, config: Option<SelfIdentConfig>) -> Vec<Row> {
    let cfg = config.unwrap_or_default();
    let window_size = cfg.window_size;
    let delta = cfg.delta;
    let k = cfg.k;
    let id_threshold = cfg.id_threshold;
    let modimizer = cfg.modimizer;
    let seed = cfg.seed;
    let kmers = read_kmers(fasta.as_ref(), k, seed);

    kmers
        .into_par_iter()
        .flat_map(|(seq, kmers)| {
            let mtx = create_self_matrix(
                kmers,
                window_size,
                delta,
                k,
                id_threshold,
                false,
                modimizer,
                seed,
            );
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
    let seed = cfg.seed;

    let kmers = generate_kmers_from_fasta(seq, k, seed);
    let mtx = create_self_matrix(
        kmers,
        window_size,
        delta,
        k,
        id_threshold,
        false,
        modimizer,
        seed,
    );
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

/// Compute the grouped sequence identity from a set of sequence self-identity matrix [`Row`]s.
/// * Works best when identity threshold set high as avoids transitioning to other non-zero regions.
///
/// # Args
/// * rows
///     * Sequence self-identity matrix rows. Each interval must be same length.
///
/// # Returns
/// * Group self-identity BED file matrix as a list of rows.
/// * Regions larger than window size are returned.
pub fn compute_group_seq_self_identity(rows: &[Row]) -> Vec<LocalRow> {
    let Some((chrom, window)) = rows
        .first()
        .map(|row| (&row.reference_name, row.query_end - row.query_start))
    else {
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

    // BFS search.
    let mut traveled = AHashSet::new();
    for x in aln_mtx.keys() {
        let y = *x;
        // Travel along the self-identity band and perform a breadth first search for any non-zero identity position.
        // Track traveled points and adjacent diagonals as exit condition.
        // * Adjacent diagonals indicate a transition into a different sequence identity group.
        // 7       * * * * +
        // 6       * * * +
        // 5       * * +
        // 4 * * * * +
        // 3 * * * +
        // 2 * * +
        // 1 * +
        // 0 +
        //   0 1 2 3 4 5 6 7
        let mut positions: VecDeque<(usize, usize)> = VecDeque::from_iter([(*x, y)]);
        let mut idents: Vec<f32> = vec![];
        let mut max_x = *x;
        while let Some(position) = positions.pop_front() {
            let (x, y) = position;

            if traveled.contains(&(x, y)) {
                continue;
            }
            // Store position traveled.
            traveled.insert((x, y));

            // Stop if both diagonal is zero.
            //   *
            //  x
            // *
            let up_right = aln_mtx
                .get(&(x + 1))
                .and_then(|col| y.checked_sub(1).and_then(|y| col.get(&y)));
            let down_left = x
                .checked_sub(1)
                .and_then(|x| aln_mtx.get(&x))
                .and_then(|col| col.get(&(y + 1)));

            let Some(col) = aln_mtx.get(&x) else {
                // Update x since we've gone into region.
                max_x = x;
                continue;
            };
            let Some(ident) = col.get(&y) else {
                continue;
            };
            idents.push(*ident);

            if up_right.is_none() && down_left.is_none() {
                max_x = x;
                continue;
            }
            // Add next positions to queue.
            // *
            // x *
            positions.push_back((x, y + 1));
            positions.push_back((x + 1, y));
        }
        if idents.is_empty() {
            continue;
        }

        let start = x * window + 1;
        let end = max_x * window + 1;
        let length = end - start;
        let n_pos = idents.len() as f32;
        // Ignore self diagonal.
        if length <= window {
            continue;
        }
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
    use flate2::read::GzDecoder;
    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    use crate::{
        compute_group_seq_self_identity, compute_local_seq_self_identity, compute_self_identity,
        Row,
    };

    fn read_self_ident(fname: &str) -> Vec<Row> {
        let reader_self_ident = BufReader::new(GzDecoder::new(File::open(fname).unwrap()));
        reader_self_ident
            .lines()
            .map_while(Result::ok)
            .map(|line| {
                let [qchrom, qst, qend, rchrom, rst, rend, ident] =
                    line.split('\t').collect::<Vec<&str>>()[..]
                else {
                    panic!("Invalid formatted line: {line}")
                };
                Row {
                    query_name: qchrom.to_owned(),
                    query_start: qst.parse().unwrap(),
                    query_end: qend.parse().unwrap(),
                    reference_name: rchrom.to_owned(),
                    reference_start: rst.parse().unwrap(),
                    reference_end: rend.parse().unwrap(),
                    perc_id_by_events: ident.parse().unwrap(),
                }
            })
            .collect::<Vec<Row>>()
    }

    #[test]
    fn test_self_ident() {
        let rows = compute_self_identity(
            "data/HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509.fa",
            None,
        );
        let reader_self_ident = BufReader::new(GzDecoder::new(
            File::open("test/self_ident/expected.bed.gz").unwrap(),
        ));
        for (line, row) in reader_self_ident.lines().map_while(Result::ok).zip(rows) {
            assert_eq!(line, row.tsv());
        }
    }

    #[test]
    fn test_group_ident() {
        let rows = read_self_ident("test/self_ident/expected.bed.gz");
        let grouped_rows = compute_group_seq_self_identity(&rows);
        let reader = BufReader::new(GzDecoder::new(
            File::open("test/group_ident/expected.bed.gz").unwrap(),
        ));
        for (line, row) in reader.lines().map_while(Result::ok).zip(grouped_rows) {
            assert_eq!(line, row.tsv());
        }
    }

    #[test]
    fn test_local_ident() {
        let rows = read_self_ident("test/self_ident/expected.bed.gz");
        let local_rows = compute_local_seq_self_identity(&rows, None);
        let reader = BufReader::new(GzDecoder::new(
            File::open("test/local_ident/expected.bed.gz").unwrap(),
        ));
        for (line, row) in reader.lines().map_while(Result::ok).zip(local_rows) {
            assert_eq!(line, row.tsv());
        }
    }
}
