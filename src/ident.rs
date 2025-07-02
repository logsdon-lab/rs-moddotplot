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
        if traveled.contains(&(*x, y)) {
            continue;
        }
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
            // *
            //  x
            //    *
            let up_left = aln_mtx
                .get(&(x + 1))
                .and_then(|col| y.checked_sub(1).and_then(|y| col.get(&y)));
            let down_right = x
                .checked_sub(1)
                .and_then(|x| aln_mtx.get(&x))
                .and_then(|col| col.get(&(y + 1)));

            if up_left.is_none() && down_right.is_none() {
                max_x = x;
                continue;
            }
            let Some(col) = aln_mtx.get(&x) else {
                // Update x since we've gone into region.
                max_x = x;
                continue;
            };
            let Some(ident) = col.get(&y) else {
                continue;
            };
            // Add next positions to queue.
            // *
            // x *
            positions.push_back((x, y + 1));
            positions.push_back((x + 1, y));
            idents.push(*ident);
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
    use std::{
        fs::File,
        io::{BufRead, BufReader, BufWriter, Write},
        path::Path,
    };

    use crate::{compute_self_identity, LocalRow, Row};

    use super::compute_group_seq_self_identity;

    #[test]
    fn test_self_ident() {
        let path_outfile = Path::new("rows.tsv");
        let rows = if let Ok(mut new_file) = File::create_new(path_outfile).map(BufWriter::new) {
            let rows = compute_self_identity(
                "data/HG00438_chr3_HG00438#1#CM089169.1_89902259-96402509.fa",
                None,
            );
            for r in rows.iter() {
                writeln!(new_file, "{}", r.tsv()).unwrap();
            }
            rows
        } else {
            let reader = BufReader::new(File::open(path_outfile).unwrap());
            let mut rows = vec![];
            for line in reader.lines() {
                let line = line.unwrap();
                let [qname, qst, qend, rname, rst, rend, ident] =
                    line.trim().split('\t').collect::<Vec<&str>>()[..]
                else {
                    panic!("Invalid columns.")
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
        let grouped_rows = compute_group_seq_self_identity(&rows);
        assert_eq!(
            vec![
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 89983,
                    end: 104980,
                    avg_perc_id_by_events: 87.974396
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 209959,
                    end: 219957,
                    avg_perc_id_by_events: 90.16794
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 394922,
                    end: 404920,
                    avg_perc_id_by_events: 95.44085
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 419917,
                    end: 429915,
                    avg_perc_id_by_events: 87.34493
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 434914,
                    end: 489903,
                    avg_perc_id_by_events: 81.49909
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 1449711,
                    end: 1469707,
                    avg_perc_id_by_events: 80.896675
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 1499701,
                    end: 1509699,
                    avg_perc_id_by_events: 80.057755
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 1549691,
                    end: 1559689,
                    avg_perc_id_by_events: 94.243454
                },
                // live alpha-satellite array.
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 1574686,
                    end: 3429315,
                    avg_perc_id_by_events: 99.29863
                },
                // hsat1a.
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 3439313,
                    end: 5168967,
                    avg_perc_id_by_events: 96.924835
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 5173966,
                    end: 5193962,
                    avg_perc_id_by_events: 94.55426
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 5198961,
                    end: 5368927,
                    avg_perc_id_by_events: 89.79428
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 5398921,
                    end: 5408919,
                    avg_perc_id_by_events: 88.830315
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 5603880,
                    end: 5643872,
                    avg_perc_id_by_events: 98.12903
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 5648871,
                    end: 5838833,
                    avg_perc_id_by_events: 92.44347
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 6118777,
                    end: 6133774,
                    avg_perc_id_by_events: 81.78596
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 6208759,
                    end: 6218757,
                    avg_perc_id_by_events: 79.344925
                },
                LocalRow {
                    chrom: "HG00438_chr3_HG00438#1#CM089169.1:89902259-96402509".to_owned(),
                    start: 6288743,
                    end: 6318737,
                    avg_perc_id_by_events: 80.197556
                }
            ],
            grouped_rows
        )
    }
}
