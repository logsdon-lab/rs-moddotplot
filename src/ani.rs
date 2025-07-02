use ahash::RandomState;
use core::str;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{common::AIndexSet, io::Row};

const BASES_TO_REMOVE: [char; 11] = ['R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N'];

fn remove_ambiguous_bases(mod_list: &mut AIndexSet<usize>, k: usize, seed: Option<u64>) {
    let rng = seed
        .map(|seed| RandomState::with_seeds(seed, seed, seed, seed))
        .unwrap_or_default();

    // https://users.rust-lang.org/t/fill-string-with-repeated-character/1121/3
    let kmers_to_remove: AIndexSet<usize> = BASES_TO_REMOVE
        .iter()
        .map(|b| rng.hash_one(std::iter::repeat_n(b, k).collect::<String>()) as usize)
        .collect();
    // Remove homopolymers of ambiguous nucleotides
    mod_list.retain(|m| !kmers_to_remove.contains(m));
}

/// Create self-identity matrix.
///
/// # Args:
/// * sequence
///     * Sequence as list of mmh3 kmers.
/// * window_size
///     * Window size.
/// * delta
///     * Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.
/// * k
///     * kmer length
/// * identity
///     * Identity cutoff threshold.
/// * ambiguous
///     * Preserve diagonal when handling strings of ambiguous homopolymers (eg. long runs of N's).
/// * modimizer
///     * Modimizer sketch size.
///     * A lower value will reduce the number of modimizers, but will increase performance.
///     * Must be less t
#[allow(clippy::too_many_arguments)]
pub(crate) fn create_self_matrix(
    sequence: Vec<usize>,
    window_size: usize,
    delta: f32,
    k: usize,
    identity: f32,
    ambiguous: bool,
    modimizer: usize,
    seed: Option<u64>,
) -> Vec<Vec<f32>> {
    // Restrict sequence sparsity to powers of 2.
    let sequence_sparsity = window_size as f32 / modimizer as f32;
    let sequence_sparsity = if window_size / modimizer <= modimizer {
        2.0f32.powf(sequence_sparsity.log2())
    } else {
        2.0f32.powf((sequence_sparsity - 1.0).log2() + 1.0)
    };
    let sketch_size = window_size / sequence_sparsity as usize;

    let no_neighbors = partition_overlaps(&sequence, window_size, 0.0, k);
    // If considering neighboring bins, need to find overlaps.
    let neighbors = if delta > 0.0 {
        partition_overlaps(&sequence, window_size, delta, k)
    } else {
        no_neighbors.clone()
    };

    let neighbors_mods = convert_to_modimizers(
        &neighbors,
        sequence_sparsity as usize,
        ambiguous,
        k,
        sketch_size,
        seed,
    );
    let no_neighbors_mods = convert_to_modimizers(
        &no_neighbors,
        sequence_sparsity as usize,
        ambiguous,
        k,
        sketch_size,
        seed,
    );
    self_containment_matrix(no_neighbors_mods, neighbors_mods, k, identity, ambiguous)
}

fn partition_overlaps(
    sequence: &[usize],
    window_size: usize,
    delta: f32,
    k: usize,
) -> Vec<&[usize]> {
    let mut kmer_list: Vec<&[usize]> = vec![];
    let kmer_to_genomic_coordinate_offset = window_size - k + 1;
    let delta_offset = window_size as f32 * delta;

    // Set the first window to contains win-k+1 kmers
    // Use round_ties even to emulate banker's rounding in python
    let starting_end_index =
        (kmer_to_genomic_coordinate_offset as f32 + delta_offset).round_ties_even() as usize;
    kmer_list.push(&sequence[0..starting_end_index]);
    let mut index = window_size - k + 1;

    // Set normal windows.
    let sequence_length = sequence.len();
    let end_of_sequence = sequence_length - window_size;
    while index <= end_of_sequence {
        let start_index = index + 1;
        let end_index = window_size + index + 1;

        let delta_start_index = (start_index as f32 - delta_offset).round_ties_even() as usize;
        let mut delta_end_index = (end_index as f32 + delta_offset).round_ties_even() as usize;
        if delta_end_index > sequence_length {
            delta_end_index = sequence_length
        };

        if let Some(kmer_set) = sequence.get(delta_start_index..delta_end_index) {
            kmer_list.push(kmer_set);
        } else {
            kmer_list.push(&sequence[delta_start_index..sequence_length]);
        }
        index += window_size;
    }

    // set the last window to get the remainder
    if index <= sequence_length - 2 {
        let final_start_index = (index as f32 + 1.0 - delta_offset).round_ties_even() as usize;
        kmer_list.push(&sequence[final_start_index..sequence_length]);
    }

    // Test that last value was added on correctly
    if let Some(Some(p)) = kmer_list.last().map(|p| p.last()) {
        assert!(*p == sequence[sequence.len() - 1])
    }
    kmer_list
}

fn populate_modimizers(
    partition: &[usize],
    sparsity: usize,
    ambiguous: bool,
    expectation: usize,
    k: usize,
    seed: Option<u64>,
) -> AIndexSet<usize> {
    let mut mod_set = AIndexSet::default();
    for kmer in partition {
        if kmer % sparsity == 0 {
            mod_set.insert(*kmer);
        }
    }
    if !ambiguous {
        remove_ambiguous_bases(&mut mod_set, k, seed);
    }
    // Decrease sparsity untile expecatition met.
    if (mod_set.len() < expectation / 2) && sparsity > 1 {
        populate_modimizers(partition, sparsity / 2, ambiguous, expectation, k, seed);
    }
    mod_set
}

fn convert_to_modimizers(
    kmer_list: &[&[usize]],
    sparsity: usize,
    ambiguous: bool,
    k: usize,
    expectation: usize,
    seed: Option<u64>,
) -> Vec<AIndexSet<usize>> {
    kmer_list
        .into_par_iter()
        .map(|prt| populate_modimizers(prt, sparsity, ambiguous, expectation, k, seed))
        .collect()
}

pub(crate) fn convert_matrix_to_bed(
    matrix: Vec<Vec<f32>>,
    window_size: usize,
    id_threshold: f32,
    query_name: &str,
    reference_name: &str,
    self_identity: bool,
) -> Vec<Row> {
    let mut bed: Vec<Row> = vec![];
    let (rows, cols) = (matrix.len(), matrix.len());
    for (x, col) in matrix.iter().enumerate().take(rows) {
        for (y, value) in col.iter().enumerate().take(cols) {
            if (self_identity && x <= y) && *value / 100.0 >= id_threshold {
                let start_x = x * window_size + 1;
                let end_x = (x + 1) * window_size;
                let start_y = y * window_size + 1;
                let end_y = (y + 1) * window_size;
                bed.push(Row {
                    query_name: query_name.to_owned(),
                    query_start: start_x,
                    query_end: end_x,
                    reference_name: reference_name.to_owned(),
                    reference_start: start_y,
                    reference_end: end_y,
                    perc_id_by_events: *value,
                });
            }
        }
    }
    bed
}

/// Calculate the binomial distance based on containment and kmer values.
///
/// # Args
/// * containment_value
///     * The containment value.
/// * kmer_value
///     * The k-mer value.
///
/// # Returns
/// * The binomial distance.
fn binomial_distance(containment_value: f32, kmer_value: usize) -> f32 {
    containment_value.powf(1.0 / kmer_value as f32)
}

/// Calculate the containment neighbors based on four sets and an identity threshold.
///
/// Args:
/// * set1:
///     * The first set.
/// * set2:
///     * The second set.
/// * set3:
///     * The third set.
/// * set4:
///     * The fourth set.
/// * identity:
///     * The identity threshold.
/// * k:
///     * Kmer value.
///
/// Returns:
/// * The containment neighbors value.
fn containment_neighbors(
    set1: &AIndexSet<usize>,
    set2: &AIndexSet<usize>,
    set3: &AIndexSet<usize>,
    set4: &AIndexSet<usize>,
    identity: f32,
    k: usize,
) -> f32 {
    let len_a = set1.len();
    let len_b = set2.len();

    let intersection_a_b_prime = set1.intersection(set4).count();
    // If len_a is zero, handle it by setting containment_a_b_prime to a default value
    let containment_a_b_prime = if len_a != 0 {
        intersection_a_b_prime as f32 / len_a as f32
    } else {
        0.0
    };

    if binomial_distance(containment_a_b_prime, k) < identity / 100.0 {
        0.0
    } else {
        let intersection_a_prime_b = set2.intersection(set3).count();
        // If len_b is zero, handle it by setting containment_a_prime_b to a default value
        let containment_a_prime_b = if len_b != 0 {
            intersection_a_prime_b as f32 / len_b as f32
        } else {
            0.0
        };

        if containment_a_b_prime > containment_a_prime_b {
            containment_a_b_prime
        } else {
            containment_a_prime_b
        }
    }
}

/// Create a self-containment matrix based on containment similarity calculations.
///
/// # Args
/// * mod_set
///     * A list of sets representing elements.
/// * mod_set_neighbors
///     * A list of sets representing neighbors for each element.
/// * k
///     * A parameter for containment similarity calculation.
///
/// # Returns
/// * An ndarray representing the self-containment matrix.
fn self_containment_matrix(
    mod_set: Vec<AIndexSet<usize>>,
    mod_set_neighbors: Vec<AIndexSet<usize>>,
    k: usize,
    identity: f32,
    ambiguous: bool,
) -> Vec<Vec<f32>> {
    let n = mod_set.len();
    let mut containment_matrix: Vec<Vec<f32>> = vec![vec![0.0; n]; n];

    for w in 0..n {
        containment_matrix[w][w] = 100.0;

        if mod_set[w].is_empty() && !ambiguous {
            containment_matrix[w][w] = 0.0;
        }

        for r in (w + 1)..n {
            let c_hat = binomial_distance(
                containment_neighbors(
                    &mod_set[w],
                    &mod_set[r],
                    &mod_set_neighbors[w],
                    &mod_set_neighbors[r],
                    identity,
                    k,
                ),
                k,
            );
            containment_matrix[r][w] = c_hat * 100.0;
            containment_matrix[w][r] = c_hat * 100.0;
        }
    }

    containment_matrix
}

#[cfg(test)]
mod test {
    use ahash::RandomState;

    use crate::{
        ani::{partition_overlaps, remove_ambiguous_bases},
        common::AIndexSet,
    };

    #[test]
    fn test_remove_ambiguous_bases() {
        let seed = 42;
        let rng = RandomState::with_seeds(seed, seed, seed, seed);

        let mut mod_list =
            AIndexSet::from_iter(["ATGC", "NNNN", "DDDD"].map(|kmer| rng.hash_one(kmer) as usize));

        let expected_mod_list =
            AIndexSet::from_iter(["ATGC"].map(|kmer: &'static str| rng.hash_one(kmer) as usize));

        remove_ambiguous_bases(&mut mod_list, 4, Some(seed));

        assert_eq!(mod_list, expected_mod_list)
    }

    #[test]
    fn test_partition_overlaps() {
        let seq = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        assert_eq!(
            partition_overlaps(&seq, 3, 0.5, 2),
            [&seq[0..4], &seq[2..8], &seq[4..10], &seq[8..]].to_vec()
        );
    }
}
