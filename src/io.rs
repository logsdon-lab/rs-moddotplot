use std::{fs::File, io::BufReader, path::Path, str};

use ahash::RandomState;
use noodles::fasta;
use rayon::iter::{ParallelBridge, ParallelIterator};

/// Sequence identity matrix row.
#[derive(Debug, PartialEq, Clone)]
pub struct Row {
    pub query_name: String,
    pub query_start: usize,
    pub query_end: usize,
    pub reference_name: String,
    pub reference_start: usize,
    pub reference_end: usize,
    /// Percent identity by events.
    pub perc_id_by_events: f32,
}

impl Row {
    pub fn header() -> &'static str {
        "#query_name\tquery_start\tquery_end\treference_name\treference_start\treference_end\tperID_by_events"
    }

    pub fn tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_start,
            self.query_end,
            self.reference_name,
            self.reference_start,
            self.reference_end,
            self.perc_id_by_events
        )
    }
}

/// Local sequence identity matrix row.
#[derive(Debug, PartialEq, Clone)]
pub struct LocalRow {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    /// Average percent identity by events.
    pub avg_perc_id_by_events: f32,
}

impl LocalRow {
    pub fn header() -> &'static str {
        "chrom\tstart\tend\tavg_perc_id_by_events"
    }

    pub fn tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}",
            self.chrom, self.start, self.end, self.avg_perc_id_by_events
        )
    }
}

/// Generate kmers from a sequence string with k length.
/// Uses ahash instead of murmurhash3.
pub(crate) fn generate_kmers_from_fasta(seq: &str, k: usize, seed: Option<u64>) -> Vec<usize> {
    let n = seq.len();
    let mut kmers = Vec::with_capacity(n - k + 1);
    // Wtf? Why does with_seed produce difference results with the same seed but with_seeds doesn't?
    // Due to runtime rng. Gets rng from operating system.
    // https://users.rust-lang.org/t/inexplicable-nondeterministic-behavior-in-scientific-computing-code/109400/13
    let rng = seed
        .map(|seed| RandomState::with_seeds(seed, seed, seed, seed))
        .unwrap_or_default();
    for i in 0..(n - k + 1) {
        let kmer = &seq[i..i + k].to_uppercase();
        let rc_kmer: String = kmer
            .chars()
            .map(|n| match n {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => n,
            })
            .rev()
            .collect();
        let fh = rng.hash_one(kmer) as usize;
        let rc = rng.hash_one(rc_kmer) as usize;

        kmers.push(if fh < rc { fh } else { rc });
    }
    kmers
}

/// Read all kmers from a fasta file.
pub(crate) fn read_kmers(
    filename: impl AsRef<Path>,
    k: usize,
    seed: Option<u64>,
) -> Vec<(String, Vec<usize>)> {
    let buf = BufReader::new(File::open(filename).unwrap());
    let mut reader = fasta::Reader::new(buf);
    reader
        .records()
        .par_bridge()
        .map(|rec| {
            let rec = rec.unwrap();
            (
                String::from_utf8(rec.name().to_vec()).unwrap(),
                generate_kmers_from_fasta(
                    str::from_utf8(rec.sequence().as_ref()).unwrap(),
                    k,
                    seed,
                ),
            )
        })
        .collect()
}
