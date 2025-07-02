#[derive(Debug, Clone)]
pub struct SelfIdentConfig {
    /// Window size in base pairs to calculate self-identity.
    pub window_size: usize,
    /// Kmer size.
    pub k: usize,
    /// Identity threshold.
    pub id_threshold: f32,
    /// Fraction of neighboring partition to include in identity estimation.
    /// * Must be between 0 and 1
    /// * Use > 0.5 is not recommended.
    pub delta: f32,
    /// Modimizer sketch size.
    ///  * A lower value will reduce the number of modimizers, but will increase performance.
    ///  * Must be less t
    pub modimizer: usize,
    /// Seed for ahash RandomState.
    pub seed: Option<u64>,
}

impl Default for SelfIdentConfig {
    fn default() -> Self {
        Self {
            window_size: 5000,
            k: 21,
            id_threshold: 0.86,
            delta: 0.5,
            modimizer: 1000,
            seed: Some(12315778783787232),
        }
    }
}
/// Configuration of local self-identity calculation.
///
/// Example:
/// * With `n_bins` of 2 and `window` size of 5.
///     * `*` is the calculated aln band
///     * `+` is self aln.
/// * This is the window calculated.
///      ```
///      // 4 * * *   +
///      // 3 * *   +
///      // 2 *   +
///      // 1   +
///      // 0 +
///      ```
///     * `*` values are averaged.
#[derive(Debug, Clone)]
pub struct LocalSelfIdentConfig {
    /// Window size in base pairs to calculate local self-identity.
    pub window_size: usize,
    /// Number of bins parallel to the self-identity diagonal to calculate local self-identify.
    /// * Larger values will consider longer-distance relationships.
    pub n_bins: usize,
    /// Number of bins parallel to the self-identity diagonal to ignore.
    pub ignore_bins: usize,
}

impl Default for LocalSelfIdentConfig {
    fn default() -> Self {
        Self {
            window_size: 5000,
            n_bins: 5,
            ignore_bins: 2,
        }
    }
}
