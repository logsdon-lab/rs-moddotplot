/// Rust implementation of ModDotPlot
/// https://github.com/marbl/ModDotPlot/commit/50ecda4eff91acd00584090afd380d4a355be7aa
mod ani;
mod cfg;
mod ident;
mod io;

pub use cfg::SelfIdentConfig;
pub use ident::{compute_self_identity, compute_seq_self_identity};
pub use io::Row;
