/// Rust implementation of ModDotPlot
/// https://github.com/marbl/ModDotPlot/commit/50ecda4eff91acd00584090afd380d4a355be7aa
mod ani;
mod cfg;
mod ident;
mod io;
mod common;

pub use cfg::{LocalSelfIdentConfig, SelfIdentConfig};
pub use ident::{
    compute_local_seq_self_identity, compute_self_identity, compute_seq_self_identity,
};
pub use io::{LocalRow, Row};
