use ahash::RandomState;
use indexmap::{IndexMap, IndexSet};

pub type AIndexMap<K, V> = IndexMap<K, V, RandomState>;
pub type AIndexSet<K> = IndexSet<K, RandomState>;
