use murmur3::murmur3_x64_128;
use std::io::Cursor;

/// mmh3 uses x64_128 variant
/// https://mmh3.readthedocs.io/en/stable/quickstart.html#basic-usage
pub fn get_hash(value: &[u8], seed: Option<u32>) -> u64 {
    murmur3_x64_128(&mut Cursor::new(value), seed.unwrap_or_default()).unwrap() as u64
}
