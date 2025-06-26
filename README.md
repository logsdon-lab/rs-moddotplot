# rs-ModDotPlot
Rust API implementation of the [`ModDotPlot`](https://github.com/marbl/ModDotPlot) ANI algorithm.

> [!NOTE]
> Uses ahash instead of MurmurHash3. Mostly due to lack of crates. Some [upside](https://github.com/tkaitchuck/aHash/blob/master/compare/readme.md#comparison-with-other-hashers)? 

## Usage
```bash
cargo add --git https://github.com/koisland/rs-moddotplot.git
```

```rust
// Compute the self-sequence identity of the CHM13 chr1 centromere.
// Write bedfile to stdout.
use std::io::{stdout, BufWriter, Write};
use rs_moddotplot::{compute_self_identity, Row};

fn main() {
    let bed = compute_self_identity("data/chm13_chr1.fa", None);
    let mut fh = BufWriter::new(stdout());
    writeln!(&mut fh, "{}", Row::header()).unwrap();
    for row in bed {
        writeln!(&mut fh, "{}", row.tsv()).unwrap();
    }
}
```

## Example
Self-identity
```bash
cargo run --example ident --release -- data/chm13_chr1.fa
```

Local self-identity
```bash
cargo run --example local_ident --release -- data/chm13_chr1.fa
```

## TODO:
* [ ] - More tests and docstrings.
* [ ] - Improve error handling and type casts.
