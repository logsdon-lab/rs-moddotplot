#!/bin/bash

set -euo pipefail

wd="test/comparison"
fa="data/chm13_chr1.fa"

mkdir -p "${wd}"

# This version
# Primary difference is colorscale (Finer-grain for visibility) and speed.
# Cannot be 1-to-1 as cannot set seed for ModDotPlot but boundaries of satellite regions are similar by eye.
cargo run --example plot_self_ident --release -- "${fa}" "${wd}/chm13_chr1_rs_moddotplot.png"

# ModDotPlot
# pip install moddotplot==0.9.9
# Cannot set colorscale. See https://github.com/marbl/ModDotPlot/issues/57
moddotplot static -f "${fa}" \
    -k 21 \
    -w 5000 \
    -o "${wd}"
