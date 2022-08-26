# Tally k-mer counts in a fasta file

## Usage
`kmer-counter --file <fasta file> --ids <output file for sequence ids> --klength <k-mer length, default 4> --out <output file> --collapse <canonicalize 0: False, 1: True; Default 0>`

Takes a fasta file of nucleotide sequences a input, and returns a plain text file with the sequence identifiers and a .npy file containing a numpy array with the k-mer counts (as 32 bit integers). By default, the counts will be canonicalized (--collapsed 1).

The counter operates directly on bytes. The order of the k-mers in the output corresponds to the cartesian products of the nucleotides A, C, G and T for k-mer size k. A list of the uncollapsed keys can be obtained in Python with:

```Python
from itertools import product
k = 4
["".join(i) for i in product("ACGT", repeat=k)]
```

## Requirements

This program is suitable for processing large read sets on an ordinary computer. From v. 0.1.2, the required memory scales with the selected k-mer size, not the number of sequence records.

A set of ~200 million HiFi reads from *Viscum album* (European mistletoe) was processed in a single batch in ~30 hours with a peak memory use of 7MB (k = 4, canonicalized).


## Installation

- Install Rust (see https://www.rust-lang.org/tools/install)
- Download the source code from this repository
- Navigate to the directory containing Cargo.toml
- Run `cargo build --release` . The resulting binary will be located in `./target/release/kmer-counter` unless otherwise specified.

