# Tally k-mer counts in a fasta file

## Usage
`kmer-counter --file <fasta file> --ids <output file for sequence ids> --klength <k-mer length, default 4> --out <output file> --collapse <canonicalize 0: False, 1: True; Default 0>`

Takes a fasta file of nucleotide sequences a input, and returns a plain text file with the sequence identifiers and a numpy array with the k-mer counts (as 32 bit integers). By default, the counts will be canonicalized (--collapsed 1).

The counter operates directly on bytes. The order of the k-mers in the output corresponds to the cartesian products of the nucleotides A, C, G and T for k-mer size k. A list of the uncollapsed keys can be obtained in Python with:

```Python
from itertools import product
k = 4
["".join(i) for i in product("ACGT", repeat=k)]
```

## Dependencies
Rust (see https://www.rust-lang.org/tools/install)

- needletail = "0.4"
- fnv = "1.0.7"
- clap = ""
- ndarray-npy = "0.8"
- ndarray = "0.15.0"
- itertools = "0.10.0"
