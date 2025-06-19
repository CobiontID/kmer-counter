# kmer-counter

This k-mer counter, based on [Needletail](https://github.com/onecodex/needletail)'s efficient FASTA parser is designed to tally k-mer counts for large sequencing read sets. It was written with downstream processing in TensorFlow or NumPy in mind, and stores the results in a [npy](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html) file. The counter operates directly on bytes, saving overhead.

For some examples of the counter in action, see https://cobiontid.github.io/

## Usage
`kmer-counter --file <fasta file> --ids <output file for sequence ids> --klength <k-mer length, default 4> --out <output file> --collapse <canonicalize 0: False, 1: True; Default 1>`

### Inputs and outputs
The program takes a fasta file of nucleotide sequences a input, and returns a plain text file with the sequence identifiers and a .npy file containing a numpy array with the k-mer counts (as 32 bit integers). The input fasta file may be gzip-ed.

#### Structure of the output

The rows in the output file correspond to sequence records, and the columns correspond to k-mer keys:

| | k1 | k2 | k3 | k4 |
|-|-|-|-|-|
| **seq 1** | | | | |
| **seq 2** | | | | |

The order of the sequence records remains the same as in the input fasta. By default, the counts will be canonicalized (--collapsed 1), so each k-mer and its reverse complement are assigned the same key. The order of the k-mers in the output corresponds to the cartesian products of the nucleotides A, C, G and T for k-mer size `k`. For reference, a list of the uncollapsed keys for a given value of `k` can be obtained in Python with:

```Python
from itertools import product
k = 4
["".join(i) for i in product("ACGT", repeat=k)]
```
### Notes
Given the structure of the output and the hashmap used while counting, this tool is not intended for counting very large k-mers that are likely to produce a sparse count matrix. Each possible k-mer is assigned a slot, regardless of whether that k-mer is subsequently observed in the input data. If you are interested in the number of distinct k-mers found in each sequence record, consider using [unique-kmer-counts](https://github.com/CobiontID/unique-kmer-counts) instead.

## Requirements

This program is suitable for processing large read sets on an ordinary computer. From v. 0.1.2, the required memory scales with the selected k-mer size (4^k for the non-canonicalized case), not the number of sequence records.

A set of ~200 million HiFi reads from *Viscum album* (European mistletoe) was processed in a single batch in ~30 hours with a peak memory use of 7MB (k = 4, canonicalized).


## Installation

- Install Rust (see https://www.rust-lang.org/tools/install)
- Download the source code from this repository
- Navigate to the directory containing Cargo.toml
- Run `cargo build --release` . The resulting binary will be located in `./target/release/kmer-counter` unless otherwise specified.

## Citation
If you use the kmer counter in your work, please cite https://doi.org/10.1093/g3journal/jkae187
