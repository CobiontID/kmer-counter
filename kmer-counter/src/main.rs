extern crate needletail;
use clap::{App, Arg};
use fnv::FnvHashMap;
use needletail::{parse_fastx_file, Sequence};
use std::error::Error;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;
use std::str;
use itertools::Itertools;

use ndarray_npy::write_npy;
use ndarray::prelude::*;

fn cartesian_product(list: &Vec<u8>, n: i32) -> Vec<Vec<u8>> {
    let mut res = vec![];
    for &i in list {
        res.push(vec![i]);
    }

    for _ in 0..n-1 {
        let mut tmp = vec![];
        for r in res {
            for &el in list {
                let mut tmp_el = r.clone();
                tmp_el.push(el);
                tmp.push(tmp_el);
            }
        }
        res = tmp;
    }
    res
}

fn rev_comp(kmer: &Vec<u8>) -> Vec<u8> {
    let comp: std::collections::HashMap<_, _> = [(65, 84), (84, 65), (67, 71), (71, 67)].iter().cloned().collect();
    let mut k = kmer.clone();
    k.reverse();
    k.iter().map(|x| comp[x]).collect::<Vec<_>>()
}

/// Converts nested `Vec`s to a 2-D array by cloning the elements.
///
/// **Panics** if the length of any axis overflows `isize`, if the
/// size in bytes of all the data overflows `isize`, or if not all the
/// rows have the same length.
fn vec_to_array<T: Clone>(v: Vec<Vec<T>>) -> Array2<T> {
    if v.is_empty() {
        return Array2::from_shape_vec((0, 0), Vec::new()).unwrap();
    }
    let nrows = v.len();
    let ncols = v[0].len();
    let mut data = Vec::with_capacity(nrows * ncols);
    for row in &v {
        assert_eq!(row.len(), ncols);
        data.extend_from_slice(&row);
    }
    Array2::from_shape_vec((nrows, ncols), data).unwrap()
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("K-mer counter")
        .version("0.1.0")
        .author("Claudia C. Weber <cw21@sanger.ac.uk>>")
        .about("Tally nucleotide counts in multi-entry fasta")
        .arg(
            Arg::with_name("file")
                .short("f")
                .long("file")
                .takes_value(true)
                .required(true)
                .help("Fasta file to tally."),
        )
        .arg(
            Arg::with_name("klength")
                .short("k")
                .long("klength")
                .takes_value(true)
                .default_value("4")
                .help("K-mer length"),
        )
        .arg(Arg::with_name("out")
                 .short("o")
                 .long("out")
                 .takes_value(true)
                 .required(true)
                 .default_value("counter_output.npy")
                 .help("Output file name."))
        .arg(Arg::with_name("ids")
                 .short("i")
                 .long("ids")
                 .takes_value(true)
                 .required(true)
                 .default_value("ids.txt")
                 .help("File to write identifiers to"))
        .arg(Arg::with_name("collapse")
                 .short("c")
                 .long("collapse")
                 .takes_value(true)
                 .required(true)
                 .default_value("1")
                 .help("Canonicalize k-mers (default 1 = True"))
        .get_matches();

    let filename = matches.value_of("file").unwrap();
    let k = matches.value_of("klength").unwrap();
    let out = matches.value_of("out").unwrap();
    let ids_out = matches.value_of("ids").unwrap();
    let canon = matches.value_of("collapse").unwrap().parse::<i32>().unwrap();

    println!("K-mer length: {:#?}", k);
    println!("File: {:#?}", filename);
    println!("Output file: {:#?}", out);

    let idfile = File::create(&ids_out).unwrap();
    let mut idfile = BufWriter::new(idfile);

    // A, C, G, T
    let acgt = [65, 67, 71, 84];
    let product = cartesian_product(&acgt.to_vec(), k.parse::<i32>().unwrap());

    let mut reader = parse_fastx_file(&filename).expect("valid path/file");

    let mut k_counts: FnvHashMap<&[u8], i32> =  FnvHashMap::default();
    //
    if canon == 1 {
        for p in &product {
            let rev = rev_comp(&p.clone());
            let check = !(k_counts.contains_key(rev.as_slice()));
            if check {
              *k_counts.entry(p.as_slice()).or_insert(0) += 0;
                }
            }
        } else {
        for p in &product {
            *k_counts.entry(p.as_slice()).or_insert(0) += 0; 
            }

        }
    
    // get keys that weren't collapsed
    let mut keys = vec![];
    for (key,_) in k_counts.iter().sorted() {
            let key_copy = key.clone();
            keys.push(key_copy);
        }
    //Uncomment line below to print keys
    //println!("Retained keys: {:?}", keys);

    let mut array = vec![];

    while let Some(record) = reader.next() {
        let mut k_counts_it = k_counts.clone();

        let seqrec = record.expect("invalid record");
        let tag = str::from_utf8(seqrec.id()).unwrap();
        writeln!(idfile, "{}", tag).unwrap();
        let norm_seq = seqrec.normalize(false);
        // normalize to make sure all the bases are consistently capitalized and
        // that we remove the newlines since this is FASTA
        let rc = norm_seq.reverse_complement();
        if canon == 1 {
            for (_, kmer, _) in norm_seq.canonical_kmers(k.parse::<u8>().unwrap(), &rc) {
                *k_counts_it.entry(kmer).or_insert(0) += 1;
            }
        } else {
            for kmer in norm_seq.kmers(k.parse::<u8>().unwrap()) {
                *k_counts_it.entry(kmer).or_insert(0) += 1;
            }
        }
        let mut line = vec![];
        for key in &keys {
            line.push(k_counts_it[key]);
        }
        array.push(line);
    }
    let arr = vec_to_array(array);
    write_npy(&out, &arr)?;
    Ok(())
}
