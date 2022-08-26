extern crate needletail;
use clap::{Command, Arg};
use fnv::FnvHashMap;
use needletail::{parse_fastx_file, Sequence};
use std::error::Error;
use std::io::prelude::*;
use std::io::{BufWriter};//,Seek,SeekFrom
use std::fs::File;
use std::str;
use itertools::Itertools;

use mktemp::Temp;

//use genawaiter::{stack::let_gen, yield_};

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

// Get reverse complement of kmers
// Unicode: A -> 65, C -> 67, G -> 71, T -> 84
fn rev_comp(kmer: &Vec<u8>) -> Vec<u8> {
    let comp: std::collections::HashMap<_, _> = [(65, 84), (84, 65), (67, 71), (71, 67)].iter().cloned().collect();
    let mut k = kmer.clone();
    k.reverse();
    k.iter().map(|x| comp[x]).collect::<Vec<_>>()
}


fn create_header(_nrows: i32, _ncols: usize) -> Vec<u8> {
    let mut header: Vec<u8> = vec![];
    header.extend(&b"{'descr': "[..]);
    header.extend("\'<i4\'".to_string().as_bytes());
    header.extend(&b", 'fortran_order': False, 'shape': "[..]);
    //let shape_pos = header.len() + 10;
    let shape = format!("{:?}", (_nrows,_ncols));
    header.extend(shape.as_bytes());
    //header.extend(FILLER);
    header.extend(&b"}"[..]);
    header
    //(header, shape_pos)
}


fn write_header(_nrows: i32, _ncols: usize, out: &str, tmp_path: mktemp::Temp) -> std::io::Result<()> {
    let mut fw = BufWriter::new(File::create(&out)?);
    fw.write_all(&[0x93u8])?;
    fw.write_all(b"NUMPY")?;
    fw.write_all(&[0x01u8, 0x00])?;

    let header = create_header(_nrows, _ncols); //let (header, shape_pos)

    let mut padding: Vec<u8> = vec![];
    padding.extend(&::std::iter::repeat(b' ').take(15 - ((header.len() + 10) % 16)).collect::<Vec<_>>());
    padding.extend(&[b'\n']);

    let len = header.len() + padding.len();
    assert! (len <= ::std::u16::MAX as usize);
    assert_eq!((len + 10) % 16, 0);

    fw.write_all(&(len as u16).to_le_bytes())?;
    fw.write_all(&header)?;
    // Padding to 8 bytes
    fw.write_all(&padding)?;
    fw.flush()?;

    let mut bin_src = File::open(&tmp_path)?;
    std::io::copy(&mut bin_src, &mut fw)?;

    Ok(())
}


fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("K-mer counter")
        .version("0.1.2")
        .author("Claudia C. Weber <cw21@sanger.ac.uk>>")
        .about("Tally nucleotide counts in multi-entry fasta")
        .arg(
            Arg::new("file")
                .short('f')
                .long("file")
                .takes_value(true)
                .required(true)
                .help("Fasta file to tally."),
        )
        .arg(
            Arg::new("klength")
                .short('k')
                .long("klength")
                .takes_value(true)
                .default_value("4")
                .help("K-mer length"),
        )
        .arg(Arg::new("out")
                 .short('o')
                 .long("out")
                 .takes_value(true)
                 .required(false)
                 .default_value("counter_output.npy")
                 .help("Output file name."))
        .arg(Arg::new("ids")
                 .short('i')
                 .long("ids")
                 .takes_value(true)
                 .required(false)
                 .default_value("ids.txt")
                 .help("File to write identifiers to"))
        .arg(Arg::new("collapse")
                 .short('c')
                 .long("collapse")
                 .takes_value(true)
                 .required(false)
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

    let _ncols = keys.len();
    println!("Number of keys: {:#?}", _ncols);

    //Uncomment line below to print keys
    //println!("Retained keys: {:?}", keys);

    let mut _nrows = 0;
    let tmp_path = Temp::new_file()?;

    {

    //let mut tmp = File::create(&tmp_path)?;
    let file = File::create(&tmp_path)?;
    let mut file = BufWriter::new(file);
    //let mut file = BufWriter::with_capacity(16384*4,file);

    while let Some(record) = reader.next() {
        _nrows += 1;
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
        //let mut line = vec![];
        for key in &keys {
            file.write_all(&k_counts_it[key].to_le_bytes())?;
            //line.push(k_counts_it[key]);
        }

    };
    file.flush()?;
    } 


    write_header(_nrows, _ncols, out, tmp_path)?;

    Ok(())
}
