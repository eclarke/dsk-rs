extern crate bio;
extern crate kmers;
extern crate clap;
#[macro_use]
extern crate error_chain;
extern crate tempdir;
extern crate byteorder;
extern crate bincode;
#[macro_use]
extern crate serde_derive;

extern crate serde;
extern crate serde_bytes;

#[macro_use]
extern crate slog;
extern crate slog_term;
extern crate slog_async;

use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::{OpenOptions, File};
use std::collections::HashMap;

use slog::Drain;
use kmers::KmerCounter;
use clap::{App, Arg, ArgMatches};
use bio::io::fastq;
use tempdir::TempDir;
use byteorder::{BigEndian, ReadBytesExt};
use serde_bytes::ByteBuf;

error_chain!{
    links {
        Kmers(kmers::errors::Error, kmers::errors::ErrorKind);
    }

    foreign_links {
        Io(::std::io::Error);
        Num(::std::num::ParseIntError);
        // Serialize(::serde::ser::Error);
    }
}



quick_main!(run);

fn parse_args<'a>() -> ArgMatches<'a> {
    let args = App::new("DSK in Rust")
        .arg(Arg::with_name("seqs")
                 .short("s")
                 .takes_value(true)
                 .required(true))
        .arg(Arg::with_name("k")
                 .help("Length of kmers to count")
                 .short("k")
                 .takes_value(true)
                 .default_value("27"))
        .arg(Arg::with_name("max_mem")
                 .help("Max RAM to use (in GB)")
                 .short("m")
                 .long("max_mem")
                 .takes_value(true)
                 .default_value("2"))
        .arg(Arg::with_name("max_disk")
                 .help("Max disk space to use (in GB)")
                 .short("d")
                 .long("max_disk")
                 .takes_value(true)
                 .default_value("3"))
        .arg(Arg::with_name("alphabet")
            .help("Nucleotide alphabet to use.\nDNA = ATGC; dna = ATGCatgc; dna+n = ATGCNatgcn; iupac = all IUPAC nucleotide symbols.\nUse the smallest applicable for your sequences; large alphabets slow down counting.\n")
            .short("a").long("alphabet")
            .possible_values(&["dna", "DNA", "dna+N", "iupac"])
            .takes_value(true)
            .default_value("DNA"))
        .get_matches();
    return args;
}

fn run() -> Result<()> {

    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let log = slog::Logger::root(drain, o!());

    let args = parse_args();
    // Kmer length
    let k: usize = args.value_of("k").unwrap().trim().parse()?;
    // Maximum memory size (in GB)
    let max_mem = args.value_of("max_mem")
        .unwrap()
        .trim()
        .parse::<f32>()
        .unwrap() * 8e9;
    // Maximum disk space (in GB)
    let max_disk = args.value_of("max_disk")
        .unwrap()
        .trim()
        .parse::<f32>()
        .unwrap() * 8e9;
    let alphabet = args.value_of("alphabet").unwrap();

    let alphabet = match alphabet {
        "dna" => bio::alphabets::dna::alphabet(),
        "DNA" => bio::alphabets::Alphabet::new(b"ATGC"),
        "dna+N" => bio::alphabets::dna::n_alphabet(),
        "iupac" => bio::alphabets::dna::iupac_alphabet(),
        _ => bail!("Invalid alphabet selection!"),
    };

    let seq_fp = args.value_of("seqs").unwrap();

    // Estimated number of kmers
    let n_kmers: usize;
    {
        debug!(log, "determining aggregate sequence length");
        let reader = fastq::Reader::from_file(seq_fp).chain_err(|| "couldn't open fastq file")?;
        // Bug: if k > sequence length, this breaks
        n_kmers = reader
            .records()
            .map(|r| r.expect("error parsing fastq record").seq().len() - (k as usize) + 1)
            .sum();
    }

    let log2k = ((2 * k) as f32).log2().ceil().exp2();
    // Number of iterations
    let iters = ((n_kmers as f32) * log2k / max_disk).ceil() as usize;
    // Number of partitions
    let parts = (((n_kmers as f32) * (log2k + 32f32)) / (0.7 * (iters as f32) * max_mem))
        .ceil() as usize;
    // Largest k for which we can use the small implementation
    let max_k = kmers::max_small_k(&alphabet);
    debug!(log, "program start"; "total kmers" => n_kmers, "partitions" => parts, "iterations" => iters, "max small k" => max_k);

    debug!(log, "creating temp files");
    let tmp_dir = TempDir::new("dsk_workspace").chain_err(|| "Couldn't create tempdir")?;
    let mut tmpfiles: Vec<std::path::PathBuf> = Vec::new();
    for p in 0..parts {
        let fp = tmp_dir.path().join(format!("kmer_bank_{}", p));
        tmpfiles.push(fp);
    }
    debug!(log, "tempfiles: {:?}", tmpfiles);
    if k <= max_k {
        let mut writers = create_writers(&tmpfiles, 10000)?;
        debug!(log, "using small kmer implementation");
        let counter = KmerCounter::for_small_k(k, &alphabet)?;
        debug!(log, "kmer size"; "bytes" => counter.bytes_per_kmer);
        for i in 0..iters {
            let reader = fastq::Reader::from_file(seq_fp)
                .chain_err(|| "Couldn't open fastq file for reading")?;
            let records = reader.records();
            for record in records {
                let record = record?;
                for kmer in counter.kmers(record.seq().iter()) {
                    let mut rdr = std::io::Cursor::new(&kmer);
                    let ukmer: usize = rdr.read_uint::<BigEndian>(kmer.len())
                        .chain_err(|| "couldn't convert to u64")? as usize;
                    if ukmer % iters == i {
                        let part = (ukmer / iters) % parts;
                        let mut file = &mut writers[part];
                        file.write(&kmer)
                            .chain_err(|| "couldn't serialize kmer")?;
                    }
                }
            }
        }
        debug!(log, "finished streaming to temp file");
        drop(writers);
        let mut readers = create_readers(&tmpfiles, 1000)?;
        let mut map: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut map_out = File::create("test.out")?;
        for part in 0..parts {
            debug!(log, "reading kmers"; "bank"=>part);
            let mut submap: HashMap<Vec<u8>, usize> = HashMap::new();
            let file = &mut readers[part];
            let mut kmer = vec![0; counter.bytes_per_kmer];
            while let Ok(bytes_read) = file.read(&mut kmer) {
                if bytes_read < counter.bytes_per_kmer {
                    break;
                }

                let count = submap.entry(kmer.clone()).or_insert(0);
                *count += 1;
            }
            debug!(log, "finished, extending primary map");
            map.extend(submap);
        }
        debug!(log, "serializing map to disk");
        bincode::serialize_into(&mut map_out, &mut map, bincode::Infinite);
        debug!(log, "finished!");

    } else {
        debug!(log, "using large kmer implementation");
        let counter = KmerCounter::for_large_k(k, &alphabet)?;
    }
    Ok(())
}

fn create_writers(paths: &Vec<std::path::PathBuf>, bufsize: usize) -> Result<Vec<BufWriter<File>>> {
    let mut writers = Vec::with_capacity(paths.len());
    for path in paths {
        let file = OpenOptions::new().read(true)
            .write(true)
            .create(true)
            .open(path)
            .chain_err(|| "couldn't create tempfile")?;
        writers.push(BufWriter::with_capacity(bufsize, file));
    }
    Ok(writers)
}

fn create_readers(paths: &Vec<std::path::PathBuf>, bufsize: usize) -> Result<Vec<BufReader<File>>> {
    let mut readers = Vec::with_capacity(paths.len());
    for path in paths {
        let file = OpenOptions::new().read(true)
            .write(true)
            .create(true)
            .open(path)
            .chain_err(|| "couldn't create tempfile")?;
        readers.push(BufReader::with_capacity(bufsize, file));
    }
    Ok(readers)
}
