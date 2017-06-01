extern crate bio;
extern crate kmers;
extern crate clap;
extern crate tempdir;
extern crate byteorder;
extern crate bincode;
extern crate num;
extern crate num_bigint;
extern crate num_traits;
extern crate fastx;
extern crate slog_stdlog;

#[macro_use] extern crate slog;
#[macro_use] extern crate error_chain;

use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::{BufWriter, BufReader, Write, Read};
use std::fmt;
use std::collections::HashMap;

use slog::Drain;

use bio::alphabets::{dna, Alphabet};
use clap::ArgMatches;
use tempdir::TempDir;

use num_bigint::{BigUint, ToBigUint};
use num_traits::ToPrimitive;
use kmers::KmerCounter;
use fastx::Records;

use bincode::{serialize_into, Infinite};

#[derive(Debug, Copy, Clone)]
pub enum SequenceFormat {
    Fasta,
    Fastq,
}

pub mod errors {
    use kmers;

    error_chain!{
        links {
            Kmers(kmers::errors::Error, kmers::errors::ErrorKind);
        }

        foreign_links {
            Io(::std::io::Error);
            ParseFloatError(::std::num::ParseFloatError);
            ParseIntError(::std::num::ParseIntError);
            Serialize(::bincode::ErrorKind);
        }
    }
}

use self::errors::*;

pub struct App {
    pub k: usize,
    pub input: String,
    pub format: SequenceFormat,
    pub iters: usize,
    pub parts: usize,
    pub alphabet: Alphabet,
    pub workspace: TempDir,
    pub output: String,
    tempfiles: Vec<PathBuf>,
    log: slog::Logger,
}

impl fmt::Display for App {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "DSK: k:{}\tinput:{}\titerations:{}\tpartitions:{}\talphabet size:{}",
            self.k, self.input, self.iters, self.parts, self.alphabet.len()
        )
    }
}

impl App {
    pub fn new<L: Into<Option<slog::Logger>>>(args: ArgMatches, logger: L) -> Result<Self> {
        let log = logger.into().unwrap_or(slog::Logger::root(slog_stdlog::StdLog.fuse(), o!()));
        let k: usize = args.value_of("k").unwrap().trim().parse()?;
        let input = args.value_of("input").unwrap();
        let format = match args.is_present("fastq") {
            true => SequenceFormat::Fastq,
            false => SequenceFormat::Fasta
        };

        // Maximum memory size (in GB)
        let max_mem = args.value_of("max_mem")
            .unwrap()
            .trim()
            .parse::<f32>()? * 8e9;
        
        // Maximum disk space (in GB)
        let max_disk = args.value_of("max_disk")
            .unwrap()
            .trim()
            .parse::<f32>()? * 8e9;

        debug!(log, "Determining total number of kmers...");
        let n_kmers: isize = records(input, format)?
            .map(|r| (r.seq.len() as isize) - (k as isize) + 1)
            .filter(|x| *x > 0).sum();
        debug!(log, ""; "kmers"=>n_kmers);
        let log2k = ((2 * k) as f32).log2().ceil().exp2();

        // Number of iterations
        let iters = ((n_kmers as f32) * log2k / max_disk).ceil() as usize;
        // Number of partitions
        let parts = (((n_kmers as f32) * (log2k + 32f32)) / (0.7 * (iters as f32) * max_mem))
            .ceil() as usize;

        let alphabet = match args.value_of("alphabet").unwrap() {
            "dna" => dna::alphabet(),
            "DNA" => Alphabet::new(b"ATGC"),
            "dna+N" => dna::n_alphabet(),
            "iupac" => dna::iupac_alphabet(),
            _ => bail!("Invalid alphabet selection!"),
        };

        let (tempdir, tempfiles) = tempfiles(parts)?;

        let outfile = args.value_of("outfile").unwrap();
        {
            // Check that we can create a writeable file at target location
            File::create(outfile).chain_err(|| format!("Couldn't create output file {}", outfile))?;
        }

        Ok(App {
            k, 
            input: String::from(input), 
            format, 
            iters, 
            parts, 
            alphabet, 
            tempfiles, 
            workspace: tempdir, 
            log, 
            output: String::from(outfile),
        })
    }

    pub fn alphabet(&self) -> &Alphabet {
        &self.alphabet
    }

    pub fn iters(&self) -> usize {
        self.iters
    }

    pub fn writers(&self) -> Result<Vec<BufWriter<File>>> {
        let mut writers = Vec::with_capacity(self.tempfiles.len());
        for path in &self.tempfiles {
            let file = File::create(path)
                .chain_err(|| "couldn't create tempfile")?;
            writers.push(BufWriter::new(file));
        }
        Ok(writers)   
    }

    pub fn readers(&self) -> Result<Vec<BufReader<File>>> {
        let mut readers = Vec::with_capacity(self.tempfiles.len());
        for path in &self.tempfiles {
            let file = File::open(path)
                .chain_err(|| "couldn't create tempfile")?;
            readers.push(BufReader::new(file));
        }
        Ok(readers)
    }


    pub fn write_small_kmers(&self) -> Result<KmerCounter<usize>> {
        let counter = KmerCounter::for_small_k(self.k, self.alphabet())?;
        let mut writers = self.writers()?;
        for i in 0..self.iters {
            debug!(self.log, "Writing kmers to disk"; "iteration"=>i+1);
            let records = records(&self.input, self.format)?;
            for record in records {
                for (bytes, kmer) in counter.kmers(record.seq.iter()) {
                    if kmer % self.iters == i {
                        let part = (kmer / self.iters) % self.parts;
                        let mut file = &mut writers[part];
                        file.write(&bytes).chain_err(|| "problem writing kmer to file")?;
                    }
                }
            }
        }
        Ok(counter)
    }

    pub fn write_large_kmers(&self) -> Result<KmerCounter<BigUint>> {
        let counter = KmerCounter::for_large_k(self.k, self.alphabet())?;
        let mut writers = self.writers()?;
        let big_iters = self.iters.to_biguint().unwrap();
        for i in 0..self.iters {
            debug!(self.log, "Writing kmers to disk"; "iteration"=>i+1);
            let records = records(&self.input, self.format)?;
            for record in records {
                for (bytes, kmer) in counter.kmers(record.seq.iter()) {
                    if &kmer % &big_iters == i.to_biguint().unwrap() {
                        let part = (&kmer / &big_iters) % self.parts.to_biguint().unwrap();
                        let part = part.to_usize().unwrap();
                        let mut file = &mut writers[part];
                        file.write(&bytes).chain_err(|| "problem writing kmer to file")?;
                    }
                }
            }
        }
        Ok(counter)
    }

    pub fn count_kmers<K>(&self, counter: &KmerCounter<K>) -> Result<HashMap<Vec<u8>, usize>> {
        let mut readers = self.readers()?;
        let mut map: HashMap<Vec<u8>, usize> = HashMap::new();

        for part in 0..self.parts {
            let infile = &mut readers[part];
            let mut kmer: Vec<u8> = vec![0; counter.bytes_per_kmer];
            while let Ok(bytes_read) = infile.read(&mut kmer) {
                if bytes_read < counter.bytes_per_kmer {
                    break;
                }
                let count = map.entry(kmer.clone()).or_insert(0);
                *count += 1;
            }
        }
        Ok(map)
    }

    pub fn write_map(&self, map: HashMap<Vec<u8>, usize>) -> Result<()> {
        let mut outfile = File::create(&self.output)?;
        serialize_into(&mut outfile, &map, Infinite).map_err(|e| *e)
            .chain_err(|| "error writing results to outfile")?;
        Ok(())
    }
    
}


fn tempfiles(n: usize) -> Result<(TempDir, Vec<PathBuf>)> {
    let tempdir = TempDir::new("dsk_workspace").chain_err(|| "couldn't create tempdir")?;
    let mut tempfiles: Vec<PathBuf> = Vec::new();
    for p in 0..n {
        let fp = tempdir.path().join(format!("kmer_bank_{}", p));
        tempfiles.push(fp);
    }
    Ok((tempdir, tempfiles))
}

fn records<P: AsRef<Path>>(path: P, format: SequenceFormat) -> Result<Records<BufReader<File>>> {
    let reader = BufReader::new(File::open(path).chain_err(|| "couldn't open sequence file")?);
    match format {
        SequenceFormat::Fasta => Ok(Records::from_fasta(reader)),
        SequenceFormat::Fastq => Ok(Records::from_fastq(reader))
    }
}