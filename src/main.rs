extern crate bio;
extern crate kmers;
extern crate clap;
#[macro_use]
extern crate error_chain;
extern crate tempdir;
extern crate byteorder;
extern crate bincode;

#[macro_use]
extern crate slog;
extern crate slog_term;
extern crate slog_async;
extern crate num;
extern crate num_bigint;
extern crate num_traits;

extern crate fastx;

pub mod dsk;

use slog::Drain;
use clap::{App, Arg, ArgMatches};

use dsk::errors::*;

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
        .arg(Arg::with_name("fastq")
                .help("input is in fastq format")
                .short("q")
                .long("fastq")
                .takes_value(false))
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

    // Set up logging
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator).build().fuse();
    let drain = slog_async::Async::new(drain).build().fuse();
    let log = slog::Logger::root(drain, o!());

    let args = parse_args();
    let app = dsk::App::new(args)?;
    info!(log, "{}", app);

    let max_k = kmers::max_small_k(app.alphabet());
    if app.k <= max_k {
        let counter = app.write_small_kmers()?;
        app.count_kmers(&counter)?;
    } else {
        let counter = app.write_large_kmers()?;
        app.count_kmers(&counter)?;
    }

    Ok(())
}

