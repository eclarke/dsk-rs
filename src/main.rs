extern crate kmers;
extern crate clap;
extern crate dsk;
extern crate env_logger;
#[macro_use] extern crate log;
#[macro_use] extern crate error_chain;

use clap::{App, Arg, ArgMatches};

error_chain! {
    links {
        App(::dsk::errors::Error, ::dsk::errors::ErrorKind);
    }
    foreign_links {
        Logging(::log::SetLoggerError);
    }
}

fn main() {
    env_logger::init().expect("Error initializing logger");
    if let Err(ref e) = run() {
        use std::io::Write;
        let stderr = &mut ::std::io::stderr();
        let errmsg = "Error writing to stderr";

        writeln!(stderr, "error: {}", e).expect(errmsg);

        for e in e.iter().skip(1) {
            writeln!(stderr, "caused by: {}", e).expect(errmsg);
        }

        if let Some(backtrace) = e.backtrace() {
            writeln!(stderr, "backtrace: {:?}", backtrace).expect(errmsg);
        }

        ::std::process::exit(1);
    }
}

fn parse_args<'a>() -> ArgMatches<'a> {
    let args = App::new("DSK in Rust")
        .arg(Arg::with_name("input")
                 .takes_value(true)
                 .required(true)
                 .index(1))
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
                .help("Input is in FASTQ format (default: FASTA)")
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

    let args = parse_args();
    let app = dsk::App::new(args)?;
    info!("\
Starting DSK with parameters: 
\t k:          {}
\t input:      {}
\t input fmt:  {:?}
\t iterations: {}
\t partitions: {}
\t workspace:  {:?}",
    app.k, app.input, app.format, app.iters, app.parts, app.workspace.path()
);

    let max_k = kmers::max_small_k(app.alphabet());
    if app.k <= max_k {
        info!("Using small kmer counter");
        let counter = app.write_small_kmers()?;
        app.count_kmers(&counter)?;
    } else {
        info!("Using large kmer counter");
        let counter = app.write_large_kmers()?;
        app.count_kmers(&counter)?;
    }

    Ok(())
}

