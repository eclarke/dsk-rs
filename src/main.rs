extern crate bio;
extern crate kmers;
extern crate clap;
extern crate log;
#[macro_use]
extern crate error_chain;

use kmers::KmerCounter;
use clap::{App, Arg, ArgMatches};

error_chain!{
    links {
        Kmers(kmers::errors::Error, kmers::errors::ErrorKind);
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
    let args = parse_args();
    // Kmer length
    let k: usize = args.value_of("k").unwrap().trim().parse().unwrap();
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
        _ => bail!("Invalid alphabet selection!")
    };

    if k <= kmers::max_small_k(&alphabet) {
        let counter = KmerCounter::for_small_k(k, &alphabet)?;
    } else {
        println!("Using large kmer counting implementation as k > {}", kmers::max_small_k(&alphabet));
        let counter = KmerCounter::for_large_k(k, &alphabet)?;
    }
    Ok(())
}

