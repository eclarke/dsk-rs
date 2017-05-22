extern crate clap;
extern crate bio;
extern crate bit_vec;
extern crate byteorder;
extern crate vec_map;
extern crate num_bigint;
extern crate num_traits;
extern crate num;
extern crate itertools;

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read};

pub mod kmers;
use kmers::KmerCounter;
//use kmers::{Kmers, Counter};

use clap::{App, Arg, ArgMatches};
use bio::alphabets::{self, RankTransform};
use bio::io::fastq;

use vec_map::VecMap;
use bit_vec::BitVec;

use num::FromPrimitive;
use num_bigint::BigUint;

fn parse_args<'a>() -> ArgMatches<'a> {
    let args = App::new("DSK in Rust")
        .arg(Arg::with_name("seqs").short("s").takes_value(true).required(true))
        .arg(Arg::with_name("k").short("k").takes_value(true).default_value("27"))
        .arg(Arg::with_name("M").short("m").takes_value(true).default_value("2"))
        .arg(Arg::with_name("D").short("d").takes_value(true).default_value("3"))
        .get_matches();
    return args
}

fn main() {
    let args = parse_args();
    // Kmer length
    let k: u32 = args.value_of("k").unwrap().trim().parse().unwrap();
    // Maximum memory size (in GB)
    let max_mem = args.value_of("M").unwrap().trim().parse::<f32>().unwrap() * 8e9;
    // Maximum disk space (in GB)
    let max_disk = args.value_of("D").unwrap().trim().parse::<f32>().unwrap() * 8e9;

    let seq_fp = args.value_of("seqs").expect("No file specified!");

    // Estimated number of kmers
    let v: usize;
    {
        let reader = fastq::Reader::from_file(seq_fp).expect("Could not open Fastq file");
        v = reader.records().map(|r| r.unwrap().seq().len() - (k as usize) + 1).sum()
    }

    // Kmer storage calculations
    let alphabet = alphabets::Alphabet::new(b"ATGC");
    let rank = alphabets::RankTransform::new(&alphabet);
    let bits = (rank.ranks.len() as f32).log2().ceil() as u32;
    // Number of bytes used to encode a k-mer
    let bytes = bits * k;

    // Number of iterations
    let log2k = ((2*k) as f32).log2().ceil().exp2();
    let iters = ((v as f32) * log2k/max_disk).ceil() as u32;
    // Number of partitions
    let parts = (((v as f32)*(log2k + 32f32))/(0.7*(iters as f32)*max_mem)).ceil() as u32;    

    println!("k:{}, v:{}, max mem:{}, max disk:{}", k, v, max_mem, max_disk);
    println!("iterations: {:.1}\npartitions: {:.1}", iters, parts);

    for i in 0..(iters) {
        write_big_kmers(k, iters, i, parts, seq_fp, &rank).expect("Error writing kmers to disk!");
    }
    for part in 0..parts {
        let map = read_big_kmers(part, bytes as usize);
        for (kmer, count) in map.iter() {
            println!("{}:{}", decode_kmer(kmer, k, &rank).unwrap_or(":(".to_owned()), count);
        }
    }
}

fn write_big_kmers(k: u32, 
        iters: u32, 
        i: u32, 
        parts: u32, 
        record_fp: &str, 
        rank: &RankTransform) -> io::Result<()> {
    let iters: num_bigint::BigUint = FromPrimitive::from_u32(iters).unwrap();
    let i: num_bigint::BigUint = FromPrimitive::from_u32(i).unwrap();
    let parts: num_bigint::BigUint = FromPrimitive::from_u32(parts).unwrap();
    let reader = fastq::Reader::from_file(record_fp)?;
    let records = reader.records();
    for record in records {
        let record = record?;
        let counter = KmerCounter::for_large_k(k as usize, &alphabets::dna::n_alphabet()); //Kmers<BigUint, _> = Kmers::count_kmers(record.seq().iter(), k, rank);
        for kmer in counter.count(record.seq().iter()) {
            let kmer = BigUint::from_bytes_be(&kmer);
            if &kmer % &iters == i {
                let j_big = (&kmer / &iters) % &parts;
                let j = num::ToPrimitive::to_usize(&j_big).unwrap();
                save_kmer(kmer.to_bytes_be(), j);
            }
        }
    }
    Ok(())
}

fn read_big_kmers(part: u32, bytes: usize) -> HashMap<BigUint, usize> {
    let mut map: HashMap<BigUint, usize> = HashMap::new();
    let mut reader = get_partition(part);
    let mut buf = Vec::with_capacity(bytes);
    while let Ok(bytes_read) = reader.read(&mut buf) {
        if bytes_read < bytes {
            break;
        }
        let kmer: BigUint = BigUint::from_bytes_be(&buf);
        let count = map.entry(kmer).or_insert(0);
        *count += 1;
    }
    return map
}

fn save_kmer(s: Vec<u8>, partition: usize) {
    unimplemented!();
}

fn get_partition(partition: u32) -> io::BufReader<File> {
    unimplemented!();
}

fn decode_kmer(s: &BigUint, k: u32, rank: &alphabets::RankTransform) -> Option<String> {
    let bits = (rank.ranks.len() as f32).log2().ceil() as u32;
    let mut b = BitVec::from_bytes(&s.to_bytes_be());
    println!("{:?}", b);
    if b.len() < (2*k as usize) {
        return None
    }
    let mut vk: VecMap<u8> = VecMap::new();
    for (k,v) in rank.ranks.iter() {
        vk.insert(*v as usize, k as u8);
    }
    let mut chars = Vec::new();
    for _ in 0..k {
        let mut x = BitVec::with_capacity(2);
        x.push(b.pop().expect("couln't get first bit"));
        x.push(b.pop().expect("couldn't get second bit"));
        let x = x.iter().rev().collect::<BitVec>().to_bytes().iter().nth(0).expect("uh oh") >> (8-bits);
        chars.push(*vk.get(x as usize).unwrap());
    }
    let chars = chars.into_iter().rev().collect();
    return Some(String::from_utf8(chars).expect("Couldn't translate to string"));
}

