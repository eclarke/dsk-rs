extern crate clap;
extern crate bio;
extern crate bit_vec;
extern crate byteorder;
extern crate vec_map;

use std::collections::HashMap;

use clap::{App, Arg, ArgMatches};
use bio::alphabets;
use bio::io::fastq;

use vec_map::VecMap;
use byteorder::{BigEndian, WriteBytesExt};
use bit_vec::BitVec;

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
    // Pseudocode for DSK algorithm
    let S: Vec<Vec<u8>>; // Set of sequences
    let k: u32; // kmer size
    let M: f64; // memory size
    let D: f64; // disk space    

    let args = parse_args();
    k = args.value_of("k").unwrap().trim().parse().unwrap();
    M = args.value_of("M").unwrap().trim().parse::<f64>().unwrap() * 8e9;
    D = args.value_of("D").unwrap().trim().parse::<f64>().unwrap() * 8e9;

    let seq_fp = args.value_of("seqs").expect("No file specified!");
    let reader = fastq::Reader::from_file(seq_fp).expect("Could not open Fastq file");
    S = reader.records().map(|r| r.unwrap().seq().to_owned()).collect();

    let alphabet = alphabets::Alphabet::new(b"ATGC");
    let rank = alphabets::RankTransform::new(&alphabet);
    // Number of k-mers
    let v: usize = S.iter().map(|s| s.len() - (k as usize) + 1).sum();
    // Number of iterations
    let n_i = (((v as f64) * ((2*k) as f64).log2().ceil().exp2())/D).ceil();
    // Number of partitions
    let n_p = (((v as f64)*(((2*k) as f64).log2().ceil().exp2() + 32f64))/(0.7*n_i*M)).ceil();

    println!("k:{}, M:{}, D:{}", k, M, D);
    // println!("Max k: {}", std::usize::BITS / alphabet.len().log2());
    println!("iterations: {:.1}\npartitions: {:.1}", n_i, n_p);

    // decode_kmer(16807957669258353, k, &rank);
    // return {};
    for i in 0..(n_i as u32) {
        let mut lists = initialize_lists(n_p as usize);
        for s in &S {
            // let seq = s.as_bytes();
            for kmer in rank.qgrams(k as u32, s) {
                let hm = kmer as u64; //hash(&kmer);
                if hm % (n_i as u64) == i as u64 {
                    let j = ((hm / (n_i as u64)) % (n_p as u64)) as usize;
                    write_kmer(kmer, &mut lists[j]);
                }
            }
        }
        for j in 0..(n_p as usize) {
            let mut map: HashMap<usize, u32> = HashMap::default();
            for kmer in &lists[j] {
                let count = map.entry(*kmer).or_insert(0);
                *count += 1;
            }
            for (kmer, count) in map.iter() {
                println!("{}:{}", decode_kmer(*kmer, k, &rank), count);
            }
        }
    }

}

fn initialize_lists(n: usize) -> Vec<Vec<usize>> {
    vec![Vec::new(); n]
}

fn write_kmer(s: usize, list: &mut Vec<usize>) {
    list.push(s)
}

fn decode_kmer(s: usize, k: u32, rank: &alphabets::RankTransform) -> String {
    let bits = (rank.ranks.len() as f32).log2().ceil() as u32;
    // let mask: usize = (1 << (k * bits)) - 1;
    let mut v = vec![];
    v.write_u64::<BigEndian>(s as u64).unwrap();
    let mut b = BitVec::from_bytes(&v);
    let mut vk: VecMap<u8> = VecMap::new();
    for (k,v) in rank.ranks.iter() {
        vk.insert(*v as usize, k as u8);
    }
    let mut chars = Vec::new();
    for _ in 0..k {
        let mut x = BitVec::with_capacity(2);
        x.push(b.pop().unwrap());
        x.push(b.pop().unwrap());
        let x = x.iter().rev().collect::<BitVec>().to_bytes().iter().nth(0).unwrap() >> (8-bits);
        chars.push(*vk.get(x as usize).unwrap());
    }
    return String::from_utf8(chars).unwrap();
}