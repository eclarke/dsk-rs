//! # Kmer counting for arbitrary sizes of k
//! 
//! Adaptation of rust-bio's bio::alphabets::RankTransform::qgrams algorithm;
//! much of the code has been lifted directly from that module and adapted
//! where necessary to provide support for larger `k`.
//! 
//! Example:
//! ```
//! let text = b"ATGCATGCATGC".into_iter();
//! let alphabet = alphabets::dna::alphabet();
//! let counter = KmerCounter::for_large_k(8, &alphabet);
//! for kmer in counter.count(text) {
//!     println!("big: {}", counter.decode(kmer));
//! }
//! ```

use std::marker::PhantomData;

use bio::alphabets::{Alphabet, RankTransform};
use bio::utils::{IntoTextIterator, TextIterator};

use num_bigint::{BigUint, ToBigUint};
use num_traits::{Zero, One};

use vec_map::VecMap;
use bit_vec::BitVec;
use byteorder::{BigEndian, WriteBytesExt};

use itertools::Itertools;

type Small = usize;
type Large = BigUint;



pub struct KmerCounter<KmerSize> {
    k: usize,
    size: PhantomData<*const KmerSize>,
    ranks: RankTransform,
    bits_per_letter: usize,
    bits2char: VecMap<u8>,
}


impl<K> KmerCounter<K> {

    fn take_bits(&self, bitvec: &mut BitVec, nbits: usize) -> BitVec {
        let mut bits = BitVec::with_capacity(nbits);
        for _  in 0..nbits {
            bits.push(bitvec.pop().unwrap_or(false));
        }
        bits.into_iter().rev().collect::<BitVec>()
    }

    pub fn decode(&self, bytes: &[u8]) -> String {
        let mut b = BitVec::from_bytes(bytes);
        let mut chars: Vec<u8> = Vec::new();
        for _ in 0..self.k {
            let chunk = self.take_bits(&mut b, self.bits_per_letter as usize);
            let byte = chunk.to_bytes()[0] as usize >> (8-self.bits_per_letter);
            chars.push(*self.bits2char.get(byte).unwrap())
        }
        chars = chars.into_iter().rev().collect_vec();
        String::from_utf8(chars).unwrap()
    }
}

impl KmerCounter<Large> {
    pub fn for_large_k(k: usize, alphabet: &Alphabet) -> KmerCounter<Large> {
        let ranks = RankTransform::new(alphabet);
        let bits_per_letter = (ranks.ranks.len() as f32).log2().ceil() as usize;
        let mut bits2char = VecMap::new();
        for (k, v) in ranks.ranks.iter() {
            bits2char.insert(*v as usize, k as u8);
        }
        KmerCounter {
            k, size: PhantomData, ranks, bits_per_letter, bits2char 
        }
    }

    pub fn count<'a, T: TextIterator<'a>>(&'a self, text: T) -> Kmers<'a, T, Large> {
        let mut kmers = Kmers::<'a, _, Large>::new(&self, text);
        for _ in 0..self.k - 1 {
            kmers.next();
        }
        kmers
    }
}

impl KmerCounter<Small> {
    pub fn for_small_k(k: usize, alphabet: &Alphabet) -> KmerCounter<Small> {
        let ranks = RankTransform::new(alphabet);
        let bits_per_letter = (ranks.ranks.len() as f32).log2().ceil() as usize;
        let mut bits2char = VecMap::new();
        for (k, v) in ranks.ranks.iter() {
            bits2char.insert(*v as usize, k as u8);
        }
        KmerCounter {
            k, size: PhantomData, ranks, bits_per_letter, bits2char 
        }
    }

    pub fn count<'a, T: TextIterator<'a>>(&'a self, text: T) -> Kmers<'a, T, Small> {
        let mut kmers = Kmers::<'a, _, Small>::new(&self, text);
        for _ in 0..self.k - 1 {
            kmers.next();
        }
        kmers
    }
}

pub struct Kmers<'a, T, K> {
    text: T,
    kmer: K,
    mask: K,
    bits_per_letter: usize,
    ranks: &'a RankTransform,
}

impl<'a, T: TextIterator<'a>> Kmers<'a, T, Large> {
    fn new(counter: &'a KmerCounter<Large>, text: T) -> Self {
        let bits_per_kmer = counter.bits_per_letter * counter.k;
        let mask  = (BigUint::one() << bits_per_kmer) - BigUint::one();
        println!("bits per kmer: {}, mask: {:b}", bits_per_kmer, mask);
        Kmers {
            text: text.into_iter(),
            kmer: BigUint::zero(),
            mask,
            bits_per_letter: counter.bits_per_letter,
            ranks: &counter.ranks,
        }
    }

    fn push(&mut self, a: u8) {
        print!("before: a={}, kmer={:b}; ", a, self.kmer);
        let a: BigUint = a.to_biguint().unwrap();
        self.kmer = &self.kmer << self.bits_per_letter;
        self.kmer = &self.kmer | a;
        self.kmer = &self.kmer & &self.mask;
        println!("after: kmer={:b}", self.kmer);
    }

    fn kmer(&self) -> Vec<u8> {
        self.kmer.to_bytes_be()
    }
}

impl<'a, T: TextIterator<'a>> Kmers<'a, T, Small> {
    fn new(counter: &'a KmerCounter<Small>, text: T) -> Self {
        let bits_per_kmer = counter.bits_per_letter * counter.k;

        Kmers {
            text: text.into_iter(),
            kmer: 0,
            mask: (1 << bits_per_kmer) - 1,
            bits_per_letter: counter.bits_per_letter,
            ranks: &counter.ranks,
        }
    }

    fn push(&mut self, a: u8) {
        self.kmer <<= self.bits_per_letter;
        self.kmer |= a as usize;
        self.kmer &= self.mask;        
    }

    fn kmer(&self) -> Vec<u8> {
        let mut bytes = vec![];
        bytes.write_uint::<BigEndian>(self.kmer as u64, 8).unwrap();
        bytes
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, T, Large> {
    type Item = Vec<u8>;
    
    fn next(&mut self) -> Option<Vec<u8>> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.push(b);
                Some(self.kmer())
            },
            None => None
        }
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, T, Small> {
    type Item = Vec<u8>;
    
    fn next(&mut self) -> Option<Vec<u8>> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.push(b);
                Some(self.kmer())
            },
            None => None
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use bio::alphabets;

    #[test]
    fn test_small_kmer() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_small_k(15, &alphabet);
        let kmers = counter.count(text).map(|kmer| counter.decode(&kmer)).collect_vec();
        assert_eq!(kmers, vec![String::from("AATTCCGGAATTCCG"), String::from("ATTCCGGAATTCCGG"), String::from("TTCCGGAATTCCGGN")]);
    }

    #[test]
    fn test_big_kmer() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_large_k(16, &alphabet);
        let kmers = counter.count(text).map(|kmer| {println!("{:?}", kmer); counter.decode(&kmer)}).collect_vec();
        assert_eq!(kmers, vec![String::from("AATTCCGGAATTCCGG"),String::from("ATTCCGGAATTCCGGN")]);
    }

    #[test]
    fn test_big_kmer2() {
        let text = b"AATTCCGGAATTCCGGN".into_iter();
        let alphabet = alphabets::dna::n_alphabet();
        let counter = KmerCounter::for_large_k(15, &alphabet);
        let kmers = counter.count(text).map(|kmer| counter.decode(&kmer)).collect_vec();
        assert_eq!(kmers, vec![String::from("AATTCCGGAATTCCG"), String::from("ATTCCGGAATTCCGG"), String::from("TTCCGGAATTCCGGN")]);
    }
}