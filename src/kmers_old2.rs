/// Adaptation of rust-bio's kmer index for arbitrary sizes of q (or k, in this case).

use std::mem;

use num_bigint::BigUint;
use num_traits::{Zero, One};
use num::FromPrimitive;

use bio::alphabets::RankTransform;
use bio::utils::TextIterator;

pub struct Kmers<'a, K, T: TextIterator<'a>> {
    text: T,
    ranks: &'a RankTransform,
    bits: u32,
    mask: K,
    kmer: K,
}

pub trait Counter<'a, K, T> {
    fn count_kmers(text: T, k: u32, rank: &'a RankTransform) -> Kmers<K, T> where T:TextIterator<'a>;
}

impl<'a, T: TextIterator<'a>> Kmers<'a, BigUint, T> {
    fn kmer_push(&mut self, a: u8) {
        let a: BigUint = FromPrimitive::from_u8(a).unwrap();
        // TODO: are these expensive clones? Maybe use mem::replace?
        let mut kmer = self.kmer.clone() << self.bits as usize;
        kmer = kmer | a;
        self.kmer = kmer & self.mask.clone();
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, BigUint, T> {
    type Item = BigUint;

    fn next(&mut self) -> Option<BigUint> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.kmer_push(b);
                Some(self.kmer.clone())
            }
            None => None,
        }
    }
}

impl<'a, T: TextIterator<'a>> Counter<'a, BigUint, T> for Kmers<'a, BigUint, T> {
    fn count_kmers(text: T, k: u32, rank: &'a RankTransform) -> Self {
        let bits = (rank.ranks.len() as f32).log2().ceil() as u32;
        // assert!((bits * q) as usize <= mem::size_of::<usize>() * 8,
                // "Expecting q to be smaller than usize / log2(|A|)");

        let mut kmers = Kmers {
            text: text.into_iter(),
            ranks: rank,
            bits: bits,
            mask: (BigUint::one() << ((k * bits) as usize)) - BigUint::one(),
            kmer: BigUint::zero(),
        };

        for _ in 0..k - 1 {
            kmers.next();
        }

        kmers
        }
}

impl<'a, T: TextIterator<'a>> Kmers<'a, usize, T> {
    fn kmer_push(&mut self, a: u8) {
        self.kmer <<= self.bits;
        self.kmer |= a as usize;
        self.kmer &= self.mask;
    }
}

impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, usize, T> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a);
                self.kmer_push(b);
                Some(self.kmer)
            }
            None => None,
        }
    }
}

impl<'a, T: TextIterator<'a>> Counter<'a, usize, T> for Kmers<'a, usize, T> {
    fn count_kmers(text: T, k: u32, rank: &'a RankTransform) -> Kmers<'a, usize, T> {
        let bits = (rank.ranks.len() as f32).log2().ceil() as u32;
        assert!((bits * k as u32) as usize <= mem::size_of::<usize>() * 8,
                "Expecting q to be smaller than usize / log2(|A|)");

        let mut kmers = Kmers {
            text: text.into_iter(),
            ranks: rank,
            bits: bits,
            mask: (1 << (k * bits)) - 1,
            kmer: 0,
        };

        for _ in 0..k - 1 {
            kmers.next();
        }

        kmers
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use bio::alphabets;

//     #[test]
//     fn test_small_kmer() {
//         let text = b"ATGCATGCATGC".into_iter();
//         let rank = RankTransform::new(&alphabets::dna::alphabet());
//         let kmer: Kmers<usize, _> = Kmers::count_kmers(text, 8, &rank);
//         for k in kmer {
//             println!("small: {}", k);
//         }
//     }

//     #[test]
//     fn test_big_kmer() {
//         let text = b"ATGCATGCATGC".into_iter();
//         let rank = RankTransform::new(&alphabets::dna::alphabet());
//         let kmer: Kmers<BigUint, _> = Kmers::count_kmers(text, 8, &rank);
//         for k in kmer {
//             println!("big: {}", k);
//         }
//     }
// }