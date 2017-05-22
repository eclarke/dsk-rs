/// Adaptation of rust-bio's kmer index for arbitrary sizes of q (or k, in this case).

use num_bigint::BigUint;
use num_traits::{Zero, One, Num};
use num::FromPrimitive;

use bio::alphabets::{RankTransform};
use bio::utils::{IntoTextIterator, TextIterator};

/// Iterator over q-grams.
pub struct Kmers<'a, T: TextIterator<'a>> {
    text: T,
    ranks: &'a RankTransform,
    bits: u32,
    mask: BigUint,
    kmer: BigUint,
}



impl<'a, T: TextIterator<'a>> Kmers<'a, T> {
    fn kmer_push(&mut self, a: u8) {
        let a: BigUint = FromPrimitive::from_u8(a).unwrap();
        // TODO: are these expensive clones? Maybe use mem::replace?
        let mut kmer = self.kmer.clone() << self.bits as usize;
        kmer = kmer | a;
        self.kmer = kmer & self.mask.clone();
    }
}


impl<'a, T: TextIterator<'a>> Iterator for Kmers<'a, T> {
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

pub fn kmers<'a, T: IntoTextIterator<'a>>(ranks: &'a RankTransform, q: u32, text: T) -> Kmers<'a, T::IntoIter> {
    let bits = (ranks.ranks.len() as f32).log2().ceil() as u32;
    // assert!((bits * q) as usize <= mem::size_of::<usize>() * 8,
            // "Expecting q to be smaller than usize / log2(|A|)");

    let mut kmers = Kmers {
        text: text.into_iter(),
        ranks: ranks,
        bits: bits,
        mask: (BigUint::one() << ((q * bits) as usize)) - BigUint::one(),
        kmer: BigUint::zero(),
    };

    for _ in 0..q - 1 {
        kmers.next();
    }

    kmers
}