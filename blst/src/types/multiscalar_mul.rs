
//! Multiscalar multiplication implementation using pippenger algorithm.
use blst::{blst_fr, blst_p1, blst_p1_add, blst_p1_double, blst_p1_in_g1};
use kzg::{Fr, G1};
use crate::types::fr::{FsFr as Scalar, FsFr};
use crate::types::g1::{FsG1 as G1Projective, FsG1};
/// Performs multiscalar multiplication reliying on Pippenger's algorithm.
/// This method was taken from `curve25519-dalek` and was originally made by
/// Oleg Andreev <oleganza@gmail.com>.

#[allow(clippy::needless_collect)]
pub fn msm_variable_base(points: &[G1Projective], scalars: &[Scalar]) -> G1Projective {
    #[cfg(feature = "parallel")]
    use rayon::prelude::*;

    let c = if scalars.len() < 32 {
        3
    } else {
        ln_without_floats(scalars.len()) + 2
    };

    let num_bits = 255usize;
    let fr_one = Scalar::one();

    let zero = G1Projective::identity();
    let window_starts: Vec<_> = (0..num_bits).step_by(c).collect();

    #[cfg(feature = "parallel")]
        let window_starts_iter = window_starts.into_par_iter();
    #[cfg(not(feature = "parallel"))]
        let window_starts_iter = window_starts.into_iter();

    // Each window is of size `c`.
    // We divide up the bits 0..num_bits into windows of size `c`, and
    // in parallel process each such window.
    let window_sums: Vec<_> = window_starts_iter
        .map(|w_start| {
            let mut res = zero;
            // We don't need the "zero" bucket, so we only have 2^c - 1 buckets
            let mut buckets = vec![zero; (1 << c) - 1];
            scalars
                .iter()
                .zip(points)
                .filter(|(s, _)| !(*s == &Scalar::zero()))
                .for_each(|(&scalar, base)| {
                    if scalar == fr_one {
                        // We only process unit scalars once in the first window.
                        if w_start == 0 {
                            res = res.add(base);
                        }
                    } else {

                        // We right-shift by w_start, thus getting rid of the
                        // lower bits.
                        let mut scalar = montgomery_reduce(
                            scalar.0.l[0],
                            scalar.0.l[1],
                            scalar.0.l[2],
                            scalar.0.l[3],
                            0,
                            0,
                            0,
                            0,
                        );

                        let scalar = divn(scalar, w_start as u32);

                        // We mod the remaining bits by the window size.
                        let inside = scalar.0;
                        let scalar = inside.l[0] % (1 << c);

                        // If the scalar is non-zero, we update the corresponding
                        // bucket.
                        // (Recall that `buckets` doesn't have a zero bucket.)
                        if scalar != 0 {
                            let space = (scalar - 1) as usize;
                            buckets[(scalar - 1) as usize] =
                                buckets[(scalar - 1) as usize].add(base);
                        }
                    }
                });
            let mut i = 0;
            let mut running_sum = G1Projective::identity();
            for b in buckets.into_iter(){
                i = i+1;
                let respo = running_sum.add(&b);
                let res = res.add(&respo);
                if( i == 429){
                    let running_sume = running_sum;
                    let rese = res;
                }
            }

            res
        })
        .collect();

    // We store the sum for the lowest window.
    let lowest = *window_sums.first().unwrap();
    // We're traversing windows from high to low.
    let other = window_sums[1..]
        .iter()
        .rev()
        .fold(zero, |mut total, sum_i| unsafe {
            total = total.add(sum_i);
            for _ in 0..c {
                let mut total_p1 = blst_p1::default();
                let raw_pointer = &mut total_p1 as *mut blst_p1;
                let total_p1_ptr = &mut total.0 as *mut blst_p1;

               blst_p1_double(raw_pointer, total_p1_ptr);
                total.0 = *raw_pointer;
                let ok = total.0;
                let ww= *raw_pointer;
                let ss = *total_p1_ptr;
            }
            total
        });
    let answer = other.add(&lowest);
    return answer;
}

fn ln_without_floats(a: usize) -> usize {
    // log2(a) * ln(2)
    (log2(a) * 69 / 100) as usize
}
fn log2(x: usize) -> u32 {
    if x <= 1 {
        return 0;
    }

    let n = x.leading_zeros();
    core::mem::size_of::<usize>() as u32 * 8 - n
}

#[inline]
pub fn divn(mut scalar: Scalar, mut n: u32) -> Scalar {
    if n >= 256 {
        return Scalar::zero();
    }

    while n >= 64 {
        let mut t = 0;
        for i in scalar.0.l.iter_mut().rev() {
            core::mem::swap(&mut t, i);
        }
        n -= 64;
    }

    if n > 0 {
        let mut t = 0;
        for i in scalar.0.l.iter_mut().rev() {
            let t2 = *i << (64 - n);
            *i >>= n;
            *i |= t;
            t = t2;
        }
    }

    scalar
}

pub const MODULUS: Scalar = Scalar(blst_fr{
    l: [
        0xffff_ffff_0000_0001,
        0x53bd_a402_fffe_5bfe,
        0x3339_d808_09a1_d805,
        0x73ed_a753_299d_7d48,
    ],
});
const INV: u64 = 0xffff_fffe_ffff_ffff;

#[inline(always)]
pub const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}
#[inline(always)]
pub const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}
#[inline(always)]
pub fn montgomery_reduce(
    r0: u64,
    r1: u64,
    r2: u64,
    r3: u64,
    r4: u64,
    r5: u64,
    r6: u64,
    r7: u64,
) -> Scalar {
    // The Montgomery reduction here is based on Algorithm 14.32 in
    // Handbook of Applied Cryptography
    // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

    let k = r0.wrapping_mul(INV);
    let (_, carry) = mac(r0, k, MODULUS.0.l[0], 0);
    let (r1, carry) = mac(r1, k, MODULUS.0.l[1], carry);
    let (r2, carry) = mac(r2, k, MODULUS.0.l[2], carry);
    let (r3, carry) = mac(r3, k, MODULUS.0.l[3], carry);
    let (r4, carry2) = adc(r4, 0, carry);

    let k = r1.wrapping_mul(INV);
    let (_, carry) = mac(r1, k, MODULUS.0.l[0], 0);
    let (r2, carry) = mac(r2, k, MODULUS.0.l[1], carry);
    let (r3, carry) = mac(r3, k, MODULUS.0.l[2], carry);
    let (r4, carry) = mac(r4, k, MODULUS.0.l[3], carry);
    let (r5, carry2) = adc(r5, carry2, carry);

    let k = r2.wrapping_mul(INV);
    let (_, carry) = mac(r2, k, MODULUS.0.l[0], 0);
    let (r3, carry) = mac(r3, k, MODULUS.0.l[1], carry);
    let (r4, carry) = mac(r4, k, MODULUS.0.l[2], carry);
    let (r5, carry) = mac(r5, k, MODULUS.0.l[3], carry);
    let (r6, carry2) = adc(r6, carry2, carry);

    let k = r3.wrapping_mul(INV);
    let (_, carry) = mac(r3, k, MODULUS.0.l[0], 0);
    let (r4, carry) = mac(r4, k, MODULUS.0.l[1], carry);
    let (r5, carry) = mac(r5, k, MODULUS.0.l[2], carry);
    let (r6, carry) = mac(r6, k, MODULUS.0.l[3], carry);
    let (r7, _) = adc(r7, carry2, carry);

    let blst_fr = blst_fr{
        l: [
            r4,
            r5,
            r6,
            r7,
        ],
    };

    let mut scalar = Scalar(blst_fr);
    let modulus = MODULUS;
    scalar.sub(&modulus);
    scalar
}

