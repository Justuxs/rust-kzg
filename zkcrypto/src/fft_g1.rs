use crate::consts::G1_GENERATOR;
use crate::kzg_proofs::FFTSettings;
use crate::kzg_types::{FsFp, FsG1Affine, FsG1ProjAddAffine, ZFr, ZG1};
use crate::multiscalar_mul::msm_variable_base;
use kzg::{Fr as KzgFr, G1Affine, G1Mul, Scalar256};
use kzg::{FFTG1, G1};
use std::ops::MulAssign;
use bls12_381 ::Scalar;
use kzg::msm::tilling_parallel_pippinger::{parallel_affine_conv, tiling_parallel_pippinger};

/*#[warn(unused_variables)]
pub fn g1_linear_combination(out: &mut ZG1, points: &[ZG1], scalars: &[ZFr], _len: usize) {
    let g1 = msm_variable_base(points, scalars);
    out.proj = g1
}*/

pub fn g1_linear_combination(out: &mut ZG1, points: &[ZG1], scalars: &[ZFr], len: usize) {
    if len < 8 {
        *out = ZG1::default();
        for i in 0..len {
            let tmp = points[i].mul(&scalars[i]);
            *out = out.add_or_dbl(&tmp);
        }
        return;
    }

    #[cfg(feature = "parallel")]
    {
        // Atleast on my machine - performance was *slightly worse* with the parallel version
        // let points = FsG1Affine::into_affines(points);
        let points = parallel_affine_conv::<ZG1, FsFp, FsG1Affine>(points);

        let scalars = scalars
            .iter()
            .map(|b| {
                let mut scalar = b.fr;
                Scalar256::from_u64(scalar.0)
            })
            .collect::<Vec<_>>();
        *out = tiling_parallel_pippinger(&points, scalars.as_slice());
    }

    #[cfg(not(feature = "parallel"))]
    {
        let ark_points = FsG1Affine::into_affines(points);
        let ark_scalars = {
            scalars
                .iter()
                .take(len)
                .map(|scalar| {
                    let bytes: [u8; 32] = scalar.fr.into();
                    Scalar256::from_u8(&bytes) //klaida
                })
                .collect::<Vec<_>>()
        };

        *out = kzg::msm::arkmsm_msm::VariableBaseMSM::multi_scalar_mul::<
            ZG1,
            FsFp,
            FsG1Affine,
            FsG1ProjAddAffine,
        >(&ark_points, &ark_scalars)
    }
}
pub fn make_data(data: usize) -> Vec<ZG1> {
    let mut vec = Vec::new();
    if data != 0 {
        vec.push(G1_GENERATOR);
        for i in 1..data as u64 {
            let res = vec[(i - 1) as usize].add_or_dbl(&G1_GENERATOR);
            vec.push(res);
        }
    }
    vec
}

impl FFTG1<ZG1> for FFTSettings {
    fn fft_g1(&self, data: &[ZG1], inverse: bool) -> Result<Vec<ZG1>, String> {
        if data.len() > self.max_width {
            return Err(String::from("data length is longer than allowed max width"));
        }
        if !data.len().is_power_of_two() {
            return Err(String::from("data length is not power of 2"));
        }

        let stride: usize = self.max_width / data.len();
        let mut ret = vec![ZG1::default(); data.len()];

        let roots = if inverse {
            &self.reverse_roots_of_unity
        } else {
            &self.expanded_roots_of_unity
        };

        fft_g1_fast(&mut ret, data, 1, roots, stride, 1);

        if inverse {
            let inv_fr_len = ZFr::from_u64(data.len() as u64).inverse();
            ret[..data.len()]
                .iter_mut()
                .for_each(|f| f.proj.mul_assign(&inv_fr_len.fr));
        }
        Ok(ret)
    }
}

pub fn fft_g1_slow(
    ret: &mut [ZG1],
    data: &[ZG1],
    stride: usize,
    roots: &[ZFr],
    roots_stride: usize,
    _width: usize,
) {
    for i in 0..data.len() {
        ret[i] = data[0].mul(&roots[0]);
        for j in 1..data.len() {
            let jv = data[j * stride];
            let r = roots[((i * j) % data.len()) * roots_stride];
            let v = jv.mul(&r);
            ret[i] = ret[i].add_or_dbl(&v);
        }
    }
}

pub fn fft_g1_fast(
    ret: &mut [ZG1],
    data: &[ZG1],
    stride: usize,
    roots: &[ZFr],
    roots_stride: usize,
    _width: usize,
) {
    let half = ret.len() / 2;
    if half > 0 {
        #[cfg(feature = "parallel")]
        {
            let (lo, hi) = ret.split_at_mut(half);
            rayon::join(
                || fft_g1_fast(hi, &data[stride..], stride * 2, roots, roots_stride * 2, 1),
                || fft_g1_fast(lo, data, stride * 2, roots, roots_stride * 2, 1),
            );
        }

        #[cfg(not(feature = "parallel"))]
        {
            fft_g1_fast(
                &mut ret[..half],
                data,
                stride * 2,
                roots,
                roots_stride * 2,
                1,
            );
            fft_g1_fast(
                &mut ret[half..],
                &data[stride..],
                stride * 2,
                roots,
                roots_stride * 2,
                1,
            );
        }

        for i in 0..half {
            let y_times_root = ret[i + half].mul(&roots[i * roots_stride]);
            ret[i + half] = ret[i].sub(&y_times_root);
            ret[i] = ret[i].add_or_dbl(&y_times_root);
        }
    } else {
        ret[0] = data[0];
    }
}
