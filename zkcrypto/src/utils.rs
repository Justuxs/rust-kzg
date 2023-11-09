use super::P1;
use crate::P2;
use bls12_381::{ G1Affine, G1Projective, G2Affine, G2Projective, Scalar};
use blst::{blst_bendian_from_scalar, BLST_ERROR, blst_fp2, blst_fr, blst_fr_from_scalar, blst_p1, blst_p1_affine, blst_p1_compress, blst_p1_from_affine, blst_p1_uncompress, blst_p2, blst_p2_affine, blst_p2_compress, blst_p2_from_affine, blst_p2_uncompress, blst_scalar, blst_scalar_fr_check, blst_scalar_from_bendian, blst_scalar_from_fr};
use kzg::eip_4844::{BYTES_PER_FIELD_ELEMENT, BYTES_PER_G1, BYTES_PER_G2};
use kzg::{G1, G2};
use crate::kzg_types::{ZG1, ZG2};

#[derive(Debug, PartialEq, Eq)]
pub struct Error;

fn form_blst_fr_to_bytes(blstfr : blst_fr) -> [u8; 32] {
    let mut scalar = blst_scalar::default();
    let mut bytes = [0u8; 32];
    unsafe {
        blst_scalar_from_fr(&mut scalar, &blstfr);
        blst_bendian_from_scalar(bytes.as_mut_ptr(), &scalar);
    }

    bytes
}

pub const fn blst_fr_into_pc_fr(fr: blst_fr) -> Scalar {
    let bytes = form_blst_fr_to_bytes(fr);
    Scalar::from_bytes(&bytes).unwrap()
}

fn from_bytes_to_blst_fr(bytes: &[u8]) -> Result<blst_fr, String> {
    bytes
        .try_into()
        .map_err(|_| {
            format!(
                "Invalid byte length. Expected {}, got {}",
                BYTES_PER_FIELD_ELEMENT,
                bytes.len()
            )
        })
        .and_then(|bytes: &[u8; BYTES_PER_FIELD_ELEMENT]| {
            let mut bls_scalar = blst_scalar::default();
            let mut fr = blst_fr::default();
            unsafe {
                blst_scalar_from_bendian(&mut bls_scalar, bytes.as_ptr());
                blst_fr_from_scalar(&mut fr, &bls_scalar);
            }
            Ok(fr)
        })
}
pub const fn pc_fr_into_blst_fr(scalar: Scalar) -> blst_fr {
    let bytes = scalar.to_bytes();
    from_bytes_to_blst_fr(bytes.as_ref()).unwrap()
}

//G1
fn blst_p1_to_bytes(p1: &P1) -> [u8; 48] {
    let mut out = [0u8; 48];
    unsafe {
        blst_p1_compress(out.as_mut_ptr(), p1);
    }
    out
}

fn from_bytes_to_blst_p1(bytes: &[u8]) -> Result<blst_p1, String> {
    bytes
        .try_into()
        .map_err(|_| {
            format!(
                "Invalid byte length. Expected {}, got {}",
                BYTES_PER_G1,
                bytes.len()
            )
        })
        .and_then(|bytes: &[u8; BYTES_PER_G1]| {
            let mut tmp = blst_p1_affine::default();
            let mut g1 = blst_p1::default();
            unsafe {
                blst_p1_from_affine(&mut g1, &tmp);
            }
            Ok(g1)
        })
}
pub const fn blst_p1_into_pc_g1projective(p1: &P1) -> G1Projective {
    let bytes = blst_p1_to_bytes(p1);
    ZG1::from_bytes(bytes.as_ref()).unwrap().proj
}
pub const fn pc_g1projective_into_blst_p1(p1: G1Projective) -> blst_p1 {
    let g1_affine = G1Affine::from(p1);
    let bytes = g1_affine.to_compressed();
    from_bytes_to_blst_p1(bytes.as_ref()).unwrap()
}

//G2

fn from_bytes_to_blst_p2(bytes: &[u8]) -> Result<blst_p2, String> {
    bytes
        .try_into()
        .map_err(|_| {
            format!(
                "Invalid byte length. Expected {}, got {}",
                BYTES_PER_G2,
                bytes.len()
            )
        })
        .and_then(|bytes: &[u8; BYTES_PER_G2]| {
            let mut tmp = blst_p2_affine::default();
            let mut g2 = blst_p2::default();
            unsafe {
                blst_p2_from_affine(&mut g2, &tmp);
            }
            Ok(g2)
        })
}

fn from_blst_p2_to_bytes(p2: &P2) -> &[u8] {
    let mut out = [0u8; BYTES_PER_G2];
    unsafe {
        blst_p2_compress(out.as_mut_ptr(), p2);
    }
    &out[..]
}

pub const fn blst_p2_into_pc_g2projective(p2: &P2) -> G2Projective {
    let bytes = from_blst_p2_to_bytes(p2);
    ZG2::from_bytes(bytes).unwrap().proj
}
pub const fn pc_g2projective_into_blst_p2(p2: G2Projective) -> blst_p2 {
    let g2_affine = G2Affine::from(p2);
    let bytes = g2_affine.to_compressed();
    from_bytes_to_blst_p2(bytes.as_ref()).unwrap()
}
