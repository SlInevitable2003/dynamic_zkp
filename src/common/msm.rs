use ark_ec::{ScalarMul, VariableBaseMSM};

use crate::common::curve::{G1, Scalar};

/// 批量归一化 bases 后做高效 MSM(Pippenger/wNAF)。
pub fn msm(bases: &[G1], scalars: &[Scalar]) -> G1 {
    assert_eq!(bases.len(), scalars.len(), "bases and scalars must have the same length");
    let affine = G1::batch_convert_to_mul_base(bases);
    G1::msm_unchecked(&affine, scalars)
}
