//! 多标量乘法(MSM)助手。

use crate::common::curve::{G1, Scalar};

/// 朴素 MSM:计算 `Σ_i scalars[i] · bases[i]`。
///
/// 学习用清晰实现;后续需要性能时换成 [`ark_ec::msm::VariableBaseMSM`]。
pub fn msm(bases: &[G1], scalars: &[Scalar]) -> G1 {
    assert_eq!(bases.len(), scalars.len(), "bases 与 scalars 长度必须一致");
    bases.iter().zip(scalars).map(|(b, s)| *b * *s).sum()
}
