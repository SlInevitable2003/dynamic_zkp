//! 配对与生成元乘法助手。

use crate::common::curve::{Engine, G1, G2, GT, Scalar};
use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;

/// `[x]_1 = x · g1`。
pub fn mul_g1(x: Scalar) -> G1 {
    G1::generator() * x
}

/// `[x]_2 = x · g2`。
pub fn mul_g2(x: Scalar) -> G2 {
    G2::generator() * x
}

/// 配对 `e(p, q) ∈ GT`。
pub fn pair(p: &G1, q: &G2) -> GT {
    Engine::pairing(*p, *q)
}
