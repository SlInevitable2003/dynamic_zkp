//! Dynamo(非通用版)的公开参数与证明结构体。

use crate::common::curve::{G1, G2};

/// 证明密钥(README §1.1.1 的 `pk`)。
#[derive(Clone, Debug)]
pub struct ProvingKey {
    pub alpha: Vec<G1>, // pk.α[i]
    pub beta: Vec<G1>,  // pk.β[i]
    pub z: Vec<G1>,     // pk.Z[i]
    pub h: Vec<G1>,     // pk.H[i]
    pub b: Vec<G1>,     // pk.B[i]
    pub lx: Vec<G1>,    // pk.LX[i]
    pub ly: Vec<G1>,    // pk.LY[i]

    pub u: G1,     // [u(τ_X, τ_Y)]_1
    pub x: G1,     // [τ_X]_1
    pub x2: G1,    // [τ_X^2]_1
    pub xm: G1,    // [τ_X^m]_1
    pub ym: G1,    // [τ_Y^m]_1
    pub xm_ym: G1, // [τ_X^m · τ_Y^m]_1
    pub x2_ym: G1, // [τ_X^2 · τ_Y^m]_1
}

impl ProvingKey {
    /// 预分配容量并置标量为单位元,由 `setup` 填充。
    pub(crate) fn with_capacity(m: usize) -> Self {
        Self {
            alpha: Vec::with_capacity(m),
            beta: Vec::with_capacity(m),
            z: Vec::with_capacity(m),
            h: Vec::with_capacity(m),
            b: Vec::with_capacity(m),
            lx: Vec::with_capacity(m),
            ly: Vec::with_capacity(m),
            u: G1::default(),
            x: G1::default(),
            x2: G1::default(),
            xm: G1::default(),
            ym: G1::default(),
            xm_ym: G1::default(),
            x2_ym: G1::default(),
        }
    }
}

/// 验证密钥(README §1.1.1 的 `vk`)。
#[derive(Clone, Debug)]
pub struct VerifyingKey {
    pub u: G2,    // [u(τ_X, τ_Y)]_2
    pub x: G2,    // [τ_X]_2
    pub x2: G2,   // [τ_X^2]_2
    pub xm: G2,   // [τ_X^m]_2
    pub ym: G2,   // [τ_Y^m]_2
    pub invm: G2, // [m^{-1}]_2
}

/// 公开实例 `x = ([z_zk]_1, [h_zk]_1)`。
#[derive(Clone, Copy, Debug)]
pub struct Instance {
    pub z: G1,
    pub h: G1,
}

/// 证明 `π`。
#[derive(Clone, Copy, Debug)]
pub struct Proof {
    pub alpha: G1, // [α_zk]_1
    pub beta: G1,  // [β_zk]_1
    pub q: G1,     // [q_zk]_1
    pub z: G1,     // [Z_zk]_1
    pub h: G1,     // [H_zk]_1
    pub b: G1,     // [B_zk]_1
}
