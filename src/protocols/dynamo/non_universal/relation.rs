//! Dynamo(非通用版)证明的见证关系。
//!
//! 由协议恒等式(README §1.1.3 第 1 条配对校验)可推出:证明通过当且仅当见证满足
//! `Σ_i z_i · y_i(τ_Y) = Σ_i h_i · L_i(τ_Y)`,即
//!
//! ```text
//! h_i = z_i - z_{σ(i)}
//! ```
//!
//! (这里 `σ` 作用于 `z` 的指标:`z_{σ(i)}`。)这刻画了 Dynamo 的“动态”关系。

use crate::common::curve::Scalar;

/// 由 `z` 与置换 `σ` 派生满足关系的 `h`:`h_i = z_i - z_{σ(i)}`。
pub fn derive_h(z: &[Scalar], sigma: &[usize]) -> Vec<Scalar> {
    assert_eq!(z.len(), sigma.len(), "z 与 σ 长度必须一致");
    (0..z.len()).map(|i| z[i] - z[sigma[i]]).collect()
}
