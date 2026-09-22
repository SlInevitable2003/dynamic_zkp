//! Dynamo 系协议共用的多项式 / 域原语:单位根、Lagrange 求值、置换、`y_i`、`u`。
//!
//! 记号见 README §0。

use crate::common::curve::Scalar;
use ark_ff::{FftField, Field};

/// 断言 `m` 是 2 的幂(返回其指数,便于复用)。
fn log2_of_power_of_two(m: usize) -> u32 {
    assert!(m.is_power_of_two(), "m must be a power of 2, got {m}");
    m.trailing_zeros()
}

/// `m` 次本原单位根 `ω`(`m` 必须为 2 的幂)。
pub fn root_of_unity(m: usize) -> Scalar {
    let _ = log2_of_power_of_two(m);
    Scalar::get_root_of_unity(m as u64)
        .expect("m does not divide the order of the multiplicative group of F")
}

/// 计算所有 `L_i(τ)`(`i ∈ [0, m)`),用闭式
///
/// ```text
/// L_i(τ) = (τ^m - 1) / (m · (ω^{-i} · τ - 1))
/// ```
///
/// 一次算满全部 `i`,O(m)。由此 `L_i(0) = 1/m` 对任意 `i` 成立。
pub fn lagrange_evals(tau: Scalar, m: usize) -> Vec<Scalar> {
    let w = root_of_unity(m);
    let num = tau.pow([m as u64]) - Scalar::ONE;
    let inv_m = Scalar::from(m as u64).inverse().expect("m is nonzero");
    let w_inv = w.inverse().expect("root of unity is nonzero");

    let mut w_pow = Scalar::ONE; // 跟踪 ω^{-i},从 i = 0 开始
    (0..m)
        .map(|_| {
            let denom = w_pow * tau - Scalar::ONE;
            let li = num * denom.inverse().expect("τ != ω^i (分母非零)") * inv_m;
            w_pow *= w_inv;
            li
        })
        .collect()
}

/// 置换 `σ` 的逆 `σ^{-1}`,满足 `sigma_inv[sigma[i]] == i`。
pub fn inverse_permutation(sigma: &[usize]) -> Vec<usize> {
    let m = sigma.len();
    let mut inv = vec![usize::MAX; m];
    for (i, &s) in sigma.iter().enumerate() {
        assert!(s < m, "sigma[{i}] = {s} 越界,不是 (0..m) 上的置换");
        assert_eq!(inv[s], usize::MAX, "sigma 不是置换:元素 {s} 重复");
        inv[s] = i;
    }
    inv
}

/// `y_i(τ_Y) = L_i(τ_Y) - L_{σ^{-1}(i)}(τ_Y)`。
pub fn y_evals(ly: &[Scalar], sigma_inv: &[usize]) -> Vec<Scalar> {
    debug_assert_eq!(ly.len(), sigma_inv.len());
    (0..ly.len()).map(|i| ly[i] - ly[sigma_inv[i]]).collect()
}

/// `u(τ_X, τ_Y) = Σ_{i} L_i(τ_X) · y_i(τ_Y)`。
pub fn u_eval(lx: &[Scalar], y: &[Scalar]) -> Scalar {
    debug_assert_eq!(lx.len(), y.len());
    lx.iter().zip(y).map(|(a, b)| *a * *b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, Field};
    use ark_std::UniformRand;

    /// 朴素几何级数 `L_i(τ) = (1/m) Σ_{j=0}^{m-1} (ω^{-i} τ)^j`,用于交叉验证闭式。
    fn lagrange_naive(tau: Scalar, m: usize) -> Vec<Scalar> {
        let w = root_of_unity(m);
        let inv_m = Scalar::from(m as u64).inverse().unwrap();
        let w_inv = w.inverse().unwrap();

        let mut out = Vec::with_capacity(m);
        let mut w_pow = Scalar::ONE; // ω^{-i}
        for _ in 0..m {
            let t = w_pow * tau;
            let mut acc = Scalar::ZERO;
            let mut t_pow = Scalar::ONE;
            for _ in 0..m {
                acc += t_pow;
                t_pow *= t;
            }
            out.push(acc * inv_m);
            w_pow *= w_inv;
        }
        out
    }

    #[test]
    fn lagrange_closed_form_matches_naive() {
        let mut rng = ark_std::test_rng();
        for m in [4usize, 8, 16] {
            let tau = Scalar::rand(&mut rng);
            assert_eq!(lagrange_evals(tau, m), lagrange_naive(tau, m), "m = {m}");
        }
    }

    #[test]
    fn lagrange_partition_of_unity() {
        let mut rng = ark_std::test_rng();
        for m in [4usize, 8, 16] {
            let tau = Scalar::rand(&mut rng);
            let sum: Scalar = lagrange_evals(tau, m).iter().sum();
            assert_eq!(sum, Scalar::ONE, "m = {m}");
        }
    }

    #[test]
    fn inverse_permutation_roundtrip() {
        let sigma = [2, 0, 3, 1];
        let inv = inverse_permutation(&sigma);
        assert_eq!(inv, [1, 3, 0, 2]);
        for i in 0..sigma.len() {
            assert_eq!(inv[sigma[i]], i);
        }
    }
}
