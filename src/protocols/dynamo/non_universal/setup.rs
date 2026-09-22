//! Dynamo(非通用版)Setup(README §1.1.1)。

use crate::common::curve::{G1, G2, Scalar};
use crate::common::poly::{inverse_permutation, lagrange_evals, u_eval, y_evals};
use crate::protocols::dynamo::non_universal::types::{ProvingKey, VerifyingKey};
use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_std::rand::Rng;
use ark_std::UniformRand;

/// 生成 `(pk, vk)`。`sigma` 是 `(0..m)` 上的置换,`m` 为 2 的幂。
pub fn setup<R: Rng>(m: usize, sigma: &[usize], rng: &mut R) -> (ProvingKey, VerifyingKey) {
    assert!(m.is_power_of_two(), "m must be a power of 2");
    assert_eq!(sigma.len(), m, "sigma must have length m");

    let tau_x = Scalar::rand(rng);
    let tau_y = Scalar::rand(rng);

    let lx = lagrange_evals(tau_x, m); // L_i(τ_X)
    let ly = lagrange_evals(tau_y, m); // L_i(τ_Y)
    let sigma_inv = inverse_permutation(sigma);
    let y = y_evals(&ly, &sigma_inv); // y_i(τ_Y)
    let u = u_eval(&lx, &y); // u(τ_X, τ_Y)

    let tau_x_m = tau_x.pow([m as u64]); // τ_X^m
    let tau_y_m = tau_y.pow([m as u64]); // τ_Y^m
    let inv_m = Scalar::from(m as u64).inverse().expect("m is nonzero");
    let inv_denom = (tau_x_m - Scalar::ONE)
        .inverse()
        .expect("τ_X^m != 1 (否则除以零)");
    let inv_tau_x = tau_x.inverse().expect("τ_X nonzero");

    let g1 = G1::generator();
    let g2 = G2::generator();

    let mut pk = ProvingKey::with_capacity(m);
    for i in 0..m {
        let li_x = lx[i];
        let yi_y = y[i];
        pk.alpha.push(g1 * (li_x * (u - yi_y) * inv_denom));
        pk.beta.push(g1 * (yi_y * (li_x - inv_m) * inv_tau_x));
        pk.z.push(g1 * (tau_y_m * li_x));
        pk.h.push(g1 * (tau_x_m * yi_y));
        pk.b.push(g1 * (yi_y * (li_x - inv_m) * tau_x));
        pk.lx.push(g1 * li_x);
        pk.ly.push(g1 * ly[i]);
    }

    pk.u = g1 * u;
    pk.x = g1 * tau_x;
    pk.x2 = g1 * (tau_x * tau_x);
    pk.xm = g1 * tau_x_m;
    pk.ym = g1 * tau_y_m;
    pk.xm_ym = g1 * (tau_x_m * tau_y_m);
    pk.x2_ym = g1 * (tau_x * tau_x * tau_y_m);

    let vk = VerifyingKey {
        u: g2 * u,
        x: g2 * tau_x,
        x2: g2 * (tau_x * tau_x),
        xm: g2 * tau_x_m,
        ym: g2 * tau_y_m,
        invm: g2 * inv_m,
    };

    (pk, vk)
}
