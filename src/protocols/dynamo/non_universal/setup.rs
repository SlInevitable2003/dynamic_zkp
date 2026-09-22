use crate::common::curve::{G1, G2, Scalar};
use crate::common::poly::{inverse_permutation, lagrange_evals, u_eval, y_evals};
use crate::protocols::dynamo::non_universal::types::{ProvingKey, VerifyingKey};

use ark_ec::scalar_mul::BatchMulPreprocessing;
use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use rayon::prelude::*;
use std::time::Instant;

pub fn setup<R: Rng>(m: usize, sigma: &[usize], rng: &mut R) -> (ProvingKey, VerifyingKey) {
    let setup_start = Instant::now();
    assert!(m.is_power_of_two(), "m must be a power of 2");
    assert_eq!(sigma.len(), m, "sigma must have length m");

    let tau_x = Scalar::rand(rng);
    let tau_y = Scalar::rand(rng);

    let t = Instant::now();
    let lx = lagrange_evals(tau_x, m); // L_i(τ_X)
    let ly = lagrange_evals(tau_y, m); // L_i(τ_Y)
    let sigma_inv = inverse_permutation(sigma);
    let y = y_evals(&ly, &sigma_inv); // y_i(τ_Y)
    let u = u_eval(&lx, &y); // u(τ_X, τ_Y)
    eprintln!("[setup] poly evals: {:?}", t.elapsed());

    let t = Instant::now();
    let tau_x_m = tau_x.pow([m as u64]); // τ_X^m
    let tau_y_m = tau_y.pow([m as u64]); // τ_Y^m
    let inv_m = Scalar::from(m as u64).inverse().expect("m is nonzero");
    let inv_denom = (tau_x_m - Scalar::ONE).inverse().expect("τ_X^m != 1");
    let inv_tau_x = tau_x.inverse().expect("τ_X nonzero");
    eprintln!("[setup] scalar constants: {:?}", t.elapsed());

    let g1 = G1::generator();
    let g2 = G2::generator();

    let t = Instant::now();
    let g1_table = BatchMulPreprocessing::new(g1, m);
    eprintln!("[setup] fixed-base table: {:?}", t.elapsed());

    let t = Instant::now();
    let alpha_scalars: Vec<Scalar> = (0..m).into_par_iter().map(|i| lx[i] * (u - y[i]) * inv_denom).collect();
    let beta_scalars: Vec<Scalar> = (0..m).into_par_iter().map(|i| y[i] * (lx[i] - inv_m) * inv_tau_x).collect();
    let z_scalars: Vec<Scalar> = (0..m).into_par_iter().map(|i| tau_y_m * lx[i]).collect();
    let h_scalars: Vec<Scalar> = (0..m).into_par_iter().map(|i| tau_x_m * y[i]).collect();
    let b_scalars: Vec<Scalar> = (0..m).into_par_iter().map(|i| y[i] * (lx[i] - inv_m) * tau_x).collect();
    eprintln!("[setup] scalar columns: {:?}", t.elapsed());

    let t = Instant::now();
    let alpha: Vec<G1> = g1_table.batch_mul(&alpha_scalars).into_iter().map(G1::from).collect();
    let beta: Vec<G1> = g1_table.batch_mul(&beta_scalars).into_iter().map(G1::from).collect();
    let z: Vec<G1> = g1_table.batch_mul(&z_scalars).into_iter().map(G1::from).collect();
    let h: Vec<G1> = g1_table.batch_mul(&h_scalars).into_iter().map(G1::from).collect();
    let b: Vec<G1> = g1_table.batch_mul(&b_scalars).into_iter().map(G1::from).collect();
    let lx_g: Vec<G1> = g1_table.batch_mul(&lx).into_iter().map(G1::from).collect();
    let ly_g: Vec<G1> = g1_table.batch_mul(&ly).into_iter().map(G1::from).collect();
    eprintln!("[setup] fixed-base msm: {:?}", t.elapsed());

    let pk = ProvingKey {
        alpha,
        beta,
        z,
        h,
        b,
        lx: lx_g,
        ly: ly_g,
        u: g1 * u,
        x: g1 * tau_x,
        x2: g1 * (tau_x * tau_x),
        xm: g1 * tau_x_m,
        ym: g1 * tau_y_m,
        xm_ym: g1 * (tau_x_m * tau_y_m),
        x2_ym: g1 * (tau_x * tau_x * tau_y_m),
    };

    let vk = VerifyingKey {
        u: g2 * u,
        x: g2 * tau_x,
        x2: g2 * (tau_x * tau_x),
        xm: g2 * tau_x_m,
        ym: g2 * tau_y_m,
        invm: g2 * inv_m,
    };

    eprintln!("[setup] total: {:?} (m = {m})", setup_start.elapsed());

    (pk, vk)
}
