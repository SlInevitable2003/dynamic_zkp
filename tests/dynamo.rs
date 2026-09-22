//! Dynamo(非通用版)集成测试。

use ark_ec::PrimeGroup;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use dynamic_zkp::common::curve::{G1, Scalar};
use dynamic_zkp::protocols::dynamo::non_universal::{derive_h, prove, setup, verify};

fn random_permutation<R: Rng>(m: usize, rng: &mut R) -> Vec<usize> {
    let mut sigma: Vec<usize> = (0..m).collect();
    for i in (1..m).rev() {
        let j = rng.gen_range(0..=i);
        sigma.swap(i, j);
    }
    sigma
}

#[test]
fn prove_verify_roundtrip() {
    let mut rng = ark_std::test_rng();
    for m in [4usize, 8, 16] {
        let sigma = random_permutation(m, &mut rng);
        let (pk, vk) = setup(m, &sigma, &mut rng);
        let z: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();
        let h = derive_h(&z, &sigma);

        let (x, pi) = prove(&pk, &z, &h, &mut rng);
        assert!(verify(&vk, &x, &pi), "roundtrip 失败:m = {m}");
    }
}

#[test]
fn tampered_proof_is_rejected() {
    let mut rng = ark_std::test_rng();
    let m = 8;
    let sigma = random_permutation(m, &mut rng);
    let (pk, vk) = setup(m, &sigma, &mut rng);
    let z: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();
    let h = derive_h(&z, &sigma);

    let (x, pi) = prove(&pk, &z, &h, &mut rng);
    assert!(verify(&vk, &x, &pi));

    let mut bad = pi;
    bad.alpha = bad.alpha + G1::generator();
    assert!(!verify(&vk, &x, &bad), "篡改 alpha 必须被拒绝");
}

#[test]
fn witness_violating_relation_is_rejected() {
    let mut rng = ark_std::test_rng();
    let m = 8;
    let sigma = random_permutation(m, &mut rng);
    let (pk, vk) = setup(m, &sigma, &mut rng);

    // 不满足 h_i = z_i - z_{σ(i)} 的随机见证应无法通过验证
    let z: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();
    let h: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();

    let (x, pi) = prove(&pk, &z, &h, &mut rng);
    assert!(!verify(&vk, &x, &pi), "违反关系的见证必须被拒绝");
}
