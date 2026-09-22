//! 端到端 demo:随机 σ / z / h,跑通 Dynamo(非通用版)setup → prove → verify。

use ark_ec::PrimeGroup;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use dynamic_zkp::common::curve::{G1, Scalar};
use dynamic_zkp::protocols::dynamo::non_universal::{derive_h, prove, setup, verify};

/// Fisher–Yates 生成 `(0..m)` 上的随机置换。
fn random_permutation<R: Rng>(m: usize, rng: &mut R) -> Vec<usize> {
    let mut sigma: Vec<usize> = (0..m).collect();
    for i in (1..m).rev() {
        let j = rng.gen_range(0..=i);
        sigma.swap(i, j);
    }
    sigma
}

fn main() {
    let m = 16; // 2 的幂
    let mut rng = ark_std::test_rng();

    let sigma = random_permutation(m, &mut rng);
    let (pk, vk) = setup(m, &sigma, &mut rng);

    let z: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();
    // 见证须满足 Dynamo 关系:h_i = z_i - z_{σ(i)}
    let h = derive_h(&z, &sigma);

    let (instance, proof) = prove(&pk, &z, &h, &mut rng);
    assert!(verify(&vk, &instance, &proof), "有效证明必须通过验证");
    println!("Dynamo (non-universal) proof verified ✓  (m = {m})");

    // 负例:篡改证明后应被拒绝
    let mut bad = proof;
    bad.alpha = bad.alpha + G1::generator();
    assert!(!verify(&vk, &instance, &bad), "篡改后的证明必须被拒绝");
    println!("tampered proof rejected ✓");
}
