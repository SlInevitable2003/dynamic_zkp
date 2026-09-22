use ark_std::rand::Rng;
use ark_std::UniformRand;
use dynamic_zkp::common::curve::Scalar;
use dynamic_zkp::protocols::dynamo::non_universal::{derive_h, prove, setup, verify};

/// Fisher–Yates 生成 `(0..m)` 上的随机置换
fn random_permutation<R: Rng>(m: usize, rng: &mut R) -> Vec<usize> {
    let mut sigma: Vec<usize> = (0..m).collect();
    for i in (1..m).rev() {
        let j = rng.gen_range(0..=i);
        sigma.swap(i, j);
    }
    sigma
}

fn main() {
    let m = 1024;
    let mut rng = ark_std::test_rng();

    let sigma = random_permutation(m, &mut rng);
    let (pk, vk) = setup(m, &sigma, &mut rng);

    let z: Vec<Scalar> = (0..m).map(|_| Scalar::rand(&mut rng)).collect();
    let h = derive_h(&z, &sigma);

    let (instance, proof) = prove(&pk, &z, &h, &mut rng);
    assert!(verify(&vk, &instance, &proof), "A valid proof can be verified");
    println!("Dynamo (non-universal) proof verified (m = {m}).");
}
