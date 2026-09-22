use crate::common::curve::{G1, Scalar};
use crate::common::msm::msm;
use crate::protocols::dynamo::non_universal::types::{Instance, Proof, ProvingKey};

use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_std::rand::Rng;
use ark_std::UniformRand;

pub fn prove<R: Rng>(
    pk: &ProvingKey,
    z: &[Scalar],
    h: &[Scalar],
    rng: &mut R,
) -> (Instance, Proof) {
    let prove_start = std::time::Instant::now();

    let m = z.len();
    assert_eq!(h.len(), m, "the length of z must match that of h");
    assert_eq!(pk.alpha.len(), m, "pk is generated for different m");

    let rho_z = Scalar::rand(rng);
    let rho_h = Scalar::rand(rng);
    let rho_qv = Scalar::rand(rng);
    let rho_qh = Scalar::rand(rng);

    let inv_m = Scalar::from(m as u64).inverse().expect("m is nonzero");

    let alpha = msm_timed("alpha", &pk.alpha, z); // [α]_1
    let beta = msm_timed("beta", &pk.beta, z); // [β]_1
    let zcomb = msm_timed("Z", &pk.z, z); // [Z]_1
    let hcomb = msm_timed("H", &pk.h, z); // [H]_1
    let bcomb = msm_timed("B", &pk.b, z); // [B]_1
    let z_lin = msm_timed("z", &pk.lx, z); // [z]_1
    let h_lin = msm_timed("h", &pk.ly, h); // [h]_1

    let g1 = G1::generator();
    let xm1 = pk.xm - g1; // G1 上的 τ_X^m - 1
    let ym1 = pk.ym - g1; // G1 上的 τ_Y^m - 1

    let z_zk = z_lin + xm1 * rho_z;
    let h_zk = h_lin + ym1 * rho_h;

    let q_zk = g1 * (-rho_h * inv_m) + xm1 * rho_qv + pk.x * rho_qh;
    let alpha_zk = alpha + pk.u * rho_z - ym1 * rho_qv;
    let beta_zk = beta - ym1 * rho_qh;
    let z_zk_p = zcomb + (pk.xm_ym - pk.ym) * rho_z; // [Z_zk]
    let h_zk_p = hcomb + (pk.xm_ym - pk.xm) * rho_h; // [H_zk]
    let b_zk = bcomb - (pk.x2_ym - pk.x2) * rho_qh; // [B_zk]

    let instance = Instance { z: z_zk, h: h_zk };
    let proof = Proof {
        alpha: alpha_zk,
        beta: beta_zk,
        q: q_zk,
        z: z_zk_p,
        h: h_zk_p,
        b: b_zk,
    };

    eprintln!("[prove] total: {:?} (m = {m})", prove_start.elapsed());

    (instance, proof)
}

fn msm_timed(label: &str, bases: &[G1], scalars: &[Scalar]) -> G1 {
    let start = std::time::Instant::now();
    let out = msm(bases, scalars);
    eprintln!("[prove] msm[{label}]: {:?}", start.elapsed());
    out
}
