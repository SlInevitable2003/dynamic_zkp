use crate::common::curve::G2;
use crate::common::pairing::pair;
use crate::protocols::dynamo::non_universal::types::{Instance, Proof, VerifyingKey};
use ark_ec::PrimeGroup;

pub fn verify(vk: &VerifyingKey, x: &Instance, pi: &Proof) -> bool {
    let g2 = G2::generator();
    let xm1 = vk.xm - g2;
    let ym1 = vk.ym - g2;

    // e(z_zk, u) == e(α_zk, Xm1) + e(β_zk, X) + e(h_zk, invm) + e(q_zk, Ym1)
    let lhs = pair(&x.z, &vk.u);
    let rhs = pair(&pi.alpha, &xm1)
        + pair(&pi.beta, &vk.x)
        + pair(&x.h, &vk.invm)
        + pair(&pi.q, &ym1);
    
    if lhs != rhs { return false; }

    // e(z_zk, Ym) == e(Z_zk, g2)
    if pair(&x.z, &vk.ym) != pair(&pi.z, &g2) { return false; }

    // e(h_zk, Xm) == e(H_zk, g2)
    if pair(&x.h, &vk.xm) != pair(&pi.h, &g2) { return false; }

    // e(β_zk, X2) == e(B_zk, g2)
    if pair(&pi.beta, &vk.x2) != pair(&pi.b, &g2) { return false; }

    true
}
