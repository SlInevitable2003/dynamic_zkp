use crate::common::curve::Scalar;
use ark_ff::{FftField, Field};

fn log2_of_power_of_two(m: usize) -> u32 {
    assert!(m.is_power_of_two(), "m must be a power of 2, got {m}");
    m.trailing_zeros()
}

pub fn root_of_unity(m: usize) -> Scalar {
    let _ = log2_of_power_of_two(m);
    Scalar::get_root_of_unity(m as u64).expect("m does not divide the order of the multiplicative group of F")
}

pub fn lagrange_evals(tau: Scalar, m: usize) -> Vec<Scalar> {
    
    let w = root_of_unity(m);
    let w_inv = w.inverse().expect("root of unity is nonzero");
    let mut w_pow = Scalar::ONE;

    let num = tau.pow([m as u64]) - Scalar::ONE;
    let inv_m = Scalar::from(m as u64).inverse().expect("m is nonzero");
    (0..m)
        .map(|_| {
            let denom = w_pow * tau - Scalar::ONE;
            let li = num * denom.inverse().expect("τ != ω^i") * inv_m;
            w_pow *= w_inv;
            li
        })
        .collect()
}

pub fn inverse_permutation(sigma: &[usize]) -> Vec<usize> {
    let m = sigma.len();
    let mut inv = vec![usize::MAX; m];
    for (i, &s) in sigma.iter().enumerate() {
        assert!(s < m, "sigma[{i}] = {s} out of bound");
        assert_eq!(inv[s], usize::MAX, "sigma has replicated element {s}");
        inv[s] = i;
    }
    inv
}

pub fn y_evals(ly: &[Scalar], sigma_inv: &[usize]) -> Vec<Scalar> {
    debug_assert_eq!(ly.len(), sigma_inv.len());
    (0..ly.len()).map(|i| ly[i] - ly[sigma_inv[i]]).collect()
}

pub fn u_eval(lx: &[Scalar], y: &[Scalar]) -> Scalar {
    debug_assert_eq!(lx.len(), y.len());
    lx.iter().zip(y).map(|(a, b)| *a * *b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{AdditiveGroup, Field};
    use ark_std::UniformRand;

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
        for i in 0..sigma.len() { assert_eq!(inv[sigma[i]], i); }
    }
}
