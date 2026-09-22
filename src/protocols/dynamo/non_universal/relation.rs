use crate::common::curve::Scalar;

pub fn derive_h(z: &[Scalar], sigma: &[usize]) -> Vec<Scalar> {
    assert_eq!(z.len(), sigma.len(), "the length of z must match that of sigma");
    (0..z.len()).map(|i| z[i] - z[sigma[i]]).collect()
}
