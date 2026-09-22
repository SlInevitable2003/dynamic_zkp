use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
use ark_ec::pairing::PairingOutput;

pub type Engine = Bls12_381;
pub type Scalar = Fr;
pub type G1 = G1Projective;
pub type G2 = G2Projective;
pub type GT = PairingOutput<Engine>;
