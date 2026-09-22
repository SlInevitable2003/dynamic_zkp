use crate::common::curve::{Engine, G1, G2, GT};
use ark_ec::pairing::Pairing;

pub fn pair(p: &G1, q: &G2) -> GT { Engine::pairing(*p, *q) }
