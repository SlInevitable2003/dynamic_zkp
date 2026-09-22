//! 集中类型别名:换曲线 / 未来抽象成泛型时只改这里。

use ark_bls12_381::{Bls12_381, Fr, G1Projective, G2Projective};
use ark_ec::pairing::PairingOutput;

/// 配对引擎(实现 [`ark_ec::pairing::Pairing`])。
pub type Engine = Bls12_381;

/// 标量域 `F`(= BLS12-381 的标量域,协议里的素数域)。
pub type Scalar = Fr;

/// `G1`(投影表示),`[x]_1 = G1::generator() * x`。
pub type G1 = G1Projective;

/// `G2`(投影表示),`[x]_2 = G2::generator() * x`。
pub type G2 = G2Projective;

/// 目标群 `G_T`。arkworks 用加法表示该群的群运算(等价于底层 `Fq12` 的乘法)。
pub type GT = PairingOutput<Engine>;
