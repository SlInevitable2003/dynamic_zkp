//! `dynamic_zkp`:基于 arkworks 实现 EUROCRYPT 2026 "Dynamic zk-SNARKs" 的学习项目。
//!
//! 结构:
//! - [`common`] —— 协议无关的共享原语(曲线类型别名、多项式、MSM、配对助手)。
//! - [`protocols`] —— 各协议实现,每个协议一个子模块(目前含 Dynamo 非通用版)。

pub mod common;
pub mod protocols;
