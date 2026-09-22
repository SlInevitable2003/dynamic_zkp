//! Dynamo 非通用版(README §1.1)。

mod prove;
mod relation;
mod setup;
mod types;
mod verify;

pub use prove::prove;
pub use relation::derive_h;
pub use setup::setup;
pub use types::{Instance, Proof, ProvingKey, VerifyingKey};
pub use verify::verify;
