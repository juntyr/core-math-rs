//! [![CI Status]][workflow] [![MSRV]][repo] [![Rust Doc Main]][docs]
//!
//! [CI Status]: https://img.shields.io/github/actions/workflow/status/juntyr/core-math-rs/ci.yml?branch=main
//! [workflow]: https://github.com/juntyr/core-math-rs/actions/workflows/ci.yml?query=branch%3Amain
//!
//! [MSRV]: https://img.shields.io/badge/MSRV-1.85.0-blue
//! [repo]: https://github.com/juntyr/core-math-rs
//!
//! [Rust Doc Main]: https://img.shields.io/badge/docs-main-blue
//! [docs]: https://juntyr.github.io/core-math-rs/core_math_rs
//!
//! # core-math-rs
//!
//! In-progress Rust port of the [CORE-MATH] project, which provides
//! on-the-shelf open-source mathematical functions with correct rounding.
//!
//! ## Cargo Features
//!
//! - `f16`: enable nightly support for `f16` functions
//!
//! ## Related Projects
//!
//! - [`core-math`]: safe Rust bindings to the [CORE-MATH] C library,
//!   which we test against
//!
//! [CORE-MATH]: https://core-math.gitlabpages.inria.fr/
//! [`core-math`]: https://github.com/jdh8/core-math

#![cfg_attr(feature = "f16", feature(f16))]
#![cfg_attr(feature = "f16", feature(f32_from_f16))]

#[cfg(feature = "f16")]
pub mod f16;
