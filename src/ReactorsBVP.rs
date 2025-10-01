//! # Reactor Boundary Value Problem (BVP) Module
//!
//! This module provides solvers for steady-state gas-phase combustion and plug-flow reactor problems
//! using boundary value problem formulations with dimensionless variables.
//!
//! ## Mathematical Model
//!
//! ### Nomenclature
//!
//! | Symbol | Description | Units |
//! |--------|-------------|-------|
//! | `l` | Length scale | m |
//! | `dT` | Temperature displacement (dimensionless) | K |
//! | `Ts` | Temperature scale | K |
//! | `Tm` | Characteristic temperature at combustion front | K |
//! | `Cp` | Heat capacity | J/(kg·K) |
//! | `λ` | Thermal conductivity | W/(m·K) |
//! | `ρ` | Gas density | kg/m³ |
//!
//! ### Model Assumptions
//!
//! - Heat capacity, thermal conductivity, thermal effects, and gas density are constants
//!   calculated at characteristic temperature `Tm`
//! - Only elementary kinetic expressions are considered (no fall-off, third-body reactions)
//!
//! ### Governing Equations
//!
//! **Dimensional form:**
//! ```text
//! m(dc_i)/dx = d/dx[ρD(dc_i)/dx] + Σ_j[ν_ij w_j]
//! mCp(dT)/dx = d/dx[λ(dT)/dx] + Σ_j[Q_j w_j]
//! ```
//!
//! **Dimensionless form:**
//!
//! With dimensionless variables `Θ = (T-dT)/Ts` and `ξ = x/l`:
//!
//! ```text
//! dΘ/dξ = q/λ
//! dq/dξ = Pe_T·q - a_T·Σ_j[Q_j w_j(Θ)]
//! dc_i/dξ = j_i/(ρD)
//! dj_i/dξ = Pe_c·j_i - a_c·Σ_j[ν_ij w_j]
//! ```
//!
//! Where:
//! - `Pe_T = lmCp/λ` (thermal Peclet number)
//! - `Pe_c = ml/(ρD)` (mass Peclet number)
//! - `a_T = l²/Ts`, `a_c = l²`
//!
//! ### Boundary Conditions
//!
//! - At `ξ = 0`: `c_i = c_i0`, `Θ = Θ0 = (T0-dT)/Ts`
//! - At `ξ = ∞`: `q = j_i = 0` (zero fluxes)
//!
//! ## Numerical Solution
//!
//! The system forms a boundary value problem with `2n+2` first-order ODEs for `n` substances.
//! Solutions are obtained using Newton-Raphson methods from the [RustedSciThe](https://crates.io/crates/RustedSciThe) 
//! package.

pub mod SimpleReactorBVP;
pub mod SimpleReactorBVP2;
pub mod SimpleReactorBVP3;
pub mod createBVP;
pub mod reactor_BVP_utils;
mod simple_reactor_bvp_tests;
mod simple_reactor_bvp_tests2;
pub mod task_parser_reactor_BVP;