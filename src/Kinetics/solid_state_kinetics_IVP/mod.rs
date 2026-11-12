//! Solid-state kinetic IVP models and solver facade.
//!
//! Models describe conversion-rate functions `f(a)` and the solver combines
//! them with an Arrhenius heating profile:
//!
//! ```text
//! da/dt = A * exp(-E / (R * (T0 + beta*t))) * f(a)
//! ```
//!
//! The scalar equation is prepared through RustedSciThe's `StiffAPI` using its
//! `solve_1d()` preset. Supported integrations are BDF, Radau, Backward Euler,
//! and LSODE2. Non-stiff methods are not exposed as choices for this workflow.
//!
//! The public module path remains `Kinetics::solid_state_kinetics_IVP`; model
//! definitions, integration lifecycle, and tests live in focused files.
//!
//! ```rust,ignore
//! use KiThe::Kinetics::solid_state_kinetics_IVP::{KineticModelIVP, KineticModelNames};
//! use RustedSciThe::numerical::StiffAPI::StiffSolverKind;
//!
//! let mut ivp = KineticModelIVP::new(StiffSolverKind::Bdf);
//! ivp.set_problem(100.0, 10.0, 298.15, 50_000.0, 1e6)?;
//! ivp.set_model(KineticModelNames::JMA, vec![2.0])?;
//! ivp.solve()?;
//! # Ok::<(), String>(())
//! ```

mod models;
mod options;
mod solution;
mod solver;

pub use models::*;
pub use options::SolidStateIvpSolverOptions;
pub use solution::IVPSolution;
pub use solver::KineticModelIVP;
