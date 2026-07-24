//! Reactor IVP module.
//!
//! This namespace contains the condensed-phase combustion IVP model used by
//! KiThe. The public surface is intentionally split into a few stable entry
//! points:
//!
//! - [`SimpleReactorIVP`] for the reactor-domain task, validation, and solve
//!   orchestration;
//! - [`solver_backend`] for the typed LSODE2 facade and solver-policy mapping;
//! - [`task_parser_reactor_IVP`] for reactor physics and task-document parsing;
//! - focused contract, solve, story, AOT, and GUI regression tests.
//!
//! The model is condensed-phase only: density is user-supplied and constant,
//! species diffusion is absent, and LSODE2 is the only production solver
//! family. The supported user-facing choices are exposed through typed enums
//! rather than free-form strings wherever possible.

pub mod SimpleReactorIVP;
pub mod reactor_ivp_aot_story_tests;
pub mod reactor_ivp_contract_tests;
pub mod reactor_ivp_quality_tests;
pub mod reactor_ivp_solve_tests;
pub mod reactor_ivp_solver_backend_tests;
pub mod reactor_ivp_story_tests;
pub mod solver_backend;
pub mod task_parser_reactor_IVP;
