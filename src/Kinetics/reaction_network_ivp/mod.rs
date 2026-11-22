//! Homogeneous, constant-volume reaction-network kinetics in physical time.
//!
//! Concentrations are mol/m^3, time is seconds, and temperature is a prescribed
//! function in kelvin. No thermal or thermodynamic balance is solved here.

mod assembly;
mod diagnostics;
mod error;
mod problem;
mod rate_law;
mod solid_state;
mod solver;
mod solver_backend;
mod temperature;

pub use assembly::AssembledNetwork;
pub use diagnostics::{NetworkIvpStage, NetworkIvpStageDiagnostics, NetworkIvpStageMeasurement};
pub use error::NetworkIvpError;
pub use problem::{NetworkProblem, ReactionSpec, SpeciesSpec, StoichTerm};
pub use rate_law::{KineticOrder, RateCoefficient, RateCoefficientBinding, RateLaw};
pub use solver::{PreparedNetworkIvp, solve_network};
pub use solver_backend::{
    NetworkIvpOptions, NetworkIvpPartialTrajectory, NetworkIvpRunMetadata, NetworkIvpSolution,
};
pub use temperature::{TemperatureBinding, TemperatureProgram};

#[cfg(test)]
mod solid_state_tests;
#[cfg(test)]
mod tests;
