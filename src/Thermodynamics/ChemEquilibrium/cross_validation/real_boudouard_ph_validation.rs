//! Test-side I4 reference-state construction for real Boudouard `P,H` cases.
//!
//! The immutable `BoudouardCarbon` fixture continues to own record lookup,
//! component identities, and reaction data. This narrow helper owns only one
//! validation concern: derive a reproducible `H_target` from an independently
//! solved fixed-temperature interior extent state before any canonical `P,H`
//! route runs.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::PurePhaseBoundaryStructuralTolerances;
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhMaterializedState, PurePhasePhProblem, PurePhasePhSolverSettings,
};
use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
    RealPurePhaseInventory, ResolvedRealPurePhaseFixture,
};

/// Boudouard I4 validation uses one explicit 1-bar convention throughout.
pub(crate) const BOUDOUARD_PH_PRESSURE_PA: f64 = 100_000.0;
/// A local interior carbon-containing state used only as a reproducible seed.
pub(crate) const BOUDOUARD_PH_REFERENCE_TEMPERATURE_K: f64 = 700.0;
const BOUDOUARD_PH_TEMPERATURE_LOWER_K: f64 = 695.0;
const BOUDOUARD_PH_TEMPERATURE_UPPER_K: f64 = 705.0;
const MIN_INTERIOR_GRAPHITE_MOLES: f64 = 1e-6;

/// Independent real-data reference for one interior Boudouard `P,H` point.
///
/// The `problem` carries a target enthalpy calculated from `state`; no
/// canonical P,T/P,H solution contributes to either value.
#[derive(Clone)]
pub(crate) struct BoudouardPhInteriorReference {
    /// Gas-only inventory whose independent extent root creates graphite.
    /// Keeping it here lets a later I3 story start from exactly the same C/O
    /// inventory rather than reconstructing a different gas composition.
    pub(crate) inventory: RealPurePhaseInventory,
    pub(crate) problem: PurePhasePhProblem,
    pub(crate) state: PurePhasePhMaterializedState,
}

impl BoudouardPhInteriorReference {
    pub(crate) fn target_enthalpy(&self) -> f64 {
        self.state.total_enthalpy
    }
}

/// Builds a physically interior Boudouard reference state from local closures.
///
/// The first temporary problem carries a harmless finite target only because
/// the generic type models P,H conditions. Its scalar inner `ln(Q)-ln(K)` root
/// is entirely independent of that target. The accepted inner state then
/// supplies the real additive target for the returned P,H problem.
pub(crate) fn build_boudouard_ph_interior_reference(
    fixture: &ResolvedRealPurePhaseFixture,
) -> Result<BoudouardPhInteriorReference, ReactionExtentError> {
    // Both gases must be present for a finite ideal-gas logarithm. This
    // gas-only inventory is deliberately CO-rich, so the independently solved
    // Boudouard extent produces an interior graphite amount at 700 K while
    // retaining exactly the elemental inventory that a later I3 activation
    // route receives.
    let inventory = RealPurePhaseInventory::new(vec![1.19, 0.005], 0.0)?;
    let bounds = TemperatureBounds::new(
        BOUDOUARD_PH_TEMPERATURE_LOWER_K,
        BOUDOUARD_PH_TEMPERATURE_UPPER_K,
    )?;
    let seed_problem = fixture.to_ph_problem(
        &inventory,
        BOUDOUARD_PH_PRESSURE_PA,
        BOUDOUARD_PH_PRESSURE_PA,
        0.0,
        bounds,
    )?;
    seed_problem
        .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())?;
    let inner = seed_problem.solve_inner_extent_at_temperature(
        BOUDOUARD_PH_REFERENCE_TEMPERATURE_K,
        PurePhasePhSolverSettings::default(),
    )?;
    let state = seed_problem
        .materialize_state_at_extent(inner.extent, BOUDOUARD_PH_REFERENCE_TEMPERATURE_K)?;
    if state.candidate_moles <= MIN_INTERIOR_GRAPHITE_MOLES {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "real_boudouard_ph_interior_reference",
            message: format!(
                "local Boudouard scalar root at {} K has non-interior graphite amount {:e}",
                BOUDOUARD_PH_REFERENCE_TEMPERATURE_K, state.candidate_moles
            ),
        });
    }
    let problem = seed_problem.with_target_enthalpy(state.total_enthalpy)?;
    Ok(BoudouardPhInteriorReference {
        inventory,
        problem,
        state,
    })
}
