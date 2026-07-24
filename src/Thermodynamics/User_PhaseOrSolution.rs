//! # Multi-Phase Thermodynamics Module
//!
//! This module provides thermodynamic calculations for systems with multiple phases or solutions.
//! It handles complex systems where substances can exist in different phases (gas, liquid, solid)
//! or different solutions simultaneously.
//!
//! ## Key Structures
//!
//! - [`PhaseOrSolution`]: Manages multiple phases, each with their own substance data
//! - [`CustomSubstance`]: Enum that unifies single-phase and multi-phase systems
//!
//! ## Hierarchy
//! ```text
//!                             MoleNumberSnapshot
//!                          /            |                \
//!                 SystemLayout     PhaseMoleNumbers   OrderedPhaseMoles
//!            /      |         \           |                   |
//!   phases   components  phase_ranges    input types          ordered vectors
//!                                        Contains total
//!                                       and component
//!                                       mole numbers for
//!                                        each phase
//! ```
//! ## Phase Key Convention
//!
//! - **Multi-phase systems**: Use `Some("phase_name")` as keys (e.g., `Some("gas")`, `Some("liquid")`)
//! - **Single-phase systems**: Use `None` as the key to represent the single phase
//!
//! ## Example Usage
//!
//! ```rust
//! use KiThe::Thermodynamics::User_PhaseOrSolution::PhaseOrSolution;
//! use std::collections::HashMap;
//!
//! // Multi-phase system with gas and liquid phases
//! let mut system = PhaseOrSolution::new();
//! // Typed builders resolve phase payloads into the shared `PhaseSystem`.
//! // Read-only inspection uses `system.phase_data_view()`.
//!
//! // Calculate Gibbs energy for all phases
//! // system.calculate_Gibbs_sym(298.15)?;
//! ```

use std::sync::Arc;

use crate::Thermodynamics::phase_layout::SystemLayout;
use RustedSciThe::symbolic::symbolic_engine::Expr;

use std::collections::HashMap;
use std::f64;
pub const R: f64 = 8.314;
#[allow(non_upper_case_globals)]
pub const R_sym: Expr = Expr::Const(R);
pub type PhaseFunction = Arc<dyn Fn(f64, Option<Vec<f64>>, Option<f64>) -> f64 + 'static>;
pub type PhaseLagrangeFunction =
    Box<dyn Fn(f64, Option<Vec<f64>>, Option<f64>, Vec<f64>) -> Vec<f64> + 'static>;

/// Legacy tuple shape retained only for `ThermodynamicsCalculatorTrait`.
///
/// New internal code reads the named `PhaseSymbolicLayout` snapshot instead.
pub(crate) type IndexedMoleVariables = (
    HashMap<Option<String>, (Option<Expr>, Option<Vec<Expr>>)>,
    Vec<Expr>,
    Vec<Expr>,
    HashMap<Option<String>, HashMap<String, Expr>>,
);

/// Legacy mole-number shape retained only at compatibility boundaries.
pub(crate) type LegacyMoleNumberSnapshot = (
    HashMap<Option<String>, (Option<f64>, Option<HashMap<String, f64>>)>,
    HashMap<Option<String>, (Option<f64>, Option<Vec<f64>>)>,
    HashMap<String, f64>,
);

/// Sparse per-phase mole input after missing declared components have been
/// normalized to zero. The map deliberately remains keyed by substance name;
/// solver-facing ordering belongs to the companion vector snapshot.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseMoleNumbers {
    total_amount: Option<f64>,
    component_amounts: Option<HashMap<String, f64>>,
}

impl PhaseMoleNumbers {
    pub fn new(total_amount: Option<f64>, component_amounts: Option<HashMap<String, f64>>) -> Self {
        Self {
            total_amount,
            component_amounts,
        }
    }

    pub fn total_amount(&self) -> Option<f64> {
        self.total_amount
    }

    pub fn component_amounts(&self) -> Option<&HashMap<String, f64>> {
        self.component_amounts.as_ref()
    }

    fn into_legacy(self) -> (Option<f64>, Option<HashMap<String, f64>>) {
        (self.total_amount, self.component_amounts)
    }
}

/// Ordered mole amounts for one phase. Its vector is aligned with the phase's
/// component range in `MoleNumberSnapshot::layout`.
#[derive(Clone, Debug, PartialEq)]
pub struct OrderedPhaseMoles {
    total_amount: Option<f64>,
    component_amounts: Vec<f64>,
}

impl OrderedPhaseMoles {
    pub fn total_amount(&self) -> Option<f64> {
        self.total_amount
    }

    pub fn component_amounts(&self) -> &[f64] {
        &self.component_amounts
    }
}

/// Typed, immutable result of normalizing sparse phase mole input.
///
/// The snapshot carries the exact layout used to build vectors, so callers do
/// not need to infer component order from a `HashMap` or a facade variant.
#[derive(Clone, Debug, PartialEq)]
pub struct MoleNumberSnapshot {
    layout: SystemLayout,
    phase_amounts: HashMap<Option<String>, PhaseMoleNumbers>,
    ordered_phase_amounts: HashMap<Option<String>, OrderedPhaseMoles>,
    total_amounts_by_substance: HashMap<String, f64>,
}

impl MoleNumberSnapshot {
    pub fn layout(&self) -> &SystemLayout {
        &self.layout
    }

    pub fn phase_amounts(&self) -> &HashMap<Option<String>, PhaseMoleNumbers> {
        &self.phase_amounts
    }

    pub fn ordered_phase_amounts(&self) -> &HashMap<Option<String>, OrderedPhaseMoles> {
        &self.ordered_phase_amounts
    }

    pub fn total_amounts_by_substance(&self) -> &HashMap<String, f64> {
        &self.total_amounts_by_substance
    }

    fn into_legacy(self) -> LegacyMoleNumberSnapshot {
        let sparse_amounts = self
            .phase_amounts
            .into_iter()
            .map(|(phase, amounts)| (phase, amounts.into_legacy()))
            .collect();
        let ordered_amounts = self
            .ordered_phase_amounts
            .into_iter()
            .map(|(phase, amounts)| {
                (
                    phase,
                    (amounts.total_amount, Some(amounts.component_amounts)),
                )
            })
            .collect();
        (
            sparse_amounts,
            ordered_amounts,
            self.total_amounts_by_substance,
        )
    }
}

#[path = "phase_moles.rs"]
mod phase_moles;
pub(crate) use phase_moles::{build_indexed_mole_variables, normalize_mole_numbers};

#[path = "phase_elements.rs"]
mod phase_elements;
pub(crate) use phase_elements::element_composition_and_molar_mass;

#[path = "phase_operations.rs"]
mod phase_operations;
pub(crate) use phase_operations::stage_phase_data_operation;

#[path = "phase_property_builders.rs"]
mod phase_property_builders;
pub(crate) use phase_property_builders::{
    stage_entropy_functions, stage_gibbs_functions, stage_symbolic_entropy, stage_symbolic_gibbs,
};

#[path = "phase_factory.rs"]
mod phase_factory;
pub use phase_factory::{
    SubstancePhaseMapping, SubstanceSystemFactory, SubstanceSystemFactoryError,
    SubstanceSystemSpec, SubstanceSystemSpecBuilder, SubstancesContainer,
};

#[path = "phase_domain.rs"]
mod phase_domain;
pub use phase_domain::{
    PhaseModel, PhasePhysicalState, PhaseResolutionSummary, PhaseSpec, ResolvedPhaseSystem,
    ResolvedPhaseSystemReport,
};

#[path = "phase_evaluation.rs"]
mod phase_evaluation;
pub(crate) use phase_evaluation::{
    NumericPhaseProperty, PhaseDataView, evaluate_numeric_phase_property,
    validate_phase_evaluation_request,
};
pub use phase_evaluation::{
    NumericThermoCacheContext, PhaseComposition, PhaseEvaluationRequest, ThermoEvaluationConditions,
};

#[path = "phase_cache.rs"]
mod phase_cache;
pub use phase_cache::{
    NestedPhaseCacheIter, NestedPhaseCacheView, PhaseThermoPropertyIter, PhaseThermoPropertyView,
    ThermoCacheSnapshot, ThermoResultSnapshot, ThermoStateSnapshot,
};
pub(crate) use phase_cache::{
    PhaseSymbolicLayout, PhaseThermoCacheBundle, build_cache_snapshot, build_multi_cache_snapshot,
    build_single_cache_snapshot, build_thermo_result_snapshot, build_thermo_state_snapshot,
    substitute_exprs,
};

#[path = "phase_system.rs"]
mod phase_system;
pub(crate) use phase_system::PhaseSystem;

#[path = "phase_facade.rs"]
mod phase_facade;
pub use phase_facade::CustomSubstance;

#[path = "phase_interfaces.rs"]
mod phase_interfaces;
pub use phase_interfaces::{
    PhaseDataPreparation, PhaseEquilibriumAssembly, PhaseLayoutAccess, PhasePropertyEvaluator,
    PhaseSymbolicPropertyBuilder,
};

#[path = "phase_system_facade.rs"]
mod phase_system_facade;
pub use phase_system_facade::PhaseOrSolution;

#[path = "phase_legacy_api.rs"]
mod phase_legacy_api;
pub use phase_legacy_api::ThermodynamicsCalculatorTrait;
