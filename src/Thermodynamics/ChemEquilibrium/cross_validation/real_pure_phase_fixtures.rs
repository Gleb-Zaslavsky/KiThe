//! Shared chemistry declarations for offline real pure-phase validation.
//!
//! A declaration describes one deliberately restricted chemistry family, not
//! a `P,T` or `P,H` scenario. Later adapters resolve it once through the
//! ordinary phase API and reuse the same record identities, provenance, and
//! thermochemical capabilities on each validation route.

use std::collections::HashMap;
use std::rc::Rc;
use std::sync::Arc;

use nalgebra::DMatrix;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_log_moles::GibbsFn;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::{
    MolarThermoFunction, ResolvedThermochemistry,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
    PurePhaseBoundaryElementComposition, PurePhaseBoundaryProblem,
    PurePhaseBoundaryStructuralTolerances,
};
use crate::Thermodynamics::ChemEquilibrium::pure_phase_ph_validation::{
    PurePhasePhConditions, PurePhasePhProblem,
};
use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseModel, ResolvedPhaseSystem, SubstanceSystemFactory, SubstanceSystemSpec,
    SubstanceSystemSpecBuilder, SubstancesContainer,
};
use crate::Thermodynamics::physical_state::PhysicalState;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// One supported restricted chemistry family for real pure-phase tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum RealPurePhaseFamily {
    WaterLiquid,
    WaterIce,
    BoudouardCarbon,
    MethaneCarbon,
}

/// Immutable chemistry-only declaration. Scenario-dependent moles,
/// temperature, and target enthalpy deliberately do not belong here.
#[derive(Debug, Clone)]
pub(crate) struct RealPurePhaseDefinition {
    pub(crate) name: &'static str,
    pub(crate) gas_phase: &'static str,
    pub(crate) candidate_phase: &'static str,
    pub(crate) gas_species: &'static [&'static str],
    pub(crate) candidate_species: &'static str,
    pub(crate) gas_stoichiometry: &'static [f64],
    pub(crate) candidate_stoichiometry: f64,
    pub(crate) elements: &'static [&'static str],
    /// Rows are gas components and columns are elements.
    pub(crate) gas_element_rows: &'static [f64],
    pub(crate) candidate_elements: &'static [f64],
    pub(crate) candidate_state: PhysicalState,
    pub(crate) gas_library: &'static str,
    pub(crate) candidate_library: &'static str,
}

/// Immutable local data resolved for one real pure-phase family.
///
/// Scenario builders borrow this payload so P,T and P,H validators cannot
/// accidentally re-resolve a different library record for the same chemistry.
#[derive(Clone)]
pub(crate) struct ResolvedRealPurePhaseFixture {
    family: RealPurePhaseFamily,
    resolved: ResolvedPhaseSystem,
    thermochemistry: ResolvedThermochemistry,
}

/// One inert gas explicitly attached to a validation scenario.
///
/// The family owns only the phase-forming reaction. An inert belongs to a
/// particular initial inventory, therefore it carries its own amount and
/// elemental row instead of being hidden in a water/carbon fixture.
#[derive(Debug, Clone)]
pub(crate) struct RealPurePhaseInertGas {
    name: String,
    moles: f64,
    elements: Vec<f64>,
}

/// Reusable gas-side inventory for independent P,T/P,H validation.
#[derive(Debug, Clone)]
pub(crate) struct RealPurePhaseGasScenario {
    reactive_gas_moles: Vec<f64>,
    inerts: Vec<RealPurePhaseInertGas>,
}

impl RealPurePhaseGasScenario {
    pub(crate) fn new(reactive_gas_moles: Vec<f64>) -> Result<Self, ReactionExtentError> {
        if reactive_gas_moles.is_empty()
            || reactive_gas_moles
                .iter()
                .any(|moles| !moles.is_finite() || *moles <= 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_reactive_gas_moles",
                message: "reactive gas moles must be finite and strictly positive".into(),
            });
        }
        Ok(Self {
            reactive_gas_moles,
            inerts: Vec::new(),
        })
    }

    /// Appends one gas component which participates in no reaction but still
    /// contributes both to ideal-gas activities and elemental totals.
    pub(crate) fn with_inert(
        mut self,
        name: impl Into<String>,
        moles: f64,
        elements: Vec<f64>,
    ) -> Result<Self, ReactionExtentError> {
        let name = name.into();
        if name.trim().is_empty() || !moles.is_finite() || moles <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_inert",
                message: "inert name must be non-empty and moles finite and strictly positive"
                    .into(),
            });
        }
        if elements.is_empty()
            || elements
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_inert_elements",
                message: "inert elemental row must be finite and non-negative".into(),
            });
        }
        if self.inerts.iter().any(|inert| inert.name == name) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_inert",
                message: format!("inert gas '{name}' was declared more than once"),
            });
        }
        self.inerts.push(RealPurePhaseInertGas {
            name,
            moles,
            elements,
        });
        Ok(self)
    }

    /// Returns the gas ordering shared by independent P,T and P,H adapters:
    /// declared reactive components followed by zero-stoichiometry inerts.
    pub(crate) fn initial_gas_moles(&self) -> Vec<f64> {
        let mut moles = self.reactive_gas_moles.clone();
        moles.extend(self.inerts.iter().map(|inert| inert.moles));
        moles
    }
}

/// Complete physical inventory shared by canonical P,T, canonical P,H, and
/// the independent pure-phase validators.
///
/// One value object for every route prevents target construction and
/// cross-validation from silently using different total material amounts.
#[derive(Debug, Clone)]
pub(crate) struct RealPurePhaseInventory {
    gas: RealPurePhaseGasScenario,
    candidate_moles: f64,
}

/// Scenario data after phase-qualified identities have been matched to the
/// resolved layout. Keeping this private prevents positional component data
/// from leaking out of the shared fixture layer.
struct MaterializedGasScenario {
    names: Vec<String>,
    moles: Vec<f64>,
    stoichiometry: Vec<f64>,
    indices: Vec<usize>,
    element_matrix: DMatrix<f64>,
    candidate_index: usize,
}

impl RealPurePhaseInventory {
    pub(crate) fn new(
        reactive_gas_moles: Vec<f64>,
        candidate_moles: f64,
    ) -> Result<Self, ReactionExtentError> {
        Self::from_gas(
            RealPurePhaseGasScenario::new(reactive_gas_moles)?,
            candidate_moles,
        )
    }

    pub(crate) fn from_gas(
        gas: RealPurePhaseGasScenario,
        candidate_moles: f64,
    ) -> Result<Self, ReactionExtentError> {
        if !candidate_moles.is_finite() || candidate_moles < 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_candidate_moles",
                message: "candidate moles must be finite and non-negative".into(),
            });
        }
        Ok(Self {
            gas,
            candidate_moles,
        })
    }

    pub(crate) fn with_inert(
        mut self,
        name: impl Into<String>,
        moles: f64,
        elements: Vec<f64>,
    ) -> Result<Self, ReactionExtentError> {
        self.gas = self.gas.with_inert(name, moles, elements)?;
        Ok(self)
    }

    pub(crate) fn gas(&self) -> &RealPurePhaseGasScenario {
        &self.gas
    }
}

impl ResolvedRealPurePhaseFixture {
    pub(crate) fn family(&self) -> RealPurePhaseFamily {
        self.family
    }

    pub(crate) fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub(crate) fn thermochemistry(&self) -> &ResolvedThermochemistry {
        &self.thermochemistry
    }

    /// Builds a physical composition in the canonical resolved layout without
    /// assuming that phase names happen to sort in a particular order.
    pub(crate) fn initial_composition(
        &self,
        scenario: &RealPurePhaseInventory,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        let materialized = self.materialize_gas_scenario(scenario.gas())?;
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        let mut moles = vec![0.0; layout.component_count()];
        for (index, value) in materialized.indices.into_iter().zip(materialized.moles) {
            moles[index] = value;
        }
        moles[materialized.candidate_index] = scenario.candidate_moles;
        MultiphaseInitialComposition::from_dense(&layout, moles)
    }

    /// Rebinds an immutable production gas-boundary state to an already
    /// validated scenario declaration.
    ///
    /// Production evidence uses the complete canonical gas-phase ordering,
    /// while independent validators keep reactive species and zero-
    /// stoichiometry inerts separate. Matching names here prevents tests from
    /// silently translating that boundary through positional assumptions.
    pub(crate) fn gas_scenario_from_canonical_boundary(
        &self,
        template: &RealPurePhaseGasScenario,
        species: &[String],
        moles: &[f64],
    ) -> Result<RealPurePhaseGasScenario, ReactionExtentError> {
        let materialized = self.materialize_gas_scenario(template)?;
        if species != materialized.names {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_boundary_species",
                message: format!(
                    "canonical gas boundary {:?} does not match fixture ordering {:?}",
                    species, materialized.names
                ),
            });
        }
        if moles.len() != materialized.names.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "canonical gas boundary has {} amounts for {} components",
                moles.len(),
                materialized.names.len()
            )));
        }
        if moles
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_boundary_moles",
                message: "canonical gas boundary moles must be finite and strictly positive".into(),
            });
        }

        let reactive_count = self.family.definition().gas_species.len();
        let mut rebound = template.clone();
        rebound
            .reactive_gas_moles
            .copy_from_slice(&moles[..reactive_count]);
        for (inert, value) in rebound.inerts.iter_mut().zip(&moles[reactive_count..]) {
            inert.moles = *value;
        }
        Ok(rebound)
    }

    /// Resolves a scenario's explicit gas ordering against this fixture's
    /// phase-qualified layout and materializes its independent element matrix.
    fn materialize_gas_scenario(
        &self,
        scenario: &RealPurePhaseGasScenario,
    ) -> Result<MaterializedGasScenario, ReactionExtentError> {
        let definition = self.family.definition();
        if scenario.reactive_gas_moles.len() != definition.gas_species.len() {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "{} expects {} reactive gas moles, got {}",
                definition.name,
                definition.gas_species.len(),
                scenario.reactive_gas_moles.len()
            )));
        }
        if scenario
            .inerts
            .iter()
            .any(|inert| inert.elements.len() != definition.elements.len())
        {
            return Err(ReactionExtentError::DimensionMismatch(format!(
                "{} inert elemental rows must have {} entries",
                definition.name,
                definition.elements.len()
            )));
        }
        let labels = self.resolved.layout().component_labels();
        let index_for = |phase: &str, species: &str| {
            let label = format!("{phase}::{species}");
            labels
                .iter()
                .position(|current| current == &label)
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "real_pure_phase_fixture_layout",
                    message: format!("{} lacks resolved component '{label}'", definition.name),
                })
        };
        let mut names = definition
            .gas_species
            .iter()
            .map(|species| (*species).to_string())
            .collect::<Vec<_>>();
        names.extend(scenario.inerts.iter().map(|inert| inert.name.clone()));
        let indices = names
            .iter()
            .map(|species| index_for(definition.gas_phase, species))
            .collect::<Result<Vec<_>, _>>()?;
        let candidate_index = index_for(definition.candidate_phase, definition.candidate_species)?;
        let mut stoichiometry = definition.gas_stoichiometry.to_vec();
        stoichiometry.extend(std::iter::repeat(0.0).take(scenario.inerts.len()));
        let base_matrix = self.family.gas_element_matrix();
        let element_matrix =
            DMatrix::from_fn(names.len(), definition.elements.len(), |row, column| {
                if row < definition.gas_species.len() {
                    base_matrix[(row, column)]
                } else {
                    scenario.inerts[row - definition.gas_species.len()].elements[column]
                }
            });
        Ok(MaterializedGasScenario {
            names,
            moles: scenario.initial_gas_moles(),
            stoichiometry,
            indices,
            element_matrix,
            candidate_index,
        })
    }

    /// Builds the independent P,T boundary problem from the exact declared
    /// chemistry. Closures are intentionally isolated behind the legacy
    /// infallible boundary API; the common fixture preflight has already
    /// established local G(T) capability over the shared temperature range.
    pub(crate) fn to_pt_boundary_problem(
        &self,
        scenario: &RealPurePhaseGasScenario,
        conditions: EquilibriumConditions,
    ) -> Result<PurePhaseBoundaryProblem, ReactionExtentError> {
        let materialized = self.materialize_gas_scenario(scenario)?;
        let definition = self.family.definition();
        let thermochemistry = Arc::new(self.thermochemistry.clone());
        let gibbs = materialized
            .indices
            .iter()
            .map(|&index| indexed_boundary_gibbs(Arc::clone(&thermochemistry), index))
            .collect();
        let problem = PurePhaseBoundaryProblem::new(
            materialized.names,
            materialized.moles,
            materialized.stoichiometry,
            definition.candidate_stoichiometry,
            gibbs,
            indexed_boundary_gibbs(thermochemistry, materialized.candidate_index),
            conditions,
            definition.candidate_species,
        )?;
        problem.with_element_composition(
            PurePhaseBoundaryElementComposition::new(
                definition
                    .elements
                    .iter()
                    .map(|name| (*name).to_string())
                    .collect(),
                materialized.element_matrix,
                definition.candidate_elements.to_vec(),
            )?,
            PurePhaseBoundaryStructuralTolerances::default(),
        )
    }

    /// Builds the independent P,H scalar problem for the exact declared
    /// chemistry, including explicitly declared zero-stoichiometry inerts.
    pub(crate) fn to_ph_problem(
        &self,
        scenario: &RealPurePhaseInventory,
        pressure: f64,
        reference_pressure: f64,
        target_enthalpy: f64,
        temperature_bounds: TemperatureBounds,
    ) -> Result<PurePhasePhProblem, ReactionExtentError> {
        let definition = self.family.definition();
        let materialized = self.materialize_gas_scenario(scenario.gas())?;
        let thermochemistry = Arc::new(self.thermochemistry.clone());
        let gibbs = materialized
            .indices
            .iter()
            .map(|&index| indexed_gibbs(Arc::clone(&thermochemistry), index))
            .collect();
        let enthalpy = materialized
            .indices
            .iter()
            .map(|&index| indexed_enthalpy(Arc::clone(&thermochemistry), index))
            .collect();
        let problem = PurePhasePhProblem::new(
            materialized.names,
            materialized.moles,
            materialized.stoichiometry,
            scenario.candidate_moles,
            definition.candidate_stoichiometry,
            definition.candidate_species,
            gibbs,
            indexed_gibbs(Arc::clone(&thermochemistry), materialized.candidate_index),
            enthalpy,
            indexed_enthalpy(thermochemistry, materialized.candidate_index),
            PurePhasePhConditions::new(
                pressure,
                reference_pressure,
                target_enthalpy,
                temperature_bounds.lower(),
                temperature_bounds.upper(),
            )?,
        )?;
        problem.with_element_composition(
            PurePhaseBoundaryElementComposition::new(
                definition
                    .elements
                    .iter()
                    .map(|name| (*name).to_string())
                    .collect(),
                materialized.element_matrix,
                definition.candidate_elements.to_vec(),
            )?,
            PurePhaseBoundaryStructuralTolerances::default(),
        )
    }
}

fn indexed_boundary_gibbs(thermochemistry: Arc<ResolvedThermochemistry>, index: usize) -> GibbsFn {
    Rc::new(move |temperature| {
        thermochemistry
            .evaluate_gibbs(temperature)
            .ok()
            .and_then(|values| values.get(index).copied())
            .unwrap_or(f64::NAN)
    })
}

fn indexed_gibbs(
    thermochemistry: Arc<ResolvedThermochemistry>,
    index: usize,
) -> MolarThermoFunction {
    Arc::new(move |temperature| {
        thermochemistry
            .evaluate_gibbs(temperature)?
            .get(index)
            .copied()
            .ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "real pure-phase Gibbs bundle has no component {index}"
                ))
            })
    })
}

fn indexed_enthalpy(
    thermochemistry: Arc<ResolvedThermochemistry>,
    index: usize,
) -> MolarThermoFunction {
    Arc::new(move |temperature| {
        thermochemistry
            .evaluate_enthalpy(temperature)?
            .get(index)
            .copied()
            .ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "real pure-phase enthalpy bundle has no component {index}"
                ))
            })
    })
}

impl RealPurePhaseFamily {
    /// Returns the declarative chemistry contract in explicit component order.
    /// Thermochemical record selection remains a separate fallible resolution
    /// step, so this declaration cannot silently depend on library iteration.
    pub(crate) fn definition(self) -> RealPurePhaseDefinition {
        match self {
            Self::WaterLiquid => RealPurePhaseDefinition {
                name: "water-liquid",
                gas_phase: "gas",
                candidate_phase: "liquid",
                gas_species: &["H2O"],
                candidate_species: "H2O",
                gas_stoichiometry: &[-1.0],
                candidate_stoichiometry: 1.0,
                elements: &["H", "O"],
                gas_element_rows: &[2.0, 1.0],
                candidate_elements: &[2.0, 1.0],
                candidate_state: PhysicalState::Liquid,
                gas_library: "NASA_gas",
                candidate_library: "NASA_cond",
            },
            Self::WaterIce => RealPurePhaseDefinition {
                name: "water-ice",
                gas_phase: "gas",
                candidate_phase: "solid",
                gas_species: &["H2O"],
                candidate_species: "H2O(s)",
                gas_stoichiometry: &[-1.0],
                candidate_stoichiometry: 1.0,
                elements: &["H", "O"],
                gas_element_rows: &[2.0, 1.0],
                candidate_elements: &[2.0, 1.0],
                candidate_state: PhysicalState::Solid,
                gas_library: "NASA_gas",
                candidate_library: "NASA_cond",
            },
            Self::BoudouardCarbon => RealPurePhaseDefinition {
                name: "boudouard-carbon",
                gas_phase: "gas",
                candidate_phase: "solid",
                gas_species: &["CO", "CO2"],
                candidate_species: "C(gr)",
                gas_stoichiometry: &[-2.0, 1.0],
                candidate_stoichiometry: 1.0,
                elements: &["C", "O"],
                gas_element_rows: &[1.0, 1.0, 1.0, 2.0],
                candidate_elements: &[1.0, 0.0],
                candidate_state: PhysicalState::Solid,
                gas_library: "NASA_gas",
                candidate_library: "NASA_cond",
            },
            Self::MethaneCarbon => RealPurePhaseDefinition {
                name: "methane-carbon",
                gas_phase: "gas",
                candidate_phase: "solid",
                gas_species: &["CH4", "H2"],
                candidate_species: "C(gr)",
                gas_stoichiometry: &[-1.0, 2.0],
                candidate_stoichiometry: 1.0,
                elements: &["C", "H"],
                gas_element_rows: &[1.0, 4.0, 0.0, 2.0],
                candidate_elements: &[1.0, 0.0],
                candidate_state: PhysicalState::Solid,
                gas_library: "NASA_gas",
                candidate_library: "NASA_cond",
            },
        }
    }

    /// Builds a state-constrained offline lookup specification. The fixture
    /// layer is intentionally explicit about libraries and disables NIST;
    /// inventory tests decide whether a local record actually satisfies it.
    pub(crate) fn offline_spec(self) -> SubstanceSystemSpec {
        self.offline_spec_with_inert(&[])
    }

    /// Builds the same restricted family with explicitly declared inert gas
    /// components. Inerts are lookup/scenario inputs with zero reaction
    /// coefficient; they never change the family reaction contract.
    pub(crate) fn offline_spec_with_inert(self, inert_gases: &[&str]) -> SubstanceSystemSpec {
        let definition = self.definition();
        let mut gas_species = definition
            .gas_species
            .iter()
            .map(|name| (*name).to_string())
            .collect::<Vec<_>>();
        gas_species.extend(inert_gases.iter().map(|name| (*name).to_string()));
        SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
            (definition.gas_phase.to_string(), gas_species),
            (
                definition.candidate_phase.to_string(),
                vec![definition.candidate_species.to_string()],
            ),
        ])))
        .with_phase_natures(Some(HashMap::from([
            (definition.gas_phase.to_string(), PhysicalState::Gas),
            (
                definition.candidate_phase.to_string(),
                definition.candidate_state,
            ),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("static real pure-phase declaration must be structurally valid")
    }

    /// Returns the independent elemental matrix in the declaration's explicit
    /// gas-species ordering. Construction code must validate conservation; the
    /// definition only prevents positional ambiguity.
    pub(crate) fn gas_element_matrix(self) -> DMatrix<f64> {
        let definition = self.definition();
        DMatrix::from_row_slice(
            definition.gas_species.len(),
            definition.elements.len(),
            definition.gas_element_rows,
        )
    }

    /// Verifies that ordinary resolution preserved this fixture's explicit
    /// state and local-library contract. This is intentionally based on the
    /// typed phase specs and lookup report, not filename conventions or an
    /// assumed JSON iteration order.
    fn validate_resolved_contract(
        self,
        resolved: &ResolvedPhaseSystem,
        inert_gases: &[&str],
    ) -> Result<(), ReactionExtentError> {
        let definition = self.definition();
        let phase = |name: &str| {
            resolved
                .phase_specs()
                .iter()
                .find(|spec| {
                    spec.id()
                        .as_option()
                        .as_ref()
                        .is_some_and(|id| id.as_str() == name)
                })
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "real_pure_phase_fixture_phase",
                    message: format!("{} lacks declared phase '{name}'", definition.name),
                })
        };
        let gas = phase(definition.gas_phase)?;
        if gas.physical_state() != PhysicalState::Gas || gas.model() != PhaseModel::IdealGas {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_fixture_gas_state",
                message: format!("{} did not resolve an ideal gas phase", definition.name),
            });
        }
        let mut expected_gas = definition
            .gas_species
            .iter()
            .map(|name| (*name).to_string())
            .collect::<Vec<_>>();
        expected_gas.extend(inert_gases.iter().map(|name| (*name).to_string()));
        if gas.components() != expected_gas {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_fixture_gas_components",
                message: format!(
                    "{} resolved gas components {:?}, expected {:?}",
                    definition.name,
                    gas.components(),
                    expected_gas
                ),
            });
        }
        let candidate = phase(definition.candidate_phase)?;
        if candidate.physical_state() != definition.candidate_state
            || candidate.model() != PhaseModel::PureCondensed
            || candidate.components() != [definition.candidate_species]
        {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_fixture_candidate_state",
                message: format!(
                    "{} did not resolve the requested pure {:?} candidate '{}'",
                    definition.name, definition.candidate_state, definition.candidate_species
                ),
            });
        }
        let report_for = |phase_name: &str| {
            resolved
                .report()
                .phases()
                .iter()
                .find(|summary| {
                    summary
                        .phase()
                        .as_option()
                        .as_ref()
                        .is_some_and(|id| id.as_str() == phase_name)
                })
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "real_pure_phase_fixture_provenance",
                    message: format!("{} lacks lookup report for '{phase_name}'", definition.name),
                })
        };
        let require_thermo_rows = |phase_name: &str, species: &[String], library: &str| {
            let report = report_for(phase_name)?;
            for substance in species {
                let selected = report.search().rows().iter().any(|row| {
                    row.property() == "Thermo"
                        && row.substance() == substance
                        && row.library() == library
                        && !row.record_key().trim().is_empty()
                });
                if !selected {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "real_pure_phase_fixture_provenance",
                        message: format!(
                            "{} lacks local Thermo provenance for {} in {}",
                            definition.name, substance, library
                        ),
                    });
                }
            }
            Ok(())
        };
        require_thermo_rows(definition.gas_phase, &expected_gas, definition.gas_library)?;
        require_thermo_rows(
            definition.candidate_phase,
            &[definition.candidate_species.to_string()],
            definition.candidate_library,
        )
    }

    /// Resolves one restricted family through the ordinary local phase API.
    /// It validates the common G/H/Cp capability boundary once, then callers
    /// may safely reuse its immutable payload across many scenarios.
    pub(crate) fn resolve_offline(
        self,
        repository: std::sync::Arc<ThermoRepository>,
        inert_gases: &[&str],
    ) -> Result<ResolvedRealPurePhaseFixture, ReactionExtentError> {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            self.offline_spec_with_inert(inert_gases),
            repository,
        )
        .map_err(|error| ReactionExtentError::InvalidProblem {
            field: "real_pure_phase_fixture_resolution",
            message: format!(
                "{} could not resolve offline: {error}",
                self.definition().name
            ),
        })?;
        if resolved.report().nist_fallback_enabled() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_fixture_nist",
                message: format!(
                    "{} unexpectedly enabled NIST fallback",
                    self.definition().name
                ),
            });
        }
        self.validate_resolved_contract(&resolved, inert_gases)?;
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        let temperature = thermochemistry.temperature_bounds().lower();
        thermochemistry.evaluate_gibbs(temperature)?;
        thermochemistry.evaluate_enthalpy(temperature)?;
        if thermochemistry
            .evaluate_heat_capacity(temperature)?
            .iter()
            .any(Option::is_none)
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "real_pure_phase_fixture_capabilities",
                message: format!(
                    "{} has no complete local Cp(T) capability",
                    self.definition().name
                ),
            });
        }
        Ok(ResolvedRealPurePhaseFixture {
            family: self,
            resolved,
            thermochemistry,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::sync::Arc;

    use super::{RealPurePhaseFamily, RealPurePhaseGasScenario, RealPurePhaseInventory};
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::TemperatureBounds;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
        PurePhaseBoundaryStructuralTolerances, PurePhaseBoundaryTolerances,
        evaluate_pure_phase_boundary,
    };
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("bundled offline repository must be available for fixture inventory")
    }

    /// Fixture resolution is read-only. Compare the complete local library
    /// payload byte-for-byte so a regression cannot hide a write API behind a
    /// convenient resolver or adapter constructor.
    fn local_library_snapshot() -> Vec<(String, Vec<u8>)> {
        let paths = with_library_manager(|manager| {
            vec![
                manager.substance_base_path().to_string(),
                manager.all_keys_substance_path().to_string(),
                manager.elements_path().to_string(),
            ]
        });
        paths
            .into_iter()
            .map(|path| {
                let bytes = fs::read(&path).unwrap_or_else(|error| {
                    panic!("must read local fixture library '{path}': {error}")
                });
                (path, bytes)
            })
            .collect()
    }

    #[test]
    fn all_declared_families_conserve_their_explicit_elements() {
        for family in [
            RealPurePhaseFamily::WaterLiquid,
            RealPurePhaseFamily::WaterIce,
            RealPurePhaseFamily::BoudouardCarbon,
            RealPurePhaseFamily::MethaneCarbon,
        ] {
            let definition = family.definition();
            let _offline_spec = family.offline_spec();
            let matrix = family.gas_element_matrix();
            for element in 0..definition.elements.len() {
                let gas = definition
                    .gas_stoichiometry
                    .iter()
                    .enumerate()
                    .map(|(species, coefficient)| coefficient * matrix[(species, element)])
                    .sum::<f64>();
                let total = gas
                    + definition.candidate_stoichiometry * definition.candidate_elements[element];
                assert!(
                    total.abs() <= 1e-12,
                    "{} does not conserve {}: {total:e}",
                    definition.name,
                    definition.elements[element]
                );
            }
        }
    }

    #[test]
    fn water_liquid_pt_and_ph_adapters_share_the_same_inert_gas_layout() {
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &["O2"])
            .expect("local water/O2 family must resolve without NIST");
        let gas = RealPurePhaseGasScenario::new(vec![0.5])
            .and_then(|scenario| scenario.with_inert("O2", 0.25, vec![0.0, 2.0]))
            .expect("water/O2 gas scenario must validate");
        let conditions = EquilibriumConditions::new(350.0, 101_325.0, 101_325.0)
            .expect("water P,T conditions must validate");
        let pt = fixture
            .to_pt_boundary_problem(&gas, conditions)
            .expect("water P,T adapter must materialize an independent problem");
        let reaction_space = pt
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect("water/O2 family must retain one strict phase-forming direction");
        assert_eq!(reaction_space.full_reaction_dimension, 1);
        assert_eq!(reaction_space.gas_only_reaction_dimension, 0);
        assert_eq!(pt.gas_species(), ["H2O", "O2"]);
        let report = evaluate_pure_phase_boundary(&pt, PurePhaseBoundaryTolerances::default())
            .expect("water P,T boundary report must remain finite");
        assert!(report.log_residual_at_absence.is_finite());

        let rebound = fixture
            .gas_scenario_from_canonical_boundary(
                &gas,
                &["H2O".to_string(), "O2".to_string()],
                &[0.4, 0.35],
            )
            .expect("canonical gas boundary must preserve reactive/inert metadata");
        assert_eq!(rebound.initial_gas_moles(), vec![0.4, 0.35]);
        assert!(matches!(
            fixture.gas_scenario_from_canonical_boundary(
                &gas,
                &["O2".to_string(), "H2O".to_string()],
                &[0.35, 0.4],
            ),
            Err(ReactionExtentError::InvalidProblem {
                field: "real_pure_phase_boundary_species",
                ..
            })
        ));

        let ph_scenario = RealPurePhaseInventory::from_gas(gas, 0.1)
            .expect("water P,H candidate inventory must validate");
        let ph = fixture
            .to_ph_problem(
                &ph_scenario,
                101_325.0,
                101_325.0,
                -161_377.821_776_069_12,
                TemperatureBounds::new(345.0, 355.0)
                    .expect("water P,H temperature bounds must validate"),
            )
            .expect("water P,H adapter must materialize an independent problem");
        assert_eq!(ph.gas_species(), pt.gas_species());
        assert_eq!(ph.gas_stoichiometry(), pt.gas_stoichiometry());
        assert_eq!(ph.candidate_name(), pt.candidate_name());
        let ph_reaction_space = ph
            .reaction_space(PurePhaseBoundaryStructuralTolerances::default())
            .expect("water P,H adapter must retain the strict reaction space");
        assert_eq!(ph_reaction_space, reaction_space);
    }

    #[test]
    fn every_inventory_family_materializes_strict_pt_and_ph_problems() {
        let before = local_library_snapshot();
        let cases = [
            (RealPurePhaseFamily::WaterLiquid, vec![0.5], 350.0),
            (RealPurePhaseFamily::WaterIce, vec![0.5], 250.0),
            (RealPurePhaseFamily::BoudouardCarbon, vec![0.5, 0.25], 700.0),
            (RealPurePhaseFamily::MethaneCarbon, vec![0.5, 0.25], 1000.0),
        ];

        for (family, reactive_gas_moles, temperature) in cases {
            let fixture = family
                .resolve_offline(local_repository(), &[])
                .unwrap_or_else(|error| {
                    panic!("{} must resolve: {error}", family.definition().name)
                });
            let gas = RealPurePhaseGasScenario::new(reactive_gas_moles)
                .expect("declared real gas inventory must validate");
            let conditions = EquilibriumConditions::new(temperature, 101_325.0, 101_325.0)
                .expect("inventory P,T conditions must validate");
            let pt = fixture
                .to_pt_boundary_problem(&gas, conditions)
                .unwrap_or_else(|error| {
                    panic!(
                        "{} P,T adapter must materialize: {error}",
                        family.definition().name
                    )
                });
            let pt_space =
                pt.validate_strict_independent_family(
                    PurePhaseBoundaryStructuralTolerances::default(),
                )
                .unwrap_or_else(|error| {
                    panic!(
                        "{} P,T structure must stay strict: {error}",
                        family.definition().name
                    )
                });
            assert_eq!(pt_space.full_reaction_dimension, 1);
            assert_eq!(pt_space.gas_only_reaction_dimension, 0);
            let pt_report =
                evaluate_pure_phase_boundary(&pt, PurePhaseBoundaryTolerances::default())
                    .unwrap_or_else(|error| {
                        panic!(
                            "{} P,T boundary must evaluate: {error}",
                            family.definition().name
                        )
                    });
            assert!(pt_report.log_residual_at_absence.is_finite());

            let ph = fixture
                .to_ph_problem(
                    &RealPurePhaseInventory::from_gas(gas, 0.1)
                        .expect("declared real P,H inventory must validate"),
                    101_325.0,
                    101_325.0,
                    0.0,
                    TemperatureBounds::new(temperature - 1.0, temperature + 1.0)
                        .expect("local P,H range must validate"),
                )
                .unwrap_or_else(|error| {
                    panic!(
                        "{} P,H adapter must materialize: {error}",
                        family.definition().name
                    )
                });
            let ph_space = ph
                .reaction_space(PurePhaseBoundaryStructuralTolerances::default())
                .unwrap_or_else(|error| {
                    panic!(
                        "{} P,H structure must stay strict: {error}",
                        family.definition().name
                    )
                });
            assert_eq!(ph_space, pt_space);
        }
        assert_eq!(before, local_library_snapshot());
    }

    /// Diagnostic inventory, deliberately ignored because it prints local
    /// catalog availability rather than enforcing that every future family is
    /// already supported. A missing capability is a typed data fact used to
    /// choose the next real validation scenario.
    #[test]
    #[ignore = "diagnostic offline real pure-phase inventory"]
    fn offline_real_pure_phase_inventory() {
        for family in [
            RealPurePhaseFamily::WaterLiquid,
            RealPurePhaseFamily::WaterIce,
            RealPurePhaseFamily::BoudouardCarbon,
            RealPurePhaseFamily::MethaneCarbon,
        ] {
            let definition = family.definition();
            match family.resolve_offline(local_repository(), &[]) {
                Ok(fixture) => {
                    let bounds = fixture.thermochemistry().temperature_bounds();
                    println!(
                        "{:18} available  T=[{:.2}, {:.2}] K  components={:?}",
                        definition.name,
                        bounds.lower(),
                        bounds.upper(),
                        fixture.resolved().layout().component_labels(),
                    );
                }
                Err(error) => println!("{:18} unavailable  {error}", definition.name),
            }
        }
    }
}
