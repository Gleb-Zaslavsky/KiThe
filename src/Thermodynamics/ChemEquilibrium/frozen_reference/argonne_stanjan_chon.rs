//! First general multicomponent I5 equilibrium characterization: Argonne/STANJAN CHON.
//!
//! This module deliberately owns only frozen-case adaptation: exact local
//! identities, an element-equivalent feed, preflight, and source-facing
//! comparison rows. The actual calculation remains the ordinary production
//! fixed-`P,T` log-mole Gibbs solve.

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use nalgebra::{DMatrix, linalg::SVD};

use crate::Kinetics::molmass::parse_formula;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::ResolvedThermochemistry;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    ArgonneStanjanChonPtReference, FrozenReferenceDataset,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::DBhandlers::thermo_api::ThermoCalculator;
use crate::Thermodynamics::User_PhaseOrSolution::element_composition_and_molar_mass;
use crate::Thermodynamics::User_PhaseOrSolution::{
    ResolvedPhaseSystem, SubstanceSystemFactory, SubstanceSystemSpec, SubstanceSystemSpecBuilder,
    SubstancesContainer,
};
use crate::Thermodynamics::User_substances::{CalculatorType, WhatIsFound};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// The sole ideal-gas phase in the published STANJAN comparison universe.
pub(crate) const STANJAN_GAS_PHASE: &str = "gas";
/// Standard pressure used by the local NASA thermochemistry activity model.
pub(crate) const STANJAN_REFERENCE_PRESSURE_PA: f64 = 101_325.0;

/// Full source-table order. `C5H12` is retained as external evidence even
/// though the locally executable element-equivalent universe does not include
/// a pentane record.
pub(crate) const STANJAN_DECLARED_SPECIES: [&str; 16] = [
    "C5H12", "CH4", "O2", "CO2", "H2O", "N2", "N", "O", "NO", "OH", "H", "N2O", "CO", "H2", "NO2",
    "HO2",
];

/// Exact local universe used by the production solve. Do not auto-expand it:
/// agreement is meaningful only against the same declared reaction space.
pub(crate) const STANJAN_LOCAL_GAS_SPECIES: [&str; 15] = [
    "CH4", "O2", "CO2", "H2O", "N2", "N", "O", "NO", "OH", "H", "N2O", "CO", "H2", "NO2", "HO2",
];

/// Local record interval and provenance used by the CHON preflight table.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct StanjanChonPreflightRow {
    /// STANJAN source identity (species formula) as published in the table.
    pub(crate) external_identity: String,
    /// Exact local gas-phase component matched by identity.
    pub(crate) component: PhaseComponentId,
    /// Library that supplied the offline thermochemistry record.
    pub(crate) library: String,
    /// Repository record key backing the local coefficients.
    pub(crate) record_key: String,
    /// Physical state required by the frozen case (always ideal gas).
    pub(crate) physical_state: PhysicalState,
    /// Lower bound of the record's valid temperature interval, K.
    pub(crate) temperature_lower_k: f64,
    /// Upper bound of the record's valid temperature interval, K.
    pub(crate) temperature_upper_k: f64,
    /// Whether the frozen target temperature lies inside this record interval.
    pub(crate) supports_target_temperature: bool,
    /// Molar Gibbs energy at the frozen temperature, J/mol.
    pub(crate) gibbs_j_mol: f64,
}

/// Rank evidence proving this is a multidimensional general equilibrium case.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct StanjanChonStructure {
    /// Number of species in the local reaction universe.
    pub(crate) component_count: usize,
    /// Number of distinct chemical elements present.
    pub(crate) element_count: usize,
    /// Rank of the element-composition matrix.
    pub(crate) element_rank: usize,
    /// Dimension of the reaction space (`component_count - element_rank`).
    pub(crate) reaction_dimension: usize,
}

/// Immutable resolved local side of the Argonne/STANJAN benchmark.
#[derive(Clone)]
pub(crate) struct ResolvedArgonneStanjanChonFixture {
    resolved: ResolvedPhaseSystem,
    preflight: Vec<StanjanChonPreflightRow>,
    structure: StanjanChonStructure,
}

/// Semantic classification for one published STANJAN component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum StanjanChonMagnitudeClass {
    /// Mole fraction at or above `1e-3`; errors reported on the relative scale.
    Major,
    /// Mole fraction between `1e-6` and `1e-3`; errors reported as `delta_log10`.
    Minor,
    /// Mole fraction below `1e-6`; errors reported as `delta_log10`.
    Trace,
    /// Published zero reactant retained only as external evidence, not solved locally.
    ExternallyZeroReactantNotSolved,
}

/// Identity-aligned local/external composition row.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct StanjanChonComparisonRow {
    /// STANJAN source identity (species formula) for this row.
    pub(crate) external_identity: String,
    /// Matched local gas-phase component, or `None` when excluded.
    pub(crate) local_component: Option<PhaseComponentId>,
    /// Semantic magnitude classification selecting the error scale.
    pub(crate) magnitude_class: StanjanChonMagnitudeClass,
    /// Mole fraction published by STANJAN.
    pub(crate) stanjan_mole_fraction: f64,
    /// Local mole fraction, or `None` when excluded from the solve.
    pub(crate) kithe_mole_fraction: Option<f64>,
    /// Absolute error `local - stanjan`, or `None` when excluded.
    pub(crate) absolute_error: Option<f64>,
    /// Relative error, or `None` when the STANJAN value is not positive.
    pub(crate) relative_error: Option<f64>,
    /// Log10 difference (`log10(local) - log10(stanjan)`), or `None` when not computable.
    pub(crate) delta_log10: Option<f64>,
    /// Human-readable reason for exclusion, or `None` when included.
    pub(crate) exclusion_reason: Option<String>,
}

/// Aggregate diagnostics for one non-zero external magnitude class.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct StanjanChonErrorSummary {
    /// Number of sampled species contributing to the aggregate.
    pub(crate) sample_count: usize,
    /// Largest absolute value among the sampled errors.
    pub(crate) maximum_absolute: f64,
    /// Root-mean-square of the sampled errors.
    pub(crate) rms: f64,
    /// Signed mean of the sampled errors (retains bias direction).
    pub(crate) mean_signed: f64,
}

/// Characterization-only comparison of the general local solver with STANJAN.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct StanjanChonComparisonReport {
    /// Identifier of the frozen STANJAN dataset that produced this report.
    pub(crate) dataset_id: String,
    /// Physical (total) pressure of the frozen case, Pa.
    pub(crate) physical_pressure_pa: f64,
    /// Standard pressure used by the local activity model, Pa.
    pub(crate) activity_reference_pressure_pa: f64,
    /// Frozen equilibrium temperature, K.
    pub(crate) temperature_k: f64,
    /// Sum of the published STANJAN mole fractions.
    pub(crate) external_sum: f64,
    /// Sum of the local (KiThe) mole fractions across non-excluded rows.
    pub(crate) kithe_sum: f64,
    /// Relative-error summary for `Major` species, when any are sampled.
    pub(crate) major_relative: Option<StanjanChonErrorSummary>,
    /// Log10-error summary for `Minor` species, when any are sampled.
    pub(crate) minor_log10: Option<StanjanChonErrorSummary>,
    /// Log10-error summary for `Trace` species, when any are sampled.
    pub(crate) trace_log10: Option<StanjanChonErrorSummary>,
    /// Identity-aligned comparison rows for every published species.
    pub(crate) species: Vec<StanjanChonComparisonRow>,
    /// Dimension of the local reaction space (11 for this case).
    pub(crate) reaction_dimension: usize,
    /// Number of phase-control transitions recorded by the production solve.
    pub(crate) phase_control_transitions: usize,
    /// Residual L2 norm reported by the accepted solution's validation.
    pub(crate) residual_l2_norm: f64,
    /// Maximum absolute element-balance error from the accepted solution.
    pub(crate) max_abs_element_balance_error: f64,
    /// Debug name of the backend that produced the accepted solution.
    pub(crate) accepted_backend: String,
}

impl ResolvedArgonneStanjanChonFixture {
    /// Resolves precisely the reviewed offline 15-species NASA-gas universe.
    pub(crate) fn resolve_offline(
        repository: Arc<ThermoRepository>,
        reference: &ArgonneStanjanChonPtReference,
    ) -> Result<Self, ReactionExtentError> {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            stanjan_gas_spec(),
            repository,
        )
        .map_err(|error| ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_resolution",
            message: format!(
                "exact 15-species CHON gas universe could not resolve offline: {error}"
            ),
        })?;
        if resolved.report().nist_fallback_enabled() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "argonne_stanjan_chon_nist",
                message: "frozen STANJAN benchmark must not enable NIST fallback".to_owned(),
            });
        }
        validate_exact_local_universe(&resolved)?;
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        if !thermochemistry
            .temperature_bounds()
            .contains(reference.temperature_k)
        {
            let bounds = thermochemistry.temperature_bounds();
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "argonne_stanjan_chon_temperature_domain",
                message: format!(
                    "common local interval [{}, {}] K excludes {} K",
                    bounds.lower(),
                    bounds.upper(),
                    reference.temperature_k
                ),
            });
        }
        let (element_matrix, _, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                    path: "argonne_stanjan_chon_elements",
                    message: format!("local component composition is unavailable: {error}"),
                })?;
        let structure = structure_from_matrix(&element_matrix)?;
        if structure
            != (StanjanChonStructure {
                component_count: 15,
                element_count: 4,
                element_rank: 4,
                reaction_dimension: 11,
            })
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "argonne_stanjan_chon_structure",
                message: format!(
                    "expected 15 components / 4 elements / rank 4 / 11 reactions, got {}/{}/{}/{}",
                    structure.component_count,
                    structure.element_count,
                    structure.element_rank,
                    structure.reaction_dimension
                ),
            });
        }
        let preflight = preflight_rows(&resolved, &thermochemistry, reference.temperature_k)?;
        if let Some(row) = preflight
            .iter()
            .find(|row| !row.supports_target_temperature)
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "argonne_stanjan_chon_component_temperature_domain",
                message: format!(
                    "{} local record interval [{}, {}] K excludes {} K",
                    row.external_identity,
                    row.temperature_lower_k,
                    row.temperature_upper_k,
                    reference.temperature_k
                ),
            });
        }
        Ok(Self {
            resolved,
            preflight,
            structure,
        })
    }

    pub(crate) fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub(crate) fn preflight(&self) -> &[StanjanChonPreflightRow] {
        &self.preflight
    }

    pub(crate) fn structure(&self) -> StanjanChonStructure {
        self.structure
    }

    /// Builds the physical element-equivalent input used by the local solve.
    /// Zero product amounts stay absent from the physical inventory; numerical
    /// log-space trace coordinates are wholly owned by the solver policy.
    pub(crate) fn initial_composition(
        &self,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        MultiphaseInitialComposition::from_sparse(
            &layout,
            vec![
                (gas_component("CH4"), 4.0),
                (gas_component("CO2"), 2.0),
                (gas_component("O2"), 8.0),
                (gas_component("N2"), 37.6),
            ],
        )
    }
}

impl StanjanChonComparisonReport {
    /// Builds an identity-based report from an accepted production `P,T` state.
    pub(crate) fn from_solution(
        dataset: &FrozenReferenceDataset<ArgonneStanjanChonPtReference>,
        fixture: &ResolvedArgonneStanjanChonFixture,
        solution: &MultiphaseEquilibriumSolution,
    ) -> Result<Self, ReactionExtentError> {
        let reference =
            dataset
                .rows()
                .first()
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "argonne_stanjan_chon_reference",
                    message: "frozen STANJAN dataset has no case".to_owned(),
                })?;
        let species = reference
            .species_mole_fractions
            .iter()
            .map(|external| {
                if external.species == "C5H12" {
                    return Ok(StanjanChonComparisonRow {
                        external_identity: external.species.clone(),
                        local_component: None,
                        magnitude_class: StanjanChonMagnitudeClass::ExternallyZeroReactantNotSolved,
                        stanjan_mole_fraction: external.mole_fraction,
                        kithe_mole_fraction: None,
                        absolute_error: None,
                        relative_error: None,
                        delta_log10: None,
                        exclusion_reason: Some(
                            "published zero reactant retained as external evidence; the local solve uses an element-equivalent 15-species universe without C5H12".to_owned(),
                        ),
                    });
                }
                let component = gas_component(&external.species);
                let local = solution.mole_fraction_for(&component).ok_or_else(|| {
                    ReactionExtentError::DimensionMismatch(format!(
                        "accepted local layout lacks STANJAN component {}",
                        external.species
                    ))
                })?;
                let absolute_error = local - external.mole_fraction;
                let relative_error = (external.mole_fraction > 0.0)
                    .then_some(absolute_error / external.mole_fraction);
                let delta_log10 = (local > 0.0 && external.mole_fraction > 0.0)
                    .then_some(local.log10() - external.mole_fraction.log10());
                Ok(StanjanChonComparisonRow {
                    external_identity: external.species.clone(),
                    local_component: Some(component),
                    magnitude_class: magnitude_class(external.mole_fraction),
                    stanjan_mole_fraction: external.mole_fraction,
                    kithe_mole_fraction: Some(local),
                    absolute_error: Some(absolute_error),
                    relative_error,
                    delta_log10,
                    exclusion_reason: None,
                })
            })
            .collect::<Result<Vec<_>, ReactionExtentError>>()?;
        let kithe_sum = species
            .iter()
            .filter_map(|row| row.kithe_mole_fraction)
            .sum::<f64>();
        let validation = solution.accepted_solution().validation();
        Ok(Self {
            dataset_id: dataset.metadata().dataset_id.clone(),
            physical_pressure_pa: reference.pressure_pa,
            activity_reference_pressure_pa: STANJAN_REFERENCE_PRESSURE_PA,
            temperature_k: reference.temperature_k,
            external_sum: reference
                .species_mole_fractions
                .iter()
                .map(|row| row.mole_fraction)
                .sum(),
            kithe_sum,
            major_relative: summary(species.iter().filter_map(|row| {
                (row.magnitude_class == StanjanChonMagnitudeClass::Major)
                    .then_some(row.relative_error)
                    .flatten()
            })),
            minor_log10: summary(species.iter().filter_map(|row| {
                (row.magnitude_class == StanjanChonMagnitudeClass::Minor)
                    .then_some(row.delta_log10)
                    .flatten()
            })),
            trace_log10: summary(species.iter().filter_map(|row| {
                (row.magnitude_class == StanjanChonMagnitudeClass::Trace)
                    .then_some(row.delta_log10)
                    .flatten()
            })),
            species,
            reaction_dimension: fixture.structure.reaction_dimension,
            phase_control_transitions: solution.phase_control_transitions(),
            residual_l2_norm: validation.residual_l2_norm,
            max_abs_element_balance_error: validation.max_abs_element_balance_error,
            accepted_backend: format!("{:?}", solution.solve_report().accepted_backend),
        })
    }
}

/// Reconstructs an elemental inventory from arbitrary molecular basis rows.
/// The formula parser makes the proof independent of the local repository and
/// intentionally permits `C5H12` even when no local thermochemistry exists.
pub(crate) fn element_totals_from_molecular_basis(
    amounts: impl IntoIterator<Item = (&'static str, f64)>,
) -> Result<BTreeMap<String, f64>, ReactionExtentError> {
    let mut totals = BTreeMap::new();
    for (formula, amount) in amounts {
        if !amount.is_finite() || amount < 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "argonne_stanjan_chon_feed",
                message: format!("amount for {formula} must be finite and non-negative"),
            });
        }
        let composition = parse_formula(formula.to_owned(), None).map_err(|error| {
            ReactionExtentError::InvalidProblem {
                field: "argonne_stanjan_chon_formula",
                message: format!("cannot parse {formula}: {error}"),
            }
        })?;
        for (element, count) in composition {
            *totals.entry(element).or_insert(0.0) += amount * count as f64;
        }
    }
    Ok(totals)
}

/// The source molecular input basis and the local executable basis must encode
/// the exact same closed C/H/O/N inventory.
pub(crate) fn verify_element_equivalent_feeds(
    reference: &ArgonneStanjanChonPtReference,
) -> Result<BTreeMap<String, f64>, ReactionExtentError> {
    let original = element_totals_from_molecular_basis(
        reference
            .original_reactants
            .iter()
            .map(|row| (row.species.as_str(), row.amount_mol))
            .map(|(species, amount)| match species {
                "C5H12" => ("C5H12", amount),
                "CH4" => ("CH4", amount),
                "O2" => ("O2", amount),
                "N2" => ("N2", amount),
                _ => unreachable!("validated frozen reactant identity"),
            }),
    )?;
    let local = element_totals_from_molecular_basis([
        ("CH4", 4.0),
        ("CO2", 2.0),
        ("O2", 8.0),
        ("N2", 37.6),
    ])?;
    let published = reference
        .element_totals
        .iter()
        .map(|row| (row.element.clone(), row.amount_mol_atoms))
        .collect::<BTreeMap<_, _>>();
    if original != local || original != published {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_element_equivalence",
            message: format!(
                "original={original:?}, local_equivalent={local:?}, published={published:?}"
            ),
        });
    }
    Ok(original)
}

/// Standard fixed `P,T` conditions for the published table.
pub(crate) fn stanjan_conditions(
    reference: &ArgonneStanjanChonPtReference,
) -> Result<EquilibriumConditions, ReactionExtentError> {
    EquilibriumConditions::new(
        reference.temperature_k,
        reference.pressure_pa,
        STANJAN_REFERENCE_PRESSURE_PA,
    )
}

/// Loads the reviewed STANJAN table through the common read-only loader.
pub(crate) fn load_argonne_stanjan_chon_dataset()
-> Result<FrozenReferenceDataset<ArgonneStanjanChonPtReference>, String> {
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/argonne_stanjan");
    FrozenReferenceDataset::load(
        directory.join("pentane_methane_air_tp_2500k_35atm.metadata.json"),
        directory.join("pentane_methane_air_tp_2500k_35atm.rows.json"),
    )
    .map_err(|error| error.to_string())
}

/// Builds the static 15-species NASA-gas substance-system spec used by the solve.
/// Search is restricted to `NASA_gas` and NIST fallback is disabled so that the
/// frozen universe resolves exactly from the local repository.
fn stanjan_gas_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([(
        STANJAN_GAS_PHASE.to_owned(),
        STANJAN_LOCAL_GAS_SPECIES
            .iter()
            .map(|name| (*name).to_owned())
            .collect(),
    )])))
    .with_phase_natures(Some(HashMap::from([(
        STANJAN_GAS_PHASE.to_owned(),
        PhysicalState::Gas,
    )])))
    .with_library_priorities(vec!["NASA_gas".to_owned()])
    .with_search_in_nist(false)
    .build()
    .expect("the static Argonne/STANJAN CHON gas universe must be structurally valid")
}

/// Verifies the resolved phase preserves the exact 15-species gas universe and
/// that every species carries an exact `NASA_gas` Thermo provenance record.
fn validate_exact_local_universe(
    resolved: &ResolvedPhaseSystem,
) -> Result<(), ReactionExtentError> {
    let phase = resolved
        .phase_specs()
        .iter()
        .find(|phase| phase.id().as_option().as_deref() == Some(STANJAN_GAS_PHASE))
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_phase",
            message: "missing declared gas phase".to_owned(),
        })?;
    if phase.physical_state() != PhysicalState::Gas
        || !phase
            .components()
            .iter()
            .map(String::as_str)
            .eq(STANJAN_LOCAL_GAS_SPECIES)
    {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_universe",
            message: "resolved phase did not preserve the exact 15-species gas universe".to_owned(),
        });
    }
    let phase_report = resolved
        .report()
        .phases()
        .iter()
        .find(|report| report.phase().as_option().as_deref() == Some(STANJAN_GAS_PHASE))
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_provenance",
            message: "missing gas lookup report".to_owned(),
        })?;
    for species in STANJAN_LOCAL_GAS_SPECIES {
        let exact_nasa = phase_report.search().rows().iter().any(|row| {
            row.property() == "Thermo"
                && row.substance() == species
                && row.library() == "NASA_gas"
                && !row.record_key().trim().is_empty()
        });
        if !exact_nasa {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "argonne_stanjan_chon_provenance",
                message: format!("{species} lacks exact NASA_gas Thermo provenance"),
            });
        }
    }
    Ok(())
}

/// Builds the per-species preflight table: temperature interval, provenance,
/// and molar Gibbs energy at the frozen target temperature.
fn preflight_rows(
    resolved: &ResolvedPhaseSystem,
    thermochemistry: &ResolvedThermochemistry,
    temperature_k: f64,
) -> Result<Vec<StanjanChonPreflightRow>, ReactionExtentError> {
    let phase_data = resolved
        .phase_data()
        .get(&Some(STANJAN_GAS_PHASE.to_owned()))
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_preflight",
            message: "missing local gas phase payload".to_owned(),
        })?;
    let gibbs = thermochemistry.evaluate_gibbs(temperature_k)?;
    STANJAN_LOCAL_GAS_SPECIES
        .iter()
        .map(|species| {
            let record = phase_data
                .get_search_result(species, WhatIsFound::Thermo)
                .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                    path: "argonne_stanjan_chon_preflight",
                    message: format!("missing Thermo record for {species}"),
                })?;
            let mut calculator = match record.calculator().cloned() {
                Some(CalculatorType::Thermo(calculator)) => calculator,
                Some(CalculatorType::Transport(_)) => {
                    return Err(ReactionExtentError::ValidationNotApplicable {
                        path: "argonne_stanjan_chon_preflight",
                        message: format!("{species} resolved transport instead of thermochemistry"),
                    });
                }
                None => {
                    return Err(ReactionExtentError::ValidationNotApplicable {
                        path: "argonne_stanjan_chon_preflight",
                        message: format!("{species} lacks a calculator"),
                    });
                }
            };
            calculator.parse_coefficients().map_err(|error| {
                ReactionExtentError::ValidationNotApplicable {
                    path: "argonne_stanjan_chon_preflight",
                    message: format!("{species} coefficients cannot be parsed: {error}"),
                }
            })?;
            let (lower, upper) = calculator.valid_temperature_interval().map_err(|error| {
                ReactionExtentError::ValidationNotApplicable {
                    path: "argonne_stanjan_chon_preflight",
                    message: format!("{species} has no valid temperature interval: {error}"),
                }
            })?;
            let index = component_index(resolved, species)?;
            let provenance = thermochemistry.provenance().get(index).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "thermochemistry provenance lacks {species}"
                ))
            })?;
            Ok(StanjanChonPreflightRow {
                external_identity: (*species).to_owned(),
                component: provenance.component().clone(),
                library: provenance.library().to_owned(),
                record_key: provenance.record_key().to_owned(),
                physical_state: PhysicalState::Gas,
                temperature_lower_k: lower,
                temperature_upper_k: upper,
                supports_target_temperature: (lower..=upper).contains(&temperature_k),
                gibbs_j_mol: gibbs[index],
            })
        })
        .collect()
}

/// Resolves the linear component index of a gas species within the resolved layout.
fn component_index(
    resolved: &ResolvedPhaseSystem,
    species: &str,
) -> Result<usize, ReactionExtentError> {
    let component = gas_component(species);
    resolved
        .layout()
        .components()
        .iter()
        .position(|candidate| candidate == &component)
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_component_layout",
            message: format!("missing exact component {}", component.label()),
        })
}

/// Computes component/element geometry and the reaction dimension from the
/// element-composition matrix, using SVD to estimate the element rank.
fn structure_from_matrix(
    element_composition: &DMatrix<f64>,
) -> Result<StanjanChonStructure, ReactionExtentError> {
    let component_count = element_composition.nrows();
    let element_count = element_composition.ncols();
    let element_rank = SVD::new(element_composition.clone(), false, false)
        .singular_values
        .iter()
        .filter(|&&value| value > 1.0e-12)
        .count();
    if component_count == 0 || element_count == 0 || element_rank == 0 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "argonne_stanjan_chon_structure",
            message: format!(
                "invalid component/element geometry: {component_count} components, {element_count} elements, rank {element_rank}"
            ),
        });
    }
    Ok(StanjanChonStructure {
        component_count,
        element_count,
        element_rank,
        reaction_dimension: component_count.saturating_sub(element_rank),
    })
}

/// Constructs the phase-qualified component id for a species in the gas phase.
fn gas_component(species: &str) -> PhaseComponentId {
    PhaseComponentId::new(PhaseId::new(Some(STANJAN_GAS_PHASE.to_owned())), species)
}

/// Classifies a published mole fraction into `Major`, `Minor`, or `Trace`,
/// selecting the error scale used by the comparison report.
fn magnitude_class(value: f64) -> StanjanChonMagnitudeClass {
    if value >= 1.0e-3 {
        StanjanChonMagnitudeClass::Major
    } else if value >= 1.0e-6 {
        StanjanChonMagnitudeClass::Minor
    } else {
        StanjanChonMagnitudeClass::Trace
    }
}

/// Aggregates an iterable of errors into a summary, or returns `None` when empty.
fn summary(values: impl Iterator<Item = f64>) -> Option<StanjanChonErrorSummary> {
    let values = values.collect::<Vec<_>>();
    (!values.is_empty()).then(|| StanjanChonErrorSummary {
        sample_count: values.len(),
        maximum_absolute: values.iter().map(|value| value.abs()).fold(0.0, f64::max),
        rms: (values.iter().map(|value| value * value).sum::<f64>() / values.len() as f64).sqrt(),
        mean_signed: values.iter().sum::<f64>() / values.len() as f64,
    })
}
