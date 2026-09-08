//! First complex I5 equilibrium characterization: NASA CEA `H2/O2` HP.
//!
//! The frozen CEA table is an external rounded answer. This module owns the
//! *local* side of the comparison: an exact declared phase universe, pinned
//! local identities, 1 kg reactant reconstruction, and a compact inventory
//! report. It does not add a CEA dependency, download anything, or expand the
//! species set selected by the publication.

use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::{DMatrix, linalg::SVD};

use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
    ReactionExtentError, ReactionExtentErrorKind,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::{
    ResolvedThermochemistry, ThermochemistryProvenance,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
    FixedPressureEnthalpySolution, PhMonolithicSeedAttemptReport, PhSolveMode, PhSolvePath,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_solver_policy::{
    SolverAttemptReport, SolverTermination,
};
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, NasaCeaH2O2HpReference,
};
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

/// Semantic phase labels used only by the frozen CEA case declaration.
pub(crate) const CEA_GAS_PHASE: &str = "gas";
pub(crate) const CEA_LIQUID_PHASE: &str = "liquid";
pub(crate) const CEA_ICE_PHASE: &str = "solid";

/// Exact CEA identities, retained in the same semantic order as the published
/// output table. The phase mapping below, rather than this presentation order,
/// determines the numerical solver layout.
pub(crate) const CEA_DECLARED_SPECIES: [&str; 11] = [
    "H", "H2", "H2O", "H2O2", "HO2", "O", "O2", "O3", "OH", "H2O(L)", "H2O(cr)",
];

const CEA_GAS_SPECIES: [&str; 9] = ["H", "H2", "H2O", "H2O2", "HO2", "O", "O2", "O3", "OH"];

/// Exact semantic mapping of one NASA CEA output identity to one local,
/// phase-qualified record. The condensed water strings deliberately reuse the
/// identities already pinned by IAPWS fixtures; they are not guessed aliases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NasaCeaComponentMapping {
    /// Identity as printed by the NASA CEA output table.
    pub(crate) cea_identity: &'static str,
    /// Local phase label for this component.
    pub(crate) phase: &'static str,
    /// Local record species name resolved against the repository.
    pub(crate) local_species: &'static str,
    /// Library expected to provide the exact local Thermo record.
    pub(crate) library: &'static str,
    /// Physical state required by the CEA case.
    pub(crate) physical_state: PhysicalState,
}

pub(crate) const CEA_COMPONENT_MAPPINGS: [NasaCeaComponentMapping; 11] = [
    NasaCeaComponentMapping {
        cea_identity: "H",
        phase: CEA_GAS_PHASE,
        local_species: "H",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "H2",
        phase: CEA_GAS_PHASE,
        local_species: "H2",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "H2O",
        phase: CEA_GAS_PHASE,
        local_species: "H2O",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "H2O2",
        phase: CEA_GAS_PHASE,
        local_species: "H2O2",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "HO2",
        phase: CEA_GAS_PHASE,
        local_species: "HO2",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "O",
        phase: CEA_GAS_PHASE,
        local_species: "O",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "O2",
        phase: CEA_GAS_PHASE,
        local_species: "O2",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "O3",
        phase: CEA_GAS_PHASE,
        local_species: "O3",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "OH",
        phase: CEA_GAS_PHASE,
        local_species: "OH",
        library: "NASA_gas",
        physical_state: PhysicalState::Gas,
    },
    NasaCeaComponentMapping {
        cea_identity: "H2O(L)",
        phase: CEA_LIQUID_PHASE,
        local_species: "H2O(L)",
        library: "NASA_cond",
        physical_state: PhysicalState::Liquid,
    },
    NasaCeaComponentMapping {
        cea_identity: "H2O(cr)",
        phase: CEA_ICE_PHASE,
        local_species: "H2O(s)",
        library: "NASA_cond",
        physical_state: PhysicalState::Solid,
    },
];

/// A resolved component inventory row suitable for a diagnostic table.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaComponentInventoryRow {
    /// CEA output identity of the component.
    pub(crate) cea_identity: String,
    /// Exact local phase-qualified component.
    pub(crate) component: PhaseComponentId,
    /// Library that supplied the record.
    pub(crate) library: String,
    /// Repository record key backing the coefficients.
    pub(crate) record_key: String,
    /// Physical state of the component.
    pub(crate) physical_state: PhysicalState,
    /// Lower bound of the common local temperature interval, K.
    pub(crate) temperature_lower_k: f64,
    /// Upper bound of the common local temperature interval, K.
    pub(crate) temperature_upper_k: f64,
    /// Molar Gibbs energy at the evaluation temperature, J/mol.
    pub(crate) gibbs_j_mol: f64,
    /// Molar enthalpy at the evaluation temperature, J/mol.
    pub(crate) enthalpy_j_mol: f64,
    /// Molar isobaric heat capacity at the evaluation temperature, J/(mol·K).
    pub(crate) heat_capacity_j_mol_k: f64,
}

/// One exact local record interval examined before the CEA solver can start.
///
/// It intentionally contains no extrapolated G/H/Cp value. A component may be
/// chemically identifiable while still being unusable at the published CEA
/// temperature; reporting that distinction is the entire point of preflight.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaComponentDomainRow {
    /// CEA output identity of the component.
    pub(crate) cea_identity: String,
    /// Exact local phase-qualified component.
    pub(crate) component: PhaseComponentId,
    /// Library that supplied the record.
    pub(crate) library: String,
    /// Repository record key backing the coefficients.
    pub(crate) record_key: String,
    /// Lower bound of the record's valid temperature interval, K.
    pub(crate) temperature_lower_k: f64,
    /// Upper bound of the record's valid temperature interval, K.
    pub(crate) temperature_upper_k: f64,
    /// Whether the reactant temperature lies inside the record interval.
    pub(crate) supports_reactant_temperature: bool,
    /// Whether the equilibrium temperature lies inside the record interval.
    pub(crate) supports_equilibrium_temperature: bool,
}

/// Input reconstructed from the CEA mass-ratio convention and local records.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaReactantReconstruction {
    /// Total reconstructed reactant mass (always one kilogram).
    pub(crate) total_mass_kg: f64,
    /// Hydrogen mass implied by the O/F mass ratio, kg.
    pub(crate) hydrogen_mass_kg: f64,
    /// Oxygen mass implied by the O/F mass ratio, kg.
    pub(crate) oxygen_mass_kg: f64,
    /// Local molar mass of H2, kg/mol.
    pub(crate) hydrogen_molar_mass_kg_mol: f64,
    /// Local molar mass of O2, kg/mol.
    pub(crate) oxygen_molar_mass_kg_mol: f64,
    /// Reconstructed H2 moles.
    pub(crate) hydrogen_moles: f64,
    /// Reconstructed O2 moles.
    pub(crate) oxygen_moles: f64,
    /// Reconstructed hydrogen atom moles (`2 * hydrogen_moles`).
    pub(crate) hydrogen_atom_moles: f64,
    /// Reconstructed oxygen atom moles (`2 * oxygen_moles`).
    pub(crate) oxygen_atom_moles: f64,
    /// Enthalpy target of the fixed-`P,H` problem, J.
    pub(crate) target_enthalpy_j: f64,
}

/// Derived structural evidence that this is a general multicomponent problem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NasaCeaMulticomponentStructure {
    /// Number of species in the declared CEA universe.
    pub(crate) component_count: usize,
    /// Number of distinct chemical elements present.
    pub(crate) element_count: usize,
    /// Rank of the element-composition matrix.
    pub(crate) element_rank: usize,
    /// Dimension of the reaction space (`component_count - element_rank`).
    pub(crate) reaction_dimension: usize,
}

/// Immutable local side of the first NASA CEA I5 benchmark.
#[derive(Clone)]
pub(crate) struct ResolvedNasaCeaH2O2HpFixture {
    resolved: ResolvedPhaseSystem,
    thermochemistry: ResolvedThermochemistry,
    component_inventory: Vec<NasaCeaComponentInventoryRow>,
    molar_masses_g_mol: HashMap<String, f64>,
    structure: NasaCeaMulticomponentStructure,
}

/// Resolved nine-component gas subset used by the first executable CEA
/// characterization. The two published zero condensed rows remain external
/// evidence, but are not falsely presented as locally TPD-validated.
#[derive(Clone)]
pub(crate) struct ResolvedNasaCeaH2O2HpGasFixture {
    resolved: ResolvedPhaseSystem,
    thermochemistry: ResolvedThermochemistry,
    component_inventory: Vec<NasaCeaComponentInventoryRow>,
    structure: NasaCeaMulticomponentStructure,
}

impl ResolvedNasaCeaH2O2HpGasFixture {
    pub(crate) fn resolve_offline(
        repository: Arc<ThermoRepository>,
        reference: &NasaCeaH2O2HpReference,
    ) -> Result<Self, ReactionExtentError> {
        ensure_gas_only_reference_eligibility(reference)?;
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            gas_benchmark_spec(),
            repository,
        )
        .map_err(|error| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_gas_inventory",
            message: format!("exact nine-component CEA gas universe could not resolve: {error}"),
        })?;
        if resolved.report().nist_fallback_enabled() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_cea_h2_o2_hp_nist",
                message: "frozen CEA benchmark must not enable NIST fallback".to_owned(),
            });
        }
        validate_mapping_subset(&resolved, gas_mappings())?;
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        let bounds = thermochemistry.temperature_bounds();
        for temperature in [
            reference.reactant_temperature_k,
            reference.equilibrium_temperature_k,
        ] {
            if !bounds.contains(temperature) {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_gas_temperature_domain",
                    message: format!(
                        "common gas interval [{}, {}] K excludes required {temperature} K",
                        bounds.lower(),
                        bounds.upper()
                    ),
                });
            }
        }
        let (matrix, _, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_gas_elements",
                    message: format!("gas component composition is unavailable: {error}"),
                })?;
        let structure = structural_evidence(&matrix)?;
        let component_inventory = component_inventory(
            &resolved,
            &thermochemistry,
            reference.equilibrium_temperature_k,
            gas_mappings(),
        )?;
        Ok(Self {
            resolved,
            thermochemistry,
            component_inventory,
            structure,
        })
    }

    pub(crate) fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub(crate) fn thermochemistry(&self) -> &ResolvedThermochemistry {
        &self.thermochemistry
    }

    pub(crate) fn component_inventory(&self) -> &[NasaCeaComponentInventoryRow] {
        &self.component_inventory
    }

    pub(crate) fn structure(&self) -> NasaCeaMulticomponentStructure {
        self.structure
    }

    pub(crate) fn initial_composition(
        &self,
        reactants: &NasaCeaReactantReconstruction,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        MultiphaseInitialComposition::from_sparse(
            &layout,
            vec![
                (
                    PhaseComponentId::new(PhaseId::new(Some(CEA_GAS_PHASE.to_owned())), "H2"),
                    reactants.hydrogen_moles,
                ),
                (
                    PhaseComponentId::new(PhaseId::new(Some(CEA_GAS_PHASE.to_owned())), "O2"),
                    reactants.oxygen_moles,
                ),
            ],
        )
    }

    /// Creates a positive, element-conserving diagnostic gas seed without
    /// borrowing the external CEA equilibrium composition. Each radical and
    /// molecular product receives the same small local amount; the remaining
    /// H/O inventory stays in H2 and O2 exactly.
    ///
    /// This is deliberately a diagnostic basin probe, not the production
    /// initial state. The normal benchmark always starts from physical H2/O2
    /// reactants and relies on the trace-seed policy for absent species.
    pub(crate) fn broad_initial_composition(
        &self,
        reactants: &NasaCeaReactantReconstruction,
        fraction: f64,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        if !fraction.is_finite() || !(0.0 < fraction && fraction < 1.0) {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_cea_h2_o2_hp_broad_seed_fraction",
                message: "broad gas-seed fraction must be finite and lie in (0, 1)".to_owned(),
            });
        }
        let hydrogen_atoms = reactants.hydrogen_atom_moles;
        let oxygen_atoms = reactants.oxygen_atom_moles;
        // The seven non-reactant species below consume 7*q H atoms and 10*q
        // O atoms. The remaining atoms reconstruct positive H2/O2 amounts.
        let q = fraction * (hydrogen_atoms / 7.0).min(oxygen_atoms / 10.0);
        let hydrogen_moles = (hydrogen_atoms - 7.0 * q) / 2.0;
        let oxygen_moles = (oxygen_atoms - 10.0 * q) / 2.0;
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        MultiphaseInitialComposition::from_sparse(
            &layout,
            vec![
                (gas_component("H"), q),
                (gas_component("H2"), hydrogen_moles),
                (gas_component("H2O"), q),
                (gas_component("H2O2"), q),
                (gas_component("HO2"), q),
                (gas_component("O"), q),
                (gas_component("O2"), oxygen_moles),
                (gas_component("O3"), q),
                (gas_component("OH"), q),
            ],
        )
    }
}

/// Diagnostic magnitude class for rounded external species amounts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NasaCeaSpeciesMagnitudeClass {
    /// External amount at or above `1e-3` kgmol/kg; compared on the relative scale.
    Major,
    /// External amount between `1e-6` and `1e-3` kgmol/kg; compared as `delta_log10`.
    Minor,
    /// External amount below `1e-6` kgmol/kg; compared as `delta_log10`.
    Trace,
    /// Published zero condensed row excluded from the local gas-only solve.
    ExternallyAbsentExcluded,
}

/// Identity-based comparison of one CEA row with the accepted KiThe state.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaSpeciesComparisonRow {
    /// CEA output identity for this row.
    pub(crate) cea_identity: String,
    /// Matched local component, or `None` when excluded from the gas-only solve.
    pub(crate) local_component: Option<PhaseComponentId>,
    /// Semantic magnitude class selecting the error scale.
    pub(crate) magnitude_class: NasaCeaSpeciesMagnitudeClass,
    /// Amount published by NASA CEA, kg-mol/kg.
    pub(crate) cea_amount_kgmol_per_kg: f64,
    /// Local (KiThe) amount, kg-mol/kg, or `None` when excluded.
    pub(crate) kithe_amount_kgmol_per_kg: Option<f64>,
    /// Absolute error `local - cea`, or `None` when excluded.
    pub(crate) absolute_error: Option<f64>,
    /// Relative error, or `None` when the CEA value is not positive.
    pub(crate) relative_error: Option<f64>,
    /// Log10 difference (`log10(local) - log10(cea)`), or `None` when not computable.
    pub(crate) delta_log10: Option<f64>,
    /// Human-readable reason for exclusion, or `None` when included.
    pub(crate) exclusion_reason: Option<String>,
}

/// Aggregate relative-error statistics for a positive external magnitude
/// class. These values are characterization data, never solver acceptance
/// gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NasaCeaRelativeErrorSummary {
    /// Number of sampled `Major` species.
    pub(crate) sample_count: usize,
    /// Largest absolute relative error among the samples.
    pub(crate) max_absolute_relative_error: f64,
    /// Root-mean-square relative error.
    pub(crate) rms_relative_error: f64,
    /// Signed mean relative error (bias indicator).
    pub(crate) mean_signed_relative_error: f64,
}

/// Aggregate logarithmic comparison statistics for minor or trace species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NasaCeaLogErrorSummary {
    /// Number of sampled minor/trace species.
    pub(crate) sample_count: usize,
    /// Largest absolute `delta_log10` among the samples.
    pub(crate) max_absolute_delta_log10: f64,
    /// Root-mean-square `delta_log10`.
    pub(crate) rms_delta_log10: f64,
    /// Signed mean `delta_log10` (bias indicator).
    pub(crate) mean_signed_delta_log10: f64,
}

/// Compact source-facing characterization metrics derived from the
/// identity-aligned species rows. Excluded condensed rows are deliberately
/// absent from every numeric aggregate.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaExternalComparisonSummary {
    /// Absolute equilibrium-temperature deviation, K.
    pub(crate) absolute_temperature_delta_k: f64,
    /// Relative equilibrium-temperature deviation.
    pub(crate) relative_temperature_delta: f64,
    /// Absolute total-amount deviation, kg-mol/kg.
    pub(crate) absolute_total_amount_delta_kgmol_per_kg: f64,
    /// Relative total-amount deviation.
    pub(crate) relative_total_amount_delta: f64,
    /// Relative-error aggregate for `Major` species, when any are sampled.
    pub(crate) major_relative: Option<NasaCeaRelativeErrorSummary>,
    /// Log-error aggregate for `Minor` species, when any are sampled.
    pub(crate) minor_log10: Option<NasaCeaLogErrorSummary>,
    /// Log-error aggregate for `Trace` species, when any are sampled.
    pub(crate) trace_log10: Option<NasaCeaLogErrorSummary>,
}

/// Typed external-characterization report for the first executable CEA case.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaEquilibriumComparisonReport {
    /// Identifier of the frozen CEA dataset that produced this report.
    pub(crate) dataset_id: String,
    /// Physical (total) pressure of the CEA case, Pa.
    pub(crate) physical_pressure_pa: f64,
    /// Standard pressure used by the local activity model, Pa.
    pub(crate) activity_reference_pressure_pa: f64,
    /// Reactant temperature of the fixed-`P,H` case, K.
    pub(crate) reactant_temperature_k: f64,
    /// Enthalpy target solved by the local engine, J.
    pub(crate) target_enthalpy_j: f64,
    /// Equilibrium temperature published by CEA, K.
    pub(crate) cea_temperature_k: f64,
    /// Equilibrium temperature returned by the local solve, K.
    pub(crate) kithe_temperature_k: f64,
    /// Temperature difference `kithe - cea`, K.
    pub(crate) temperature_delta_k: f64,
    /// Total equilibrium amount published by CEA, kg-mol/kg.
    pub(crate) cea_total_kgmol_per_kg: f64,
    /// Total equilibrium amount returned by the local solve, kg-mol/kg.
    pub(crate) kithe_total_kgmol_per_kg: f64,
    /// Compact aggregate external comparison metrics.
    pub(crate) external_summary: NasaCeaExternalComparisonSummary,
    /// Identity-aligned per-species comparison rows.
    pub(crate) species: Vec<NasaCeaSpeciesComparisonRow>,
    /// Solve path taken by the fixed-`P,H` engine.
    pub(crate) solve_path: PhSolvePath,
    /// Human-readable fallback reason, if a fallback occurred.
    pub(crate) fallback_reason: Option<String>,
    /// Phase-control transitions in the accepted equilibrium state.
    pub(crate) final_phase_transitions: usize,
    /// Phase-control events recorded on the solution path.
    pub(crate) trial_phase_events: usize,
    /// Absolute enthalpy closure error of the accepted solution, J.
    pub(crate) enthalpy_error_j: f64,
    /// Scaled (dimensionless) enthalpy closure error.
    pub(crate) scaled_enthalpy_error: f64,
    /// Residual L2 norm from the accepted solution's validation.
    pub(crate) residual_l2_norm: f64,
    /// Maximum absolute element-balance error from the accepted solution.
    pub(crate) max_abs_element_balance_error: f64,
}

/// One compact backend trace entry retained for a direct P,H route probe.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct NasaCeaBackendAttemptDiagnostic {
    /// Debug name of the nonlinear backend.
    pub(crate) backend: String,
    /// Debug representation of the attempt outcome.
    pub(crate) outcome: String,
    /// Failure kind, when the attempt failed.
    pub(crate) failure_kind: Option<String>,
    /// Failure reason text, when available.
    pub(crate) reason: Option<String>,
    /// Solver termination state, when metrics were reported.
    pub(crate) termination: Option<SolverTermination>,
    /// Number of nonlinear iterations, when metrics were reported.
    pub(crate) iterations: Option<usize>,
    /// Number of residual evaluations, when metrics were reported.
    pub(crate) residual_evaluations: Option<usize>,
    /// Number of Jacobian evaluations, when metrics were reported.
    pub(crate) jacobian_evaluations: Option<usize>,
}

/// A direct route outcome for the exact CEA gas benchmark.
///
/// It preserves only one compact solver-cascade layer. Full backend traces
/// remain available on the original typed error and do not belong in the
/// external-comparison report.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaPhRouteDiagnostic {
    /// The solver mode explicitly requested by the caller.
    pub(crate) requested_mode: PhSolveMode,
    /// Temperature seed used to start the direct route, K.
    pub(crate) temperature_seed_k: f64,
    /// Trace mole floor applied to absent species.
    pub(crate) trace_floor: f64,
    /// Path actually taken, or `None` when the route failed.
    pub(crate) solve_path: Option<PhSolvePath>,
    /// Human-readable fallback reason, when a fallback occurred.
    pub(crate) fallback_reason: Option<String>,
    /// Error kind when the route failed, or `None` on success.
    pub(crate) failure_kind: Option<ReactionExtentErrorKind>,
    /// Error message when the route failed, or `None` on success.
    pub(crate) failure_message: Option<String>,
    /// Compact backend-attempt trace.
    pub(crate) backend_attempts: Vec<NasaCeaBackendAttemptDiagnostic>,
    /// Temperature-seed recovery attempts, when applicable.
    pub(crate) temperature_seed_attempts: Vec<PhMonolithicSeedAttemptReport>,
    /// Inner backend attempts when a nested route was used.
    pub(crate) inner_backend_attempts: Option<usize>,
    /// Inner nonlinear iterations when a nested route was used.
    pub(crate) inner_nonlinear_iterations: Option<usize>,
    /// Number of temperature trials recorded, when available.
    pub(crate) temperature_trials: Option<usize>,
    /// External comparison report, present only when the route succeeded.
    pub(crate) comparison: Option<NasaCeaEquilibriumComparisonReport>,
}

impl NasaCeaEquilibriumComparisonReport {
    pub(crate) fn from_solution(
        dataset: &FrozenReferenceDataset<NasaCeaH2O2HpReference>,
        reactants: &NasaCeaReactantReconstruction,
        activity_reference_pressure_pa: f64,
        solution: &FixedPressureEnthalpySolution,
    ) -> Result<Self, ReactionExtentError> {
        let reference =
            dataset
                .rows()
                .first()
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "nasa_cea_h2_o2_hp_reference",
                    message: "frozen CEA dataset has no case".to_owned(),
                })?;
        ensure_gas_only_reference_eligibility(reference)?;
        let accepted = solution.equilibrium();
        let species = reference
            .species_amounts
            .iter()
            .map(|external| {
                let mapping = CEA_COMPONENT_MAPPINGS
                    .iter()
                    .find(|mapping| mapping.cea_identity == external.species)
                    .expect("validated frozen species must have a static mapping");
                if mapping.physical_state != PhysicalState::Gas {
                    return NasaCeaSpeciesComparisonRow {
                        cea_identity: external.species.clone(),
                        local_component: None,
                        magnitude_class:
                            NasaCeaSpeciesMagnitudeClass::ExternallyAbsentExcluded,
                        cea_amount_kgmol_per_kg: external.amount_kgmol_per_kg,
                        kithe_amount_kgmol_per_kg: None,
                        absolute_error: None,
                        relative_error: None,
                        delta_log10: None,
                        exclusion_reason: Some(format!(
                            "published zero component excluded from local solve: {} record domain does not cover the CEA HP temperature",
                            mapping.local_species
                        )),
                    };
                }
                let component = PhaseComponentId::new(
                    PhaseId::new(Some(mapping.phase.to_owned())),
                    mapping.local_species,
                );
                let local = accepted
                    .moles_for(&component)
                    .expect("gas fixture and static mapping must have identical identities")
                    / 1_000.0;
                let absolute_error = local - external.amount_kgmol_per_kg;
                let relative_error = (external.amount_kgmol_per_kg > 0.0)
                    .then_some(absolute_error / external.amount_kgmol_per_kg);
                let delta_log10 = (local > 0.0 && external.amount_kgmol_per_kg > 0.0)
                    .then_some(local.log10() - external.amount_kgmol_per_kg.log10());
                let magnitude_class = if external.amount_kgmol_per_kg >= 1.0e-3 {
                    NasaCeaSpeciesMagnitudeClass::Major
                } else if external.amount_kgmol_per_kg >= 1.0e-6 {
                    NasaCeaSpeciesMagnitudeClass::Minor
                } else {
                    NasaCeaSpeciesMagnitudeClass::Trace
                };
                NasaCeaSpeciesComparisonRow {
                    cea_identity: external.species.clone(),
                    local_component: Some(component),
                    magnitude_class,
                    cea_amount_kgmol_per_kg: external.amount_kgmol_per_kg,
                    kithe_amount_kgmol_per_kg: Some(local),
                    absolute_error: Some(absolute_error),
                    relative_error,
                    delta_log10,
                    exclusion_reason: None,
                }
            })
            .collect::<Vec<_>>();
        let kithe_total_kgmol_per_kg = accepted.component_moles().iter().sum::<f64>() / 1_000.0;
        let external_summary = external_comparison_summary(
            reference,
            solution.temperature(),
            kithe_total_kgmol_per_kg,
            &species,
        );
        let validation = accepted.accepted_solution().validation();
        Ok(Self {
            dataset_id: dataset.metadata().dataset_id.clone(),
            physical_pressure_pa: reference.pressure_pa,
            activity_reference_pressure_pa,
            reactant_temperature_k: reference.reactant_temperature_k,
            target_enthalpy_j: reactants.target_enthalpy_j,
            cea_temperature_k: reference.equilibrium_temperature_k,
            kithe_temperature_k: solution.temperature(),
            temperature_delta_k: solution.temperature() - reference.equilibrium_temperature_k,
            cea_total_kgmol_per_kg: reference.total_amount_kgmol_per_kg,
            kithe_total_kgmol_per_kg,
            external_summary,
            species,
            solve_path: solution.report().solve_path(),
            fallback_reason: solution
                .report()
                .fallback_reason()
                .map(|reason| format!("{:?}: {}", reason.error_kind(), reason.message())),
            final_phase_transitions: accepted.phase_control_transitions(),
            trial_phase_events: solution.report().phase_control_transitions(),
            enthalpy_error_j: solution.enthalpy_error(),
            scaled_enthalpy_error: solution.scaled_enthalpy_error(),
            residual_l2_norm: validation.residual_l2_norm,
            max_abs_element_balance_error: validation.max_abs_element_balance_error,
        })
    }

    /// Count of condensed rows excluded from the gas-only external comparison.
    pub(crate) fn excluded_condensed_count(&self) -> usize {
        self.species
            .iter()
            .filter(|row| {
                row.magnitude_class == NasaCeaSpeciesMagnitudeClass::ExternallyAbsentExcluded
            })
            .count()
    }
}

impl NasaCeaPhRouteDiagnostic {
    /// Converts one explicit route result into compact, comparison-safe
    /// diagnostics. The caller owns the requested mode and seed, so an `Auto`
    /// recovery never stands in for a direct monolithic experiment.
    pub(crate) fn from_result(
        dataset: &FrozenReferenceDataset<NasaCeaH2O2HpReference>,
        reactants: &NasaCeaReactantReconstruction,
        activity_reference_pressure_pa: f64,
        requested_mode: PhSolveMode,
        temperature_seed_k: f64,
        trace_floor: f64,
        result: Result<FixedPressureEnthalpySolution, ReactionExtentError>,
    ) -> Result<Self, ReactionExtentError> {
        match result {
            Ok(solution) => {
                let report = solution.report();
                let backend_attempts = report
                    .monolithic_evidence()
                    .map(|evidence| compact_attempts(&evidence.solve_report().attempts))
                    .unwrap_or_default();
                let temperature_seed_attempts = report
                    .monolithic_evidence()
                    .and_then(|evidence| evidence.temperature_seed_report())
                    .map(|seed_report| seed_report.attempts().to_vec())
                    .unwrap_or_default();
                Ok(Self {
                    requested_mode,
                    temperature_seed_k,
                    trace_floor,
                    solve_path: Some(report.solve_path()),
                    fallback_reason: report
                        .fallback_reason()
                        .map(|reason| format!("{:?}: {}", reason.error_kind(), reason.message())),
                    failure_kind: None,
                    failure_message: None,
                    backend_attempts,
                    temperature_seed_attempts,
                    inner_backend_attempts: Some(report.inner_backend_attempts()),
                    inner_nonlinear_iterations: Some(report.inner_nonlinear_iterations()),
                    temperature_trials: Some(report.trials().len()),
                    comparison: Some(NasaCeaEquilibriumComparisonReport::from_solution(
                        dataset,
                        reactants,
                        activity_reference_pressure_pa,
                        &solution,
                    )?),
                })
            }
            Err(error) => Ok(Self {
                requested_mode,
                temperature_seed_k,
                trace_floor,
                solve_path: None,
                fallback_reason: None,
                failure_kind: Some(error.kind()),
                failure_message: Some(error.to_string()),
                backend_attempts: compact_error_attempts(&error),
                temperature_seed_attempts: temperature_seed_attempts_from_error(&error),
                inner_backend_attempts: None,
                inner_nonlinear_iterations: None,
                temperature_trials: None,
                comparison: None,
            }),
        }
    }
}

/// Verifies every published condensed CEA component has exactly zero amount,
/// so a reduced gas-only local solve is a faithful representation.
fn ensure_gas_only_reference_eligibility(
    reference: &NasaCeaH2O2HpReference,
) -> Result<(), ReactionExtentError> {
    for mapping in CEA_COMPONENT_MAPPINGS
        .iter()
        .filter(|mapping| mapping.physical_state != PhysicalState::Gas)
    {
        let amount = reference
            .species_amounts
            .iter()
            .find(|row| row.species == mapping.cea_identity)
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_gas_only_eligibility",
                message: format!(
                    "external CEA row '{}' is required before a gas-only comparison",
                    mapping.cea_identity
                ),
            })?;
        if amount.amount_kgmol_per_kg != 0.0 {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_gas_only_eligibility",
                message: format!(
                    "cannot exclude external condensed component '{}' with positive amount {} kgmol/kg",
                    mapping.cea_identity, amount.amount_kgmol_per_kg
                ),
            });
        }
    }
    Ok(())
}

/// Builds the compact temperature, total-amount, and per-class error summary
/// from identity-aligned species rows.
fn external_comparison_summary(
    reference: &NasaCeaH2O2HpReference,
    kithe_temperature_k: f64,
    kithe_total_kgmol_per_kg: f64,
    species: &[NasaCeaSpeciesComparisonRow],
) -> NasaCeaExternalComparisonSummary {
    let temperature_delta = kithe_temperature_k - reference.equilibrium_temperature_k;
    let total_delta = kithe_total_kgmol_per_kg - reference.total_amount_kgmol_per_kg;
    NasaCeaExternalComparisonSummary {
        absolute_temperature_delta_k: temperature_delta.abs(),
        relative_temperature_delta: temperature_delta / reference.equilibrium_temperature_k,
        absolute_total_amount_delta_kgmol_per_kg: total_delta.abs(),
        relative_total_amount_delta: total_delta / reference.total_amount_kgmol_per_kg,
        major_relative: relative_error_summary(species.iter().filter_map(|row| {
            (row.magnitude_class == NasaCeaSpeciesMagnitudeClass::Major)
                .then_some(row.relative_error)
                .flatten()
        })),
        minor_log10: log_error_summary(species.iter().filter_map(|row| {
            (row.magnitude_class == NasaCeaSpeciesMagnitudeClass::Minor)
                .then_some(row.delta_log10)
                .flatten()
        })),
        trace_log10: log_error_summary(species.iter().filter_map(|row| {
            (row.magnitude_class == NasaCeaSpeciesMagnitudeClass::Trace)
                .then_some(row.delta_log10)
                .flatten()
        })),
    }
}

/// Aggregates `Major`-class relative errors, or returns `None` when empty.
fn relative_error_summary(
    values: impl Iterator<Item = f64>,
) -> Option<NasaCeaRelativeErrorSummary> {
    let values = values.collect::<Vec<_>>();
    (!values.is_empty()).then(|| NasaCeaRelativeErrorSummary {
        sample_count: values.len(),
        max_absolute_relative_error: values.iter().map(|value| value.abs()).fold(0.0, f64::max),
        rms_relative_error: (values.iter().map(|value| value * value).sum::<f64>()
            / values.len() as f64)
            .sqrt(),
        mean_signed_relative_error: values.iter().sum::<f64>() / values.len() as f64,
    })
}

/// Aggregates `Minor`/`Trace`-class `delta_log10` errors, or `None` when empty.
fn log_error_summary(values: impl Iterator<Item = f64>) -> Option<NasaCeaLogErrorSummary> {
    let values = values.collect::<Vec<_>>();
    (!values.is_empty()).then(|| NasaCeaLogErrorSummary {
        sample_count: values.len(),
        max_absolute_delta_log10: values.iter().map(|value| value.abs()).fold(0.0, f64::max),
        rms_delta_log10: (values.iter().map(|value| value * value).sum::<f64>()
            / values.len() as f64)
            .sqrt(),
        mean_signed_delta_log10: values.iter().sum::<f64>() / values.len() as f64,
    })
}

/// Flattens a typed solver error into compact backend-attempt diagnostics.
fn compact_error_attempts(error: &ReactionExtentError) -> Vec<NasaCeaBackendAttemptDiagnostic> {
    let mut attempts = Vec::new();
    append_error_attempts(error, &mut attempts);
    attempts
}

/// Recursively appends backend attempts from a nested typed solver error tree.
fn append_error_attempts(
    error: &ReactionExtentError,
    output: &mut Vec<NasaCeaBackendAttemptDiagnostic>,
) {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => {
            output.extend(compact_attempts(attempts));
        }
        ReactionExtentError::TemperatureTrialFailed { cause, .. }
        | ReactionExtentError::TemperatureRangePointFailed { cause, .. } => {
            append_error_attempts(cause, output);
        }
        ReactionExtentError::PhAutoFallbackFailed { monolithic, nested } => {
            append_error_attempts(monolithic, output);
            append_error_attempts(nested, output);
        }
        ReactionExtentError::PhMonolithicSeedRecoveryFailed { attempts } => {
            for attempt in attempts {
                append_error_attempts(attempt.cause(), output);
            }
        }
        _ => {}
    }
}

/// Reconstructs rejected temperature-seed attempt reports from an error.
fn temperature_seed_attempts_from_error(
    error: &ReactionExtentError,
) -> Vec<PhMonolithicSeedAttemptReport> {
    match error {
        ReactionExtentError::PhMonolithicSeedRecoveryFailed { attempts } => attempts
            .iter()
            .map(|attempt| PhMonolithicSeedAttemptReport {
                temperature_seed_k: attempt.temperature_seed_k(),
                accepted: false,
                failure_kind: Some(attempt.cause().kind()),
                error: Some(attempt.cause().to_string()),
                started_backend_attempts: started_backend_attempts_for_error(attempt.cause()),
                nonlinear_iterations: nonlinear_iterations_for_error(attempt.cause()),
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Counts backends that were actually started within a cascade error.
fn started_backend_attempts_for_error(error: &ReactionExtentError) -> usize {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts
            .iter()
            .filter(|attempt| attempt.is_started())
            .count(),
        _ => 0,
    }
}

/// Sums reported nonlinear iterations across a cascade error's attempts.
fn nonlinear_iterations_for_error(error: &ReactionExtentError) -> usize {
    match error {
        ReactionExtentError::AllBackendsFailed { attempts }
        | ReactionExtentError::CascadeAborted { attempts, .. } => attempts
            .iter()
            .filter_map(|attempt| attempt.metrics.as_ref())
            .map(|metrics| metrics.iterations)
            .sum(),
        _ => 0,
    }
}

/// Converts a slice of solver-attempt reports into compact diagnostics.
fn compact_attempts(attempts: &[SolverAttemptReport]) -> Vec<NasaCeaBackendAttemptDiagnostic> {
    attempts
        .iter()
        .map(|attempt| NasaCeaBackendAttemptDiagnostic {
            backend: format!("{:?}", attempt.backend),
            outcome: format!("{:?}", attempt.outcome),
            failure_kind: attempt
                .outcome
                .failure_kind()
                .map(|kind| format!("{kind:?}")),
            reason: attempt.outcome.reason().map(ToOwned::to_owned),
            termination: attempt.metrics.as_ref().map(|metrics| metrics.termination),
            iterations: attempt.metrics.as_ref().map(|metrics| metrics.iterations),
            residual_evaluations: attempt
                .metrics
                .as_ref()
                .map(|metrics| metrics.residual_evaluations),
            jacobian_evaluations: attempt
                .metrics
                .as_ref()
                .map(|metrics| metrics.jacobian_evaluations),
        })
        .collect()
}

/// Constructs a phase-qualified component id in the CEA gas phase.
fn gas_component(substance: &str) -> PhaseComponentId {
    PhaseComponentId::new(PhaseId::new(Some(CEA_GAS_PHASE.to_owned())), substance)
}

/// Outcome of the exact local CEA component inventory before a solver is
/// permitted to construct equations. This report is available even when the
/// selected records do not cover the full benchmark temperature domain.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaCeaH2O2HpPreflightReport {
    rows: Vec<NasaCeaComponentDomainRow>,
    common_temperature_lower_k: f64,
    common_temperature_upper_k: f64,
}

impl NasaCeaH2O2HpPreflightReport {
    pub(crate) fn rows(&self) -> &[NasaCeaComponentDomainRow] {
        &self.rows
    }

    /// Whether every component record covers both required temperatures.
    pub(crate) fn supports_required_temperatures(&self) -> bool {
        self.rows
            .iter()
            .all(|row| row.supports_reactant_temperature && row.supports_equilibrium_temperature)
    }

    /// Returns the common (lower, upper) temperature interval, K.
    pub(crate) fn common_temperature_interval(&self) -> (f64, f64) {
        (
            self.common_temperature_lower_k,
            self.common_temperature_upper_k,
        )
    }

    /// Builds a diagnostic message listing records that do not cover the case.
    fn unavailable_message(&self, reference: &NasaCeaH2O2HpReference) -> String {
        let unavailable = self
            .rows
            .iter()
            .filter(|row| {
                !row.supports_reactant_temperature || !row.supports_equilibrium_temperature
            })
            .map(|row| {
                format!(
                    "{} -> {} [{} , {}] K",
                    row.cea_identity,
                    row.component.label(),
                    row.temperature_lower_k,
                    row.temperature_upper_k
                )
            })
            .collect::<Vec<_>>();
        format!(
            "required temperatures {} K and {} K are outside exact local record intervals: {}; common interval is [{}, {}] K",
            reference.reactant_temperature_k,
            reference.equilibrium_temperature_k,
            unavailable.join("; "),
            self.common_temperature_lower_k,
            self.common_temperature_upper_k,
        )
    }
}

impl ResolvedNasaCeaH2O2HpFixture {
    /// Resolves only the eleven species printed by the NASA CEA example.
    ///
    /// The test stays offline by policy. A missing exact component is reported
    /// as `ValidationNotApplicable`; callers must not replace it by an alias,
    /// a different water polymorph, or a broader element search.
    pub(crate) fn resolve_offline(
        repository: Arc<ThermoRepository>,
        reference: &NasaCeaH2O2HpReference,
    ) -> Result<Self, ReactionExtentError> {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            declared_spec(),
            repository,
        )
        .map_err(|error| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_exact_component_inventory",
            message: format!("exact offline CEA universe could not resolve: {error}"),
        })?;
        if resolved.report().nist_fallback_enabled() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_cea_h2_o2_hp_nist",
                message: "frozen CEA benchmark must not enable NIST fallback".to_owned(),
            });
        }
        validate_resolved_mapping(&resolved)?;
        let preflight = preflight_temperature_domain(&resolved, reference)?;
        if !preflight.supports_required_temperatures() {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_temperature_domain",
                message: preflight.unavailable_message(reference),
            });
        }
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        let bounds = thermochemistry.temperature_bounds();
        for temperature in [
            reference.reactant_temperature_k,
            reference.equilibrium_temperature_k,
        ] {
            if !bounds.contains(temperature) {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_temperature_domain",
                    message: format!(
                        "common local interval [{}, {}] K excludes required {temperature} K",
                        bounds.lower(),
                        bounds.upper()
                    ),
                });
            }
        }
        let (matrix, molar_masses, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_molar_masses",
                    message: format!(
                        "local exact CEA records lack composition/molar mass: {error}"
                    ),
                })?;
        let structure = structural_evidence(&matrix)?;
        let component_inventory = component_inventory(
            &resolved,
            &thermochemistry,
            reference.equilibrium_temperature_k,
            &CEA_COMPONENT_MAPPINGS,
        )?;
        Ok(Self {
            resolved,
            thermochemistry,
            component_inventory,
            molar_masses_g_mol: molar_masses,
            structure,
        })
    }

    pub(crate) fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    pub(crate) fn thermochemistry(&self) -> &ResolvedThermochemistry {
        &self.thermochemistry
    }

    pub(crate) fn component_inventory(&self) -> &[NasaCeaComponentInventoryRow] {
        &self.component_inventory
    }

    pub(crate) fn structure(&self) -> NasaCeaMulticomponentStructure {
        self.structure
    }

    /// Reconstructs exactly one kilogram of reactants from CEA's O/F *mass*
    /// ratio. Local phase composition machinery supplies molar masses in
    /// g/mol, which are converted explicitly to kg/mol before division.
    pub(crate) fn reconstruct_reactants(
        &self,
        reference: &NasaCeaH2O2HpReference,
    ) -> Result<NasaCeaReactantReconstruction, ReactionExtentError> {
        reconstruct_reactants_from_resolved(
            &self.resolved,
            &self.thermochemistry,
            &self.molar_masses_g_mol,
            reference,
        )
    }

    /// Creates physical H2/O2 inventory with all product gases and condensed
    /// candidates absent. The numerical trace floor is selected later by the
    /// production solver; this value object remains physically truthful.
    pub(crate) fn initial_composition(
        &self,
        reactants: &NasaCeaReactantReconstruction,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        MultiphaseInitialComposition::from_sparse(
            &layout,
            vec![
                (
                    PhaseComponentId::new(PhaseId::new(Some(CEA_GAS_PHASE.to_owned())), "H2"),
                    reactants.hydrogen_moles,
                ),
                (
                    PhaseComponentId::new(PhaseId::new(Some(CEA_GAS_PHASE.to_owned())), "O2"),
                    reactants.oxygen_moles,
                ),
            ],
        )
    }
}

/// Resolves only identity and native coefficient intervals for the exact CEA
/// universe. It is intentionally separate from [`ResolvedNasaCeaH2O2HpFixture`]
/// so diagnostics can describe an unavailable record without evaluating it
/// beyond its declared domain.
pub(crate) fn preflight_nasa_cea_h2_o2_hp(
    repository: Arc<ThermoRepository>,
    reference: &NasaCeaH2O2HpReference,
) -> Result<NasaCeaH2O2HpPreflightReport, ReactionExtentError> {
    let resolved =
        SubstanceSystemFactory::resolve_phase_system_with_repository(declared_spec(), repository)
            .map_err(|error| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_exact_component_inventory",
            message: format!("exact offline CEA universe could not resolve: {error}"),
        })?;
    if resolved.report().nist_fallback_enabled() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "nasa_cea_h2_o2_hp_nist",
            message: "frozen CEA benchmark must not enable NIST fallback".to_owned(),
        });
    }
    validate_resolved_mapping(&resolved)?;
    preflight_temperature_domain(&resolved, reference)
}

/// Reconstructs CEA's one-kilogram H2/O2 input from exact local reactant
/// records only. This deliberately succeeds independently of whether an
/// optional condensed candidate later covers the final equilibrium temperature:
/// input reconstruction is not permission to solve a reduced product universe.
pub(crate) fn reconstruct_nasa_cea_h2_o2_reactants(
    repository: Arc<ThermoRepository>,
    reference: &NasaCeaH2O2HpReference,
) -> Result<NasaCeaReactantReconstruction, ReactionExtentError> {
    let resolved =
        SubstanceSystemFactory::resolve_phase_system_with_repository(reactant_spec(), repository)
            .map_err(|error| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_reactant_resolution",
            message: format!("exact local H2/O2 reactants could not resolve offline: {error}"),
        })?;
    let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
    if !thermochemistry
        .temperature_bounds()
        .contains(reference.reactant_temperature_k)
    {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_reactant_temperature_domain",
            message: format!(
                "local H2/O2 interval [{}, {}] K excludes {} K",
                thermochemistry.temperature_bounds().lower(),
                thermochemistry.temperature_bounds().upper(),
                reference.reactant_temperature_k,
            ),
        });
    }
    let (_, molar_masses, _) =
        element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
            .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_reactant_molar_masses",
                message: format!("local H2/O2 molar masses are unavailable: {error}"),
            })?;
    reconstruct_reactants_from_resolved(&resolved, &thermochemistry, &molar_masses, reference)
}

/// Loads the reviewed external CEA result using the generic read-only loader.
pub(crate) fn load_nasa_cea_h2_o2_hp_dataset()
-> Result<FrozenReferenceDataset<NasaCeaH2O2HpReference>, String> {
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_cea");
    FrozenReferenceDataset::load(
        directory.join("h2_o2_hp_scitech_2025.metadata.json"),
        directory.join("h2_o2_hp_scitech_2025.rows.json"),
    )
    .map_err(|error| error.to_string())
}

/// Builds the full eleven-component multiphase CEA universe spec.
fn declared_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
        (
            CEA_GAS_PHASE.to_owned(),
            CEA_GAS_SPECIES
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
        ),
        (CEA_LIQUID_PHASE.to_owned(), vec!["H2O(L)".to_owned()]),
        (CEA_ICE_PHASE.to_owned(), vec!["H2O(s)".to_owned()]),
    ])))
    .with_phase_natures(Some(HashMap::from([
        (CEA_GAS_PHASE.to_owned(), PhysicalState::Gas),
        (CEA_LIQUID_PHASE.to_owned(), PhysicalState::Liquid),
        (CEA_ICE_PHASE.to_owned(), PhysicalState::Solid),
    ])))
    .with_library_priorities(vec!["NASA_gas".to_owned(), "NASA_cond".to_owned()])
    .with_search_in_nist(false)
    .build()
    .expect("the static NASA CEA H2/O2 universe must be structurally valid")
}

/// Builds the nine-component gas-only benchmark universe spec.
fn gas_benchmark_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([(
        CEA_GAS_PHASE.to_owned(),
        CEA_GAS_SPECIES
            .iter()
            .map(|name| (*name).to_owned())
            .collect(),
    )])))
    .with_phase_natures(Some(HashMap::from([(
        CEA_GAS_PHASE.to_owned(),
        PhysicalState::Gas,
    )])))
    .with_library_priorities(vec!["NASA_gas".to_owned()])
    .with_search_in_nist(false)
    .build()
    .expect("the static nine-component NASA CEA gas universe must be valid")
}

/// Builds the two-species H2/O2 reactant-only universe spec.
fn reactant_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([(
        CEA_GAS_PHASE.to_owned(),
        vec!["H2".to_owned(), "O2".to_owned()],
    )])))
    .with_phase_natures(Some(HashMap::from([(
        CEA_GAS_PHASE.to_owned(),
        PhysicalState::Gas,
    )])))
    .with_library_priorities(vec!["NASA_gas".to_owned()])
    .with_search_in_nist(false)
    .build()
    .expect("the static NASA CEA H2/O2 reactant declaration must be valid")
}

/// Validates the resolved universe against the full static component mappings.
fn validate_resolved_mapping(resolved: &ResolvedPhaseSystem) -> Result<(), ReactionExtentError> {
    validate_mapping_subset(resolved, &CEA_COMPONENT_MAPPINGS)
}

/// Returns the component mappings that belong to the gas-only universe.
fn gas_mappings() -> &'static [NasaCeaComponentMapping] {
    &CEA_COMPONENT_MAPPINGS[..CEA_GAS_SPECIES.len()]
}

/// Verifies each mapping's phase, physical state, component identity, and exact
/// local Thermo provenance are preserved by the resolved universe.
fn validate_mapping_subset(
    resolved: &ResolvedPhaseSystem,
    mappings: &[NasaCeaComponentMapping],
) -> Result<(), ReactionExtentError> {
    for mapping in mappings {
        let phase = resolved
            .phase_specs()
            .iter()
            .find(|phase| phase.id().as_option().as_deref() == Some(mapping.phase))
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_phase_identity",
                message: format!("missing declared '{}' phase", mapping.phase),
            })?;
        if phase.physical_state() != mapping.physical_state
            || phase
                .components()
                .iter()
                .all(|name| name != mapping.local_species)
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_component_identity",
                message: format!(
                    "CEA '{}' did not retain pinned {}::{} identity",
                    mapping.cea_identity, mapping.phase, mapping.local_species
                ),
            });
        }
        let phase_report = resolved
            .report()
            .phases()
            .iter()
            .find(|report| report.phase().as_option().as_deref() == Some(mapping.phase))
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_lookup_provenance",
                message: format!("missing lookup report for '{}' phase", mapping.phase),
            })?;
        let selected = phase_report.search().rows().iter().any(|row| {
            row.property() == "Thermo"
                && row.substance() == mapping.local_species
                && row.library() == mapping.library
                && !row.record_key().trim().is_empty()
        });
        if !selected {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_lookup_provenance",
                message: format!(
                    "CEA '{}' lacks exact local Thermo provenance in {}",
                    mapping.cea_identity, mapping.library
                ),
            });
        }
    }
    Ok(())
}

/// Evaluates G/H/Cp and provenance for every mapped component at a temperature.
fn component_inventory(
    resolved: &ResolvedPhaseSystem,
    thermochemistry: &ResolvedThermochemistry,
    temperature_k: f64,
    mappings: &[NasaCeaComponentMapping],
) -> Result<Vec<NasaCeaComponentInventoryRow>, ReactionExtentError> {
    let gibbs = thermochemistry.evaluate_gibbs(temperature_k)?;
    let enthalpy = thermochemistry.evaluate_enthalpy(temperature_k)?;
    let heat_capacity = thermochemistry.evaluate_heat_capacity(temperature_k)?;
    mappings
        .iter()
        .map(|mapping| {
            let index = component_index(resolved, mapping.phase, mapping.local_species)?;
            let provenance = thermochemistry.provenance().get(index).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "thermochemistry provenance lacks CEA component {}",
                    mapping.cea_identity
                ))
            })?;
            inventory_row(
                mapping,
                provenance,
                thermochemistry,
                index,
                &gibbs,
                &enthalpy,
                &heat_capacity,
            )
        })
        .collect()
}

/// Builds the per-component record-interval preflight and the common interval.
fn preflight_temperature_domain(
    resolved: &ResolvedPhaseSystem,
    reference: &NasaCeaH2O2HpReference,
) -> Result<NasaCeaH2O2HpPreflightReport, ReactionExtentError> {
    let mut common_lower = f64::NEG_INFINITY;
    let mut common_upper = f64::INFINITY;
    let mut rows = Vec::with_capacity(CEA_COMPONENT_MAPPINGS.len());
    for mapping in CEA_COMPONENT_MAPPINGS {
        let phase_data = resolved
            .phase_data()
            .get(&Some(mapping.phase.to_owned()))
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_component_interval",
                message: format!("missing local phase payload for '{}'", mapping.phase),
            })?;
        let record = phase_data
            .get_search_result(mapping.local_species, WhatIsFound::Thermo)
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_component_interval",
                message: format!("missing Thermo record for '{}'", mapping.cea_identity),
            })?;
        let mut calculator = match record.calculator().cloned() {
            Some(CalculatorType::Thermo(calculator)) => calculator,
            Some(CalculatorType::Transport(_)) => {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_component_interval",
                    message: format!("'{}' resolved a transport record", mapping.cea_identity),
                });
            }
            None => {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_cea_h2_o2_hp_component_interval",
                    message: format!("'{}' has no calculator", mapping.cea_identity),
                });
            }
        };
        calculator.parse_coefficients().map_err(|error| {
            ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_component_interval",
                message: format!(
                    "'{}' coefficients cannot be parsed: {error}",
                    mapping.cea_identity
                ),
            }
        })?;
        let (lower, upper) = calculator.valid_temperature_interval().map_err(|error| {
            ReactionExtentError::ValidationNotApplicable {
                path: "nasa_cea_h2_o2_hp_component_interval",
                message: format!(
                    "'{}' has no valid temperature interval: {error}",
                    mapping.cea_identity
                ),
            }
        })?;
        common_lower = common_lower.max(lower);
        common_upper = common_upper.min(upper);
        rows.push(NasaCeaComponentDomainRow {
            cea_identity: mapping.cea_identity.to_owned(),
            component: PhaseComponentId::new(
                PhaseId::new(Some(mapping.phase.to_owned())),
                mapping.local_species,
            ),
            library: mapping.library.to_owned(),
            record_key: record.record_key().to_owned(),
            temperature_lower_k: lower,
            temperature_upper_k: upper,
            supports_reactant_temperature: (lower..=upper)
                .contains(&reference.reactant_temperature_k),
            supports_equilibrium_temperature: (lower..=upper)
                .contains(&reference.equilibrium_temperature_k),
        });
    }
    Ok(NasaCeaH2O2HpPreflightReport {
        rows,
        common_temperature_lower_k: common_lower,
        common_temperature_upper_k: common_upper,
    })
}

/// Assembles one inventory row from provenance and evaluated thermochemistry.
fn inventory_row(
    mapping: &NasaCeaComponentMapping,
    provenance: &ThermochemistryProvenance,
    thermochemistry: &ResolvedThermochemistry,
    index: usize,
    gibbs: &[f64],
    enthalpy: &[f64],
    heat_capacity: &[Option<f64>],
) -> Result<NasaCeaComponentInventoryRow, ReactionExtentError> {
    let heat_capacity_j_mol_k = heat_capacity
        .get(index)
        .and_then(|value| *value)
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_heat_capacity",
            message: format!("CEA '{}' lacks Cp capability", mapping.cea_identity),
        })?;
    let bounds = thermochemistry.temperature_bounds();
    Ok(NasaCeaComponentInventoryRow {
        cea_identity: mapping.cea_identity.to_owned(),
        component: provenance.component().clone(),
        library: provenance.library().to_owned(),
        record_key: provenance.record_key().to_owned(),
        physical_state: mapping.physical_state,
        temperature_lower_k: bounds.lower(),
        temperature_upper_k: bounds.upper(),
        gibbs_j_mol: gibbs[index],
        enthalpy_j_mol: enthalpy[index],
        heat_capacity_j_mol_k,
    })
}

/// Resolves the linear layout index of a phase-qualified component.
fn component_index(
    resolved: &ResolvedPhaseSystem,
    phase: &str,
    substance: &str,
) -> Result<usize, ReactionExtentError> {
    let component = PhaseComponentId::new(PhaseId::new(Some(phase.to_owned())), substance);
    resolved
        .layout()
        .components()
        .iter()
        .position(|candidate| candidate == &component)
        .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_component_layout",
            message: format!("missing exact component '{}'", component.label()),
        })
}

/// Reads a local molar mass and converts it from g/mol to kg/mol.
fn molar_mass_kg_mol(
    molar_masses_g_mol: &HashMap<String, f64>,
    substance: &str,
) -> Result<f64, ReactionExtentError> {
    let value = molar_masses_g_mol.get(substance).copied().ok_or_else(|| {
        ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_molar_mass",
            message: format!("exact local molar mass for '{substance}' is unavailable"),
        }
    })?;
    let converted = value / 1_000.0;
    if !converted.is_finite() || converted <= 0.0 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_molar_mass",
            message: format!("local molar mass for '{substance}' is invalid: {value} g/mol"),
        });
    }
    Ok(converted)
}

/// Reconstructs one kilogram of H2/O2 reactants from the CEA O/F mass ratio and
/// local records, computing the fixed-`P,H` enthalpy target.
fn reconstruct_reactants_from_resolved(
    resolved: &ResolvedPhaseSystem,
    thermochemistry: &ResolvedThermochemistry,
    molar_masses_g_mol: &HashMap<String, f64>,
    reference: &NasaCeaH2O2HpReference,
) -> Result<NasaCeaReactantReconstruction, ReactionExtentError> {
    let hydrogen_molar_mass_kg_mol = molar_mass_kg_mol(molar_masses_g_mol, "H2")?;
    let oxygen_molar_mass_kg_mol = molar_mass_kg_mol(molar_masses_g_mol, "O2")?;
    let total_mass_kg = 1.0;
    let hydrogen_mass_kg = total_mass_kg / (1.0 + reference.oxidizer_fuel_mass_ratio);
    let oxygen_mass_kg = total_mass_kg - hydrogen_mass_kg;
    let hydrogen_moles = hydrogen_mass_kg / hydrogen_molar_mass_kg_mol;
    let oxygen_moles = oxygen_mass_kg / oxygen_molar_mass_kg_mol;
    let enthalpy = thermochemistry.evaluate_enthalpy(reference.reactant_temperature_k)?;
    let h2 = component_index(resolved, CEA_GAS_PHASE, "H2")?;
    let o2 = component_index(resolved, CEA_GAS_PHASE, "O2")?;
    let target_enthalpy_j = hydrogen_moles * enthalpy[h2] + oxygen_moles * enthalpy[o2];
    if !target_enthalpy_j.is_finite() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "nasa_cea_h2_o2_hp_target_enthalpy",
            message: "local H2/O2 reactant enthalpy is non-finite".to_owned(),
        });
    }
    Ok(NasaCeaReactantReconstruction {
        total_mass_kg,
        hydrogen_mass_kg,
        oxygen_mass_kg,
        hydrogen_molar_mass_kg_mol,
        oxygen_molar_mass_kg_mol,
        hydrogen_moles,
        oxygen_moles,
        hydrogen_atom_moles: 2.0 * hydrogen_moles,
        oxygen_atom_moles: 2.0 * oxygen_moles,
        target_enthalpy_j,
    })
}

/// Computes component/element geometry and the reaction dimension via SVD.
fn structural_evidence(
    element_composition: &DMatrix<f64>,
) -> Result<NasaCeaMulticomponentStructure, ReactionExtentError> {
    let component_count = element_composition.nrows();
    let element_count = element_composition.ncols();
    let element_rank = SVD::new(element_composition.clone(), false, false)
        .singular_values
        .iter()
        .filter(|&&value| value > 1.0e-12)
        .count();
    if component_count == 0 || element_count == 0 || element_rank == 0 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_cea_h2_o2_hp_structure",
            message: format!(
                "invalid component/element geometry: {component_count} components, {element_count} elements, rank {element_rank}"
            ),
        });
    }
    Ok(NasaCeaMulticomponentStructure {
        component_count,
        element_count,
        element_rank,
        reaction_dimension: component_count.saturating_sub(element_rank),
    })
}
