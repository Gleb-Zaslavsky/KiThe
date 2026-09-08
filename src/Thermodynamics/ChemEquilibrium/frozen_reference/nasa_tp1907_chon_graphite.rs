//! NASA TP-1907 Table 11.3E CHON + graphite I5 benchmark adapter.
//!
//! This module owns a reviewed reduced local universe and source-facing
//! evidence only.  It never implements a special carbon solver: all accepted
//! states come from the ordinary resolved `P,T` phase-control workflow.

use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

use nalgebra::{DMatrix, linalg::SVD};

use crate::Kinetics::molmass::parse_formula;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::additive_total_enthalpy;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
    MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_thermochemistry::ResolvedThermochemistry;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, NasaTp1907MultiphaseReference,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::User_PhaseOrSolution::{
    ResolvedPhaseSystem, SubstanceSystemFactory, SubstanceSystemSpec, SubstanceSystemSpecBuilder,
    SubstancesContainer, element_composition_and_molar_mass,
};
use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
use crate::Thermodynamics::physical_state::PhysicalState;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// Semantic gas phase used by the local production problem.
pub(crate) const TP1907_GAS_PHASE: &str = "gas";
/// Pure solid graphite candidate used by the local production problem.
pub(crate) const TP1907_GRAPHITE_PHASE: &str = "graphite";
/// Activity reference pressure of local NASA records.
pub(crate) const TP1907_REFERENCE_PRESSURE_PA: f64 = 101_325.0;

/// Reviewed reduced gas universe. It is deliberately broader than the columns
/// printed in the selected rows, but deliberately not presented as NASA's
/// original 55-species universe.
pub(crate) const TP1907_LOCAL_GAS_SPECIES: [&str; 17] = [
    "Ar", "CH4", "CO", "CO2", "H", "HO2", "H2", "H2O", "N", "NH3", "NO", "NO2", "N2", "N2O", "O",
    "OH", "O2",
];

/// Structural evidence for the concrete local reaction space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct NasaTp1907Structure {
    /// Number of species in the local reaction universe.
    pub(crate) component_count: usize,
    /// Number of distinct chemical elements present.
    pub(crate) element_count: usize,
    /// Rank of the element-composition matrix.
    pub(crate) element_rank: usize,
    /// Dimension of the reaction space (`component_count - element_rank`).
    pub(crate) reaction_dimension: usize,
}

/// Local record and interval provenance captured before any solve.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaTp1907PreflightRow {
    /// Exact local phase-qualified component.
    pub(crate) component: PhaseComponentId,
    /// Library that supplied the exact Thermo record.
    pub(crate) library: String,
    /// Repository record key backing the coefficients.
    pub(crate) record_key: String,
    /// Lower bound of the common local temperature interval, K.
    pub(crate) lower_temperature_k: f64,
    /// Upper bound of the common local temperature interval, K.
    pub(crate) upper_temperature_k: f64,
}

/// One identity-aligned source/local composition comparison value.
///
/// Water ice and liquid are deliberately represented as excluded source-zero
/// values here because their local records do not enter this 700--720 K
/// graphite benchmark. They are not silently converted into gas components.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaTp1907CompositionComparisonRow {
    /// Source identity (gas species or condensed phase).
    pub(crate) identity: String,
    /// System-total mole fraction published by NASA.
    pub(crate) source_system_fraction: f64,
    /// Local system-total mole fraction, or `None` when excluded.
    pub(crate) kithe_system_fraction: Option<f64>,
    /// Absolute error `local - source`, or `None` when excluded.
    pub(crate) absolute_error: Option<f64>,
    /// Relative error, or `None` when the source fraction is not positive.
    pub(crate) relative_error: Option<f64>,
    /// Log10 difference, or `None` when not computable.
    pub(crate) delta_log10: Option<f64>,
    /// Human-readable reason for exclusion, or `None` when included.
    pub(crate) exclusion_reason: Option<String>,
}

/// Resolved immutable side of the NASA Table 11.3E benchmark.
#[derive(Clone)]
pub(crate) struct ResolvedNasaTp1907ChonGraphiteFixture {
    resolved: ResolvedPhaseSystem,
    thermochemistry: ResolvedThermochemistry,
    structure: NasaTp1907Structure,
    preflight: Vec<NasaTp1907PreflightRow>,
}

impl ResolvedNasaTp1907ChonGraphiteFixture {
    /// Resolves the reviewed local NASA gas plus NASA condensed graphite set.
    pub(crate) fn resolve_offline(
        repository: Arc<ThermoRepository>,
    ) -> Result<Self, ReactionExtentError> {
        let resolved =
            SubstanceSystemFactory::resolve_phase_system_with_repository(tp1907_spec(), repository)
                .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_tp1907_chon_graphite_resolution",
                    message: format!(
                        "reviewed CHON + graphite universe did not resolve offline: {error}"
                    ),
                })?;
        if resolved.report().nist_fallback_enabled() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_tp1907_chon_graphite_nist",
                message: "frozen NASA TP-1907 evidence must not enable NIST fallback".to_owned(),
            });
        }
        validate_local_contract(&resolved)?;
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)?;
        let (element_matrix, _, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .map_err(|error| ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_tp1907_chon_graphite_elements",
                    message: format!("local elemental composition is unavailable: {error}"),
                })?;
        let structure = structure_from_matrix(&element_matrix)?;
        if structure.component_count != TP1907_LOCAL_GAS_SPECIES.len() + 1
            || structure.element_count != 5
            || structure.element_rank != 5
            || structure.reaction_dimension != structure.component_count - 5
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_tp1907_chon_graphite_structure",
                message: format!(
                    "expected 18 components / 5 elements / rank 5 / 13 reactions, got {}/{}/{}/{}",
                    structure.component_count,
                    structure.element_count,
                    structure.element_rank,
                    structure.reaction_dimension,
                ),
            });
        }
        let preflight = preflight_rows(&resolved, &thermochemistry)?;
        Ok(Self {
            resolved,
            thermochemistry,
            structure,
            preflight,
        })
    }

    pub(crate) fn resolved(&self) -> &ResolvedPhaseSystem {
        &self.resolved
    }

    /// Returns the immutable local thermochemistry that belongs to this exact
    /// phase-qualified fixture.  P,H evidence must use this bundle rather than
    /// rebuilding one from a separately resolved system.
    pub(crate) fn thermochemistry(&self) -> &ResolvedThermochemistry {
        &self.thermochemistry
    }

    /// Computes the physical mass of the executable closed inventory in grams.
    ///
    /// This is intentionally derived from the ordinary `CH4/CO2/O2/N2/Ar`
    /// basis actually passed to the solver.  The value is also checked against
    /// the conceptual `CH2 + dry air` source construction, so a future change
    /// cannot silently alter the `J/g -> J` conversion used by P,H fixtures.
    pub(crate) fn initial_inventory_mass_g(&self) -> Result<f64, ReactionExtentError> {
        executable_inventory_mass_g()
    }

    /// Evaluates the currently supported additive total enthalpy of an
    /// accepted local state.  The helper is a read-only diagnostic boundary;
    /// it neither updates a thermochemical cache nor changes solver state.
    pub(crate) fn total_enthalpy_j(
        &self,
        solution: &MultiphaseEquilibriumSolution,
    ) -> Result<f64, ReactionExtentError> {
        let enthalpy = self
            .thermochemistry
            .evaluate_enthalpy(solution.conditions().temperature())?;
        additive_total_enthalpy(solution.component_moles(), &enthalpy)
    }
    pub(crate) fn structure(&self) -> NasaTp1907Structure {
        self.structure
    }
    pub(crate) fn preflight(&self) -> &[NasaTp1907PreflightRow] {
        &self.preflight
    }

    /// Builds a gas-only physical inventory. Graphite is intentionally absent
    /// so both source rows exercise production phase appearance or stability.
    pub(crate) fn initial_composition(
        &self,
    ) -> Result<MultiphaseInitialComposition, ReactionExtentError> {
        let feed = reconstructed_source_feed()?;
        let layout = MultiphaseEquilibriumLayout::new(self.resolved.phase_specs().to_vec())?;
        MultiphaseInitialComposition::from_sparse(
            &layout,
            vec![
                (gas_component("CH4"), 0.5),
                (gas_component("CO2"), 0.5 + feed.co2_moles),
                (gas_component("O2"), feed.o2_moles - 0.5),
                (gas_component("N2"), feed.n2_moles),
                (gas_component("Ar"), feed.ar_moles),
            ],
        )
    }

    pub(crate) fn conditions(
        &self,
        reference: &NasaTp1907MultiphaseReference,
    ) -> Result<EquilibriumConditions, ReactionExtentError> {
        self.conditions_at(reference.temperature_k, reference.pressure_pa)
    }

    /// Builds validated local conditions for a source row or an internal TPD
    /// probe without broadening the frozen external dataset.
    pub(crate) fn conditions_at(
        &self,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<EquilibriumConditions, ReactionExtentError> {
        if !self
            .thermochemistry
            .temperature_bounds()
            .contains(temperature_k)
        {
            let bounds = self.thermochemistry.temperature_bounds();
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nasa_tp1907_chon_graphite_temperature_domain",
                message: format!(
                    "common local interval [{}, {}] K excludes {} K",
                    bounds.lower(),
                    bounds.upper(),
                    temperature_k
                ),
            });
        }
        EquilibriumConditions::new(temperature_k, pressure_pa, TP1907_REFERENCE_PRESSURE_PA)
    }

    /// Aligns printed NASA values with accepted local system fractions.
    ///
    /// The production snapshot's `mole_fraction_for` is phase-local; it would
    /// return one for a pure graphite phase. This method therefore derives
    /// `n_i / sum(n_all_phases)` explicitly before comparing the source's
    /// heterogeneous system-total fractions.
    pub(crate) fn compare_system_composition(
        &self,
        reference: &NasaTp1907MultiphaseReference,
        solution: &MultiphaseEquilibriumSolution,
    ) -> Result<Vec<NasaTp1907CompositionComparisonRow>, ReactionExtentError> {
        let system_total_moles = solution.component_moles().iter().sum::<f64>();
        if !system_total_moles.is_finite() || system_total_moles <= 0.0 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_tp1907_chon_graphite_system_total",
                message: "accepted solution has no finite positive total physical amount"
                    .to_owned(),
            });
        }
        let mut rows =
            Vec::with_capacity(reference.gas_species.len() + reference.condensed_species.len());
        for external in &reference.gas_species {
            let component = gas_component(&external.species);
            let local = solution.moles_for(&component).ok_or_else(|| {
                ReactionExtentError::DimensionMismatch(format!(
                    "accepted solution lacks NASA TP-1907 gas component {}",
                    component.label()
                ))
            })? / system_total_moles;
            rows.push(comparison_row(
                external.species.clone(),
                external.system_mole_fraction,
                Some(local),
                None,
            ));
        }
        for external in &reference.condensed_species {
            if external.phase == "C(gr)" {
                let local = solution.moles_for(&graphite_component()).ok_or_else(|| {
                    ReactionExtentError::DimensionMismatch(
                        "accepted solution lacks graphite component".to_owned(),
                    )
                })? / system_total_moles;
                rows.push(comparison_row(
                    external.phase.clone(),
                    external.system_mole_fraction,
                    Some(local),
                    None,
                ));
            } else {
                rows.push(comparison_row(
                    external.phase.clone(),
                    external.system_mole_fraction,
                    None,
                    Some(
                        "published zero condensed water is outside the reviewed local graphite universe; no extrapolation was performed".to_owned(),
                    ),
                ));
            }
        }
        Ok(rows)
    }
}

/// Source-faithful NASA TP-1906 dry-air molar basis, normalized to one `O2`.
///
/// TP-1907 refers to the dry-air convention defined by TP-1906. Atmospheric
/// `CO2` is part of that convention and must remain in the element inventory:
/// treating the printed F/A as an equation for inferred `N2` loses this carbon
/// and makes a near-boundary graphite comparison ambiguous.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NasaTp1906DryAir {
    /// Moles of N2 per mole of O2 in dry air.
    pub(crate) n2_per_o2: f64,
    /// Moles of Ar per mole of O2 in dry air.
    pub(crate) ar_per_o2: f64,
    /// Moles of CO2 per mole of O2 in dry air (part of the element inventory).
    pub(crate) co2_per_o2: f64,
}

impl NasaTp1906DryAir {
    pub(crate) const fn source_faithful() -> Self {
        Self {
            n2_per_o2: 3.727_587,
            ar_per_o2: 0.044_706_8,
            co2_per_o2: 0.001_522_8,
        }
    }
}

/// Explicit, reproducible TP-1906/1907 feed reconstructed for the local
/// executable basis.
///
/// The source publishes an atom-ratio model fuel rather than a molecular
/// compound. This fixture uses one carbon mole-atom of conceptual `CH2`, then
/// constructs its oxidizer with the ordinary published equivalence ratio and
/// exact TP-1906 dry-air ratios. The executable basis below is ordinary
/// molecules only; `CH2` is never offered to the thermochemical repository.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct NasaTp1907ReconstructedFeed {
    /// Fuel hydrogen-to-carbon atom ratio (conceptual `CH2` fuel).
    pub(crate) h_to_c: f64,
    /// Ordinary equivalence ratio used to size the oxidizer.
    pub(crate) equivalence_ratio: f64,
    /// Oxidizer scaling factor for one carbon mole-atom of fuel.
    pub(crate) alpha: f64,
    /// Reconstructed O2 moles (per one carbon mole-atom).
    pub(crate) o2_moles: f64,
    /// Reconstructed N2 moles from the dry-air ratio.
    pub(crate) n2_moles: f64,
    /// Reconstructed Ar moles from the dry-air ratio.
    pub(crate) ar_moles: f64,
    /// Reconstructed CO2 moles from the dry-air ratio.
    pub(crate) co2_moles: f64,
}

const MODEL_FUEL_MOLAR_MASS_G_PER_MOL: f64 = 12.011 + 2.0 * 1.008;
const CH4_MOLAR_MASS_G_PER_MOL: f64 = 12.011 + 4.0 * 1.008;
const O2_MOLAR_MASS_G_PER_MOL: f64 = 2.0 * 15.999;
const N2_MOLAR_MASS_G_PER_MOL: f64 = 2.0 * 14.007;
const AR_MOLAR_MASS_G_PER_MOL: f64 = 39.948;
const CO2_MOLAR_MASS_G_PER_MOL: f64 = 12.011 + 2.0 * 15.999;
const TP1907_H_TO_C: f64 = 2.0;
const TP1907_EQUIVALENCE_RATIO: f64 = 1.25;
const TP1907_PUBLISHED_FUEL_AIR_RATIO: f64 = 0.084_535;
const TP1907_PUBLISHED_CHEMICAL_EQUIVALENCE_RATIO: f64 = 1.249_6;

impl NasaTp1907ReconstructedFeed {
    /// Fuel mass divided by the complete TP-1906 dry-air mass.
    pub(crate) fn fuel_air_mass_ratio(self) -> f64 {
        MODEL_FUEL_MOLAR_MASS_G_PER_MOL
            / (self.o2_moles * O2_MOLAR_MASS_G_PER_MOL
                + self.n2_moles * N2_MOLAR_MASS_G_PER_MOL
                + self.ar_moles * AR_MOLAR_MASS_G_PER_MOL
                + self.co2_moles * CO2_MOLAR_MASS_G_PER_MOL)
    }

    /// Chemical-equivalence diagnostic on the source oxygen-reactant basis.
    ///
    /// Atmospheric `CO2` belongs to dry-air mass and elemental inventory, but
    /// is not counted as fresh oxidizing `O2`. The result is compared to the
    /// independently rounded TP-1907 value and never selects feed coefficients.
    pub(crate) fn chemical_equivalence_ratio(self) -> f64 {
        let stoichiometric_o2 = (4.0 + self.h_to_c) * 0.25;
        stoichiometric_o2 / self.o2_moles
    }

    /// Elemental inventory of the conceptual model-fuel source basis.
    pub(crate) fn source_element_totals(self) -> BTreeMap<String, f64> {
        BTreeMap::from([
            ("C".to_owned(), 1.0 + self.co2_moles),
            ("H".to_owned(), self.h_to_c),
            ("O".to_owned(), 2.0 * (self.o2_moles + self.co2_moles)),
            ("N".to_owned(), 2.0 * self.n2_moles),
            ("Ar".to_owned(), self.ar_moles),
        ])
    }
}

/// Reconstructs the closed C/H/O/N/Ar inventory from the TP-1906 dry-air
/// formula and TP-1907's ordinary equivalence ratio.
pub(crate) fn reconstructed_source_feed() -> Result<NasaTp1907ReconstructedFeed, ReactionExtentError>
{
    let dry_air = NasaTp1906DryAir::source_faithful();
    let alpha = (4.0 + TP1907_H_TO_C) / (4.0 * TP1907_EQUIVALENCE_RATIO);
    let feed = NasaTp1907ReconstructedFeed {
        h_to_c: TP1907_H_TO_C,
        equivalence_ratio: TP1907_EQUIVALENCE_RATIO,
        alpha,
        o2_moles: alpha,
        n2_moles: alpha * dry_air.n2_per_o2,
        ar_moles: alpha * dry_air.ar_per_o2,
        co2_moles: alpha * dry_air.co2_per_o2,
    };
    if !feed.alpha.is_finite()
        || feed.alpha <= 0.0
        || !feed.o2_moles.is_finite()
        || feed.o2_moles <= 0.5
        || !feed.n2_moles.is_finite()
        || feed.n2_moles <= 0.0
        || !feed.ar_moles.is_finite()
        || feed.ar_moles <= 0.0
        || !feed.co2_moles.is_finite()
        || feed.co2_moles <= 0.0
    {
        return Err(ReactionExtentError::InvalidProblem {
            field: "nasa_tp1907_reconstructed_feed",
            message: "TP-1906/1907 source formula produced an invalid dry-air feed".to_owned(),
        });
    }
    Ok(feed)
}

/// Verifies that the conceptual model-fuel basis and executable ordinary-molecule basis preserve every element.
pub(crate) fn verify_element_equivalent_feed() -> Result<BTreeMap<String, f64>, ReactionExtentError>
{
    let feed = reconstructed_source_feed()?;
    let source = feed.source_element_totals();
    let local = element_totals_from_molecular_basis([
        ("CH4", 0.5),
        ("CO2", 0.5 + feed.co2_moles),
        ("O2", feed.o2_moles - 0.5),
        ("N2", feed.n2_moles),
        ("Ar", feed.ar_moles),
    ])?;
    let equivalent = source.len() == local.len()
        && source.iter().all(|(element, source_amount)| {
            local
                .get(element)
                .is_some_and(|local_amount| (source_amount - local_amount).abs() <= 1.0e-12)
        });
    if !equivalent {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_element_equivalence",
            message: format!("source={source:?}, executable={local:?}"),
        });
    }
    Ok(source)
}

/// Returns the exact mass of the ordinary-molecule inventory used by
/// [`ResolvedNasaTp1907ChonGraphiteFixture::initial_composition`].
///
/// A TP-1906 heterogeneous specific enthalpy is extensive only after this
/// conversion.  Keeping it beside the feed reconstruction makes the chosen
/// molecular-weight convention visible and prevents a guessed total mass from
/// entering a frozen P,H comparison.
pub(crate) fn executable_inventory_mass_g() -> Result<f64, ReactionExtentError> {
    let feed = reconstructed_source_feed()?;
    let local_mass = 0.5 * CH4_MOLAR_MASS_G_PER_MOL
        + (0.5 + feed.co2_moles) * CO2_MOLAR_MASS_G_PER_MOL
        + (feed.o2_moles - 0.5) * O2_MOLAR_MASS_G_PER_MOL
        + feed.n2_moles * N2_MOLAR_MASS_G_PER_MOL
        + feed.ar_moles * AR_MOLAR_MASS_G_PER_MOL;
    let source_mass = MODEL_FUEL_MOLAR_MASS_G_PER_MOL
        + feed.o2_moles * O2_MOLAR_MASS_G_PER_MOL
        + feed.n2_moles * N2_MOLAR_MASS_G_PER_MOL
        + feed.ar_moles * AR_MOLAR_MASS_G_PER_MOL
        + feed.co2_moles * CO2_MOLAR_MASS_G_PER_MOL;
    if !local_mass.is_finite() || local_mass <= 0.0 {
        return Err(ReactionExtentError::InvalidProblem {
            field: "nasa_tp1907_chon_graphite_inventory_mass",
            message: "executable TP-1906/1907 inventory mass must be finite and positive"
                .to_owned(),
        });
    }
    if (local_mass - source_mass).abs() > 1.0e-12 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_inventory_mass",
            message: format!(
                "ordinary-molecule inventory mass {local_mass:e} g differs from source basis {source_mass:e} g"
            ),
        });
    }
    Ok(local_mass)
}

/// Validates rounded TP-1907 diagnostics without making either quantity an
/// input to the source reconstruction.
pub(crate) fn validate_published_feed_diagnostics(
    feed: NasaTp1907ReconstructedFeed,
) -> Result<(), ReactionExtentError> {
    let fuel_air_delta = (feed.fuel_air_mass_ratio() - TP1907_PUBLISHED_FUEL_AIR_RATIO).abs();
    if fuel_air_delta > 5.0e-7 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_fuel_air_ratio",
            message: format!(
                "source-faithful F/A={} differs from published {} by {}",
                feed.fuel_air_mass_ratio(),
                TP1907_PUBLISHED_FUEL_AIR_RATIO,
                fuel_air_delta
            ),
        });
    }
    let chemical_delta =
        (feed.chemical_equivalence_ratio() - TP1907_PUBLISHED_CHEMICAL_EQUIVALENCE_RATIO).abs();
    if chemical_delta > 5.0e-4 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_chemical_equivalence_ratio",
            message: format!(
                "oxygen-reactant chemical ER={} differs from published {} by {}",
                feed.chemical_equivalence_ratio(),
                TP1907_PUBLISHED_CHEMICAL_EQUIVALENCE_RATIO,
                chemical_delta
            ),
        });
    }
    Ok(())
}

/// Loads the read-only Table 11.3E evidence through the common frozen loader.
pub(crate) fn load_nasa_tp1907_chon_graphite_dataset()
-> Result<FrozenReferenceDataset<NasaTp1907MultiphaseReference>, String> {
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_tp1907");
    FrozenReferenceDataset::load(
        directory.join("chon_graphite_er125_1atm.metadata.json"),
        directory.join("chon_graphite_er125_1atm.rows.json"),
    )
    .map_err(|error| error.to_string())
}

/// Builds the static CHON + graphite multiphase universe spec for the solve.
fn tp1907_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
        (
            TP1907_GAS_PHASE.to_owned(),
            TP1907_LOCAL_GAS_SPECIES
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
        ),
        (TP1907_GRAPHITE_PHASE.to_owned(), vec!["C(gr)".to_owned()]),
    ])))
    .with_phase_natures(Some(HashMap::from([
        (TP1907_GAS_PHASE.to_owned(), PhysicalState::Gas),
        (TP1907_GRAPHITE_PHASE.to_owned(), PhysicalState::Solid),
    ])))
    .with_library_priorities(vec!["NASA_gas".to_owned(), "NASA_cond".to_owned()])
    .with_search_in_nist(false)
    .build()
    .expect("static NASA TP-1907 CHON + graphite declaration must be valid")
}

/// Verifies the resolved system preserves the exact gas and graphite declarations.
fn validate_local_contract(resolved: &ResolvedPhaseSystem) -> Result<(), ReactionExtentError> {
    let gas = resolved
        .phase_specs()
        .iter()
        .find(|phase| phase.id().as_option().as_deref() == Some(TP1907_GAS_PHASE));
    let graphite = resolved
        .phase_specs()
        .iter()
        .find(|phase| phase.id().as_option().as_deref() == Some(TP1907_GRAPHITE_PHASE));
    let exact_gas = gas.is_some_and(|phase| {
        phase.physical_state() == PhysicalState::Gas
            && phase
                .components()
                .iter()
                .map(String::as_str)
                .eq(TP1907_LOCAL_GAS_SPECIES)
    });
    let exact_graphite = graphite.is_some_and(|phase| {
        phase.physical_state() == PhysicalState::Solid
            && phase.components().iter().map(String::as_str).eq(["C(gr)"])
    });
    if !exact_gas || !exact_graphite {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_universe",
            message: "resolved phase system did not preserve exact gas plus graphite declarations"
                .to_owned(),
        });
    }
    Ok(())
}

/// Builds the provenance and interval preflight table for every component.
fn preflight_rows(
    resolved: &ResolvedPhaseSystem,
    thermochemistry: &ResolvedThermochemistry,
) -> Result<Vec<NasaTp1907PreflightRow>, ReactionExtentError> {
    thermochemistry
        .provenance()
        .iter()
        .map(|row| {
            let expected_library =
                if row.component().phase.as_option().as_deref() == Some(TP1907_GAS_PHASE) {
                    "NASA_gas"
                } else {
                    "NASA_cond"
                };
            if row.library() != expected_library || row.record_key().trim().is_empty() {
                return Err(ReactionExtentError::ValidationNotApplicable {
                    path: "nasa_tp1907_chon_graphite_provenance",
                    message: format!(
                        "{} lacks exact {expected_library} provenance",
                        row.component().label()
                    ),
                });
            }
            if !resolved
                .layout()
                .components()
                .iter()
                .any(|component| component == row.component())
            {
                return Err(ReactionExtentError::DimensionMismatch(format!(
                    "provenance component {} is absent from layout",
                    row.component().label()
                )));
            }
            // `ResolvedThermochemistry` publishes an intersection interval.  It
            // is the actual executable contract for this all-component fixture.
            let bounds = thermochemistry.temperature_bounds();
            Ok(NasaTp1907PreflightRow {
                component: row.component().clone(),
                library: row.library().to_owned(),
                record_key: row.record_key().to_owned(),
                lower_temperature_k: bounds.lower(),
                upper_temperature_k: bounds.upper(),
            })
        })
        .collect()
}

/// Computes component/element geometry and the reaction dimension via SVD.
fn structure_from_matrix(
    matrix: &DMatrix<f64>,
) -> Result<NasaTp1907Structure, ReactionExtentError> {
    let rank = SVD::new(matrix.clone(), false, false)
        .singular_values
        .iter()
        .filter(|&&value| value > 1.0e-12)
        .count();
    if matrix.nrows() == 0 || matrix.ncols() == 0 || rank == 0 {
        return Err(ReactionExtentError::ValidationNotApplicable {
            path: "nasa_tp1907_chon_graphite_structure",
            message: "empty or rank-zero element matrix".to_owned(),
        });
    }
    Ok(NasaTp1907Structure {
        component_count: matrix.nrows(),
        element_count: matrix.ncols(),
        element_rank: rank,
        reaction_dimension: matrix.nrows().saturating_sub(rank),
    })
}

/// Reconstructs an elemental inventory from arbitrary molecular basis rows.
fn element_totals_from_molecular_basis(
    amounts: impl IntoIterator<Item = (&'static str, f64)>,
) -> Result<BTreeMap<String, f64>, ReactionExtentError> {
    let mut totals = BTreeMap::new();
    for (formula, amount) in amounts {
        let composition = parse_formula(formula.to_owned(), None).map_err(|error| {
            ReactionExtentError::InvalidProblem {
                field: "nasa_tp1907_chon_graphite_formula",
                message: format!("cannot parse {formula}: {error}"),
            }
        })?;
        for (element, count) in composition {
            *totals.entry(element).or_insert(0.0) += amount * count as f64;
        }
    }
    Ok(totals)
}

/// Constructs a phase-qualified component id in the TP-1907 gas phase.
pub(crate) fn gas_component(species: &str) -> PhaseComponentId {
    PhaseComponentId::new(PhaseId::new(Some(TP1907_GAS_PHASE.to_owned())), species)
}

/// Constructs the phase-qualified graphite (`C(gr)`) component id.
pub(crate) fn graphite_component() -> PhaseComponentId {
    PhaseComponentId::new(
        PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned())),
        "C(gr)",
    )
}

/// Assembles one composition comparison row with derived error metrics.
fn comparison_row(
    identity: String,
    source_system_fraction: f64,
    kithe_system_fraction: Option<f64>,
    exclusion_reason: Option<String>,
) -> NasaTp1907CompositionComparisonRow {
    let absolute_error = kithe_system_fraction.map(|local| local - source_system_fraction);
    let relative_error = absolute_error
        .filter(|_| source_system_fraction > 0.0)
        .map(|error| error / source_system_fraction);
    let delta_log10 = kithe_system_fraction
        .filter(|local| *local > 0.0 && source_system_fraction > 0.0)
        .map(|local| local.log10() - source_system_fraction.log10());
    NasaTp1907CompositionComparisonRow {
        identity,
        source_system_fraction,
        kithe_system_fraction,
        absolute_error,
        relative_error,
        delta_log10,
        exclusion_reason,
    }
}
