//! Frozen NASA TP-1906/1907 CHON + graphite `P,H` benchmark adapter.
//!
//! TP-1906 owns the external heterogeneous-mixture specific enthalpy, while
//! TP-1907 owns the separate composition and graphite-topology evidence.  This
//! module joins only source-compatible rows and converts the former to an
//! extensive target through the exact local executable inventory.  It does not
//! add a benchmark-specific equilibrium algorithm or enthalpy correction.

use std::path::PathBuf;
use std::sync::Arc;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
    EquilibriumConstraint, TemperatureBounds,
};
use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedPhaseEnthalpyRequest;
use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, NasaTp1906HeterogeneousEnthalpyReference, NasaTp1907MultiphaseReference,
};
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nasa_tp1907_chon_graphite::{
    ResolvedNasaTp1907ChonGraphiteFixture, TP1907_GRAPHITE_PHASE,
    load_nasa_tp1907_chon_graphite_dataset,
};
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
    PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::thermo_lib_api::ThermoRepository;

/// One source-compatible TP-1906 enthalpy and TP-1907 composition pair.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaTp1906Tp1907PhCase {
    /// TP-1906 heterogeneous-equilibrium enthalpy source row.
    pub(crate) enthalpy: NasaTp1906HeterogeneousEnthalpyReference,
    /// Matched TP-1907 multiphase composition source row.
    pub(crate) composition: NasaTp1907MultiphaseReference,
    /// Exact mass of the ordinary-molecule closed inventory used by KiThe.
    pub(crate) inventory_mass_g: f64,
    /// Extensive production input derived only from `specific_enthalpy_j_g`
    /// and `inventory_mass_g`.
    pub(crate) target_enthalpy_j: f64,
}

impl NasaTp1906Tp1907PhCase {
    /// Source expectation for the graphite phase, owned by TP-1907 rather
    /// than inferred from the TP-1906 enthalpy value.
    pub(crate) fn expected_graphite_active(&self) -> bool {
        self.composition
            .condensed_species
            .iter()
            .find(|entry| entry.phase == "C(gr)")
            .is_some_and(|entry| entry.system_mole_fraction > 0.0)
    }

    /// Reports the original specific source value without concealing its unit.
    pub(crate) fn source_specific_enthalpy_j_g(&self) -> f64 {
        self.enthalpy.specific_enthalpy_j_g
    }
}

/// Reference-state preflight for one external temperature.
///
/// It deliberately describes a discrepancy instead of classifying it as a
/// solver error.  A later review may establish a physical conversion, but no
/// empirical offset may be introduced from these rows.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NasaTp1906EnthalpyPreflight {
    /// Temperature of the joined source pair, K.
    pub(crate) temperature_k: f64,
    /// TP-1906 published specific enthalpy, J/g.
    pub(crate) source_specific_enthalpy_j_g: f64,
    /// Local P,T witness specific enthalpy, J/g.
    pub(crate) local_specific_enthalpy_j_g: f64,
    /// Difference `local - source` specific enthalpy, J/g.
    pub(crate) delta_specific_enthalpy_j_g: f64,
    /// Extensive local total enthalpy, J.
    pub(crate) local_total_enthalpy_j: f64,
    /// Exact inventory mass used for the `J/g -> J` conversion, g.
    pub(crate) inventory_mass_g: f64,
}

/// Immutable local side of the external heterogeneous `P,H` benchmark.
#[derive(Clone)]
pub(crate) struct ResolvedNasaTp1906ChonGraphitePhFixture {
    tp1907: ResolvedNasaTp1907ChonGraphiteFixture,
    tp1906: FrozenReferenceDataset<NasaTp1906HeterogeneousEnthalpyReference>,
    tp1907_rows: FrozenReferenceDataset<NasaTp1907MultiphaseReference>,
    inventory_mass_g: f64,
}

impl ResolvedNasaTp1906ChonGraphitePhFixture {
    /// Resolves the existing local TP-1907 universe and loads both independent
    /// frozen source tables.  No network fallback or mutable repository path
    /// belongs to this evidence adapter.
    pub(crate) fn resolve_offline(
        repository: Arc<ThermoRepository>,
    ) -> Result<Self, ReactionExtentError> {
        let tp1907 = ResolvedNasaTp1907ChonGraphiteFixture::resolve_offline(repository)?;
        let tp1906 = load_nasa_tp1906_chon_graphite_enthalpy_dataset().map_err(|message| {
            ReactionExtentError::ValidationNotApplicable {
                path: "nasa_tp1906_chon_graphite_enthalpy_load",
                message,
            }
        })?;
        let tp1907_rows = load_nasa_tp1907_chon_graphite_dataset().map_err(|message| {
            ReactionExtentError::ValidationNotApplicable {
                path: "nasa_tp1906_tp1907_composition_load",
                message,
            }
        })?;
        let inventory_mass_g = tp1907.initial_inventory_mass_g()?;
        Ok(Self {
            tp1907,
            tp1906,
            tp1907_rows,
            inventory_mass_g,
        })
    }

    /// Underlying source-faithful resolved phase system shared with the P,T
    /// TP-1907 benchmark.
    pub(crate) fn tp1907(&self) -> &ResolvedNasaTp1907ChonGraphiteFixture {
        &self.tp1907
    }

    /// Exact inventory mass used for every source `J/g -> J` conversion.
    pub(crate) fn inventory_mass_g(&self) -> f64 {
        self.inventory_mass_g
    }

    /// Builds semantic source joins without relying on row order.
    pub(crate) fn cases(&self) -> Result<Vec<NasaTp1906Tp1907PhCase>, ReactionExtentError> {
        self.tp1906
            .rows()
            .iter()
            .map(|enthalpy| {
                let composition = self
                    .tp1907_rows
                    .rows()
                    .iter()
                    .find(|composition| same_source_family(enthalpy, composition))
                    .cloned()
                    .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                        path: "nasa_tp1906_tp1907_semantic_join",
                        message: format!(
                            "no TP-1907 composition row matches TP-1906 enthalpy at {} K",
                            enthalpy.temperature_k
                        ),
                    })?;
                let target_enthalpy_j = enthalpy.specific_enthalpy_j_g * self.inventory_mass_g;
                if !target_enthalpy_j.is_finite() {
                    return Err(ReactionExtentError::InvalidProblem {
                        field: "nasa_tp1906_target_enthalpy",
                        message: format!(
                            "TP-1906 target at {} K became non-finite after J/g -> J conversion",
                            enthalpy.temperature_k
                        ),
                    });
                }
                Ok(NasaTp1906Tp1907PhCase {
                    enthalpy: enthalpy.clone(),
                    composition,
                    inventory_mass_g: self.inventory_mass_g,
                    target_enthalpy_j,
                })
            })
            .collect()
    }

    /// Solves an independent canonical P,T witness at a source temperature.
    /// This is the mandatory enthalpy-reference preflight, not a replacement
    /// for the later P,H request.
    pub(crate) fn solve_pt_witness(
        &self,
        temperature_k: f64,
    ) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
        solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                self.tp1907.resolved(),
                self.tp1907.conditions_at(temperature_k, 101_325.0)?,
                self.tp1907.initial_composition()?,
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
    }

    /// Evaluates one P,T witness against the separate TP-1906 source value.
    pub(crate) fn enthalpy_preflight(
        &self,
        case: &NasaTp1906Tp1907PhCase,
        solution: &MultiphaseEquilibriumSolution,
    ) -> Result<NasaTp1906EnthalpyPreflight, ReactionExtentError> {
        if (solution.conditions().temperature() - case.enthalpy.temperature_k).abs() > 1.0e-12 {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_tp1906_preflight_temperature",
                message: "P,T witness temperature does not match the joined source row".to_owned(),
            });
        }
        let local_total_enthalpy_j = self.tp1907.total_enthalpy_j(solution)?;
        let local_specific_enthalpy_j_g = local_total_enthalpy_j / self.inventory_mass_g;
        if !local_specific_enthalpy_j_g.is_finite() {
            return Err(ReactionExtentError::InvalidProblem {
                field: "nasa_tp1906_local_specific_enthalpy",
                message: "local P,T witness has non-finite specific enthalpy".to_owned(),
            });
        }
        Ok(NasaTp1906EnthalpyPreflight {
            temperature_k: case.enthalpy.temperature_k,
            source_specific_enthalpy_j_g: case.source_specific_enthalpy_j_g(),
            local_specific_enthalpy_j_g,
            delta_specific_enthalpy_j_g: local_specific_enthalpy_j_g
                - case.source_specific_enthalpy_j_g(),
            local_total_enthalpy_j,
            inventory_mass_g: self.inventory_mass_g,
        })
    }

    /// Builds one canonical production P,H request.  `initial_temperature_k`
    /// is solely a numerical seed chosen by the caller; it is never read from
    /// the TP-1906/1907 source row.
    pub(crate) fn ph_request(
        &self,
        case: &NasaTp1906Tp1907PhCase,
        initial_temperature_k: f64,
    ) -> Result<ResolvedPhaseEnthalpyRequest<'_>, ReactionExtentError> {
        let bounds = self.tp1907.thermochemistry().temperature_bounds();
        ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
            self.tp1907.resolved(),
            self.tp1907.initial_composition()?,
            EquilibriumConstraint::ph(
                case.enthalpy.pressure_pa,
                case.enthalpy.pressure_pa,
                case.target_enthalpy_j,
                initial_temperature_k,
            )?,
            TemperatureBounds::new(bounds.lower(), bounds.upper())?,
            self.tp1907.thermochemistry().clone(),
        )
        .map(|request| request.with_phase_control_policy(PhaseControlPolicy::default()))
    }

    /// Returns whether graphite is active in an accepted production result.
    pub(crate) fn graphite_is_active(&self, solution: &MultiphaseEquilibriumSolution) -> bool {
        matches!(
            solution.phase_status(&PhaseId::new(Some(TP1907_GRAPHITE_PHASE.to_owned()))),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        )
    }
}

/// Loads the TP-1906 enthalpy source without merging its provenance with
/// TP-1907 composition evidence.
pub(crate) fn load_nasa_tp1906_chon_graphite_enthalpy_dataset()
-> Result<FrozenReferenceDataset<NasaTp1906HeterogeneousEnthalpyReference>, String> {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_tp1906");
    FrozenReferenceDataset::load(
        directory.join("chon_graphite_er125_1atm_heterogeneous_enthalpy.metadata.json"),
        directory.join("chon_graphite_er125_1atm_heterogeneous_enthalpy.rows.json"),
    )
    .map_err(|error| error.to_string())
}

/// Returns whether a TP-1906 enthalpy row and a TP-1907 composition row share
/// the exact source physical family (temperature, pressure, and ratios), so
/// they can be joined as one benchmark case.
fn same_source_family(
    enthalpy: &NasaTp1906HeterogeneousEnthalpyReference,
    composition: &NasaTp1907MultiphaseReference,
) -> bool {
    enthalpy.temperature_k == composition.temperature_k
        && enthalpy.pressure_pa == composition.pressure_pa
        && enthalpy.fuel_h_to_c_atom_ratio == composition.fuel_h_to_c_atom_ratio
        && enthalpy.fuel_air_mass_ratio == composition.fuel_air_mass_ratio
        && enthalpy.equivalence_ratio == composition.equivalence_ratio
        && enthalpy.chemical_equivalence_ratio == composition.chemical_equivalence_ratio
        && enthalpy.dry_air == composition.dry_air
}
