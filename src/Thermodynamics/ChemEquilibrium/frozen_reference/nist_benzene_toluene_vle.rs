//! Independent I1 Raoult-law evidence for the frozen benzene/toluene P-x case.
//!
//! This module owns the I1 Raoult oracle and the I4 local-record declaration.
//! It intentionally knows nothing about Gibbs residuals, TPD, or active-set
//! mutation. Bubble pressure and vapor composition come only from frozen pure
//! endpoints, leaving the oracle independent from later I3 candidate tests.

use std::path::PathBuf;

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, NistBenzeneTolueneVleReference,
};
use crate::Thermodynamics::User_PhaseOrSolution::{
    PhaseModel, PhaseSpec, ResolvedPhaseSystem, SubstanceSystemFactory, SubstanceSystemSpec,
};
use crate::Thermodynamics::phase_layout::PhaseId;
use crate::Thermodynamics::physical_state::PhysicalState;

/// Isotherm temperature of the frozen benzene/toluene P-x dataset, K.
pub(crate) const BENZENE_TOLUENE_TEMPERATURE_K: f64 = 353.15;
/// Canonical gas-phase identity for this binary VLE fixture.
pub(crate) const BENZENE_TOLUENE_GAS_PHASE: &str = "gas";
/// Canonical liquid-solution identity for this binary VLE fixture.
pub(crate) const BENZENE_TOLUENE_LIQUID_PHASE: &str = "liquid";

const BENZENE_GAS: &str = "C6H6";
const TOLUENE_GAS: &str = "C7H8";
const BENZENE_LIQUID: &str = "C6H6(L)";
const TOLUENE_LIQUID: &str = "C7H8(L)";

/// Raoult-law result corresponding to one frozen P-x measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct BenzeneTolueneRaoultReference {
    /// Liquid-phase benzene mole fraction.
    pub(crate) liquid_benzene_mole_fraction: f64,
    /// Liquid-phase toluene mole fraction (`1 - x_benzene`).
    pub(crate) liquid_toluene_mole_fraction: f64,
    /// Independent Raoult bubble pressure, Pa.
    pub(crate) bubble_pressure_pa: f64,
    /// Vapor-phase benzene mole fraction.
    pub(crate) vapor_benzene_mole_fraction: f64,
    /// Vapor-phase toluene mole fraction.
    pub(crate) vapor_toluene_mole_fraction: f64,
}

/// Declares the exact local phase universe for the binary VLE story.
///
/// The liquid is explicitly an `IdealSolution`; merely labelling it liquid
/// would select the wrong `PureCondensed` activity contract. Lookup remains
/// offline and state-specific, with gas and liquid standard states supplied by
/// the local NASA gas and condensed libraries respectively.
pub(crate) fn benzene_toluene_vle_spec() -> SubstanceSystemSpec {
    SubstanceSystemSpec::from_phases(vec![
        PhaseSpec::ideal_gas(
            PhaseId::new(Some(BENZENE_TOLUENE_GAS_PHASE.to_owned())),
            vec![BENZENE_GAS.to_owned(), TOLUENE_GAS.to_owned()],
        )
        .expect("static benzene/toluene gas declaration must be valid"),
        PhaseSpec::ideal_solution(
            PhaseId::new(Some(BENZENE_TOLUENE_LIQUID_PHASE.to_owned())),
            vec![BENZENE_LIQUID.to_owned(), TOLUENE_LIQUID.to_owned()],
            PhysicalState::Liquid,
        )
        .expect("static benzene/toluene liquid declaration must be valid"),
    ])
    .expect("static benzene/toluene VLE declaration must be valid")
    .with_lookup_policy(
        vec!["NASA_gas".to_owned(), "NASA_cond".to_owned()],
        vec!["NASA_gas".to_owned(), "NASA_cond".to_owned()],
        None,
        false,
    )
}

/// Resolves the reviewed binary VLE universe from the local repository only.
///
/// This boundary intentionally rejects online NIST fallback. The frozen NIST
/// table is I5 evidence, while the thermochemistry used by the solver remains
/// independent I4 data from the local repository.
pub(crate) fn resolve_offline_benzene_toluene_vle()
-> Result<ResolvedPhaseSystem, ReactionExtentError> {
    let resolved =
        SubstanceSystemFactory::resolve_spec(benzene_toluene_vle_spec()).map_err(|error| {
            ReactionExtentError::ValidationNotApplicable {
                path: "nist_benzene_toluene_vle_resolution",
                message: format!(
                    "reviewed benzene/toluene VLE system did not resolve offline: {error}"
                ),
            }
        })?;
    if resolved.report().nist_fallback_enabled() {
        return Err(ReactionExtentError::InvalidProblem {
            field: "nist_benzene_toluene_vle_nist",
            message: "frozen NIST VLE evidence must not query NIST at test runtime".to_owned(),
        });
    }
    let expected = [
        (
            BENZENE_TOLUENE_GAS_PHASE,
            PhysicalState::Gas,
            PhaseModel::IdealGas,
            [BENZENE_GAS, TOLUENE_GAS],
        ),
        (
            BENZENE_TOLUENE_LIQUID_PHASE,
            PhysicalState::Liquid,
            PhaseModel::IdealSolution,
            [BENZENE_LIQUID, TOLUENE_LIQUID],
        ),
    ];
    for (phase_id, state, model, components) in expected {
        let phase = resolved
            .phase_specs()
            .iter()
            .find(|phase| phase.id().as_option().as_deref() == Some(phase_id))
            .ok_or_else(|| ReactionExtentError::ValidationNotApplicable {
                path: "nist_benzene_toluene_vle_phase_layout",
                message: format!("resolved system lacks declared {phase_id} phase"),
            })?;
        if phase.physical_state() != state
            || phase.model() != model
            || !phase.components().iter().map(String::as_str).eq(components)
        {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "nist_benzene_toluene_vle_phase_layout",
                message: format!(
                    "resolved {phase_id} phase does not preserve its declared state, model, and component order"
                ),
            });
        }
    }
    Ok(resolved)
}

/// Loads the reviewed offline ThermoML subset. No network or mutable cache is
/// involved in test execution.
pub(crate) fn load_nist_benzene_toluene_vle_dataset()
-> Result<FrozenReferenceDataset<NistBenzeneTolueneVleReference>, String> {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src")
        .join("Thermodynamics")
        .join("ChemEquilibrium")
        .join("frozen_reference")
        .join("data")
        .join("nist_thermoml");
    FrozenReferenceDataset::load(
        directory.join("benzene_toluene_vle_353_15k.metadata.json"),
        directory.join("benzene_toluene_vle_353_15k.rows.json"),
    )
    .map_err(|error| error.to_string())
}

/// Computes the independent ideal-liquid bubble and vapor composition route.
pub(crate) fn raoult_reference(
    dataset: &FrozenReferenceDataset<NistBenzeneTolueneVleReference>,
    row: &NistBenzeneTolueneVleReference,
) -> Result<BenzeneTolueneRaoultReference, String> {
    let pure_toluene = dataset
        .rows()
        .first()
        .ok_or_else(|| "benzene/toluene frozen P-x dataset has no toluene endpoint".to_owned())?;
    let pure_benzene = dataset
        .rows()
        .last()
        .ok_or_else(|| "benzene/toluene frozen P-x dataset has no benzene endpoint".to_owned())?;
    if row.temperature_k != BENZENE_TOLUENE_TEMPERATURE_K {
        return Err(format!(
            "benzene/toluene Raoult oracle only owns {BENZENE_TOLUENE_TEMPERATURE_K} K, got {} K",
            row.temperature_k
        ));
    }
    let x_benzene = row.liquid_benzene_mole_fraction;
    let x_toluene = 1.0 - x_benzene;
    let bubble_pressure_pa = x_benzene * pure_benzene.experimental_pressure_pa
        + x_toluene * pure_toluene.experimental_pressure_pa;
    if !bubble_pressure_pa.is_finite() || bubble_pressure_pa <= 0.0 {
        return Err("Raoult bubble pressure must be finite and positive".to_owned());
    }
    let y_benzene = x_benzene * pure_benzene.experimental_pressure_pa / bubble_pressure_pa;
    let y_toluene = x_toluene * pure_toluene.experimental_pressure_pa / bubble_pressure_pa;
    if !y_benzene.is_finite()
        || !y_toluene.is_finite()
        || y_benzene < 0.0
        || y_toluene < 0.0
        || (y_benzene + y_toluene - 1.0).abs() > 1.0e-12
    {
        return Err("Raoult vapor composition must be a finite binary simplex point".to_owned());
    }
    Ok(BenzeneTolueneRaoultReference {
        liquid_benzene_mole_fraction: x_benzene,
        liquid_toluene_mole_fraction: x_toluene,
        bubble_pressure_pa,
        vapor_benzene_mole_fraction: y_benzene,
        vapor_toluene_mole_fraction: y_toluene,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        BENZENE_TOLUENE_TEMPERATURE_K, load_nist_benzene_toluene_vle_dataset, raoult_reference,
    };

    #[test]
    fn frozen_benzene_toluene_raoult_oracle_preserves_bubble_dew_round_trip() {
        let dataset = load_nist_benzene_toluene_vle_dataset().unwrap();
        for row in dataset.rows() {
            let reference = raoult_reference(&dataset, row).unwrap();
            let dew_pressure_pa = 1.0
                / (reference.vapor_benzene_mole_fraction
                    / dataset.rows().last().unwrap().experimental_pressure_pa
                    + reference.vapor_toluene_mole_fraction
                        / dataset.rows().first().unwrap().experimental_pressure_pa);
            let recovered_x_benzene = reference.vapor_benzene_mole_fraction * dew_pressure_pa
                / dataset.rows().last().unwrap().experimental_pressure_pa;
            assert!((dew_pressure_pa - reference.bubble_pressure_pa).abs() <= 1.0e-10);
            assert!(
                (recovered_x_benzene - reference.liquid_benzene_mole_fraction).abs() <= 1.0e-12
            );
            assert!(
                (reference.liquid_benzene_mole_fraction + reference.liquid_toluene_mole_fraction
                    - 1.0)
                    .abs()
                    <= 1.0e-12
            );
            assert!(
                (reference.vapor_benzene_mole_fraction + reference.vapor_toluene_mole_fraction
                    - 1.0)
                    .abs()
                    <= 1.0e-12
            );
            assert_eq!(row.temperature_k, BENZENE_TOLUENE_TEMPERATURE_K);
        }
    }
}
