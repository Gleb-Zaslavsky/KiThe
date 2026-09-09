//! Frozen NASA CEA RP-1311 Example 3 source data and preflight contract.
//!
//! This module intentionally stops at source validation for the first pass.
//! The production H,P route will be added only after the immutable external
//! rows and the independently reconstructed feed have been reviewed.

use serde::Deserialize;
use std::collections::BTreeMap;

use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
    FrozenReferenceDataset, FrozenReferenceError, FrozenReferenceRow, ReferenceColumnSchema,
};

pub(crate) const DATASET_ID: &str = "nasa_cea.rp1311.example3.hydrocarbon_air.hp.v1";
pub(crate) const TRACE_THRESHOLD: f64 = 1e-15;
pub(crate) const OXIDIZER_FUEL_MASS_RATIO: f64 = 17.0;
pub(crate) const FUEL_C7H8_MASS_FRACTION: f64 = 0.4;
pub(crate) const FUEL_C8H18_MASS_FRACTION: f64 = 0.6;
pub(crate) const SOURCE_SPECIFIC_ENTHALPY_J_KG: f64 = 317_838.0;

/// Test-side decomposition of NASA CEA's pseudo-reactant Air.
///
/// This is intentionally not exposed as a production `Air` alias. It exists
/// only to reconstruct the source enthalpy from ordinary local gas records.
pub(crate) const CEA_AIR_SPECIES: [(&str, f64); 4] = [
    ("N2", 0.780840),
    ("O2", 0.209476),
    ("Ar", 0.009365),
    ("CO2", 0.000319),
];
pub(crate) const CEA_SPECIES: [&str; 40] = [
    "Ar",
    "CN",
    "CO",
    "CO2",
    "COOH",
    "H",
    "H2",
    "H2O",
    "H2O2",
    "HCHO,formaldehy",
    "HCN",
    "HCO",
    "HCOOH",
    "HNC",
    "HNCO",
    "HNO",
    "HNO2",
    "HNO3",
    "HO2",
    "N",
    "N2",
    "N2H2",
    "N2O",
    "N2O3",
    "N2O4",
    "N3",
    "N3H",
    "NCO",
    "NH",
    "NH2",
    "NH2NO2",
    "NH2OH",
    "NH3",
    "NO",
    "NO2",
    "NO3",
    "O",
    "O2",
    "O3",
    "OH",
];

pub(crate) const TRACE_ONLY_SPECIES: [&str; 39] = [
    "(HCOOH)2",
    "C",
    "C10H8,naphthale",
    "C2",
    "C2H",
    "C2H2,acetylene",
    "C2H2,vinylidene",
    "C2H3,vinyl",
    "C2H4",
    "C2H4O,ethylen-o",
    "C2H5",
    "C2H5OH",
    "C2H6",
    "C3H3,1-propynl",
    "C3H3,2-propynl",
    "C3H6O,acetone",
    "C3H6O,propanal",
    "C3H6O,propylox",
    "C4H2,butadiyne",
    "C4H6,1butyne",
    "C4H6,2butyne",
    "C6H14,n-hexane",
    "C7H16,2-methylh",
    "CH",
    "CH2",
    "CH2CO,ketene",
    "CH2OH",
    "CH3",
    "CH3CHO,ethanal",
    "CH3CN",
    "CH3CO,acetyl",
    "CH3COOH",
    "CH3N2CH3",
    "CH3O",
    "CH3O2CH3",
    "CH3OCH3",
    "CH3OH",
    "CH3OOH",
    "CH4",
];

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CeaMoleFraction {
    pub(crate) species: String,
    pub(crate) mole_fraction: f64,
}

#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NasaCeaRp1311Example3Reference {
    pub(crate) problem_type: String,
    pub(crate) pressure_pa: f64,
    pub(crate) equilibrium_temperature_k: f64,
    pub(crate) source_specific_enthalpy_j_kg: f64,
    pub(crate) species_mole_fractions: Vec<CeaMoleFraction>,
}

/// Explicit source-side feed on a one-kilogram fuel basis.
///
/// Air is represented by the elemental composition published by CEA rather
/// than by a substituted O2/N2/Ar gas recipe. The resulting elemental totals
/// are the only quantities passed toward the later production formulation.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ReconstructedFeed {
    pub(crate) total_mass_kg: f64,
    pub(crate) fuel_mass_kg: f64,
    pub(crate) air_mass_kg: f64,
    /// Pure-liquid fuel amounts on the physical one-kilogram fuel basis.
    pub(crate) c7h8_moles: f64,
    pub(crate) c8h18_moles: f64,
    /// Nominal moles of the source-defined CEA Air pseudo-reactant.
    pub(crate) nominal_air_moles: f64,
    pub(crate) elemental_moles: BTreeMap<&'static str, f64>,
    pub(crate) target_enthalpy_j: f64,
}

pub(crate) fn reconstruct_feed() -> ReconstructedFeed {
    // SI values.  The source molecular masses are commonly tabulated in
    // g/mol; this source-side reconstruction must use kg/mol because its
    // masses are in kg and its amounts are in mol.
    const C7H8_MOLAR_MASS_KG_MOL: f64 = 0.092_141_04;
    const C8H18_MOLAR_MASS_KG_MOL: f64 = 0.114_232_0;
    reconstruct_feed_with_fuel_molar_masses(C7H8_MOLAR_MASS_KG_MOL, C8H18_MOLAR_MASS_KG_MOL)
}

/// Reconstructs the same physical CEA source basis using explicitly supplied
/// liquid-fuel molar masses in kg/mol.  Test fixtures use this variant with
/// masses resolved from local records, so their seed and source inventory are
/// tied to identical data provenance.
pub(crate) fn reconstruct_feed_with_fuel_molar_masses(
    c7h8_molar_mass_kg_mol: f64,
    c8h18_molar_mass_kg_mol: f64,
) -> ReconstructedFeed {
    const ATOMIC_MASS: [(&str, f64); 5] = [
        ("C", 12.011),
        ("H", 1.008),
        ("O", 15.999),
        ("N", 14.007),
        ("Ar", 39.948),
    ];
    const AIR_PER_NOMINAL_MOLE: [(&str, f64); 4] = [
        ("N", 1.561680),
        ("O", 0.419590),
        ("Ar", 0.009365),
        ("C", 0.000319),
    ];

    let fuel_mass_kg = 1.0;
    let air_mass_kg = OXIDIZER_FUEL_MASS_RATIO * fuel_mass_kg;
    let total_mass_kg = fuel_mass_kg + air_mass_kg;
    assert!(c7h8_molar_mass_kg_mol.is_finite() && c7h8_molar_mass_kg_mol > 0.0);
    assert!(c8h18_molar_mass_kg_mol.is_finite() && c8h18_molar_mass_kg_mol > 0.0);
    let c7h8_moles = FUEL_C7H8_MASS_FRACTION * fuel_mass_kg / c7h8_molar_mass_kg_mol;
    let c8h18_moles = FUEL_C8H18_MASS_FRACTION * fuel_mass_kg / c8h18_molar_mass_kg_mol;

    let air_molar_mass_kg_mol = AIR_PER_NOMINAL_MOLE
        .iter()
        .map(|(element, count)| {
            let atomic_mass = ATOMIC_MASS
                .iter()
                .find(|(name, _)| name == element)
                .map(|(_, value)| *value)
                .expect("CEA Air element must have an atomic mass");
            count * atomic_mass
        })
        .sum::<f64>()
        / 1000.0;
    let nominal_air_moles = air_mass_kg / air_molar_mass_kg_mol;

    let mut elemental_moles = BTreeMap::new();
    for (element, count) in AIR_PER_NOMINAL_MOLE {
        elemental_moles.insert(element, count * nominal_air_moles);
    }
    *elemental_moles.entry("C").or_default() += 7.0 * c7h8_moles + 8.0 * c8h18_moles;
    *elemental_moles.entry("H").or_default() += 8.0 * c7h8_moles + 18.0 * c8h18_moles;

    ReconstructedFeed {
        total_mass_kg,
        fuel_mass_kg,
        air_mass_kg,
        c7h8_moles,
        c8h18_moles,
        nominal_air_moles,
        elemental_moles,
        target_enthalpy_j: total_mass_kg * SOURCE_SPECIFIC_ENTHALPY_J_KG,
    }
}

const SCHEMA: [ReferenceColumnSchema; 5] = [
    ReferenceColumnSchema {
        column: "problem_type",
        unit: "text",
    },
    ReferenceColumnSchema {
        column: "pressure_pa",
        unit: "Pa",
    },
    ReferenceColumnSchema {
        column: "equilibrium_temperature_k",
        unit: "K",
    },
    ReferenceColumnSchema {
        column: "source_specific_enthalpy_j_kg",
        unit: "J/kg",
    },
    ReferenceColumnSchema {
        column: "species_mole_fractions",
        unit: "mole fraction",
    },
];

impl FrozenReferenceRow for NasaCeaRp1311Example3Reference {
    fn schema() -> &'static [ReferenceColumnSchema] {
        &SCHEMA
    }

    fn validate_rows(dataset_id: &str, rows: &[Self]) -> Result<(), FrozenReferenceError> {
        if rows.len() != 3 {
            return Err(FrozenReferenceError::InvalidRows {
                dataset_id: dataset_id.to_owned(),
                row: None,
                field: "rows".to_owned(),
                message: "RP-1311 Example 3 must contain exactly three pressure rows".to_owned(),
            });
        }
        for (index, row) in rows.iter().enumerate() {
            if row.problem_type != "HP" {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "problem_type".to_owned(),
                    message: format!("expected HP, got {}", row.problem_type),
                });
            }
            if !row.pressure_pa.is_finite()
                || row.pressure_pa <= 0.0
                || !row.equilibrium_temperature_k.is_finite()
                || !row.source_specific_enthalpy_j_kg.is_finite()
            {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "scalar".to_owned(),
                    message: "pressure, temperature and enthalpy must be finite".to_owned(),
                });
            }
            if row.species_mole_fractions.len() != CEA_SPECIES.len() {
                return Err(FrozenReferenceError::InvalidRows {
                    dataset_id: dataset_id.to_owned(),
                    row: Some(index),
                    field: "species_mole_fractions".to_owned(),
                    message: format!("expected {} published species", CEA_SPECIES.len()),
                });
            }
            for (position, (expected, value)) in CEA_SPECIES
                .iter()
                .zip(row.species_mole_fractions.iter())
                .enumerate()
            {
                if value.species != *expected
                    || !value.mole_fraction.is_finite()
                    || value.mole_fraction < 0.0
                    || value.mole_fraction > 1.0
                {
                    return Err(FrozenReferenceError::InvalidRows {
                        dataset_id: dataset_id.to_owned(),
                        row: Some(index),
                        field: format!("species_mole_fractions[{position}]"),
                        message: format!("expected non-negative finite {} mole fraction", expected),
                    });
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn load_dataset()
-> Result<FrozenReferenceDataset<NasaCeaRp1311Example3Reference>, String> {
    let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("src/Thermodynamics/ChemEquilibrium/frozen_reference/data/nasa_cea");
    FrozenReferenceDataset::load(
        directory.join("rp1311_example3_hp.metadata.json"),
        directory.join("rp1311_example3_hp.rows.json"),
    )
    .map_err(|error| error.to_string())
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_constraints::{
        EquilibriumConstraint, TemperatureBounds,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_multiphase_domain::{
        MultiphaseEquilibriumLayout, MultiphaseInitialComposition,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedThermochemistry;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::{
        ResolvedPhaseEnthalpyRequest, solve_resolved_ph,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::TraceSpeciesSeedPolicy;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::EquilibriumSolveOptions;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemFactory, SubstanceSystemSpecBuilder, SubstancesContainer,
        element_composition_and_molar_mass,
    };
    use crate::Thermodynamics::User_substances::Phases;
    use crate::Thermodynamics::phase_layout::{PhaseComponentId, PhaseId};
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};

    use super::*;

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("bundled offline thermochemistry repository must be available")
    }

    /// Resolves the actual local liquid records used to audit the source fuel
    /// basis.  The conversion is deliberately visible at this boundary:
    /// `element_composition_and_molar_mass` returns g/mol, while the feed
    /// reconstruction uses kg and mol.
    fn local_liquid_fuel_molar_masses_kg_mol(repository: Arc<ThermoRepository>) -> (f64, f64) {
        let spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
                ("fuel_c7h8".to_owned(), vec!["C7H8(L)".to_owned()]),
                ("fuel_c8h18".to_owned(), vec!["C8H18(L),n-octa".to_owned()]),
            ])))
            .with_phase_natures(Some(HashMap::from([
                ("fuel_c7h8".to_owned(), Phases::Liquid),
                ("fuel_c8h18".to_owned(), Phases::Liquid),
            ])))
            .with_library_priorities(vec!["NASA_cond".to_owned()])
            .with_search_in_nist(false)
            .build()
            .expect("fuel-mass audit spec must be valid");
        let resolved =
            SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)
                .expect("local liquid fuel records must resolve offline");
        assert!(!resolved.report().nist_fallback_enabled());
        let (_, masses_g_mol, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .expect("local liquid fuel molar masses must be available");
        (
            masses_g_mol["C7H8(L)"] / 1000.0,
            masses_g_mol["C8H18(L),n-octa"] / 1000.0,
        )
    }

    #[test]
    fn rp1311_example3_source_dataset_has_complete_semantic_universe() {
        let dataset = load_dataset().expect("RP-1311 frozen source must load");
        assert_eq!(dataset.metadata().dataset_id, DATASET_ID);
        assert_eq!(dataset.rows().len(), 3);
        assert_eq!(CEA_SPECIES.len(), 40);
        assert_eq!(TRACE_ONLY_SPECIES.len(), 39);
        assert_eq!(TRACE_THRESHOLD, 1e-15);
        assert!(dataset.rows().iter().all(|row| row.problem_type == "HP"));
        assert_eq!(dataset.rows()[0].pressure_pa, 100.0e5);
        assert_eq!(dataset.rows()[1].pressure_pa, 10.0e5);
        assert_eq!(dataset.rows()[2].pressure_pa, 1.0e5);
        assert!(dataset.rows().iter().all(|row| {
            row.species_mole_fractions
                .iter()
                .all(|value| value.mole_fraction >= 0.0)
        }));
        println!(
            "RP-1311 Example 3 preflight: rows=3 non_trace_species=40 trace_only_species=39 verdict=ExternalFixtureSourceComplete"
        );
    }

    #[test]
    fn rp1311_example3_feed_reconstruction_preserves_source_conventions() {
        let feed = reconstruct_feed();
        assert_eq!(feed.fuel_mass_kg, 1.0);
        assert_eq!(
            feed.air_mass_kg / feed.fuel_mass_kg,
            OXIDIZER_FUEL_MASS_RATIO
        );
        assert!((FUEL_C7H8_MASS_FRACTION + FUEL_C8H18_MASS_FRACTION - 1.0).abs() < 1e-15);
        assert_eq!(feed.total_mass_kg, 18.0);
        assert_eq!(feed.target_enthalpy_j, 18.0 * SOURCE_SPECIFIC_ENTHALPY_J_KG);
        assert!((feed.c7h8_moles - 0.4 / 0.092_141_04).abs() < 1.0e-12);
        assert!((feed.c8h18_moles - 0.6 / 0.114_232_0).abs() < 1.0e-12);
        assert!(feed.nominal_air_moles > 500.0);
        for element in ["C", "H", "O", "N", "Ar"] {
            assert!(
                feed.elemental_moles
                    .get(element)
                    .copied()
                    .unwrap_or_default()
                    > 0.0
            );
        }
        println!(
            "RP-1311 feed preflight: n_C7H8={:.9e} mol n_C8H18={:.9e} mol nominal_Air={:.9e} mol total_mass={:.6} kg O/F={:.6} H_target={:.6e} J elements={:?}",
            feed.c7h8_moles,
            feed.c8h18_moles,
            feed.nominal_air_moles,
            feed.total_mass_kg,
            feed.air_mass_kg / feed.fuel_mass_kg,
            feed.target_enthalpy_j,
            feed.elemental_moles
        );
    }

    #[test]
    fn rp1311_example3_local_liquid_mass_audit_preserves_one_kg_fuel_basis() {
        let (c7_mass_kg_mol, c8_mass_kg_mol) =
            local_liquid_fuel_molar_masses_kg_mol(local_repository());
        let feed = reconstruct_feed_with_fuel_molar_masses(c7_mass_kg_mol, c8_mass_kg_mol);
        let reconstructed_fuel_mass =
            feed.c7h8_moles * c7_mass_kg_mol + feed.c8h18_moles * c8_mass_kg_mol;

        // This follows resolved record masses rather than the source-side
        // constants used by `reconstruct_feed`, preventing a self-consistent
        // g/mol-versus-kg/mol mistake from passing unnoticed.
        assert!((reconstructed_fuel_mass - 1.0).abs() < 1.0e-12);
        assert!((feed.air_mass_kg - 17.0).abs() < 1.0e-12);
        assert!((reconstructed_fuel_mass + feed.air_mass_kg - 18.0).abs() < 1.0e-12);
        assert!(feed.c7h8_moles > 4.0 && feed.c8h18_moles > 5.0);
        println!(
            "RP-1311 local mass audit: n_C7H8={:.9e} mol n_C8H18={:.9e} mol fuel_mass={:.9e} kg air_mass={:.9e} kg total_mass={:.9e} kg",
            feed.c7h8_moles,
            feed.c8h18_moles,
            reconstructed_fuel_mass,
            feed.air_mass_kg,
            reconstructed_fuel_mass + feed.air_mass_kg,
        );
    }

    #[test]
    fn rp1311_example3_liquid_fuel_enthalpy_preflight_is_offline_and_explicit() {
        let spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
                ("liquid_c7h8".to_owned(), vec!["C7H8(L)".to_owned()]),
                (
                    "liquid_c8h18".to_owned(),
                    vec!["C8H18(L),n-octa".to_owned()],
                ),
            ])))
            .with_phase_natures(Some(HashMap::from([
                ("liquid_c7h8".to_owned(), Phases::Liquid),
                ("liquid_c8h18".to_owned(), Phases::Liquid),
            ])))
            .with_library_priorities(vec!["NASA_cond".to_owned()])
            .with_search_in_nist(false)
            .build()
            .expect("liquid fuel preflight spec must be valid");
        let resolved =
            SubstanceSystemFactory::resolve_phase_system_with_repository(spec, local_repository())
                .expect("both liquid fuel records must resolve offline");
        assert!(!resolved.report().nist_fallback_enabled());

        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("liquid fuel records must provide enthalpy closures");
        let bounds = thermochemistry.temperature_bounds();
        assert!(bounds.contains(298.15));
        let enthalpies = thermochemistry
            .evaluate_enthalpy(298.15)
            .expect("local liquid fuel enthalpies must be evaluable at 298.15 K");
        let (_, molar_masses, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .expect("liquid fuel elemental composition must be available");
        let c7_moles = 0.4 / (molar_masses["C7H8(L)"] / 1000.0);
        let c8_moles = 0.6 / (molar_masses["C8H18(L),n-octa"] / 1000.0);
        let fuel_enthalpy_j = c7_moles * enthalpies[0] + c8_moles * enthalpies[1];
        assert!(fuel_enthalpy_j.is_finite());
        println!(
            "RP-1311 liquid-fuel enthalpy preflight: bounds={:.3}..{:.3} K h_C7H8={:.6e} J/mol h_C8H18={:.6e} J/mol fuel_H={:.6e} J",
            bounds.lower(),
            bounds.upper(),
            enthalpies[0],
            enthalpies[1],
            fuel_enthalpy_j
        );
        println!(
            "CEA Air enthalpy comparison: deferred because the source defines Air by elemental composition, not by a single local thermochemical record"
        );
    }

    #[test]
    fn rp1311_example3_cea_air_reactant_enthalpy_characterization() {
        let air_spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([(
                "air_source_reconstruction".to_owned(),
                CEA_AIR_SPECIES
                    .iter()
                    .map(|(name, _)| (*name).to_owned())
                    .collect(),
            )])))
            .with_phase_natures(Some(HashMap::from([(
                "air_source_reconstruction".to_owned(),
                Phases::Gas,
            )])))
            .with_library_priorities(vec!["NASA_gas".to_owned()])
            .with_search_in_nist(false)
            .build()
            .expect("CEA Air source reconstruction spec must be valid");
        let resolved_air = SubstanceSystemFactory::resolve_phase_system_with_repository(
            air_spec,
            local_repository(),
        )
        .expect("N2/O2/Ar/CO2 must resolve from local NASA gas data");
        assert!(!resolved_air.report().nist_fallback_enabled());
        assert_eq!(resolved_air.phase_specs().len(), 1);
        assert_eq!(resolved_air.phase_specs()[0].components().len(), 4);

        let air_thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved_air)
            .expect("CEA Air reconstruction species must expose thermochemistry");
        assert!(air_thermochemistry.temperature_bounds().contains(700.0));
        let air_enthalpies = air_thermochemistry
            .evaluate_enthalpy(700.0)
            .expect("local NASA gas enthalpies must be evaluable at 700 K");
        let (_, air_molar_masses, _) = element_composition_and_molar_mass(
            resolved_air.phase_data(),
            resolved_air.layout(),
            None,
        )
        .expect("CEA Air reconstruction composition must be available");

        let air_fraction_sum: f64 = CEA_AIR_SPECIES.iter().map(|(_, fraction)| *fraction).sum();
        assert!((2.0 * CEA_AIR_SPECIES[0].1 - 1.561680).abs() < 1e-12);
        assert!((2.0 * CEA_AIR_SPECIES[1].1 + 2.0 * CEA_AIR_SPECIES[3].1 - 0.419590).abs() < 1e-12);
        assert!((CEA_AIR_SPECIES[2].1 - 0.009365).abs() < 1e-15);
        assert!((CEA_AIR_SPECIES[3].1 - 0.000319).abs() < 1e-15);
        assert!((air_fraction_sum - 1.0).abs() < 1e-12);

        let air_molar_mass_kg_mol: f64 = CEA_AIR_SPECIES
            .iter()
            .map(|(species, fraction)| fraction * air_molar_masses[*species] / 1000.0)
            .sum();
        let nominal_air_moles = 17.0 / air_molar_mass_kg_mol;
        let h_air_molar: f64 = CEA_AIR_SPECIES
            .iter()
            .zip(air_enthalpies.iter())
            .map(|((_, fraction), enthalpy)| fraction * enthalpy)
            .sum();
        let h_air = nominal_air_moles * h_air_molar;

        let liquid_spec =
            SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
                ("fuel_c7h8".to_owned(), vec!["C7H8(L)".to_owned()]),
                ("fuel_c8h18".to_owned(), vec!["C8H18(L),n-octa".to_owned()]),
            ])))
            .with_phase_natures(Some(HashMap::from([
                ("fuel_c7h8".to_owned(), Phases::Liquid),
                ("fuel_c8h18".to_owned(), Phases::Liquid),
            ])))
            .with_library_priorities(vec!["NASA_cond".to_owned()])
            .with_search_in_nist(false)
            .build()
            .expect("fuel enthalpy spec must be valid");
        let resolved_fuel = SubstanceSystemFactory::resolve_phase_system_with_repository(
            liquid_spec,
            local_repository(),
        )
        .expect("fuel liquid records must resolve without NIST fallback");
        assert!(!resolved_fuel.report().nist_fallback_enabled());
        let fuel_thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved_fuel)
            .expect("fuel liquid records must expose thermochemistry");
        assert!(fuel_thermochemistry.temperature_bounds().contains(298.15));
        let fuel_enthalpies = fuel_thermochemistry
            .evaluate_enthalpy(298.15)
            .expect("fuel liquid enthalpies must be evaluable");
        let (_, fuel_molar_masses, _) = element_composition_and_molar_mass(
            resolved_fuel.phase_data(),
            resolved_fuel.layout(),
            None,
        )
        .expect("fuel composition must be available");
        let h_fuel = 0.4 / (fuel_molar_masses["C7H8(L)"] / 1000.0) * fuel_enthalpies[0]
            + 0.6 / (fuel_molar_masses["C8H18(L),n-octa"] / 1000.0) * fuel_enthalpies[1];
        let h_local = h_air + h_fuel;
        let h_cea = 18.0 * SOURCE_SPECIFIC_ENTHALPY_J_KG;
        let specific_h_local = h_local / 18.0;
        let relative_delta = (specific_h_local - SOURCE_SPECIFIC_ENTHALPY_J_KG).abs()
            / SOURCE_SPECIFIC_ENTHALPY_J_KG;
        let verdict = if relative_delta < 0.01 {
            "ReactantEnthalpyConventionAligned"
        } else {
            "ReactantEnthalpyConventionOffset"
        };
        println!("RP-1311 CEA Air enthalpy characterization");
        println!(
            "  Air fractions: N2={:.6} O2={:.6} Ar={:.6} CO2={:.6}",
            CEA_AIR_SPECIES[0].1, CEA_AIR_SPECIES[1].1, CEA_AIR_SPECIES[2].1, CEA_AIR_SPECIES[3].1
        );
        println!(
            "  Air molar mass={:.9e} kg/mol nominal_air_moles={:.9e}",
            air_molar_mass_kg_mol, nominal_air_moles
        );
        println!(
            "  h_N2={:.9e} h_O2={:.9e} h_Ar={:.9e} h_CO2={:.9e} J/mol",
            air_enthalpies[0], air_enthalpies[1], air_enthalpies[2], air_enthalpies[3]
        );
        println!(
            "  h_air_molar={:.9e} J/mol H_air={:.9e} J H_fuel={:.9e} J",
            h_air_molar, h_air, h_fuel
        );
        println!(
            "  H_local={:.9e} J H_CEA={:.9e} J h_local={:.9e} J/kg h_CEA={:.9e} J/kg delta={:.9e} J/kg relative={:.9e} verdict={}",
            h_local,
            h_cea,
            specific_h_local,
            SOURCE_SPECIFIC_ENTHALPY_J_KG,
            specific_h_local - SOURCE_SPECIFIC_ENTHALPY_J_KG,
            relative_delta,
            verdict
        );
        assert!(h_local.is_finite() && specific_h_local.is_finite());
    }

    /// Runs the canonical production P,H route against all three frozen CEA
    /// pressure rows.  The liquid fuel records are used to reconstruct and
    /// audit the source enthalpy; the high-temperature equilibrium universe
    /// is deliberately gas-only, matching the published CEA species domain.
    ///
    /// This is an external characterization test, not a thermochemistry
    /// fitting test.  Its assertions guard identity, temperature, species
    /// comparison, conservation, and offline lookup provenance separately.
    #[test]
    #[ignore = "release characterization of the RP-1311 Example 3 production P,H route"]
    fn i5_nasa_cea_rp1311_example3_production_ph_comparison() {
        // Recorded release evidence uses the local NASA gas catalog against a
        // separately frozen CEA table. These are deliberately broad regression
        // envelopes: they detect a broken feed/reference convention without
        // pretending two database revisions are bitwise-identical.
        const MAX_TEMPERATURE_DELTA_K: f64 = 2.0;
        const MAX_MAJOR_SPECIES_RELATIVE_ERROR: f64 = 0.12;
        const MAX_MINOR_TRACE_LOG10_ERROR: f64 = 1.1;
        const MAX_ALL_SPECIES_RELATIVE_ERROR: f64 = 1.0;
        const STANDARD_STATE_PRESSURE_PA: f64 = 100_000.0;
        let dataset = load_dataset().expect("RP-1311 source dataset must load");
        let repository = local_repository();
        let (c7h8_molar_mass_kg_mol, c8h18_molar_mass_kg_mol) =
            local_liquid_fuel_molar_masses_kg_mol(repository.clone());
        let feed = reconstruct_feed_with_fuel_molar_masses(
            c7h8_molar_mass_kg_mol,
            c8h18_molar_mass_kg_mol,
        );
        let gas_species: Vec<String> = CEA_SPECIES.iter().map(|name| (*name).to_owned()).collect();
        // The local NASA gas catalog has no octane gas record.  The initial
        // seed therefore uses an element-equivalent CO/H2 representation;
        // the independently reconstructed CEA source enthalpy remains the
        // P,H target and is not replaced by this computational seed.
        let spec = SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from(
            [("gas".to_owned(), gas_species)],
        )))
        .with_phase_natures(Some(HashMap::from([("gas".to_owned(), Phases::Gas)])))
        .with_library_priorities(vec!["NASA_gas".to_owned()])
        .with_search_in_nist(false)
        .build()
        .expect("RP-1311 gas production spec must be valid");
        let resolved =
            SubstanceSystemFactory::resolve_phase_system_with_repository(spec, repository)
                .expect("all RP-1311 gas and fuel records must resolve offline");
        assert!(!resolved.report().nist_fallback_enabled());
        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("RP-1311 gas universe must expose thermochemistry");
        let bounds = thermochemistry.temperature_bounds();
        assert!(bounds.contains(700.0));
        assert!(bounds.contains(2400.0));
        let (_, molar_masses, _) =
            element_composition_and_molar_mass(resolved.phase_data(), resolved.layout(), None)
                .expect("RP-1311 gas composition must be available");
        let air_molar_mass = CEA_AIR_SPECIES
            .iter()
            .map(|(name, fraction)| fraction * molar_masses[*name] / 1000.0)
            .sum::<f64>();
        let nominal_air_moles = feed.air_mass_kg / air_molar_mass;
        let c7_moles = feed.c7h8_moles;
        let c8_moles = feed.c8h18_moles;
        let reconstructed_fuel_mass =
            c7_moles * c7h8_molar_mass_kg_mol + c8_moles * c8h18_molar_mass_kg_mol;
        let reconstructed_air_mass = nominal_air_moles * air_molar_mass;
        assert!((reconstructed_fuel_mass - feed.fuel_mass_kg).abs() < 1.0e-12);
        assert!((reconstructed_air_mass - feed.air_mass_kg).abs() < 1.0e-10);
        assert!(
            (reconstructed_fuel_mass + reconstructed_air_mass - feed.total_mass_kg).abs() < 1.0e-10
        );
        let fuel_carbon_moles = 7.0 * c7_moles + 8.0 * c8_moles;
        let fuel_hydrogen_moles = 8.0 * c7_moles + 18.0 * c8_moles;
        let carbon_monoxide_moles = fuel_carbon_moles;
        let hydrogen_moles = fuel_hydrogen_moles / 2.0;
        let oxygen_adjustment_moles = carbon_monoxide_moles / 2.0;
        assert!(hydrogen_moles > 0.0 && oxygen_adjustment_moles > 0.0);
        let layout = MultiphaseEquilibriumLayout::new(resolved.phase_specs().to_vec())
            .expect("production layout must be constructible");
        let entries = CEA_AIR_SPECIES
            .iter()
            .map(|(name, fraction)| {
                (
                    PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), *name),
                    nominal_air_moles * fraction
                        - if *name == "O2" {
                            oxygen_adjustment_moles
                        } else {
                            0.0
                        },
                )
            })
            .chain([
                (
                    PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "CO"),
                    carbon_monoxide_moles,
                ),
                (
                    PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), "H2"),
                    hydrogen_moles,
                ),
            ])
            .collect();
        let initial = MultiphaseInitialComposition::from_sparse(&layout, entries)
            .expect("source-side gas inventory must fit the resolved layout");

        // Independent element check for the computational CO/H2 seed.  This
        // guards the seed transformation without turning it into a production
        // pseudo-species or an implicit Air model.
        let air_n2 = nominal_air_moles * CEA_AIR_SPECIES[0].1;
        let air_o2 = nominal_air_moles * CEA_AIR_SPECIES[1].1;
        let air_ar = nominal_air_moles * CEA_AIR_SPECIES[2].1;
        let air_co2 = nominal_air_moles * CEA_AIR_SPECIES[3].1;
        let seed_elements = [
            ("C", carbon_monoxide_moles + air_co2),
            ("H", 2.0 * hydrogen_moles),
            (
                "O",
                2.0 * (air_o2 - oxygen_adjustment_moles) + 2.0 * air_co2 + carbon_monoxide_moles,
            ),
            ("N", 2.0 * air_n2),
            ("Ar", air_ar),
        ];
        for (element, seed_total) in seed_elements {
            let source_total = *feed
                .elemental_moles
                .get(element)
                .expect("reconstructed source feed must contain every seed element");
            assert!(
                (seed_total - source_total).abs() < 1.0e-8 * source_total.abs().max(1.0),
                "CO/H2 seed changed {element} total: seed={seed_total:e} source={source_total:e}"
            );
        }
        println!(
            "feed audit: n_C7H8={c7_moles:.9e} mol n_C8H18={c8_moles:.9e} mol fuel_C={fuel_carbon_moles:.9e} mol fuel_H={fuel_hydrogen_moles:.9e} mol nominal_Air={nominal_air_moles:.9e} mol reconstructed_mass={:.9e} kg target_H={:.9e} J",
            reconstructed_fuel_mass + reconstructed_air_mass,
            feed.target_enthalpy_j,
        );
        println!("seed element totals: {:?}", seed_elements);

        println!("RP-1311 Example 3 production P,H comparison");
        println!(
            "pressure_Pa | P/P0 | CEA_T_K | KiThe_T_K | delta_T_K | max_species_rel | residual | balance | backend/iterations | provenance"
        );
        for reference in dataset.rows() {
            let pressure_ratio = reference.pressure_pa / STANDARD_STATE_PRESSURE_PA;
            let expected_ratio = reference.pressure_pa / 100_000.0;
            assert!((pressure_ratio - expected_ratio).abs() < 1.0e-12);
            let request = ResolvedPhaseEnthalpyRequest::from_resolved_thermochemistry(
                &resolved,
                initial.clone(),
                EquilibriumConstraint::ph(
                    reference.pressure_pa,
                    STANDARD_STATE_PRESSURE_PA,
                    feed.target_enthalpy_j,
                    700.0,
                )
                .expect("CEA P,H constraint must be valid"),
                TemperatureBounds::new(bounds.lower(), bounds.upper()).unwrap(),
                thermochemistry.clone(),
            )
            .unwrap()
            .with_solve_options(
                EquilibriumSolveOptions::default()
                    .with_trace_seed_policy(TraceSpeciesSeedPolicy::Absolute { floor: 1.0e-30 }),
            );
            let solution = solve_resolved_ph(request)
                .expect("canonical production P,H route must solve RP-1311 row");
            let temperature = solution.temperature();
            let max_species_rel = CEA_SPECIES
                .iter()
                .zip(reference.species_mole_fractions.iter())
                .map(|(name, external)| {
                    let component =
                        PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), *name);
                    let actual = solution
                        .equilibrium()
                        .mole_fraction_for(&component)
                        .expect("published solution must contain every CEA gas species");
                    (actual - external.mole_fraction).abs() / external.mole_fraction.max(1.0e-12)
                })
                .fold(0.0_f64, f64::max);
            let max_major_rel = CEA_SPECIES
                .iter()
                .zip(reference.species_mole_fractions.iter())
                .filter(|(_, external)| external.mole_fraction >= 1.0e-3)
                .map(|(name, external)| {
                    let component =
                        PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), *name);
                    let actual = solution
                        .equilibrium()
                        .mole_fraction_for(&component)
                        .expect("published solution must contain every major gas species");
                    (actual - external.mole_fraction).abs() / external.mole_fraction
                })
                .fold(0.0_f64, f64::max);
            let max_minor_trace_log10 = CEA_SPECIES
                .iter()
                .zip(reference.species_mole_fractions.iter())
                .filter(|(_, external)| {
                    external.mole_fraction > 0.0 && external.mole_fraction < 1.0e-3
                })
                .map(|(name, external)| {
                    let component =
                        PhaseComponentId::new(PhaseId::new(Some("gas".to_owned())), *name);
                    let actual = solution
                        .equilibrium()
                        .mole_fraction_for(&component)
                        .expect("published solution must contain every minor gas species");
                    actual.max(1.0e-300).log10() - external.mole_fraction.log10()
                })
                .map(f64::abs)
                .fold(0.0_f64, f64::max);
            let fraction_sum: f64 = CEA_SPECIES
                .iter()
                .map(|name| {
                    solution
                        .equilibrium()
                        .mole_fraction_for(&PhaseComponentId::new(
                            PhaseId::new(Some("gas".to_owned())),
                            *name,
                        ))
                        .unwrap_or_default()
                })
                .sum();
            let validation = solution.equilibrium().accepted_solution().validation();
            let delta_t = temperature - reference.equilibrium_temperature_k;
            println!(
                "{:.0} | {:.0} | {:.3} | {:.3} | {:+.3} | all={:.3e} major={:.3e} minor_trace_log10={:.3e} | {:.3e} | {:.3e} | {:?}/{} | NASA_gas/local/offline",
                reference.pressure_pa,
                pressure_ratio,
                reference.equilibrium_temperature_k,
                temperature,
                delta_t,
                max_species_rel,
                max_major_rel,
                max_minor_trace_log10,
                validation.residual_l2_norm,
                validation.max_abs_element_balance_error,
                solution.equilibrium().solve_report().accepted_backend,
                solution.report().inner_nonlinear_iterations(),
            );
            assert!(temperature.is_finite());
            assert!((200.0..=6000.0).contains(&temperature));
            assert!(delta_t.is_finite() && max_major_rel.is_finite());
            assert!(
                delta_t.abs() <= MAX_TEMPERATURE_DELTA_K,
                "CEA temperature regression at P={:.0} Pa: delta_T={delta_t:.6} K exceeds {MAX_TEMPERATURE_DELTA_K} K",
                reference.pressure_pa,
            );
            assert!(
                max_major_rel <= MAX_MAJOR_SPECIES_RELATIVE_ERROR,
                "CEA major-species regression at P={:.0} Pa: max_relative={max_major_rel:.6e}",
                reference.pressure_pa,
            );
            assert!(
                max_minor_trace_log10 <= MAX_MINOR_TRACE_LOG10_ERROR,
                "CEA minor/trace regression at P={:.0} Pa: max_log10={max_minor_trace_log10:.6e}",
                reference.pressure_pa,
            );
            assert!(
                max_species_rel <= MAX_ALL_SPECIES_RELATIVE_ERROR,
                "CEA all-species regression at P={:.0} Pa: max_relative={max_species_rel:.6e}",
                reference.pressure_pa,
            );
            assert!((fraction_sum - 1.0).abs() < 1.0e-8);
            assert!(
                validation.residual_l2_norm.is_finite() && validation.residual_l2_norm < 1.0e-5
            );
            assert!(
                validation.max_abs_element_balance_error.is_finite()
                    && validation.max_abs_element_balance_error < 1.0e-5
            );
        }
    }
}
