//! Shared frozen fixtures and facade setup for equilibrium benchmarks.
//!
//! Keep benchmark-specific timing/report models in their own targets. This
//! module owns only the repeated physical fixtures so different dashboards
//! cannot silently drift to different species lists or inventories.

#![allow(dead_code)]

use KiThe::Thermodynamics::ChemEquilibrium::prelude::*;

pub const PRESSURE_PA: f64 = 101_325.0;
pub const REFERENCE_PRESSURE_PA: f64 = 100_000.0;
pub const REFERENCE_TEMPERATURE_K: f64 = 2_500.0;

pub const STANJAN_15: [&str; 15] = [
    "CH4", "O2", "CO2", "H2O", "N2", "N", "O", "NO", "OH", "H", "N2O", "CO", "H2", "NO2", "HO2",
];

pub const TP1907_GAS_17: [&str; 17] = [
    "Ar", "CH4", "CO", "CO2", "H", "HO2", "H2", "H2O", "N", "NH3", "NO", "NO2", "N2", "N2O", "O",
    "OH", "O2",
];

pub const RP1311_40: [&str; 40] = [
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

pub fn component(phase: Option<&str>, species: &str) -> PhaseComponentId {
    PhaseComponentId::new(PhaseId::new(phase.map(str::to_owned)), species)
}

pub fn tp1907_element_seed() -> MultiphaseInitialComposition {
    let phases = vec![
        PhaseSpec::ideal_gas(
            PhaseId::new(Some("gas".to_owned())),
            TP1907_GAS_17
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
        )
        .expect("TP-1907 gas seed phase must be valid"),
        PhaseSpec::pure_condensed(
            PhaseId::new(Some("graphite".to_owned())),
            vec!["C(gr)".to_owned()],
            PhysicalState::Solid,
        )
        .expect("TP-1907 graphite seed phase must be valid"),
    ];
    let layout =
        MultiphaseEquilibriumLayout::new(phases).expect("TP-1907 seed layout must be valid");
    MultiphaseInitialComposition::from_sparse(
        &layout,
        [
            (component(Some("gas"), "CH4"), 0.5),
            (component(Some("gas"), "CO2"), 0.50182736),
            (component(Some("gas"), "O2"), 0.7),
            (component(Some("gas"), "N2"), 4.4731044),
            (component(Some("gas"), "Ar"), 0.05364816),
        ]
        .to_vec(),
    )
    .expect("TP-1907 seed must be valid")
}
