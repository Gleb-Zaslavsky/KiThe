//! Negative feasibility evidence for the ideal Ar/H2O/CO2 three-phase lifecycle.
//!
//! The offline catalog currently has no phase-qualified condensed `CO2`
//! thermochemistry. In addition, the existing local water-ice pair begins at
//! `200 K`, above the supplied CO2 sublimation correlation range. These tests
//! protect both facts from being hidden behind an accidental online fallback
//! or a guessed record name. They intentionally do not construct a reaction
//! basis, compute TPD, or invoke phase control.

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedThermochemistry;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemFactory, SubstanceSystemFactoryError, SubstanceSystemSpec,
        SubstanceSystemSpecBuilder, SubstancesContainer,
    };
    use crate::Thermodynamics::User_substances::Phases;
    use crate::Thermodynamics::User_substances_error::SubsDataError;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};

    const EXTERNAL_CO2_SUBLIMATION_LOWER_K: f64 = 154.26;
    const EXTERNAL_CO2_SUBLIMATION_UPPER_K: f64 = 195.89;
    const PREFERRED_LOWER_K: f64 = 180.0;
    const PREFERRED_UPPER_K: f64 = 195.0;

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("bundled offline thermochemistry repository must be available")
    }

    fn three_phase_spec() -> SubstanceSystemSpec {
        SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
            (
                "gas".to_string(),
                vec!["Ar".to_string(), "H2O".to_string(), "CO2".to_string()],
            ),
            ("ice".to_string(), vec!["H2O(s)".to_string()]),
            ("dry_ice".to_string(), vec!["CO2(s)".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            ("gas".to_string(), Phases::Gas),
            ("ice".to_string(), Phases::Solid),
            ("dry_ice".to_string(), Phases::Solid),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("Ar/H2O/CO2 three-phase preflight declaration must be valid")
    }

    fn water_ice_spec() -> SubstanceSystemSpec {
        SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["Ar".to_string(), "H2O".to_string()]),
            ("ice".to_string(), vec!["H2O(s)".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            ("gas".to_string(), Phases::Gas),
            ("ice".to_string(), Phases::Solid),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("Ar/H2O gas-ice preflight declaration must be valid")
    }

    fn co2_gas_spec() -> SubstanceSystemSpec {
        SubstanceSystemSpecBuilder::new(SubstancesContainer::SinglePhase(vec!["CO2".to_string()]))
            .with_phase_natures(Some(HashMap::from([("gas".to_string(), Phases::Gas)])))
            .with_library_priorities(vec!["NASA_gas".to_string()])
            .with_search_in_nist(false)
            .build()
            .expect("CO2 gas preflight declaration must be valid")
    }

    /// The existing gas record must be part of the feasibility gate as well:
    /// a solid closure cannot be paired with a gas standard state outside the
    /// latter's declared validity interval.
    #[test]
    fn local_co2_gas_interval_does_not_overlap_giauque_egan_pressure_range() {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            co2_gas_spec(),
            local_repository(),
        )
        .expect("CO2 gas record must resolve from the offline catalog");
        let bounds = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("local CO2 gas record must expose thermochemistry")
            .temperature_bounds();

        println!("CO2(g) and Giauque-Egan range compatibility preflight");
        println!(
            "  local CO2(g) interval: {:.2}..{:.2} K",
            bounds.lower(),
            bounds.upper()
        );
        println!(
            "  Giauque-Egan CO2(s) pressure interval: {EXTERNAL_CO2_SUBLIMATION_LOWER_K:.2}..{EXTERNAL_CO2_SUBLIMATION_UPPER_K:.2} K"
        );

        assert!(bounds.lower() > EXTERNAL_CO2_SUBLIMATION_UPPER_K);
    }

    /// The existing local water pair has a real interval, but it does not
    /// overlap the initial CO2 Antoine range. This is distinct from the
    /// missing dry-ice record and must remain visible in feasibility reports.
    #[test]
    fn ar_water_local_ice_interval_does_not_overlap_supplied_co2_oracle() {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            water_ice_spec(),
            local_repository(),
        )
        .expect("Ar/H2O gas-ice records must resolve from the offline catalog");
        let bounds = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("local water records must expose thermochemistry")
            .temperature_bounds();

        println!("Ar/H2O and CO2(s) external-range compatibility preflight");
        println!(
            "  local H2O(g)/H2O(s) interval: {:.2}..{:.2} K",
            bounds.lower(),
            bounds.upper()
        );
        println!(
            "  supplied CO2(s) external interval: {EXTERNAL_CO2_SUBLIMATION_LOWER_K:.2}..{EXTERNAL_CO2_SUBLIMATION_UPPER_K:.2} K"
        );

        assert!(bounds.lower() > EXTERNAL_CO2_SUBLIMATION_UPPER_K);
    }

    /// The requested dry-ice record is deliberately unavailable in the bundled
    /// offline catalog. A future fixture must introduce a reviewed local record
    /// with provenance and validity bounds, rather than silently using NIST.
    #[test]
    fn ar_water_co2_offline_catalog_lacks_condensed_co2_feasibility_preflight() {
        let error = SubstanceSystemFactory::resolve_phase_system_with_repository(
            three_phase_spec(),
            local_repository(),
        )
        .expect_err("the offline catalog must not pretend that CO2(s) is resolved");

        println!("Ar/H2O/CO2 ideal three-phase local preflight");
        println!(
            "  external CO2(s) sublimation domain: {EXTERNAL_CO2_SUBLIMATION_LOWER_K:.2}..{EXTERNAL_CO2_SUBLIMATION_UPPER_K:.2} K"
        );
        println!(
            "  requested local design window: {PREFERRED_LOWER_K:.2}..{PREFERRED_UPPER_K:.2} K"
        );
        println!("  result: {error}");

        assert!(matches!(
            error,
            SubstanceSystemFactoryError::Resolution(SubsDataError::SubstanceNotFound(name))
                if name == "CO2(s)"
        ));
    }
}
