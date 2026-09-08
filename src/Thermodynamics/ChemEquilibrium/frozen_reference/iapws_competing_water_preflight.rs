//! Feasibility evidence for the competing-water-candidate phase-order study.
//!
//! This module deliberately stops before phase-control. It answers the
//! prerequisite question first: can the bundled local catalog resolve Ar(g),
//! H2O(g), H2O(l), and ice Ih together, and what temperature domain do those
//! exact records expose? The later order-invariance story must not silently
//! extrapolate a record or replace this lookup with synthetic data.

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_ph_workflow::ResolvedThermochemistry;
    use crate::Thermodynamics::User_PhaseOrSolution::{
        SubstanceSystemFactory, SubstanceSystemSpec, SubstanceSystemSpecBuilder,
        SubstancesContainer,
    };
    use crate::Thermodynamics::User_substances::Phases;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("bundled offline thermochemistry repository must be available")
    }

    fn competing_water_spec() -> SubstanceSystemSpec {
        SubstanceSystemSpecBuilder::new(SubstancesContainer::MultiPhase(HashMap::from([
            ("gas".to_string(), vec!["Ar".to_string(), "H2O".to_string()]),
            ("liquid".to_string(), vec!["H2O".to_string()]),
            ("ice".to_string(), vec!["H2O(s)".to_string()]),
        ])))
        .with_phase_natures(Some(HashMap::from([
            ("gas".to_string(), Phases::Gas),
            ("liquid".to_string(), Phases::Liquid),
            ("ice".to_string(), Phases::Solid),
        ])))
        .with_library_priorities(vec!["NASA_gas".to_string(), "NASA_cond".to_string()])
        .with_search_in_nist(false)
        .build()
        .expect("competing water preflight declaration must be structurally valid")
    }

    #[test]
    fn competing_water_local_records_resolve_without_nist() {
        let resolved = SubstanceSystemFactory::resolve_phase_system_with_repository(
            competing_water_spec(),
            local_repository(),
        )
        .expect("Ar/H2O gas-liquid-ice records must resolve from the local catalog");

        assert!(!resolved.report().nist_fallback_enabled());
        assert_eq!(resolved.phase_specs().len(), 3);
        let phase = |name: &str| {
            resolved
                .phase_specs()
                .iter()
                .find(|spec| spec.id().as_option().as_deref() == Some(name))
                .unwrap_or_else(|| panic!("resolved preflight lacks phase '{name}'"))
        };
        assert_eq!(phase("gas").components(), ["Ar", "H2O"]);
        assert_eq!(phase("liquid").components(), ["H2O"]);
        assert_eq!(phase("ice").components(), ["H2O(s)"]);

        let thermochemistry = ResolvedThermochemistry::from_resolved_system(&resolved)
            .expect("all competing water records must expose common thermochemistry");
        let bounds = thermochemistry.temperature_bounds();

        println!("Ar + H2O competing-phase local preflight");
        println!("  phases: {:?}", resolved.phase_specs());
        println!(
            "  common local temperature bounds: {:.6}..{:.6} K",
            bounds.lower(),
            bounds.upper()
        );
        for provenance in thermochemistry.provenance() {
            println!(
                "  component={} phase={:?} library={} record={} state={}",
                provenance.component().substance,
                provenance.component().phase,
                provenance.library(),
                provenance.record_key(),
                provenance.state()
            );
            println!(
                "    standard-state pressure: {:?}",
                provenance.standard_state_pressure()
            );
        }

        // This preflight is intentionally only a catalog/contract check. The
        // triple-point working temperature and competing TPD assertions are
        // selected only after this evidence has been reviewed.
        assert!(bounds.lower().is_finite() && bounds.upper().is_finite());
        // The current local NASA gas/condensed records meet only at the
        // declared 273.15 K endpoint. This is useful feasibility evidence,
        // but it is not yet a robust below-triple-point interval.
        assert!((bounds.lower() - 273.15).abs() < 1e-9);
        assert!((bounds.upper() - 273.15).abs() < 1e-9);
    }
}
