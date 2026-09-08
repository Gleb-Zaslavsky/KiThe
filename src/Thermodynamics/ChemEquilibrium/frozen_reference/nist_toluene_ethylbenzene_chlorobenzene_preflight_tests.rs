//! Regression tests for the ternary VLE capability gate.

use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::ReactionExtentError;
use crate::Thermodynamics::ChemEquilibrium::frozen_reference_nist_toluene_ethylbenzene_chlorobenzene_preflight::{
    TernaryVlePreflightDecision, nist_ternary_thermoml_audit, preflight_primary_ternary_vle,
};

#[test]
fn primary_ternary_source_audit_preserves_actual_thermoml_semantics() {
    let source = nist_ternary_thermoml_audit();
    assert_eq!(source.doi, "10.1021/je020186c");
    assert_eq!(
        source.component_order,
        ["toluene", "ethylbenzene", "chlorobenzene"]
    );
    assert_eq!(source.property, "boiling temperature at pressure P, K");
    assert_eq!(source.pressure_levels_kpa, [26.66, 53.33, 79.99, 101.32]);
    assert_eq!(source.ternary_row_count, 48);
    assert!(source.liquid_composition_published);
    assert!(!source.vapor_composition_published);
    assert!(source.property_uncertainty_published);
    assert!(!source.same_payload_has_pure_saturation_support);
}

#[test]
fn primary_ternary_preflight_refuses_missing_exact_local_states_without_substitution() {
    let report = preflight_primary_ternary_vle();
    println!("{}", report.summary());
    assert_eq!(report.local_records.len(), 6);

    let toluene_gas = report
        .local_records
        .iter()
        .find(|record| record.request.record_key == "C7H8")
        .unwrap();
    let toluene_liquid = report
        .local_records
        .iter()
        .find(|record| record.request.record_key == "C7H8(L)")
        .unwrap();
    let ethylbenzene_gas = report
        .local_records
        .iter()
        .find(|record| record.request.record_key == "C8H10,ethylbenz")
        .unwrap();
    assert!(toluene_gas.is_available());
    assert!(toluene_liquid.is_available());
    assert!(ethylbenzene_gas.is_available());
    assert_eq!(toluene_gas.selected_library.as_deref(), Some("NASA_gas"));
    assert_eq!(toluene_gas.selected_record_key.as_deref(), Some("C7H8"));
    assert_eq!(
        toluene_liquid.selected_library.as_deref(),
        Some("NASA_cond")
    );
    assert_eq!(
        toluene_liquid.selected_record_key.as_deref(),
        Some("C7H8(L)")
    );
    assert_eq!(
        ethylbenzene_gas.selected_record_key.as_deref(),
        Some("C8H10,ethylbenz")
    );
    assert!(toluene_gas.supports_standard_gibbs && toluene_gas.supports_enthalpy);
    assert!(toluene_liquid.supports_standard_gibbs && toluene_liquid.supports_enthalpy);
    assert!(toluene_gas.temperature_interval_k.is_some());
    assert!(toluene_liquid.temperature_interval_k.is_some());
    assert!(toluene_gas.standard_state_pressure_provenance.is_some());
    assert!(ethylbenzene_gas.phase_model_compatible);

    for missing_key in ["C8H10(L),ethylbenz", "C6H5Cl", "C6H5Cl(L)"] {
        let missing = report
            .local_records
            .iter()
            .find(|record| record.request.record_key == missing_key)
            .unwrap();
        assert!(!missing.is_available());
        assert!(missing.temperature_interval_k.is_none());
        assert!(!missing.supports_standard_gibbs && !missing.supports_enthalpy);
    }

    let TernaryVlePreflightDecision::ValidationNotApplicable {
        missing_or_incompatible_records,
        ..
    } = &report.decision
    else {
        panic!("the current local repository must not claim this ternary fixture is executable");
    };
    assert!(
        missing_or_incompatible_records
            .iter()
            .any(|record| record.contains("ethylbenzene") && record.contains("C8H10(L),ethylbenz"))
    );
    assert!(
        missing_or_incompatible_records
            .iter()
            .any(|record| record.contains("chlorobenzene") && record.contains("C6H5Cl"))
    );

    let error = report.require_executable().unwrap_err();
    assert!(matches!(
        error,
        ReactionExtentError::ValidationNotApplicable { .. }
    ));
    assert!(report.summary().contains("C8H10(L),ethylbenz"));
    assert!(report.summary().contains("C6H5Cl(L)"));
}

#[test]
fn primary_ternary_preflight_uses_no_online_nist_fallback() {
    let report = preflight_primary_ternary_vle();
    for record in &report.local_records {
        if let Some(reason) = &record.unavailable_reason {
            assert!(reason.contains("NIST fallback disabled"));
        }
    }
}
