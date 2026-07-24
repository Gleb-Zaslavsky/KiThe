//! Regression tests for the NASA-7 Levenberg-Marquardt fitting strategies.
//!
//! These tests protect two distinct contracts:
//! - sequential LM may fit entropy and enthalpy in stages, but the published
//!   polynomial must remain accurate for all three thermodynamic properties;
//! - simultaneous LM must retain independent `Cp`, `dH`, and `dS` residual
//!   blocks while preserving the existing property-wise acceptance contract.
//!
//! A successful result is finite and keeps the maximum pointwise relative error
//! below two percent for every property over a real adjacent NASA interval.
//! The assertions rebuild raw property samples from the source ranges and call
//! `reports_for_coefficients` for the coefficient vector actually published by
//! the fitting method; coefficient finiteness alone is not considered success.

#[cfg(test)]
mod tests {
    use crate::Thermodynamics::DBhandlers::NASAdata::{Coeffs, NASAError, NASAdata};
    use crate::Thermodynamics::DBhandlers::NASAdata_fitting::FittingReport;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use std::collections::HashMap;

    const MAX_PROPERTY_RELATIVE_ERROR: f64 = 2.0e-2;

    fn load_adjacent_ranges(species: &str) -> (NASAdata, Coeffs, Coeffs) {
        let thermo_data = ThermoData::new();
        let record = thermo_data
            .LibThermoData
            .get("NASA_gas")
            .and_then(|library| library.get(species))
            .unwrap_or_else(|| panic!("missing NASA_gas fixture for {species}"));
        let mut nasa = NASAdata::new();
        nasa.from_serde(record.clone())
            .unwrap_or_else(|err| panic!("failed to parse {species}: {err}"));
        nasa.parse_coefficients()
            .unwrap_or_else(|err| panic!("failed to extract {species} coefficients: {err}"));

        let lower = nasa
            .interval_for_this_T(900.0)
            .unwrap_or_else(|err| panic!("missing lower range for {species}: {err}"))
            .1;
        let upper = nasa
            .interval_for_this_T(1100.0)
            .unwrap_or_else(|err| panic!("missing upper range for {species}: {err}"))
            .1;
        assert_ne!(
            lower.coeff, upper.coeff,
            "{species} fixture must cross a real boundary"
        );
        (nasa, lower, upper)
    }

    fn assert_property_reports(
        strategy: &str,
        species: &str,
        reports: &HashMap<String, FittingReport>,
    ) {
        for property in [
            "Cp fitting report",
            "dh fitting report",
            "ds fitting report",
        ] {
            let report = reports
                .get(property)
                .unwrap_or_else(|| panic!("{strategy} omitted {property} for {species}"));
            assert!(
                report.l2_norm.is_finite(),
                "{strategy} {species} {property} has non-finite L2"
            );
            assert!(
                report.max_norm.is_finite(),
                "{strategy} {species} {property} has non-finite max error"
            );
            assert!(
                report.max_norm <= MAX_PROPERTY_RELATIVE_ERROR,
                "{strategy} {species} {property} max relative error {} exceeds {}",
                report.max_norm,
                MAX_PROPERTY_RELATIVE_ERROR
            );
            assert_eq!(
                report.number_of_significant_points, 0,
                "{strategy} {species} {property} contains >10% residuals"
            );
        }
    }

    fn independently_compare_published_coefficients(
        nasa: &NASAdata,
        lower: Coeffs,
        upper: Coeffs,
    ) -> Result<HashMap<String, FittingReport>, NASAError> {
        let (temps, cp, h, s) =
            NASAdata::prepare_fit_data_for_2_range_fitting(lower, upper, 900.0, 1100.0)
                .expect("real NASA fixtures must produce raw adjacent-range samples");

        nasa.reports_for_coefficients(&temps, &cp, &h, &s, nasa.coeffs)
    }

    #[test]
    fn sequential_lm_keeps_all_properties_accurate_for_real_adjacent_ranges() {
        for species in ["CO", "CO2", "NO", "N2"] {
            let (mut nasa, lower, upper) = load_adjacent_ranges(species);
            let reports = nasa
                .fitting_adjacent2(lower.clone(), upper.clone(), 900.0, 1100.0)
                .unwrap_or_else(|err| panic!("sequential LM failed for {species}: {err}"));
            let independent_reports = independently_compare_published_coefficients(
                &nasa, lower, upper,
            )
            .unwrap_or_else(|err| panic!("sequential LM comparison failed for {species}: {err}"));

            assert_property_reports("sequential LM returned", species, &reports);
            assert_property_reports("sequential LM published", species, &independent_reports);
        }
    }

    #[test]
    fn simultaneous_lm_preserves_independent_property_accuracy() {
        for species in ["CO", "CO2", "NO", "N2"] {
            let (mut nasa, lower, upper) = load_adjacent_ranges(species);
            let reports = nasa
                .fitting_adjacent3(lower.clone(), upper.clone(), 900.0, 1100.0)
                .unwrap_or_else(|err| panic!("simultaneous LM failed for {species}: {err}"));
            let independent_reports = independently_compare_published_coefficients(
                &nasa, lower, upper,
            )
            .unwrap_or_else(|err| panic!("simultaneous LM comparison failed for {species}: {err}"));

            assert_property_reports("simultaneous LM returned", species, &reports);
            assert_property_reports("simultaneous LM published", species, &independent_reports);
        }
    }
}
