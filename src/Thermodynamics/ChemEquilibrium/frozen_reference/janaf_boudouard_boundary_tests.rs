//! Frozen-JANAF Boudouard `P,T` boundary and lifecycle evidence.
//!
//! Unlike `janaf_boudouard_tests`, this module owns pressure roots and bounded
//! phase-control stories. It reuses the same primary frozen species data but
//! keeps the external analytical oracle, independent I1/I2 root, canonical
//! TPD root, and production lifecycle as explicitly separate routes.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        ReactionExtentError, SolveError,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::PhaseStatus;
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, JanafBoudouardReference,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_boundary::{
        BOUDOUARD_BOUNDARY_Y_CO, BOUDOUARD_BOUNDARY_Y_CO2, JanafBoudouardBoundaryReference,
        JanafBoudouardBoundaryReport, JanafBoudouardBoundaryRow, JanafBoudouardPressureRoot,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_janaf_boudouard_thermochemistry::JANAF_STANDARD_PRESSURE_PA;
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
        PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
        stable_inactive_evidence_from_solution,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
        PurePhaseBoundaryStructuralTolerances, PurePhaseBoundaryTolerances,
        evaluate_pure_phase_boundary,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_solution::MultiphaseEquilibriumSolution;
    use crate::Thermodynamics::ChemEquilibrium::phase_equilibrium_workflow::{
        PhaseControlPolicy, ResolvedPhaseEquilibriumRequest, solve_resolved_pt,
    };
    use crate::Thermodynamics::ChemEquilibrium::real_pure_phase_fixtures::{
        RealPurePhaseFamily, RealPurePhaseGasScenario, RealPurePhaseInventory,
        ResolvedRealPurePhaseFixture,
    };
    use crate::Thermodynamics::phase_layout::PhaseId;
    use crate::Thermodynamics::thermo_lib_api::{ThermoData, ThermoRepository};
    use crate::library_manager::with_library_manager;

    const ROOT_LOWER_PRESSURE_PA: f64 = 100.0;
    const ROOT_UPPER_PRESSURE_PA: f64 = 1_000_000.0;
    const ROOT_MAX_ITERATIONS: usize = 96;
    const ROOT_RELATIVE_BRACKET_TOLERANCE: f64 = 1e-9;
    const MAX_INTERNAL_ROOT_RELATIVE_DELTA: f64 = 3e-6;
    const BOUNDARY_TEMPERATURES_K: [f64; 3] = [800.0, 900.0, 1_000.0];
    const ACCEPTED_BALANCE_ABSOLUTE_TOLERANCE: f64 = 1e-6;
    const ACCEPTED_BALANCE_RELATIVE_TOLERANCE: f64 = 1e-6;
    // Regression guards selected after repeated debug/release characterization.
    // They are not JANAF source uncertainty or physical acceptance tolerances.
    const MAX_EXTERNAL_BOUNDARY_RELATIVE_ERROR_GUARD: f64 = 0.015;
    const MAX_EXTERNAL_BOUNDARY_RMS_ERROR_GUARD: f64 = 0.015;

    fn janaf_directory() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("janaf")
    }

    fn janaf_dataset() -> FrozenReferenceDataset<JanafBoudouardReference> {
        let directory = janaf_directory();
        FrozenReferenceDataset::load(
            directory.join("boudouard_reaction_thermodynamics.metadata.json"),
            directory.join("boudouard_reaction_thermodynamics.rows.json"),
        )
        .expect("reviewed frozen JANAF Boudouard dataset must load")
    }

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("the bundled offline thermochemistry repository must be available")
    }

    fn boudouard_fixture() -> ResolvedRealPurePhaseFixture {
        RealPurePhaseFamily::BoudouardCarbon
            .resolve_offline(local_repository(), &[])
            .expect("the pinned local Boudouard fixture must resolve without NIST")
    }

    fn conditions(temperature_k: f64, pressure_pa: f64) -> EquilibriumConditions {
        EquilibriumConditions::new(temperature_k, pressure_pa, JANAF_STANDARD_PRESSURE_PA)
            .expect("Boudouard boundary conditions must be finite and positive")
    }

    fn gas_scenario() -> RealPurePhaseGasScenario {
        RealPurePhaseGasScenario::new(vec![BOUDOUARD_BOUNDARY_Y_CO, BOUDOUARD_BOUNDARY_Y_CO2])
            .expect("controlled 50/50 Boudouard gas must be valid")
    }

    fn evidence_request() -> PurePhaseProductionEvidenceRequest {
        PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some("gas".to_string())),
            PhaseId::new(Some("solid".to_string())),
        )
    }

    fn selected_external_references() -> Vec<JanafBoudouardBoundaryReference> {
        let dataset = janaf_dataset();
        BOUNDARY_TEMPERATURES_K
            .into_iter()
            .map(|temperature| {
                let row = dataset
                    .rows()
                    .iter()
                    .copied()
                    .find(|row| row.temperature_k == temperature)
                    .unwrap_or_else(|| {
                        panic!("frozen JANAF dataset lacks requested {temperature} K row")
                    });
                JanafBoudouardBoundaryReference::from_frozen_reference(row).unwrap_or_else(
                    |error| {
                        panic!("cannot derive JANAF Boudouard boundary at {temperature} K: {error}")
                    },
                )
            })
            .collect()
    }

    fn finite_root_value(
        label: &str,
        pressure_pa: f64,
        result: Result<f64, ReactionExtentError>,
    ) -> Result<f64, ReactionExtentError> {
        let value = result?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(ReactionExtentError::InvalidProblem {
                field: "janaf_boudouard_pressure_root",
                message: format!("{label} returned non-finite value at {pressure_pa:e} Pa"),
            })
        }
    }

    fn bisect_log_pressure<F>(
        label: &str,
        mut evaluate: F,
    ) -> Result<JanafBoudouardPressureRoot, ReactionExtentError>
    where
        F: FnMut(f64) -> Result<f64, ReactionExtentError>,
    {
        let mut lower = ROOT_LOWER_PRESSURE_PA;
        let mut upper = ROOT_UPPER_PRESSURE_PA;
        let mut lower_value = finite_root_value(label, lower, evaluate(lower))?;
        let upper_value = finite_root_value(label, upper, evaluate(upper))?;
        if lower_value == 0.0 {
            return Ok(JanafBoudouardPressureRoot {
                pressure_pa: lower,
                residual: lower_value,
                iterations: 0,
            });
        }
        if upper_value == 0.0 {
            return Ok(JanafBoudouardPressureRoot {
                pressure_pa: upper,
                residual: upper_value,
                iterations: 0,
            });
        }
        if lower_value.signum() == upper_value.signum() {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "janaf_boudouard_pressure_root",
                message: format!(
                    "{label} does not bracket a root on [{lower:e}, {upper:e}] Pa: values {lower_value:e}, {upper_value:e}"
                ),
            });
        }

        for iteration in 1..=ROOT_MAX_ITERATIONS {
            // The expected first-pass boundaries span kilopascals to bars.
            // Bisection in ln(P) preserves one relative contract everywhere.
            let pressure_pa = ((lower.ln() + upper.ln()) * 0.5).exp();
            let residual = finite_root_value(label, pressure_pa, evaluate(pressure_pa))?;
            if residual == 0.0 || (upper - lower) / pressure_pa <= ROOT_RELATIVE_BRACKET_TOLERANCE {
                return Ok(JanafBoudouardPressureRoot {
                    pressure_pa,
                    residual,
                    iterations: iteration,
                });
            }
            if lower_value.signum() != residual.signum() {
                upper = pressure_pa;
            } else {
                lower = pressure_pa;
                lower_value = residual;
            }
        }
        Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
    }

    fn validate_boundary_structure(fixture: &ResolvedRealPurePhaseFixture, temperature_k: f64) {
        let problem = fixture
            .to_pt_boundary_problem(&gas_scenario(), conditions(temperature_k, 100_000.0))
            .expect("Boudouard I1/I2 problem must materialize");
        let space = problem
            .validate_strict_independent_family(PurePhaseBoundaryStructuralTolerances::default())
            .expect("Boudouard family must retain exactly one phase-forming direction");
        assert_eq!(space.gas_only_reaction_dimension, 0);
        assert_eq!(space.full_reaction_dimension, 1);
    }

    fn independent_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature_k: f64,
    ) -> Result<JanafBoudouardPressureRoot, ReactionExtentError> {
        validate_boundary_structure(fixture, temperature_k);
        bisect_log_pressure("I1/I2 ln(Q)-ln(K)", |pressure_pa| {
            let problem = fixture
                .to_pt_boundary_problem(&gas_scenario(), conditions(temperature_k, pressure_pa))?;
            Ok(
                evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default())?
                    .log_residual_at_absence,
            )
        })
    }

    fn graphite_tpd_at_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<f64, ReactionExtentError> {
        let inventory = RealPurePhaseInventory::from_gas(gas_scenario(), 0.0)?;
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(temperature_k, pressure_pa),
                fixture.initial_composition(&inventory)?,
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )?;
        let request = evidence_request();
        let evidence = match solution.phase_status(&request.candidate_phase) {
            Some(PhaseStatus::Inactive) => {
                stable_inactive_evidence_from_solution(&solution, &request)?
            }
            Some(PhaseStatus::Active | PhaseStatus::Appeared) => {
                activation_evidence_from_solution(&solution, &request)?
            }
            status => {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "janaf_boudouard_tpd_root",
                    message: format!(
                        "graphite has unsupported phase status {status:?} during TPD root search"
                    ),
                });
            }
        };
        evidence
            .boundary_minimum_tpd
            .filter(|value| value.is_finite())
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "janaf_boudouard_tpd_root",
                message: "canonical graphite evidence lacks a finite boundary TPD".to_string(),
            })
    }

    fn canonical_tpd_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature_k: f64,
    ) -> Result<JanafBoudouardPressureRoot, ReactionExtentError> {
        bisect_log_pressure("canonical TPD(C | gas)", |pressure_pa| {
            graphite_tpd_at_pressure(fixture, temperature_k, pressure_pa)
        })
    }

    fn solve_lifecycle(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> MultiphaseEquilibriumSolution {
        let inventory = RealPurePhaseInventory::from_gas(gas_scenario(), 0.0)
            .expect("gas-only Boudouard lifecycle inventory must validate");
        solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(temperature_k, pressure_pa),
                fixture
                    .initial_composition(&inventory)
                    .expect("lifecycle inventory must match Boudouard layout"),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .expect("Boudouard phase-control lifecycle must solve")
    }

    fn local_library_snapshot() -> Vec<(String, Vec<u8>)> {
        with_library_manager(|manager| {
            vec![
                manager.substance_base_path().to_string(),
                manager.all_keys_substance_path().to_string(),
                manager.elements_path().to_string(),
            ]
        })
        .into_iter()
        .map(|path| {
            (
                path.clone(),
                fs::read(&path).expect("local library must be readable"),
            )
        })
        .collect()
    }

    fn frozen_janaf_snapshot() -> (Vec<u8>, Vec<u8>) {
        let directory = janaf_directory();
        (
            fs::read(directory.join("boudouard_reaction_thermodynamics.metadata.json"))
                .expect("frozen JANAF metadata must be readable"),
            fs::read(directory.join("boudouard_reaction_thermodynamics.rows.json"))
                .expect("frozen JANAF rows must be readable"),
        )
    }

    fn assert_janaf_boundary_external_regression_envelope(report: &JanafBoudouardBoundaryReport) {
        assert!(
            report.max_external_relative_error() <= MAX_EXTERNAL_BOUNDARY_RELATIVE_ERROR_GUARD,
            "JANAF Boudouard boundary external quality regression: max relative error={:e} exceeds reviewed guard={MAX_EXTERNAL_BOUNDARY_RELATIVE_ERROR_GUARD:e}",
            report.max_external_relative_error(),
        );
        let rms = report.rms_independent_external_error();
        assert!(
            rms <= MAX_EXTERNAL_BOUNDARY_RMS_ERROR_GUARD,
            "JANAF Boudouard boundary external quality regression: I1/I2 RMS relative error={rms:e} exceeds reviewed guard={MAX_EXTERNAL_BOUNDARY_RMS_ERROR_GUARD:e}",
        );
    }

    #[test]
    fn i5_janaf_boudouard_analytical_boundary_uses_selected_primary_rows() {
        let boundaries = selected_external_references();
        assert_eq!(
            boundaries
                .iter()
                .map(|row| row.temperature_k)
                .collect::<Vec<_>>(),
            BOUNDARY_TEMPERATURES_K
        );
        for boundary in boundaries {
            boundary.validate().unwrap();
            assert_eq!(boundary.gas_y_co, 0.5);
            assert_eq!(boundary.gas_y_co2, 0.5);
        }
    }

    #[test]
    fn i5_janaf_boudouard_independent_pressure_obeys_delta_nu_minus_one() {
        let fixture = boudouard_fixture();
        let temperature_k = 900.0;
        let pressure_1_pa = 10_000.0;
        let pressure_2_pa = 250_000.0;
        let residual = |pressure_pa| {
            let problem = fixture
                .to_pt_boundary_problem(&gas_scenario(), conditions(temperature_k, pressure_pa))
                .unwrap();
            evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default())
                .unwrap()
                .log_residual_at_absence
        };
        let delta = residual(pressure_2_pa) - residual(pressure_1_pa);
        assert!(
            (delta + (pressure_2_pa / pressure_1_pa).ln()).abs() <= 1e-11,
            "Boudouard lnQ pressure invariant failed: delta={delta:e}"
        );
    }

    #[test]
    fn i5_janaf_boudouard_internal_pressure_roots_agree() {
        let before = local_library_snapshot();
        let frozen_before = frozen_janaf_snapshot();
        let fixture = boudouard_fixture();
        let rows = selected_external_references()
            .into_iter()
            .map(|external| {
                let independent = independent_boundary_pressure(&fixture, external.temperature_k)
                    .unwrap_or_else(|error| {
                        panic!(
                            "I1/I2 Boudouard root failed at {} K: {error}",
                            external.temperature_k
                        )
                    });
                let canonical_tpd =
                    canonical_tpd_boundary_pressure(&fixture, external.temperature_k)
                        .unwrap_or_else(|error| {
                            panic!(
                                "TPD Boudouard root failed at {} K: {error}",
                                external.temperature_k
                            )
                        });
                JanafBoudouardBoundaryRow {
                    external,
                    independent,
                    canonical_tpd,
                }
            })
            .collect();
        let report = JanafBoudouardBoundaryReport::new(rows).unwrap();
        assert_eq!(report.rows().len(), BOUNDARY_TEMPERATURES_K.len());
        report
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap_or_else(|error| panic!("Boudouard internal root mismatch: {error}"));
        assert_janaf_boundary_external_regression_envelope(&report);
        assert_eq!(before, local_library_snapshot());
        assert_eq!(frozen_before, frozen_janaf_snapshot());
    }

    #[test]
    fn i5_janaf_boudouard_local_lifecycle_follows_local_pressure_boundary() {
        let before = local_library_snapshot();
        let frozen_before = frozen_janaf_snapshot();
        let fixture = boudouard_fixture();
        let boundary = independent_boundary_pressure(&fixture, 900.0)
            .expect("local 900 K Boudouard boundary must be available");
        let request = evidence_request();

        // Delta nu_gas = -1: higher pressure favours CO2 + graphite. This
        // uses the local I1/I2 boundary only; JANAF remains external evidence.
        let appeared = solve_lifecycle(&fixture, 900.0, boundary.pressure_pa * 3.0);
        assert!(matches!(
            appeared.phase_status(&request.candidate_phase),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        ));
        let activation = activation_evidence_from_solution(&appeared, &request)
            .expect("high-pressure graphite must retain activation evidence");
        assert!(activation.boundary_minimum_tpd.unwrap() < 0.0);

        let inactive = solve_lifecycle(&fixture, 900.0, boundary.pressure_pa / 3.0);
        assert_eq!(
            inactive.phase_status(&request.candidate_phase),
            Some(PhaseStatus::Inactive)
        );
        let absence = stable_inactive_evidence_from_solution(&inactive, &request)
            .expect("low-pressure graphite must retain stable-inactive evidence");
        assert!(absence.boundary_minimum_tpd.unwrap() > 0.0);

        for solution in [&appeared, &inactive] {
            let acceptance = solution
                .acceptance_report()
                .expect("bounded phase-control solution must publish acceptance");
            // Match the accepted candidate contract: conservation is absolute
            // plus relative to the original element scale, never an arbitrary
            // tighter assertion introduced by this test fixture.
            let initial_element_scale = 1.5;
            let accepted_limit = ACCEPTED_BALANCE_ABSOLUTE_TOLERANCE
                + ACCEPTED_BALANCE_RELATIVE_TOLERANCE * initial_element_scale;
            assert!(
                acceptance.final_validation.max_abs_element_balance_error <= accepted_limit,
                "Boudouard lifecycle conservation exceeded the production-scale contract: error={:e}, limit={accepted_limit:e}",
                acceptance.final_validation.max_abs_element_balance_error,
            );
        }
        assert_eq!(before, local_library_snapshot());
        assert_eq!(frozen_before, frozen_janaf_snapshot());
    }

    #[test]
    #[ignore = "release I5 JANAF Boudouard pressure characterization with conservative software-regression envelope"]
    fn i5_janaf_boudouard_pressure_boundary_diagnostic() {
        let frozen_before = frozen_janaf_snapshot();
        let libraries_before = local_library_snapshot();
        let fixture = boudouard_fixture();
        println!("KiThe standard-state pressure provenance:");
        for row in fixture.thermochemistry().provenance() {
            println!(
                "  {:<16} {}",
                row.component().label(),
                row.standard_state_pressure()
            );
        }
        let rows = selected_external_references()
            .into_iter()
            .map(|external| JanafBoudouardBoundaryRow {
                independent: independent_boundary_pressure(&fixture, external.temperature_k)
                    .unwrap_or_else(|error| {
                        panic!("I1/I2 failed at {} K: {error}", external.temperature_k)
                    }),
                canonical_tpd: canonical_tpd_boundary_pressure(&fixture, external.temperature_k)
                    .unwrap_or_else(|error| {
                        panic!("TPD failed at {} K: {error}", external.temperature_k)
                    }),
                external,
            })
            .collect();
        let report = JanafBoudouardBoundaryReport::new(rows).unwrap();
        report
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap();
        assert_janaf_boundary_external_regression_envelope(&report);
        println!("{report}");
        assert_eq!(frozen_janaf_snapshot(), frozen_before);
        assert_eq!(local_library_snapshot(), libraries_before);
    }
}
