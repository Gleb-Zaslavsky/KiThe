//! Frozen-IAPWS ice-Ih sublimation evidence on the real local water records.
//!
//! The I1/I2 route, gas-reference TPD route, and ice-reference TPD route are
//! intentionally separate. Their strict agreement validates KiThe's internal
//! formulation; the IAPWS comparison is a visible characterization of the
//! current NASA/ideal-phase model, not a synthetic acceptance tolerance.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        ReactionExtentError, SolveError,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        PhaseStatus, PhaseTransitionReason,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, WaterIceSublimationPressureReference,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_iapws_ice_sublimation::{
        IapwsIceSublimationComparisonContract, IceSublimationComparisonReport,
        IceSublimationComparisonRow, IceSublimationRoot,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_production_adapter::{
        PurePhaseProductionEvidenceRequest, activation_evidence_from_solution,
        stable_inactive_evidence_from_solution,
    };
    use crate::Thermodynamics::ChemEquilibrium::phase_boundary_validation::{
        PurePhaseBoundaryTolerances, evaluate_pure_phase_boundary,
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

    const REFERENCE_PRESSURE_PA: f64 = 101_325.0;
    const ROOT_LOWER_PRESSURE_PA: f64 = 0.01;
    const ROOT_UPPER_PRESSURE_PA: f64 = 1_000.0;
    const ROOT_MAX_ITERATIONS: usize = 112;
    const ROOT_RELATIVE_BRACKET_TOLERANCE: f64 = 1e-10;
    const MAX_INTERNAL_ROOT_RELATIVE_DELTA: f64 = 2e-7;
    // Regression guards selected after repeated debug/release characterization.
    // They are not IAPWS source uncertainty or physical acceptance tolerances.
    const MAX_EXTERNAL_RELATIVE_ERROR_GUARD: f64 = 0.025;
    const MAX_EXTERNAL_RMS_ERROR_GUARD: f64 = 0.025;

    fn iapws_directory() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("iapws")
    }

    fn iapws_dataset() -> FrozenReferenceDataset<WaterIceSublimationPressureReference> {
        let directory = iapws_directory();
        FrozenReferenceDataset::load(
            directory.join("water_ice_sublimation_low_pressure.metadata.json"),
            directory.join("water_ice_sublimation_low_pressure.rows.json"),
        )
        .expect("reviewed IAPWS ice-Ih frozen dataset must load")
    }

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("the bundled offline thermochemistry repository must be available")
    }

    fn conditions(temperature: f64, pressure_pa: f64) -> EquilibriumConditions {
        EquilibriumConditions::new(temperature, pressure_pa, REFERENCE_PRESSURE_PA)
            .expect("IAPWS ice diagnostic conditions must be finite and positive")
    }

    fn water_gas_scenario() -> RealPurePhaseGasScenario {
        RealPurePhaseGasScenario::new(vec![1.0])
            .expect("direct pure-water gas boundary scenario must be valid")
    }

    fn bisect_log_pressure<F>(
        label: &str,
        mut evaluate: F,
    ) -> Result<(f64, f64, usize), ReactionExtentError>
    where
        F: FnMut(f64) -> Result<f64, ReactionExtentError>,
    {
        let mut lower = ROOT_LOWER_PRESSURE_PA;
        let mut upper = ROOT_UPPER_PRESSURE_PA;
        let mut lower_value = finite_root_value(label, lower, evaluate(lower))?;
        let upper_value = finite_root_value(label, upper, evaluate(upper))?;
        if lower_value == 0.0 {
            return Ok((lower, lower_value, 0));
        }
        if upper_value == 0.0 {
            return Ok((upper, upper_value, 0));
        }
        if lower_value.signum() == upper_value.signum() {
            return Err(ReactionExtentError::ValidationNotApplicable {
                path: "iapws_ice_pressure_root",
                message: format!(
                    "{label} does not bracket a root on [{lower:e}, {upper:e}] Pa: values {lower_value:e}, {upper_value:e}"
                ),
            });
        }

        for iteration in 1..=ROOT_MAX_ITERATIONS {
            // The table spans four orders of magnitude; bisection in ln(P)
            // gives the same relative stopping semantics at every row.
            let pressure = ((lower.ln() + upper.ln()) * 0.5).exp();
            let value = finite_root_value(label, pressure, evaluate(pressure))?;
            if value == 0.0 || (upper - lower) / pressure <= ROOT_RELATIVE_BRACKET_TOLERANCE {
                return Ok((pressure, value, iteration));
            }
            if lower_value.signum() != value.signum() {
                upper = pressure;
            } else {
                lower = pressure;
                lower_value = value;
            }
        }

        Err(ReactionExtentError::SolveError(SolveError::MaxIterations))
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
                field: "iapws_ice_pressure_root",
                message: format!("{label} returned non-finite value at {pressure_pa:e} Pa"),
            })
        }
    }

    fn independent_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
    ) -> Result<IceSublimationRoot, ReactionExtentError> {
        let gas = water_gas_scenario();
        let (pressure_pa, residual, iterations) = bisect_log_pressure("I1/I2 ln(Q)-ln(K)", |p| {
            let problem = fixture.to_pt_boundary_problem(&gas, conditions(temperature, p))?;
            Ok(
                evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default())?
                    .log_residual_at_absence,
            )
        })?;
        Ok(IceSublimationRoot {
            pressure_pa,
            residual,
            iterations,
        })
    }

    fn ice_tpd_at_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
        pressure_pa: f64,
    ) -> Result<f64, ReactionExtentError> {
        let inventory = RealPurePhaseInventory::from_gas(water_gas_scenario(), 0.0)?;
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(temperature, pressure_pa),
                fixture.initial_composition(&inventory)?,
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )?;
        let request = PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some("gas".to_owned())),
            PhaseId::new(Some("solid".to_owned())),
        );
        let ice = request.candidate_phase.clone();
        let evidence = match solution.phase_status(&ice) {
            Some(PhaseStatus::Inactive) => {
                stable_inactive_evidence_from_solution(&solution, &request)?
            }
            Some(PhaseStatus::Active | PhaseStatus::Appeared) => {
                activation_evidence_from_solution(&solution, &request)?
            }
            status => {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "iapws_ice_tpd_root",
                    message: format!(
                        "ice has unsupported phase status {status:?} during TPD root search"
                    ),
                });
            }
        };
        evidence
            .boundary_minimum_tpd
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "iapws_ice_tpd_root",
                message: "canonical ice evidence lacks a boundary TPD".to_owned(),
            })
    }

    fn ice_from_gas_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
    ) -> Result<IceSublimationRoot, ReactionExtentError> {
        let (pressure_pa, residual, iterations) = bisect_log_pressure("TPD(ice | gas)", |p| {
            ice_tpd_at_pressure(fixture, temperature, p)
        })?;
        Ok(IceSublimationRoot {
            pressure_pa,
            residual,
            iterations,
        })
    }

    fn solve_from_ice_reference(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
        pressure_pa: f64,
    ) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
        // The gas trace is solely a finite log-coordinate seed. The physical
        // accepted inventory starts as pure ice, so gas TPD is evaluated from
        // the correct condensed reference when no gas assemblage is active.
        let inventory = RealPurePhaseInventory::new(vec![1e-300], 1.0)?;
        solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(temperature, pressure_pa),
                fixture.initial_composition(&inventory)?,
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
    }

    fn gas_tpd_from_ice_reference(
        solution: &MultiphaseEquilibriumSolution,
    ) -> Result<f64, ReactionExtentError> {
        let gas = PhaseId::new(Some("gas".to_owned()));
        let gas_index = phase_index(solution, &gas)?;
        match solution.phase_status(&gas) {
            Some(PhaseStatus::Inactive) => solution
                .acceptance_report()
                .and_then(|report| report.phase_stability.get(gas_index))
                .and_then(|report| report.minimum_tpd)
                .filter(|value| value.is_finite())
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "iapws_ice_gas_tpd",
                    message: "stable absent gas lacks final finite TPD evidence".to_owned(),
                }),
            Some(PhaseStatus::Active | PhaseStatus::Appeared) => solution
                .phase_control_report()
                .and_then(|report| {
                    report
                        .transitions
                        .iter()
                        .find(|transition| {
                            transition
                                .activated
                                .iter()
                                .any(|phase| phase.index() == gas_index)
                        })
                        .and_then(|transition| transition.minimum_tpds.get(gas_index))
                        .and_then(|value| *value)
                })
                .filter(|value| value.is_finite())
                .ok_or_else(|| ReactionExtentError::InvalidProblem {
                    field: "iapws_ice_gas_tpd",
                    message: "appearing gas lacks finite pre-activation TPD evidence".to_owned(),
                }),
            status => Err(ReactionExtentError::InvalidProblem {
                field: "iapws_ice_gas_tpd",
                message: format!("gas has unsupported final phase status {status:?}"),
            }),
        }
    }

    fn gas_from_ice_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
    ) -> Result<IceSublimationRoot, ReactionExtentError> {
        let (pressure_pa, residual, iterations) = bisect_log_pressure("TPD(gas | ice)", |p| {
            let solution = solve_from_ice_reference(fixture, temperature, p)?;
            gas_tpd_from_ice_reference(&solution)
        })?;
        Ok(IceSublimationRoot {
            pressure_pa,
            residual,
            iterations,
        })
    }

    fn phase_index(
        solution: &MultiphaseEquilibriumSolution,
        phase_id: &PhaseId,
    ) -> Result<usize, ReactionExtentError> {
        solution
            .phases()
            .iter()
            .position(|phase| phase.id() == phase_id)
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "iapws_ice_phase_layout",
                message: format!("accepted solution lacks phase {phase_id:?}"),
            })
    }

    #[test]
    fn i5_iapws_ice_sublimation_dataset_matches_local_validity_interval() {
        let dataset = iapws_dataset();
        let fixture = RealPurePhaseFamily::WaterIce
            .resolve_offline(local_repository(), &[])
            .expect("local H2O(g)/H2O(s) records must resolve without NIST");
        let bounds = fixture.thermochemistry().temperature_bounds();
        let applicable = dataset
            .rows()
            .iter()
            .filter(|row| bounds.contains(row.temperature_k))
            .count();

        IapwsIceSublimationComparisonContract::characterization_only()
            .validate_dataset(&dataset)
            .unwrap();
        assert_eq!(dataset.rows().len(), 8);
        assert_eq!(
            applicable,
            dataset.rows().len(),
            "local bounds are [{}, {}] K",
            bounds.lower(),
            bounds.upper()
        );
        assert_eq!(dataset.rows().first().unwrap().temperature_k, 200.0);
        assert_eq!(dataset.rows().last().unwrap().temperature_k, 270.0);
    }

    #[test]
    fn i5_iapws_ice_trace_pressure_root_stays_finite() {
        let fixture = RealPurePhaseFamily::WaterIce
            .resolve_offline(local_repository(), &[])
            .expect("local ice fixture must resolve");
        let independent = independent_boundary_pressure(&fixture, 200.0)
            .expect("200 K independent ice boundary must be finite");
        let ice_from_gas = ice_from_gas_boundary_pressure(&fixture, 200.0)
            .expect("200 K ice TPD boundary must be finite");
        let gas_from_ice = gas_from_ice_boundary_pressure(&fixture, 200.0)
            .expect("200 K gas TPD boundary must be finite");
        for root in [independent, ice_from_gas, gas_from_ice] {
            assert!(root.pressure_pa.is_finite() && root.pressure_pa > 0.0);
            assert!(root.residual.is_finite());
        }
    }

    #[test]
    fn i5_direct_water_ice_lifecycle_tracks_internal_boundary() {
        let fixture = RealPurePhaseFamily::WaterIce
            .resolve_offline(local_repository(), &[])
            .expect("local ice fixture must resolve");
        let boundary = independent_boundary_pressure(&fixture, 250.0)
            .expect("internal 250 K ice boundary must be available");
        let deposited = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(250.0, boundary.pressure_pa * 3.0),
                fixture
                    .initial_composition(
                        &RealPurePhaseInventory::from_gas(water_gas_scenario(), 0.0).unwrap(),
                    )
                    .unwrap(),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .expect("pressure above the local sublimation boundary must deposit ice");
        let gas = PhaseId::new(Some("gas".to_owned()));
        let ice = PhaseId::new(Some("solid".to_owned()));
        assert_eq!(deposited.phase_status(&gas), Some(PhaseStatus::Inactive));
        assert!(matches!(
            deposited.phase_status(&ice),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        ));
        assert!(deposited.phase_control_report().unwrap().transitions.iter().any(|transition| {
            transition.deactivated.iter().any(|phase| phase.index() == 0)
                && matches!(transition.reason, PhaseTransitionReason::BoundaryUnstableActivePhase { minimum_tpd, .. } if minimum_tpd > 0.0)
        }));

        let sublimated = solve_from_ice_reference(&fixture, 250.0, boundary.pressure_pa / 3.0)
            .expect("pressure below the local sublimation boundary must sublimate ice");
        assert!(matches!(
            sublimated.phase_status(&gas),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        ));
        assert_eq!(sublimated.phase_status(&ice), Some(PhaseStatus::Inactive));
        assert!(gas_tpd_from_ice_reference(&sublimated).unwrap() < 0.0);
        for solution in [&deposited, &sublimated] {
            assert!((solution.component_moles().iter().sum::<f64>() - 1.0).abs() <= 1e-10);
        }
    }

    fn assert_iapws_ice_external_regression_envelope(report: &IceSublimationComparisonReport) {
        assert!(
            report.max_external_relative_error() <= MAX_EXTERNAL_RELATIVE_ERROR_GUARD,
            "IAPWS ice-Ih external quality regression: max relative error={:e} exceeds reviewed guard={MAX_EXTERNAL_RELATIVE_ERROR_GUARD:e}",
            report.max_external_relative_error(),
        );
        let rms = report.independent_external_rms();
        assert!(
            rms <= MAX_EXTERNAL_RMS_ERROR_GUARD,
            "IAPWS ice-Ih external quality regression: independent RMS relative error={rms:e} exceeds reviewed guard={MAX_EXTERNAL_RMS_ERROR_GUARD:e}",
        );
    }

    #[test]
    #[ignore = "release I5 IAPWS ice-Ih characterization with conservative software-regression envelope"]
    fn i5_iapws_water_ice_ih_low_pressure_boundary_diagnostic() {
        let directory = iapws_directory();
        let metadata_path = directory.join("water_ice_sublimation_low_pressure.metadata.json");
        let rows_path = directory.join("water_ice_sublimation_low_pressure.rows.json");
        let metadata_before = fs::read(&metadata_path).unwrap();
        let rows_before = fs::read(&rows_path).unwrap();
        let dataset = iapws_dataset();
        let fixture = RealPurePhaseFamily::WaterIce
            .resolve_offline(local_repository(), &[])
            .expect("IAPWS ice diagnostic must resolve local NASA data without NIST");
        let bounds = fixture.thermochemistry().temperature_bounds();
        let excluded = dataset
            .rows()
            .iter()
            .filter(|row| !bounds.contains(row.temperature_k))
            .collect::<Vec<_>>();
        assert!(
            excluded.is_empty(),
            "local validity [{}, {}] K excludes frozen IAPWS rows {:?}",
            bounds.lower(),
            bounds.upper(),
            excluded
                .iter()
                .map(|row| row.temperature_k)
                .collect::<Vec<_>>()
        );

        let contract = IapwsIceSublimationComparisonContract::characterization_only();
        contract.validate_dataset(&dataset).unwrap();
        let mut rows = Vec::with_capacity(dataset.rows().len());
        for row in dataset.rows() {
            let independent = independent_boundary_pressure(&fixture, row.temperature_k)
                .unwrap_or_else(|error| panic!("I1/I2 failed at {} K: {error}", row.temperature_k));
            let ice_from_gas = ice_from_gas_boundary_pressure(&fixture, row.temperature_k)
                .unwrap_or_else(|error| {
                    panic!("TPD(ice|gas) failed at {} K: {error}", row.temperature_k)
                });
            let gas_from_ice = gas_from_ice_boundary_pressure(&fixture, row.temperature_k)
                .unwrap_or_else(|error| {
                    panic!("TPD(gas|ice) failed at {} K: {error}", row.temperature_k)
                });
            rows.push(IceSublimationComparisonRow {
                temperature_k: row.temperature_k,
                iapws_pressure_pa: row.pressure_pa,
                independent,
                ice_from_gas,
                gas_from_ice,
            });
        }
        let report = IceSublimationComparisonReport::new(&dataset, rows).unwrap();
        contract.validate_report(&report).unwrap();
        assert_eq!(report.rows().len(), dataset.rows().len());
        report
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap_or_else(|error| {
                panic!("strict ice three-way internal contract failed: {error}")
            });
        assert_iapws_ice_external_regression_envelope(&report);
        println!("{report}");
        assert_eq!(fs::read(metadata_path).unwrap(), metadata_before);
        assert_eq!(fs::read(rows_path).unwrap(), rows_before);
    }
}
