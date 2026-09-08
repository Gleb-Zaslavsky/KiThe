//! Frozen-IAPWS evidence for the low-pressure liquid/vapor boundary of water.
//!
//! The normal test below checks only the frozen external dataset. The ignored
//! diagnostic intentionally performs no external pass/fail comparison yet:
//! it first characterizes the difference between IAPWS values and the current
//! KiThe ideal-gas/pure-condensed model. By contrast, it keeps I1/I2 and TPD
//! roots on a strict internal agreement contract.

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::sync::Arc;

    use serde_json::{Value, json};

    use crate::Thermodynamics::ChemEquilibrium::equilibrium_nonlinear::{
        ReactionExtentError, SolveError,
    };
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_problem::EquilibriumConditions;
    use crate::Thermodynamics::ChemEquilibrium::equilibrium_workflows::{
        PhaseStabilityStatus, PhaseStatus, PhaseTransitionReason,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference::{
        FrozenReferenceDataset, WaterSaturationPressureReference,
    };
    use crate::Thermodynamics::ChemEquilibrium::frozen_reference_iapws_water_saturation::{
        ExternalComparisonMode, IapwsWaterSaturationComparisonContract,
        WaterSaturationComparisonReport, WaterSaturationComparisonRow, WaterSaturationRoot,
        WaterSaturationSymmetryReport, WaterSaturationSymmetryRow,
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
    const ROOT_LOWER_PRESSURE_PA: f64 = 1.0;
    const ROOT_UPPER_PRESSURE_PA: f64 = 1_000_000.0;
    const ROOT_MAX_ITERATIONS: usize = 96;
    const ROOT_RELATIVE_BRACKET_TOLERANCE: f64 = 1e-10;
    const MAX_INTERNAL_ROOT_RELATIVE_DELTA: f64 = 2e-7;
    // Regression guards selected after repeated debug/release characterization.
    // They are not IAPWS source uncertainty or physical acceptance tolerances.
    const MAX_EXTERNAL_RELATIVE_ERROR_GUARD: f64 = 0.02;
    const MAX_EXTERNAL_RMS_ERROR_GUARD: f64 = 0.02;

    fn iapws_directory() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("Thermodynamics")
            .join("ChemEquilibrium")
            .join("frozen_reference")
            .join("data")
            .join("iapws")
    }

    fn iapws_dataset() -> FrozenReferenceDataset<WaterSaturationPressureReference> {
        let directory = iapws_directory();
        FrozenReferenceDataset::load(
            directory.join("water_liquid_saturation_low_pressure.metadata.json"),
            directory.join("water_liquid_saturation_low_pressure.rows.json"),
        )
        .expect("reviewed IAPWS frozen dataset must load")
    }

    fn parse_iapws_rows(
        rows: &Value,
    ) -> Result<FrozenReferenceDataset<WaterSaturationPressureReference>, String> {
        let directory = iapws_directory();
        let metadata = fs::read_to_string(
            directory.join("water_liquid_saturation_low_pressure.metadata.json"),
        )
        .map_err(|error| error.to_string())?;
        FrozenReferenceDataset::from_json_strs(
            &metadata,
            &serde_json::to_string(rows).map_err(|error| error.to_string())?,
            "water_liquid_saturation_low_pressure.rows.json",
        )
        .map_err(|error| error.to_string())
    }

    fn local_repository() -> Arc<ThermoRepository> {
        ThermoData::try_default_repository()
            .expect("the bundled offline thermochemistry repository must be available")
    }

    fn conditions(temperature: f64, pressure_pa: f64) -> EquilibriumConditions {
        EquilibriumConditions::new(temperature, pressure_pa, REFERENCE_PRESSURE_PA)
            .expect("IAPWS diagnostic conditions must be finite and positive")
    }

    fn water_gas_scenario() -> RealPurePhaseGasScenario {
        RealPurePhaseGasScenario::new(vec![1.0])
            .expect("direct pure-water gas boundary scenario must be valid")
    }

    fn carrier_water_gas_scenario() -> RealPurePhaseGasScenario {
        water_gas_scenario()
            // This non-reacting carrier exists only for the metamorphic check
            // below. The production IAPWS diagnostic uses direct H2O; adding
            // O2 must not move the independently computed p_H2O boundary.
            .with_inert("O2", 1.0, vec![0.0, 2.0])
            .expect("pure-water gas boundary scenario must be valid")
    }

    fn water_partial_pressure(pressure_pa: f64, gas: &RealPurePhaseGasScenario) -> f64 {
        let moles = gas.initial_gas_moles();
        pressure_pa * moles[0] / moles.iter().sum::<f64>()
    }

    fn bisect_pressure<F>(
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
                path: "iapws_water_pressure_root",
                message: format!(
                    "{label} does not bracket a root on [{lower:e}, {upper:e}] Pa: values {lower_value:e}, {upper_value:e}"
                ),
            });
        }

        for iteration in 1..=ROOT_MAX_ITERATIONS {
            // Saturation pressure spans orders of magnitude, so bisection in
            // log-pressure resolves every table row uniformly.
            let pressure = (lower.ln() + upper.ln()).mul_add(0.5, 0.0).exp();
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
                field: "iapws_water_pressure_root",
                message: format!("{label} returned non-finite value at {pressure_pa:e} Pa"),
            })
        }
    }

    fn independent_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        gas: &RealPurePhaseGasScenario,
        temperature: f64,
    ) -> Result<WaterSaturationRoot, ReactionExtentError> {
        let (total_pressure_pa, residual, iterations) =
            bisect_pressure("I1/I2 ln(Q)-ln(K)", |pressure_pa| {
                let problem =
                    fixture.to_pt_boundary_problem(gas, conditions(temperature, pressure_pa))?;
                Ok(
                    evaluate_pure_phase_boundary(&problem, PurePhaseBoundaryTolerances::default())?
                        .log_residual_at_absence,
                )
            })?;
        Ok(WaterSaturationRoot {
            total_pressure_pa,
            partial_water_pressure_pa: water_partial_pressure(total_pressure_pa, gas),
            residual,
            iterations,
        })
    }

    fn canonical_tpd_at_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        inventory: &RealPurePhaseInventory,
        temperature: f64,
        pressure_pa: f64,
    ) -> Result<f64, ReactionExtentError> {
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(temperature, pressure_pa),
                fixture.initial_composition(inventory)?,
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )?;
        let request = PurePhaseProductionEvidenceRequest::new(
            PhaseId::new(Some("gas".to_owned())),
            PhaseId::new(Some("liquid".to_owned())),
        );
        let liquid = request.candidate_phase.clone();
        let evidence = match solution.phase_status(&liquid) {
            Some(PhaseStatus::Inactive) => {
                stable_inactive_evidence_from_solution(&solution, &request)?
            }
            Some(PhaseStatus::Active | PhaseStatus::Appeared) => {
                activation_evidence_from_solution(&solution, &request)?
            }
            status => {
                return Err(ReactionExtentError::InvalidProblem {
                    field: "iapws_water_tpd_root",
                    message: format!(
                        "liquid water has unsupported phase status {status:?} during TPD root search"
                    ),
                });
            }
        };
        evidence
            .boundary_minimum_tpd
            .ok_or_else(|| ReactionExtentError::InvalidProblem {
                field: "iapws_water_tpd_root",
                message: "canonical phase-control evidence lacks a boundary TPD".to_owned(),
            })
    }

    /// Solves from a condensed-only physical inventory. The trace gas amount
    /// is a log-coordinate seed, not an active phase inventory.
    fn solve_from_liquid_reference(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
        pressure_pa: f64,
    ) -> Result<MultiphaseEquilibriumSolution, ReactionExtentError> {
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

    /// Extracts the gas driving force evaluated against the liquid-only
    /// reference state. A stable absent gas uses final TPD evidence; an
    /// evaporating gas uses the immutable pre-activation transition evidence.
    fn gas_tpd_from_liquid_reference(
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
                    field: "iapws_water_gas_tpd",
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
                    field: "iapws_water_gas_tpd",
                    message: "appearing gas lacks finite pre-activation TPD evidence".to_owned(),
                }),
            status => Err(ReactionExtentError::InvalidProblem {
                field: "iapws_water_gas_tpd",
                message: format!("gas has unsupported final phase status {status:?}"),
            }),
        }
    }

    fn condensed_reference_gas_tpd_at_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
        pressure_pa: f64,
    ) -> Result<f64, ReactionExtentError> {
        let solution = solve_from_liquid_reference(fixture, temperature, pressure_pa)?;
        gas_tpd_from_liquid_reference(&solution)
    }

    fn condensed_side_tpd_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        temperature: f64,
    ) -> Result<WaterSaturationRoot, ReactionExtentError> {
        let (total_pressure_pa, residual, iterations) =
            bisect_pressure("gas TPD from liquid reference", |pressure_pa| {
                condensed_reference_gas_tpd_at_pressure(fixture, temperature, pressure_pa)
            })?;
        Ok(WaterSaturationRoot {
            total_pressure_pa,
            partial_water_pressure_pa: total_pressure_pa,
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
                field: "iapws_water_phase_layout",
                message: format!("accepted solution lacks phase {phase_id:?}"),
            })
    }

    fn canonical_tpd_boundary_pressure(
        fixture: &ResolvedRealPurePhaseFixture,
        inventory: &RealPurePhaseInventory,
        gas: &RealPurePhaseGasScenario,
        temperature: f64,
    ) -> Result<WaterSaturationRoot, ReactionExtentError> {
        let (total_pressure_pa, residual, iterations) =
            bisect_pressure("canonical TPD", |pressure_pa| {
                canonical_tpd_at_pressure(fixture, inventory, temperature, pressure_pa)
            })?;
        Ok(WaterSaturationRoot {
            total_pressure_pa,
            partial_water_pressure_pa: water_partial_pressure(total_pressure_pa, gas),
            residual,
            iterations,
        })
    }

    #[test]
    fn i5_iapws_low_pressure_water_saturation_dataset_is_valid() {
        let dataset = iapws_dataset();
        let contract = IapwsWaterSaturationComparisonContract::characterization_only();
        assert!(dataset.is_external_evidence());
        assert_eq!(dataset.rows().len(), 10);
        assert_eq!(dataset.rows().first().unwrap().temperature_k, 275.0);
        assert_eq!(dataset.rows().last().unwrap().temperature_k, 320.0);
        assert_eq!(
            dataset.metadata().source.organization,
            "International Association for the Properties of Water and Steam (IAPWS)"
        );
        assert_eq!(
            dataset.metadata().source.version.as_deref(),
            Some("SR1-86(1992)")
        );
        assert_eq!(
            contract.mode(),
            ExternalComparisonMode::CharacterizationOnly
        );
        assert!(contract.model_scope().contains("NASA"));
        contract.validate_dataset(&dataset).unwrap();
    }

    #[test]
    fn i5_iapws_water_saturation_contract_rejects_non_liquid_and_non_monotonic_rows() {
        let directory = iapws_directory();
        let rows_path = directory.join("water_liquid_saturation_low_pressure.rows.json");
        let original: Value = serde_json::from_slice(&fs::read(rows_path).unwrap()).unwrap();

        let mut triple_point = original.clone();
        triple_point["rows"][0]["temperature_k"] = json!(273.16);
        let error = parse_iapws_rows(&triple_point).unwrap_err();
        assert!(error.contains("above 273.16 K"), "{error}");

        let mut non_monotonic_pressure = original;
        non_monotonic_pressure["rows"][1]["pressure_pa"] = json!(500.0);
        let error = parse_iapws_rows(&non_monotonic_pressure).unwrap_err();
        assert!(
            error.contains("pressures must be strictly increasing"),
            "{error}"
        );
    }

    fn assert_iapws_external_regression_envelope(report: &WaterSaturationComparisonReport) {
        assert!(
            report.max_external_relative_error() <= MAX_EXTERNAL_RELATIVE_ERROR_GUARD,
            "IAPWS liquid-water external quality regression: max relative error={:e} exceeds reviewed guard={MAX_EXTERNAL_RELATIVE_ERROR_GUARD:e}",
            report.max_external_relative_error(),
        );
        for (route, rms) in [
            ("independent I1/I2", report.independent_external_rms()),
            ("canonical TPD", report.canonical_external_rms()),
        ] {
            assert!(
                rms <= MAX_EXTERNAL_RMS_ERROR_GUARD,
                "IAPWS liquid-water external quality regression: {route} RMS relative error={rms:e} exceeds reviewed guard={MAX_EXTERNAL_RMS_ERROR_GUARD:e}",
            );
        }
    }

    #[test]
    #[ignore = "release I5 IAPWS liquid-water characterization with conservative software-regression envelope"]
    fn i5_iapws_water_liquid_low_pressure_boundary_diagnostic() {
        // 275..320 K keeps p_sat near 0.7..10.5 kPa, where the ideal-gas
        // approximation f ~= p is materially better than near 1 atm. It does
        // not trespass below the triple point, where ice is the stable phase.
        let directory = iapws_directory();
        let metadata_path = directory.join("water_liquid_saturation_low_pressure.metadata.json");
        let rows_path = directory.join("water_liquid_saturation_low_pressure.rows.json");
        let metadata_before = fs::read(&metadata_path).unwrap();
        let rows_before = fs::read(&rows_path).unwrap();
        let dataset = iapws_dataset();
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("IAPWS water/liquid diagnostic must resolve local NASA data without NIST");
        let gas = water_gas_scenario();
        let inventory = RealPurePhaseInventory::from_gas(gas.clone(), 0.0)
            .expect("absent liquid boundary inventory must validate");

        let contract = IapwsWaterSaturationComparisonContract::characterization_only();
        contract.validate_dataset(&dataset).unwrap();
        let mut rows = Vec::with_capacity(dataset.rows().len());
        let mut symmetry_rows = Vec::with_capacity(dataset.rows().len());
        for row in dataset.rows() {
            let independent = independent_boundary_pressure(&fixture, &gas, row.temperature_k)
                .unwrap_or_else(|error| {
                    panic!(
                        "{} failed to locate independent boundary at {} K: {error}",
                        dataset.provenance_context(rows.len()),
                        row.temperature_k
                    )
                });
            let canonical =
                canonical_tpd_boundary_pressure(&fixture, &inventory, &gas, row.temperature_k)
                    .unwrap_or_else(|error| {
                        panic!(
                            "{} failed to locate canonical TPD boundary at {} K: {error}",
                            dataset.provenance_context(rows.len()),
                            row.temperature_k
                        )
                    });
            let gas_from_liquid = condensed_side_tpd_boundary_pressure(&fixture, row.temperature_k)
                .unwrap_or_else(|error| {
                    panic!(
                        "{} failed to locate condensed-side TPD boundary at {} K: {error}",
                        dataset.provenance_context(rows.len()),
                        row.temperature_k
                    )
                });
            rows.push(WaterSaturationComparisonRow {
                temperature_k: row.temperature_k,
                iapws_pressure_pa: row.pressure_pa,
                independent,
                canonical_tpd: canonical,
            });
            symmetry_rows.push(WaterSaturationSymmetryRow {
                temperature_k: row.temperature_k,
                independent,
                liquid_from_gas: canonical,
                gas_from_liquid,
            });
        }
        let report = WaterSaturationComparisonReport::new(&dataset, rows).unwrap();
        let symmetry = WaterSaturationSymmetryReport::new(symmetry_rows).unwrap();
        contract.validate_report(&report).unwrap();
        assert_eq!(report.dataset_id(), dataset.metadata().dataset_id);
        assert_eq!(report.rows().len(), dataset.rows().len());
        report
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap_or_else(|error| panic!("{}: {error}", dataset.provenance_context(0)));
        symmetry
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap_or_else(|error| panic!("{}: {error}", dataset.provenance_context(0)));
        assert_iapws_external_regression_envelope(&report);
        println!("{report}\n{symmetry}");
        assert_eq!(fs::read(metadata_path).unwrap(), metadata_before);
        assert_eq!(fs::read(rows_path).unwrap(), rows_before);
    }

    #[test]
    fn i5_direct_water_boundary_is_invariant_to_zero_stoichiometry_carrier() {
        let direct_fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("direct water fixture must resolve locally");
        let carrier_fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &["O2"])
            .expect("carrier water fixture must resolve locally");
        let direct_gas = water_gas_scenario();
        let carrier_gas = carrier_water_gas_scenario();

        let direct = independent_boundary_pressure(&direct_fixture, &direct_gas, 300.0)
            .expect("direct H2O boundary must be evaluable");
        let with_carrier = independent_boundary_pressure(&carrier_fixture, &carrier_gas, 300.0)
            .expect("carrier H2O boundary must be evaluable");
        let relative_delta = (with_carrier.partial_water_pressure_pa
            - direct.partial_water_pressure_pa)
            / direct.partial_water_pressure_pa;

        assert!(
            relative_delta.abs() <= 1e-10,
            "zero-stoichiometry carrier changed p_H2O boundary: direct={} Pa, carrier={} Pa, relative delta={relative_delta:e}",
            direct.partial_water_pressure_pa,
            with_carrier.partial_water_pressure_pa,
        );
    }

    #[test]
    fn direct_water_phase_control_accepts_boundary_pressure_search() {
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("direct water fixture must resolve locally");
        let gas = water_gas_scenario();
        let inventory = RealPurePhaseInventory::from_gas(gas, 0.0)
            .expect("direct absent-liquid inventory must validate");

        let high_pressure = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(275.0, ROOT_UPPER_PRESSURE_PA),
                fixture.initial_composition(&inventory).unwrap(),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .unwrap_or_else(|error| panic!("direct water high-pressure solve failed: {error:#?}"));
        let gas_id = PhaseId::new(Some("gas".to_owned()));
        let liquid_id = PhaseId::new(Some("liquid".to_owned()));
        assert_eq!(
            high_pressure.phase_status(&gas_id),
            Some(PhaseStatus::Inactive)
        );
        assert!(matches!(
            high_pressure.phase_status(&liquid_id),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        ));
        let phase_control = high_pressure
            .phase_control_report()
            .expect("bounded solve must publish lifecycle evidence");
        assert!(phase_control.transitions.iter().any(|transition| {
            transition
                .deactivated
                .iter()
                .any(|phase| phase.index() == 0)
                && matches!(
                    transition.reason,
                    PhaseTransitionReason::BoundaryUnstableActivePhase {
                        minimum_tpd,
                        ..
                    } if minimum_tpd > 0.0
                )
        }));
        assert!((high_pressure.component_moles().iter().sum::<f64>() - 1.0).abs() <= 1e-10);
        assert!(
            condensed_reference_gas_tpd_at_pressure(&fixture, 275.0, ROOT_UPPER_PRESSURE_PA)
                .expect("high-pressure gas TPD must be available")
                > 0.0
        );

        let canonical =
            canonical_tpd_boundary_pressure(&fixture, &inventory, &water_gas_scenario(), 275.0)
                .unwrap_or_else(|error| panic!("direct water boundary search failed: {error:#?}"));
        let independent = independent_boundary_pressure(&fixture, &water_gas_scenario(), 275.0)
            .expect("independent direct-water boundary must be available");
        let relative_delta = (canonical.partial_water_pressure_pa
            - independent.partial_water_pressure_pa)
            / independent.partial_water_pressure_pa;
        assert!(
            relative_delta.abs() <= MAX_INTERNAL_ROOT_RELATIVE_DELTA,
            "direct canonical and independent roots differ: canonical={} Pa, independent={} Pa, relative delta={relative_delta:e}",
            canonical.partial_water_pressure_pa,
            independent.partial_water_pressure_pa,
        );
    }

    #[test]
    fn direct_water_three_way_boundary_roots_agree() {
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("direct water fixture must resolve locally");
        let gas = water_gas_scenario();
        let inventory = RealPurePhaseInventory::from_gas(gas.clone(), 0.0)
            .expect("direct absent-liquid inventory must validate");
        let independent = independent_boundary_pressure(&fixture, &gas, 300.0)
            .expect("independent water boundary must be available");
        let liquid_from_gas = canonical_tpd_boundary_pressure(&fixture, &inventory, &gas, 300.0)
            .expect("gas-side liquid TPD boundary must be available");
        let gas_from_liquid = condensed_side_tpd_boundary_pressure(&fixture, 300.0)
            .expect("liquid-side gas TPD boundary must be available");
        let report = WaterSaturationSymmetryReport::new(vec![WaterSaturationSymmetryRow {
            temperature_k: 300.0,
            independent,
            liquid_from_gas,
            gas_from_liquid,
        }])
        .unwrap();

        report
            .validate_internal_roots(MAX_INTERNAL_ROOT_RELATIVE_DELTA)
            .unwrap_or_else(|error| panic!("direct water three-way boundary mismatch: {error}"));
    }

    #[test]
    fn direct_water_gas_reference_keeps_gas_fixed_and_liquid_evaluated() {
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("direct water fixture must resolve locally");
        let gas = water_gas_scenario();
        let inventory = RealPurePhaseInventory::from_gas(gas, 0.0)
            .expect("direct absent-liquid inventory must validate");
        let solution = solve_resolved_pt(
            ResolvedPhaseEquilibriumRequest::new(
                fixture.resolved(),
                conditions(275.0, ROOT_LOWER_PRESSURE_PA),
                fixture.initial_composition(&inventory).unwrap(),
            )
            .with_phase_control_policy(PhaseControlPolicy::default()),
        )
        .expect("low-pressure gas reference must solve");
        let gas_id = PhaseId::new(Some("gas".to_owned()));
        let liquid_id = PhaseId::new(Some("liquid".to_owned()));
        let acceptance = solution
            .acceptance_report()
            .expect("bounded solve must publish acceptance evidence");
        let gas_stability = &acceptance.phase_stability[phase_index(&solution, &gas_id).unwrap()];
        let liquid_stability =
            &acceptance.phase_stability[phase_index(&solution, &liquid_id).unwrap()];

        assert_eq!(
            gas_stability.status,
            PhaseStabilityStatus::FixedGasAssemblage
        );
        assert_eq!(gas_stability.minimum_tpd, None);
        assert_eq!(liquid_stability.status, PhaseStabilityStatus::Evaluated);
        assert!(
            liquid_stability
                .minimum_tpd
                .is_some_and(|value| value > 0.0)
        );
    }

    #[test]
    fn direct_water_phase_control_evaporates_from_liquid_only_reference() {
        let fixture = RealPurePhaseFamily::WaterLiquid
            .resolve_offline(local_repository(), &[])
            .expect("direct water fixture must resolve locally");
        let low_pressure = solve_from_liquid_reference(&fixture, 275.0, ROOT_LOWER_PRESSURE_PA)
            .unwrap_or_else(|error| {
                panic!("direct water low-pressure evaporation failed: {error:#?}")
            });
        let gas_id = PhaseId::new(Some("gas".to_owned()));
        let liquid_id = PhaseId::new(Some("liquid".to_owned()));
        assert!(matches!(
            low_pressure.phase_status(&gas_id),
            Some(PhaseStatus::Active | PhaseStatus::Appeared)
        ));
        assert_eq!(
            low_pressure.phase_status(&liquid_id),
            Some(PhaseStatus::Inactive)
        );
        let gas_tpd = gas_tpd_from_liquid_reference(&low_pressure)
            .expect("evaporation transition must retain gas TPD evidence");
        assert!(gas_tpd < 0.0, "evaporating gas TPD must be negative");
        assert!(
            low_pressure
                .phase_control_report()
                .expect("bounded solve must publish lifecycle evidence")
                .transitions
                .iter()
                .any(|transition| transition.activated.iter().any(|phase| phase.index() == 0))
        );
        assert!((low_pressure.component_moles().iter().sum::<f64>() - 1.0).abs() <= 1e-10);
    }
}
