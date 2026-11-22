use RustedSciThe::numerical::StiffAPI::{
    PreparedStiffSolver, StiffSolveOutcome, StiffTelemetrySnapshot,
};
use RustedSciThe::numerical::{CancellationToken, SolverTelemetryMode};
use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;
use std::time::Instant;

use super::diagnostics::{NetworkIvpStage, NetworkIvpStageDiagnostics};
use super::error::{NetworkIvpError, invalid};
use super::problem::NetworkProblem;
use super::rate_law::{RateCoefficient, RateCoefficientBinding, RateLaw};
use super::solid_state::validate_conversion;
use super::solver_backend::{
    NetworkIvpOptions, NetworkIvpPartialTrajectory, NetworkIvpSolution, prepare_native,
};
use super::temperature::{TemperatureBinding, check_temperature};

/// A prepared symbolic network and native solver. Repeated runs may change the
/// initial concentrations and numeric temperature parameters, but not the
/// equations, temperature law shape, or interval.
/// Failed runs never replace the last accepted solution.
pub struct PreparedNetworkIvp<'a> {
    problem: &'a NetworkProblem,
    options: &'a NetworkIvpOptions,
    native: PreparedStiffSolver,
    attempted: bool,
    accepted_temperature: Option<TemperatureBinding>,
    native_temperature: Option<TemperatureBinding>,
    parameter_count: usize,
    temperature_parameter_count: usize,
    rate_parameter_schema: Vec<RateCoefficientBinding>,
    rate_binding_indices: HashMap<String, usize>,
    accepted_rate_bindings: Vec<RateCoefficientBinding>,
    native_bindings_match_accepted: bool,
    last_solution: Option<NetworkIvpSolution>,
    last_attempt_telemetry: Option<StiffTelemetrySnapshot>,
    last_failed_trajectory: Option<NetworkIvpPartialTrajectory>,
    stage_diagnostics: Option<NetworkIvpStageDiagnostics>,
}

impl<'a> PreparedNetworkIvp<'a> {
    pub fn new(
        problem: &'a NetworkProblem,
        options: &'a NetworkIvpOptions,
    ) -> Result<Self, NetworkIvpError> {
        options.validate()?;
        let mut stage_diagnostics = NetworkIvpStageDiagnostics::new(options.telemetry);
        let assembly_started = stage_diagnostics
            .as_ref()
            .filter(|diagnostics| diagnostics.records_timings())
            .map(|_| Instant::now());
        let assembled = problem.assemble();
        if let Some(diagnostics) = &mut stage_diagnostics {
            diagnostics.record(
                NetworkIvpStage::EquationAssembly,
                assembly_started.map(|started| started.elapsed()),
            );
        }
        let rate_parameter_schema = rate_parameter_schema(problem);
        let parameter_count = assembled.parameter_names.len();
        let rate_parameter_count = rate_parameter_schema
            .iter()
            .map(|binding| binding.coefficient.parameter_count())
            .sum::<usize>();
        let rate_binding_indices = rate_parameter_schema
            .iter()
            .enumerate()
            .map(|(index, binding)| (binding.reaction_id.clone(), index))
            .collect();

        let preparation_started = stage_diagnostics
            .as_ref()
            .filter(|diagnostics| diagnostics.records_timings())
            .map(|_| Instant::now());
        let native = prepare_native(problem, assembled, options)?;
        if let Some(diagnostics) = &mut stage_diagnostics {
            diagnostics.record(
                NetworkIvpStage::BackendPreparation,
                preparation_started.map(|started| started.elapsed()),
            );
        }
        Ok(Self {
            problem,
            options,
            native,
            attempted: false,
            accepted_temperature: problem.temperature().initial_binding(),
            native_temperature: problem.temperature().initial_binding(),
            parameter_count,
            temperature_parameter_count: parameter_count - rate_parameter_count,
            rate_parameter_schema: rate_parameter_schema.clone(),
            accepted_rate_bindings: rate_parameter_schema,
            native_bindings_match_accepted: true,
            rate_binding_indices,
            last_solution: None,
            last_attempt_telemetry: None,
            last_failed_trajectory: None,
            stage_diagnostics,
        })
    }

    pub fn solve(&mut self) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        let initial = self.problem.initial_concentrations();
        let temperature = self.accepted_temperature;
        self.solve_request(&initial, temperature, None)
    }

    pub fn solve_with_initial_concentrations(
        &mut self,
        initial: &[f64],
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        let temperature = self.accepted_temperature;
        self.solve_request(initial, temperature, None)
    }

    /// Rebind only numeric values of the already prepared temperature law.
    /// This is a fresh solve, not continuation from the previous endpoint.
    pub fn solve_with_temperature_binding(
        &mut self,
        initial: &[f64],
        binding: TemperatureBinding,
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        self.solve_request(initial, Some(binding), None)
    }

    /// Rebind the complete set of numeric reaction coefficients by reaction
    /// ID. Input order is irrelevant; omitted or incompatible entries fail
    /// before the native solver is touched.
    pub fn solve_with_rate_coefficients(
        &mut self,
        initial: &[f64],
        bindings: &[RateCoefficientBinding],
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        let temperature = self.accepted_temperature;
        self.solve_request(initial, temperature, Some(bindings))
    }

    /// Rebind temperature and all numeric reaction coefficients atomically
    /// for one fresh solve while reusing the prepared native callbacks.
    pub fn solve_with_temperature_and_rate_coefficients(
        &mut self,
        initial: &[f64],
        temperature: TemperatureBinding,
        bindings: &[RateCoefficientBinding],
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        self.solve_request(initial, Some(temperature), Some(bindings))
    }

    fn solve_request(
        &mut self,
        initial: &[f64],
        temperature: Option<TemperatureBinding>,
        bindings: Option<&[RateCoefficientBinding]>,
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        self.last_attempt_telemetry = None;
        self.last_failed_trajectory = None;
        self.validate_initial(initial)?;
        self.validate_binding(temperature)?;
        let (rate_bindings_to_commit, parameter_values) =
            self.parameter_values_for(bindings, temperature)?;
        self.native_temperature = temperature;
        self.native_bindings_match_accepted = temperature == self.accepted_temperature
            && rate_bindings_to_commit
                .as_ref()
                .is_none_or(|rates| *rates == self.accepted_rate_bindings);
        let original_initial = !self.attempted && self.initial_matches_problem(initial);
        let result = if self.parameter_count > 0 {
            let initial = DVector::from_vec(self.problem.initial_state(initial));
            let parameters = DVector::from_vec(parameter_values);
            self.run_integration(|native| {
                native.try_solve_with_parameter_values(initial, parameters)
            })
        } else if original_initial {
            self.run_integration(PreparedStiffSolver::try_solve)
        } else {
            let initial = DVector::from_vec(self.problem.initial_state(initial));
            self.run_integration(|native| native.try_solve_with_initial_state(initial))
        };
        self.attempted = true;
        self.capture_attempt_diagnostics(result.is_err(), true);
        let outcome = result.map_err(NetworkIvpError::Integration)?;
        self.publish_outcome(outcome, initial, temperature, rate_bindings_to_commit)
    }

    pub fn solve_with_cancellation_token(
        &mut self,
        initial: &[f64],
        cancellation: CancellationToken,
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        self.last_attempt_telemetry = None;
        self.last_failed_trajectory = None;
        self.validate_initial(initial)?;
        let cancelled_before_run = cancellation.is_cancelled();
        if self.accepted_temperature != self.native_temperature
            || !self.native_bindings_match_accepted
        {
            return Err(NetworkIvpError::UnsupportedCapability {
                family: self.options.family,
                capability: "cancellation after a failed parameter rebind requires a successful restoration of the accepted bind",
            });
        }
        if self.attempted || !self.initial_matches_problem(initial) {
            let (t0, tf) = self.problem.time_interval();
            self.native
                .try_restart_with_initial_state(
                    t0,
                    DVector::from_vec(self.problem.initial_state(initial)),
                    tf,
                )
                .map_err(NetworkIvpError::Restart)?;
        }
        self.attempted = true;
        let result =
            self.run_integration(|native| native.try_solve_with_cancellation_token(cancellation));
        self.capture_attempt_diagnostics(result.is_err(), !cancelled_before_run);
        let outcome = result.map_err(NetworkIvpError::CancellationAwareIntegration)?;
        self.publish_outcome(outcome, initial, self.accepted_temperature, None)
    }

    fn publish_outcome(
        &mut self,
        outcome: StiffSolveOutcome,
        initial: &[f64],
        binding: Option<TemperatureBinding>,
        rate_bindings: Option<Vec<RateCoefficientBinding>>,
    ) -> Result<&NetworkIvpSolution, NetworkIvpError> {
        let postprocessing_started = self.timer_start();
        let result = make_candidate(self.problem, self.options, &self.native, outcome, binding);
        self.record_stage(
            NetworkIvpStage::OutputPostprocessing,
            postprocessing_started,
        );
        let mut candidate = result?;
        candidate.run_metadata = super::solver_backend::NetworkIvpRunMetadata {
            initial_concentrations: initial.to_vec(),
            time_interval_seconds: self.problem.time_interval(),
            temperature_binding: binding,
            rate_coefficients: rate_bindings
                .as_deref()
                .unwrap_or(&self.accepted_rate_bindings)
                .to_vec(),
        };
        candidate.stage_diagnostics = self.stage_diagnostics.clone();
        self.last_attempt_telemetry = candidate.telemetry.clone();
        self.last_solution = Some(candidate);
        self.accepted_temperature = binding;
        if let Some(rate_bindings) = rate_bindings {
            self.accepted_rate_bindings = rate_bindings;
        }
        self.native_bindings_match_accepted = true;
        Ok(self.last_solution.as_ref().unwrap())
    }

    fn parameter_values_for(
        &self,
        bindings: Option<&[RateCoefficientBinding]>,
        temperature: Option<TemperatureBinding>,
    ) -> Result<(Option<Vec<RateCoefficientBinding>>, Vec<f64>), NetworkIvpError> {
        let mut parameter_values = Vec::with_capacity(self.parameter_count);
        if self.temperature_parameter_count > 0 {
            let Some(binding) = temperature else {
                return Err(NetworkIvpError::UnbindableTemperatureProgram);
            };
            let (values, count) = binding.parameter_values();
            if count != self.temperature_parameter_count {
                return Err(NetworkIvpError::TemperatureBindingMismatch {
                    expected: if self.temperature_parameter_count == 1 {
                        "constant"
                    } else {
                        "linear"
                    },
                    actual: binding.kind(),
                });
            }
            parameter_values.extend_from_slice(&values[..count]);
        }
        let rate_bindings_to_commit = if let Some(bindings) = bindings {
            let normalized = self.validate_rate_bindings(bindings)?;
            for binding in &normalized {
                binding
                    .coefficient
                    .append_parameter_values(&mut parameter_values);
            }
            (normalized != self.accepted_rate_bindings).then_some(normalized)
        } else {
            for binding in &self.accepted_rate_bindings {
                binding
                    .coefficient
                    .append_parameter_values(&mut parameter_values);
            }
            None
        };
        debug_assert_eq!(parameter_values.len(), self.parameter_count);
        Ok((rate_bindings_to_commit, parameter_values))
    }

    fn validate_rate_bindings(
        &self,
        bindings: &[RateCoefficientBinding],
    ) -> Result<Vec<RateCoefficientBinding>, NetworkIvpError> {
        let mut ordered: Vec<Option<&RateCoefficientBinding>> =
            vec![None; self.rate_parameter_schema.len()];
        for binding in bindings {
            let Some(&index) = self.rate_binding_indices.get(&binding.reaction_id) else {
                return Err(NetworkIvpError::UnknownRateBinding {
                    reaction_id: binding.reaction_id.clone(),
                });
            };
            if ordered[index].replace(binding).is_some() {
                return Err(NetworkIvpError::DuplicateRateBinding {
                    reaction_id: binding.reaction_id.clone(),
                });
            }
            let expected = &self.rate_parameter_schema[index].coefficient;
            if expected.kind() != binding.coefficient.kind() {
                return Err(NetworkIvpError::RateBindingKindMismatch {
                    reaction_id: binding.reaction_id.clone(),
                    expected: expected.kind(),
                    actual: binding.coefficient.kind(),
                });
            }
            validate_rate_coefficient(&binding.reaction_id, &binding.coefficient)?;
        }
        if bindings.len() != self.rate_parameter_schema.len() {
            return Err(NetworkIvpError::RateBindingCountMismatch {
                expected: self.rate_parameter_schema.len(),
                actual: bindings.len(),
            });
        }
        Ok(ordered
            .into_iter()
            .map(|binding| {
                binding
                    .expect("complete rate bindings were checked by count and IDs")
                    .clone()
            })
            .collect())
    }

    fn validate_binding(&self, binding: Option<TemperatureBinding>) -> Result<(), NetworkIvpError> {
        match (self.problem.temperature().initial_binding(), binding) {
            (None, None) => Ok(()),
            (None, Some(_)) => Err(NetworkIvpError::UnbindableTemperatureProgram),
            (Some(expected), Some(actual)) if expected.kind() == actual.kind() => {
                let (t0, tf) = self.problem.time_interval();
                actual.to_program().validate(t0, tf)
            }
            (Some(expected), Some(actual)) => Err(NetworkIvpError::TemperatureBindingMismatch {
                expected: expected.kind(),
                actual: actual.kind(),
            }),
            (Some(expected), None) => Err(NetworkIvpError::TemperatureBindingMismatch {
                expected: expected.kind(),
                actual: "none",
            }),
        }
    }

    fn validate_initial(&self, initial: &[f64]) -> Result<(), NetworkIvpError> {
        if initial.len() != self.problem.species().len() {
            return Err(invalid(
                "initial_concentrations",
                "state length differs from species count",
            ));
        }
        if initial
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "initial_concentrations",
                "physical state must be finite and non-negative",
            ));
        }
        self.problem.validate_solid_state_resources(initial)?;
        Ok(())
    }

    fn initial_matches_problem(&self, initial: &[f64]) -> bool {
        initial.len() == self.problem.species().len()
            && initial
                .iter()
                .zip(self.problem.species())
                .all(|(value, species)| *value == species.initial_concentration)
    }

    fn run_integration<T>(&mut self, run: impl FnOnce(&mut PreparedStiffSolver) -> T) -> T {
        let started = self.timer_start();
        let result = run(&mut self.native);
        self.record_stage(NetworkIvpStage::Integration, started);
        result
    }

    fn timer_start(&self) -> Option<Instant> {
        self.stage_diagnostics
            .as_ref()
            .filter(|diagnostics| diagnostics.records_timings())
            .map(|_| Instant::now())
    }

    fn record_stage(&mut self, stage: NetworkIvpStage, started: Option<Instant>) {
        if let Some(diagnostics) = &mut self.stage_diagnostics {
            diagnostics.record(stage, started.map(|started| started.elapsed()));
        }
    }

    pub fn last_solution(&self) -> Option<&NetworkIvpSolution> {
        self.last_solution.as_ref()
    }

    /// Configured numeric binding, initially from the problem and updated only
    /// after an accepted rebind. None means a fixed symbolic law.
    pub fn current_temperature_binding(&self) -> Option<TemperatureBinding> {
        self.accepted_temperature
    }

    /// Complete accepted set of numeric coefficients for rebindable rate
    /// laws, ordered by the original reaction list.
    pub fn current_rate_coefficients(&self) -> &[RateCoefficientBinding] {
        &self.accepted_rate_bindings
    }

    /// Available native diagnostics from the most recent run request. None
    /// may mean Off, pre-validation rejection, or no attributable failure data.
    pub fn last_attempt_telemetry(&self) -> Option<&StiffTelemetrySnapshot> {
        self.last_attempt_telemetry.as_ref()
    }

    /// Bounded raw native samples observed after the last failed integration,
    /// if the trajectory advanced. RST does not yet expose a solve ID to
    /// prove provenance across repeated failures. Never a published solution.
    pub fn last_failed_trajectory(&self) -> Option<&NetworkIvpPartialTrajectory> {
        self.last_failed_trajectory.as_ref()
    }

    /// KiThe-owned cumulative stage data for this prepared session. It is
    /// available after a failed solve as well as after an accepted one.
    pub fn stage_diagnostics(&self) -> Option<&NetworkIvpStageDiagnostics> {
        self.stage_diagnostics.as_ref()
    }

    fn native_snapshot(&self) -> Option<StiffTelemetrySnapshot> {
        if self.options.telemetry == SolverTelemetryMode::Off {
            return None;
        }
        match &self.native {
            PreparedStiffSolver::Radau(prepared) => prepared
                .solver()
                .last_solve_telemetry()
                .cloned()
                .map(StiffTelemetrySnapshot::Radau),
            _ => self.native.telemetry_snapshot().ok(),
        }
    }

    fn capture_attempt_diagnostics(&mut self, failed: bool, permit_partial: bool) {
        if failed {
            self.last_failed_trajectory = if permit_partial {
                self.failed_trajectory_snapshot()
            } else {
                None
            };
            self.last_attempt_telemetry = self
                .last_failed_trajectory
                .as_ref()
                .and_then(|_| self.native_snapshot());
        } else {
            self.last_attempt_telemetry = self.native_snapshot();
        }
    }

    fn failed_trajectory_snapshot(&self) -> Option<NetworkIvpPartialTrajectory> {
        if self.options.telemetry == SolverTelemetryMode::Off {
            return None;
        }
        let view = self.native.trajectory_view().ok()?;
        let (start, end) = self.problem.time_interval();
        if view.axis.len() < 2
            || view.values.nrows() != view.axis.len()
            || view.values.ncols()
                != self.problem.species().len() + self.problem.solid_state_count()
            || (view.axis[0] - start).abs() > 1e-9 * (end - start).abs().max(1.0)
            || view.axis[view.axis.len() - 1] >= end
            || view.axis[view.axis.len() - 1] <= view.axis[0]
        {
            return None;
        }
        let first = view.axis.len().saturating_sub(32);
        Some(NetworkIvpPartialTrajectory {
            time_seconds: view.axis.iter().skip(first).copied().collect(),
            states: (first..view.axis.len())
                .map(|row| view.values.row(row).iter().copied().collect())
                .collect(),
            omitted_samples: first,
        })
    }
}

fn rate_parameter_schema(problem: &NetworkProblem) -> Vec<RateCoefficientBinding> {
    problem
        .reactions()
        .iter()
        .filter_map(|reaction| {
            let coefficient = match &reaction.rate_law {
                RateLaw::MassAction { coefficient }
                | RateLaw::EmpiricalMassAction { coefficient, .. }
                | RateLaw::SolidState { coefficient, .. } => coefficient.clone(),
                RateLaw::CustomExpr(_) => return None,
            };
            Some(RateCoefficientBinding::new(
                reaction.id.clone(),
                coefficient,
            ))
        })
        .collect()
}

fn validate_rate_coefficient(
    reaction_id: &str,
    coefficient: &RateCoefficient,
) -> Result<(), NetworkIvpError> {
    let valid = match coefficient {
        RateCoefficient::Constant(value) => value.is_finite() && *value >= 0.0,
        RateCoefficient::Arrhenius {
            prefactor,
            activation_energy,
        } => prefactor.is_finite() && *prefactor >= 0.0 && activation_energy.is_finite(),
    };
    if !valid {
        return Err(invalid(
            format!("rate_bindings[{reaction_id}]"),
            "rate parameters must be finite and prefactor non-negative",
        ));
    }
    Ok(())
}

/// Prepare and solve one physical network without retaining its native solver.
pub fn solve_network(
    problem: &NetworkProblem,
    options: &NetworkIvpOptions,
) -> Result<NetworkIvpSolution, NetworkIvpError> {
    let mut prepared = PreparedNetworkIvp::new(problem, options)?;
    prepared.solve()?;
    Ok(prepared.last_solution.take().unwrap())
}

fn make_candidate(
    problem: &NetworkProblem,
    options: &NetworkIvpOptions,
    native: &PreparedStiffSolver,
    outcome: StiffSolveOutcome,
    binding: Option<TemperatureBinding>,
) -> Result<NetworkIvpSolution, NetworkIvpError> {
    match outcome {
        StiffSolveOutcome::Radau(solution) => {
            let trajectory = solution
                .output_trajectory()
                .map_err(NetworkIvpError::RadauExtraction)?;
            let telemetry = (options.telemetry != SolverTelemetryMode::Off)
                .then(|| StiffTelemetrySnapshot::from_radau_solution(&solution));
            validate_candidate(problem, binding, &trajectory.0, &trajectory.1, telemetry)
        }
        _ => {
            let view = native
                .trajectory_view()
                .map_err(NetworkIvpError::Postprocessing)?;
            let candidate = validate_candidate(problem, binding, view.axis, view.values, None)?;
            let telemetry = if options.telemetry == SolverTelemetryMode::Off {
                None
            } else {
                Some(
                    native
                        .telemetry_snapshot()
                        .map_err(NetworkIvpError::Telemetry)?,
                )
            };
            Ok(NetworkIvpSolution {
                telemetry,
                ..candidate
            })
        }
    }
}

pub(super) fn validate_candidate(
    problem: &NetworkProblem,
    binding: Option<TemperatureBinding>,
    axis: &DVector<f64>,
    values: &DMatrix<f64>,
    telemetry: Option<StiffTelemetrySnapshot>,
) -> Result<NetworkIvpSolution, NetworkIvpError> {
    if axis.is_empty()
        || values.nrows() != axis.len()
        || values.ncols() != problem.species().len() + problem.solid_state_count()
    {
        return Err(NetworkIvpError::InvalidOutput {
            reason: "trajectory shape differs from problem",
        });
    }
    let mut time_seconds = Vec::with_capacity(axis.len());
    let mut concentrations = Vec::with_capacity(axis.len());
    let mut conversions = Vec::with_capacity(axis.len());
    let mut temperature_kelvin = Vec::with_capacity(axis.len());
    let rebound_temperature = binding.map(TemperatureBinding::to_program);
    let temperature_program = rebound_temperature
        .as_ref()
        .unwrap_or(problem.temperature());
    for (row_index, time) in axis.iter().copied().enumerate() {
        let previous_time = time_seconds.last().copied();
        if !time.is_finite() || previous_time.is_some_and(|previous| time <= previous) {
            return Err(NetworkIvpError::InvalidTimeSample {
                index: row_index,
                time,
                previous_time,
            });
        }
        let temperature = temperature_program.at(time);
        check_temperature(time, temperature)?;
        for (species_index, value) in values
            .row(row_index)
            .columns(0, problem.species().len())
            .iter()
            .copied()
            .enumerate()
        {
            if !value.is_finite() || value < 0.0 {
                return Err(NetworkIvpError::InvalidConcentration {
                    species_id: problem.species()[species_index].id.clone(),
                    time,
                    value,
                });
            }
        }
        let mut alpha_row = Vec::with_capacity(problem.solid_state_count());
        for (conversion_index, (reaction, value)) in problem
            .reactions()
            .iter()
            .filter(|reaction| matches!(reaction.rate_law, RateLaw::SolidState { .. }))
            .zip(
                values
                    .row(row_index)
                    .columns(problem.species().len(), problem.solid_state_count())
                    .iter()
                    .copied(),
            )
            .enumerate()
        {
            let model = match &reaction.rate_law {
                RateLaw::SolidState { model, .. } => model,
                _ => unreachable!(),
            };
            validate_conversion(&reaction.id, model, value).map_err(|_| {
                NetworkIvpError::InvalidConversion {
                    reaction_id: reaction.id.clone(),
                    time,
                    value,
                }
            })?;
            if conversions
                .last()
                .is_some_and(|previous: &Vec<f64>| value < previous[conversion_index])
            {
                return Err(NetworkIvpError::InvalidConversion {
                    reaction_id: reaction.id.clone(),
                    time,
                    value,
                });
            }
            alpha_row.push(value);
        }
        time_seconds.push(time);
        concentrations.push(
            values
                .row(row_index)
                .columns(0, problem.species().len())
                .iter()
                .copied()
                .collect(),
        );
        conversions.push(alpha_row);
        temperature_kelvin.push(temperature);
    }
    let (t0, tf) = problem.time_interval();
    let coverage_tolerance = 1e-9 * (tf - t0).abs().max(1.0);
    if (time_seconds[0] - t0).abs() > coverage_tolerance
        || (time_seconds.last().unwrap() - tf).abs() > coverage_tolerance
    {
        return Err(NetworkIvpError::InvalidOutput {
            reason: "trajectory does not cover requested time interval",
        });
    }
    Ok(NetworkIvpSolution {
        run_metadata: super::solver_backend::NetworkIvpRunMetadata {
            initial_concentrations: problem.initial_concentrations(),
            time_interval_seconds: (t0, tf),
            temperature_binding: binding,
            rate_coefficients: rate_parameter_schema(problem),
        },
        species_ids: problem.species().iter().map(|s| s.id.clone()).collect(),
        solid_state_reaction_ids: problem
            .reactions()
            .iter()
            .filter(|reaction| matches!(reaction.rate_law, RateLaw::SolidState { .. }))
            .map(|reaction| reaction.id.clone())
            .collect(),
        time_seconds,
        concentrations,
        conversions,
        temperature_kelvin,
        telemetry,
        stage_diagnostics: None,
    })
}
