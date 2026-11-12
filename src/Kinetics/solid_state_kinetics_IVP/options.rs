use RustedSciThe::numerical::SolverTelemetryMode;

/// Numerical controls shared by the scalar solid-state IVP adapter.
///
/// `max_steps` limits integration attempts where the selected native solver
/// exposes that control. `newton_tolerance` and `max_newton_iterations` apply
/// specifically to Backward Euler.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidStateIvpSolverOptions {
    pub max_step: f64,
    pub rtol: f64,
    pub atol: f64,
    pub first_step: Option<f64>,
    pub max_steps: Option<usize>,
    pub newton_tolerance: f64,
    pub max_newton_iterations: usize,
    pub telemetry: SolverTelemetryMode,
}

impl Default for SolidStateIvpSolverOptions {
    fn default() -> Self {
        Self {
            max_step: 0.1,
            rtol: 1e-3,
            atol: 1e-3,
            first_step: None,
            max_steps: None,
            newton_tolerance: 1e-3,
            max_newton_iterations: 100_000,
            telemetry: SolverTelemetryMode::Off,
        }
    }
}

impl SolidStateIvpSolverOptions {
    pub(crate) fn validate(&self) -> Result<(), String> {
        for (field, value) in [
            ("max_step", self.max_step),
            ("rtol", self.rtol),
            ("atol", self.atol),
            ("newton_tolerance", self.newton_tolerance),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(format!("{field} must be finite and positive"));
            }
        }
        if let Some(first_step) = self.first_step {
            if !first_step.is_finite() || first_step <= 0.0 || first_step > self.max_step {
                return Err(
                    "first_step must be finite, positive, and no larger than max_step".into(),
                );
            }
        }
        if self.max_steps == Some(0) {
            return Err("max_steps must be positive when specified".into());
        }
        if self.max_newton_iterations == 0 {
            return Err("max_newton_iterations must be positive".into());
        }
        Ok(())
    }
}
