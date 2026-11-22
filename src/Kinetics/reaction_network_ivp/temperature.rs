use RustedSciThe::symbolic::symbolic_engine::Expr;

use super::error::{NetworkIvpError, invalid};
use super::problem::validate_expr;

#[derive(Debug, Clone)]
pub enum TemperatureProgram {
    Constant(f64),
    Linear {
        initial_kelvin: f64,
        kelvin_per_second: f64,
    },
    Symbolic(Expr),
}

/// Numeric values that may be rebound without changing the symbolic form of
/// a prepared temperature program. For a linear program, T(t)=T0+beta*t.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TemperatureBinding {
    Constant(f64),
    Linear {
        initial_kelvin: f64,
        kelvin_per_second: f64,
    },
}

impl TemperatureBinding {
    pub(crate) fn parameter_values(self) -> ([f64; 2], usize) {
        match self {
            Self::Constant(kelvin) => ([kelvin, 0.0], 1),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => ([initial_kelvin, kelvin_per_second], 2),
        }
    }

    pub(crate) fn to_program(self) -> TemperatureProgram {
        match self {
            Self::Constant(kelvin) => TemperatureProgram::Constant(kelvin),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => TemperatureProgram::Linear {
                initial_kelvin,
                kelvin_per_second,
            },
        }
    }

    pub(crate) fn kind(self) -> &'static str {
        match self {
            Self::Constant(_) => "constant",
            Self::Linear { .. } => "linear",
        }
    }
}

impl TemperatureProgram {
    pub(crate) fn initial_binding(&self) -> Option<TemperatureBinding> {
        match self {
            Self::Constant(kelvin) => Some(TemperatureBinding::Constant(*kelvin)),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => Some(TemperatureBinding::Linear {
                initial_kelvin: *initial_kelvin,
                kelvin_per_second: *kelvin_per_second,
            }),
            Self::Symbolic(_) => None,
        }
    }

    pub(crate) fn prepared_expression(&self) -> (Expr, Vec<String>, Vec<f64>) {
        match self {
            Self::Constant(kelvin) => (
                Expr::Var("p_temperature".into()),
                vec!["p_temperature".into()],
                vec![*kelvin],
            ),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => (
                Expr::Var("p_temperature".into())
                    + Expr::Var("p_heating_rate".into()) * Expr::Var("t".into()),
                vec!["p_temperature".into(), "p_heating_rate".into()],
                vec![*initial_kelvin, *kelvin_per_second],
            ),
            Self::Symbolic(expr) => (expr.clone(), Vec::new(), Vec::new()),
        }
    }

    pub(crate) fn validate(&self, start: f64, end: f64) -> Result<(), NetworkIvpError> {
        match self {
            Self::Constant(value) => check_temperature(start, *value),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => {
                if !kelvin_per_second.is_finite() {
                    return Err(invalid("temperature.kelvin_per_second", "must be finite"));
                }
                check_temperature(start, *initial_kelvin + *kelvin_per_second * start)?;
                check_temperature(end, *initial_kelvin + *kelvin_per_second * end)
            }
            Self::Symbolic(expr) => {
                validate_expr(expr, "temperature.expression", &mut |symbol| {
                    (symbol == "t").then_some(()).ok_or_else(|| {
                        NetworkIvpError::TemperatureSymbol {
                            symbol: symbol.to_owned(),
                        }
                    })
                })?;
                // Endpoint checks do not establish positivity between samples.
                check_temperature(start, self.at(start))?;
                check_temperature(end, self.at(end))
            }
        }
    }

    pub fn expression(&self) -> Expr {
        match self {
            Self::Constant(value) => Expr::Const(*value),
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => {
                Expr::Const(*initial_kelvin)
                    + Expr::Const(*kelvin_per_second) * Expr::Var("t".into())
            }
            Self::Symbolic(expr) => expr.clone(),
        }
    }

    pub fn at(&self, time: f64) -> f64 {
        match self {
            Self::Constant(value) => *value,
            Self::Linear {
                initial_kelvin,
                kelvin_per_second,
            } => *initial_kelvin + *kelvin_per_second * time,
            Self::Symbolic(expr) => expr.eval_expression(&["t"], &[time]),
        }
    }
}

pub(crate) fn check_temperature(time: f64, value: f64) -> Result<(), NetworkIvpError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(NetworkIvpError::TemperatureDomain { time, value })
    }
}
