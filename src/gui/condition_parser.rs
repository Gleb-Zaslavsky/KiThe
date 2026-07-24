//! Typed parsing helpers for GUI numeric condition fields.
//!
//! The GUI accepts a small set of user-editable physical conditions:
//! temperature, pressure, and temperature ranges.
//! These helpers reject malformed, non-finite, and non-positive values before the
//! calculation layer sees them.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConditionField {
    Temperature,
    Pressure,
    T0,
    Tend,
}

impl ConditionField {
    fn display_name(self) -> &'static str {
        match self {
            Self::Temperature => "Temperature",
            Self::Pressure => "Pressure",
            Self::T0 => "T0",
            Self::Tend => "Tend",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParsedRange {
    pub t0: f64,
    pub tend: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConditionParseError {
    InvalidNumber {
        field: ConditionField,
        input: String,
    },
    NotFinite {
        field: ConditionField,
        value: f64,
    },
    NotPositive {
        field: ConditionField,
        value: f64,
    },
    InvalidRange {
        t0: f64,
        tend: f64,
    },
}

impl std::fmt::Display for ConditionParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidNumber { field, .. } => {
                write!(f, "{} must be a valid number", field.display_name())
            }
            Self::NotFinite { field, .. } => write!(f, "{} must be finite", field.display_name()),
            Self::NotPositive { field, .. } => {
                write!(f, "{} must be positive", field.display_name())
            }
            Self::InvalidRange { .. } => {
                write!(f, "T0 must be smaller than Tend")
            }
        }
    }
}

impl std::error::Error for ConditionParseError {}

fn parse_positive_finite(raw: &str, field: ConditionField) -> Result<f64, ConditionParseError> {
    let trimmed = raw.trim();
    let value = trimmed
        .parse::<f64>()
        .map_err(|_| ConditionParseError::InvalidNumber {
            field,
            input: raw.to_string(),
        })?;

    if !value.is_finite() {
        return Err(ConditionParseError::NotFinite { field, value });
    }

    if value <= 0.0 {
        return Err(ConditionParseError::NotPositive { field, value });
    }

    Ok(value)
}

pub fn parse_temperature(raw: &str) -> Result<f64, ConditionParseError> {
    parse_positive_finite(raw, ConditionField::Temperature)
}

pub fn parse_pressure(raw: &str) -> Result<f64, ConditionParseError> {
    parse_positive_finite(raw, ConditionField::Pressure)
}

pub fn parse_temperature_range(
    t0_raw: &str,
    tend_raw: &str,
) -> Result<ParsedRange, ConditionParseError> {
    let t0 = parse_positive_finite(t0_raw, ConditionField::T0)?;
    let tend = parse_positive_finite(tend_raw, ConditionField::Tend)?;

    if t0 >= tend {
        return Err(ConditionParseError::InvalidRange { t0, tend });
    }

    Ok(ParsedRange { t0, tend })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_scalar_conditions_accept_valid_positive_values() {
        assert_eq!(parse_temperature("300.0").unwrap(), 300.0);
        assert_eq!(parse_pressure("101325.0").unwrap(), 101325.0);
    }

    #[test]
    fn parse_scalar_conditions_reject_bad_inputs() {
        let cases = [
            (
                "Temperature",
                parse_temperature("abc").unwrap_err().to_string(),
                "Temperature must be a valid number",
            ),
            (
                "Temperature",
                parse_temperature("0").unwrap_err().to_string(),
                "Temperature must be positive",
            ),
            (
                "Temperature",
                parse_temperature("-1").unwrap_err().to_string(),
                "Temperature must be positive",
            ),
            (
                "Temperature",
                parse_temperature("NaN").unwrap_err().to_string(),
                "Temperature must be finite",
            ),
            (
                "Temperature",
                parse_temperature("inf").unwrap_err().to_string(),
                "Temperature must be finite",
            ),
            (
                "Pressure",
                parse_pressure("abc").unwrap_err().to_string(),
                "Pressure must be a valid number",
            ),
            (
                "Pressure",
                parse_pressure("0").unwrap_err().to_string(),
                "Pressure must be positive",
            ),
            (
                "Pressure",
                parse_pressure("-1").unwrap_err().to_string(),
                "Pressure must be positive",
            ),
            (
                "Pressure",
                parse_pressure("NaN").unwrap_err().to_string(),
                "Pressure must be finite",
            ),
            (
                "Pressure",
                parse_pressure("-inf").unwrap_err().to_string(),
                "Pressure must be finite",
            ),
        ];

        for (field, actual, expected) in cases {
            assert_eq!(actual, expected, "field {field}");
        }
    }

    #[test]
    fn parse_temperature_range_rejects_bad_inputs() {
        assert_eq!(
            parse_temperature_range("abc", "350.0")
                .unwrap_err()
                .to_string(),
            "T0 must be a valid number"
        );
        assert_eq!(
            parse_temperature_range("300.0", "abc")
                .unwrap_err()
                .to_string(),
            "Tend must be a valid number"
        );
        assert_eq!(
            parse_temperature_range("0", "350.0")
                .unwrap_err()
                .to_string(),
            "T0 must be positive"
        );
    }

    #[test]
    fn parse_temperature_range_rejects_reversed_range() {
        let err = parse_temperature_range("400.0", "350.0").unwrap_err();
        assert_eq!(err.to_string(), "T0 must be smaller than Tend");
    }

    #[test]
    fn parse_temperature_range_accepts_valid_bounds() {
        let range = parse_temperature_range("300.0", "350.0").unwrap();
        assert_eq!(range.t0, 300.0);
        assert_eq!(range.tend, 350.0);
    }
}
