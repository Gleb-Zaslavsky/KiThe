//! Strict scalar conversions at the reactor task-document boundary.
//!
//! Task documents can contain compatibility scalar variants produced by older
//! parsers and GUI versions. The conversion rules live here so every caller
//! agrees on what is lossless. In particular, floating-point counters are
//! accepted only when they are finite, non-negative, integral, and representable
//! as `usize` without changing their value.

use RustedSciThe::command_interpreter::task_parser::Value;
use std::fmt;

/// Explains why a task-document scalar cannot represent a `usize` exactly.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum UsizeConversionError {
    NegativeInteger(i64),
    OutOfRangeInteger(i64),
    NegativeFloat(f64),
    NonFiniteFloat(f64),
    FractionalFloat(f64),
    OutOfRangeFloat(f64),
    InvalidString(String),
    MissingOptionalValue,
    UnsupportedType(&'static str),
}

impl fmt::Display for UsizeConversionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NegativeInteger(value) => {
                write!(f, "negative integer `{value}` cannot represent usize")
            }
            Self::OutOfRangeInteger(value) => {
                write!(f, "integer `{value}` is outside the usize range")
            }
            Self::NegativeFloat(value) => {
                write!(f, "negative float `{value}` cannot represent usize")
            }
            Self::NonFiniteFloat(value) => {
                write!(f, "non-finite float `{value}` cannot represent usize")
            }
            Self::FractionalFloat(value) => {
                write!(
                    f,
                    "fractional float `{value}` cannot represent usize without truncation"
                )
            }
            Self::OutOfRangeFloat(value) => {
                write!(f, "float `{value}` is outside the exact usize range")
            }
            Self::InvalidString(value) => {
                write!(f, "string `{value}` is not a valid usize")
            }
            Self::MissingOptionalValue => write!(f, "optional integer is None"),
            Self::UnsupportedType(kind) => {
                write!(f, "value type `{kind}` cannot represent usize")
            }
        }
    }
}

/// Converts all supported compatibility scalar shapes without lossy casts.
pub(crate) fn try_value_as_usize(value: &Value) -> Result<usize, UsizeConversionError> {
    match value {
        Value::Usize(value) => Ok(*value),
        Value::Integer(value) => {
            if *value < 0 {
                return Err(UsizeConversionError::NegativeInteger(*value));
            }
            usize::try_from(*value).map_err(|_| UsizeConversionError::OutOfRangeInteger(*value))
        }
        Value::Float(value) => {
            if !value.is_finite() {
                return Err(UsizeConversionError::NonFiniteFloat(*value));
            }
            if *value < 0.0 {
                return Err(UsizeConversionError::NegativeFloat(*value));
            }
            if value.fract() != 0.0 {
                return Err(UsizeConversionError::FractionalFloat(*value));
            }

            // Use an exclusive power-of-two boundary instead of
            // `usize::MAX as f64`: on 64-bit targets that cast rounds up to
            // 2^64 and would otherwise admit one out-of-range value.
            let upper_exclusive = (usize::MAX as u128 + 1) as f64;
            if *value >= upper_exclusive {
                return Err(UsizeConversionError::OutOfRangeFloat(*value));
            }

            let converted = *value as usize;
            if converted as f64 != *value {
                return Err(UsizeConversionError::OutOfRangeFloat(*value));
            }
            Ok(converted)
        }
        Value::String(value) => value
            .trim()
            .parse::<usize>()
            .map_err(|_| UsizeConversionError::InvalidString(value.clone())),
        Value::Optional(Some(inner)) => try_value_as_usize(inner),
        Value::Optional(None) => Err(UsizeConversionError::MissingOptionalValue),
        Value::Vector(_) => Err(UsizeConversionError::UnsupportedType("Vector")),
        Value::Boolean(_) => Err(UsizeConversionError::UnsupportedType("Boolean")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_only_lossless_usize_compatibility_shapes() {
        assert_eq!(try_value_as_usize(&Value::Usize(0)), Ok(0));
        assert_eq!(try_value_as_usize(&Value::Integer(42)), Ok(42));
        assert_eq!(try_value_as_usize(&Value::Float(42.0)), Ok(42));
        assert_eq!(
            try_value_as_usize(&Value::String(" 42 ".to_string())),
            Ok(42)
        );
        assert_eq!(
            try_value_as_usize(&Value::Optional(Some(Box::new(Value::Usize(42))))),
            Ok(42)
        );
        assert_eq!(
            try_value_as_usize(&Value::Usize(usize::MAX)),
            Ok(usize::MAX)
        );
    }

    #[test]
    fn rejects_lossy_or_invalid_usize_shapes() {
        assert!(matches!(
            try_value_as_usize(&Value::Integer(-1)),
            Err(UsizeConversionError::NegativeInteger(-1))
        ));
        assert!(matches!(
            try_value_as_usize(&Value::Float(-1.0)),
            Err(UsizeConversionError::NegativeFloat(value)) if value == -1.0
        ));
        assert!(matches!(
            try_value_as_usize(&Value::Float(1.9)),
            Err(UsizeConversionError::FractionalFloat(value)) if value == 1.9
        ));
        assert!(matches!(
            try_value_as_usize(&Value::Float(f64::NAN)),
            Err(UsizeConversionError::NonFiniteFloat(value)) if value.is_nan()
        ));
        assert!(matches!(
            try_value_as_usize(&Value::Float(f64::INFINITY)),
            Err(UsizeConversionError::NonFiniteFloat(value)) if value.is_infinite()
        ));
        let float_upper_exclusive = (usize::MAX as u128 + 1) as f64;
        assert!(matches!(
            try_value_as_usize(&Value::Float(float_upper_exclusive)),
            Err(UsizeConversionError::OutOfRangeFloat(value))
                if value == float_upper_exclusive
        ));
        assert!(matches!(
            try_value_as_usize(&Value::String("184467440737095516160".to_string())),
            Err(UsizeConversionError::InvalidString(_))
        ));
        assert!(matches!(
            try_value_as_usize(&Value::Optional(None)),
            Err(UsizeConversionError::MissingOptionalValue)
        ));
    }
}
