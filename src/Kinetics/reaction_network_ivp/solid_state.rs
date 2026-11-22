use crate::Kinetics::solid_state_kinetics_IVP::KineticModel;

use super::error::NetworkIvpError;

pub(crate) fn validate_model(
    reaction_id: &str,
    model: &KineticModel,
    alpha: f64,
) -> Result<(), NetworkIvpError> {
    let valid_parameters = match model {
        KineticModel::SB(m, n, p) => [m, n, p].iter().all(|v| v.is_finite() && **v >= 0.0),
        KineticModel::JMA(m) | KineticModel::Ac(m) => m.is_finite() && *m > 0.0,
        KineticModel::Dec(m) => m.is_finite() && *m >= 0.0,
        KineticModel::PTe(m, n) => [m, n].iter().all(|v| v.is_finite() && **v >= 0.0),
        KineticModel::SBtp(m, n, c) => [m, n, c].iter().all(|v| v.is_finite() && **v >= 0.0),
        _ => true,
    };
    if !valid_parameters {
        return Err(NetworkIvpError::SolidStateInput {
            reaction_id: reaction_id.to_owned(),
            field: "model parameters",
            reason: "must be finite and non-negative (JMA/Ac exponent strictly positive)",
        });
    }
    validate_conversion(reaction_id, model, alpha)
}

pub(crate) fn validate_conversion(
    reaction_id: &str,
    model: &KineticModel,
    alpha: f64,
) -> Result<(), NetworkIvpError> {
    // An endpoint needs a usable state Jacobian and must not drive alpha past 1.
    let domain = match model {
        KineticModel::F1 | KineticModel::F2 | KineticModel::F3 => (true, true),
        KineticModel::R2 | KineticModel::R3 => (true, false),
        _ => (false, false),
    };
    let inside = alpha.is_finite()
        && (alpha > 0.0 || (domain.0 && alpha == 0.0))
        && (alpha < 1.0 || (domain.1 && alpha == 1.0));
    if !inside {
        return Err(domain_error(
            reaction_id,
            model,
            alpha,
            "outside the solver-admissible conversion domain",
        ));
    }
    let rate = model.get_expr().eval_expression(&["a"], &[alpha]);
    if !rate.is_finite() || rate < 0.0 {
        return Err(domain_error(
            reaction_id,
            model,
            alpha,
            "model rate must be finite and non-negative",
        ));
    }
    Ok(())
}

fn domain_error(
    reaction_id: &str,
    model: &KineticModel,
    alpha: f64,
    reason: &'static str,
) -> NetworkIvpError {
    NetworkIvpError::SolidStateDomain {
        reaction_id: reaction_id.to_owned(),
        model: format!("{model:?}"),
        alpha,
        reason,
    }
}
