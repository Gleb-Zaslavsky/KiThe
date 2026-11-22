use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;

use super::error::{NetworkIvpError, invalid};
use super::problem::NetworkProblem;
use super::rate_law::RateLaw;
use super::solid_state::validate_conversion;
use super::temperature::check_temperature;

/// Ordered symbolic IVP passed to RST. Concentrations `c0..cN` come first;
/// SolidState conversion coordinates `alpha0..` follow in reaction order.
#[derive(Debug, Clone)]
pub struct AssembledNetwork {
    pub equations: Vec<Expr>,
    pub variables: Vec<String>,
    pub parameter_names: Vec<String>,
    pub parameter_values: Vec<f64>,
    /// Reaction-by-species net coefficients, in input reaction order.
    pub stoichiometry: Vec<Vec<f64>>,
}

impl NetworkProblem {
    pub fn assemble(&self) -> AssembledNetwork {
        let mut variables: Vec<String> =
            (0..self.species().len()).map(|i| format!("c{i}")).collect();
        variables.extend((0..self.solid_state_count()).map(|i| format!("alpha{i}")));
        let species_index: HashMap<&str, usize> = self
            .species()
            .iter()
            .enumerate()
            .map(|(index, species)| (species.id.as_str(), index))
            .collect();
        let (temperature, mut parameter_names, mut parameter_values) =
            self.temperature().prepared_expression();
        let mut equations = vec![Expr::Const(0.0); variables.len()];
        let mut stoichiometry = Vec::with_capacity(self.reactions().len());
        let mut conversion_index = 0;
        for (reaction_index, reaction) in self.reactions().iter().enumerate() {
            let mut row = vec![0.0; self.species().len()];
            for term in &reaction.reactants {
                let index = species_index[term.species_id.as_str()];
                row[index] -= term.coefficient;
            }
            for term in &reaction.products {
                let index = species_index[term.species_id.as_str()];
                row[index] += term.coefficient;
            }
            let rate = match &reaction.rate_law {
                RateLaw::MassAction { coefficient } => multiply_orders(
                    {
                        let (expression, names, values) =
                            coefficient.parameterized(reaction_index, &temperature);
                        parameter_names.extend(names);
                        parameter_values.extend(values);
                        expression
                    },
                    reaction
                        .reactants
                        .iter()
                        .map(|term| (term.species_id.as_str(), term.coefficient)),
                    &species_index,
                    &variables,
                ),
                RateLaw::EmpiricalMassAction {
                    coefficient,
                    orders,
                } => multiply_orders(
                    {
                        let (expression, names, values) =
                            coefficient.parameterized(reaction_index, &temperature);
                        parameter_names.extend(names);
                        parameter_values.extend(values);
                        expression
                    },
                    orders
                        .iter()
                        .map(|order| (order.species_id.as_str(), order.exponent)),
                    &species_index,
                    &variables,
                ),
                RateLaw::CustomExpr(expr) => expr.substitute_variable("T", &temperature),
                RateLaw::SolidState {
                    coefficient,
                    model,
                    capacity,
                    ..
                } => {
                    let (coefficient_expr, names, values) =
                        coefficient.parameterized(reaction_index, &temperature);
                    parameter_names.extend(names);
                    parameter_values.extend(values);
                    let alpha_name = &variables[self.species().len() + conversion_index];
                    let conversion_rhs = coefficient_expr
                        * model
                            .get_expr()
                            .substitute_variable("a", &Expr::Var(alpha_name.clone()));
                    equations[self.species().len() + conversion_index] = conversion_rhs.clone();
                    conversion_index += 1;
                    Expr::Const(*capacity) * conversion_rhs
                }
            };
            for (index, coefficient) in row.iter().copied().enumerate() {
                if coefficient != 0.0 {
                    equations[index] =
                        equations[index].clone() + Expr::Const(coefficient) * rate.clone();
                }
            }
            stoichiometry.push(row);
        }
        AssembledNetwork {
            equations,
            variables,
            parameter_names,
            parameter_values,
            stoichiometry,
        }
    }

    /// Independent RHS probe using the problem's initial conversion values.
    /// The returned vector includes conversion equations after species rows.
    pub fn evaluate_rhs(
        &self,
        time: f64,
        concentrations: &[f64],
    ) -> Result<Vec<f64>, NetworkIvpError> {
        if concentrations.len() != self.species().len() {
            return Err(invalid(
                "concentrations",
                "state length differs from species count",
            ));
        }
        if concentrations
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "concentrations",
                "physical state must be finite and non-negative",
            ));
        }
        let conversions: Vec<f64> = self
            .reactions()
            .iter()
            .filter_map(|reaction| match &reaction.rate_law {
                RateLaw::SolidState {
                    initial_conversion, ..
                } => Some(*initial_conversion),
                _ => None,
            })
            .collect();
        self.evaluate_rhs_with_conversions(time, concentrations, &conversions)
    }

    pub fn evaluate_rhs_with_conversions(
        &self,
        time: f64,
        concentrations: &[f64],
        conversions: &[f64],
    ) -> Result<Vec<f64>, NetworkIvpError> {
        if concentrations.len() != self.species().len()
            || conversions.len() != self.solid_state_count()
        {
            return Err(invalid(
                "state",
                "concentration or conversion count differs from problem",
            ));
        }
        if concentrations
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(invalid(
                "concentrations",
                "physical state must be finite and non-negative",
            ));
        }
        for (reaction, &alpha) in self
            .reactions()
            .iter()
            .filter(|reaction| matches!(reaction.rate_law, RateLaw::SolidState { .. }))
            .zip(conversions)
        {
            if let RateLaw::SolidState { model, .. } = &reaction.rate_law {
                validate_conversion(&reaction.id, model, alpha)?;
            }
        }
        check_temperature(time, self.temperature().at(time))?;
        let assembly = self.assemble();
        let mut names: Vec<&str> = assembly.variables.iter().map(String::as_str).collect();
        names.push("t");
        names.extend(assembly.parameter_names.iter().map(String::as_str));
        let mut values = concentrations.to_vec();
        values.extend_from_slice(conversions);
        values.push(time);
        values.extend_from_slice(&assembly.parameter_values);
        assembly
            .equations
            .iter()
            .enumerate()
            .map(|(component, expr)| {
                let value = expr.eval_expression(&names, &values);
                if value.is_finite() {
                    Ok(value)
                } else {
                    if component < self.species().len() {
                        Err(NetworkIvpError::RateDomain {
                            time,
                            component,
                            species_id: self.species()[component].id.clone(),
                        })
                    } else {
                        let reaction = self
                            .reactions()
                            .iter()
                            .filter(|reaction| {
                                matches!(reaction.rate_law, RateLaw::SolidState { .. })
                            })
                            .nth(component - self.species().len())
                            .unwrap();
                        Err(NetworkIvpError::SolidStateDomain {
                            reaction_id: reaction.id.clone(),
                            model: format!(
                                "{:?}",
                                match &reaction.rate_law {
                                    RateLaw::SolidState { model, .. } => model,
                                    _ => unreachable!(),
                                }
                            ),
                            alpha: conversions[component - self.species().len()],
                            reason: "evaluated conversion rate is nonfinite",
                        })
                    }
                }
            })
            .collect()
    }
}

fn multiply_orders<'a>(
    initial: Expr,
    orders: impl Iterator<Item = (&'a str, f64)>,
    species_index: &HashMap<&str, usize>,
    variables: &[String],
) -> Expr {
    orders.fold(initial, |rate, (species_id, exponent)| {
        if exponent == 0.0 {
            return rate;
        }
        let concentration = Expr::Var(variables[species_index[species_id]].clone());
        let factor = if exponent == 1.0 {
            concentration
        } else {
            concentration.pow(Expr::Const(exponent))
        };
        rate * factor
    })
}
