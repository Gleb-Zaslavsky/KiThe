use std::collections::HashSet;

use RustedSciThe::symbolic::symbolic_engine::Expr;

use super::error::{NetworkIvpError, invalid};
use super::rate_law::{RateCoefficient, RateLaw};
use super::solid_state::validate_model;
use super::temperature::TemperatureProgram;

/// A species label need not be a chemical formula. Its initial concentration
/// is in mol/m^3; the vector order fixes `c0`, `c1`, ... symbolic coordinates.
#[derive(Debug, Clone)]
pub struct SpeciesSpec {
    pub id: String,
    pub initial_concentration: f64,
}

#[derive(Debug, Clone)]
pub struct StoichTerm {
    pub species_id: String,
    pub coefficient: f64,
}

impl StoichTerm {
    pub fn new(species_id: impl Into<String>, coefficient: f64) -> Self {
        Self {
            species_id: species_id.into(),
            coefficient,
        }
    }
}

/// Forward reaction with a complete rate law. Reverse processes are separate
/// reaction stages; no reverse rate is inferred from thermodynamics.
#[derive(Debug, Clone)]
pub struct ReactionSpec {
    pub id: String,
    pub reactants: Vec<StoichTerm>,
    pub products: Vec<StoichTerm>,
    pub rate_law: RateLaw,
}

#[derive(Debug, Clone)]
pub struct NetworkProblem {
    species: Vec<SpeciesSpec>,
    reactions: Vec<ReactionSpec>,
    temperature: TemperatureProgram,
    start_time: f64,
    end_time: f64,
}

impl NetworkProblem {
    pub fn new(
        species: Vec<SpeciesSpec>,
        reactions: Vec<ReactionSpec>,
        temperature: TemperatureProgram,
        start_time: f64,
        end_time: f64,
    ) -> Result<Self, NetworkIvpError> {
        if species.is_empty() {
            return Err(invalid("species", "at least one species is required"));
        }
        if !start_time.is_finite() || !end_time.is_finite() || end_time <= start_time {
            return Err(invalid("time", "finite end time must exceed start time"));
        }
        let mut species_ids = HashSet::new();
        for (index, entry) in species.iter().enumerate() {
            if entry.id.trim().is_empty() {
                return Err(invalid(format!("species[{index}].id"), "must not be blank"));
            }
            if !species_ids.insert(entry.id.as_str()) {
                return Err(NetworkIvpError::DuplicateSpecies {
                    id: entry.id.clone(),
                });
            }
            if !entry.initial_concentration.is_finite() || entry.initial_concentration < 0.0 {
                return Err(invalid(
                    format!("species[{index}].initial_concentration"),
                    "must be finite and non-negative",
                ));
            }
        }
        temperature.validate(start_time, end_time)?;
        let mut reaction_ids = HashSet::new();
        for (index, reaction) in reactions.iter().enumerate() {
            if reaction.id.trim().is_empty() {
                return Err(invalid(
                    format!("reactions[{index}].id"),
                    "must not be blank",
                ));
            }
            if !reaction_ids.insert(reaction.id.as_str()) {
                return Err(NetworkIvpError::DuplicateReaction {
                    id: reaction.id.clone(),
                });
            }
            if reaction.reactants.is_empty() && reaction.products.is_empty() {
                return Err(invalid(
                    format!("reactions[{index}]"),
                    "reaction must contain at least one species",
                ));
            }
            for (side, terms) in [
                (&"reactants", &reaction.reactants),
                (&"products", &reaction.products),
            ] {
                let mut seen = HashSet::new();
                for (term_index, term) in terms.iter().enumerate() {
                    if !species_ids.contains(term.species_id.as_str()) {
                        return Err(NetworkIvpError::UnknownSpecies {
                            reaction_id: reaction.id.clone(),
                            species_id: term.species_id.clone(),
                        });
                    }
                    if !seen.insert(term.species_id.as_str()) {
                        return Err(invalid(
                            format!("reactions[{index}].{side}[{term_index}]"),
                            "duplicate species on one side",
                        ));
                    }
                    if !term.coefficient.is_finite() || term.coefficient <= 0.0 {
                        return Err(invalid(
                            format!("reactions[{index}].{side}[{term_index}].coefficient"),
                            "must be finite and positive",
                        ));
                    }
                }
            }
            match &reaction.rate_law {
                RateLaw::MassAction { coefficient }
                | RateLaw::EmpiricalMassAction { coefficient, .. }
                | RateLaw::SolidState { coefficient, .. } => {
                    let values = match coefficient {
                        RateCoefficient::Constant(value) => vec![*value],
                        RateCoefficient::Arrhenius {
                            prefactor,
                            activation_energy,
                        } => {
                            vec![*prefactor, *activation_energy]
                        }
                    };
                    if values.iter().any(|value| !value.is_finite()) || values[0] < 0.0 {
                        return Err(invalid(
                            format!("reactions[{index}].rate_law"),
                            "rate parameters must be finite and prefactor non-negative",
                        ));
                    }
                    if let RateLaw::EmpiricalMassAction { orders, .. } = &reaction.rate_law {
                        let mut seen = HashSet::new();
                        for (order_index, order) in orders.iter().enumerate() {
                            if !species_ids.contains(order.species_id.as_str()) {
                                return Err(NetworkIvpError::UnknownSpecies {
                                    reaction_id: reaction.id.clone(),
                                    species_id: order.species_id.clone(),
                                });
                            }
                            if !seen.insert(order.species_id.as_str()) {
                                return Err(invalid(
                                    format!("reactions[{index}].rate_law.orders[{order_index}]"),
                                    "duplicate kinetic order for species",
                                ));
                            }
                            if !order.exponent.is_finite() || order.exponent < 0.0 {
                                return Err(invalid(
                                    format!(
                                        "reactions[{index}].rate_law.orders[{order_index}].exponent"
                                    ),
                                    "must be finite and non-negative",
                                ));
                            }
                        }
                    }
                    if let RateLaw::SolidState {
                        model,
                        capacity,
                        initial_conversion,
                        ..
                    } = &reaction.rate_law
                    {
                        if !capacity.is_finite() || *capacity <= 0.0 {
                            return Err(NetworkIvpError::SolidStateInput {
                                reaction_id: reaction.id.clone(),
                                field: "capacity",
                                reason: "must be finite and positive",
                            });
                        }
                        validate_model(&reaction.id, model, *initial_conversion)?;
                    }
                }
                RateLaw::CustomExpr(expr) => validate_expr(
                    expr,
                    &format!("reactions[{index}].rate_law.expression"),
                    &mut |symbol| {
                        let canonical_species = symbol
                            .strip_prefix('c')
                            .and_then(|index| index.parse::<usize>().ok())
                            .filter(|index| {
                                *index < species.len() && symbol == format!("c{index}")
                            });
                        if symbol == "t" || symbol == "T" || canonical_species.is_some() {
                            Ok(())
                        } else {
                            Err(NetworkIvpError::RateSymbol {
                                reaction_id: reaction.id.clone(),
                                symbol: symbol.to_owned(),
                            })
                        }
                    },
                )?,
            }
        }
        let problem = Self {
            species,
            reactions,
            temperature,
            start_time,
            end_time,
        };
        problem.validate_solid_state_resources(&problem.initial_concentrations())?;
        Ok(problem)
    }

    pub fn species(&self) -> &[SpeciesSpec] {
        &self.species
    }
    pub fn reactions(&self) -> &[ReactionSpec] {
        &self.reactions
    }
    pub fn temperature(&self) -> &TemperatureProgram {
        &self.temperature
    }
    pub fn time_interval(&self) -> (f64, f64) {
        (self.start_time, self.end_time)
    }
    pub fn initial_concentrations(&self) -> Vec<f64> {
        self.species
            .iter()
            .map(|entry| entry.initial_concentration)
            .collect()
    }

    pub(crate) fn initial_state(&self, concentrations: &[f64]) -> Vec<f64> {
        let mut state = Vec::with_capacity(concentrations.len() + self.solid_state_count());
        state.extend_from_slice(concentrations);
        state.extend(
            self.reactions
                .iter()
                .filter_map(|reaction| match &reaction.rate_law {
                    RateLaw::SolidState {
                        initial_conversion, ..
                    } => Some(*initial_conversion),
                    _ => None,
                }),
        );
        state
    }

    pub(crate) fn solid_state_count(&self) -> usize {
        self.reactions
            .iter()
            .filter(|reaction| matches!(reaction.rate_law, RateLaw::SolidState { .. }))
            .count()
    }

    pub(crate) fn validate_solid_state_resources(
        &self,
        concentrations: &[f64],
    ) -> Result<(), NetworkIvpError> {
        let mut owners = vec![None::<usize>; self.species.len()];
        for (reaction_index, reaction) in self.reactions.iter().enumerate() {
            for (species_index, species) in self.species.iter().enumerate() {
                let products = reaction
                    .products
                    .iter()
                    .find(|term| term.species_id == species.id)
                    .map_or(0.0, |term| term.coefficient);
                let reactants = reaction
                    .reactants
                    .iter()
                    .find(|term| term.species_id == species.id)
                    .map_or(0.0, |term| term.coefficient);
                let net = products - reactants;
                if net >= 0.0 {
                    continue;
                }
                if let Some(other) = owners[species_index] {
                    if matches!(reaction.rate_law, RateLaw::SolidState { .. })
                        || matches!(self.reactions[other].rate_law, RateLaw::SolidState { .. })
                    {
                        return Err(NetworkIvpError::SharedSolidStateResource {
                            reaction_id: reaction.id.clone(),
                            other_reaction_id: self.reactions[other].id.clone(),
                            species_id: species.id.clone(),
                        });
                    }
                } else {
                    owners[species_index] = Some(reaction_index);
                }
                if let RateLaw::SolidState {
                    capacity,
                    initial_conversion,
                    ..
                } = &reaction.rate_law
                {
                    let required = -net * capacity * (1.0 - initial_conversion);
                    let available = concentrations[species_index];
                    if !required.is_finite() || available < required {
                        return Err(NetworkIvpError::SolidStateResourceShortfall {
                            reaction_id: reaction.id.clone(),
                            species_id: species.id.clone(),
                            required,
                            available,
                        });
                    }
                }
            }
        }
        Ok(())
    }
}

pub(crate) fn validate_expr(
    expr: &Expr,
    field: &str,
    check_symbol: &mut impl FnMut(&str) -> Result<(), NetworkIvpError>,
) -> Result<(), NetworkIvpError> {
    match expr {
        Expr::Var(symbol) => check_symbol(symbol),
        Expr::Const(value) if value.is_finite() => Ok(()),
        Expr::Const(_) => Err(invalid(field, "contains a non-finite constant")),
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) | Expr::Pow(a, b) => {
            validate_expr(a, field, check_symbol)?;
            validate_expr(b, field, check_symbol)
        }
        Expr::Exp(a)
        | Expr::Ln(a)
        | Expr::sin(a)
        | Expr::cos(a)
        | Expr::tg(a)
        | Expr::ctg(a)
        | Expr::arcsin(a)
        | Expr::arccos(a)
        | Expr::arctg(a)
        | Expr::arcctg(a) => validate_expr(a, field, check_symbol),
    }
}
