use crate::Kinetics::solid_state_kinetics_IVP::KineticModel;
use RustedSciThe::symbolic::symbolic_engine::Expr;

/// A concentration exponent in an empirical rate law. It is not a
/// stoichiometric coefficient and never changes the species balance directly.
#[derive(Debug, Clone)]
pub struct KineticOrder {
    pub species_id: String,
    pub exponent: f64,
}

impl KineticOrder {
    pub fn new(species_id: impl Into<String>, exponent: f64) -> Self {
        Self {
            species_id: species_id.into(),
            exponent,
        }
    }
}

/// A rate constant in units appropriate to the reaction order. For an
/// Arrhenius coefficient, activation energy is J/mol and T is kelvin.
#[derive(Debug, Clone, PartialEq)]
pub enum RateCoefficient {
    Constant(f64),
    Arrhenius {
        prefactor: f64,
        activation_energy: f64,
    },
}

/// Complete numeric coefficient replacement for one prepared reaction.
/// The coefficient variant must match the prepared rate law.
#[derive(Debug, Clone, PartialEq)]
pub struct RateCoefficientBinding {
    pub reaction_id: String,
    pub coefficient: RateCoefficient,
}

impl RateCoefficientBinding {
    pub fn new(reaction_id: impl Into<String>, coefficient: RateCoefficient) -> Self {
        Self {
            reaction_id: reaction_id.into(),
            coefficient,
        }
    }
}

/// The forward rate is mol reaction/(m^3*s). `CustomExpr` is already the
/// whole rate; SolidState multiplies its conversion rate by explicit capacity.
/// Custom expressions use `c0`, `c1`, ... in species order, `t`, and `T`.
#[derive(Debug, Clone)]
pub enum RateLaw {
    MassAction {
        coefficient: RateCoefficient,
    },
    EmpiricalMassAction {
        coefficient: RateCoefficient,
        orders: Vec<KineticOrder>,
    },
    CustomExpr(Expr),
    /// Independent conversion extent: alpha' = k(T) f(alpha),
    /// reaction rate = capacity * alpha'.
    SolidState {
        coefficient: RateCoefficient,
        model: KineticModel,
        capacity: f64,
        initial_conversion: f64,
    },
}

impl RateCoefficient {
    pub(crate) fn parameterized(
        &self,
        reaction_index: usize,
        temperature: &Expr,
    ) -> (Expr, Vec<String>, Vec<f64>) {
        match self {
            Self::Constant(value) => {
                let name = format!("p_rate_{reaction_index}_k");
                (Expr::Var(name.clone()), vec![name], vec![*value])
            }
            Self::Arrhenius {
                prefactor,
                activation_energy,
            } => {
                let gas_constant = 8.314_462_618_153_24;
                let prefactor_name = format!("p_rate_{reaction_index}_a");
                let energy_name = format!("p_rate_{reaction_index}_e");
                let expression = Expr::Var(prefactor_name.clone())
                    * ((Expr::Const(-1.0 / gas_constant) * Expr::Var(energy_name.clone()))
                        / temperature.clone())
                    .exp();
                (
                    expression,
                    vec![prefactor_name, energy_name],
                    vec![*prefactor, *activation_energy],
                )
            }
        }
    }

    pub(crate) fn parameter_count(&self) -> usize {
        match self {
            Self::Constant(_) => 1,
            Self::Arrhenius { .. } => 2,
        }
    }

    pub(crate) fn append_parameter_values(&self, output: &mut Vec<f64>) {
        match self {
            Self::Constant(value) => output.push(*value),
            Self::Arrhenius {
                prefactor,
                activation_energy,
            } => {
                output.push(*prefactor);
                output.push(*activation_energy);
            }
        }
    }

    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Self::Constant(_) => "constant",
            Self::Arrhenius { .. } => "Arrhenius",
        }
    }
}
