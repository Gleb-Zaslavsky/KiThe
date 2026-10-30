pub mod kinetics;
mod mechfinder;

use crate::Kinetics::error::{KineticsError, KineticsResult};
use RustedSciThe::symbolic::symbolic_engine::Expr;
/// The module is equipped with a library of kinetic parameters of chemical reactions obtained as a result of parsing publicly available databases
/// The module takes as input the name of the library and the vector of substances and then produces the following data:
/// 1) all reactions of starting substances with each other, and all their possible products with each other.
/// 2) HashMap with kinetic data of all found reactions
/// ru
/// Модуль снабжен библиотекой кинетических параметров химических реакций, полученной в результате парсинга общедоступныхбаз данных
/// Модуль берет на вход название библиотеки и вектор веществ а затем выдает следующие данные:
/// 1) все реакции исходных веществ между собой, и всех их возможных продуктов между собой.
/// 2) HashMap с кинетическими данными всех найденных реакций
/// ----------------------------------------------------------------
// git config --global pack.windowmemory 80000000
use kinetics::{ElementaryStruct, FalloffStruct, PressureStruct, ThreeBodyStruct};
use log::info;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::collections::{HashMap, HashSet};
use std::f64;
use std::fmt;

fn expr_contains_non_finite_constant(expr: &Expr) -> bool {
    match expr {
        Expr::Var(_) => false,
        Expr::Const(value) => !value.is_finite(),
        Expr::Add(lhs, rhs)
        | Expr::Sub(lhs, rhs)
        | Expr::Mul(lhs, rhs)
        | Expr::Div(lhs, rhs)
        | Expr::Pow(lhs, rhs) => {
            expr_contains_non_finite_constant(lhs) || expr_contains_non_finite_constant(rhs)
        }
        Expr::Exp(value)
        | Expr::Ln(value)
        | Expr::sin(value)
        | Expr::cos(value)
        | Expr::tg(value)
        | Expr::ctg(value)
        | Expr::arcsin(value)
        | Expr::arccos(value)
        | Expr::arctg(value)
        | Expr::arcctg(value) => expr_contains_non_finite_constant(value),
    }
}

/// Maximum number of materialized samples in a temperature sweep.
pub const MAX_TEMPERATURE_SWEEP_POINTS: usize = 100_000;
/// enum for types of chemical kinetics rate contant functions
#[derive(Debug, PartialEq, Serialize, Clone)]
#[serde(rename_all = "lowercase")]
pub enum ReactionType {
    Elem,
    Falloff,
    #[serde(rename = "pres")]
    Pressure,
    #[serde(rename = "three-body")]
    ThreeBody,
    Empirical,
}
impl<'de> Deserialize<'de> for ReactionType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        match s.as_str() {
            "elem" => Ok(ReactionType::Elem),
            "falloff" => Ok(ReactionType::Falloff),
            "pressure" | "pres" => Ok(ReactionType::Pressure),
            "three-body" | "threebody" => Ok(ReactionType::ThreeBody),
            "empirical" => Ok(ReactionType::Empirical),
            _ => Err(serde::de::Error::custom(format!(
                "Unknown reaction type: {}",
                s
            ))),
        }
    }
}
/// struct for reaction data
#[derive(Debug, Deserialize, Serialize, Clone)]
pub struct ReactionData {
    #[serde(rename = "type")]
    pub reaction_type: ReactionType,
    pub eq: String,
    pub react: Option<HashMap<String, f64>>,
    #[serde(flatten)]
    pub data: ReactionKinetics,
}
/// enum for structs of different types of kinetics
#[derive(Debug, Serialize, Deserialize, Clone)]
#[serde(untagged)]
pub enum ReactionKinetics {
    Falloff(FalloffStruct),
    Pressure(PressureStruct),
    ThreeBody(ThreeBodyStruct),
    Elementary(ElementaryStruct),
}
impl ReactionData {
    fn invalid_payload(
        &self,
        field: impl Into<String>,
        reason: impl fmt::Display,
    ) -> KineticsError {
        KineticsError::InvalidReactionPayload {
            reaction: self.eq.clone(),
            field: field.into(),
            reason: reason.to_string(),
        }
    }

    fn invalid_input(&self, field: &'static str, reason: impl fmt::Display) -> KineticsError {
        KineticsError::InvalidKineticInput {
            reaction: self.eq.clone(),
            field,
            reason: reason.to_string(),
        }
    }

    fn validate_arrhenius(&self, values: &[f64], field: &'static str) -> KineticsResult<()> {
        if values.len() != 3 {
            return Err(self.invalid_payload(
                field,
                format!(
                    "expected Arrhenius [A, n, E] with 3 values, got {}",
                    values.len()
                ),
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(self.invalid_payload(field, "Arrhenius values must all be finite"));
        }
        if values[0] < 0.0 {
            return Err(self.invalid_payload(
                field,
                format!(
                    "pre-exponential factor A must be non-negative, got {}",
                    values[0]
                ),
            ));
        }
        Ok(())
    }

    fn validate_efficiencies(
        &self,
        efficiencies: &HashMap<String, f64>,
        field: &'static str,
    ) -> KineticsResult<()> {
        if efficiencies
            .iter()
            .any(|(name, value)| name.trim().is_empty() || !value.is_finite() || *value < 0.0)
        {
            return Err(self.invalid_payload(
                field,
                "efficiencies require nonempty species names and finite non-negative values",
            ));
        }
        Ok(())
    }

    fn validate_payload(&self) -> KineticsResult<()> {
        match &self.data {
            ReactionKinetics::Elementary(data) => {
                self.validate_arrhenius(&data.Arrenius, "Arrenius")?;
            }
            ReactionKinetics::Falloff(data) => {
                self.validate_arrhenius(&data.low_rate, "low_rate")?;
                self.validate_arrhenius(&data.high_rate, "high_rate")?;
                if let Some(efficiencies) = &data.eff {
                    self.validate_efficiencies(efficiencies, "eff")?;
                }
                if let Some(troe) = &data.troe {
                    if !matches!(troe.len(), 3 | 4)
                        || troe.iter().any(|value| !value.is_finite())
                        || !(0.0..=1.0).contains(&troe[0])
                        || troe[1..].iter().any(|value| *value <= 0.0)
                    {
                        return Err(self.invalid_payload(
                            "troe",
                            "Troe parameters must be [A, T3, T1] or [A, T3, T1, T2] with finite values and positive characteristic temperatures",
                        ));
                    }
                }
            }
            ReactionKinetics::Pressure(data) => {
                if data.Arrenius.is_empty() {
                    return Err(self.invalid_payload(
                        "Arrenius",
                        "pressure-dependent model must contain at least one pressure node",
                    ));
                }
                let mut seen_pressures = HashSet::new();
                for (key, parameters) in &data.Arrenius {
                    let pressure = key.parse::<f64>().map_err(|_| {
                        self.invalid_payload("Arrenius", format!("invalid pressure node `{key}`"))
                    })?;
                    if !pressure.is_finite() || pressure <= 0.0 {
                        return Err(self.invalid_payload(
                            "Arrenius",
                            format!("pressure node must be finite and positive, got `{key}`"),
                        ));
                    }
                    if !seen_pressures.insert(pressure.to_bits()) {
                        return Err(self.invalid_payload(
                            "Arrenius",
                            format!("ambiguous duplicate numeric pressure node `{key}`"),
                        ));
                    }
                    self.validate_arrhenius(parameters, "Arrenius node")?;
                }
            }
            ReactionKinetics::ThreeBody(data) => {
                self.validate_arrhenius(&data.Arrenius, "Arrenius")?;
                self.validate_efficiencies(&data.eff, "eff")?;
            }
        }
        Ok(())
    }

    /// Validate the declared model and its serialized numerical payload.
    pub fn validate(&self) -> KineticsResult<()> {
        self.validate_reaction_type()?;
        self.validate_payload()
    }

    fn validate_temperature(&self, temp: f64) -> KineticsResult<()> {
        if !temp.is_finite() || temp <= 0.0 {
            return Err(self.invalid_input(
                "temperature",
                format!("must be finite and positive, got {temp}"),
            ));
        }
        Ok(())
    }

    fn validate_concentrations(&self, concentrations: &HashMap<String, f64>) -> KineticsResult<()> {
        if concentrations
            .iter()
            .any(|(name, value)| name.trim().is_empty() || !value.is_finite() || *value < 0.0)
        {
            return Err(self.invalid_input(
                "concentrations",
                "species names must be nonempty and values finite and non-negative",
            ));
        }
        Ok(())
    }

    fn validate_symbolic_concentrations(
        &self,
        concentrations: &HashMap<String, Expr>,
    ) -> KineticsResult<()> {
        if concentrations
            .iter()
            .any(|(name, value)| name.trim().is_empty() || expr_contains_non_finite_constant(value))
        {
            return Err(self.invalid_input(
                "symbolic concentrations",
                "species names must be nonempty and constant values finite",
            ));
        }
        Ok(())
    }

    fn validate_pressure(&self, pressure: f64) -> KineticsResult<()> {
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(self.invalid_input(
                "pressure",
                format!("must be finite and positive, got {pressure}"),
            ));
        }
        Ok(())
    }

    fn validate_finite_rate(&self, rate: f64) -> KineticsResult<f64> {
        if !rate.is_finite() {
            return Err(KineticsError::NonFiniteRate {
                reaction: self.eq.clone(),
                value: rate.to_string(),
            });
        }
        Ok(rate)
    }

    fn missing_concentration_error(&self) -> KineticsError {
        KineticsError::MissingKineticInput {
            reaction: self.eq.clone(),
            field: "concentrations",
        }
    }

    fn missing_pressure_error(&self) -> KineticsError {
        KineticsError::MissingKineticInput {
            reaction: self.eq.clone(),
            field: "pressure",
        }
    }

    /// Build an elementary reaction without validation for checked constructors.
    ///
    /// Prefer [`ReactionData::try_new_elementary`] for untrusted or user-authored data.
    pub fn new_elementary(
        eq: String,
        arrenius: Vec<f64>,
        react: Option<HashMap<String, f64>>,
    ) -> Self {
        Self {
            reaction_type: ReactionType::Elem,
            eq,
            react: react,
            data: ReactionKinetics::Elementary(ElementaryStruct::new(arrenius)),
        }
    }

    /// Build a falloff reaction without validation for checked constructors.
    ///
    /// Prefer [`ReactionData::try_new_falloff`] for untrusted or user-authored data.
    pub fn new_falloff(
        eq: String,
        low_rate: Vec<f64>,
        high_rate: Vec<f64>,
        eff: Option<HashMap<String, f64>>,
        troe: Option<Vec<f64>>,
        react: Option<HashMap<String, f64>>,
    ) -> Self {
        Self {
            reaction_type: ReactionType::Falloff,
            eq,
            react: react,
            data: ReactionKinetics::Falloff(FalloffStruct::new(low_rate, high_rate, eff, troe)),
        }
    }

    /// Build a three-body reaction without validation for checked constructors.
    ///
    /// Prefer [`ReactionData::try_new_three_body`] for untrusted or user-authored data.
    pub fn new_three_body(
        eq: String,
        arrenius: Vec<f64>,
        eff: HashMap<String, f64>,
        react: Option<HashMap<String, f64>>,
    ) -> Self {
        Self {
            reaction_type: ReactionType::ThreeBody,
            eq,
            react: react,
            data: ReactionKinetics::ThreeBody(ThreeBodyStruct::new(arrenius, eff)),
        }
    }

    /// Build a pressure-dependent reaction without validation for checked constructors.
    ///
    /// Prefer [`ReactionData::try_new_pressure`] for untrusted or user-authored data.
    pub fn new_pressure(
        eq: String,
        arrenius: HashMap<String, Vec<f64>>,
        react: Option<HashMap<String, f64>>,
    ) -> Self {
        Self {
            reaction_type: ReactionType::Pressure,
            eq,
            react: react,
            data: ReactionKinetics::Pressure(PressureStruct::new(arrenius)),
        }
    }

    /// Construct an elementary reaction and validate its payload immediately.
    pub fn try_new_elementary(
        eq: String,
        arrenius: Vec<f64>,
        react: Option<HashMap<String, f64>>,
    ) -> KineticsResult<Self> {
        let reaction = Self::new_elementary(eq, arrenius, react);
        reaction.validate()?;
        Ok(reaction)
    }

    /// Construct a falloff reaction and validate its payload immediately.
    pub fn try_new_falloff(
        eq: String,
        low_rate: Vec<f64>,
        high_rate: Vec<f64>,
        eff: Option<HashMap<String, f64>>,
        troe: Option<Vec<f64>>,
        react: Option<HashMap<String, f64>>,
    ) -> KineticsResult<Self> {
        let reaction = Self::new_falloff(eq, low_rate, high_rate, eff, troe, react);
        reaction.validate()?;
        Ok(reaction)
    }

    /// Construct a three-body reaction and validate its payload immediately.
    pub fn try_new_three_body(
        eq: String,
        arrenius: Vec<f64>,
        eff: HashMap<String, f64>,
        react: Option<HashMap<String, f64>>,
    ) -> KineticsResult<Self> {
        let reaction = Self::new_three_body(eq, arrenius, eff, react);
        reaction.validate()?;
        Ok(reaction)
    }

    /// Construct a pressure-dependent reaction and validate its payload immediately.
    pub fn try_new_pressure(
        eq: String,
        arrenius: HashMap<String, Vec<f64>>,
        react: Option<HashMap<String, f64>>,
    ) -> KineticsResult<Self> {
        let reaction = Self::new_pressure(eq, arrenius, react);
        reaction.validate()?;
        Ok(reaction)
    }

    ///function that checks if the reaction_type field matches the variant of the data field in a ReactionData instance.
    /// Return an error when the declared reaction type does not match its payload.
    pub fn validate_reaction_type(&self) -> KineticsResult<()> {
        match (&self.reaction_type, &self.data) {
            (ReactionType::Elem, ReactionKinetics::Elementary(_)) => Ok(()),
            (ReactionType::Falloff, ReactionKinetics::Falloff(_)) => Ok(()),
            (ReactionType::Pressure, ReactionKinetics::Pressure(_)) => Ok(()),
            (ReactionType::ThreeBody, ReactionKinetics::ThreeBody(_)) => Ok(()),
            (ReactionType::Empirical, _) => Err(KineticsError::UnsupportedReactionModel(
                "empirical symbolic rate expressions are not implemented yet".to_string(),
            )),
            _ => Err(KineticsError::InvalidReactionData(format!(
                "Mismatch between reaction_type ({:?}) and data variant ({:?})",
                self.reaction_type,
                std::mem::discriminant(&self.data)
            ))),
        }
    }
    /// Evaluate the rate coefficient at an absolute temperature in kelvin.
    ///
    /// Negative activation energies and temperature exponents are supported.
    /// A finite result of `0.0` is valid when produced by the model, for example
    /// from zero effective collider concentration or floating-point underflow;
    /// it is never substituted for a non-finite result. Non-finite results are
    /// returned as `KineticsError::NonFiniteRate`.
    /// For falloff reactions, concentrations must use the source mechanism's
    /// convention (commonly mol/cm^3 for CHEMKIN data); the facade performs no
    /// unit conversion, and efficiencies are dimensionless.
    pub fn K_const(
        &self,
        temp: f64,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> KineticsResult<f64> {
        self.validate()?;
        self.validate_temperature(temp)?;
        let rate = match &self.data {
            ReactionKinetics::Elementary(data) => data.K_const(temp),
            ReactionKinetics::Falloff(data) => {
                let concentrations =
                    concentrations.ok_or_else(|| self.missing_concentration_error())?;
                self.validate_concentrations(&concentrations)?;
                data.K_const(temp, concentrations)
            }
            ReactionKinetics::Pressure(data) => {
                let pressure = pres.ok_or_else(|| self.missing_pressure_error())?;
                self.validate_pressure(pressure)?;
                data.K_const(temp, pressure)
            }
            ReactionKinetics::ThreeBody(data) => {
                let concentrations =
                    concentrations.ok_or_else(|| self.missing_concentration_error())?;
                self.validate_concentrations(&concentrations)?;
                data.K_const(temp, concentrations)
            }
        };
        self.validate_finite_rate(rate)
    }
    ///function to calculate reaction rate constant for the range of temperatures: from T0 to Tend, number of points is n
    /// ATTENTION! don't forget to use absolute temperature in Kelvin!
    pub fn K_const_for_T_range(
        &self,
        T0: f64,
        Tend: f64,
        n: usize,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, f64>>,
    ) -> KineticsResult<Vec<f64>> {
        let mut K_const_values = Vec::new();
        let T = Self::create_uniform_vector(T0, Tend, n)?;
        for Ti in T {
            // println!("T, {}", Ti);
            let k = self.K_const(Ti, pres, concentrations.clone())?;
            K_const_values.push(k);
        }
        Ok(K_const_values)
    }
    /// Creates a vector of n uniformly distributed numbers from T0 to Tend (inclusive)
    pub fn create_uniform_vector(T0: f64, Tend: f64, n: usize) -> KineticsResult<Vec<f64>> {
        if n < 2 {
            return Err(KineticsError::InvalidReactionData(
                "n must be at least 2 to include both boundary values".to_string(),
            ));
        }
        if n > MAX_TEMPERATURE_SWEEP_POINTS {
            return Err(KineticsError::ResourceLimitExceeded {
                resource: "temperature sweep",
                requested: n,
                limit: MAX_TEMPERATURE_SWEEP_POINTS,
            });
        }
        if !T0.is_finite() || !Tend.is_finite() || T0 <= 0.0 || Tend <= 0.0 {
            return Err(KineticsError::InvalidReactionData(
                "temperature sweep boundaries must be finite and positive".to_string(),
            ));
        }
        let step = (Tend - T0) / (n - 1) as f64;
        if !step.is_finite() {
            return Err(KineticsError::InvalidReactionData(
                "temperature sweep step must be finite".to_string(),
            ));
        }
        Ok((0..n).map(|i| T0 + i as f64 * step).collect())
    }
    /// generate the symbolic representation of the reaction rate constant
    pub fn K_sym(
        &self,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, Expr>>,
    ) -> KineticsResult<Expr> {
        self.validate()?;
        match &self.data {
            ReactionKinetics::Elementary(data) => Ok(data.K_expr()),
            ReactionKinetics::Falloff(data) => {
                let concentrations =
                    concentrations.ok_or_else(|| self.missing_concentration_error())?;
                self.validate_symbolic_concentrations(&concentrations)?;
                Ok(data.K_expr(concentrations))
            }
            ReactionKinetics::Pressure(data) => {
                let pressure = pres.ok_or_else(|| self.missing_pressure_error())?;
                self.validate_pressure(pressure)?;
                Ok(data.K_expr(pressure))
            }
            ReactionKinetics::ThreeBody(data) => {
                let concentrations =
                    concentrations.ok_or_else(|| self.missing_concentration_error())?;
                self.validate_symbolic_concentrations(&concentrations)?;
                Ok(data.K_expr(concentrations))
            }
        }
    }
    /// it is often prefered to use a scaled temperature, so this function returns the symbolic representation of the reaction rate constant with scaled temperature
    pub fn K_sym_with_scaled_T(
        &self,
        pres: Option<f64>,
        concentrations: Option<HashMap<String, Expr>>,
        scaling_T: Option<Expr>,
    ) -> KineticsResult<Expr> {
        let mut K_sym = self.K_sym(pres, concentrations)?;
        if let Some(scaling_T) = scaling_T {
            K_sym = K_sym.substitute_variable("T", &scaling_T);
        }
        Ok(K_sym)
    }

    /// Convert ReactionData to serde_json::Value
    pub fn to_serde_value(&self) -> KineticsResult<Value> {
        Ok(serde_json::to_value(self)?)
    }
}
///parses each reaction data into a ReactionData struct, validates the reaction type, and stores the parsed data in a hash map (reaction_data_hash) with a unique code (react_code) as the key.
///  It also collects the reaction equations in a vector (equations).
pub fn parse_kinetic_data(
    big_mech: &str,
    vec_of_reactions: &[String],
    vec_of_reaction_value: Vec<Value>,
) -> KineticsResult<(Map<String, Value>, Vec<String>)> {
    if vec_of_reaction_value.len() != vec_of_reactions.len() {
        return Err(KineticsError::LengthMismatch {
            context: "reaction IDs and reaction payloads",
            left: vec_of_reactions.len(),
            right: vec_of_reaction_value.len(),
        });
    }

    let mut reaction_data_hash = Map::new();
    let mut equations = Vec::new();
    info!("______________PARCING REACTION DATA INTO STRUCTS________");
    for (index, (reaction_record, reaction_id)) in vec_of_reaction_value
        .iter()
        .zip(vec_of_reactions)
        .enumerate()
    {
        let react_code = format!("{}_{}", big_mech, reaction_id);
        if reaction_data_hash.contains_key(&react_code) {
            return Err(KineticsError::DuplicateReactionId {
                mechanism: big_mech.to_string(),
                reaction_id: reaction_id.clone(),
            });
        }
        let reactiondata = serde_json::from_value::<ReactionData>(reaction_record.clone())
            .map_err(KineticsError::from)
            .map_err(|error| {
                error.with_context(format!(
                    "failed to deserialize reaction `{reaction_id}` at index {index} in mechanism `{big_mech}`"
                ))
            })?;
        reactiondata.validate().map_err(|error| {
            error.with_context(format!(
                "invalid reaction `{reaction_id}` at index {index} in mechanism `{big_mech}`"
            ))
        })?;
        equations.push(reactiondata.eq.clone());
        let value = serde_json::to_value(&reactiondata).map_err(KineticsError::from).map_err(
            |error| {
                error.with_context(format!(
                    "failed to serialize reaction `{reaction_id}` at index {index} in mechanism `{big_mech}`"
                ))
            },
        )?;
        reaction_data_hash.insert(react_code, value);
    }
    info!("______________PARCING REACTION DATA INTO STRUCTS ENDED________");
    Ok((reaction_data_hash, equations))
}
/// parse Vec of serde Values with reaction data
pub fn parse_kinetic_data_vec(
    vec_of_reaction_value: Vec<Value>,
) -> KineticsResult<(Vec<ReactionData>, Vec<String>)> {
    parse_kinetic_data_vec_inner(vec_of_reaction_value, None)
}

pub(crate) fn parse_kinetic_data_vec_with_context(
    vec_of_reaction_value: Vec<Value>,
    reaction_contexts: &[String],
) -> KineticsResult<(Vec<ReactionData>, Vec<String>)> {
    if reaction_contexts.len() != vec_of_reaction_value.len() {
        return Err(KineticsError::LengthMismatch {
            context: "reaction payloads and diagnostic identities",
            left: vec_of_reaction_value.len(),
            right: reaction_contexts.len(),
        });
    }
    parse_kinetic_data_vec_inner(vec_of_reaction_value, Some(reaction_contexts))
}

fn parse_kinetic_data_vec_inner(
    vec_of_reaction_value: Vec<Value>,
    reaction_contexts: Option<&[String]>,
) -> KineticsResult<(Vec<ReactionData>, Vec<String>)> {
    let mut reaction_dat = Vec::new();
    let mut equations = Vec::new();

    for (index, reaction_record) in vec_of_reaction_value.into_iter().enumerate() {
        let parsed_reaction = (|| {
            let reactiondata = serde_json::from_value::<ReactionData>(reaction_record)
                .map_err(KineticsError::from)
                .map_err(|error| {
                    error.with_context(format!(
                        "failed to deserialize reaction payload at index {index}"
                    ))
                })?;
            reactiondata.validate().map_err(|error| {
                error.with_context(format!(
                    "invalid reaction `{}` at index {index}",
                    reactiondata.eq
                ))
            })?;
            Ok::<_, KineticsError>(reactiondata)
        })();
        let reactiondata = parsed_reaction.map_err(|error| {
            if let Some(context) = reaction_contexts.and_then(|contexts| contexts.get(index)) {
                error.with_context(format!("{context} at batch index {index}"))
            } else {
                error
            }
        })?;
        equations.push(reactiondata.eq.clone());
        reaction_dat.push(reactiondata);
    }
    Ok((reaction_dat, equations))
}
/// struct for chemical mechanism construction
#[allow(non_camel_case_types)] // Preserve the historical public type name.
#[derive(Debug)]
pub struct Mechanism_search {
    pub task_substances: Vec<String>,
    pub task_library: String,
    pub mechanism: Vec<String>,
    pub reactants: Vec<String>,
    pub vec_of_reactions: Vec<String>,
    pub reactdata: Vec<ReactionData>,
}
// set the task to construct mechanism
impl Mechanism_search {
    pub fn new(task_substances: Vec<String>, task_library: String) -> Self {
        Self {
            task_substances,
            task_library,
            mechanism: Vec::new(),
            reactants: Vec::new(),
            vec_of_reactions: Vec::new(),
            reactdata: Vec::new(),
        }
    }

    pub fn default() -> Self {
        Self {
            task_substances: Vec::new(),
            task_library: String::new(),
            mechanism: Vec::new(),
            reactants: Vec::new(),
            vec_of_reactions: Vec::new(),
            reactdata: Vec::new(),
        }
    }
    /// find chemical mechanism using mechfinder API
    pub fn mechfinder_api(&mut self) -> KineticsResult<(Vec<String>, Vec<String>, Vec<String>)> {
        /*
        let tuple = [ "O", "NH3", "NO", "O2", "N2", "N2O", "CO", "C"];
        O,NH3,NO,O2,N2,N2O,CO,C
        let big_mech = "NUIG".to_string();
            let vec: Vec<&str> =  tuple.into_iter().collect();
        */

        let vec: Vec<&str> = self.task_substances.iter().map(|s| s.as_str()).collect();
        let big_mech = self.task_library.clone();
        info!(
            "Р·Р°РґР°РЅРёРµ {:?}, Р±РёР±Р»РёРѕС‚РµРєР° {:?}",
            &big_mech, &vec
        );
        let (mechanism, reactants, vec_of_reactions, vec_of_reaction_value) =
            mechfinder::mechfinder(&big_mech, vec)?;
        // info!("mechanism {:?}", &mechanism);
        let (reactdata, _) = parse_kinetic_data_vec(vec_of_reaction_value)?;
        self.mechanism = mechanism;
        self.reactants = reactants;
        self.vec_of_reactions = vec_of_reactions;
        self.reactdata = reactdata;
        Ok((
            self.mechanism.to_owned(),
            self.reactants.to_owned(),
            self.vec_of_reactions.to_owned(),
        ))
    }
}
//tests
#[cfg(test)]
const ELEM_TESTING_JSON: &str = r#"{"type": "elem",
                 "eq": "NAPH+C2H3<=>NAPHV+C2H4",
                  "Arrenius": [0.408, 4.02, 36822.949]}"#;
#[cfg(test)]
const FALOFF_TESTING_JSON: &str = r#" {"type": "falloff",
                 "eq": "C4H71-3+CH3(+M)<=>C5H10-2(+M)",
                 "low_rate": [3.91e+60, -12.81, 26143.75],
                "high_rate": [100000000000000.0, -0.32, -1097.2009],
                 "eff": {"H2": 2.0, "H2O": 6.0, "CH4": 2.0, "CO": 1.5, "CO2": 2.0, "C2H6": 3.0, "AR": 0.7},
                 "troe": [0.104, 1606.0, 60000.0, 6118.0]} "#;
#[cfg(test)]
const THREE_BODY_TESTING_JSON: &str = r#"{"type": "threebody",
              "eq": "H2+M<=>H+H+M",
              "Arrenius": [4.577e+19, -1.4, 436705.19999999995],
             "eff": {"H2": 2.5, "H2O": 12.0, "CO": 1.9, "CO2": 3.8, "HE": 0.83, "CH4": 2.0, "C2H6": 3.0} }"#;
#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn root_cause(error: &KineticsError) -> &KineticsError {
        match error {
            KineticsError::Context { source, .. } => root_cause(source),
            other => other,
        }
    }

    #[test]
    fn test_mechfinder_api() {
        let mut mech_search = Mechanism_search::new(
            vec!["O".to_string(), "NH3".to_string(), "NO".to_string()],
            "NUIG".to_string(),
        );

        let mech_result = mech_search.mechfinder_api();
        assert!(mech_result.is_ok());
        let (mechanism, reactants, vec_of_reactions) = mech_result.unwrap();

        assert!(!mechanism.is_empty());
        assert!(!reactants.is_empty());
        assert!(!vec_of_reactions.is_empty());
    }

    #[test]
    #[ignore = "recursive mechanism stress story; included in the Kinetics campaign"]
    fn story_cantera_recursive_mechanism_construction_reaches_a_consistent_closure() {
        use crate::Kinetics::User_reactions::KinData;

        let seeds = vec!["H", "H2", "O", "O2", "OH", "H2O", "CO", "CO2", "CH4"];
        let mut kinetics = KinData::new();
        kinetics
            .construct_mechanism(
                seeds.iter().map(|seed| (*seed).to_string()).collect(),
                "Cantera".to_string(),
            )
            .expect("recursive Cantera mechanism construction should complete");

        let mechanism = kinetics.canonical_reaction_ids().unwrap();
        let equations = kinetics.equations();
        let reactions = kinetics.reaction_data().unwrap();
        let discovered_species = kinetics.substances();
        assert!(!mechanism.is_empty());
        assert!(
            seeds
                .iter()
                .all(|seed| discovered_species.iter().any(|species| species == seed))
        );
        assert_eq!(mechanism.len(), equations.len());
        assert_eq!(mechanism.len(), reactions.len());
        assert!(
            reactions
                .iter()
                .zip(equations.iter())
                .all(|(reaction, equation)| reaction.eq == *equation)
        );
        println!(
            "Cantera recursive mechanism closure: reactions={}, species={}, seeds={}",
            mechanism.len(),
            discovered_species.len(),
            seeds.len()
        );
    }

    #[test]
    fn test_default_values() {
        let mech_search = Mechanism_search::default();

        assert!(mech_search.task_substances.is_empty());
        assert!(mech_search.task_library.is_empty());
        assert!(mech_search.mechanism.is_empty());
        assert!(mech_search.reactants.is_empty());
        assert!(mech_search.vec_of_reactions.is_empty());
    }

    #[test]

    fn test_ELEM_parse_kinetic_data() {
        let big_mech: &str = "NUIG";
        let vec_of_reactions = vec!["1".to_string()];
        let reaction = ELEM_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let parsed = parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);
        assert!(parsed.is_ok());
        let (ReactionDataHash, _) = parsed.unwrap();

        assert!(!ReactionDataHash.is_empty());
        //   let elem_saved_to_hash = ReactionDataHash[test_reactions_numbers[0]];
        let key = format!("{}_{}", big_mech, &vec_of_reactions[0]);
        let elem_react_testing_instance: ElementaryStruct =
            serde_json::from_value::<ElementaryStruct>(ReactionDataHash[&key].clone()).unwrap();
        assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }
    #[test]

    fn test_THREEBODY_parse_kinetic_data() {
        let ThreeBodyStruct_test_data: &str = r#"{"Arrenius": [4.577e+19, -1.4, 436705.19999999995],
       "eff": {"H2": 2.5, "H2O": 12.0, "CO": 1.9, "CO2": 3.8, "HE": 0.83, "CH4": 2.0, "C2H6": 3.0} }"#;
        let big_mech: &str = "NUIG";
        let vec_of_reactions = vec!["2".to_string()];
        let reaction = THREE_BODY_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let parsed = parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);
        assert!(parsed.is_ok());
        let (ReactionDataHash, _) = parsed.unwrap();

        assert!(!ReactionDataHash.is_empty());
        let threebody_react_testing_instance: ThreeBodyStruct =
            serde_json::from_value::<ThreeBodyStruct>(
                serde_json::from_str(ThreeBodyStruct_test_data).unwrap(),
            )
            .unwrap();
        let mut Concentrations: HashMap<String, f64> = HashMap::new();
        Concentrations.insert("H".to_string(), 0.5);
        Concentrations.insert("O".to_string(), 0.5);
        assert!(threebody_react_testing_instance.K_const(298.15, Concentrations) > 0.0);
        // assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }
    #[test]
    fn test_THREEBODY_from_lib() {
        use crate::Kinetics::User_reactions::KinData;
        let mut kinetics = KinData::new();
        let C1_react = Some(vec!["C1".to_string()]);
        kinetics.shortcut_reactions = C1_react.clone();

        kinetics
            .get_reactions_from_shortcuts()
            .expect("fixture shortcut should resolve");

        kinetics
            .reactdata_parsing()
            .expect("fixture reaction data should parse");
        assert!(kinetics.vec_of_reaction_data.iter().len() > 0);
    }
    #[test]

    fn test_FALOFF_parse_kinetic_data() {
        let big_mech: &str = "NUIG";
        let vec_of_reactions = vec!["3".to_string()];
        let reaction = FALOFF_TESTING_JSON;
        let vec_of_reaction_value: Vec<Value> = vec![serde_json::from_str(reaction).unwrap()];
        let parsed = parse_kinetic_data(big_mech, &vec_of_reactions, vec_of_reaction_value);
        assert!(parsed.is_ok());
        let (ReactionDataHash, _) = parsed.unwrap();
        assert!(!ReactionDataHash.is_empty());
        let key = format!("{}_{}", big_mech, &vec_of_reactions[0]);
        let falloff_react_testing_instance: FalloffStruct =
            serde_json::from_value::<FalloffStruct>(ReactionDataHash[&key].clone()).unwrap();
        let mut Concentrations: HashMap<String, f64> = HashMap::new();
        Concentrations.insert("H".to_string(), 0.5);
        Concentrations.insert("O".to_string(), 0.5);
        assert!(falloff_react_testing_instance.K_const(298.15, Concentrations) > 0.0);
        // assert!(elem_react_testing_instance.K_const(298.15) > 0.0);
    }

    #[test]
    fn test_three_body_deserialization() {
        // The parser should map the JSON payload to the three-body variant.
        let three_body_json = json!({
            "type": "three-body",
            "eq": "H2+M<=>H+H+M",
            "Arrenius": [4.577e+19, -1.4, 436705.19999999995],
            "eff": {
                "H2": 2.5,
                "H2O": 12.0,
                "CO": 1.9,
                "CO2": 3.8,
                "HE": 0.83,
                "CH4": 2.0,
                "C2H6": 3.0
            }
        });

        let reaction_data: ReactionData = serde_json::from_value(three_body_json).unwrap();
        assert_eq!(
            reaction_data.reaction_type,
            ReactionType::ThreeBody,
            "wrong reaction type!"
        );
        assert!(matches!(reaction_data.data, ReactionKinetics::ThreeBody(_)));
    }

    #[test]
    fn test_parse_kinetic_data_vec_rejects_missing_equation() {
        // Missing `eq` should fail fast and surface a domain error.
        let payload = json!({
            "type": "elem",
            "Arrenius": [1.0, 0.0, 0.0]
        });

        let err = parse_kinetic_data_vec(vec![payload]).unwrap_err();
        assert!(matches!(root_cause(&err), KineticsError::Json(_)));
    }

    #[test]
    fn test_reaction_data_deserialization_rejects_unknown_reaction_type() {
        // Invalid type tags should be rejected before any derived state is built.
        let payload = json!({
            "type": "not-a-real-type",
            "eq": "A=>B",
            "Arrenius": [1.0, 0.0, 0.0]
        });

        let err = serde_json::from_value::<ReactionData>(payload).unwrap_err();
        assert!(err.to_string().contains("Unknown reaction type"));
    }

    #[test]
    fn test_parse_kinetic_data_vec_rejects_missing_rate_parameters() {
        // Keep a malformed reaction record around so missing numeric payloads fail loudly.
        let payload = json!({
            "type": "elem",
            "eq": "A=>B"
        });

        let err = parse_kinetic_data_vec(vec![payload]).unwrap_err();
        assert!(matches!(root_cause(&err), KineticsError::Json(_)));
    }

    #[test]
    fn reaction_payload_validation_runs_for_both_mechanism_parsers() {
        let payload = json!({
            "type": "elem",
            "eq": "A=>B",
            "Arrenius": [1.0, 0.0]
        });

        let decoded = serde_json::from_value::<ReactionData>(payload.clone())
            .expect("serde shape is valid even though the domain payload is incomplete");
        assert!(matches!(
            decoded.validate(),
            Err(KineticsError::InvalidReactionPayload { reaction, field, reason })
                if reaction == "A=>B" && field == "Arrenius" && reason.contains("3 values")
        ));

        let vector_error = parse_kinetic_data_vec(vec![payload.clone()]).unwrap_err();
        assert!(matches!(
            root_cause(&vector_error),
            KineticsError::InvalidReactionPayload { reaction, field, reason }
                if reaction == "A=>B" && field == "Arrenius" && reason.contains("3 values")
        ));
        assert!(matches!(
            &vector_error,
            KineticsError::Context { context, .. } if context.contains("index 0")
        ));

        let mechanism_error =
            parse_kinetic_data("test-mechanism", &["reaction-1".into()], vec![payload])
                .unwrap_err();
        assert!(matches!(
            root_cause(&mechanism_error),
            KineticsError::InvalidReactionPayload { reaction, field, reason }
                if reaction == "A=>B" && field == "Arrenius" && reason.contains("3 values")
        ));
        assert!(matches!(
            &mechanism_error,
            KineticsError::Context { context, .. }
                if context.contains("reaction-1") && context.contains("test-mechanism")
        ));
    }

    #[test]
    fn parser_errors_preserve_serde_source_and_reaction_context() {
        let payload = json!({
            "type": "elem",
            "eq": "A=>B",
            "Arrenius": "not-an-array"
        });
        let err = parse_kinetic_data_vec(vec![payload]).unwrap_err();

        assert!(matches!(
            &err,
            KineticsError::Context { context, .. } if context.contains("index 0")
        ));
        assert!(matches!(root_cause(&err), KineticsError::Json(_)));
        assert!(std::error::Error::source(root_cause(&err)).is_some());
    }

    #[test]
    fn mechanism_parser_rejects_id_payload_length_mismatch() {
        let err = parse_kinetic_data("test", &["reaction-1".into()], Vec::new()).unwrap_err();
        assert!(matches!(
            err,
            KineticsError::LengthMismatch {
                context: "reaction IDs and reaction payloads",
                left: 1,
                right: 0
            }
        ));
    }

    #[test]
    fn try_reaction_constructors_validate_payloads_immediately() {
        assert!(ReactionData::try_new_elementary("A=>B".into(), vec![1.0], None).is_err());
        assert!(
            ReactionData::try_new_elementary("A=>B".into(), vec![1.0, 0.0, f64::NAN], None)
                .is_err()
        );
        assert!(ReactionData::try_new_elementary("A=>B".into(), vec![1.0, 0.0, 0.0], None).is_ok());

        assert!(
            ReactionData::try_new_falloff(
                "A(+M)=>B(+M)".into(),
                vec![1.0, 0.0, 0.0],
                vec![1.0, 0.0, 0.0],
                None,
                Some(vec![0.5, 1000.0]),
                None,
            )
            .is_err()
        );

        assert!(
            ReactionData::try_new_three_body(
                "A+M=>B+M".into(),
                vec![1.0, 0.0, 0.0],
                HashMap::from([("M".into(), -1.0)]),
                None,
            )
            .is_err()
        );

        assert!(
            ReactionData::try_new_pressure(
                "A=>B".into(),
                HashMap::from([
                    ("1".into(), vec![1.0, 0.0, 0.0]),
                    ("1.0".into(), vec![2.0, 0.0, 0.0])
                ]),
                None,
            )
            .is_err()
        );
    }

    #[test]
    fn mechanism_parser_rejects_duplicate_reaction_ids_without_overwriting() {
        let payload = json!({
            "type": "elem",
            "eq": "A=>B",
            "Arrenius": [1.0, 0.0, 0.0]
        });
        let err = parse_kinetic_data(
            "fixture",
            &["R1".into(), "R1".into()],
            vec![payload.clone(), payload],
        )
        .unwrap_err();

        assert!(matches!(
            err,
            KineticsError::DuplicateReactionId { mechanism, reaction_id }
                if mechanism == "fixture" && reaction_id == "R1"
        ));
    }

    #[test]
    fn test_pres_deserialization() {
        // This regression test verifies pressure-dependent JSON deserialization.
        let pres_json = json!({
            "type": "pres",
            "eq": "SC4H9<=>C3H6+CH3",
            "Arrenius": {
                "0.01": [2.89e+40, -9.76, 140552.983],
                "0.1": [1.8e+44, -10.5, 154800.281]
            }
        });

        let reaction_data: ReactionData = serde_json::from_value(pres_json).unwrap();
        assert_eq!(
            reaction_data.reaction_type,
            ReactionType::Pressure,
            "wrong reaction type!"
        );
        assert!(matches!(reaction_data.data, ReactionKinetics::Pressure(_)));
    }
    #[test]
    fn test_pres_data_deserialization() {
        const PRES_TESTING_JSON: &str = r#"{
                        "Arrenius":{"0.01": [2.89e+40, -9.76, 140552.983],
                                    "0.1": [1.8e+44, -10.5, 154800.281]
                                    }}"#;

        let pres_val = serde_json::from_str(PRES_TESTING_JSON).expect("valid pressure JSON");
        let pres_data = serde_json::from_value::<PressureStruct>(pres_val)
            .expect("pressure payload should deserialize");
        assert!(!pres_data.Arrenius.is_empty());
    }
}
