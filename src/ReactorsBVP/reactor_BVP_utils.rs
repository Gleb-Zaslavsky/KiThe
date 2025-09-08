//! # Reactor BVP Utilities Module
//!
//! This module provides essential utility structures and functions for configuring and solving
//! boundary value problems (BVPs) in chemical reactor modeling. It acts as a configuration
//! layer that simplifies the setup of complex numerical solvers.
//!
//! ## Purpose
//!
//! The module addresses the challenge of configuring BVP solvers with dozens of parameters
//! by providing high-level configuration structs that automatically expand to full parameter
//! maps. Instead of manually specifying tolerances for each variable (C0, C1, J0, J1, etc.),
//! users can set tolerances by variable type ("C", "J", "Teta", "q").
//!
//! ## Main Structures
//!
//! - **`InitialConfig`**: Generates physically meaningful initial guess matrices using templates
//!   - `InitialTemplate` enum: Linear, exponential, sigmoid, and constant profiles
//!   - Automatically maps variable types to appropriate initial conditions
//!
//! - **`ToleranceConfig`**: Manages solver convergence tolerances
//!   - Expands 4 tolerance values to full maps for all variables
//!   - Provides sensible defaults for different variable types
//!
//! - **`BoundsConfig`**: Handles variable bounds for constrained optimization
//!   - Physical bounds (concentrations ∈ [0,1], fluxes unbounded)
//!   - Automatic expansion to all solver variables
//!
//! - **`ScalingConfig`**: Dimensionless scaling parameters
//!   - Temperature and length scales for equation normalization
//!   - Validation and conversion utilities
//!
//! ## Key Functions
//!
//! - `create_tolerance_map()`: Alternative HashMap-based tolerance configuration
//! - `create_bounds_map()`: Alternative HashMap-based bounds configuration
//!
//! ## Interesting Code Features
//!
//! ### Smart Variable Type Detection
//! The `get_variable_type()` method uses string pattern matching to automatically
//! determine variable types from solver unknowns ("C0" → "C", "J1" → "J").
//!
//! ### Template-Based Initial Conditions
//! The `InitialTemplate` enum provides mathematical functions for generating
//! realistic initial profiles. The sigmoid template is particularly clever for
//! modeling sharp transitions in reaction zones.
//!
//! ### Automatic Expansion Pattern
//! All config structs follow the same pattern: store 4 base values, then expand
//! to n×4 full maps. This reduces configuration complexity from O(n) to O(1).
//!
//! ### DMatrix Construction Trick
//! The `generate_initial_guess()` method builds the matrix column-wise by extending
//! a flat vector, then uses `DMatrix::from_vec()` for efficient construction.
//!
//! ### Functional Template Generation
//! Each template type implements its mathematical function inline using match
//! expressions, making the code both readable and efficient.

use super::SimpleReactorBVP::ReactorError;
use nalgebra::DMatrix;
use std::collections::HashMap;

/// Template types for initial guess generation
#[derive(Debug, Clone)]
pub enum InitialTemplate {
    /// Linear interpolation between start and end values
    Linear { start: f64, end: f64 },
    /// Decreasing exponential: start * exp(-decay * z)
    DecreasingExp { start: f64, decay: f64 },
    /// Increasing exponential: end * (1 - exp(-growth * z))
    IncreasingExp { end: f64, growth: f64 },
    /// Constant value throughout domain
    Constant { value: f64 },
    /// Sigmoid transition: start + (end-start) / (1 + exp(-steepness*(z-center)))
    Sigmoid {
        start: f64,
        end: f64,
        steepness: f64,
        center: f64,
    },
}

impl InitialTemplate {
    /// Generate values for n_steps grid points
    pub fn generate(&self, n_steps: usize) -> Vec<f64> {
        let mut values = Vec::with_capacity(n_steps);

        for i in 0..n_steps {
            let z = i as f64 / (n_steps - 1) as f64; // z ∈ [0, 1]

            let value = match self {
                InitialTemplate::Linear { start, end } => start + (end - start) * z,
                InitialTemplate::DecreasingExp { start, decay } => start * (-decay * z).exp(),
                InitialTemplate::IncreasingExp { end, growth } => end * (1.0 - (-growth * z).exp()),
                InitialTemplate::Constant { value } => *value,
                InitialTemplate::Sigmoid {
                    start,
                    end,
                    steepness,
                    center,
                } => start + (end - start) / (1.0 + (-steepness * (z - center)).exp()),
            };

            values.push(value);
        }

        values
    }
}

/// Configuration for generating initial guess matrix
///  struct for creating initial guess. There are 4 types of variables "C", "Teta", "J" and
///  "q" and enum with several types of behaviour like Linear(starting point, end point),
///  DecreaingExp(arametes...), IncreaingExp(parameters..) and some other templates. User must provide
///  pairs "type of parameter -template", unknowns say C1, C2,...J1, J2, q, Teta
#[derive(Debug, Clone)]
pub struct InitialConfig {
    /// Templates for each variable type
    pub templates: HashMap<String, InitialTemplate>,
}

impl InitialConfig {
    pub fn new() -> Self {
        Self {
            templates: HashMap::new(),
        }
    }

    /// Set template for variable type ("C", "J", "Teta", "q")
    pub fn set_template(&mut self, var_type: &str, template: InitialTemplate) {
        self.templates.insert(var_type.to_string(), template);
    }

    /// Generate initial guess matrix
    pub fn generate_initial_guess(
        &self,
        unknowns: &[String],
        n_steps: usize,
    ) -> Result<DMatrix<f64>, ReactorError> {
        let n_vars = unknowns.len();
        let mut data = Vec::with_capacity(n_vars * n_steps);

        for unknown in unknowns {
            let var_type = self.get_variable_type(unknown)?;
            let template = self.templates.get(&var_type).ok_or_else(|| {
                ReactorError::MissingData(format!(
                    "No template found for variable type: {}",
                    var_type
                ))
            })?;

            let values = template.generate(n_steps);
            data.extend(values);
        }

        Ok(DMatrix::from_vec(n_vars, n_steps, data))
    }

    /// Determine variable type from unknown name
    fn get_variable_type(&self, unknown: &str) -> Result<String, ReactorError> {
        if unknown == "Teta" {
            Ok("Teta".to_string())
        } else if unknown == "q" {
            Ok("q".to_string())
        } else if unknown.starts_with('C') {
            Ok("C".to_string())
        } else if unknown.starts_with('J') {
            Ok("J".to_string())
        } else {
            Err(ReactorError::InvalidConfiguration(format!(
                "Unknown variable type for: {}",
                unknown
            )))
        }
    }
}

impl Default for InitialConfig {
    fn default() -> Self {
        let mut config = Self::new();

        // Default templates
        config.set_template(
            "C",
            InitialTemplate::Linear {
                start: 0.5,
                end: 0.1,
            },
        );
        config.set_template("J", InitialTemplate::Constant { value: 0.0 });
        config.set_template(
            "Teta",
            InitialTemplate::Linear {
                start: 0.0,
                end: 1.0,
            },
        );
        config.set_template("q", InitialTemplate::Constant { value: 0.0 });

        config
    }
}
/// Configuration for solver relative tolerances
///
/// Automatically expands to full tolerance map for all concentration (C0,C1,...)
/// and flux (J0,J1,...) variables
#[derive(Debug, Clone)]
pub struct ToleranceConfig {
    /// Relative tolerance for all concentration variables (Ci)
    pub C: f64,
    /// Relative tolerance for all flux variables (Ji)
    pub J: f64,
    /// Relative tolerance for dimensionless temperature (Teta)
    pub Teta: f64,
    /// Relative tolerance for heat flux (q)
    pub q: f64,
}

impl Default for ToleranceConfig {
    fn default() -> Self {
        Self {
            C: 1e-4,
            J: 1e-4,
            Teta: 1e-5,
            q: 1e-4,
        }
    }
}

impl ToleranceConfig {
    /// Create new tolerance configuration
    pub fn new(C: f64, J: f64, Teta: f64, q: f64) -> Self {
        Self { C, J, Teta, q }
    }

    /// Convert to full tolerance map with entries for all variables
    ///
    /// Creates entries: Teta, q, C0, C1, ..., J0, J1, ...
    pub fn to_full_tolerance_map(&self, substances: &[String]) -> HashMap<String, f64> {
        let mut tolerance_map = HashMap::new();

        tolerance_map.insert("Teta".to_string(), self.Teta);
        tolerance_map.insert("q".to_string(), self.q);

        for (i, _) in substances.iter().enumerate() {
            tolerance_map.insert(format!("C{}", i), self.C);
            tolerance_map.insert(format!("J{}", i), self.J);
        }

        tolerance_map
    }
}

/// Create tolerance map from HashMap configuration
///
/// Alternative to ToleranceConfig struct - takes HashMap with keys "C", "J", "Teta", "q"
pub fn create_tolerance_map(
    tolerance_config: HashMap<String, f64>,
    substances: &[String],
) -> HashMap<String, f64> {
    let mut full_tolerance_map = HashMap::new();

    let default_C = tolerance_config.get("C").copied().unwrap_or(1e-4);
    let default_J = tolerance_config.get("J").copied().unwrap_or(1e-4);
    let default_Teta = tolerance_config.get("Teta").copied().unwrap_or(1e-5);
    let default_q = tolerance_config.get("q").copied().unwrap_or(1e-4);

    full_tolerance_map.insert("Teta".to_string(), default_Teta);
    full_tolerance_map.insert("q".to_string(), default_q);

    for (i, _) in substances.iter().enumerate() {
        full_tolerance_map.insert(format!("C{}", i), default_C);
        full_tolerance_map.insert(format!("J{}", i), default_J);
    }

    full_tolerance_map
}

/// Configuration for solver variable bounds
///
/// Automatically expands to full bounds map for all concentration and flux variables
#[derive(Debug, Clone)]
pub struct BoundsConfig {
    /// Bounds for all concentration variables (Ci) - typically (0.0, 1.0)
    pub C: (f64, f64),
    /// Bounds for all flux variables (Ji) - typically (-∞, ∞)
    pub J: (f64, f64),
    /// Bounds for dimensionless temperature (Teta)
    pub Teta: (f64, f64),
    /// Bounds for heat flux (q)
    pub q: (f64, f64),
}

impl Default for BoundsConfig {
    fn default() -> Self {
        Self {
            C: (0.0, 1.0),
            J: (-1e20, 1e20),
            Teta: (-10.0, 10.0),
            q: (-1e20, 1e20),
        }
    }
}

impl BoundsConfig {
    /// Create new bounds configuration
    pub fn new(C: (f64, f64), J: (f64, f64), Teta: (f64, f64), q: (f64, f64)) -> Self {
        Self { C, J, Teta, q }
    }

    /// Convert to full bounds map with entries for all variables
    ///
    /// Creates entries: Teta, q, C0, C1, ..., J0, J1, ...
    pub fn to_full_bounds_map(&self, substances: &[String]) -> HashMap<String, (f64, f64)> {
        let mut bounds_map = HashMap::new();

        bounds_map.insert("Teta".to_string(), self.Teta);
        bounds_map.insert("q".to_string(), self.q);

        for (i, _) in substances.iter().enumerate() {
            bounds_map.insert(format!("C{}", i), self.C);
            bounds_map.insert(format!("J{}", i), self.J);
        }

        bounds_map
    }
}

/// Create bounds map from HashMap configuration
///
/// Alternative to BoundsConfig struct - takes HashMap with keys "C", "J", "Teta", "q"
pub fn create_bounds_map(
    bounds_config: HashMap<String, (f64, f64)>,
    substances: &[String],
) -> HashMap<String, (f64, f64)> {
    let mut full_bounds_map = HashMap::new();

    let default_C = bounds_config.get("C").copied().unwrap_or((0.0, 1.0));
    let default_J = bounds_config.get("J").copied().unwrap_or((-1e20, 1e20));
    let default_Teta = bounds_config.get("Teta").copied().unwrap_or((-10.0, 10.0));
    let default_q = bounds_config.get("q").copied().unwrap_or((-1e20, 1e20));

    full_bounds_map.insert("Teta".to_string(), default_Teta);
    full_bounds_map.insert("q".to_string(), default_q);

    for (i, _) in substances.iter().enumerate() {
        full_bounds_map.insert(format!("C{}", i), default_C);
        full_bounds_map.insert(format!("J{}", i), default_J);
    }

    full_bounds_map
}

/// Configuration for dimensionless scaling parameters
///
/// Defines the characteristic scales used to convert dimensional equations to dimensionless form
#[derive(Debug, Clone)]
pub struct ScalingConfig {
    /// Temperature scaling parameter: Θ = (T - dT)/dT
    pub dT: f64,
    /// Length scaling parameter: z = x/L
    pub L: f64,
    ///
    pub T_scale: f64,
}

impl Default for ScalingConfig {
    fn default() -> Self {
        Self {
            dT: 100.0, // K
            L: 0.01,   // m
            T_scale: 100.0,
        }
    }
}

impl ScalingConfig {
    /// Create new scaling configuration
    pub fn new(dT: f64, L: f64, T_scale: f64) -> Self {
        Self { dT, L, T_scale }
    }

    /// Validate scaling parameters
    pub fn validate(&self) -> Result<(), ReactorError> {
        if self.dT <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Temperature scaling dT must be positive".to_string(),
            ));
        }
        if self.L <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Length scaling L must be positive".to_string(),
            ));
        }
        Ok(())
    }

    /// Convert to HashMap for backward compatibility
    pub fn to_hashmap(&self) -> HashMap<String, f64> {
        let mut map = HashMap::new();
        map.insert("dT".to_string(), self.dT);
        map.insert("L".to_string(), self.L);
        map
    }

    /// Create from HashMap
    pub fn from_hashmap(map: &HashMap<String, f64>) -> Result<Self, ReactorError> {
        let dT = map.get("dT").copied().ok_or_else(|| {
            ReactorError::MissingData("Missing dT in scaling parameters".to_string())
        })?;
        let L = map.get("L").copied().ok_or_else(|| {
            ReactorError::MissingData("Missing L in scaling parameters".to_string())
        })?;
        let T_scale = map.get("T_scale").copied().ok_or_else(|| {
            ReactorError::MissingData("Missing L in scaling parameters".to_string())
        })?;
        let config = Self::new(dT, L, T_scale);
        config.validate()?;
        Ok(config)
    }
}
