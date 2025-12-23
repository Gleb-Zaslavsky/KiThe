//! # Pair Diffusion Coefficient Calculation Module
//!
//! ## Overview
//! This module provides comprehensive tools for calculating binary diffusion coefficients
//! between gas pairs using kinetic theory and collision integrals. It supports both
//! individual pair calculations and batch processing of multiple substances.
//!
//! ## Key Features
//! - Binary diffusion coefficient calculations using Chapman-Enskog theory
//! - Support for polar and non-polar molecules with dipole moment corrections
//! - Numerical, symbolic, and closure-based calculations
//! - Taylor series expansions and temperature integration
//! - Multi-substance batch processing with matrix output
//! - Unit conversion support (m²/s, cm²/s)
//!
//! ## Main Structures
//! - [`PairDiffusion`]: Individual pair diffusion calculations
//! - [`MultiSubstanceDiffusion`]: Batch processing for multiple substances
//!
//! ## Usage Example
//! ```rust,ignore
//! // Single pair calculation
//! let mut pair_diff = PairDiffusion::new(input_a, input_b, 300.0, 101325.0, 28.0, 30.0);
//! pair_diff.calc_coefficients();
//! pair_diff.calculate_D();
//!
//! // Multi-substance calculation
//! let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);
//! multi_diff.add_substances_from_library(vec!["CO", "N2", "O2"], "Aramco_transport")?;
//! multi_diff.calculate_all_pairs();
//! let matrix = multi_diff.create_diffusion_matrix();
//! ```

use crate::Thermodynamics::DBhandlers::TRANSPORTdata::{
    TransportData, TransportInput, omega_11_calc, omega_11_calc_sym,
};
use crate::Thermodynamics::DBhandlers::transport_api::TransportError;
use crate::Thermodynamics::thermo_lib_api::ThermoData;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;
use std::fmt;

/// Error types for pair diffusion calculations
#[derive(Debug)]
pub enum PairDiffusionError {
    /// Temperature range is invalid (T_min >= T_max or negative values)
    InvalidTemperatureRange,
    /// Error during symbolic computation or integration
    SymbolicError(String),
    /// Unsupported unit specification
    InvalidUnit(String),
    /// Requested substance not found in library
    SubstanceNotFound(String),
}
impl fmt::Display for PairDiffusionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairDiffusionError::InvalidTemperatureRange => {
                write!(f, "Invalid temperature range specified.")
            }
            PairDiffusionError::SymbolicError(msg) => {
                write!(f, "Symbolic computation error: {}", msg)
            }
            PairDiffusionError::InvalidUnit(unit) => write!(f, "Invalid unit specified: {}", unit),
            PairDiffusionError::SubstanceNotFound(name) => {
                write!(f, "Substance not found: {}", name)
            }
        }
    }
}
/// Binary diffusion coefficient calculator for gas pairs
///
/// This structure handles the calculation of diffusion coefficients between two gas species
/// using kinetic theory. It supports numerical calculations, symbolic expressions, and
/// temperature-dependent closures.
///
/// # Theory
/// Uses Chapman-Enskog kinetic theory with collision integrals (Ω₁₁) to calculate
/// binary diffusion coefficients. Handles both polar and non-polar molecules with
/// appropriate mixing rules for Lennard-Jones parameters.
///
/// # Formula
/// D₁₂ = 1.858×10⁻³ × T^1.5 × √(1/M₁ + 1/M₂) / (P × σ₁₂² × Ω₁₁)
pub struct PairDiffusion {
    /// Transport parameters for substance A (σ, ε/k, μ, α)
    pub substance_A_input: TransportInput,
    /// Transport parameters for substance B (σ, ε/k, μ, α)
    pub substance_B_input: TransportInput,
    /// Temperature range for integration calculations (T_min, T_max)
    pub T_interval: Option<(f64, f64)>,
    /// Current temperature [K]
    pub T: f64,
    /// System pressure [Pa]
    pub P: f64,
    /// Molar mass of substance A [kg/mol]
    pub M_A: f64,
    /// Molar mass of substance B [kg/mol]
    pub M_B: f64,
    /// Mixed well depth parameter ε₁₂/k [K]
    pub e_k: f64,
    /// Mixed collision diameter σ₁₂ [Å]
    pub sigma: f64,
    /// Mixed dipole moment parameter [Debye]
    pub mu: f64,
    /// Calculated diffusion coefficient [m²/s or cm²/s]
    pub D: f64,
    /// Unit specification for diffusion coefficient
    pub D_unit: Option<String>,
    /// Unit conversion multiplier
    pub unit_multiplier: f64,
    /// Temperature-dependent diffusion coefficient function
    pub D_closure: Box<dyn Fn(f64) -> f64 + 'static>,
    /// Symbolic expression for diffusion coefficient D(T)
    pub D_symbolic: Expr,
}
impl PairDiffusion {
    pub fn default() -> Self {
        PairDiffusion {
            substance_A_input: TransportInput::default(),
            substance_B_input: TransportInput::default(),
            T_interval: None,
            T: 300.0,
            P: 101325.0,
            M_A: 0.0,
            M_B: 0.0,
            e_k: 0.0,
            sigma: 0.0,
            mu: 0.0,
            D: 0.0,
            D_unit: None,
            unit_multiplier: 1.0,
            D_closure: Box::new(|_T| 0.0),
            D_symbolic: Expr::Const(0.0),
        }
    }
    pub fn new(
        substance_A_input: TransportInput,
        substance_B_input: TransportInput,
        T: f64,
        P: f64,
        M_A: f64,
        M_B: f64,
    ) -> Self {
        PairDiffusion {
            substance_A_input,
            substance_B_input,
            T_interval: None,
            T,
            P,
            M_A,
            M_B,
            e_k: 0.0,
            sigma: 0.0,
            mu: 0.0,
            D: 0.0,
            D_unit: None,
            unit_multiplier: 1.0,
            D_closure: Box::new(|_T| 0.0),
            D_symbolic: Expr::Const(0.0),
        }
    }

    pub fn set_T_interval(&mut self, T_min: f64, T_max: f64) {
        self.T_interval = Some((T_min, T_max));
    }

    pub fn set_D_unit(&mut self, unit: Option<String>) -> Result<(), PairDiffusionError> {
        if let Some(unit) = unit {
            self.D_unit = Some(unit.clone());
            match unit.as_str() {
                "m2/s" => self.unit_multiplier = 1.0,
                "cm2/s" => self.unit_multiplier = 1e4,
                _ => return Err(PairDiffusionError::InvalidUnit(unit)),
            }
        }
        Ok(())
    }

    pub fn calc_coefficients(&mut self) {
        // e_k = 'well_depth'
        //  sigma_k = 'diam'
        //  mu_A= 'dipole'
        //  alpha = 'polar'
        let TransportInput {
            Altname: _,
            polar: alpha_A,
            Form: _,
            diam: sigma_k_A,
            well_depth: e_k_A,
            rot_relax: _,
            dipole: mu_A,
        } = self.substance_A_input;
        let TransportInput {
            Altname: _,
            polar: alpha_B,
            Form: _,
            diam: sigma_k_B,
            well_depth: e_k_B,
            rot_relax: _,
            dipole: mu_B,
        } = self.substance_B_input;
        const K_B: f64 = 1.38e-23;

        let (e_k, sigma, mu) = if (mu_A == 0.0 && mu_B == 0.0) || (mu_A != 0.0 && mu_B != 0.0) {
            let e_k = (e_k_A * e_k_B).powf(0.5);
            let sigma = (sigma_k_A + sigma_k_B) / 2.0;
            let mu = (mu_A * mu_B).powf(0.5);
            (e_k, sigma, mu)
        } else {
            let (mu, e_k_p, sigma_p, e_k_n, sigma_n, alpha_n) = if mu_A != 0.0 {
                (mu_A, e_k_A, sigma_k_A, e_k_B, sigma_k_B, alpha_B)
            } else {
                (mu_B, e_k_B, sigma_k_B, e_k_A, sigma_k_A, alpha_A)
            };

            let alpha_1 = alpha_n / sigma_n.powi(3);
            let mu_1 = 3.162e-10 * mu / (e_k_p * K_B * sigma_p.powi(3)).powf(0.5);
            let xi = 1.0 + (1.0 / 4.0) * alpha_1 * mu_1 * (e_k_p / e_k_n).powf(0.5);
            let sigma = xi.powf(-1.0 / 6.0) * (sigma_n + sigma_p) / 2.0;
            let e_k = xi.powi(2) * (e_k_n * e_k_p).powf(0.5);
            (e_k, sigma, 0.0)
        };
        self.e_k = e_k;
        self.sigma = sigma;
        self.mu = mu;
    }

    pub fn calculate_D(&mut self) {
        let P = self.P;
        let T = self.T;
        let e_k = self.e_k;
        let sigma = self.sigma;
        let mu = self.mu;
        let M_a = self.M_A;
        let M_b = self.M_B;
        let um = self.unit_multiplier;
        let omega = omega_11_calc(T, e_k, mu, sigma);
        let m1 = 1.0 / M_a + 1.0 / M_b;
        let D = um * 1.858e-3 * T.powf(1.5) * m1.powf(0.5) / (P * sigma.powi(2) * omega);

        self.D = D;
    }

    pub fn calculate_D_closure(&mut self) {
        let P = self.P;

        let e_k = self.e_k.clone();
        let sigma = self.sigma.clone();
        let mu = self.mu.clone();
        let M_a = self.M_A;
        let M_b = self.M_B;
        let um = self.unit_multiplier;
        let m1 = 1.0 / M_a + 1.0 / M_b;
        let D_closure = move |T| {
            let omega = omega_11_calc(T, e_k, mu, sigma);
            um * 1.858e-3 * T.powf(1.5) * m1.powf(0.5) / (P * sigma.powi(2) * omega)
        };
        self.D_closure = Box::new(D_closure);
    }

    pub fn calculate_D_symbolic(&mut self) {
        let P = self.P;

        let e_k = self.e_k.clone();
        let sigma = self.sigma.clone();
        let mu = self.mu.clone();
        let M_a = self.M_A;
        let M_b = self.M_B;
        let um = self.unit_multiplier;
        let m1 = 1.0 / M_a + 1.0 / M_b;
        let omega = omega_11_calc_sym(e_k, mu, sigma);
        let T = Expr::Var("T".to_string());
        let P = Expr::Const(P);
        let k = Expr::Const(um * 1.858e-3);
        let sigma = Expr::Const(sigma);
        let m1 = Expr::Const(m1);
        let D_symbolic = k * T.pow(Expr::Const(1.5)) * m1.pow(Expr::Const(0.5))
            / (P * sigma.pow(Expr::Const(2.0)) * omega);
        self.D_symbolic = D_symbolic;
    }

    pub fn taylor_series(&mut self, T0: f64, n: usize) {
        let D_sym = self.D_symbolic.clone();
        let series = D_sym.taylor_series1D("T", T0, n);
        self.D_symbolic = series;
    }

    pub fn integr_mean(&mut self) -> Result<(), PairDiffusionError> {
        use RustedSciThe::symbolic::symbolic_integration::QuadMethod;
        let (T_min, T_max) = self
            .T_interval
            .ok_or(PairDiffusionError::InvalidTemperatureRange)?;

        if T_min >= T_max {
            return Err(PairDiffusionError::InvalidTemperatureRange);
        }
        let T_range = T_max - T_min;
        let inv_T_range = 1.0 / T_range;
        let D_sym = self.D_symbolic.clone();
        let I = inv_T_range
            * D_sym
                .quad(QuadMethod::GaussLegendre, 3, T_min, T_max, None)
                .map_err(|e| {
                    PairDiffusionError::SymbolicError(format!("Failed to integrate entropy: {}", e))
                })?;
        self.D = I;
        Ok(())
    }
}

impl Clone for PairDiffusion {
    fn clone(&self) -> Self {
        PairDiffusion {
            substance_A_input: self.substance_A_input.clone(),
            substance_B_input: self.substance_B_input.clone(),
            T_interval: self.T_interval,
            T: self.T,
            P: self.P,
            M_A: self.M_A,
            M_B: self.M_B,
            e_k: self.e_k,
            sigma: self.sigma,
            mu: self.mu,
            D: self.D,
            D_unit: self.D_unit.clone(),
            unit_multiplier: self.unit_multiplier,
            D_closure: Box::new(|_T| 0.0),
            D_symbolic: self.D_symbolic.clone(),
        }
    }
}
use std::fmt::Debug;
impl Debug for PairDiffusion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PairDiffusion")
            .field("substance_A_input", &self.substance_A_input)
            .field("substance_B_input", &self.substance_B_input)
            .field("T_interval", &self.T_interval)
            .field("T", &self.T)
            .field("P", &self.P)
            .field("M_A", &self.M_A)
            .field("M_B", &self.M_B)
            .field("e_k", &self.e_k)
            .field("sigma", &self.sigma)
            .field("mu", &self.mu)
            .field("D", &self.D)
            .field("D_unit", &self.D_unit)
            .field("unit_multiplier", &self.unit_multiplier)
            .field("D_symbolic", &self.D_symbolic)
            .finish()
    }
}
#[derive(Debug)]
/// Multi-substance diffusion coefficient calculator
///
/// Handles batch calculation of binary diffusion coefficients for multiple gas species.
/// Automatically loads transport data from libraries and calculates all possible pairs
/// using double-loop iteration. Results can be accessed as individual pairs or as
/// a symmetric diffusion matrix.
///
/// # Features
/// - Batch loading from thermodynamic databases
/// - Automatic calculation of all N×N pairs (with N(N+1)/2 unique calculations)
/// - Matrix representation for numerical algorithms
/// - Common unit and parameter management
/// - Error handling for missing substances or libraries
///
/// # Usage Pattern
/// ```rust,ignore
/// let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);
/// multi_diff.add_substances_from_library(vec!["CO", "N2", "O2"], "Aramco_transport")?;
/// multi_diff.set_molar_masses(masses);
/// multi_diff.calculate_all_pairs();
/// let matrix = multi_diff.create_diffusion_matrix();
/// ```
#[derive(Clone)]
pub struct MultiSubstanceDiffusion {
    /// Storage for individual substance transport data
    pub substances: HashMap<String, TransportData>,
    /// Ordered list of substance names for matrix indexing
    pub substance_names: Vec<String>,
    /// Calculated pair diffusion data, keyed by (substance_A, substance_B)
    pub pair_diffusions: HashMap<(String, String), PairDiffusion>,
    /// Symmetric N×N diffusion coefficient matrix [m²/s or cm²/s]
    pub diffusion_matrix: Option<Vec<Vec<f64>>>,
    /// System temperature [K]
    pub T: f64,
    /// System pressure [Pa]
    pub P: f64,
    /// Common unit for all diffusion coefficients
    pub D_unit: Option<String>,
}

impl MultiSubstanceDiffusion {
    pub fn default() -> Self {
        Self {
            substances: HashMap::new(),
            substance_names: Vec::new(),
            pair_diffusions: HashMap::new(),
            diffusion_matrix: None,
            T: 300.0,
            P: 101325.0,
            D_unit: None,
        }
    }
    /// Creates new multi-substance diffusion calculator
    ///
    /// # Arguments
    /// * `T` - System temperature [K]
    /// * `P` - System pressure [Pa]
    pub fn new(T: f64, P: f64) -> Self {
        Self {
            substances: HashMap::new(),
            substance_names: Vec::new(),
            pair_diffusions: HashMap::new(),
            diffusion_matrix: None,
            T,
            P,
            D_unit: None,
        }
    }

    /// Adds individual substance with transport data
    ///
    /// # Arguments
    /// * `name` - Substance identifier (e.g., "CO", "N2")
    /// * `transport_data` - Pre-configured TransportData instance
    pub fn add_substance(&mut self, name: String, mut transport_data: TransportData) {
        transport_data.P = self.P;
        self.substances.insert(name.clone(), transport_data);
        if !self.substance_names.contains(&name) {
            self.substance_names.push(name);
        }
    }

    /// Batch loads substances from thermodynamic database
    ///
    /// # Arguments
    /// * `substance_names` - List of substance identifiers
    /// * `library_name` - Database name (e.g., "Aramco_transport")
    ///
    /// # Returns
    /// * `Ok(())` - All substances loaded successfully
    /// * `Err(TransportError)` - Library or substance not found
    pub fn add_substances_from_library(
        &mut self,
        substance_names: Vec<String>,
        library_name: &str,
    ) -> Result<(), TransportError> {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get(library_name).ok_or_else(|| {
            TransportError::MissingData(format!("Library {} not found", library_name))
        })?;

        for name in substance_names {
            let data = sublib.get(&name).ok_or_else(|| {
                TransportError::MissingData(format!("Substance {} not found", name.to_string()))
            })?;

            let mut transport_data = TransportData::new();
            transport_data.from_serde(data.clone())?;
            transport_data.P = self.P;

            self.add_substance(name, transport_data);
        }
        Ok(())
    }

    /// Calculates diffusion coefficients for all substance pairs
    ///
    /// Uses double loop to compute N(N+1)/2 unique pairs.
    /// Must set molar masses first via `set_molar_masses()`.
    pub fn calculate_all_pairs(&mut self) {
        self.pair_diffusions.clear();

        for (i, name_a) in self.substance_names.iter().enumerate() {
            for (j, name_b) in self.substance_names.iter().enumerate() {
                if i <= j {
                    let transport_a = &self.substances[name_a];
                    let transport_b = &self.substances[name_b];

                    let mut pair_diff = PairDiffusion::new(
                        transport_a.input.clone(),
                        transport_b.input.clone(),
                        self.T,
                        self.P,
                        transport_a.M,
                        transport_b.M,
                    );

                    if let Some(ref unit) = self.D_unit {
                        let _ = pair_diff.set_D_unit(Some(unit.clone()));
                    }

                    pair_diff.calc_coefficients();
                    pair_diff.calculate_D();

                    self.pair_diffusions
                        .insert((name_a.clone(), name_b.clone()), pair_diff);
                }
            }
        }
    }

    /// Creates temperature-dependent closures for all pairs
    ///
    /// Must call `calculate_all_pairs()` first.
    pub fn calculate_all_closures(&mut self) {
        for pair_diff in self.pair_diffusions.values_mut() {
            pair_diff.calculate_D_closure();
        }
    }

    /// Creates symbolic expressions for all pairs
    ///
    /// Must call `calculate_all_pairs()` first.
    pub fn calculate_all_symbolic(&mut self) {
        for pair_diff in self.pair_diffusions.values_mut() {
            pair_diff.calculate_D_symbolic();
        }
    }

    /// Creates symmetric N×N diffusion coefficient matrix
    ///
    /// # Returns
    /// Symmetric matrix where matrix[i][j] = D_ij [m²/s or cm²/s]
    pub fn create_diffusion_matrix(&mut self) -> Vec<Vec<f64>> {
        let n = self.substance_names.len();
        let mut matrix = vec![vec![0.0; n]; n];

        for (i, name_a) in self.substance_names.iter().enumerate() {
            for (j, name_b) in self.substance_names.iter().enumerate() {
                let key = if i <= j {
                    (name_a.clone(), name_b.clone())
                } else {
                    (name_b.clone(), name_a.clone())
                };

                if let Some(pair_diff) = self.pair_diffusions.get(&key) {
                    matrix[i][j] = pair_diff.D;
                    matrix[j][i] = pair_diff.D;
                }
            }
        }

        self.diffusion_matrix = Some(matrix.clone());
        matrix
    }

    /// Retrieves PairDiffusion data for specific substance pair
    ///
    /// # Arguments
    /// * `sub_a` - First substance name
    /// * `sub_b` - Second substance name
    ///
    /// # Returns
    /// * `Some(&PairDiffusion)` - Pair data found
    /// * `None` - Pair not calculated
    pub fn get_pair_diffusion(&self, sub_a: &str, sub_b: &str) -> Option<&PairDiffusion> {
        let key1 = (sub_a.to_string(), sub_b.to_string());
        let key2 = (sub_b.to_string(), sub_a.to_string());

        self.pair_diffusions
            .get(&key1)
            .or_else(|| self.pair_diffusions.get(&key2))
    }

    /// Gets diffusion coefficient value for specific pair
    ///
    /// # Arguments
    /// * `sub_a` - First substance name
    /// * `sub_b` - Second substance name
    ///
    /// # Returns
    /// * `Some(f64)` - Diffusion coefficient [m²/s or cm²/s]
    /// * `None` - Pair not found
    pub fn get_diffusion_coefficient(&self, sub_a: &str, sub_b: &str) -> Option<f64> {
        self.get_pair_diffusion(sub_a, sub_b).map(|pd| pd.D)
    }

    /// Sets common unit for all diffusion coefficients
    ///
    /// # Arguments
    /// * `unit` - Unit string ("m2/s" or "cm2/s")
    pub fn set_common_D_unit(&mut self, unit: String) {
        self.D_unit = Some(unit);
    }

    /// Prints formatted diffusion matrix to console
    ///
    /// Must call `create_diffusion_matrix()` first.
    pub fn print_diffusion_matrix(&self) {
        if let Some(ref matrix) = self.diffusion_matrix {
            println!("Diffusion Matrix:");
            print!("\t");
            for name in &self.substance_names {
                print!("{:>12}", name);
            }
            println!();

            for (i, name) in self.substance_names.iter().enumerate() {
                print!("{:>8}", name);
                for j in 0..self.substance_names.len() {
                    print!("{:>12.2e}", matrix[i][j]);
                }
                println!();
            }
        }
    }

    /// Sets molar masses for multiple substances
    ///
    /// # Arguments
    /// * `masses` - HashMap mapping substance names to molar masses [kg/mol]
    pub fn set_molar_masses(&mut self, masses: HashMap<String, f64>) {
        for (name, mass) in masses {
            if let Some(transport_data) = self.substances.get_mut(&name) {
                transport_data.M = mass;
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Thermodynamics::DBhandlers::TRANSPORTdata::TransportData;
    use crate::Thermodynamics::thermo_lib_api::ThermoData;
    use approx::assert_relative_eq;
    #[test]
    fn test_pair_diffusion() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let N2_data = sublib.get("N2").unwrap();
        let mut tr1 = TransportData::new();
        tr1.from_serde(CO_data.clone()).unwrap();
        let CO_input = tr1.input;
        let mut tr2 = TransportData::new();
        tr2.from_serde(N2_data.clone()).unwrap();
        let N2_input = tr2.input;
        let mut pair_diff = PairDiffusion::new(
            CO_input,
            N2_input,
            300.0,
            101325.0,
            28.01 / 1000.0,
            30.01 / 1000.0,
        );
        let _ = pair_diff.set_D_unit(Some("m2/s".to_string()));
        pair_diff.calc_coefficients();
        pair_diff.calculate_D();
        println!("D_CO_NO at 300K: {} m^2/s* 1e6", pair_diff.D * 1e6);
    }

    #[test]
    fn test_closure() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let N2_data = sublib.get("N2").unwrap();
        let mut tr1 = TransportData::new();
        tr1.from_serde(CO_data.clone()).unwrap();
        let CO_input = tr1.input;
        let mut tr2 = TransportData::new();
        tr2.from_serde(N2_data.clone()).unwrap();
        let N2_input = tr2.input;
        let mut pair_diff = PairDiffusion::new(
            CO_input,
            N2_input,
            300.0,
            101325.0,
            28.01 / 1000.0,
            30.01 / 1000.0,
        );
        pair_diff.calc_coefficients();
        pair_diff.calculate_D_closure();

        let D_300 = (pair_diff.D_closure)(300.0);
        let D_400 = (pair_diff.D_closure)(400.0);

        assert!(D_400 > D_300);
        println!("D at 300K: {}, D at 400K: {}", D_300, D_400);
    }

    #[test]
    fn test_symbolic() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let N2_data = sublib.get("N2").unwrap();
        let mut tr1 = TransportData::new();
        tr1.from_serde(CO_data.clone()).unwrap();
        let CO_input = tr1.input;
        let mut tr2 = TransportData::new();
        tr2.from_serde(N2_data.clone()).unwrap();
        let N2_input = tr2.input;
        let mut pair_diff = PairDiffusion::new(
            CO_input,
            N2_input,
            300.0,
            101325.0,
            28.01 / 1000.0,
            30.01 / 1000.0,
        );
        pair_diff.calc_coefficients();
        pair_diff.calculate_D_symbolic();

        let D_sym_func = pair_diff.D_symbolic.lambdify1D();
        let D_300 = D_sym_func(300.0);
        let D_400 = D_sym_func(400.0);

        assert!(D_400 > D_300);
        println!("Symbolic D at 300K: {}, at 400K: {}", D_300, D_400);
    }

    #[test]
    fn test_taylor_series() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let N2_data = sublib.get("N2").unwrap();
        let mut tr1 = TransportData::new();
        tr1.from_serde(CO_data.clone()).unwrap();
        let CO_input = tr1.input;
        let mut tr2 = TransportData::new();
        tr2.from_serde(N2_data.clone()).unwrap();
        let N2_input = tr2.input;
        let mut pair_diff = PairDiffusion::new(
            CO_input,
            N2_input,
            300.0,
            101325.0,
            28.01 / 1000.0,
            30.01 / 1000.0,
        );
        pair_diff.calc_coefficients();
        pair_diff.calculate_D_symbolic();

        let original_func = pair_diff.D_symbolic.lambdify1D();
        let original_300 = original_func(300.0);

        pair_diff.taylor_series(300.0, 3);
        let taylor_func = pair_diff.D_symbolic.lambdify1D();
        let taylor_300 = taylor_func(300.0);

        assert_relative_eq!(original_300, taylor_300, epsilon = 1e-6);
        println!("Original: {}, Taylor: {}", original_300, taylor_300);
    }

    #[test]
    fn test_integr_mean() {
        let thermo_data = ThermoData::new();
        let sublib = thermo_data.LibThermoData.get("Aramco_transport").unwrap();
        let CO_data = sublib.get("CO").unwrap();
        let N2_data = sublib.get("N2").unwrap();
        let mut tr1 = TransportData::new();
        tr1.from_serde(CO_data.clone()).unwrap();
        let CO_input = tr1.input;
        let mut tr2 = TransportData::new();
        tr2.from_serde(N2_data.clone()).unwrap();
        let N2_input = tr2.input;
        let mut pair_diff = PairDiffusion::new(
            CO_input,
            N2_input,
            300.0,
            101325.0,
            28.01 / 1000.0,
            30.01 / 1000.0,
        );
        pair_diff.calc_coefficients();
        pair_diff.calculate_D_symbolic();
        pair_diff.set_T_interval(300.0, 400.0);

        let result = pair_diff.integr_mean();
        assert!(result.is_ok());
        println!("Mean D over 300-400K: {}", pair_diff.D);
    }

    #[test]
    fn test_unit_conversion() {
        let mut pair_diff = PairDiffusion::default();

        assert!(pair_diff.set_D_unit(Some("cm2/s".to_string())).is_ok());
        assert_eq!(pair_diff.unit_multiplier, 1e4);

        assert!(pair_diff.set_D_unit(Some("m2/s".to_string())).is_ok());
        assert_eq!(pair_diff.unit_multiplier, 1.0);

        assert!(pair_diff.set_D_unit(Some("invalid".to_string())).is_err());
    }
    /////////////////////////////////////
    #[test]
    fn test_multi_substance_diffusion() {
        let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);

        let result = multi_diff.add_substances_from_library(
            vec!["CO".to_string(), "N2".to_string(), "O2".to_string()],
            "Aramco_transport",
        );
        assert!(result.is_ok());

        let masses = [
            ("CO", 28.01 / 1000.0),
            ("N2", 28.014 / 1000.0),
            ("O2", 31.998 / 1000.0),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), *v))
        .collect();
        multi_diff.set_molar_masses(masses);

        multi_diff.set_common_D_unit("cm2/s".to_string());
        multi_diff.calculate_all_pairs();

        let matrix = multi_diff.create_diffusion_matrix();
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);

        let d_co_n2 = multi_diff.get_diffusion_coefficient("CO", "N2");
        assert!(d_co_n2.is_some());

        multi_diff.print_diffusion_matrix();
        println!("CO-N2 diffusion: {:?}", d_co_n2);
    }
    #[test]
    fn test_multi_substance_closures() {
        let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);

        let _ = multi_diff.add_substances_from_library(
            vec!["CO".to_string(), "N2".to_string()],
            "Aramco_transport",
        );

        let masses = [("CO", 28.01 / 1000.0), ("N2", 28.014 / 1000.0)]
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        multi_diff.set_molar_masses(masses);

        multi_diff.calculate_all_pairs();
        multi_diff.calculate_all_closures();

        let pair_diff = multi_diff.get_pair_diffusion("CO", "N2").unwrap();
        let D_300 = (pair_diff.D_closure)(300.0);
        let D_400 = (pair_diff.D_closure)(400.0);

        assert!(D_400 > D_300);
    }

    #[test]
    fn test_multi_substance_symbolic() {
        let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);

        let _ = multi_diff.add_substances_from_library(
            vec!["CO".to_string(), "N2".to_string()],
            "Aramco_transport",
        );

        let masses = [("CO", 28.01 / 1000.0), ("N2", 28.014 / 1000.0)]
            .iter()
            .map(|(k, v)| (k.to_string(), *v))
            .collect();
        multi_diff.set_molar_masses(masses);

        multi_diff.calculate_all_pairs();
        multi_diff.calculate_all_symbolic();

        let pair_diff = multi_diff.get_pair_diffusion("CO", "N2").unwrap();
        let D_sym_func = pair_diff.D_symbolic.lambdify1D();
        let D_300 = D_sym_func(300.0);
        let D_400 = D_sym_func(400.0);

        assert!(D_400 > D_300);
    }

    #[test]
    fn test_matrix_symmetry() {
        let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);

        let _ = multi_diff.add_substances_from_library(
            vec!["CO".to_string(), "N2".to_string(), "O2".to_string()],
            "Aramco_transport",
        );

        let masses = [
            ("CO", 28.01 / 1000.0),
            ("N2", 28.014 / 1000.0),
            ("O2", 31.998 / 1000.0),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), *v))
        .collect();
        multi_diff.set_molar_masses(masses);

        multi_diff.calculate_all_pairs();
        let matrix = multi_diff.create_diffusion_matrix();

        for i in 0..matrix.len() {
            for j in 0..matrix[i].len() {
                assert_relative_eq!(matrix[i][j], matrix[j][i], epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_error_handling() {
        let mut multi_diff = MultiSubstanceDiffusion::new(300.0, 101325.0);

        let result = multi_diff
            .add_substances_from_library(vec!["INVALID_SUBSTANCE".to_string()], "Aramco_transport");
        assert!(result.is_err());

        let result =
            multi_diff.add_substances_from_library(vec!["CO".to_string()], "INVALID_LIBRARY");
        assert!(result.is_err());
    }
}
