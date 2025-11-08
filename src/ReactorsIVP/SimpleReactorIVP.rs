//! NOT READY!

//! # Simple Reactor IVP Module
//!
//! This module provides a comprehensive framework for modeling chemical reactors using boundary value problems (IVP).
//! It implements dimensionless reactor equations for mass and heat transfer with chemical reactions.
//!
//! ## Main Structures
//!
//! - **`SimpleReactorTask`**: Main reactor modeling structure that aggregates kinetics, thermodynamics, and transport properties
//! - **`IVPSolver`**: Wrapper for different IVP solvers
//! - **`ToleranceConfig`**: Helper for setting solver tolerances across all variables (C, J, Teta, q)
//! - **`BoundsConfig`**: Helper for setting variable bounds across all variables
//! - **`FastElemReact`**: Simple structure for elementary reactions with Arrhenius parameters
//!
//! ## Key Features
//!
//! - **Automatic dimensionless scaling**: Converts dimensional equations to dimensionless form for numerical stability
//! - **Flexible tolerance/bounds setup**: Automatically expands simple configs to full variable maps (C0,C1,... J0,J1,...)
//! - **Multiple solver backends**: Supports both damped Newton-Raphson and adaptive mesh refinement solvers
//! - **Transport property calculation**: Automatic calculation of Peclet numbers and Lewis numbers
//! - **Post-processing**: Converts dimensionless results back to dimensional form
//!
//! ## Mathematical Model
//!
//! The module solves dimensionless reactor equations:
//! - **Mass balance**: dCi/dz = Ji/(D*ρ), dJi/dz = Pe_D*Ji - l²*Gi
//! - **Energy balance**: dTeta/dz = q/λ, dq/dz = Pe_q*q - l²*Q/dT
//!
//! Where z = x/L (dimensionless coordinate), Teta = (T-T₀)/dT (dimensionless temperature)
//! PAY ATTENTION TO THE DIMENSION OF INPUT PARAMETERS
use crate::Kinetics::User_reactions::KinData;
use crate::Kinetics::mechfinder_api::ReactionData;
use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
use crate::ReactorsBVP::reactor_BVP_utils::{
    BoundsConfig, ScalingConfig, ToleranceConfig, create_bounds_map, create_tolerance_map,
};

use RustedSciThe::symbolic::symbolic_engine::Expr;
use log::info;

use nalgebra::{DMatrix, DVector};

use std::collections::HashMap;
use std::fmt;

/// Universal gas constant in J/(mol·K)
pub const R_G: f64 = 8.314;

#[derive(Debug, Clone)]
pub struct SolutionQuality {
    pub energy_balane_error_abs: f64,
    pub energy_balane_error_rel: f64,

    /// steps where sum of molar fractions is larger then threshhold
    pub sum_of_mass_fractions: Vec<(usize, f64)>,
    pub atomic_mass_balance_error: Vec<(usize, f64)>,
}

impl Default for SolutionQuality {
    fn default() -> Self {
        Self {
            energy_balane_error_abs: 0.0,
            energy_balane_error_rel: 0.0,

            sum_of_mass_fractions: Vec::new(),
            atomic_mass_balance_error: Vec::new(),
        }
    }
}

/*
       (1/l)*  d(Lambda*( (1/l)* dT/d(x/l)  )  )/d(x/l) - c*m*( (1/l)* dT/d(x/l) + Q =0
        define: z = x/l
        d(Lambda* dT/dz   )/dz - c*m*l* dT/dz + l^2*Q =0
        define  Teta = (T-dT)/dT
         d(Lambda* dTeta/dz   )/dz - c*m*l* dTeta/dz + l^2*Q/dT =0
        define: q = Lambda* dTeta/dz
         d(Lambda* dTeta/dz   )/dz - c*m*l* dT/dz + l^2*Q/dT =0
        got equations:
        dT/dz = q/Lambda
        dq/dz - (c*m*l)/Lambda *q + (l^2)*Q/dT =0
        and finally:
        dT/dz = q/Lambda
        dq/dz - Pe_q*q + (l^2)*Q/dT =0

        (1/l)*  d(D*ro*( (1/l)* dCi/d(x/l)  )  )/d(x/l) - m*( (1/l)* dCi/d(x/l) + Gi =0
        define: z = x/l
        d(D*ro* dCi/dz   )/dz - m*l* dT/dz + l^2*Gi =0
        define: Ji = D*ro dC/dz
         d(Ji)/dz - m*l* dС/dz + l^2*Gi =0
        got equations:
        dCi/dz = Ji/D*ro;
        dJi/dz - (m*l)/ro*D *Ji + (l^2)*Gi =0

*/
/// Simple structure for elementary chemical reactions with Arrhenius kinetics
///
/// Rate = A * T^n * exp(-E/(R*T)) * ∏[Ci]^νi
#[derive(Debug, Clone)]
pub struct FastElemReact {
    /// Chemical equation (e.g., "A + B => C + D")
    pub eq: String,
    /// Pre-exponential factor (units depend on reaction order)
    pub A: f64,
    /// Temperature exponent (dimensionless)
    pub n: f64,
    /// Activation energy (J/mol)
    pub E: f64,
    /// Heat of reaction (J/kg)
    pub Q: f64,
}

/// Main reactor modeling structure that aggregates all reactor properties and methods
///
/// This structure handles the complete workflow from kinetic data to IVP solution:
/// 1. Kinetic preprocessing (stoichiometry, rate expressions)
/// 2. Transport property calculations (Peclet numbers, diffusion)
/// 3. Dimensionless scaling and equation setup
/// 4. IVP solving with multiple solver options
/// 5. Post-processing and visualization
#[derive(Debug, Clone)]
pub struct SimpleReactorTask {
    /// Optional problem identifier
    pub problem_name: Option<String>,
    /// Optional problem description
    pub problem_description: Option<String>,
    /// Kinetic data (reactions, substances, rate constants)
    pub kindata: KinData,
    /// Heat effects for each reaction (J/kg)
    pub thermal_effects: Vec<f64>,
    /// Pressure (Pa)
    pub P: f64,
    /// Mean temperature (K)
    pub Tm: f64,
    /// Heat capacity (J/kg·K)
    pub Cp: f64,
    /// Boundary conditions for substances and temperature
    pub boundary_condition: HashMap<String, f64>,
    /// Thermal conductivity (W/m·K)
    pub Lambda: f64,
    pub ro: f64,
    /// Mass flow rate (kg/s)
    pub m: f64,
    /// Scaling parameters for dimensionless transformation
    pub scaling: ScalingConfig,
    /// Characteristic length (m)
    pub L: f64,
    /// Temperature scaling expression: T = dT*(Teta + 1)
    pub T_scaling: Expr,
    /// Mean molar mass (kg/mol)
    pub M: f64,

    /// Thermal Peclet number: Pe_q = L*m*Cp/λ
    pub Pe_q: f64,

    /// Reaction rate expressions for each reaction
    pub map_eq_rate: HashMap<String, Expr>,
    /// System of differential equations (substance -> (variable, equation))
    pub map_of_equations: HashMap<String, (String, Expr)>,
    /// heat release function
    pub heat_release: Expr,
    /// IVP solver instance
    pub solver: IVPSolver,
}
/// Boundary Value Problem solver wrapper
///
/// Supports multiple solver backends:

#[derive(Debug, Clone)]
pub struct IVPSolver {
    /// Independent variable name (typically "x" or "z")
    pub arg_name: String,
    /// Domain range (start, end) - typically (0.0, 1.0) for dimensionless
    pub x_range: (f64, f64),
    /// Names of unknown variables ["Teta", "q", "C0", "J0", "C1", "J1", ...]
    pub unknowns: Vec<String>,
    /// System of differential equations dy/dx = f(x,y)
    pub eq_system: Vec<Expr>,
    /// Boundary conditions: variable -> (boundary_index, value)
    pub BorderConditions: HashMap<String, (usize, f64)>,
    /// Solution matrix (variables × mesh_points)
    pub solution: Option<DMatrix<f64>>,
    /// Spatial mesh points
    pub x_mesh: Option<DVector<f64>>,
    /// struct that stores balance errors - filled at the end of solution
    pub quality: SolutionQuality,
}

impl Default for IVPSolver {
    fn default() -> Self {
        Self {
            arg_name: "x".to_string(),
            x_range: (0.0, 1.0),
            unknowns: Vec::new(),
            eq_system: Vec::new(),
            BorderConditions: HashMap::new(),
            solution: None,
            x_mesh: None,
            quality: SolutionQuality::default(),
        }
    }
}
impl IVPSolver {
    /// Get reference to solution matrix
    ///
    /// Returns None if solve hasn't been called yet
    pub fn get_solution(&self) -> Option<&DMatrix<f64>> {
        self.solution.as_ref()
    }

    /// Debug print solution summary
    ///
    /// Prints solution matrix dimensions and sample values for each variable
    pub fn debug_solution(&self) {
        if let Some(solution) = &self.solution {
            println!("\n=== SOLUTION DEBUG ===");
            println!(
                "Solution matrix shape: {} x {}",
                solution.nrows(),
                solution.ncols()
            );
            println!("Unknowns: {:?}", self.unknowns);

            // Print first and lust few values for each variable
            for (i, var_name) in self.unknowns.iter().enumerate() {
                if i < solution.ncols() {
                    let col = solution.column(i);
                    println!(
                        "{}: [first...{:.6}, {:.6}, {:.6}, ... last: {:.6}, {:.6}, {:.6}]",
                        var_name,
                        col[0],
                        col[1.min(col.len() - 1)],
                        col[2.min(col.len() - 1)],
                        col[col.len() - 3],
                        col[col.len() - 2],
                        col[col.len() - 1]
                    );
                }
            }
            println!("=== END DEBUG ===\n");
        }
    }
}

impl SimpleReactorTask {
    /// Create new reactor task with default values
    pub fn new() -> Self {
        Self {
            problem_name: None,
            problem_description: None,
            kindata: KinData::new(),
            thermal_effects: Vec::new(),
            P: 0.0,
            Tm: 0.0,
            Cp: 0.0,
            boundary_condition: HashMap::new(),
            Lambda: 0.0,
            ro: 0.0,

            m: 0.0,
            scaling: ScalingConfig::default(),
            L: 1.0,
            T_scaling: Expr::Const(0.0),
            M: 0.0,

            Pe_q: 0.0,
            map_eq_rate: HashMap::new(),
            map_of_equations: HashMap::new(),
            heat_release: Expr::Const(0.0),
            solver: IVPSolver::default(),
        }
    }

    /// Create tolerance map from simplified config for this reactor's substances
    pub fn create_tolerance_map_for_system(
        &self,
        tolerance_config: HashMap<String, f64>,
    ) -> HashMap<String, f64> {
        create_tolerance_map(tolerance_config, &self.kindata.substances)
    }

    /// Create bounds map from simplified config for this reactor's substances
    pub fn create_bounds_map_for_system(
        &self,
        bounds_config: HashMap<String, (f64, f64)>,
    ) -> HashMap<String, (f64, f64)> {
        create_bounds_map(bounds_config, &self.kindata.substances)
    }

    /// Set scaling parameters using ScalingConfig
    pub fn set_scaling(&mut self, scaling: ScalingConfig) -> Result<(), ReactorError> {
        scaling.validate()?;
        self.scaling = scaling;
        Ok(())
    }

    /// Set scaling parameters from individual values
    pub fn set_scaling_values(
        &mut self,
        dT: f64,
        L: f64,
        T_scale: f64,
    ) -> Result<(), ReactorError> {
        let scaling = ScalingConfig::new(dT, L, T_scale);
        self.set_scaling(scaling)
    }
    /////////////////////////////////SETTERS////////////////////////////////////////////////////////////////////////////////
    /// Set problem name for identification
    pub fn set_problem_name(&mut self, name: &str) {
        self.problem_name = Some(name.to_string());
    }

    /// Set problem description
    pub fn set_problem_description(&mut self, description: &str) {
        self.problem_description = Some(description.to_string());
    }
    /// Set boundary conditions for substances and temperature
    ///
    /// Keys should include substance names and "T" for temperature
    pub fn set_boundary_conditions(&mut self, conditions: HashMap<String, f64>) {
        self.boundary_condition = conditions;
    }
    /// Set all reactor parameters at once
    ///
    /// Convenient method to set all physical and transport properties
    pub fn set_parameters(
        &mut self,
        thermal_effects: Vec<f64>,
        P: f64,
        Tm: f64,
        Cp: f64,
        boundary_condition: HashMap<String, f64>,
        Lambda: f64,
        Diffusion: HashMap<String, f64>,
        m: f64,

        scaling: ScalingConfig,
    ) {
        self.thermal_effects = thermal_effects;
        self.P = P;
        self.Tm = Tm;
        self.Cp = Cp;
        self.boundary_condition = boundary_condition;
        self.Lambda = Lambda;

        self.m = m;

        self.scaling = scaling;
    }
    /// Complete IVP setup workflow
    ///
    /// Orchestrates the entire setup process:
    /// 1. Kinetic preprocessing
    /// 2. Scaling and transport calculations  
    /// 3. Equation system assembly
    /// 4. Boundary condition setup
    pub fn setup_IVP(&mut self) -> Result<(), ReactorError> {
        //
        self.check_task()?;
        info!("task checked!");
        // Process
        self.scaling_processing()?;
        info!("scaling processed!");
        //
        // Process kinetics
        self.kinetic_processing()?;
        info!("kinetics processed!");

        // Calculate mean molar mass
        self.mean_molar_mass()?;
        info!("mean molar mass calculated");

        // dbg!("here");
        // Calculate Peclet numbers
        self.peclet_numbers()?;
        info!("Peclet numbers calculated");
        // Create IVP equations
        self.create_IVP_equations()?;
        info!("IVP equations created");
        self.set_solver_BC()?;
        info!("Boundary conditions created!");
        self.check_before_solution()?;
        info!("IVP setup completed!");
        self.pretty_print_task();
        self.pretty_print_equations();
        self.pretty_print_reaction_rates();

        Ok(())
    }

    /// Set transport properties
    pub fn set_transport_properties(&mut self, lambda: f64, cp: f64) {
        self.Lambda = lambda;
        self.Cp = cp;
    }

    /// Set operating conditions
    pub fn set_operating_conditions(&mut self, pressure: f64, temperature: f64, mass_flow: f64) {
        self.P = pressure;
        self.Tm = temperature;
        self.m = mass_flow;
    }

    /// Set thermal effects (heat of reaction) for each reaction
    pub fn set_thermal_effects(&mut self, thermal_effects: Vec<f64>) {
        self.thermal_effects = thermal_effects;
    }

    /// Set elementary reactions from FastElemReact structures
    ///
    /// Convenient method to quickly set up elementary reactions with Arrhenius kinetics
    pub fn fast_react_set(&mut self, vec_of_maps: Vec<FastElemReact>) -> Result<(), ReactorError> {
        let mut eq_vec: Vec<String> = Vec::new();
        let mut elementary_reaction_vec = Vec::new();
        let mut Q_vec = Vec::new();
        for (idx, map_of_reactiondata) in vec_of_maps.iter().enumerate() {
            // Check equation string
            let eq = if !map_of_reactiondata.eq.is_empty() {
                map_of_reactiondata.eq.clone()
            } else {
                return Err(ReactorError::MissingData(format!(
                    "No equation in input hashmap at index {}",
                    idx
                )));
            };

            // Check Arrhenius parameters
            let A = map_of_reactiondata.A;
            let n = map_of_reactiondata.n;
            let E = map_of_reactiondata.E;
            let Q = map_of_reactiondata.Q;
            // Check for NaN (in case of uninitialized f64)
            if A.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'A' in input hashmap at index {}",
                    idx
                )));
            }
            if n.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'n' in input hashmap at index {}",
                    idx
                )));
            }
            if E.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'E' in input hashmap at index {}",
                    idx
                )));
            }
            if Q.is_nan() {
                return Err(ReactorError::MissingData(format!(
                    "Missing Arrhenius parameter 'Q' in input hashmap at index {}",
                    idx
                )));
            }
            let arrenius = vec![A, n, E];
            let reactdata = ReactionData::new_elementary(eq.clone(), arrenius, None);
            eq_vec.push(eq);
            Q_vec.push(Q);
            elementary_reaction_vec.push(reactdata);
        }

        let mut kindata = KinData::new();
        kindata.vec_of_equations = eq_vec;
        kindata.vec_of_reaction_data = Some(elementary_reaction_vec);
        self.kindata = kindata;
        self.thermal_effects = Q_vec;
        Ok(())
    }
    ///////////////////////////////////////////VALIDATION////////////////////////////////////////////////
    /// Validate reactor task configuration
    ///
    /// Checks:
    /// - All physical properties are positive
    /// - Diffusion coefficients exist for all substances
    /// - Thermal effects match number of reactions
    /// - Scaling parameters are valid
    /// - Boundary conditions are complete
    pub fn check_task(&self) -> Result<(), ReactorError> {
        // Check basic properties
        if self.P <= 0.0 {
            return Err(ReactorError::MissingData("P must be positive".to_string()));
        }
        if self.Tm <= 0.0 {
            return Err(ReactorError::MissingData("Tm must be positive".to_string()));
        }
        if self.Cp <= 0.0 {
            return Err(ReactorError::MissingData("Cp must be positive".to_string()));
        }
        if self.Lambda <= 0.0 {
            return Err(ReactorError::MissingData(
                "Lambda must be positive".to_string(),
            ));
        }
        if self.m <= 0.0 {
            return Err(ReactorError::MissingData("m must be positive".to_string()));
        }

        // Check thermal effects length
        if self.thermal_effects.len() != self.kindata.vec_of_equations.len() {
            return Err(ReactorError::InvalidConfiguration(
                "Thermal effects length must match number of reactions".to_string(),
            ));
        }

        // Validate scaling parameters
        self.scaling.validate()?;

        // Check boundary conditions
        if !self.boundary_condition.contains_key("T") {
            return Err(ReactorError::MissingData(
                "Missing T in boundary conditions".to_string(),
            ));
        }
        for substance in &self.kindata.substances {
            if !self.boundary_condition.contains_key(substance) {
                return Err(ReactorError::MissingData(format!(
                    "Missing boundary condition for {}",
                    substance
                )));
            }
        }

        Ok(())
    }
    /// Validate system before solving
    ///
    /// Checks that all arrays have consistent dimensions:
    /// - Equation system length = 2*n_substances + 2
    /// - Unknown variables, equations, and boundary conditions match
    /// - Peclet numbers are calculated and positive
    pub fn check_before_solution(&self) -> Result<(), ReactorError> {
        let n_substances = self.kindata.substances.len();
        let expected_len = n_substances + 2;

        if self.solver.eq_system.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "eq_system length {} != expected {}",
                self.solver.eq_system.len(),
                expected_len
            )));
        }
        if self.solver.unknowns.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "unknowns length {} != expected {}",
                self.solver.unknowns.len(),
                expected_len
            )));
        }
        if self.map_of_equations.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "map_of_equations length {} != expected {}",
                self.map_of_equations.len(),
                expected_len
            )));
        }
        if self.solver.BorderConditions.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "BorderConditions length {} != expected {}",
                self.solver.BorderConditions.len(),
                expected_len
            )));
        }

        if self.M <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "M must be positive".to_string(),
            ));
        }
        if self.Pe_q <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Pe_q must be positive".to_string(),
            ));
        }
        Ok(())
    }
    ///////////////////////////////////////////KINETICS AND THERMAL PREPROCESSING////////////////////////////////////////////////

    ///
    pub fn kinetic_processing(&mut self) -> Result<(), ReactorError> {
        let kd = &mut self.kindata;
        // stoichiometry and element matrix
        kd.analyze_reactions();
        // in elementary reactions there are only Arrhenius parameters - no concentration or pressure dependencies
        kd.calc_sym_constants(None, None, Some(self.T_scaling.clone()));
        Ok(())
    }

    ///
    pub fn mean_molar_mass(&mut self) -> Result<(), ReactorError> {
        println!("DEBUG mean_molar_mass: Entering function");
        println!(
            "DEBUG mean_molar_mass: Current molar masses: {:?}",
            self.kindata.stecheodata.vec_of_molmasses
        );
        let mut mean_mass_inv = 0.0;
        let molar_masses = self
            .kindata
            .stecheodata
            .vec_of_molmasses
            .as_ref()
            .ok_or_else(|| ReactorError::MissingData("Molar masses not calculated".to_string()))?;
        println!(
            "DEBUG mean_molar_mass: Using molar masses: {:?}",
            molar_masses
        );

        // M_mean = sum_i(xi*Mi)
        //   x(i) = (ω(i) / M(i)) / ∑(ω(j) / M(j)) where ω(j) - is mass fruction
        // so M_men = ∑(xi*Mi) = ∑ω(i)  / ∑(ω(j) / M(j))= 1/∑(ω(j) / M(j))
        if self.M == 0.0 || self.M.is_nan() {
            for (i, substance) in self.kindata.substances.iter().enumerate() {
                if let Some(conc) = self.boundary_condition.get(substance) {
                    let mol_mass = molar_masses.get(i).ok_or_else(|| {
                        ReactorError::IndexOutOfBounds(format!(
                            "Molar mass index {} out of bounds",
                            i
                        ))
                    })?;
                    mean_mass_inv += conc / mol_mass;
                }
            }
            let mean_mass = 1.0 / mean_mass_inv;
            self.M = mean_mass / 1000.0; // from g/mol to kg/mol
        }
        Ok(())
    }
    ///
    /// Process scaling parameters for dimensionless equations
    /// Teta = (T - dT)/T_scaling
    /// Creates temperature scaling: T = (Teta*T_scaling + dT)
    /// Sets characteristic length L for spatial scaling: z = x/L
    pub fn scaling_processing(&mut self) -> Result<(), ReactorError> {
        // Validate scaling parameters
        self.scaling.validate()?;

        // Create temperature scaling expression: T = dT*(Teta + 1)
        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        let Teta = Expr::Var("Teta".to_owned());
        self.T_scaling = Teta.clone() * Expr::Const(T_scale) + Expr::Const(dT);

        // Set characteristic length
        self.L = self.scaling.L;
        Ok(())
    }
    pub fn peclet_numbers(&mut self) -> Result<(), ReactorError> {
        if self.Lambda <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Lambda must be positive".to_string(),
            ));
        }
        if self.m <= 0.0 {
            return Err(ReactorError::InvalidConfiguration(
                "Mass flow rate must be positive".to_string(),
            ));
        }
        if self.L.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing L".to_string()));
        }
        if self.m.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing m".to_string()));
        }
        if self.Cp.is_nan() {
            return Err(ReactorError::InvalidConfiguration("missing Cp".to_string()));
        }
        if self.Lambda.is_nan() {
            return Err(ReactorError::InvalidConfiguration(
                "missing Lambda".to_string(),
            ));
        }

        let Pe_q = (self.L * self.m * self.Cp) / self.Lambda;

        self.Pe_q = Pe_q;

        Ok(())
    }
    pub fn ideal_gas_density(&self) -> f64 {
        self.M * self.P / (R_G * self.Tm)
    }
}

/////////////////////////////////////////TESTS/////////////////////////////////////////////////////////
