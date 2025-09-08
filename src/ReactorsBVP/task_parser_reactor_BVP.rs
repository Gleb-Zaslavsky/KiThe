//! # Task Parser for Reactor BVP Module
//!
//! This module provides comprehensive parsing capabilities for chemical reactor boundary value problems (BVPs)
//! from structured text files. It bridges the gap between human-readable configuration files and the complex
//! data structures required by numerical BVP solvers.
//!
//! ## Purpose
//!
//! The module addresses the challenge of configuring complex reactor simulations through text files rather
//! than hardcoded parameters. It enables users to define complete reactor problems including:
//! - Chemical kinetics (reactions, Arrhenius parameters)
//! - Thermodynamic properties (heat capacity, thermal conductivity)
//! - Transport properties (diffusion coefficients)
//! - Boundary conditions and initial conditions
//! - Solver settings (tolerances, bounds, numerical methods)
//!
//! ## Main Methods
//!
//! - **`set_reactor_params_from_hashmap()`**: Core parsing method that extracts reactor parameters from DocumentMap
//! - **`parse_toleranse_and_bounds()`**: Specialized parser for solver tolerance and bounds configurations
//! - **`parse_file()`**: File I/O wrapper that creates DocumentParser from file path
//! - **`parse_parameters_with_exact_names()`**: High-level parsing orchestrator
//! - **`solve_from_map()`**: Complete workflow from parsed data to solved BVP
//! - **`solve_from_file()`**: One-shot method: file → parsing → solving
//!
//! ## Non-Obvious Code Features & Tips
//!
//! ### Variable Shadowing Prevention
//! The code carefully avoids variable name conflicts when parsing nested structures.
//! Example: `problem_name_val` instead of reusing `process_conditions` in if-let statements.
//!
//! ### DocumentMap Navigation Pattern
//! Uses consistent `.get().expect().clone().unwrap()[0].as_type()` pattern for extracting
//! values from the nested HashMap structure returned by the document parser.
//!
//! ### Groups Parsing Logic
//! The groups parsing implements a two-stage process: first check if groups are enabled,
//! then iterate through substances to build the nested HashMap structure.
//!
//! ### Vector Extraction Trick
//! For arrays like thermal_effects and reaction parameters, uses `.as_vector().unwrap().clone()`
//! to extract Vec<f64> directly from the parser's Value enum.
//!
//! ### Bounds Tuple Construction
//! Bounds parsing extracts 2-element vectors and converts them to tuples using array indexing:
//! `(bounds_vec[0], bounds_vec[1])` for efficient tuple creation.
//!
//! ### Error Handling Strategy
//! Uses `.expect()` for required parameters and `if let Some()` for optional ones,
//! providing clear error messages for missing required configuration.

use super::SimpleReactorBVP::{FastElemReact, SimpleReactorTask};
use crate::ReactorsBVP::reactor_BVP_utils::{BoundsConfig, ScalingConfig, ToleranceConfig};
use RustedSciThe::Utils::task_parser::{DocumentMap, DocumentParser};
use RustedSciThe::numerical::BVP_Damp::NR_Damp_solver_damped::{NRBVP, SolverParams};
use nalgebra::DMatrix;
use std::collections::HashMap;
impl SimpleReactorTask {
    /// Parses reactor parameters from a DocumentMap and populates the SimpleReactorTask structure.
    ///
    /// This is the core parsing method that extracts all reactor configuration from the parsed
    /// document structure. It handles reactions, process conditions, boundary conditions,
    /// diffusion coefficients, and chemical groups.
    ///
    /// # Arguments
    /// * `task_hashmap` - The parsed document structure containing all configuration sections
    ///
    /// # Sections Parsed
    /// - `reactions`: Arrhenius parameters (A, n, E, Q) for each reaction
    /// - `process_conditions`: Physical properties (Tm, P, Cp, Lambda, m, etc.)
    /// - `boundary_condition`: Initial/boundary values for all variables
    /// - `diffusion_coefficients`: Transport properties for each substance
    /// - Chemical groups: Atomic composition for custom substance definitions
    pub fn set_reactor_params_from_hashmap(&mut self, task_hashmap: DocumentMap) {
        // reactions
        let reactions = task_hashmap
            .get("reactions")
            .expect("Failed to get solver settings");

        let reactions: HashMap<String, Vec<f64>> = reactions
            .iter()
            .map(|(key, value)| {
                let binding = value.clone().unwrap();
                let value = binding[0]
                    .as_vector()
                    .unwrap()
                    .clone()
                    .into_iter()
                    .map(|val| val)
                    .collect::<Vec<f64>>();

                (key.to_owned(), value)
            })
            .collect();
        let mut vec_of_struct: Vec<FastElemReact> = Vec::new();
        for (key, value) in reactions {
            let A = value[0];
            let n = value[1];
            let E = value[2];
            let Q = value[3];
            let react = FastElemReact {
                eq: key,
                A,
                n,
                E,
                Q,
            };

            vec_of_struct.push(react);
        }

        let _ = self.fast_react_set(vec_of_struct);

        //
        let process_conditions = task_hashmap
            .get("process_conditions")
            .expect("Failed to get solver settings");
        // String and f64 params
        let problem_name = if let Some(problem_name_val) = process_conditions.get("problem_name") {
            problem_name_val.clone().unwrap()[0]
                .as_option_string()
                .cloned()
        } else {
            None
        };
        // dbg!(&roblem_name_val);
        self.problem_name = problem_name;
        let problem_description =
            if let Some(problem_description) = process_conditions.get("problem_description") {
                problem_description.clone().unwrap()[0]
                    .as_option_string()
                    .cloned()
            } else {
                None
            };
        let substances: Vec<String> = process_conditions
            .get("substances")
            .expect("Failed to get substances")
            .clone()
            .unwrap()
            .iter()
            .map(|val| val.as_string().unwrap().clone())
            .collect();

        self.kindata.substances = substances;

        self.problem_description = problem_description;
        let Tm = process_conditions
            .get("Tm")
            .expect("Failed to get Tm")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        self.Tm = Tm;
        let L = process_conditions
            .get("L")
            .expect("Failed to get L")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let dT = process_conditions
            .get("dT")
            .expect("Failed to get dT")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let T_scale = process_conditions
            .get("T_scale")
            .expect("Failed to get T_scale")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let scaling = ScalingConfig::new(dT, L, T_scale);
        self.scaling = scaling;
        let P = process_conditions
            .get("P")
            .expect("Failed to get P")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        self.P = P;
        let Cp = process_conditions
            .get("Cp")
            .expect("Failed to get Cp")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        self.Cp = Cp;
        let Lambda = process_conditions
            .get("Lambda")
            .expect("Failed to get Lambda")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        self.Lambda = Lambda;
        let m = process_conditions
            .get("m")
            .expect("Failed to get m")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        self.m = m;
        if let Some(M) = process_conditions.get("M") {
            let M = M.clone().unwrap()[0].as_float().unwrap();
            self.M = M;
        }
        // parsing vectors and hashmaps
        let thermal_effects = process_conditions
            .get("thermal_effects")
            .expect("Failed to get thermal_effects")
            .clone()
            .unwrap()[0]
            .as_vector()
            .unwrap()
            .clone();

        self.thermal_effects = thermal_effects;
        let boundary_condition = task_hashmap
            .get("boundary_condition")
            .expect("Failed to get solver settings");
        let boundary_condition: HashMap<String, f64> = boundary_condition
            .iter()
            .map(|(key, value)| {
                let binding = value.clone().unwrap();

                let value = binding[0].as_float().unwrap();

                (key.to_owned(), value)
            })
            .collect();
        self.boundary_condition = boundary_condition;
        let diffusion_coefficients = task_hashmap
            .get("diffusion_coefficients")
            .expect("Failed to get solver settings");
        let diffusion_coefficients: HashMap<String, f64> = diffusion_coefficients
            .iter()
            .map(|(key, value)| {
                let binding = value.clone().unwrap();
                let value = binding[0].as_float().unwrap();

                (key.to_owned(), value)
            })
            .collect();
        self.Diffusion = diffusion_coefficients;
        // groups
        if let Some(groups_val) = process_conditions.get("groups") {
            let groups_enabled = groups_val.clone().unwrap()[0].as_boolean().unwrap();
            if groups_enabled {
                let mut groups: HashMap<String, HashMap<String, usize>> = HashMap::new();
                for sub_i in self.kindata.substances.iter() {
                    if let Some(sub_i_groups) = task_hashmap.get(sub_i.as_str()) {
                        let sub_i_groups_hashmap: HashMap<String, usize> = sub_i_groups
                            .iter()
                            .map(|(key, value)| {
                                let binding = value.clone().unwrap();
                                let value = binding[0].as_usize().unwrap();
                                (key.to_owned(), value)
                            })
                            .collect();
                        groups.insert(sub_i.clone(), sub_i_groups_hashmap);
                    }
                }
                self.kindata.groups = Some(groups);
            }
        }
    }
    /// Parses solver tolerance and bounds configurations from the document.
    ///
    /// Extracts numerical solver settings including relative tolerances and variable bounds
    /// for the four main variable types: C (concentrations), J (fluxes), Teta (temperature), q (heat flux).
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the parsed configuration
    ///
    /// # Returns
    /// * `Ok((ToleranceConfig, BoundsConfig))` - Parsed tolerance and bounds configurations
    /// * `Err(String)` - Error message if required sections are missing
    ///
    /// # Expected Document Structure
    /// ```text
    /// rel_tolerance
    /// C: 1e-5
    /// J: 1e-5
    /// Teta: 1e-5
    /// q: 1e-5
    /// bounds
    /// C: -10.0, 10.0
    /// J: -1e20, 1e20
    /// Teta: -100.0, 100.0
    /// q: -1e20, 1e20
    /// ```
    fn parse_toleranse_and_bounds(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(ToleranceConfig, BoundsConfig), String> {
        let task_hashmap = parser.get_result().unwrap();
        // Parse tolerance configuration
        let rel_tolerance = task_hashmap.get("rel_tolerance").unwrap();
        let C = rel_tolerance
            .get("C")
            .ok_or("Missing C in rel_tolerance")?
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let J = rel_tolerance
            .get("J")
            .ok_or("Missing J in rel_tolerance")?
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let Teta = rel_tolerance
            .get("Teta")
            .ok_or("Missing Teta in rel_tolerance")?
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let q = rel_tolerance
            .get("q")
            .ok_or("Missing q in rel_tolerance")?
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let tolerance_config = ToleranceConfig::new(C, J, Teta, q);

        // Parse bounds configuration
        let bounds = task_hashmap.get("bounds").unwrap();
        let C_bounds = bounds
            .get("C")
            .expect("Missing C in bounds")
            .clone()
            .unwrap();

        let C = (
            C_bounds[0].as_float().unwrap(),
            C_bounds[1].as_float().unwrap(),
        );

        let J_bounds = bounds
            .get("J")
            .expect("Missing J in bounds")
            .clone()
            .unwrap();
        let J = (
            J_bounds[0].as_float().unwrap(),
            J_bounds[1].as_float().unwrap(),
        );

        let Teta_bounds = bounds
            .get("Teta")
            .expect("Missing Teta in bounds")
            .clone()
            .unwrap();
        let Teta = (
            Teta_bounds[0].as_float().unwrap(),
            Teta_bounds[1].as_float().unwrap(),
        );

        let q_bounds = bounds
            .get("q")
            .expect("Missing q in bounds")
            .clone()
            .unwrap();
        let q = (
            q_bounds[0].as_float().unwrap(),
            q_bounds[1].as_float().unwrap(),
        );

        let bounds_config = BoundsConfig::new(C, J, Teta, q);

        Ok((tolerance_config, bounds_config))
    }

    /// Extracts basic solver settings (time domain and grid size) from the parsed document.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the configuration
    ///
    /// # Returns
    /// * `(t0, t_end, n_steps)` - Start time, end time, and number of grid steps
    pub fn parse_basic_settings(&mut self, parser: DocumentParser) -> (f64, f64, usize, String) {
        let task_hashmap = parser.get_result().unwrap();
        let process_conditions = task_hashmap
            .get("process_conditions")
            .expect("Failed to get solver settings");
        let t0 = process_conditions
            .get("t0")
            .expect("Failed to get t0")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let t_end = process_conditions
            .get("t_end")
            .expect("Failed to get t_end")
            .clone()
            .unwrap()[0]
            .as_float()
            .unwrap();
        let n_steps = process_conditions
            .get("n_steps")
            .expect("Failed to get n_steps")
            .clone()
            .unwrap()[0]
            .as_usize()
            .unwrap();
        let arg = if let Some(arg) = process_conditions.get("arg") {
            arg.clone().unwrap()[0].as_string().unwrap().clone()
        } else {
            "x".to_string()
        };
        (t0, t_end, n_steps, arg)
    }

    pub fn set_postpocessing_from_hashmap(&mut self, parser: &mut DocumentParser) {
        let result: DocumentMap = parser.get_result().unwrap().clone();
        let solver_settings = result
            .get("postprocessing")
            .expect("Failed to get postpocessing");
        // flag to create plot via rust native crate
        let plot_flag = if let Some(plot) = solver_settings.get("plot") {
            plot.clone().unwrap()[0]
                .as_boolean()
                .expect("Failed to get plot as bool")
        } else {
            false
        };
        // flag to create plot via GNU plot library
        let gnuplot_flag = if let Some(gnuplotflag) = solver_settings.get("gnuplot") {
            gnuplotflag.clone().unwrap()[0]
                .as_boolean()
                .expect("Failed to get gnuplot as bool")
        } else {
            false
        };
        // flag to save solution to txt
        let save_flag = if let Some(save) = solver_settings.get("save") {
            save.clone().unwrap()[0]
                .as_boolean()
                .expect("Failed to get save as bool")
        } else {
            false
        };
        let save_to_csv = if let Some(save) = solver_settings.get("save_to_csv") {
            save.clone().unwrap()[0]
                .as_boolean()
                .expect("Failed to get save_to_csv as bool")
        } else {
            false
        };

        let name = if let Some(name) = solver_settings.get("filename") {
            name.clone().unwrap()[0].as_string().cloned()
        } else {
            None
        };
        // return from dimensionless to dimensioned unknowns
        let return_to_dimension = if let Some(return_to_dimension) = solver_settings.get("return_to_dimension") {
            return_to_dimension.clone().unwrap()[0].as_boolean().unwrap()
        } else {
            true
        };
        if return_to_dimension {
            self.postprocessing();
        }
        if plot_flag {
            self.plot();
        }
        if gnuplot_flag {
            self.gnuplot();
        };
        if save_flag {
            self.save_to_file(name.clone())
        };
        if save_to_csv {
            self.save_to_csv(name);
        };

    }

    /// Complete workflow: parses configuration, sets up BVP, configures solver, and solves the problem.
    ///
    /// This is a high-level method that orchestrates the entire process from parsed document
    /// to solved BVP. It handles all intermediate steps including BVP setup, solver configuration,
    /// tolerance/bounds parsing, and postprocessing.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser containing the complete problem configuration
    ///
    /// # Process Flow
    /// 1. Parse reactor parameters
    /// 2. Setup BVP equations and boundary conditions
    /// 3. Configure NRBVP solver with parsed settings
    /// 4. Parse tolerances and bounds
    /// 5. Set initial guess and solve
    /// 6. Apply postprocessing (plotting, data export)
    pub fn solve_from_map(&mut self, mut parser: DocumentParser) {
        self.parse_parameters_with_exact_names(&mut parser)
            .expect("Failed to parse parameters");
        let _ = self.setup_bvp();
        assert!(!self.solver.unknowns.is_empty());
        assert!(!self.solver.eq_system.is_empty());

        let mut nr = NRBVP::default();
        let _ = nr.parse_settings(&mut parser);
        let solver = &self.solver;
        nr.eq_system = solver.eq_system.clone();
        nr.values = solver.unknowns.clone();
        let BC = solver.BorderConditions.clone();
        let BC: HashMap<String, Vec<(usize, f64)>> =
            BC.iter().map(|(k, v)| (k.clone(), vec![*v])).collect();
        nr.BorderConditions = BC;
        let (tolerance_config, bounds_config) =
            self.parse_toleranse_and_bounds(&mut parser).unwrap();
        let rel_tolerance = Some(tolerance_config.to_full_tolerance_map(&self.kindata.substances));
        let bounds = Some(bounds_config.to_full_bounds_map(&self.kindata.substances));
        nr.rel_tolerance = rel_tolerance;
        nr.Bounds = bounds;
        let (t0, t_end, n_steps, arg) = self.parse_basic_settings(parser.clone());
        nr.t0 = t0;
        nr.t_end = t_end;
        nr.arg = arg.clone();
        self.solver.arg_name = arg;
        let ig = vec![1e-2; n_steps * self.solver.unknowns.len()];
        let initial_guess = DMatrix::from_vec(self.solver.unknowns.len(), n_steps, ig);
        nr.initial_guess = initial_guess;
        nr.n_steps = n_steps;
        nr.before_solve_preprocessing();
        nr.solve();
        self.solver.x_mesh = Some(nr.x_mesh.clone());
        let y_result = nr.get_result();
        self.solver.solution = y_result;
        self.set_postpocessing_from_hashmap(&mut parser);
    }

    /// One-shot method: reads configuration file, parses it, and solves the BVP problem.
    ///
    /// Convenience method that combines file reading, parsing, and solving in a single call.
    /// Ideal for batch processing or when the entire workflow can be automated.
    ///
    /// # Arguments
    /// * `path` - Path to the configuration file
    pub fn solve_from_file(&mut self, path: std::path::PathBuf) {
        let mut reactor = SimpleReactorTask::new();
        // Test parsing
        let mut parser = reactor
            .parse_file(Some(path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        reactor.solve_from_map(parser);
    }

    /// Creates a DocumentParser from a file path.
    ///
    /// Simple wrapper around DocumentParser creation that handles file I/O.
    ///
    /// # Arguments
    /// * `path` - Optional path to configuration file
    ///
    /// # Returns
    /// * `Ok(DocumentParser)` - Configured parser ready for document parsing
    /// * `Err(String)` - File I/O error message
    pub fn parse_file(
        &mut self,
        path: Option<std::path::PathBuf>,
    ) -> Result<DocumentParser, String> {
        let mut parser = DocumentParser::new(String::new());
        parser.setting_from_file(path)?;
        Ok(parser)
    }

    /// High-level parsing orchestrator that processes the document and populates reactor parameters.
    ///
    /// Coordinates the document parsing process and delegates to `set_reactor_params_from_hashmap`
    /// for actual parameter extraction.
    ///
    /// # Arguments
    /// * `parser` - DocumentParser with loaded configuration
    ///
    /// # Returns
    /// * `Ok(())` - Successful parsing
    /// * `Err(String)` - Parsing error message
    pub fn parse_parameters_with_exact_names(
        &mut self,
        parser: &mut DocumentParser,
    ) -> Result<(), String> {
        let _ = parser.parse_document();

        let result: DocumentMap = parser
            .get_result()
            .ok_or("No result after parsing")?
            .clone();
        self.set_reactor_params_from_hashmap(result);
        Ok(())
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a template configuration file for reactor BVP tasks.
///
/// Generates a commented template file that users can fill with their specific values.
/// The template includes all required sections with example values and explanatory comments.
pub fn create_template() {
    use std::fs::File;
    use std::io::Write;

    let template_content = r#"# Reactor BVP Configuration Template
# Fill in the values below for your specific problem

# Process conditions - main problem parameters
process_conditions
# Optional problem identification
problem_name: Some(YourProblemName)
problem_description: Some(YourProblemDescription)
# List of chemical substances (comma-separated)
substances: Substance1, Substance2
# Start position/time
t0: 0.0
# End position/time
t_end: 1.0
# Number of grid points
n_steps: 200
# Independent variable name (x, t, etc.)
arg: x
# Maximum temperature [K]
Tm: 1500.0
# Characteristic length [m]
L: 9e-4
# Temperature difference [K]
dT: 600.0
# Temperature scaling [K]
T_scale: 600.0
# Pressure [Pa]
P: 1e6
# Heat capacity [J/kg/K]
Cp: 1464.4
# Thermal conductivity [W/m/K]
Lambda: 0.07
# Mass parameter [kg]
m: 0.0043
# Molar mass [kg/mol]
M: 0.0342
# Thermal effects for each reaction [J/mol]
thermal_effects: [102000.0]
# Enable custom atomic groups (true/false)
groups: true

# Boundary conditions - initial values for all variables
boundary_condition
# Initial concentration/mole fraction
Substance1: 0.999
# Initial concentration/mole fraction
Substance2: 0.001
# Initial temperature [K]
T: 800.0

# Diffusion coefficients for each substance [m²/s]
diffusion_coefficients
Substance1: 0.000009296
Substance2: 0.000009296

# Atomic composition (only if groups: true)
# Define atomic composition for each substance
Substance1
# Number of hydrogen atoms
H: 4
# Number of nitrogen atoms
N: 8
# Number of carbon atoms
C: 8
# Number of oxygen atoms
O: 8

Substance2
H: 6
C: 1
O: 1

# Chemical reactions with Arrhenius parameters
# Format: Reactant=>Product: [A, n, E, Q]
# A: pre-exponential factor, n: temperature exponent
# E: activation energy [J/mol], Q: heat of reaction [J/mol]
reactions
Substance1=>Substance2: [130000.0, 0.0, 20920.0, 102000.0]

# Numerical solver settings
solver_settings
# Discretization scheme
scheme: forward
# Solution method
method: Sparse
# Convergence strategy
strategy: Damped
# Linear system solver
linear_sys_method: None
# Absolute tolerance
abs_tolerance: 1e-5
# Maximum iterations
max_iterations: 100
# Logging level
loglevel: Some(info)
# Disable log saving
dont_save_logs: true

# Variable bounds for solver
bounds
# Concentration bounds
C: -10.0, 10.0
# Flux bounds
J: -1e20, 1e20
# Temperature bounds
Teta: -100.0, 100.0
# Heat flux bounds
q: -1e20, 1e20

# Relative tolerances for each variable type
rel_tolerance
# Concentration tolerance
C: 1e-5
# Flux tolerance
J: 1e-5
# Temperature tolerance
Teta: 1e-5
# Heat flux tolerance
q: 1e-5

# Advanced solver strategy parameters
strategy_params
# Maximum Jacobian updates
max_jac: Some(3)
# Maximum damping iterations
max_damp_iter: Some(10)
# Damping factor
damp_factor: Some(0.5)
# Adaptive strategy
adaptive: None

# Output and visualization options
postprocessing
# Generate gnuplot output
gnuplot: true
# Save results to CSV
save_to_csv: false
# Output filename prefix
filename: output_name
"#;

    let mut file = File::create("template.txt").expect("Failed to create template.txt");
    file.write_all(template_content.as_bytes())
        .expect("Failed to write template");
    println!("Template created: template.txt");
}
////////////////////////////////////////////////////TESTS///////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;
    const task_content: &str = r#"
        process_conditions
        problem_name: Some(HMXTest)
        problem_description: Some(HMXdecompositiontest)
        substances: HMX, HMXprod
        t0: 0.0
        t_end: 1.0
        n_steps: 200
        arg:x
        Tm: 1500.0
        L: 9e-4
        dT: 600.0
        T_scale: 600.0
        P: 1e6
        Cp: 1464.4
        Lambda: 0.07
        m: 0.0043
        M: 0.0342
        thermal_effects: [102000.0]
        groups:true
        boundary_condition
        HMX: 0.999
        HMXprod: 0.001
        T: 800.0
        diffusion_coefficients
        HMX: 0.000009296
        HMXprod: 0.000009296
        HMX
        H: 4
        N: 8
        C: 8
        O: 8
        HMXprod
        H: 6
        C: 1
        O: 1
        reactions
        HMX=>HMXprod: [130000.0, 0.0, 20920.0, 102000.0]
        solver_settings
        scheme: forward
        method: Sparse
        strategy: Damped
        linear_sys_method: None
        abs_tolerance: 1e-5
        max_iterations: 100
        loglevel: Some(info)
        dont_save_logs: true
        bounds
        C: -10.0, 10.0
        J:  -1e20, 1e20
        Teta:-100.0, 100.0
        q: -1e20, 1e20
        rel_tolerance
        C: 1e-5
        J: 1e-5
        Teta: 1e-5
        q:  1e-5
        strategy_params
        max_jac: Some(3)
        max_damp_iter: Some(10)
        damp_factor: Some(0.5)
        adaptive: None
        postprocessing
        gnuplot:true
        save_to_csv:false
        filename: meow
        "#;
    #[test]
    fn parse_directly() {
        let mut parser = DocumentParser::new(task_content.to_string());
        let result: DocumentMap = parser.parse_document().unwrap().clone();
        println!("result {:?}", result);
    }

    #[test]
    fn parsng_task_elementary() {
        let mut reactor = SimpleReactorTask::new();

        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");

        // HMX data from create_hmx function

        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        println!("parser {:?}", parser);

        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("Failed to parse parameters");

        // Verify parsed data
        assert_eq!(reactor.problem_name, Some("HMXTest".to_string()));
        assert_eq!(
            reactor.kindata.substances,
            vec!["HMX".to_string(), "HMXprod".to_string()]
        );

        assert_eq!(reactor.Tm, 1500.0);
        assert_eq!(reactor.P, 1e6);
        assert_eq!(reactor.Cp, 1464.4);
        assert_eq!(reactor.Lambda, 0.07);
        assert_eq!(reactor.m, 0.0043);
        assert_eq!(reactor.M, 0.0342);
        assert_eq!(reactor.thermal_effects, vec![102000.0]);

        // Verify boundary conditions
        assert_eq!(reactor.boundary_condition.get("HMX"), Some(&0.999));
        assert_eq!(reactor.boundary_condition.get("HMXprod"), Some(&0.001));
        assert_eq!(reactor.boundary_condition.get("T"), Some(&800.0));

        // Verify diffusion coefficients
        assert_eq!(reactor.Diffusion.get("HMX"), Some(&0.000009296));
        assert_eq!(reactor.Diffusion.get("HMXprod"), Some(&0.000009296));

        // Verify groups
        assert!(reactor.kindata.groups.is_some());
        let groups = reactor.kindata.groups.as_ref().unwrap();
        assert_eq!(groups.get("HMX").unwrap().get("H"), Some(&4));
        assert_eq!(groups.get("HMX").unwrap().get("N"), Some(&8));
        assert_eq!(groups.get("HMXprod").unwrap().get("H"), Some(&6));
        assert_eq!(groups.get("HMXprod").unwrap().get("C"), Some(&1));

        // Verify reactions were parsed
        assert_eq!(reactor.kindata.vec_of_equations.len(), 1);
        assert_eq!(reactor.kindata.vec_of_equations[0], "HMX=>HMXprod");
    }

    #[test]
    fn test_parse_toleranse_and_bounds() {
        let mut reactor = SimpleReactorTask::new();

        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("test_bounds.txt");

        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Parse file
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();

        // Test parsing tolerances and bounds
        let (tolerance_config, bounds_config) = reactor
            .parse_toleranse_and_bounds(&mut parser)
            .expect("Failed to parse tolerances and bounds");

        // Verify tolerance config
        assert_eq!(tolerance_config.C, 1e-5);
        assert_eq!(tolerance_config.J, 1e-5);
        assert_eq!(tolerance_config.Teta, 1e-5);
        assert_eq!(tolerance_config.q, 1e-5);

        // Verify bounds config
        assert_eq!(bounds_config.C, (-10.0, 10.0));
        assert_eq!(bounds_config.J, (-1e20, 1e20));
        assert_eq!(bounds_config.Teta, (-100.0, 100.0));
        assert_eq!(bounds_config.q, (-1e20, 1e20));
    }

    #[test]
    fn test_setup_bvp_from_file() {
        let mut reactor = SimpleReactorTask::new();
        // Create temporary directory
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");
        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        // println!("parser {:?}", parser);

        reactor
            .parse_parameters_with_exact_names(&mut parser)
            .expect("Failed to parse parameters");

        let res = reactor.setup_bvp();
        match res {
            Ok(_) => println!("ok"),
            Err(e) => println!("error {:?}", e),
        }
        let rates = reactor.map_eq_rate.clone();
        for (eq, rate) in rates {
            println!("reaction {} rate {}", eq, rate);
        }
        println!("\n \n");
        let system = reactor.map_of_equations.clone();
        for (subs, (variable, eq)) in system {
            println!("subs: {} | variable: {} | eq: {} | \n", subs, variable, eq);
        }
        let bc = &reactor.solver.BorderConditions;
        println!("bc {:?}", bc);
        println!(" unknowns{:?}", reactor.solver.unknowns);
    }

    #[test]
    fn test_solve_from_map() {
        let mut reactor = SimpleReactorTask::new();
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let file_path = temp_dir.path().join("hmx_task.txt");
        // Write to temporary file
        let mut file = File::create(&file_path).expect("Failed to create temp file");
        file.write_all(task_content.as_bytes())
            .expect("Failed to write to temp file");

        // Test parsing
        let mut parser = reactor
            .parse_file(Some(file_path))
            .expect("Failed to parse file");
        let _ = parser.parse_document();
        reactor.solve_from_map(parser);
        // println!("parser {:?}", parser);
    }
}
