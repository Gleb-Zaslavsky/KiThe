//! # BVP Equation Creation Module
//!
//! This module contains the core mathematical engine for creating boundary value problem equations
//! from reactor physics and chemistry. It transforms dimensional reactor equations into dimensionless
//! symbolic form suitable for numerical solution.
//!
//! ## Key Mathematical Transformations
//!
//! ### Energy Balance
//! - **Dimensional**: d(λ∇T)/dx - ṁCp(dT/dx) + Q = 0
//! - **Dimensionless**: dΘ/dz = q/λ, dq/dz = Pe_q·q - (L²/dT)·∑(Qⱼ·Rⱼ)
//!
//! ### Mass Balance  
//! - **Dimensional**: d(D·ρ∇Cᵢ)/dx - ṁ(dCᵢ/dx) + Wᵢ = 0
//! - **Dimensionless**: dCᵢ/dz = Jᵢ/(D·ρ), dJᵢ/dz = Pe_D·Jᵢ - L²·Wᵢ
//!
//! Where:
//! - z = x/L (dimensionless coordinate)
//! - Θ = (T-T₀)/dT (dimensionless temperature)
//! - Pe_q = L·ṁ·Cp/λ (thermal Peclet number)
//! - Pe_D = ṁ·L/(D·ρ) (mass Peclet number)
//! - Rⱼ = kⱼ·∏[Cᵢ]^νᵢ (reaction rates)
//! - Wᵢ = ∑(νᵢⱼ·Rⱼ·Mᵢ) (mass production rates)
//!
//! ## Interesting Features
//!
//! - **Automatic symbolic differentiation**: Creates exact analytical Jacobians
//! - **Dimensionless scaling**: Ensures numerical stability across different scales
//! - **Flexible reaction kinetics**: Supports arbitrary Arrhenius expressions
//! - **Mass-weighted rates**: Properly accounts for molecular weight differences
//! - **Boundary condition mapping**: Automatically maps physical BCs to mathematical variables
use super::SimpleReactorIVP::SimpleReactorTask;
use crate::ReactorsBVP::SimpleReactorBVP::ReactorError;
use RustedSciThe::symbolic::symbolic_engine::Expr;
use std::collections::HashMap;

impl SimpleReactorTask {
    /////////////////////////CREATING SYMBOLIC EQUATIONS///////////////////////////////////

    /// Create the complete system of dimensionless BVP equations
    ///
    /// This is the most critical method in the entire reactor modeling framework.
    /// It transforms the dimensional reactor physics into a mathematically tractable
    /// system of ordinary differential equations.
    ///
    /// ## Mathematical Framework
    ///
    /// The method creates a system of 2n+2 equations for n substances:
    /// 1. **Energy equations** (2): dΘ/dz, dq/dz
    /// 2. **Mass equations** (2n): dCᵢ/dz, dJᵢ/dz for each substance i
    ///
    /// ## Algorithm Steps
    ///
    /// 1. **Validation**: Check dimensions and data consistency
    /// 2. **Rate Expression Assembly**: Build symbolic reaction rates Rⱼ = kⱼ·∏[Cᵢ]^νᵢⱼ
    /// 3. **Energy Balance Construction**:
    ///    - dΘ/dz = q/λ (Fourier's law)
    ///    - dq/dz = Pe_q·q - (L²/dT)·∑(Qⱼ·Rⱼ) (energy conservation)
    /// 4. **Mass Balance Construction**:
    ///    - dCᵢ/dz = Jᵢ/(D·ρ) (Fick's law)
    ///    - dJᵢ/dz = Pe_D·Jᵢ - L²·Wᵢ (mass conservation)
    /// 5. **Production Rate Calculation**: Wᵢ = ∑(νᵢⱼ·Rⱼ·Mᵢ)
    ///
    /// ## Key Insights
    ///
    /// - **Dimensionless scaling** prevents numerical issues from vastly different scales
    /// - **Symbolic expressions** enable exact Jacobian computation for Newton methods
    /// - **Mass-weighted rates** ensure proper stoichiometric balance
    /// - **Peclet numbers** capture the physics of convection vs diffusion/conduction
    ///
    /// ## Returns
    ///
    /// Populates:
    /// - `self.solver.eq_system`: Vector of RHS expressions
    /// - `self.solver.unknowns`: Variable names ["Teta", "q", "C0", "J0", "C1", "J1", ...]
    /// - `self.map_of_equations`: Mapping from physical quantities to mathematical variables
    /// - `self.map_eq_rate`: Reaction rate expressions for each reaction
    pub fn create_IVP_equations(&mut self) -> Result<(), ReactorError> {
        let substances = &self.kindata.substances;
        let n = substances.len();

        if n == 0 {
            return Err(ReactorError::MissingData(
                "No substances found in kinetic data".to_string(),
            ));
        }

        let equations = self.kindata.vec_of_equations.clone();
        let k = equations.len();

        if k == 0 {
            return Err(ReactorError::MissingData(
                "No reactions found in kinetic data".to_string(),
            ));
        }

        // Extract scaling parameters for dimensionless transformation
        let T_scale = self.scaling.T_scale;
        // Scaling factors for energy and mass equations
        let qm = Expr::Const(self.L.powi(2) / T_scale); // Energy source term scaling: L²/dT
        let qc = Expr::Const(self.L.powi(2)); // Mass source term scaling: L²

        let (Ci_expr, Ci) = Expr::IndexedVars(n, "C");

        let Ci: Vec<String> = Ci.iter().map(|var| var.replace("_", "")).collect();

        let k_sym_vec = self.kindata.K_sym_vec.as_ref().ok_or_else(|| {
            ReactorError::MissingData("Symbolic rate constants not calculated".to_string())
        })?;

        if k_sym_vec.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Mismatch between number of reactions and rate constants".to_string(),
            ));
        }
        // vector pf vectors for powers of reagents in each reaction
        let G_react = &self.kindata.stecheodata.stecheo_reags;
        let stoich_matrix = &self.kindata.stecheodata.stecheo_matrx;

        if stoich_matrix.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Stoichiometric matrix size mismatch".to_string(),
            ));
        }

        let Q: Vec<f64> = self.thermal_effects.clone();
        if Q.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Thermal effects vector size mismatch".to_string(),
            ));
        }

        let Pe_q = self.Pe_q;

        let Mi = self.kindata.stecheodata.vec_of_molmasses.clone().unwrap();
        let Mi: Vec<f64> = Mi.iter().map(|m| *m / 1000.0).collect();

        let ro_m = Expr::Const(self.ro);
        let mut map_eq_rate: HashMap<String, Expr> = HashMap::new();
        let mut Rates: Vec<Expr> = Vec::new();

        for j in 0..k {
            let K_j = k_sym_vec[j].clone();
            let mut rate_expr = K_j.clone();
            //  println!("rate expr {}", &rate_expr);
            for i in 0..n {
                let powers_ji = G_react.get(j).and_then(|row| row.get(i)).ok_or_else(|| {
                    ReactorError::IndexOutOfBounds(format!("G_react[{}][{}] out of bounds", j, i))
                })?;
                let Powers_ji = Expr::Const(*powers_ji);
                let ci_expr = Ci_expr.get(i).ok_or_else(|| {
                    ReactorError::IndexOutOfBounds(format!("Ci_expr[{}] out of bounds", i))
                })?;
                let Mi = Expr::Const(Mi[i].clone());
                rate_expr = rate_expr * (ro_m.clone() * ci_expr.clone() / Mi).pow(Powers_ji);
            }
            rate_expr = rate_expr.simplify_();
            map_eq_rate.insert(equations[j].clone(), rate_expr.clone());
            Rates.push(rate_expr);
        }

        self.map_eq_rate = map_eq_rate;
        let mut vec_of_unknowns: Vec<String> = Vec::new();
        let mut vec_of_equations = Vec::new();
        let mut map_of_equations: HashMap<String, (String, Expr)> = HashMap::new();
        // creating RHS for energy balance equations
        // unknown q (heat fux) in the LHS = dq/dx
        let q = Expr::Var("q".to_string());
        // unknown Teta in the left side = dTeta/dx...
        vec_of_unknowns.push("Teta".to_string());
        // is equal to RHS q/Lambda
        let lambda = Expr::Const(self.Lambda);

        let RHS_Teta = q.clone() / lambda;
        vec_of_equations.push(RHS_Teta.clone());
        map_of_equations.insert("Teta".to_string(), ("Teta".to_string(), RHS_Teta));

        vec_of_unknowns.push("q".to_string());
        let mut Heat_prod = Expr::Const(0.0);
        for j in 0..k {
            let Qj = Expr::Const(Q[j]); //  
            Heat_prod = (Heat_prod + Qj * Rates[j].clone()).simplify_();
        }
        self.heat_release = Heat_prod.clone();
        let Pe_q_expr = Expr::Const(Pe_q);
        let RHS_q = q * Pe_q_expr - Heat_prod * qm;
        vec_of_equations.push(RHS_q.clone());
        map_of_equations.insert("q".to_string(), ("q".to_string(), RHS_q));
        // STEP 3: Create mass balance equations for each substance
        // For each substance i: dCᵢ/dz = - L²·Wᵢ
        for i in 0..n {
            let ci = Ci.get(i).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!("Ci[{}] out of bounds", i))
            })?;

            // First mass equation: dCᵢ/dz = Jᵢ/(D·ρ) (Fick's law in dimensionless form)
            vec_of_unknowns.push(ci.clone());

            // Calculate mass production rate: Wᵢ = ∑ⱼ(νᵢⱼ·Rⱼ·Mᵢ)
            let mut Wi = Expr::Const(0.0); // Mass production rate for substance i
            for j in 0..k {
                let stoich_coeff_val =
                    stoich_matrix
                        .get(j)
                        .and_then(|row| row.get(i))
                        .ok_or_else(|| {
                            ReactorError::IndexOutOfBounds(format!(
                                "stoich_matrix[{}][{}] out of bounds",
                                j, i
                            ))
                        })?; // Stoichiometric coefficient νᵢⱼ

                let stoich_coeff = Expr::Const(*stoich_coeff_val);
                let rate_j = Rates[j].clone(); // Reaction rate Rⱼ
                let Mi = Expr::Const(Mi[i]); // Molecular weight Mᵢ

                // Add contribution from reaction j: νᵢⱼ·Rⱼ·Mᵢ
                Wi = Wi + rate_j * stoich_coeff * Mi;
            }

            // Second mass equation: dJᵢ/dz = Pe_D·Jᵢ - L²·Wᵢ (mass conservation)

            let RHS_i = (-Wi * qc.clone()).simplify_(); // Convection - production
            vec_of_equations.push(RHS_i.clone());
            map_of_equations.insert(format!("{}_flux", substances[i]), (ci.clone(), RHS_i));
        }

        self.solver.eq_system = vec_of_equations;
        self.solver.unknowns = vec_of_unknowns;
        self.map_of_equations = map_of_equations;
        Ok(())
    }
    /////////////////////////////////CREATING BC///////////////////////////////////////////////////////////////

    /// Set boundary conditions for the BVP solver
    ///
    /// Maps physical boundary conditions to mathematical variables and determines
    /// whether each condition applies at the left (inlet, z=0) or right (outlet, z=1) boundary.
    ///
    /// ## Boundary Condition Logic
    ///
    /// - **Concentrations (Cᵢ)**: Set at inlet (left boundary, z=0)
    /// - **Fluxes (Jᵢ)**: Set at outlet (right boundary, z=1) - typically zero
    /// - **Temperature (Θ)**: Set at inlet (left boundary, z=0)
    /// - **Heat flux (q)**: Set at outlet (right boundary, z=1) - typically zero
    ///
    /// ## Temperature Scaling
    ///
    /// Physical temperature T is converted to dimensionless form:
    /// Θ = (T - dT)/dT, where dT is the temperature scaling parameter
    pub fn set_solver_BC(&mut self) -> Result<(), ReactorError> {
        let border_conditions: HashMap<String, f64> = self.boundary_condition.clone();
        let map_of_equstions: HashMap<String, (String, Expr)> = self.map_of_equations.clone();

        // Boundary condition format: {variable → (boundary_index, value)}
        // boundary_index: 0 = left boundary (inlet, z=0), 1 = right boundary (outlet, z=1)
        // variable: solver variable name ("Teta", "q", "C0", "J0", etc.)
        let mut BorderConditions: HashMap<String, (usize, f64)> = HashMap::new();

        // Process substance boundary conditions
        // Concentrations → left boundary (inlet conditions)
        // Fluxes → right boundary (outlet conditions, typically zero)
        for (key, (variable, _)) in map_of_equstions.iter() {
            if key.ends_with("_flux") {
                // Flux variable: set at outlet (right boundary)
                let substance = key.strip_suffix("_flux").unwrap();
                match border_conditions.get(&format!("{}_flux", substance)) {
                    Some(condition_value) => {
                        BorderConditions.insert(variable.clone(), (1, *condition_value));
                    }
                    None => {
                        // Default: zero flux at outlet (closed boundary)
                        BorderConditions.insert(variable.clone(), (1, 1e-10));
                    }
                }
            } else if !key.starts_with("Teta") && !key.starts_with("q") {
                // Concentration variable: set at inlet (left boundary)
                match border_conditions.get(key) {
                    Some(condition_value) => {
                        BorderConditions.insert(variable.clone(), (0, *condition_value));
                    }
                    None => {}
                }
            }
        }

        // Process temperature boundary condition with dimensionless scaling
        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        // Convert physical temperature to dimensionless form: Θ = (T - dT)/dT
        match border_conditions.get("T") {
            Some(T0) => {
                let Teta = (*T0 - dT) / T_scale; // Dimensionless temperature at inlet
                BorderConditions.insert("Teta".to_string(), (0, Teta)); // Set at left boundary
            }
            None => {
                return Err(ReactorError::MissingData(
                    "Initial temperature boundary condition 'T' not found".to_string(),
                ));
            }
        }

        // Process heat flux boundary condition
        match border_conditions.get("q") {
            Some(q_final) => {
                BorderConditions.insert("q".to_string(), (1, *q_final)); // Set at outlet
            }
            None => {
                // Default: zero heat flux at outlet (adiabatic boundary)
                BorderConditions.insert("q".to_string(), (1, 1e-10));
            }
        }

        // Store boundary conditions in solver
        self.solver.BorderConditions = BorderConditions;
        Ok(())
    }
}
