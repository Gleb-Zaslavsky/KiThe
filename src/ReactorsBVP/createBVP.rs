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

use super::SimpleReactorBVP::{ReactorError, SimpleReactorTask, SymbolicRhsAssemblyPolicy};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use rayon::prelude::*;
use std::collections::HashMap;

const AUTO_PARALLEL_ASSEMBLY_THRESHOLD: usize = 20;

#[derive(Debug)]
struct BvpEquationAssemblyContext<'a> {
    substances: &'a [String],
    equations: &'a [String],
    k_sym_vec: &'a [Expr],
    stoich_reags: &'a [Vec<f64>],
    stoich_matrix: &'a [Vec<f64>],
    thermal_effects: &'a [f64],
    pe_d: &'a [f64],
    d_ro_values: Vec<f64>,
    molar_masses_kg: Vec<f64>,
    ci_expr: Vec<Expr>,
    ci_names: Vec<String>,
    ji_expr: Vec<Expr>,
    ji_names: Vec<String>,
    ro_m: Expr,
    qm: Expr,
    qc: Expr,
    pe_q: f64,
}

#[derive(Debug)]
struct SpeciesEquationAssembly {
    index: usize,
    substance: String,
    ci_name: String,
    ji_name: String,
    ci_rhs: Expr,
    ji_rhs: Expr,
}

impl SimpleReactorTask {
    /// Build a normalized snapshot of the kinetic state used by BVP assembly.
    fn build_bvp_equation_context(&self) -> Result<BvpEquationAssemblyContext<'_>, ReactorError> {
        let substances = &self.kindata.substances;
        let n = substances.len();
        if n == 0 {
            return Err(ReactorError::MissingData(
                "No substances found in kinetic data".to_string(),
            ));
        }

        let equations = &self.kindata.vec_of_equations;
        let k = equations.len();
        if k == 0 {
            return Err(ReactorError::MissingData(
                "No reactions found in kinetic data".to_string(),
            ));
        }

        let k_sym_vec = self.kindata.K_sym_vec.as_ref().ok_or_else(|| {
            ReactorError::MissingData("Symbolic rate constants not calculated".to_string())
        })?;
        if k_sym_vec.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Mismatch between number of reactions and rate constants".to_string(),
            ));
        }

        let stoich_reags = &self.kindata.stecheodata.stecheo_reags;
        let stoich_matrix = &self.kindata.stecheodata.stecheo_matrx;
        if stoich_reags.len() != k || stoich_matrix.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Stoichiometric reaction data size mismatch".to_string(),
            ));
        }
        for (idx, row) in stoich_reags.iter().enumerate() {
            if row.len() != n {
                return Err(ReactorError::InvalidConfiguration(format!(
                    "stoicheodata.stecheo_reags[{}] length {} != substances {}",
                    idx,
                    row.len(),
                    n
                )));
            }
        }
        for (idx, row) in stoich_matrix.iter().enumerate() {
            if row.len() != n {
                return Err(ReactorError::InvalidConfiguration(format!(
                    "stoicheodata.stecheo_matrx[{}] length {} != substances {}",
                    idx,
                    row.len(),
                    n
                )));
            }
        }

        let thermal_effects = self.thermal_effects.as_slice();
        if thermal_effects.len() != k {
            return Err(ReactorError::InvalidConfiguration(
                "Thermal effects vector size mismatch".to_string(),
            ));
        }

        let pe_d = self.Pe_D.as_slice();
        if pe_d.len() != n {
            return Err(ReactorError::InvalidConfiguration(
                "Peclet numbers vector size mismatch".to_string(),
            ));
        }

        for substance in substances.iter() {
            if !self.D_ro_map.contains_key(substance) {
                return Err(ReactorError::MissingData(format!(
                    "Transport coefficient for `{}` is missing from the normalized snapshot",
                    substance
                )));
            }
        }
        let mut d_ro_values = Vec::with_capacity(n);
        for substance in substances.iter() {
            let d_ro = *self.D_ro_map.get(substance).ok_or_else(|| {
                ReactorError::MissingData(format!(
                    "Transport coefficient for `{}` is missing from the normalized snapshot",
                    substance
                ))
            })?;
            if !d_ro.is_finite() || d_ro <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Transport coefficient for `{}` must be finite and positive, got {}",
                    substance, d_ro
                )));
            }
            d_ro_values.push(d_ro);
        }

        let molar_masses = self
            .kindata
            .stecheodata
            .vec_of_molmasses
            .as_ref()
            .ok_or_else(|| ReactorError::MissingData("Molar masses not calculated".to_string()))?;
        if molar_masses.len() != n {
            return Err(ReactorError::InvalidConfiguration(
                "Molar masses vector size mismatch".to_string(),
            ));
        }

        let mut molar_masses_kg = Vec::with_capacity(n);
        for (idx, mass) in molar_masses.iter().enumerate() {
            if !mass.is_finite() || *mass <= 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Molar mass at index {} must be finite and positive, got {}",
                    idx, mass
                )));
            }
            molar_masses_kg.push(*mass / 1000.0);
        }

        let (ci_expr, ci_names_raw) = Expr::IndexedVars(n, "C");
        let (ji_expr, ji_names_raw) = Expr::IndexedVars(n, "J");
        let ci_names: Vec<String> = ci_names_raw
            .iter()
            .map(|var| var.replace("_", ""))
            .collect();
        let ji_names: Vec<String> = ji_names_raw
            .iter()
            .map(|var| var.replace("_", ""))
            .collect();

        Ok(BvpEquationAssemblyContext {
            substances,
            equations,
            k_sym_vec,
            stoich_reags,
            stoich_matrix,
            thermal_effects,
            pe_d,
            d_ro_values,
            molar_masses_kg,
            ci_expr,
            ci_names,
            ji_expr,
            ji_names,
            ro_m: Expr::Const(self.ideal_gas_density()?),
            qm: Expr::Const(self.L.powi(2) / self.scaling.T_scale),
            qc: Expr::Const(self.L.powi(2)),
            pe_q: self.Pe_q,
        })
    }

    /// Resolve the effective symbolic assembly policy.
    fn resolve_symbolic_rhs_assembly_policy(
        &self,
        equation_count: usize,
    ) -> SymbolicRhsAssemblyPolicy {
        match self.symbolic_rhs_assembly_policy {
            SymbolicRhsAssemblyPolicy::Auto => {
                if equation_count < AUTO_PARALLEL_ASSEMBLY_THRESHOLD {
                    SymbolicRhsAssemblyPolicy::Sequential
                } else {
                    SymbolicRhsAssemblyPolicy::Parallel
                }
            }
            policy => policy,
        }
    }

    /// Build the symbolic reaction rate for a single reaction.
    fn build_rate_expression(
        context: &BvpEquationAssemblyContext<'_>,
        reaction_index: usize,
    ) -> Result<Expr, ReactorError> {
        let mut rate_expr = context
            .k_sym_vec
            .get(reaction_index)
            .ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "k_sym_vec[{}] out of bounds",
                    reaction_index
                ))
            })?
            .clone();

        let reaction_powers = context.stoich_reags.get(reaction_index).ok_or_else(|| {
            ReactorError::IndexOutOfBounds(format!(
                "stoich_reags[{}] out of bounds",
                reaction_index
            ))
        })?;

        for species_index in 0..context.substances.len() {
            let power = reaction_powers.get(species_index).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "stoich_reags[{}][{}] out of bounds",
                    reaction_index, species_index
                ))
            })?;
            let ci_expr = context.ci_expr.get(species_index).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!("ci_expr[{}] out of bounds", species_index))
            })?;
            let molar_mass = context.molar_masses_kg.get(species_index).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "molar_masses_kg[{}] out of bounds",
                    species_index
                ))
            })?;
            rate_expr = rate_expr
                * (context.ro_m.clone() * ci_expr.clone() / Expr::Const(*molar_mass))
                    .pow(Expr::Const(*power));
        }

        Ok(rate_expr.simplify_())
    }

    /// Build ordered rate expressions and the rate lookup map.
    fn build_rate_expressions(
        &self,
        context: &BvpEquationAssemblyContext<'_>,
        policy: SymbolicRhsAssemblyPolicy,
    ) -> Result<(HashMap<String, Expr>, Vec<Expr>), ReactorError> {
        let reaction_count = context.equations.len();
        let mut rate_entries: Vec<(usize, String, Expr)> = match policy {
            SymbolicRhsAssemblyPolicy::Sequential => {
                let mut entries = Vec::with_capacity(reaction_count);
                for reaction_index in 0..reaction_count {
                    let rate_expr = Self::build_rate_expression(context, reaction_index)?;
                    entries.push((
                        reaction_index,
                        context.equations[reaction_index].clone(),
                        rate_expr,
                    ));
                }
                entries
            }
            SymbolicRhsAssemblyPolicy::Parallel => (0..reaction_count)
                .into_par_iter()
                .map(|reaction_index| {
                    Self::build_rate_expression(context, reaction_index).map(|rate_expr| {
                        (
                            reaction_index,
                            context.equations[reaction_index].clone(),
                            rate_expr,
                        )
                    })
                })
                .collect::<Result<Vec<_>, ReactorError>>()?,
            SymbolicRhsAssemblyPolicy::Auto => {
                unreachable!("Auto policy should be resolved before assembly")
            }
        };

        rate_entries.sort_by_key(|entry| entry.0);
        let mut map_eq_rate = HashMap::with_capacity(rate_entries.len());
        let mut rates = Vec::with_capacity(rate_entries.len());
        for (_, equation, rate_expr) in rate_entries {
            map_eq_rate.insert(equation, rate_expr.clone());
            rates.push(rate_expr);
        }
        Ok((map_eq_rate, rates))
    }

    /// Build the symbolic equation block for one species.
    fn build_species_equation(
        context: &BvpEquationAssemblyContext<'_>,
        rates: &[Expr],
        species_index: usize,
    ) -> Result<SpeciesEquationAssembly, ReactorError> {
        let substance = context
            .substances
            .get(species_index)
            .ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "substances[{}] out of bounds",
                    species_index
                ))
            })?
            .clone();
        let ci_name = context
            .ci_names
            .get(species_index)
            .ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!("ci_names[{}] out of bounds", species_index))
            })?
            .clone();
        let ji_name = context
            .ji_names
            .get(species_index)
            .ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!("ji_names[{}] out of bounds", species_index))
            })?
            .clone();
        let ji_expr = context
            .ji_expr
            .get(species_index)
            .ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!("ji_expr[{}] out of bounds", species_index))
            })?
            .clone();

        let d_ro = *context.d_ro_values.get(species_index).ok_or_else(|| {
            ReactorError::IndexOutOfBounds(format!("d_ro_values[{}] out of bounds", species_index))
        })?;
        let ci_rhs = ji_expr.clone() / Expr::Const(d_ro);

        let mut wi = Expr::Const(0.0);
        for reaction_index in 0..rates.len() {
            let stoich_coeff = context
                .stoich_matrix
                .get(reaction_index)
                .and_then(|row| row.get(species_index))
                .ok_or_else(|| {
                    ReactorError::IndexOutOfBounds(format!(
                        "stoich_matrix[{}][{}] out of bounds",
                        reaction_index, species_index
                    ))
                })?;
            let molar_mass = context.molar_masses_kg.get(species_index).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "molar_masses_kg[{}] out of bounds",
                    species_index
                ))
            })?;
            wi = wi
                + rates[reaction_index].clone()
                    * Expr::Const(*stoich_coeff)
                    * Expr::Const(*molar_mass);
        }

        let ji_rhs = (ji_expr.clone() * Expr::Const(context.pe_d[species_index])
            - wi * context.qc.clone())
        .simplify_();

        Ok(SpeciesEquationAssembly {
            index: species_index,
            substance,
            ci_name,
            ji_name,
            ci_rhs,
            ji_rhs,
        })
    }

    /// Build all species equation blocks using the selected policy.
    fn build_species_equations(
        &self,
        context: &BvpEquationAssemblyContext<'_>,
        rates: &[Expr],
        policy: SymbolicRhsAssemblyPolicy,
    ) -> Result<Vec<SpeciesEquationAssembly>, ReactorError> {
        let species_count = context.substances.len();
        let mut blocks: Vec<SpeciesEquationAssembly> = match policy {
            SymbolicRhsAssemblyPolicy::Sequential => {
                let mut entries = Vec::with_capacity(species_count);
                for species_index in 0..species_count {
                    entries.push(Self::build_species_equation(context, rates, species_index)?);
                }
                entries
            }
            SymbolicRhsAssemblyPolicy::Parallel => (0..species_count)
                .into_par_iter()
                .map(|species_index| Self::build_species_equation(context, rates, species_index))
                .collect::<Result<Vec<_>, ReactorError>>()?,
            SymbolicRhsAssemblyPolicy::Auto => {
                unreachable!("Auto policy should be resolved before assembly")
            }
        };
        blocks.sort_by_key(|block| block.index);
        Ok(blocks)
    }

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
    pub fn create_bvp_equations(&mut self) -> Result<(), ReactorError> {
        let context = self.build_bvp_equation_context()?;
        let policy = self.resolve_symbolic_rhs_assembly_policy(context.equations.len());
        let (map_eq_rate, rates) = self.build_rate_expressions(&context, policy)?;
        let species_blocks = self.build_species_equations(&context, &rates, policy)?;
        let substance_count = context.substances.len();
        // Keep the normalized snapshot borrowed here so we do not duplicate
        // reaction metadata before the backend assembly starts.
        let thermal_effects = context.thermal_effects;
        let pe_q = context.pe_q;
        let qm = &context.qm;

        let mut vec_of_unknowns: Vec<String> = Vec::with_capacity(2 + 2 * substance_count);
        let mut vec_of_equations: Vec<Expr> = Vec::with_capacity(2 + 2 * substance_count);
        let mut map_of_equations: HashMap<String, (String, Expr)> =
            HashMap::with_capacity(2 + 2 * substance_count);

        let q = Expr::Var("q".to_string());
        let rhs_teta = q.clone() / Expr::Const(self.Lambda);
        vec_of_unknowns.push("Teta".to_string());
        vec_of_equations.push(rhs_teta.clone());
        map_of_equations.insert("Teta".to_string(), ("Teta".to_string(), rhs_teta));

        let mut heat_release = Expr::Const(0.0);
        for (reaction_index, rate_expr) in rates.iter().enumerate() {
            heat_release = (heat_release
                + Expr::Const(thermal_effects[reaction_index]) * rate_expr.clone())
            .simplify_();
        }
        self.heat_release = heat_release.clone();

        let rhs_q = q * Expr::Const(pe_q) - heat_release * qm.clone();
        vec_of_unknowns.push("q".to_string());
        vec_of_equations.push(rhs_q.clone());
        map_of_equations.insert("q".to_string(), ("q".to_string(), rhs_q));

        for block in species_blocks {
            vec_of_unknowns.push(block.ci_name.clone());
            vec_of_equations.push(block.ci_rhs.clone());
            map_of_equations.insert(block.substance.clone(), (block.ci_name, block.ci_rhs));

            vec_of_unknowns.push(block.ji_name.clone());
            vec_of_equations.push(block.ji_rhs.clone());
            map_of_equations.insert(
                format!("{}_flux", block.substance),
                (block.ji_name, block.ji_rhs),
            );
        }

        self.map_eq_rate = map_eq_rate;
        self.solver.eq_system = vec_of_equations;
        self.solver.unknowns = vec_of_unknowns;
        self.map_of_equations = map_of_equations;
        return Ok(());
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
        let border_conditions = &self.boundary_condition;
        let substances = &self.kindata.substances;
        let expected_len = 2 * substances.len() + 2;
        if self.solver.unknowns.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "Unknowns length {} does not match expected canonical size {}",
                self.solver.unknowns.len(),
                expected_len
            )));
        }
        if self.map_of_equations.len() != expected_len {
            return Err(ReactorError::InvalidConfiguration(format!(
                "map_of_equations length {} does not match expected canonical size {}",
                self.map_of_equations.len(),
                expected_len
            )));
        }

        // Boundary condition format: {variable → (boundary_index, value)}
        // boundary_index: 0 = left boundary (inlet, z=0), 1 = right boundary (outlet, z=1)
        // variable: solver variable name ("Teta", "q", "C0", "J0", etc.)
        let mut border_conditions_map: HashMap<String, (usize, f64)> =
            HashMap::with_capacity(expected_len);

        // Process temperature boundary condition with dimensionless scaling
        let dT = self.scaling.dT;
        let T_scale = self.scaling.T_scale;
        // Convert physical temperature to dimensionless form: Θ = (T - dT)/dT
        match border_conditions.get("T") {
            Some(T0) => {
                if !T0.is_finite() || *T0 <= 0.0 {
                    return Err(ReactorError::InvalidNumericValue(format!(
                        "Boundary temperature must be finite and positive, got {}",
                        T0
                    )));
                }
                let Teta = (*T0 - dT) / T_scale; // Dimensionless temperature at inlet
                border_conditions_map.insert("Teta".to_string(), (0, Teta)); // Set at left boundary
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
                if !q_final.is_finite() {
                    return Err(ReactorError::InvalidNumericValue(format!(
                        "Boundary heat flux must be finite, got {}",
                        q_final
                    )));
                }
                border_conditions_map.insert("q".to_string(), (1, *q_final)); // Set at outlet
            }
            None => {
                // Default: zero heat flux at outlet (adiabatic boundary)
                border_conditions_map.insert("q".to_string(), (1, 1e-10));
            }
        }

        for (index, substance) in substances.iter().enumerate() {
            let ci_name = self.solver.unknowns.get(2 + 2 * index).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "unknowns[{}] out of bounds while building concentration boundary conditions",
                    2 + 2 * index
                ))
            })?;
            let ji_name = self.solver.unknowns.get(2 + 2 * index + 1).ok_or_else(|| {
                ReactorError::IndexOutOfBounds(format!(
                    "unknowns[{}] out of bounds while building flux boundary conditions",
                    2 + 2 * index + 1
                ))
            })?;

            let concentration = border_conditions.get(substance).ok_or_else(|| {
                ReactorError::MissingData(format!("Missing boundary condition for `{}`", substance))
            })?;
            if !concentration.is_finite() || *concentration < 0.0 {
                return Err(ReactorError::InvalidNumericValue(format!(
                    "Boundary fraction for `{}` must be finite and non-negative, got {}",
                    substance, concentration
                )));
            }
            border_conditions_map.insert(ci_name.clone(), (0, *concentration));

            let mut flux_key = String::with_capacity(substance.len() + 5);
            flux_key.push_str(substance);
            flux_key.push_str("_flux");
            match border_conditions.get(&flux_key) {
                Some(condition_value) => {
                    if !condition_value.is_finite() {
                        return Err(ReactorError::InvalidNumericValue(format!(
                            "Boundary flux for `{}` must be finite, got {}",
                            substance, condition_value
                        )));
                    }
                    border_conditions_map.insert(ji_name.clone(), (1, *condition_value));
                }
                None => {
                    // Default: zero flux at outlet (closed boundary).
                    border_conditions_map.insert(ji_name.clone(), (1, 1e-10));
                }
            }
        }

        // Store boundary conditions in solver
        self.solver.BorderConditions = border_conditions_map;
        Ok(())
    }
}
