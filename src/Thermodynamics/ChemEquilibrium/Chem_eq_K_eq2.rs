use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq::{
    EquilibriumLogMoles, GibbsFn, Phase, PhaseKind, R, Solvers, reaction_phase_stoichiometry,
    species_to_phase_map,
};
use crate::Thermodynamics::ChemEquilibrium::Chem_eq_K_eq3::ReactionExtentError;
use crate::Thermodynamics::User_substances::{LibraryPriority, Phases, SubsData};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use log::{error, info};
use nalgebra::DMatrix;
use std::collections::HashMap;
use std::default::Default;
use std::rc::Rc;
//////////////////////////////////////////NEW PHASES/////////////////////////////////////////////////////////////
///
///
///
/// For each phase P not currently active, we compute:
/// ΔG = ∑_i μ_i*n_i
/// If ΔG < 0, the phase is thermodynamically stable → must appear.

pub fn compute_phase_creation_dg(
    y: &[f64], // log-moles (used only for μ reference)
    gibbs: &[GibbsFn],
    phases: &[Phase],
    species_phase: &[usize],
    temperature: f64,
    pressure: f64,
    p0: f64,
) -> Vec<f64> {
    let m = y.len();
    let n_phases = phases.len();
    let rt = R * temperature;

    // --- phase activity offsets ---
    let mut phi = vec![0.0; n_phases];
    for (j, phase) in phases.iter().enumerate() {
        phi[j] = match phase.kind {
            PhaseKind::IdealGas => (pressure / p0).ln(),
            PhaseKind::IdealSolution => 0.0,
        };
    }

    let mut dg = vec![0.0; n_phases];

    for i in 0..m {
        let p = species_phase[i];
        let mu_i = gibbs[i](temperature) + rt * phi[p];
        dg[p] += mu_i;
    }

    dg
}

/// A phase is inactive if:
/// n_phase<phase_epsn

pub fn detect_inactive_phases(n_phase: &[f64], phase_eps: f64) -> Vec<usize> {
    let mut inactive = Vec::new();

    for (j, &np) in n_phase.iter().enumerate() {
        if np < phase_eps {
            inactive.push(j);
        }
    }

    inactive
}

pub fn activate_phase(
    y: &mut Vec<f64>, // log-moles (modified in place)
    phase_id: usize,
    species_phase: &[usize],
    elements: &DMatrix<f64>, // m × E
    element_totals: &[f64],  // E
    seed_fraction: f64,      // e.g. 1e-8
) {
    let m = y.len();
    let e = elements.ncols();

    // --- identify species in the phase ---
    let species: Vec<usize> = (0..m).filter(|&i| species_phase[i] == phase_id).collect();

    if species.is_empty() {
        return;
    }

    // --- compute current species moles ---
    let mut n = vec![0.0; m];
    for i in 0..m {
        n[i] = y[i].exp();
    }

    // --- compute elemental pool available ---
    let mut element_available = element_totals.to_vec();
    for el in 0..e {
        for i in 0..m {
            element_available[el] -= elements[(i, el)] * n[i];
        }
    }

    // --- distribute a small neutral seed ---
    for &i in &species {
        let mut max_seed = f64::INFINITY;

        for el in 0..e {
            let a = elements[(i, el)];
            if a > 0.0 {
                max_seed = max_seed.min(element_available[el] / a);
            }
        }

        let seed = seed_fraction * max_seed.max(0.0);
        if seed > 0.0 {
            n[i] += seed;
        }
    }

    // --- update log-moles ---
    for i in 0..m {
        y[i] = n[i].max(1e-300).ln();
    }
}

pub fn deactivate_phases(
    y: &mut Vec<f64>, // log-moles (modified in place)
    phases_to_deactivate: &[usize],
    species_phase: &[usize],
    log_floor: f64, // e.g. -700.0
) {
    let m = y.len();

    for &phase_id in phases_to_deactivate {
        for i in 0..m {
            if species_phase[i] == phase_id {
                y[i] = log_floor;
            }
        }
    }
}

pub fn deactivate_phases_mass_conserving(
    y: &mut Vec<f64>,
    phases_to_deactivate: &[usize],
    species_phase: &[usize],
    reference_phase: usize,
    log_floor: f64,
) {
    let m = y.len();
    let mut n = vec![0.0; m];

    for i in 0..m {
        n[i] = y[i].exp();
    }

    for &p in phases_to_deactivate {
        for i in 0..m {
            if species_phase[i] == p {
                let amount = n[i];
                n[i] = 0.0;

                // move to reference phase species (first match)
                for j in 0..m {
                    if species_phase[j] == reference_phase {
                        n[j] += amount;
                        break;
                    }
                }
            }
        }
    }

    for i in 0..m {
        y[i] = n[i].max(1e-300).ln();
    }
}

pub fn detect_inactive_species(
    y: &[f64],
    species_phase: &[usize],
    phase_active: &[bool],
    subs_eps: f64,
) -> Vec<bool> {
    let m = y.len();
    let mut active = vec![true; m];

    for i in 0..m {
        let n = y[i].exp();
        let p = species_phase[i];

        if !phase_active[p] && n < subs_eps {
            active[i] = false;
        }
    }

    active
}

pub struct PhaseManager {
    pub phase_eps: f64, // destruction threshold
    pub subs_eps: f64,
    pub dg_create: f64, // e.g. -1e-6 * RT
    pub dg_keep: f64,   // e.g. +1e-8 * RT
}

impl PhaseManager {
    pub fn new(phase_eps: f64, subs_eps: f64, dg_create: f64, dg_keep: f64) -> Self {
        Self {
            phase_eps,
            subs_eps,
            dg_create,
            dg_keep,
        }
    }
    pub fn set_thresholds(&mut self, T: f64) {
        self.dg_create = -1e-6 * R * T;
        self.dg_keep = 1e-8 * R * T;
    }
    /// Detect phases that must be created (ΔG < 0)
    pub fn detect_phase_creation(
        &self,
        y: &[f64],
        gibbs: &[GibbsFn],
        phases: &[Phase],
        species_phase: &[usize],
        temperature: f64,
        pressure: f64,
        p0: f64,
        active: &[bool],
    ) -> Vec<usize> {
        let dg =
            compute_phase_creation_dg(y, gibbs, phases, species_phase, temperature, pressure, p0);

        dg.iter()
            .enumerate()
            .filter(|&(p, &v)| !active[p] && v < self.dg_create)
            .map(|(p, _)| p)
            .collect()
    }

    /// Detect phases that must be removed (n_phase < eps)
    pub fn detect_phase_destruction(&self, n_phase: &[f64], active: &[bool]) -> Vec<usize> {
        n_phase
            .iter()
            .enumerate()
            .filter(|&(p, &n)| active[p] && n < self.phase_eps)
            .map(|(p, _)| p)
            .collect()
    }
}
impl Default for PhaseManager {
    fn default() -> Self {
        Self {
            phase_eps: 1e-30,
            subs_eps: 1e-30,
            dg_create: -1e-6 * 298.0 * R,
            dg_keep: 1e-8 * 298.0 * R,
        }
    }
}
pub fn compute_phase_totals(
    y: &[f64],               // log-moles
    species_phase: &[usize], // m
) -> Vec<f64> {
    let m = y.len();
    let n_phases = species_phase.iter().copied().max().unwrap_or(0) + 1;

    let mut n_phase = vec![0.0; n_phases];

    for i in 0..m {
        let p = species_phase[i];
        n_phase[p] += y[i].exp();
    }

    n_phase
}

/*

///
///  loop:
///     solve equilibrium (NR / trust region)
///
///     compute n_phase
///
///     if inactive phases detected:
///         deactivate phases
///         restart solver
///
///     if ΔG < 0 for missing phase:
///         activate phase
///         restart solver
///
///     break
///
///



*/
impl EquilibriumLogMoles {
    pub fn solve_with_phase_control(&mut self) -> Result<(), ReactionExtentError> {
        let elements = self.elem_composition.clone();
        let element_totals = self.elements_vector.clone();
        let species_phase = self.species_phase.clone();
        let n = self.phases.len();
        let mut phase_active = vec![true; self.phases.len()];

        loop {
            // 1. Solve equilibrium for current active phases
            self.solve()?;
            let phase_manager = &self.phase_manager;
            let mut y = self.solution.clone();
            // 2. Compute phase totals
            let n_phase = compute_phase_totals(&y, &species_phase);

            // 3. Detect phase destruction
            let destroyed = phase_manager.detect_phase_destruction(&n_phase, &phase_active);

            if !destroyed.is_empty() {
                info!("phases {:?} destroyed", destroyed);
                deactivate_phases(&mut y, &destroyed, &species_phase, -700.0);

                for &p in &destroyed {
                    phase_active[p] = false;
                }

                continue; // restart solver
            }

            // 4. Detect phase creation
            let created = phase_manager.detect_phase_creation(
                &y,
                &self.gibbs,
                &self.phases,
                &species_phase,
                self.T,
                self.P,
                self.p0,
                &phase_active,
            );

            if !created.is_empty() {
                info!("phases {:?} created", created);
                for &p in &created {
                    activate_phase(&mut y, p, &species_phase, &elements, &element_totals, 1e-8);
                    phase_active[p] = true;
                }

                continue; // restart solver
            }

            // 5. Nothing changed → converged equilibrium
            break;
        }
        Ok(())
    }
}
///////////////////////////SYMBOLIC//////////////////////////////////////////////////////////////////////
/// Generates symbolic residual expressions for equilibrium equations
///
/// Creates symbolic expressions for residuals that can be used for
/// analytical differentiation or symbolic manipulation.
pub fn multiphase_equilibrium_residual_generator_sym(
    reactions: DMatrix<f64>,  // m × r
    elements: DMatrix<f64>,   // m × E
    element_totals: Vec<f64>, // E
    gibbs_sym: Vec<Expr>,     // m
    phases: Vec<Phase>,

    pressure: f64,
    p0: f64,
) -> Result<Vec<Expr>, ReactionExtentError> {
    let m = reactions.nrows();
    let r = reactions.ncols();
    let e = elements.ncols();

    if gibbs_sym.len() != m || element_totals.len() != e {
        return Err(ReactionExtentError::Other(
            "Dimension mismatch in residual generator".to_string(),
        ));
    }

    let species_phase = species_to_phase_map(&phases, m)?;
    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);
    let rt = Expr::Const(R) * Expr::Var("T".to_string());

    // log-mole symbolic variables (y_i = ln n_i)
    let y = Expr::IndexedVars(m, "y").0;
    if y.len() != m {
        return Err(ReactionExtentError::ResidualEvaluation(
            "y length mismatch".to_string(),
        ));
    }

    // --- species mole numbers as n_i = exp(y_i) ---
    let mut n: Vec<Expr> = Vec::with_capacity(m);
    for i in 0..m {
        n.push(y[i].clone().exp());
    }

    // --- phase mole totals ---
    let mut n_phase: Vec<Expr> = vec![Expr::Const(0.0); phases.len()];
    for i in 0..m {
        let j = species_phase[i];
        n_phase[j] += n[i].clone();
    }

    // --- phase offsets ---
    let mut phi: Vec<Expr> = vec![Expr::Const(0.0); phases.len()];
    for (j, phase) in phases.iter().enumerate() {
        phi[j] = match phase.kind {
            PhaseKind::IdealGas => Expr::Const((pressure / p0).ln()),
            PhaseKind::IdealSolution => Expr::Const(0.0),
        };
    }

    // Build residuals: r reaction eqs followed by e element balance eqs
    let mut f: Vec<Expr> = vec![Expr::Const(0.0); r + e];

    for k in 0..r {
        let mut sum_ln_n = Expr::Const(0.0);
        let mut sum_ln_n_phase = Expr::Const(0.0);
        let mut sum_phi = Expr::Const(0.0);
        let mut dg0 = Expr::Const(0.0);

        // species contributions: ν_{ik} * y_i (since y_i = ln n_i) and ν_{ik} * g0_i(T)
        for i in 0..m {
            let nu = reactions[(i, k)];
            if nu == 0.0 {
                continue;
            }
            let nu_expr = Expr::Const(nu);
            sum_ln_n += nu_expr.clone() * y[i].clone();
            dg0 += nu_expr.clone() * gibbs_sym[i].clone();
        }

        // phase contributions: Δn_{kj} * ln(N_j) and Δn_{kj} * φ_j
        for j in 0..phases.len() {
            let dnk = delta_n[k][j];
            if dnk == 0.0 {
                continue;
            }
            let dnk_expr = Expr::Const(dnk);
            let nj = n_phase[j].clone();
            sum_ln_n_phase += dnk_expr.clone() * nj.ln();
            sum_phi += dnk_expr.clone() * phi[j].clone();
        }

        let ln_k = -dg0 / rt.clone();
        f[k] = (sum_ln_n - sum_ln_n_phase + sum_phi - ln_k).simplify();
    }

    // element balance equations: sum_i a_{i,el} * n_i - element_totals[el]
    for el in 0..e {
        let mut sum_el = Expr::Const(0.0);
        for i in 0..m {
            sum_el += Expr::Const(elements[(i, el)]) * n[i].clone();
        }
        f[r + el] = (sum_el - Expr::Const(element_totals[el])).simplify();
    }

    Ok(f)
}

////////////////////////CONVENIENCE FUNCTIONS////////////////////////////////
/// Convenience function for gas-phase equilibrium calculations
///
/// Sets up complete equilibrium solver for gas-phase systems with
/// automatic database searching and Gibbs function generation.
///
/// # Arguments
/// * `subs` - List of substance names
/// * `T` - Temperature in K
/// * `P` - Pressure in Pa
/// * `solver` - Choice of numerical solver
/// * `loglevel` - Logging level (None to disable)
/// * `scaling` - Whether to enable equation scaling
/// the function don't run calculation by itself only prepare data
/// for calculation
pub fn gas_solver(
    subs: Vec<String>,
    T: f64,
    P: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    scaling: bool,
) -> EquilibriumLogMoles {
    let mut user_subs = SubsData::new();

    user_subs.substances = subs;
    let map_of_phases = user_subs
        .substances
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;
    user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
    user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

    let _ = user_subs.search_substances().unwrap();
    let _ = user_subs.parse_all_thermal_coeffs();
    let _ = user_subs.extract_all_thermal_coeffs(T);
    let _ = user_subs.calculate_elem_composition_and_molar_mass(None);

    let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();
    let mut g_vec: Vec<GibbsFn> = Vec::new();

    for sub in &user_subs.substances {
        let boxed_fn = map_of_functions.remove(sub).unwrap();
        let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
        g_vec.push(gibbs_fn);
    }
    let mut instance = EquilibriumLogMoles::new();

    instance.elem_composition = user_subs.elem_composition_matrix.clone().unwrap().clone();
    instance.subs_data = user_subs;
    instance.with_loglevel(loglevel);
    instance.solver = solver;
    instance.scaling_flag = scaling;
    instance.gibbs = g_vec;

    instance.P = P;
    instance.T = T;
    instance
}

pub fn gas_solver_from_elements(
    elements: Vec<String>,
    map_of_nonzero_moles: HashMap<String, f64>,
    T: f64,
    P: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    scaling: bool,
) -> EquilibriumLogMoles {
    let mut user_subs = SubsData::new();

    let _ = user_subs.search_by_elements_only(elements);
    let subs = user_subs.substances.clone();
    let _ = user_subs.parse_all_thermal_coeffs();
    let _ = user_subs.extract_all_thermal_coeffs(T);
    let _ = user_subs.calculate_elem_composition_and_molar_mass(None);

    let mut map_of_functions = user_subs.calculate_dG0_fun_one_phase();
    let mut g_vec: Vec<GibbsFn> = Vec::new();

    for sub in &subs {
        let boxed_fn = map_of_functions.remove(sub).unwrap();
        let gibbs_fn: GibbsFn = Rc::new(move |T: f64| boxed_fn(T));
        g_vec.push(gibbs_fn);
    }
    let mut instance = EquilibriumLogMoles::new();
    instance.with_loglevel(loglevel);
    instance.elem_composition = user_subs.elem_composition_matrix.clone().unwrap().clone();
    instance.subs_data = user_subs;
    info!("calculating n0");
    instance.set_n0_from_non_zero_map(map_of_nonzero_moles);
    let n0 = instance.n0.clone();

    instance.solver = solver;
    instance.scaling_flag = scaling;
    instance.gibbs = g_vec;
    let species: Vec<usize> = (0..subs.len()).map(|x| x as usize).collect();

    instance.set_problem(
        n0,
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: species,
        }],
        P,
    );
    info!(" phase data{:?}", instance.phases);
    instance.P = P;
    instance.T = T;

    instance
}

pub fn gas_solver_for_T_range(
    subs: Vec<String>,
    n0: Vec<f64>,
    P: f64,
    T_start: f64,
    T_end: f64,
    T_step: f64,
    solver: Solvers,
    loglevel: Option<&str>,
    parallel: bool,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    let mut user_subs = SubsData::new();

    user_subs.substances = subs.clone();
    let map_of_phases = user_subs
        .substances
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;

    user_subs.set_library_priority("NASA_gas".to_string(), LibraryPriority::Priority);
    user_subs.set_library_priority("NUIG_thermo".to_string(), LibraryPriority::Permitted);

    let mut instance = EquilibriumLogMoles::new();
    instance.scaling_flag = false;
    instance.with_loglevel(loglevel);
    let species: Vec<usize> = (0..n0.len()).map(|x| x as usize).collect();
    assert_eq!(subs.len(), species.len());
    instance.set_problem(
        n0,
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: species,
        }],
        P,
    );
    instance.subs_data = user_subs;
    instance.solver = solver;
    instance.create_equilibrium_system()?;
    if parallel {
        instance.solve_for_T_range_par2(T_start, T_end, T_step)?
    } else {
        instance.solve_for_T_range(T_start, T_end, T_step)?
    };
    instance.create_moles_table();
    Ok(instance)
}

pub fn gas_solver_for_T_range_for_elements(
    elements: Vec<String>,
    map_of_nonzero_moles: HashMap<String, f64>,
    P: f64,
    T_start: f64,
    T_end: f64,
    T_step: f64,
    solver: Solvers,
    loglevel: Option<&str>,
) -> Result<EquilibriumLogMoles, ReactionExtentError> {
    let mut user_subs = SubsData::new();
    let _ = user_subs.search_by_elements_only(elements);
    let _ = user_subs.parse_all_thermal_coeffs();
    let _ = user_subs.calculate_elem_composition_and_molar_mass(None);
    let subs = user_subs.substances.clone();
    let map_of_phases = subs
        .clone()
        .iter()
        .map(|s| (s.clone(), Some(Phases::Gas)))
        .collect();
    user_subs.map_of_phases = map_of_phases;

    let mut instance = EquilibriumLogMoles::new();
    instance.elem_composition = user_subs.elem_composition_matrix.clone().unwrap().clone();
    instance.subs_data = user_subs;
    instance.set_n0_from_non_zero_map(map_of_nonzero_moles);
    let n0 = instance.n0.clone();
    instance.scaling_flag = false;
    instance.with_loglevel(loglevel);
    let species: Vec<usize> = (0..n0.len()).map(|x| x as usize).collect();

    instance.set_problem(
        n0,
        vec![Phase {
            kind: PhaseKind::IdealGas,
            species: species,
        }],
        P,
    );

    instance.solver = solver;
    let _ = instance.create_stoich_matrix();
    //  instance.create_equilibrium_system()?;
    instance.solve_for_T_range(T_start, T_end, T_step)?;
    instance.create_moles_table();
    Ok(instance)
}

////////////////////////misc///////////////////////////////
/// Computes finite difference approximation of Jacobian matrix
///
/// Used for testing or when analytical Jacobian is not available.
/// Less accurate and slower than analytical Jacobian.
pub fn finite_difference_jacobian<F>(
    f: &F,
    x: &[f64],
    eps: f64,
) -> Result<DMatrix<f64>, ReactionExtentError>
where
    F: Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>,
{
    let f0 = f(x)?;
    let n = x.len();
    let m = f0.len();

    let mut j = DMatrix::<f64>::zeros(m, n);

    for k in 0..n {
        let mut x_pert = x.to_vec();
        x_pert[k] += eps;

        let f1 = f(&x_pert)?;

        for i in 0..m {
            j[(i, k)] = (f1[i] - f0[i]) / eps;
        }
    }

    Ok(j)
}
