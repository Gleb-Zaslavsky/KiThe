use nalgebra::DMatrix;
#[derive(Debug, Clone, Copy)]
pub enum PhaseKind {
    IdealGas,
    IdealSolution,
}

#[derive(Debug, Clone)]
pub struct Phase {
    pub kind: PhaseKind,
    /// Indices of species belonging to this phase
    pub species: Vec<usize>,
}

use std::sync::Arc;

pub type GibbsFn = Arc<dyn Fn(f64) -> f64 + Send + Sync>;

pub fn species_to_phase_map(phases: &[Phase], num_species: usize) -> Vec<usize> {
    let mut map = vec![None; num_species];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            if map[i].is_some() {
                panic!(
                    "Species {} appears in more than one phase. \
                     Duplicate species per phase instead (e.g. O2_g, O2_l).",
                    i
                );
            }
            map[i] = Some(j);
        }
    }

    map.into_iter()
        .map(|x| x.expect("Species not assigned to any phase"))
        .collect()
}

/// Δn[k][j] = sum of stoichiometric coefficients of reaction k
///             for species belonging to phase j
pub fn reaction_phase_stoichiometry(
    reactions: &DMatrix<f64>, // m × r
    phases: &[Phase],
) -> Vec<Vec<f64>> {
    let r = reactions.ncols();

    let mut delta_n = vec![vec![0.0; phases.len()]; r];

    for (j, phase) in phases.iter().enumerate() {
        for &i in &phase.species {
            for k in 0..r {
                delta_n[k][j] += reactions[(i, k)];
            }
        }
    }

    delta_n
}

pub fn multiphase_equilibrium_residual_generator(
    n0: Vec<f64>,
    reactions: DMatrix<f64>, // m × r
    gibbs: Vec<GibbsFn>,     // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    r_gas: f64,
) -> impl Fn(&[f64]) -> Vec<f64> {
    let m = n0.len();
    let r = reactions.ncols();

    assert_eq!(reactions.nrows(), m);
    assert_eq!(gibbs.len(), m);

    let species_phase = species_to_phase_map(&phases, m);
    let delta_n = reaction_phase_stoichiometry(&reactions, &phases);

    move |xi: &[f64]| {
        assert_eq!(xi.len(), r);

        // --- species mole numbers ---
        let mut n = n0.clone();
        for i in 0..m {
            for k in 0..r {
                n[i] += reactions[(i, k)] * xi[k];
            }
        }

        // --- phase mole totals ---
        let mut n_phase = vec![0.0; phases.len()];
        for i in 0..m {
            let j = species_phase[i];
            n_phase[j] += n[i];
        }

        // --- phase log-activity offsets φ_j ---
        let mut phi = vec![0.0; phases.len()];
        for (j, phase) in phases.iter().enumerate() {
            phi[j] = match phase.kind {
                PhaseKind::IdealGas => (pressure / p0).ln(),
                PhaseKind::IdealSolution => 0.0,
            };
        }

        // --- residuals ---
        let mut f = vec![0.0; r];

        for k in 0..r {
            let mut sum_ln_n = 0.0;
            let mut dg0 = 0.0;

            for i in 0..m {
                sum_ln_n += reactions[(i, k)] * n[i].ln();
                dg0 += reactions[(i, k)] * gibbs[i](temperature);
            }

            let mut sum_ln_n_phase = 0.0;
            let mut sum_phi = 0.0;

            for j in 0..phases.len() {
                let dnk = delta_n[k][j];
                if dnk != 0.0 {
                    sum_ln_n_phase += dnk * n_phase[j].ln();
                    sum_phi += dnk * phi[j];
                }
            }

            let ln_k = -dg0 / (r_gas * temperature);

            f[k] = sum_ln_n - sum_ln_n_phase + sum_phi - ln_k;
        }

        f
    }
}
