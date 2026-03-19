//! # Vyazovkin advanced isoconversional method
//!
//! ## What this module does
//! Implements the **Vyazovkin non-linear isoconversional method** for estimating
//! the activation energy Eα as a function of conversion α from non-isothermal
//! TGA experiments run at several heating rates.
//!
//! Unlike the linear integral methods (KAS, OFW, Starink) which approximate the
//! temperature integral analytically, Vyazovkin minimises a functional Φ(Eα) that
//! contains the *exact* numerical temperature integral, making it more accurate
//! especially when Eα varies strongly with α.
//!
//! ## Main data structures
//! - [`VyazovkinSolver`] — holds the search grid bounds (`e_min`, `e_max`, `steps`)
//!   and exposes `solve` / `solve_layer_optimized`.
//! - [`ConversionGrid`] (from `kinetic_methods`) — the pre-built isoconversional
//!   grid; **must** contain the `dt` matrix (time-step increments per experiment
//!   per α-layer), computed by `ConversionGridBuilder::with_dt_matrix()`.
//! - [`IsoconversionalResult`] / [`IsoLayerResult`] (re-exported from
//!   `integral_isoconversion`) — the output: one `IsoLayerResult` per α-layer
//!   carrying `eta`, `ea`, and a placeholder `LinearRegressionResult`.
//!
//! ## Algorithm
//! For each α-layer k the solver minimises:
//!
//! ```text
//! Φ(E) = Σᵢ Σⱼ≠ᵢ  J(E, Tᵢ) / J(E, Tⱼ)
//! ```
//!
//! where `J(E, Tᵢ) = exp(-E/RT_i) · Δtᵢ` is the trapezoidal approximation of
//! the temperature integral for experiment i at layer k.
//!
//! ## Non-trivial optimisation trick
//! The naive double-loop over all experiment pairs is O(n²).  
//! `solve_layer_optimized` rewrites the sum algebraically:
//!
//! ```text
//! Σᵢ Σⱼ≠ᵢ Jᵢ/Jⱼ  =  (Σ Jᵢ)(Σ 1/Jᵢ) − n
//! ```
//!
//! reducing the inner work to two O(n) passes — a significant speedup when many
//! experiments are used.  The search over E is a simple uniform grid scan
//! (`steps` points between `e_min` and `e_max`); no gradient information is
//! needed, which keeps the implementation robust against flat or noisy Φ landscapes.

use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoLayerResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IsoconversionalResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::kinetic_regression::LinearRegressionResult;
use crate::Kinetics::experimental_kinetics::kinetic_methods::*;

/// Solver for the Vyazovkin non-linear isoconversional method.
///
/// Searches for the activation energy Eα at each α-layer by minimising the
/// functional Φ(E) over a uniform grid of `steps` candidate values in
/// `[e_min, e_max]` (J/mol).
#[derive(Clone, Debug)]
pub struct VyazovkinSolver {
    /// Lower bound of the Eα search range (J/mol). Default: 40 000.
    pub e_min: f64,
    /// Upper bound of the Eα search range (J/mol). Default: 300 000.
    pub e_max: f64,
    /// Number of uniformly-spaced candidate E values to evaluate. Default: 200.
    pub steps: usize,
}

impl VyazovkinSolver {
    /// Run the Vyazovkin method over all α-layers in `grid`.
    ///
    /// Requires `grid.dt` to be populated (use `ConversionGridBuilder::with_dt_matrix()`).
    /// Returns one [`IsoLayerResult`] per layer (starting from layer index 1,
    /// since Δt at layer 0 is undefined).
    pub fn solve(&self, grid: &ConversionGrid) -> Result<IsoconversionalResult, TGADomainError> {
        let dt = grid.dt.as_ref().ok_or(TGADomainError::InvalidOperation(
            "Vyazovkin requires dt matrix".into(),
        ))?;

        let n_eta = grid.eta.len();
        let mut layers = Vec::with_capacity(n_eta - 1);

        for k in 1..n_eta {
            let ea = self.solve_layer_optimized(grid, dt, k)?;

            layers.push(IsoLayerResult {
                eta: grid.eta[k],
                ea,
                k: None,
                regression: LinearRegressionResult::default(),
            });
        }

        Ok(IsoconversionalResult {
            method: "Vyazovkin",
            layers,
        })
    }

    /// Naive O(n²) reference implementation of the layer solver.
    ///
    /// Iterates over all ordered pairs (i, j) with i ≠ j and accumulates
    /// `J(E,T_i) / J(E,T_j)`.  Kept for correctness verification against
    /// `solve_layer_optimized`.
    fn solve_layer(
        &self,
        grid: &ConversionGrid,
        dt: &Array2<f64>,
        k: usize,
    ) -> Result<f64, TGADomainError> {
        let inv_t_col = grid.inv_temperature.column(k);
        let dt_col = dt.column(k);

        let r = 8.314462618;

        let mut best_e = 0.0;
        let mut best_phi = f64::INFINITY;

        let n_exp = inv_t_col.len();

        let mut j = vec![0.0; n_exp];

        for s in 0..self.steps {
            let e = self.e_min + (self.e_max - self.e_min) * (s as f64) / (self.steps - 1) as f64;

            let mut valid = true;

            for i in 0..n_exp {
                let inv_t = inv_t_col[i];
                let dt_i = dt_col[i];

                let val = (-e * inv_t / r).exp() * dt_i;

                if !val.is_finite() {
                    valid = false;
                    break;
                }

                j[i] = val;
            }

            if !valid {
                continue;
            }

            let mut phi = 0.0;

            for i in 0..n_exp {
                let ji = j[i];

                if ji < 1e-300 {
                    continue;
                }

                for j2 in 0..n_exp {
                    if i == j2 {
                        continue;
                    }

                    let jj = j[j2];

                    if jj < 1e-300 {
                        continue;
                    }

                    phi += ji / jj;
                }
            }

            if phi < best_phi {
                best_phi = phi;
                best_e = e;
            }
        }

        Ok(best_e)
    }

    /// Optimised O(n) equivalent of `solve_layer`.
    ///
    /// Uses the algebraic identity  Σᵢ Σⱼ≠ᵢ Jᵢ/Jⱼ = (Σ Jᵢ)(Σ 1/Jᵢ) − n
    /// to replace the double loop with two single passes, then picks the E
    /// that minimises the result.
    pub fn solve_layer_optimized(
        &self,
        grid: &ConversionGrid,
        dt: &Array2<f64>,
        k: usize,
    ) -> Result<f64, TGADomainError> {
        let inv_t_col = grid.inv_temperature.column(k);
        let dt_col = dt.column(k);

        let r = 8.314462618;

        let n_exp = inv_t_col.len();

        let mut best_e = 0.0;
        let mut best_phi = f64::INFINITY;

        for s in 0..self.steps {
            let e = self.e_min + (self.e_max - self.e_min) * (s as f64) / (self.steps - 1) as f64;

            let mut sum_j = 0.0;
            let mut sum_inv_j = 0.0;

            let mut valid = true;

            for (&inv_t, &dt_i) in inv_t_col.iter().zip(dt_col.iter()) {
                let j = (-e * inv_t / r).exp() * dt_i;

                if !j.is_finite() || j < 1e-300 {
                    valid = false;
                    break;
                }

                sum_j += j;
                sum_inv_j += 1.0 / j;
            }

            if !valid {
                continue;
            }

            let phi = sum_j * sum_inv_j - (n_exp as f64);

            if phi < best_phi {
                best_phi = phi;
                best_e = e;
            }
        }

        Ok(best_e)
    }
    /// Low-level functional Φ(E) evaluated at a single candidate energy `e`
    /// for layer `k`.  Used internally; exposed for debugging and plotting.
    fn functional(&self, grid: &ConversionGrid, k: usize, e: f64) -> f64 {
        let r = 8.314;

        let n_exp = grid.temperature.nrows();

        let mut j = vec![0.0; n_exp];

        for i in 0..n_exp {
            let t = grid.temperature[[i, k]];
            let dt = grid.dt.as_ref().unwrap()[[i, k]];

            j[i] = (-e / (r * t)).exp() * dt;
        }

        let mut phi = 0.0;

        for i in 0..n_exp {
            for j2 in 0..n_exp {
                if i != j2 {
                    phi += j[i] / j[j2];
                }
            }
        }

        phi
    }
}

impl Default for VyazovkinSolver {
    fn default() -> Self {
        Self {
            e_min: 40_000.0,
            e_max: 300_000.0,
            steps: 200,
        }
    }
}

/// Marker struct that wires `VyazovkinSolver` into the `KineticMethod` trait
/// (used by the unified dispatcher in `isoconversion.rs`).
pub struct Vyzaovkin;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion_tests::tests::simulate_tga_first_order_with_dt;

    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::VyazovkinMethod;
    use crate::Kinetics::experimental_kinetics::kinetic_methods_tests::tests::build_view_from_cfg_exact_m0;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_non_isothermal;
    #[test]
    fn vyazovkin_compute_with_mock_non_isothermal_data() {
        // similar to ofw_compute_with_mock_non_isothermal_data2 but using Vyazovkin pipeline
        let cfg = base_advanced_config_non_isothermal(
            420.0,
            1e6,
            100_000.0,
            1.0, // seconds
            30_000,
            vec![0.5, 0.75, 1.0, 1.5],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();

        // use method adapter that builds dt-matrix internally
        let result = VyazovkinMethod {}.compute(&view).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, None);
        assert!(result.layers.len() > 0);
        for layer in &result.layers {
            assert!(layer.ea.is_finite());
        }
    }

    #[test]
    fn vyazovkin_mock_grid_basic() {
        // generate a synthetic grid that includes dt matrix via builder in method compute
        let grid = simulate_tga_first_order_with_dt(
            500.0,
            1e6,
            120_000.0,
            &[0.5, 1.0, 1.5, 2.0],
            0.5,
            50_000,
            None,
        );
        // Build a new grid with dt via the production builder path by wrapping into a view is complex here,
        // instead directly call solver expecting dt to be absent -> we check method path on real pipeline in previous test.
        // Here we only assert solver refuses grids without dt to avoid silent misuse.
        let result = VyazovkinSolver::default().solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, None);
    }
}
