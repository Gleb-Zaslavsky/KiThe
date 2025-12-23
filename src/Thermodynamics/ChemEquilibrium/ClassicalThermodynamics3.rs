use std::f64;

use nalgebra::{DMatrix, SVD};

/// Result of reaction discovery from element matrix
#[derive(Debug, Clone)]
pub struct ReactionBasis {
    /// Rank of element matrix A
    pub rank: usize,
    /// Number of independent reactions (nullspace dimension)
    pub num_reactions: usize,
    /// Reaction stoichiometry matrix N (m x r)
    /// Each column is one independent reaction
    pub reactions: DMatrix<f64>,
}

/// Computes a basis of independent chemical reactions from an element
/// composition matrix.
///
/// # Mathematical meaning
///
/// Let `A` be the element composition matrix with:
/// - rows corresponding to chemical species
/// - columns corresponding to chemical elements
/// - `A[i, e]` equal to the number of atoms of element `e` in species `i`
///
/// Any physically admissible chemical reaction must conserve each element:
///
/// ```text
/// Aᵀ · ν = 0
/// ```
///
/// where `ν` is the vector of stoichiometric coefficients of a reaction.
/// Therefore, all independent reactions lie in the nullspace of `Aᵀ`.
///
/// This function computes a basis of that nullspace using singular value
/// decomposition (SVD).
///
/// # Physical meaning
///
/// - Each column of the returned matrix represents one independent chemical
///   reaction.
/// - The sign and scaling of each reaction vector are arbitrary and have no
///   physical significance.
/// - Together, these reactions span all possible composition changes that
///   conserve elements.
///
/// The number of independent reactions is:
///
/// ```text
/// r = number_of_species − rank(A)
/// ```
///
/// This construction guarantees that:
/// - Element conservation is satisfied automatically
/// - No redundant or dependent reactions are introduced
/// - Equilibrium systems formulated in reaction extents are well-conditioned
///
/// # Arguments
///
/// * `a`   – Element composition matrix `(m × E)`
/// * `tol` – Singular value threshold used to determine numerical rank
///
/// # Returns
///
/// A [`ReactionBasis`] containing:
/// - the rank of the element matrix
/// - the number of independent reactions
/// - a reaction stoichiometry matrix `(m × r)`
///
/// # Notes
///
/// This function does **not** classify species as reactants or products.
/// That distinction emerges naturally from the solution of the equilibrium
/// equations.
///
/// # Panics
///
/// Panics if SVD decomposition fails (which should not happen for finite
/// matrices).

pub fn compute_reaction_basis(a: &DMatrix<f64>, tol: f64) -> ReactionBasis {
    let m = a.nrows(); // species
    let _e = a.ncols(); // elements

    // Compute SVD of A^T (E x m)
    let at = a.transpose();
    let svd = SVD::new(at.clone(), true, true);

    let _u = svd.u.expect("SVD failed: U");
    let v_t = svd.v_t.expect("SVD failed: Vt");
    let singular_values = svd.singular_values;

    // Rank = number of singular values above tolerance
    let rank = singular_values.iter().filter(|&&s| s > tol).count();

    // Nullspace dimension
    let num_reactions = if rank <= m { m - rank } else { 0 };

    // Nullspace basis vectors are columns of V corresponding
    // to zero singular values
    //
    // V is m x m, stored as V^T => rows of V^T
    let mut reactions = DMatrix::<f64>::zeros(m, num_reactions);

    // Only proceed if we have reactions and valid dimensions
    if num_reactions > 0 && rank < v_t.nrows() {
        let mut col = 0;
        for i in rank..std::cmp::min(rank + num_reactions, v_t.nrows()) {
            // ith row of V^T is ith column of V
            let v_col = v_t.row(i).transpose();
            if col < reactions.ncols() {
                reactions.set_column(col, &v_col);
                col += 1;
            }
        }
    }

    ReactionBasis {
        rank,
        num_reactions,
        reactions,
    }
}

pub fn mole_numbers(xi: &[f64], n0: &[f64], reactions: &DMatrix<f64>) -> Vec<f64> {
    let m = n0.len();
    let r = xi.len();

    let mut n = n0.to_vec();

    for i in 0..m {
        for j in 0..r {
            n[i] += reactions[(i, j)] * xi[j];
        }
    }

    n
}
pub fn reaction_delta_g(reactions: &DMatrix<f64>, gibbs: &[GibbsFn]) -> Vec<GibbsFn> {
    let m = reactions.nrows();
    let r = reactions.ncols();

    assert_eq!(gibbs.len(), m);

    (0..r)
        .map(|j| {
            // Copy stoichiometry column (cheap, numeric)
            let nu = reactions.column(j).clone_owned();
            let g = gibbs.to_vec(); // Arc is Clone — cheap and safe

            Arc::new(move |t: f64| {
                let mut dg = 0.0;
                for i in 0..m {
                    dg += nu[i] * g[i](t);
                }
                dg
            }) as GibbsFn
        })
        .collect()
}
/// Standard-state Gibbs free energy [J/mol] as a function of temperature
use std::sync::Arc;

pub type GibbsFn = Arc<dyn Fn(f64) -> f64 + Send + Sync>;

#[allow(non_upper_case_globals)]
const R_g: f64 = 8.314;
pub fn equilibrium_residual_generator(
    n0: Vec<f64>,
    reactions: DMatrix<f64>,
    delta_g: Vec<GibbsFn>,
    temperature: f64,
    pressure: f64,
) -> impl Fn(&[f64]) -> Vec<f64> {
    let r = reactions.ncols();
    let m = reactions.nrows();

    move |xi: &[f64]| {
        let n = mole_numbers(xi, &n0, &reactions);
        let N_tot: f64 = n.iter().sum();

        let mut f = vec![0.0; r];

        for j in 0..r {
            let mut sum = 0.0;
            let mut dN = 0.0;

            for i in 0..m {
                sum += reactions[(i, j)] * n[i].ln();
                dN += reactions[(i, j)];
            }

            let ln_Kp = -delta_g[j](temperature) / (R_g * temperature);

            f[j] = sum - dN * N_tot.ln() + dN * (pressure / 101325.0).ln() - ln_Kp;
        }

        f
    }
}

////////////////////////////////////TESTS////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    const TOL: f64 = 1e-10;

    #[test]
    fn diatomic_dissociation_o2_o() {
        // Species: O2, O
        // Elements: O
        //
        // Reaction: O2 ⇌ 2 O

        let a = DMatrix::from_row_slice(
            2,
            1,
            &[
                2.0, // O2
                1.0, // O
            ],
        );

        let basis = compute_reaction_basis(&a, TOL);
        println!("basis: {:?}", basis);
        assert_eq!(basis.rank, 1);
        assert_eq!(basis.num_reactions, 1);

        // Check element conservation: Aᵀ * ν = 0
        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }

    #[test]
    fn nitrogen_oxygen_system() {
        // Species: N2, O2, NO, NO2
        // Elements: N, O

        let a = DMatrix::from_row_slice(
            4,
            2,
            &[
                2.0, 0.0, // N2
                0.0, 2.0, // O2
                1.0, 1.0, // NO
                1.0, 2.0, // NO2
            ],
        );

        let basis = compute_reaction_basis(&a, TOL);

        assert_eq!(basis.rank, 2);
        assert_eq!(basis.num_reactions, 2);

        // Check nullspace condition
        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }

    #[test]
    fn methane_combustion_species() {
        // Species: CH4, O2, CO2, H2O
        // Elements: C, H, O

        let a = DMatrix::from_row_slice(
            4,
            3,
            &[
                1.0, 4.0, 0.0, // CH4
                0.0, 0.0, 2.0, // O2
                1.0, 0.0, 2.0, // CO2
                0.0, 2.0, 1.0, // H2O
            ],
        );

        let basis = compute_reaction_basis(&a, TOL);

        assert_eq!(basis.rank, 3);
        assert_eq!(basis.num_reactions, 1);

        let residual = a.transpose() * &basis.reactions;
        assert!(residual.norm() < 1e-12);
    }
}
