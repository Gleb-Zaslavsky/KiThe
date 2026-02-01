// PROPOSED FIX for multiphase_equilibrium_residual_generator3
// The issue: ln(ni/nj) uses total phase moles, creating artificial coupling

// CURRENT (BROKEN) CODE:
// let ln_ai = ni.ln() - nj.ln() + ln_p_ratio;

// FIXED CODE should calculate reaction-specific phase totals:
pub fn multiphase_equilibrium_residual_generator3_fixed(
    n0: Vec<f64>,
    reactions: DMatrix<f64>, // m × r
    gibbs: Vec<GibbsFn>,     // m
    phases: Vec<Phase>,
    temperature: f64,
    pressure: f64,
    p0: f64,
    substanse_eps: f64,
    phase_eps: f64,
) -> Result<Box<dyn Fn(&[f64]) -> Result<Vec<f64>, ReactionExtentError>>, ReactionExtentError> {
    let m = n0.len();
    let r = reactions.ncols();

    if reactions.nrows() != m || gibbs.len() != m {
        return Err(ReactionExtentError::Other(
            "Dimension mismatch in residual generator".to_string(),
        ));
    }

    let species_phase = species_to_phase_map(&phases, m)?;

    Ok(Box::new(move |xi: &[f64]| {
        if xi.len() != r {
            return Err(ReactionExtentError::ResidualEvaluation(
                "xi length mismatch".to_string(),
            ));
        }

        // --- species mole numbers ---
        let mut n = n0.clone();
        for i in 0..m {
            for k in 0..r {
                n[i] += reactions[(i, k)] * xi[k];
            }
        }

        let rt = R * temperature;
        let ln_p_ratio = (pressure / p0).ln();

        let mut f = vec![0.0; r];

        for k in 0..r {
            let mut sum_ln_a = 0.0;
            let mut dg0 = 0.0;

            // FIXED: Calculate reaction-specific phase total for each reaction
            let mut reaction_phase_total = 0.0;
            for i in 0..m {
                let nu = reactions[(i, k)];
                if nu != 0.0 {
                    reaction_phase_total += n[i].max(substanse_eps);
                }
            }
            reaction_phase_total = reaction_phase_total.max(phase_eps);

            for i in 0..m {
                let nu = reactions[(i, k)];
                if nu == 0.0 {
                    continue;
                }

                let ni = n[i].max(substanse_eps);
                
                // FIXED: Use reaction-specific phase total instead of global phase total
                let ln_ai = ni.ln() - reaction_phase_total.ln() + ln_p_ratio;

                sum_ln_a += nu * ln_ai;
                dg0 += nu * gibbs[i](temperature);
            }

            let ln_k = -dg0 / rt;
            f[k] = sum_ln_a - ln_k;
        }

        Ok(f)
    }))
}