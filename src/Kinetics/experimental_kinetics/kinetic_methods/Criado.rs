use super::super::super::solid_state_kinetics_IVP::KineticModelNames;

use crate::Kinetics::experimental_kinetics::kinetic_methods::*;
/// ранжирование моделей.
#[derive(Clone, Debug)]
pub struct CriadoModelScore {
    pub model: KineticModelNames,
    pub error: f64,
}

#[derive(Clone, Debug)]
pub struct CriadoResult {
    pub experiment_id: String,
    pub ranking: Vec<CriadoModelScore>,
}

//========================================================================================================
//  SOLVER
//===============================================================================================================

#[derive(Clone, Copy, Debug)]
pub enum CriadoMode {
    Differential,
    IntegralIsothermal,
    IntegralNonIsothermal,
}

pub struct CriadoSolver {
    pub mode: CriadoMode,
    pub alpha_min: f64,
    pub alpha_max: f64,
    pub reference_alpha: f64,

    pub alpha_grid: Vec<f64>,

    pub models: Vec<KineticModelNames>,

    pub master_curves: Vec<Vec<f64>>,
    pub master_curves_iso: Vec<Vec<f64>>,

    pub master_curves_noniso: Vec<Vec<f64>>,
    pub result: Option<Vec<CriadoResult>>,
}

impl CriadoSolver {
    pub fn new() -> Self {
        let alpha_grid = (1..99).map(|i| i as f64 / 100.0).collect::<Vec<_>>();

        let models = criado_models();

        Self {
            mode: CriadoMode::IntegralIsothermal,

            alpha_min: 0.05,
            alpha_max: 0.95,
            reference_alpha: 0.5,
            alpha_grid,
            models,
            master_curves: Vec::new(),
            master_curves_iso: Vec::new(),
            master_curves_noniso: Vec::new(),
            result: None,
        }
    }
    pub fn reference_alpha(mut self, alpha: f64) -> Self {
        self.reference_alpha = alpha;
        self
    }
    //=========================================================================
    // ISOTHERMAL CASE
    //=========================================================================
    pub fn precompute_master_curves(&mut self) -> Result<(), TGADomainError> {
        self.master_curves.clear();

        for model in &self.models {
            let g = integrate_g2(model, &self.alpha_grid);

            let g_ref = interpolate(&self.alpha_grid, &g, self.reference_alpha)?;
            if !g_ref.is_finite() || g_ref == 0.0 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Criado reference integral is invalid for model {model:?}"
                )));
            }

            let curve: Vec<f64> = g.iter().map(|v| v / g_ref).collect();
            if curve.iter().any(|value| !value.is_finite()) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Criado produced non-finite master-curve values for model {model:?}"
                )));
            }

            self.master_curves.push(curve);
        }
        Ok(())
    }
    fn experimental_curve_isothermal(
        &self,
        exp: &ExperimentData,
    ) -> Result<Vec<(f64, f64)>, TGADomainError> {
        let alpha = &exp.conversion;
        let time = &exp.time;
        let ref_alpha = self.reference_alpha;
        let t_ref = interpolate(alpha, time, ref_alpha)?;
        if t_ref == 0.0 {
            return Err(TGADomainError::InvalidOperation(format!(
                "Criado reference time is zero for experiment '{}'",
                exp.meta.id
            )));
        }

        let mut curve = Vec::new();

        for (&a, &t) in alpha.iter().zip(time) {
            if a < self.alpha_min || a > self.alpha_max {
                continue;
            }

            let normalized_time = t / t_ref;
            if !normalized_time.is_finite() {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Criado produced a non-finite normalized time for experiment '{}'",
                    exp.meta.id
                )));
            }
            curve.push((a, normalized_time));
        }

        Ok(curve)
    }

    fn solve_experiment_isothermal(
        &self,
        exp: &ExperimentData,
    ) -> Result<Vec<CriadoModelScore>, TGADomainError> {
        let exp_curve = self.experimental_curve_isothermal(exp)?;
        if exp_curve.is_empty() {
            return Err(TGADomainError::InvalidConversionRange);
        }

        let alpha: Vec<f64> = exp_curve.iter().map(|(a, _)| *a).collect();

        let projected = project_models_to_experiment(self, &alpha)?;

        let mut scores = evaluate_models_fast(&exp_curve, &projected, &self.models);
        sort_model_scores(&mut scores, &exp.meta.id)?;
        Ok(scores)
    }

    pub fn solve_integral_isothermal(
        &mut self,
        data: &KineticDataView,
    ) -> Result<Vec<CriadoResult>, TGADomainError> {
        if self.master_curves.is_empty() {
            self.precompute_master_curves()?;
        }

        let mut results = Vec::new();

        for exp in &data.experiments {
            let ranking = self.solve_experiment_isothermal(exp)?;

            results.push(CriadoResult {
                experiment_id: exp.meta.id.clone(),
                ranking,
            });
        }
        self.result = Some(results.clone());
        Ok(results)
    }
    //===========================================================================================
    // INTEGRAL NONISOTHERMAL
    //==============================================================================================
    pub fn precompute_nonisothermal_curves(&mut self) -> Result<(), TGADomainError> {
        self.master_curves_noniso.clear();

        for model in &self.models {
            let g = integrate_g2(model, &self.alpha_grid);

            let mut z = Vec::new();

            let model_fn = model.get_fn();
            let model_fn = |a: f64| model_fn(a, vec![]);

            // The normalized Criado product g(a)*f(a) is identical for F1 and
            // the A2-A4 Avrami models, so this method cannot distinguish them.
            for (i, a) in self.alpha_grid.iter().enumerate() {
                let val = g[i] * model_fn(*a);

                z.push(val);
            }

            let z_ref = interpolate(&self.alpha_grid, &z, self.reference_alpha)?;
            if !z_ref.is_finite() || z_ref == 0.0 {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Criado reference value is invalid for model {model:?}"
                )));
            }

            for v in &mut z {
                *v /= z_ref;
            }

            self.master_curves_noniso.push(z);
        }
        Ok(())
    }

    fn solve_experiment_nonisothermal(
        &self,
        exp: &ExperimentData,
    ) -> Result<Vec<CriadoModelScore>, TGADomainError> {
        let exp_curve = experimental_z_curve(exp, self)?;
        if exp_curve.is_empty() {
            return Err(TGADomainError::InvalidConversionRange);
        }

        let mut scores = Vec::new();

        for (i, model) in self.models.iter().enumerate() {
            let err = evaluate_model_noniso(&exp_curve, self, i)?;

            scores.push(CriadoModelScore {
                model: model.clone(),
                error: err,
            });
        }

        if scores.iter().any(|score| !score.error.is_finite()) {
            return Err(TGADomainError::InvalidOperation(format!(
                "Criado produced a non-finite model score for experiment '{}'",
                exp.meta.id
            )));
        }
        scores.sort_by(|a, b| a.error.total_cmp(&b.error));

        Ok(scores)
    }

    fn solve_integral_nonisothermal(
        &mut self,
        data: &KineticDataView,
    ) -> Result<Vec<CriadoResult>, TGADomainError> {
        if self.master_curves_noniso.is_empty() {
            self.precompute_nonisothermal_curves()?;
        }

        let mut results = Vec::new();

        for exp in &data.experiments {
            let ranking = self.solve_experiment_nonisothermal(exp)?;

            results.push(CriadoResult {
                experiment_id: exp.meta.id.clone(),
                ranking,
            });
        }

        self.result = Some(results.clone());
        Ok(results)
    }

    //=============================================================================================
    // MAIN SOLVER
    //====================================================================================================
    pub fn solve(&mut self, data: &KineticDataView) -> Result<Vec<CriadoResult>, TGADomainError> {
        match self.mode {
            CriadoMode::IntegralIsothermal => self.solve_integral_isothermal(data),

            CriadoMode::IntegralNonIsothermal => self.solve_integral_nonisothermal(data),

            CriadoMode::Differential => self.solve_differential(data),
        }
    }

    //==========================================================================================

    //===================================================================================
    pub fn solve_differential(
        &mut self,
        _data: &KineticDataView,
    ) -> Result<Vec<CriadoResult>, TGADomainError> {
        Err(TGADomainError::InvalidOperation(
            "Criado differential mode not implemented yet".into(),
        ))
    }

    //==================================================================
    // INTERPRETATION
    //=================================================================
    pub fn global_criado_interpretation(
        &self,
    ) -> Result<CriadoGlobalInterpretation, TGADomainError> {
        if let Some(result) = self.result.clone() {
            let interp = global_criado_interpretation_local(&result);
            return Ok(interp);
        } else {
            return Err(TGADomainError::InvalidOperation("no result".to_string()));
        }
    }
}

pub fn criado_models() -> Vec<KineticModelNames> {
    vec![
        KineticModelNames::F1,
        KineticModelNames::F2,
        KineticModelNames::F3,
        KineticModelNames::R2,
        KineticModelNames::R3,
        KineticModelNames::D1,
        KineticModelNames::D2,
        KineticModelNames::D3,
        KineticModelNames::D4,
        KineticModelNames::A2,
        KineticModelNames::A3,
        KineticModelNames::A4,
    ]
}

fn integrate_g(f: &dyn Fn(f64) -> f64, alpha: &[f64]) -> Vec<f64> {
    // Composite Simpson integration of g(a) = \int_0^a da' / f(a')
    // This provides much better accuracy than a simple midpoint or trapezoidal rule
    // while keeping the implementation easy and robust for non-uniform grids.

    if alpha.is_empty() {
        return Vec::new();
    }

    let mut g = Vec::with_capacity(alpha.len());
    g.push(0.0);
    debug_assert!(alpha[0] >= 0.0);
    debug_assert!(alpha.windows(2).all(|w| w[0] < w[1]));
    let safe_inv = |a: f64| {
        let fa = f(a);
        if !fa.is_finite() || fa.abs() < 1e-15 {
            // Avoid division by zero / non-finite values by clamping.
            // In solid-state kinetics the integrand can be singular at a=0 for some models.
            1.0 / 1e-15
        } else {
            1.0 / fa
        }
    };

    for i in 1..alpha.len() {
        let a0 = alpha[i - 1];
        let a1 = alpha[i];
        let h = a1 - a0;
        let amid = 0.5 * (a0 + a1);

        let i0 = safe_inv(a0);
        let im = safe_inv(amid);
        let i1 = safe_inv(a1);

        // Simpson's rule over each interval
        let delta = h * (i0 + 4.0 * im + i1) / 6.0;
        g.push(g[i - 1] + delta);
    }

    g
}

fn integrate_g2(model: &KineticModelNames, alpha: &[f64]) -> Vec<f64> {
    if let Ok(analytical_integral) = model.get_g() {
        let g_result = alpha.iter().map(|&a| analytical_integral(a)).collect();
        g_result
    } else {
        // fallback to

        integrate_simpson(model, alpha)
    }
}

fn integrate_simpson(model: &KineticModelNames, alpha: &[f64]) -> Vec<f64> {
    // Fall back to the general numerical integrator (Simpson) using the model's function.
    let model_fn = model.get_fn();
    let f = |a: f64| model_fn(a, vec![]);
    integrate_g(&f, alpha)
}
#[allow(dead_code)]
fn evaluate_model_monotonic(
    exp_curve: &[(f64, f64)],
    solver: &CriadoSolver,
    model_idx: usize,
) -> f64 {
    let model_curve = &solver.master_curves[model_idx];

    let mut err = 0.0;

    let mut i = 0;

    for (a, v) in exp_curve {
        while i + 1 < solver.alpha_grid.len() && solver.alpha_grid[i + 1] < *a {
            i += 1;
        }

        let x0 = solver.alpha_grid[i];
        let x1 = solver.alpha_grid[i + 1];

        let y0 = model_curve[i];
        let y1 = model_curve[i + 1];

        let w = (a - x0) / (x1 - x0);

        let m = y0 + w * (y1 - y0);

        let d = v - m;

        err += d * d;
    }

    err.sqrt()
}

fn evaluate_models_fast(
    exp_curve: &[(f64, f64)],

    projected_models: &[Vec<f64>],

    models: &[KineticModelNames],
) -> Vec<CriadoModelScore> {
    let mut scores = Vec::new();

    for (i, model_curve) in projected_models.iter().enumerate() {
        let mut err = 0.0;

        for (j, (_, exp_v)) in exp_curve.iter().enumerate() {
            let d = exp_v - model_curve[j];

            err += d * d;
        }

        scores.push(CriadoModelScore {
            model: models[i].clone(),
            error: err.sqrt(),
        });
    }

    scores
}
pub fn interpolate(x: &[f64], y: &[f64], xq: f64) -> Result<f64, TGADomainError> {
    if x.len() != y.len() || x.len() < 2 {
        return Err(TGADomainError::InvalidOperation(format!(
            "Criado interpolation requires equally sized arrays with at least two points (x={}, y={})",
            x.len(),
            y.len()
        )));
    }
    if !xq.is_finite() || x.iter().chain(y.iter()).any(|value| !value.is_finite()) {
        return Err(TGADomainError::InvalidOperation(
            "Criado interpolation requires finite values".into(),
        ));
    }
    if let Some((index, pair)) = x.windows(2).enumerate().find(|(_, pair)| pair[0] > pair[1]) {
        return Err(TGADomainError::InvalidOperation(format!(
            "Criado interpolation requires non-decreasing X; x[{}]={} and x[{}]={}",
            index,
            pair[0],
            index + 1,
            pair[1]
        )));
    }

    // bounds
    if xq <= x[0] {
        return Ok(y[0]);
    }

    if xq >= x[x.len() - 1] {
        return Ok(y[y.len() - 1]);
    }

    // binary search
    let mut lo = 0;
    let mut hi = x.len() - 1;

    while hi - lo > 1 {
        let mid = (lo + hi) / 2;

        if x[mid] > xq {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    let x0 = x[lo];
    let x1 = x[hi];

    let y0 = y[lo];
    let y1 = y[hi];

    let w = (xq - x0) / (x1 - x0);

    let value = y0 + w * (y1 - y0);
    if !value.is_finite() {
        return Err(TGADomainError::InvalidOperation(
            "Criado interpolation produced a non-finite value".into(),
        ));
    }
    Ok(value)
}

pub fn find_alpha_index(alpha: &[f64], a: f64) -> usize {
    alpha
        .windows(2)
        .position(|w| w[0] <= a && w[1] >= a)
        .unwrap_or(alpha.len() - 2)
}

fn project_models_to_experiment(
    solver: &CriadoSolver,
    exp_alpha: &[f64],
) -> Result<Vec<Vec<f64>>, TGADomainError> {
    let mut projected = Vec::new();

    for curve in &solver.master_curves {
        let mut v = Vec::with_capacity(exp_alpha.len());

        for &a in exp_alpha {
            let m = interpolate(&solver.alpha_grid, curve, a)?;

            v.push(m);
        }

        projected.push(v);
    }

    Ok(projected)
}

fn evaluate_model_noniso(
    exp_curve: &[(f64, f64)],
    solver: &CriadoSolver,
    model_idx: usize,
) -> Result<f64, TGADomainError> {
    let model_curve = &solver.master_curves_noniso[model_idx];

    let mut err = 0.0;

    for (a, z_exp) in exp_curve {
        let z_model = interpolate(&solver.alpha_grid, model_curve, *a)?;

        let d = z_exp - z_model;

        err += d * d;
    }

    let score = err.sqrt();
    if !score.is_finite() {
        return Err(TGADomainError::InvalidOperation(
            "Criado produced a non-finite model score".into(),
        ));
    }
    Ok(score)
}

//=======================================================================================================================
//  NONISOTHERMAL
//==============================================================================================================

//Экспериментальный master plot

fn experimental_z_curve(
    exp: &ExperimentData,
    solver: &CriadoSolver,
) -> Result<Vec<(f64, f64)>, TGADomainError> {
    let alpha = &exp.conversion;
    if alpha.is_empty() {
        return Err(TGADomainError::EmptyColumn("conversion".into()));
    }

    let exp_min = alpha[0];
    let exp_max = alpha[alpha.len() - 1];

    if solver.reference_alpha < exp_min || solver.reference_alpha > exp_max {
        return Err(TGADomainError::InvalidConversionRange);
    }

    let eta_grid: Vec<f64> = solver
        .alpha_grid
        .iter()
        .copied()
        .filter(|&a| a >= solver.alpha_min && a <= solver.alpha_max && a >= exp_min && a <= exp_max)
        .collect();

    if eta_grid.is_empty() {
        return Err(TGADomainError::InvalidConversionRange);
    }

    let (_, temp, rate) = interpolate_linear(exp, &eta_grid).map_err(|error| {
        TGADomainError::InvalidOperation(format!(
            "Criado interpolation failed for experiment '{}': {error}",
            exp.meta.id
        ))
    })?;

    let ref_rate = interpolate(alpha, &exp.conversion_rate, solver.reference_alpha)?;
    let ref_temp = interpolate(alpha, &exp.temperature, solver.reference_alpha)?;
    let z_ref = ref_rate * ref_temp * ref_temp;

    if !z_ref.is_finite() || z_ref == 0.0 {
        return Err(TGADomainError::InvalidOperation(format!(
            "Criado reference rate/temperature is invalid for experiment '{}'",
            exp.meta.id
        )));
    }

    let mut curve = Vec::with_capacity(eta_grid.len());

    for i in 0..eta_grid.len() {
        let z = rate[i] * temp[i] * temp[i] / z_ref;
        if !z.is_finite() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Criado produced a non-finite point for experiment '{}' at conversion {}",
                exp.meta.id, eta_grid[i]
            )));
        }
        curve.push((eta_grid[i], z));
    }

    Ok(curve)
}

//=======================================================================================================================
//  AUTOMATED INTERPRETATION
//=========================================================================================================================
#[derive(Clone, Debug)]
pub struct CriadoGlobalInterpretation {
    pub mechanism: ReactionMechanism,

    pub best_model: KineticModelNames,

    pub confidence: f64,

    pub mechanism_votes: Vec<(ReactionMechanism, usize)>,
}

pub fn aggregate_mechanisms(
    results: &[CriadoResult],

    top_k: usize,
) -> Vec<(ReactionMechanism, usize)> {
    use std::collections::HashMap;

    let mut counts = HashMap::new();

    for res in results {
        for m in res.ranking.iter().take(top_k) {
            let mech = mechanism_of_model(m.model.clone());

            *counts.entry(mech).or_insert(0) += 1;
        }
    }

    let mut v: Vec<_> = counts.into_iter().collect();

    v.sort_by(|a, b| b.1.cmp(&a.1));

    v
}

pub fn global_criado_interpretation_local(results: &[CriadoResult]) -> CriadoGlobalInterpretation {
    let votes = aggregate_mechanisms(results, 3);

    let mechanism = votes
        .first()
        .map(|(m, _)| *m)
        .unwrap_or(ReactionMechanism::Unknown);

    let best_model = results[0].ranking[0].model.clone();

    let confidence = votes[0].1 as f64 / votes.iter().map(|v| v.1).sum::<usize>() as f64;

    CriadoGlobalInterpretation {
        mechanism,

        best_model,

        confidence,

        mechanism_votes: votes,
    }
}

pub fn best_model_frequency(results: &[CriadoResult]) -> Vec<(KineticModelNames, usize)> {
    use std::collections::HashMap;

    let mut counts = HashMap::new();

    for r in results {
        let Some(best) = r.ranking.first() else {
            continue;
        };
        let m = best.model.clone();

        *counts.entry(m).or_insert(0) += 1;
    }

    let mut v: Vec<_> = counts.into_iter().collect();

    v.sort_by(|a, b| b.1.cmp(&a.1));

    v
}

fn sort_model_scores(
    scores: &mut [CriadoModelScore],
    experiment_id: &str,
) -> Result<(), TGADomainError> {
    if scores.is_empty() {
        return Err(TGADomainError::InvalidOperation(format!(
            "Criado produced no model scores for experiment '{experiment_id}'"
        )));
    }
    if scores.iter().any(|score| !score.error.is_finite()) {
        return Err(TGADomainError::InvalidOperation(format!(
            "Criado produced a non-finite model score for experiment '{experiment_id}'"
        )));
    }
    scores.sort_by(|a, b| a.error.total_cmp(&b.error));
    Ok(())
}
/*
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReactionMechanism {

    ReactionOrder,

    ContractingGeometry,

    Diffusion,

    NucleationGrowth,

    PowerLaw,

    Unknown,
}

*/

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ReactionMechanism {
    ReactionOrder,

    ContractingGeometry,

    Diffusion,

    NucleationGrowth,

    PowerLaw,

    Unknown,
}

pub fn mechanism_of_model(model: KineticModelNames) -> ReactionMechanism {
    match model {
        KineticModelNames::F1 | KineticModelNames::F2 | KineticModelNames::F3 => {
            ReactionMechanism::ReactionOrder
        }

        KineticModelNames::R2 | KineticModelNames::R3 => ReactionMechanism::ContractingGeometry,

        KineticModelNames::D1
        | KineticModelNames::D2
        | KineticModelNames::D3
        | KineticModelNames::D4 => ReactionMechanism::Diffusion,

        KineticModelNames::A2 | KineticModelNames::A3 | KineticModelNames::A4 => {
            ReactionMechanism::NucleationGrowth
        }

        KineticModelNames::P2 | KineticModelNames::P3 => ReactionMechanism::PowerLaw,

        _ => ReactionMechanism::Unknown,
    }
}
/*
Идея простая:
confidence =
error_best / error_second
Если
error_second ≫ error_best

*/

pub struct CriadoInterpretation {
    pub best_model: KineticModelNames,

    pub mechanism: ReactionMechanism,

    pub confidence: f64,
}

pub fn interpret_criado(ranking: &[CriadoModelScore]) -> CriadoInterpretation {
    let best = &ranking[0];

    let second = &ranking[1];

    let mechanism = mechanism_of_model(best.model.clone());

    let confidence = second.error / best.error;

    CriadoInterpretation {
        best_model: best.model.clone(),

        mechanism,

        confidence,
    }
}

pub fn mechanism_vote(ranking: &[CriadoModelScore], top_k: usize) -> ReactionMechanism {
    use std::collections::HashMap;

    let mut counts = HashMap::new();

    for r in ranking.iter().take(top_k) {
        let mech = mechanism_of_model(r.model.clone());

        *counts.entry(mech).or_insert(0) += 1;
    }

    counts
        .into_iter()
        .max_by_key(|(_, v)| *v)
        .map(|(m, _)| m)
        .unwrap_or(ReactionMechanism::Unknown)
}

pub fn interpret_criado_robust(ranking: &[CriadoModelScore]) -> CriadoInterpretation {
    let best = &ranking[0].model;

    let mechanism = mechanism_vote(ranking, 3);

    let confidence = ranking[1].error / ranking[0].error;

    CriadoInterpretation {
        best_model: best.clone(),

        mechanism,

        confidence,
    }
}
//==============================================================================================================
//TESTING
//===========================================================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn criado_score_sort_rejects_non_finite_scores_without_panicking() {
        let mut scores = vec![CriadoModelScore {
            model: KineticModelNames::F1,
            error: f64::NAN,
        }];
        assert!(sort_model_scores(&mut scores, "invalid").is_err());
    }

    #[test]
    fn best_model_frequency_ignores_empty_rankings() {
        let results = vec![CriadoResult {
            experiment_id: "empty".to_string(),
            ranking: Vec::new(),
        }];
        assert!(best_model_frequency(&results).is_empty());
    }

    #[test]
    fn criado_interpolation_validates_inputs_and_preserves_clamped_bounds() {
        assert_eq!(interpolate(&[0.0, 1.0], &[2.0, 4.0], 0.5).unwrap(), 3.0);
        assert_eq!(interpolate(&[0.0, 1.0], &[2.0, 4.0], -1.0).unwrap(), 2.0);
        assert!(interpolate(&[0.0], &[2.0], 0.5).is_err());
        assert_eq!(interpolate(&[0.0, 0.0], &[2.0, 4.0], 0.5).unwrap(), 4.0);
        assert!(interpolate(&[0.0, 1.0], &[2.0, f64::INFINITY], 0.5).is_err());
    }

    use crate::Kinetics::experimental_kinetics::kinetic_methods_tests::tests::build_view_from_cfg_exact_m0;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_isothermal;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::base_advanced_config_non_isothermal;
    use std::time::Instant;

    // ── helpers ──────────────────────────────────────────────────────────────

    /// Temperatures used across most tests (K).
    fn iso_temps() -> Vec<f64> {
        vec![520.0, 540.0, 560.0, 580.0, 600.0, 620.0]
    }

    fn assert_f1_in_nonisothermal_best_tie(ranking: &[CriadoModelScore]) {
        let best_error = ranking.first().expect("ranking must not be empty").error;
        let f1_error = ranking
            .iter()
            .find(|score| score.model == KineticModelNames::F1)
            .expect("F1 must be included in Criado candidates")
            .error;
        let tolerance = 1e-10 * best_error.abs().max(f1_error.abs()).max(1e-12) + 1e-14;

        assert!(
            f1_error <= best_error + tolerance,
            "F1 should be tied with the best non-isothermal model: best={best_error:.16e}, F1={f1_error:.16e}"
        );
        for equivalent_model in [
            KineticModelNames::A2,
            KineticModelNames::A3,
            KineticModelNames::A4,
        ] {
            let error = ranking
                .iter()
                .find(|score| score.model == equivalent_model)
                .expect("all Avrami candidates must be included")
                .error;
            assert!(
                (error - f1_error).abs() <= tolerance,
                "{equivalent_model:?} should be indistinguishable from F1 in the normalized non-isothermal Criado plot: F1={f1_error:.16e}, model={error:.16e}"
            );
        }
    }

    #[test]
    fn criado_integral_compute_with_mock_isothermal_data() {
        let now = Instant::now();
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, iso_temps(), 0.1, 10_000);
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = CriadoSolver::new()
            .solve_integral_isothermal(&view)
            .unwrap();
        for res in &result {
            let winner = res.ranking.first().unwrap();
            let winner = winner.model.clone();
            let expected_winner = KineticModelNames::F1;
            assert!(
                winner == expected_winner,
                "Expected winner: {:?}, got: {:?}",
                expected_winner,
                winner
            );

            let interp = interpret_criado(&res.ranking);

            println!(
                "Best model: {:?}, mechanism: {:?}, confidence: {:.2}",
                interp.best_model, interp.mechanism, interp.confidence,
            );
        }

        // println!("{:#?}", result);
        println!("Time elapsed: {:?}", now.elapsed());
    }

    #[test]
    fn criado_integral_compute_with_mock_isothermal_data2() {
        let now = Instant::now();
        let cfg = base_advanced_config_isothermal(
            1e6,
            100_000.0,
            vec![520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0],
            1.0,
            30_000,
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let result = CriadoSolver::new()
            .solve_integral_isothermal(&view)
            .unwrap();
        for res in &result {
            let winner = res.ranking.first().unwrap();
            let winner = winner.model.clone();
            let expected_winner = KineticModelNames::F1;
            assert!(
                winner == expected_winner,
                "Expected winner: {:?}, got: {:?}",
                expected_winner,
                winner
            );

            let interp = interpret_criado(&res.ranking);

            println!(
                "Best model: {:?}, mechanism: {:?}, confidence: {:.2}",
                interp.best_model, interp.mechanism, interp.confidence,
            );
        }
        //println!("{:#?}", result);
        println!("Time elapsed: {:?}", now.elapsed());
    }

    #[test]
    fn criado_integral_compute_with_mock_isothermal_data3() {
        let now = Instant::now();
        let cfg = base_advanced_config_isothermal(
            1e7,
            130_000.0,
            vec![600.0, 625.0, 650.0, 675.0, 700.0, 725.0],
            0.5,
            150_000,
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let mut solver = CriadoSolver::new();
        let result = solver.solve_integral_isothermal(&view).unwrap();
        // println!("{:#?}", result);
        println!("Time elapsed: {:?}", now.elapsed());
        for res in &result {
            let winner = res.ranking.first().unwrap();
            let winner = winner.model.clone();
            let expected_winner = KineticModelNames::F1;
            assert!(
                winner == expected_winner,
                "Expected winner: {:?}, got: {:?}",
                expected_winner,
                winner
            );

            let interp = interpret_criado(&res.ranking);

            println!(
                "Best model: {:?}, mechanism: {:?}, confidence: {:.2}",
                interp.best_model, interp.mechanism, interp.confidence,
            );
        }
        let interp = solver.global_criado_interpretation().unwrap();
        println!("interp {:?}", interp);
    }

    //=================================================================================
    // NON-ISOTHERMAL
    //=================================================================================
    #[test]
    fn criado_reports_experiment_outside_reference_conversion_range() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            140_000.0,
            0.1,
            100_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();

        let mut solver = CriadoSolver::new();
        let result = solver.solve_integral_nonisothermal(&view);
        assert!(matches!(
            result,
            Err(TGADomainError::InvalidConversionRange)
        ));
    }

    #[test]
    fn criado_compute_with_mock_non_isothermal_data2() {
        let now = Instant::now();
        let cfg = base_advanced_config_non_isothermal(
            220.0,
            1e6,
            70_000.0,
            1.0,
            30_000,
            vec![0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let mut solver = CriadoSolver::new();
        let result = solver.solve_integral_nonisothermal(&view).unwrap();
        println!("Time elapsed: {:?}", now.elapsed());
        println!("{:#?}", result);
        for res in &result {
            assert_f1_in_nonisothermal_best_tie(&res.ranking);

            let interp = interpret_criado(&res.ranking);

            println!(
                "Best model: {:?}, mechanism: {:?}, confidence: {:.2}",
                interp.best_model, interp.mechanism, interp.confidence,
            );
        }
        let _interp = solver.global_criado_interpretation().unwrap();
    }

    #[test]
    fn criado_compute_with_mock_non_isothermal_data3() {
        let now = Instant::now();
        let cfg = base_advanced_config_non_isothermal(
            600.0,
            1e7,
            150_000.0,
            0.5,
            50_000,
            vec![0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        let solver = CriadoSolver::new();
        let mut solver = solver.reference_alpha(0.15);
        let result = solver.solve_integral_nonisothermal(&view).unwrap();
        println!("Time elapsed: {:?}", now.elapsed());
        println!("{:#?}", result);
        for res in &result {
            assert_f1_in_nonisothermal_best_tie(&res.ranking);

            let interp = interpret_criado(&res.ranking);

            println!(
                "Best model: {:?}, mechanism: {:?}, confidence: {:.2}",
                interp.best_model, interp.mechanism, interp.confidence,
            );
        }
        let _interp = solver.global_criado_interpretation().unwrap();
    }
}
