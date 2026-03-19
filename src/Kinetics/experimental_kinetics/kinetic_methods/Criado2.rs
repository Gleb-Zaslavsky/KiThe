use super::super::super::solid_state_kinetics_IVP::KineticModelNames;
use crate::Kinetics::experimental_kinetics::kinetic_methods::Criado::*;
use crate::Kinetics::experimental_kinetics::kinetic_methods::*;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnNature;

/*



*/

pub struct CKAGuess {
    pub n: f64,
    pub m: f64,
    pub bounds_n: (f64, f64),
    pub bounds_m: (f64, f64),
}

pub fn guess_from_model(model: KineticModelNames) -> CKAGuess {
    match model {
        KineticModelNames::F1 => CKAGuess {
            n: 1.0,
            m: 0.0,
            bounds_n: (0.5, 1.5),
            bounds_m: (0.0, 0.2),
        },

        KineticModelNames::F2 => CKAGuess {
            n: 1.0,
            m: 1.0,
            bounds_n: (0.5, 1.5),
            bounds_m: (0.5, 1.5),
        },

        KineticModelNames::A2 => CKAGuess {
            n: 1.0,
            m: 0.5,
            bounds_n: (0.5, 1.5),
            bounds_m: (0.2, 1.0),
        },

        KineticModelNames::A3 => CKAGuess {
            n: 1.0,
            m: 0.66,
            bounds_n: (0.5, 1.5),
            bounds_m: (0.3, 1.2),
        },

        KineticModelNames::R2 => CKAGuess {
            n: 0.5,
            m: 0.0,
            bounds_n: (0.2, 1.0),
            bounds_m: (0.0, 0.2),
        },

        KineticModelNames::R3 => CKAGuess {
            n: 0.66,
            m: 0.0,
            bounds_n: (0.3, 1.2),
            bounds_m: (0.0, 0.2),
        },

        _ => CKAGuess {
            n: 1.0,
            m: 1.0,
            bounds_n: (0.0, 3.0),
            bounds_m: (0.0, 3.0),
        },
    }
}

pub fn guess_from_criado(criado: &CriadoResult) -> CKAGuess {
    let best_model = criado.ranking[0].model;

    guess_from_model(best_model)
}

pub fn blended_guess(criado: &CriadoResult) -> CKAGuess {
    let mut n = 0.0;
    let mut m = 0.0;
    let mut w_sum = 0.0;

    for (i, model_score) in criado.ranking.iter().take(3).enumerate() {
        let w = 1.0 / (i as f64 + 1.0);

        let g = guess_from_model(model_score.model);

        n += g.n * w;
        m += g.m * w;

        w_sum += w;
    }

    CKAGuess {
        n: n / w_sum,
        m: m / w_sum,
        bounds_n: (0.0, 3.0),
        bounds_m: (0.0, 3.0),
    }
}

//==================================================================================================
// KineticMethod impl (wrapper)
//==================================================================================================
#[derive(Clone, Debug)]
pub struct CriadoKineticMethod {
    pub mode: CriadoMode,
    pub alpha_min: f64,
    pub alpha_max: f64,
    pub reference_alpha: f64,
}

impl Default for CriadoKineticMethod {
    fn default() -> Self {
        let solver = CriadoSolver::new();
        Self {
            mode: solver.mode,
            alpha_min: solver.alpha_min,
            alpha_max: solver.alpha_max,
            reference_alpha: solver.reference_alpha,
        }
    }
}

impl KineticMethod for CriadoKineticMethod {
    type Output = Vec<CriadoResult>;

    fn name(&self) -> &'static str {
        "Criado"
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let mut solver = CriadoSolver::new();
        solver.mode = self.mode;
        solver.alpha_min = self.alpha_min;
        solver.alpha_max = self.alpha_max;
        solver.reference_alpha = self.reference_alpha;
        solver.solve(data)
    }

    fn requirements(&self) -> KineticRequirements {
        let needs_rate = matches!(
            self.mode,
            CriadoMode::IntegralNonIsothermal | CriadoMode::Differential
        );
        let needs_isothermal = matches!(self.mode, CriadoMode::IntegralIsothermal);
        KineticRequirements {
            min_experiments: 1,
            needs_conversion: true,
            needs_conversion_rate: needs_rate,
            needs_temperature: needs_isothermal,
            needs_heating_rate: false,
        }
    }

    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        match self.mode {
            CriadoMode::IntegralIsothermal => vec![ColumnNature::Conversion, ColumnNature::Time],
            CriadoMode::IntegralNonIsothermal | CriadoMode::Differential => vec![
                ColumnNature::Conversion,
                ColumnNature::Temperature,
                ColumnNature::Time,
                ColumnNature::ConversionRate,
            ],
        }
    }
}
