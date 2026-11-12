use RustedSciThe::symbolic::symbolic_engine::Expr;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
use tabled::{Table, Tabled};

#[allow(non_upper_case_globals)]
pub const two: Expr = Expr::Const(2.0);
#[allow(non_upper_case_globals)]
pub const three: Expr = Expr::Const(3.0);
#[allow(non_upper_case_globals)]
pub const four: Expr = Expr::Const(4.0);
#[allow(non_upper_case_globals)]
pub const one: Expr = Expr::Const(1.0);
#[allow(non_upper_case_globals)]
pub const half: Expr = Expr::Const(0.5);
#[allow(non_upper_case_globals)]
pub const neg: Expr = Expr::Const(-1.0);
/// 2*(1-a)*(-ln(1-a))^(1/2)
pub fn A2() -> Expr {
    let a = Expr::Var("a".to_owned());
    let f = two * (one - a.clone()) * Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(half));
    return f;
}
/// 3*(1-a)*(-ln(1-a))^(2/3)
pub fn A3() -> Expr {
    let a = Expr::Var("a".to_owned());
    let f = three
        * (one - a.clone())
        * Expr::Pow(
            Box::new(-Expr::ln(one - a)),
            Box::new(Expr::Const(2.0 / 3.0)),
        );
    return f;
}
///  4*(1-a)*(-lm(1-a))^(3/4)
pub fn A4() -> Expr {
    let a = Expr::Var("a".to_owned());
    let f = four
        * (one - a.clone())
        * Expr::Pow(
            Box::new(-Expr::ln(one - a)),
            Box::new(Expr::Const(3.0 / 4.0)),
        );
    return f;
}
/// (1-a)^0.5
pub fn R2() -> Expr {
    let a = Expr::Var("a".to_owned());

    return Expr::Pow(Box::new(one - a), Box::new(half));
}

/// (1-a)^(2/3)
pub fn R3() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(one - a), Box::new(Expr::Const(2.0 / 3.0)));
}
/// a^-1
pub fn D1() -> Expr {
    let a = Expr::Var("a".to_owned());

    return Expr::Pow(Box::new(a), Box::new(neg));
}
/// (-ln(1-a))^-1
pub fn D2() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(neg));
}
/// (3/2)*(1-a)^(2/3)*(1-(1-a)^(1/3))**-1  
pub fn D3() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Const(1.5)
        * Expr::Pow(Box::new(one - a.clone()), Box::new(Expr::Const(2.0 / 3.0)))
        * Expr::Pow(
            Box::new(one - Expr::Pow(Box::new(one - a), Box::new(Expr::Const(1.0 / 3.0)))),
            Box::new(neg),
        );
}
/// (3/2)*((1-a)^(-1/3)-1)^-1
pub fn D4() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Const(1.5)
        * Expr::Pow(
            Box::new(Expr::Pow(Box::new(one - a), Box::new(Expr::Const(-1.0 / 3.0))) - one),
            Box::new(neg),
        );
}
///  a^(-1/2)
pub fn P2_3() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(a), Box::new(Expr::Const(-0.5)));
}
/// a^(1/2)
pub fn P2() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(a), Box::new(half));
}
/// a^(3/2)
pub fn P3() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(a), Box::new(Expr::Const(0.75)));
}
/// 1-a
pub fn F1() -> Expr {
    let a = Expr::Var("a".to_owned());
    return one - a;
}
/// (1 - a)^2
pub fn F2() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(one - a), Box::new(two));
}
/// (1-a)^3
pub fn F3() -> Expr {
    let a = Expr::Var("a".to_owned());
    return Expr::Pow(Box::new(one - a), Box::new(three));
}

/// models with parameters
///  (1-a)^m*a^n*(-ln(1-a))^p
pub fn Sestak_Berggen(m: f64, n: f64, p: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    let p = Expr::Const(p);
    return Expr::Pow(Box::new(one - a.clone()), Box::new(m))
        * Expr::Pow(Box::new(a.clone()), Box::new(n))
        * Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(p));
}
/// sestak bergren model jacobian  log(1-a)^p*(n*(1-a)^m*a^(n-1) - m*(1-a)^(m-1)*a^n) - p*log(1-a)^(p-1)*(1-a)^(m-1)*a^n
pub fn SB_j(m: f64, n: f64, p: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    let p = Expr::Const(p);
    let remaining = one - a.clone();
    let logarithmic_conversion = -Expr::ln(remaining.clone());
    return Expr::Pow(
        Box::new(logarithmic_conversion.clone()),
        Box::new(p.clone()),
    ) * (n.clone()
        * Expr::Pow(Box::new(remaining.clone()), Box::new(m.clone()))
        * Expr::Pow(Box::new(a.clone()), Box::new(n.clone() - one))
        - m.clone()
            * Expr::Pow(Box::new(remaining.clone()), Box::new(m.clone() - one))
            * Expr::Pow(Box::new(a.clone()), Box::new(n.clone())))
        + p.clone()
            * Expr::Pow(Box::new(logarithmic_conversion), Box::new(p - one))
            * Expr::Pow(Box::new(remaining), Box::new(m - one))
            * Expr::Pow(Box::new(a), Box::new(n));
}
///  m*(1-a)*(-ln(1-a))^(1-1/m)
pub fn Johnson_Mehl_Avrami(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    return m.clone()
        * (one - a.clone())
        * Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(one - one / m));
}
/// m*(-ln(1-a)^(-1/m)/m - log(1-a)^(-1/m) - (-ln(1-a))**(1-1/m))
pub fn JMA_jac(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let logarithmic_conversion = -Expr::ln(one - a);
    return (m.clone() - one)
        * Expr::Pow(
            Box::new(logarithmic_conversion.clone()),
            Box::new(-one / m.clone()),
        )
        - m.clone() * Expr::Pow(Box::new(logarithmic_conversion), Box::new(one - one / m));
}
///  m*a^(1-1/m)
pub fn Accelerating_model(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    return m.clone() * Expr::Pow(Box::new(a), Box::new(one - one / m));
}
/// (1-1/m)*m*a**(-1/m)
pub fn Ac_jac(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    return (one - one / m.clone()) * m.clone() * Expr::Pow(Box::new(a), Box::new(-one / m));
}
/// (1-a)^m
pub fn Decelerating_model(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    return Expr::Pow(Box::new(one - a), Box::new(m));
}
///  -m*(1-a)^(m-1)
pub fn Dec_jac(m: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    return -m.clone() * Expr::Pow(Box::new(one - a), Box::new(m - one));
}
/// (1-a)^m * a^n
pub fn Prout_Tompkins_ext(m: f64, n: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    return Expr::Pow(Box::new(one - a.clone()), Box::new(m)) * Expr::Pow(Box::new(a), Box::new(n));
}
/// n*(1-a)^m * a^(n-1) - m*(1-a)^(m-1) * a^n
pub fn PTe_jac(m: f64, n: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    return n.clone()
        * Expr::Pow(Box::new(one - a.clone()), Box::new(m.clone()))
        * Expr::Pow(Box::new(a.clone()), Box::new(n.clone() - one))
        - m.clone()
            * Expr::Pow(Box::new(one - a.clone()), Box::new(m - one))
            * Expr::Pow(Box::new(a), Box::new(n));
}
/// c*(1-a)**m*a**n
pub fn Sestak_Berggren_trunc_with_param(m: f64, n: f64, c: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    let c = Expr::Const(c);
    return c
        * Expr::Pow(Box::new(one - a.clone()), Box::new(m))
        * Expr::Pow(Box::new(a), Box::new(n));
}
/// c*(n*(1-a)**m * a**(n-1) - m*(1-a)**(m-1) * a**n)
pub fn SBtp_jac(m: f64, n: f64, c: f64) -> Expr {
    let a = Expr::Var("a".to_owned());
    let m = Expr::Const(m);
    let n = Expr::Const(n);
    let c = Expr::Const(c);
    return c
        * (n.clone()
            * Expr::Pow(Box::new(one - a.clone()), Box::new(m.clone()))
            * Expr::Pow(Box::new(a.clone()), Box::new(n.clone() - one))
            - m.clone()
                * Expr::Pow(Box::new(one - a.clone()), Box::new(m - one))
                * Expr::Pow(Box::new(a), Box::new(n)));
}

/// Enumeration of available solid-state kinetic model names.
///
/// This enum represents different kinetic models used in solid-state reactions,
/// including nucleation and growth models (A2-A4), contracting geometry models (R2-R3),
/// diffusion models (D1-D4), power law models (P2-P3), and reaction order models (F1-F3).
///
/// # Parameter Requirements
/// - Parameter-free models: A2, A3, A4, R2, R3, D1, D2, D3, D4, P2_3, P2, P3, F1, F2, F3
/// - Parameterized models:
///   - SB: 3 parameters (m, n, p)
///   - JMA: 1 parameter (m)
///   - Ac: 1 parameter (m)
///   - Dec: 1 parameter (m)
///   - PTe: 2 parameters (m, n)
///   - SBtp: 3 parameters (m, n, c)
#[derive(Debug, Clone, PartialEq, EnumIter, Eq, Hash, Copy)]
pub enum KineticModelNames {
    A2,
    A3,
    A4,
    R2,
    R3,
    D1,
    D2,
    D3,
    D4,
    P2_3,
    P2,
    P3,
    F1,
    F2,
    F3,
    SB,
    JMA,
    Ac,
    Dec,
    PTe,
    SBtp,
}

impl KineticModelNames {
    pub fn map_of_names_and_formulas(&self) -> String {
        match self {
            KineticModelNames::A2 => "2*(1-a)*(-ln(1-a))^(1/2)".to_string(),
            KineticModelNames::A3 => "3*(1-a)*(-ln(1-a))^(2/3)".to_string(),
            KineticModelNames::A4 => "4*(1-a)*(-ln(1-a))^(3/4)".to_string(),
            KineticModelNames::R2 => "(1-a)^(1/2)".to_string(),
            KineticModelNames::R3 => "(1-a)^(2/3)".to_string(),
            KineticModelNames::D1 => "a^(-1)".to_string(),
            KineticModelNames::D2 => "(-ln(1-a))^(-1)".to_string(),
            KineticModelNames::D3 => "1.5*(1-a)^(2/3)*(1-(1 -a)^(1/3))^(-1)".to_string(),
            KineticModelNames::D4 => "(3/2)*((1-a)^(-1/3)-1)^-1".to_string(),
            KineticModelNames::P2_3 => "a^(-0.5)".to_string(),
            KineticModelNames::P2 => "a^(1/2)".to_string(),
            KineticModelNames::P3 => "a^(0.75)".to_string(),
            KineticModelNames::F1 => "(1-a)".to_string(),
            KineticModelNames::F2 => "(1-a)^(2)".to_string(),
            KineticModelNames::F3 => "(1-a)^(3)".to_string(),
            KineticModelNames::SB => " (1-a)^m*a^n*(-ln(1-a))^p".to_string(),
            KineticModelNames::JMA => "m*(1-a)*(-ln(1-a))^(1-1/m)".to_string(),
            KineticModelNames::Ac => "m*a^(1-1/m)".to_string(),
            KineticModelNames::Dec => "(1-a)^m".to_string(),
            KineticModelNames::PTe => "(1-a)^m * a^n".to_string(),
            KineticModelNames::SBtp => "c*(1-a)^m*a^n".to_string(),
        }
    }
    pub fn required_params(&self) -> Vec<String> {
        match self {
            KineticModelNames::A2
            | KineticModelNames::A3
            | KineticModelNames::A4
            | KineticModelNames::R2
            | KineticModelNames::R3
            | KineticModelNames::D1
            | KineticModelNames::D2
            | KineticModelNames::D3
            | KineticModelNames::D4
            | KineticModelNames::P2_3
            | KineticModelNames::P2
            | KineticModelNames::P3
            | KineticModelNames::F1
            | KineticModelNames::F2
            | KineticModelNames::F3 => vec![],
            KineticModelNames::SB => vec!["m".to_string(), "n".to_string(), "p".to_string()],
            KineticModelNames::JMA => vec!["m".to_string()],
            KineticModelNames::Ac => vec!["m".to_string()],
            KineticModelNames::Dec => vec!["m".to_string()],
            KineticModelNames::PTe => vec!["m".to_string(), "n".to_string()],
            KineticModelNames::SBtp => vec!["m".to_string(), "n".to_string(), "c".to_string()],
        }
    }

    pub fn from_names_to_expr(&self, vec_of_params: Vec<f64>) -> KineticModel {
        match self {
            KineticModelNames::A2 => KineticModel::A2,
            KineticModelNames::A3 => KineticModel::A3,
            KineticModelNames::A4 => KineticModel::A4,
            KineticModelNames::R2 => KineticModel::R2,
            KineticModelNames::R3 => KineticModel::R3,
            KineticModelNames::D1 => KineticModel::D1,
            KineticModelNames::D2 => KineticModel::D2,
            KineticModelNames::D3 => KineticModel::D3,
            KineticModelNames::D4 => KineticModel::D4,
            KineticModelNames::P2_3 => KineticModel::P2_3,
            KineticModelNames::P2 => KineticModel::P2,
            KineticModelNames::P3 => KineticModel::P3,
            KineticModelNames::F1 => KineticModel::F1,
            KineticModelNames::F2 => KineticModel::F2,
            KineticModelNames::F3 => KineticModel::F3,
            KineticModelNames::SB => {
                let (m, n, p) = (vec_of_params[0], vec_of_params[1], vec_of_params[2]);
                KineticModel::SB(m, n, p)
            }
            KineticModelNames::JMA => {
                let m = vec_of_params[0];
                KineticModel::JMA(m)
            }
            KineticModelNames::Ac => {
                let m = vec_of_params[0];
                KineticModel::Ac(m)
            }
            KineticModelNames::Dec => {
                let m = vec_of_params[0];
                KineticModel::Dec(m)
            }
            KineticModelNames::PTe => {
                let (m, n) = (vec_of_params[0], vec_of_params[1]);
                KineticModel::PTe(m, n)
            }
            KineticModelNames::SBtp => {
                let (m, n, c) = (vec_of_params[0], vec_of_params[1], vec_of_params[2]);
                KineticModel::SBtp(m, n, c)
            }
        }
    }
    pub fn get_fn(&self) -> Box<dyn Fn(f64, Vec<f64>) -> f64 + '_> {
        let closure_for_expr = |vec_of_params: Vec<f64>| self.from_names_to_expr(vec_of_params);
        let func = move |a: f64, vec_of_params: Vec<f64>| {
            let this_is_expression = closure_for_expr(vec_of_params);
            let res = this_is_expression.get_fn();
            res(a)
        };
        Box::new(func)
    }
    pub fn get_g(&self) -> Result<Box<dyn Fn(f64) -> f64>, String> {
        let closure: Box<dyn Fn(f64) -> f64> = match self {
            KineticModelNames::A2 => Box::new(|a: f64| (-(1.0 - a).ln()).powf(0.5_f64)),
            KineticModelNames::A3 => Box::new(|a: f64| (-(1.0 - a).ln()).powf(1.0 / 3.0)),
            KineticModelNames::A4 => Box::new(|a: f64| (-(1.0 - a).ln()).powf(1.0 / 4.0)),
            KineticModelNames::R2 => Box::new(|a: f64| 2.0 - 2.0 * (1.0 - a).powf(0.5)),
            KineticModelNames::R3 => Box::new(|a: f64| 3.0 * (1.0 - (1.0 - a).powf(1.0 / 3.0))),
            KineticModelNames::D1 => Box::new(|a: f64| 0.5 * a.powf(2.0)),
            KineticModelNames::D2 => Box::new(|a: f64| a + (1.0 - a) * (1.0 - a).ln()),
            KineticModelNames::D3 => Box::new(|a: f64| (1.0 - (1.0 - a).powf(1.0 / 3.0)).powi(2)),
            KineticModelNames::D4 => {
                Box::new(|a: f64| -(2.0 / 3.0) * a - (1.0 - a).powf(2.0 / 3.0))
            }
            KineticModelNames::P2_3 => Box::new(|a: f64| (2.0 / 3.0) * a.powf(1.5)),
            KineticModelNames::P2 => Box::new(|a: f64| 2.0 * a.sqrt()),
            KineticModelNames::P3 => Box::new(|a: f64| 4.0 * a.powf(0.25)),
            KineticModelNames::F1 => Box::new(|a: f64| -(1.0 - a).ln()),
            KineticModelNames::F2 => Box::new(|a: f64| (1.0 / (1.0 - a)) - 1.0),

            KineticModelNames::F3 => {
                Box::new(|a: f64| (1.0 / (2.0 * (1.0 - a).powf(2.0))) - 1.0 / 2.0)
            }

            _ => return Err("no analutical integral".to_string()),
        };
        Ok(closure)
    }
    pub fn pretty_print() {
        #[derive(Tabled)]
        struct ModelRow {
            code: String,
            expression: String,
        }

        let data: Vec<ModelRow> = KineticModelNames::iter()
            .map(|model| ModelRow {
                code: format!("{:?}", model),
                expression: model.map_of_names_and_formulas(),
            })
            .collect();

        let mut binding = Table::new(data);

        let table = binding.with(tabled::settings::Style::rounded());
        println!("{}", table);
    }
}

/// Concrete kinetic model instances with their parameters.
///
/// This enum represents instantiated kinetic models with specific parameter values.
/// Each variant corresponds to a `KineticModelNames` variant but includes the actual
/// parameter values needed for computation.
///
/// # Examples
/// ```rust, ignore
/// use your_crate::KineticModel;
///
/// let model = KineticModel::JMA(2.0); // Johnson-Mehl-Avrami with m=2.0
/// let expr = model.get_expr(); // Get symbolic expression
/// ```
#[derive(Debug, Clone, PartialEq)]
pub enum KineticModel {
    A2,
    A3,
    A4,
    R2,
    R3,
    D1,
    D2,
    D3,
    D4,
    P2_3,
    P2,
    P3,
    F1,
    F2,
    F3,
    SB(f64, f64, f64),
    JMA(f64),
    Ac(f64),
    Dec(f64),
    PTe(f64, f64),
    SBtp(f64, f64, f64),
}

impl KineticModel {
    pub fn required_params(&self) -> Vec<String> {
        match self {
            KineticModel::A2
            | KineticModel::A3
            | KineticModel::A4
            | KineticModel::R2
            | KineticModel::R3
            | KineticModel::D1
            | KineticModel::D2
            | KineticModel::D3
            | KineticModel::D4
            | KineticModel::P2_3
            | KineticModel::P2
            | KineticModel::P3
            | KineticModel::F1
            | KineticModel::F2
            | KineticModel::F3 => vec![],
            KineticModel::SB(_, _, _) => vec!["m".to_string(), "n".to_string(), "p".to_string()],
            KineticModel::JMA(_) => vec!["m".to_string()],
            KineticModel::Ac(_) => vec!["m".to_string()],
            KineticModel::Dec(_) => vec!["m".to_string()],
            KineticModel::PTe(_, _) => vec!["m".to_string(), "n".to_string()],
            KineticModel::SBtp(_, _, _) => vec!["m".to_string(), "n".to_string(), "c".to_string()],
        }
    }
    pub fn from_name_and_params(name: KineticModelNames, params: Vec<f64>) -> Result<Self, String> {
        match name {
            KineticModelNames::A2 => Ok(KineticModel::A2),
            KineticModelNames::A3 => Ok(KineticModel::A3),
            KineticModelNames::A4 => Ok(KineticModel::A4),
            KineticModelNames::R2 => Ok(KineticModel::R2),
            KineticModelNames::R3 => Ok(KineticModel::R3),
            KineticModelNames::D1 => Ok(KineticModel::D1),
            KineticModelNames::D2 => Ok(KineticModel::D2),
            KineticModelNames::D3 => Ok(KineticModel::D3),
            KineticModelNames::D4 => Ok(KineticModel::D4),
            KineticModelNames::P2_3 => Ok(KineticModel::P2_3),
            KineticModelNames::P2 => Ok(KineticModel::P2),
            KineticModelNames::P3 => Ok(KineticModel::P3),
            KineticModelNames::F1 => Ok(KineticModel::F1),
            KineticModelNames::F2 => Ok(KineticModel::F2),
            KineticModelNames::F3 => Ok(KineticModel::F3),
            KineticModelNames::SB => {
                if params.len() != 3 {
                    return Err("SB model requires 3 parameters (m, n, p)".to_string());
                }
                Ok(KineticModel::SB(params[0], params[1], params[2]))
            }
            KineticModelNames::JMA => {
                if params.len() != 1 {
                    return Err("JMA model requires 1 parameter (m)".to_string());
                }
                Ok(KineticModel::JMA(params[0]))
            }
            KineticModelNames::Ac => {
                if params.len() != 1 {
                    return Err("Ac model requires 1 parameter (m)".to_string());
                }
                Ok(KineticModel::Ac(params[0]))
            }
            KineticModelNames::Dec => {
                if params.len() != 1 {
                    return Err("Dec model requires 1 parameter (m)".to_string());
                }
                Ok(KineticModel::Dec(params[0]))
            }
            KineticModelNames::PTe => {
                if params.len() != 2 {
                    return Err("PTe model requires 2 parameters (m, n)".to_string());
                }
                Ok(KineticModel::PTe(params[0], params[1]))
            }
            KineticModelNames::SBtp => {
                if params.len() != 3 {
                    return Err("SBtp model requires 3 parameters (m, n, c)".to_string());
                }
                Ok(KineticModel::SBtp(params[0], params[1], params[2]))
            }
        }
    }

    pub fn get_expr(&self) -> Expr {
        match self {
            KineticModel::A2 => A2(),
            KineticModel::A3 => A3(),
            KineticModel::A4 => A4(),
            KineticModel::R2 => R2(),
            KineticModel::R3 => R3(),
            KineticModel::D1 => D1(),
            KineticModel::D2 => D2(),
            KineticModel::D3 => D3(),
            KineticModel::D4 => D4(),
            KineticModel::P2_3 => P2_3(),
            KineticModel::P2 => P2(),
            KineticModel::P3 => P3(),
            KineticModel::F1 => F1(),
            KineticModel::F2 => F2(),
            KineticModel::F3 => F3(),
            KineticModel::SB(m, n, p) => Sestak_Berggen(*m, *n, *p),
            KineticModel::JMA(m) => Johnson_Mehl_Avrami(*m),
            KineticModel::Ac(m) => Accelerating_model(*m),
            KineticModel::Dec(m) => Decelerating_model(*m),
            KineticModel::PTe(m, n) => Prout_Tompkins_ext(*m, *n),
            KineticModel::SBtp(m, n, c) => Sestak_Berggren_trunc_with_param(*m, *n, *c),
        }
    }

    fn get_jac(&self) -> Expr {
        match self {
            KineticModel::SB(m, n, p) => SB_j(*m, *n, *p),
            KineticModel::JMA(m) => JMA_jac(*m),
            KineticModel::Ac(m) => Ac_jac(*m),
            KineticModel::Dec(m) => Dec_jac(*m),
            KineticModel::PTe(m, n) => PTe_jac(*m, *n),
            KineticModel::SBtp(m, n, c) => SBtp_jac(*m, *n, *c),
            _ => self.get_expr().diff("a"),
        }
    }
    pub fn get_fn(&self) -> Box<dyn Fn(f64) -> f64> {
        let expr = self.get_expr();
        let closure = expr.lambdify1D();
        closure
    }
}

/// Creates a kinetic model instance and its jacobian from a model name and parameters.
///
/// This convenience function combines model instantiation and jacobian computation
/// in a single call with proper error handling.
///
/// # Arguments
/// * `name` - The kinetic model type to create
/// * `params` - Parameters required by the model (empty for parameter-free models)
///
/// # Returns
/// * `Ok((KineticModel, Expr))` - The model instance and its jacobian expression
/// * `Err(String)` - Error message if parameters are invalid
///
/// # Examples
/// ```rust, ignore
/// let (model, jac) = create_kinetic_model(KineticModelNames::JMA, vec![2.0])?;
/// let (model, jac) = create_kinetic_model(KineticModelNames::A2, vec![])?;
/// ```
pub fn create_kinetic_model(
    name: KineticModelNames,
    params: Vec<f64>,
) -> Result<(KineticModel, Expr), String> {
    match KineticModel::from_name_and_params(name, params) {
        Ok(kinmodel) => {
            let jac = kinmodel.get_jac();
            Ok((kinmodel, jac))
        }
        Err(e) => Err(e),
    }
}
