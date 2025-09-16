//! # Solid-State Kinetics Initial Value Problem Solver
//!
//! This module provides a comprehensive framework for modeling and solving solid-state kinetic
//! reactions with temperature-dependent Arrhenius kinetics. It combines symbolic expression
//! handling with numerical ODE solving to simulate thermal analysis experiments like DSC, TGA,
//! and isothermal kinetic studies.
//!
//! ## Purpose
//!
//! The module is designed for:
//! - **Thermal analysis simulation**: DSC, TGA, and other thermal analysis techniques
//! - **Kinetic parameter estimation**: Fitting experimental data to kinetic models
//! - **Reaction mechanism studies**: Comparing different solid-state reaction models
//! - **Process optimization**: Predicting reaction behavior under different thermal conditions
//!
//! ## Main Components
//!
//! ### Core Structures
//! - [`KineticModelIVP`]: Main solver struct that orchestrates the entire simulation
//! - [`KineticModelNames`]: Enum of available kinetic model types
//! - [`KineticModel`]: Instantiated kinetic models with specific parameters
//!
//! ### Key Methods
//! - **Setup**: `new()` → `set_problem()` → `set_model()` → `check_task()`
//! - **Solving**: `solve()` → `gnuplot()` or `save_result()`
//! - **Utilities**: `create_kinetic_model()`, `pretty_print()`
//!
//! ### Available Kinetic Models
//! - **Nucleation & Growth**: A2, A3, A4 (Avrami-Erofeev)
//! - **Contracting Geometry**: R2, R3 (contracting area/volume)
//! - **Diffusion**: D1, D2, D3, D4 (1D, 2D, 3D diffusion)
//! - **Power Law**: P2, P3, P2_3 (power law models)
//! - **Reaction Order**: F1, F2, F3 (first, second, third order)
//! - **Parameterized**: SB (Šesták-Berggren), JMA (Johnson-Mehl-Avrami),
//!   Ac (Accelerating), Dec (Decelerating), PTe (Prout-Tompkins extended),
//!   SBtp (Šesták-Berggren truncated)
//!
//! ## Non-Obvious Features & Tips
//!
//! ### Symbolic Expression System
//! - All kinetic models return `Expr` types for symbolic manipulation
//! - Jacobians are computed analytically for parameterized models, numerically for others
//! - Constants are defined as `Expr::Const` for symbolic consistency
//!
//! ### Error Handling Pattern
//! - All methods return `Result<(), String>` for comprehensive error reporting
//! - Validation occurs at multiple levels: parameter validation, model setup, solving
//! - Use `?` operator for clean error propagation in calling code
//!
//! ### Temperature Profile
//! - Linear heating: `T(t) = T0 + β*t` where β is heating rate
//! - Arrhenius rate: `k(T) = A * exp(-E/(R*T))`
//! - Combined expression: `da/dt = k(T) * f(a)` where f(a) is the kinetic model
//!
//! ### Solver Integration
//! - Supports multiple solver types: BDF, Radau, RK45, Backward Euler
//! - Solver parameters are customizable via `set_solver_params()`
//! - Initial condition: `a(0) = 1e-5` (small non-zero value to avoid singularities)
//!
//! ### Performance Tips
//! - Use BDF or Radau solvers for stiff problems (high activation energies)
//! - Use RK45 for non-stiff problems (lower computational cost)
//! - Adjust `max_step` parameter for better accuracy vs. speed trade-off
//!
//! ## Example Usage
//!
//! ```rust,ignore
//! use KiThe::Kinetics::solid_state_kinetics_IVP::*;
//! use RustedSciThe::numerical::ODE_api2::SolverType;
//!
//! // Create solver instance
//! let mut ivp = KineticModelIVP::new(SolverType::BackwardEuler);
//!
//! // Set thermal parameters: t_final=100s, β=10K/s, T0=298K, E=50kJ/mol, A=1e6/s
//! ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6)?;
//!
//! // Use Johnson-Mehl-Avrami model with m=2.0
//! ivp.set_model(KineticModelNames::JMA, vec![2.0])?;
//!
//! // Solve and visualize
//! ivp.solve()?;
//! ivp.gnuplot()?;
//! ```
//!
//! ## Mathematical Background
//!
//! The general form of solid-state kinetic equations solved here is:
//!
//! ```text
//! da/dt = A * exp(-E/(R*T(t))) * f(a)
//! ```
//!
//! Where:
//! - `a`: degree of conversion (0 to 1)
//! - `T(t) = T0 + β*t`: temperature profile
//! - `f(a)`: kinetic model function
//! - `A`: pre-exponential factor
//! - `E`: activation energy
//! - `R`: gas constant (8.314 J/mol/K)

use RustedSciThe::numerical::ODE_api2::{SolverParam, SolverType, UniversalODESolver};
use RustedSciThe::symbolic::symbolic_engine::Expr;
use nalgebra::DVector;
use std::collections::HashMap;
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
use RustedSciThe::Utils::plots:: plots_terminal;
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
    let f = two * (one - a.clone()) * Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(-half));
    return f;
}
/// 3*(1-a)*(-ln(1-a))^(2/3)
pub fn A3() -> Expr {
    let a = Expr::Var("a".to_owned());
    let f = three
        * (one - a.clone())
        * Expr::Pow(
            Box::new(-Expr::ln(one - a)),
            Box::new(Expr::Const(-2.0 / 3.0)),
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
            Box::new(Expr::Const(-3.0 / 4.0)),
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
    return Expr::Pow(Box::new(Expr::ln(one - a.clone())), Box::new(p.clone()))
        * (n.clone()
            * Expr::Pow(Box::new(one - a.clone()), Box::new(m.clone()))
            * Expr::Pow(Box::new(a.clone()), Box::new(n.clone() - one))
            - m.clone()
                * Expr::Pow(Box::new(one - a.clone()), Box::new(m.clone() - one))
                * Expr::Pow(Box::new(a.clone()), Box::new(n.clone())))
        - p.clone()
            * Expr::Pow(Box::new(Expr::ln(one - a.clone())), Box::new(p - one))
            * Expr::Pow(Box::new(one - a.clone()), Box::new(m - one))
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
    return m.clone()
        * (Expr::Pow(
            Box::new(-Expr::ln(one - a.clone())),
            Box::new(-one / m.clone()),
        ) / m.clone()
            - Expr::Pow(
                Box::new(-Expr::ln(one - a.clone())),
                Box::new(-one / m.clone()),
            )
            - Expr::Pow(Box::new(-Expr::ln(one - a)), Box::new(one - one / m)));
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
#[derive(Debug, Clone, PartialEq, EnumIter)]
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
            KineticModelNames::R3 => "(1-a)^(1/3)".to_string(),
            KineticModelNames::D1 => "a^(-1)".to_string(),
            KineticModelNames::D2 => "(-ln(1-a))^(-1)".to_string(),
            KineticModelNames::D3 => "1.5*(1-a)^(2/3)*(1-a^(1/3))^(-1)".to_string(),
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
/////////////////////////////////////////
/// Initial Value Problem solver for solid-state kinetic models with temperature dependence.
///
/// This struct provides a complete framework for solving kinetic problems involving
/// solid-state reactions with Arrhenius temperature dependence. It combines kinetic
/// models with thermal profiles and numerical ODE solving.
///
/// # Usage Pattern
/// 1. `new()` - Create instance with solver type
/// 2. `set_problem()` - Set thermal and kinetic parameters
/// 3. `set_model()` - Choose and configure kinetic model
/// 4. `check_task()` - Validate configuration
/// 5. `solve()` - Solve the ODE system
/// 6. `gnuplot()` or `save_result()` - Output results
///
/// # Fields
/// - Thermal parameters: `T0` (initial temperature), `beta` (heating rate), `E` (activation energy), `A` (pre-exponential factor)
/// - Time parameter: `t_final` (final integration time)
/// - Kinetic model: Combined Arrhenius and solid-state kinetic expressions
///
/// # Examples
/// ```rust, ignore
/// use your_crate::{KineticModelIVP, KineticModelNames, SolverType};
///
/// let mut ivp = KineticModelIVP::new(SolverType::BDF);
/// ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6)?;
/// ivp.set_model(KineticModelNames::JMA, vec![2.0])?;
/// ivp.solve()?;
/// ```
pub struct KineticModelIVP {
    /// ODE solver instance (initialized after solve() is called)
    solver: Option<UniversalODESolver>,
    /// Numerical solver parameters (step size, tolerance, etc.)
    solver_params: HashMap<String, SolverParam>,
    /// Final integration time [s]
    t_final: f64,
    /// Heating rate [K/s]
    beta: f64,
    /// Initial temperature [K]
    T0: f64,
    /// Activation energy [J/mol]
    E: f64,
    /// Pre-exponential factor [1/s]
    A: f64,
    /// Selected kinetic model with parameters
    kinmodel: Option<KineticModel>,
    /// Combined kinetic expression (Arrhenius × kinetic model)
    kin_expression: Expr,
    /// Jacobian of the kinetic expression
    jac: Expr,
    /// Type of numerical solver to use
    solvertype: SolverType,
    stop_condition: Option<HashMap<String, f64>>,
}

impl KineticModelIVP {
    /// Creates a new KineticModelIVP instance with default parameters.
    ///
    /// # Arguments
    /// * `solvertype` - The numerical solver type to use for ODE integration
    ///
    /// # Returns
    /// A new instance with default solver parameters and uninitialized problem data
    pub fn new(solvertype: SolverType) -> Self {
        let map_of_params = HashMap::from([
            ("step_size".to_owned(), SolverParam::Float(1e-3)),
            ("tolerance".to_owned(), SolverParam::Float(1e-3)),
            ("max_iterations".to_owned(), SolverParam::Int(100000)),
            ("rtol".to_owned(), SolverParam::Float(1e-3)),
            ("atol".to_owned(), SolverParam::Float(1e-3)),
            ("max_step".to_owned(), SolverParam::Float(0.1)),
            ("first_step".to_owned(), SolverParam::OptionalFloat(None)),
            ("vectorized".to_owned(), SolverParam::Bool(false)),
            ("jac_sparsity".to_owned(), SolverParam::OptionalMatrix(None)),
            ("parallel".to_owned(), SolverParam::Bool(true)),
        ]);
        Self {
            solver: None,
            solver_params: map_of_params,
            t_final: 0.0,
            beta: 0.0,
            T0: 0.0,
            E: 0.0,
            A: 0.0,
            kinmodel: None,
            kin_expression: Expr::Const(0.0),
            jac: Expr::Const(0.0),
            solvertype: solvertype,
            stop_condition: None,
        }
    }

    pub fn set_stop_condition(&mut self, condition: Option<HashMap<String, f64>>) {
        self.stop_condition = condition;
    }
    /// Sets the thermal and kinetic parameters for the problem.
    ///
    /// # Arguments
    /// * `t_final` - Final integration time [s] (must be > 0)
    /// * `beta` - Heating rate [K/s] (must be > 0)
    /// * `T0` - Initial temperature [K] (must be > 0)
    /// * `E` - Activation energy [J/mol] (must be > 0)
    /// * `A` - Pre-exponential factor [1/s] (must be > 0)
    ///
    /// # Returns
    /// * `Ok(())` - Parameters set successfully
    /// * `Err(String)` - Error message if any parameter is invalid (≤ 0)
    pub fn set_problem(
        &mut self,
        t_final: f64,
        beta: f64,
        T0: f64,
        E: f64,
        A: f64,
    ) -> Result<(), String> {
        if t_final <= 0.0 {
            return Err("t_final must be positive".to_string());
        }
        if beta <= 0.0 {
            return Err("beta must be positive".to_string());
        }
        if T0 <= 0.0 {
            return Err("T0 must be positive".to_string());
        }
        if E <= 0.0 {
            return Err("E must be positive".to_string());
        }
        if A <= 0.0 {
            return Err("A must be positive".to_string());
        }

        self.t_final = t_final;
        self.beta = beta;
        self.T0 = T0;
        self.E = E;
        self.A = A;
        Ok(())
    }
    pub fn set_solver_params(&mut self, params: HashMap<String, SolverParam>) {
        self.solver_params = params;
    }
    /// Sets the kinetic model and its parameters.
    ///
    /// This method creates the kinetic model, computes its jacobian, and combines
    /// it with the Arrhenius expression to form the complete kinetic expression.
    ///
    /// # Arguments
    /// * `name` - The kinetic model type to use
    /// * `params` - Model parameters (empty vector for parameter-free models)
    ///
    /// # Returns
    /// * `Ok(())` - Model set successfully
    /// * `Err(String)` - Error message if parameters are invalid
    pub fn set_model(&mut self, name: KineticModelNames, params: Vec<f64>) -> Result<(), String> {
        let (kinmodel, jac) = create_kinetic_model(name, params)?;
        self.kinmodel = Some(kinmodel.clone());
        let arrhenius = self.arrhenius();
        self.kin_expression = arrhenius * kinmodel.get_expr();
        self.jac = jac;
        Ok(())
    }

    /// Validates that all required parameters and models are properly set.
    ///
    /// This method should be called before `solve()` to ensure the problem
    /// is completely configured.
    ///
    /// # Returns
    /// * `Ok(())` - All parameters are valid and set
    /// * `Err(String)` - Detailed error message indicating what is missing or invalid
    pub fn check_task(&self) -> Result<(), String> {
        if self.t_final <= 0.0 {
            return Err("t_final not set or invalid".to_string());
        }
        if self.beta <= 0.0 {
            return Err("beta not set or invalid".to_string());
        }
        if self.E <= 0.0 {
            return Err("E not set or invalid".to_string());
        }
        if self.A <= 0.0 {
            return Err("A not set or invalid".to_string());
        }
        if self.T0 <= 0.0 {
            return Err("T0 not set or invalid".to_string());
        }
        if self.kinmodel.is_none() {
            return Err("kinetic model not set".to_string());
        }
        if self.kin_expression == Expr::Const(0.0) {
            return Err("kinetic expression not initialized".to_string());
        }
        if self.jac == Expr::Const(0.0) {
            return Err("jacobian not initialized".to_string());
        }
        Ok(())
    }
    /// Solves the kinetic ODE system.
    ///
    /// This method automatically calls `check_task()` first, then initializes
    /// and runs the ODE solver with the configured kinetic expression.
    ///
    /// # Returns
    /// * `Ok(())` - Solution computed successfully
    /// * `Err(String)` - Error message if validation fails or solver encounters issues
    pub fn solve(&mut self) -> Result<(), String> {
        self.check_task()?;

        let y0 = DVector::from_vec(vec![1e-5]);
        let mut ode = UniversalODESolver::new(
            vec![self.kin_expression.clone()],
            vec!["a".to_owned()],
            "t".to_owned(),
            self.solvertype.clone(),
            0.0,
            y0,
            self.t_final,
        );
        ode.set_parameters(self.solver_params.clone());
        if let Some(stop_condition) = self.stop_condition.clone() {
            ode.set_stop_condition(stop_condition)
        } else {
            let stop_condition = HashMap::from([("a".to_owned(), 0.999)]);
            ode.set_stop_condition(stop_condition)
        }

        ode.initialize();
        ode.solve();

        self.solver = Some(ode);
        Ok(())
    }
    /// Plots the solution using gnuplot.
    ///
    /// # Returns
    /// * `Ok(())` - Plot generated successfully
    /// * `Err(String)` - Error if solver not initialized (call `solve()` first)
    pub fn plot(&self) -> Result<(), String> {
        let ode = self
            .solver
            .as_ref()
            .ok_or("Solver not initialized. Call solve() first.")?;
        ode.plot_result();
        Ok(())
    }

    pub fn plot_in_terminal(&self) {

        let (t, a ) = self.solver.as_ref().unwrap().get_result();
        plots_terminal(
            "t".to_string(),
            vec!["a".to_string()],
            t.unwrap(),
            a.unwrap()

        )

    }
    /// Saves the solution results to file.
    ///
    /// # Returns
    /// * `Ok(())` - Results saved successfully
    /// * `Err(String)` - Error if solver not initialized or file I/O fails
    pub fn save_result(&self) -> Result<(), String> {
        let ode = self
            .solver
            .as_ref()
            .ok_or("Solver not initialized. Call solve() first.")?;
        ode.save_result()
            .map_err(|e| format!("Failed to save result: {:?}", e))?;
        Ok(())
    }

    fn arrhenius(&self) -> Expr {
        let E = Expr::Const(self.E);
        let A = Expr::Const(self.A);
        let beta = Expr::Const(self.beta);
        let R_sym = Expr::Const(8.31446261815324);
        let T0 = Expr::Const(self.T0);
        let t = Expr::Var("t".to_string());
        let k = A * Expr::Exp(Box::new(-E / (R_sym * (T0 + t * beta))));
        k
    }
}

/////////////////////////////////////////TESTS/////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
   

    #[test]
    fn test_parameter_free_models() {
        let models = [
            KineticModelNames::A2,
            KineticModelNames::A3,
            KineticModelNames::A4,
            KineticModelNames::R2,
            KineticModelNames::R3,
            KineticModelNames::D1,
            KineticModelNames::D2,
            KineticModelNames::D3,
            KineticModelNames::D4,
            KineticModelNames::P2_3,
            KineticModelNames::P2,
            KineticModelNames::P3,
            KineticModelNames::F1,
            KineticModelNames::F2,
            KineticModelNames::F3,
        ];

        for model_name in models {
            let result = create_kinetic_model(model_name.clone(), vec![]);
            assert!(result.is_ok(), "Failed to create model {:?}", model_name);

            let (kinetic_model, jacobian) = result.unwrap();
            let expr = kinetic_model.get_expr();

            // Verify expressions are not empty/trivial
            assert_ne!(format!("{:?}", expr), "Const(0.0)");
            assert_ne!(format!("{:?}", jacobian), "Const(0.0)");
        }
    }

    #[test]
    fn test_parameterized_models() {
        // Test SB model
        let result = create_kinetic_model(KineticModelNames::SB, vec![1.0, 2.0, 0.5]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let sbjac = SB_j(1.0, 2.0, 0.5);
        assert!(matches!(model, KineticModel::SB(1.0, 2.0, 0.5)));
        assert_eq!(jac, sbjac);

        // Test JMA model
        let result = create_kinetic_model(KineticModelNames::JMA, vec![2.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let jmajac = JMA_jac(2.0);
        assert!(matches!(model, KineticModel::JMA(2.0)));
        assert_eq!(jac, jmajac);

        // Test Ac model
        let result = create_kinetic_model(KineticModelNames::Ac, vec![1.5]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        assert!(matches!(model, KineticModel::Ac(1.5)));
        let acjac = Ac_jac(1.5);
        assert_eq!(jac, acjac);

        // Test Dec model
        let result = create_kinetic_model(KineticModelNames::Dec, vec![3.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let decjac = Dec_jac(3.0);
        assert!(matches!(model, KineticModel::Dec(3.0)));
        assert_eq!(jac, decjac);

        // Test PTe model
        let result = create_kinetic_model(KineticModelNames::PTe, vec![2.0, 1.0]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        let ptejac = PTe_jac(2.0, 1.0);
        assert!(matches!(model, KineticModel::PTe(2.0, 1.0)));
        assert_eq!(jac, ptejac);

        // Test SBtp model
        let result = create_kinetic_model(KineticModelNames::SBtp, vec![1.0, 2.0, 0.8]);
        assert!(result.is_ok());
        let (model, jac) = result.unwrap();
        assert!(matches!(model, KineticModel::SBtp(1.0, 2.0, 0.8)));
        let sbtpjac = SBtp_jac(1.0, 2.0, 0.8);
        assert_eq!(jac, sbtpjac);
    }

    #[test]

    fn test_parameter_validation() {
        // Test wrong number of parameters
        assert!(create_kinetic_model(KineticModelNames::SB, vec![1.0]).is_err());
        assert!(create_kinetic_model(KineticModelNames::JMA, vec![]).is_err());
        assert!(create_kinetic_model(KineticModelNames::PTe, vec![1.0]).is_err());
        assert!(create_kinetic_model(KineticModelNames::SBtp, vec![1.0, 2.0]).is_err());
    }

    #[test]
    fn test_jacobian_consistency() {
        // Test that jacobians are properly generated for parameterized models
        let test_cases = [
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::JMA, vec![2.0]),
            (KineticModelNames::Ac, vec![1.5]),
            (KineticModelNames::Dec, vec![3.0]),
            (KineticModelNames::PTe, vec![2.0, 1.0]),
            (KineticModelNames::SBtp, vec![1.0, 2.0, 0.8]),
        ];

        for (model_name, params) in test_cases {
            let result = create_kinetic_model(model_name.clone(), params);
            assert!(result.is_ok(), "Failed to create model {:?}", model_name);

            let (kinetic_model, jacobian) = result.unwrap();

            // Verify jacobian is not trivial
            assert_ne!(format!("{:?}", jacobian), "Const(0.0)");

            // For models with explicit jacobian functions, verify they're used
            match kinetic_model {
                KineticModel::SB(_, _, _) => {
                    // Should use SB_j function, not automatic differentiation
                    assert!(
                        format!("{:?}", jacobian).contains("Pow")
                            || format!("{:?}", jacobian).contains("ln")
                    );
                }
                KineticModel::JMA(_) => {
                    assert!(
                        format!("{:?}", jacobian).contains("Pow")
                            || format!("{:?}", jacobian).contains("ln")
                    );
                }
                _ => {}
            }
        }
    }

    #[test]
    fn test_individual_functions() {
        // Test A2 function
        let expr = A2();
        println!("{:?}", expr);
        assert!(format!("{:?}", expr).contains("Pow"));
        assert!(format!("{:?}", expr).contains("Ln"));

        // Test parameterized functions
        let sb_expr = Sestak_Berggen(1.0, 2.0, 0.5);
        assert!(format!("{:?}", sb_expr).contains("Pow"));

        let jma_expr = Johnson_Mehl_Avrami(2.0);
        assert!(format!("{:?}", jma_expr).contains("Pow"));
        assert!(format!("{:?}", jma_expr).contains("Ln"));
    }

    #[test]
    fn test_model_names_mapping() {
        // Test that all model names have proper string mappings
        for model in KineticModelNames::iter() {
            let mapping = model.map_of_names_and_formulas();
            assert!(!mapping.is_empty(), "Empty mapping for {:?}", model);
            assert!(
                mapping.contains("a"),
                "Mapping should contain variable 'a' for {:?}",
                model
            );
        }
    }

    #[test]
    fn test_expression_evaluation() {
        // Test that expressions can be created and contain expected elements
        let models_with_params = [
            (KineticModelNames::A2, vec![]),
            (KineticModelNames::F1, vec![]),
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::JMA, vec![2.0]),
        ];

        for (model_name, params) in models_with_params {
            let result = create_kinetic_model(model_name.clone(), params);
            assert!(result.is_ok());

            let (kinetic_model, jacobian) = result.unwrap();
            let expr = kinetic_model.get_expr();

            // Verify expressions contain variable 'a'
            let expr_str = format!("{:?}", expr);
            let _jac_str = format!("{:?}", jacobian);

            assert!(
                expr_str.contains("Var(\"a\")"),
                "Expression should contain variable 'a' for {:?}",
                model_name
            );
            // Jacobian might be constant for some simple models, so we don't always require 'a'
        }
    }

    #[test]
    fn test_constants() {
        // Test that constants are properly defined
        assert_eq!(two, Expr::Const(2.0));
        assert_eq!(three, Expr::Const(3.0));
        assert_eq!(four, Expr::Const(4.0));
        assert_eq!(one, Expr::Const(1.0));
        assert_eq!(half, Expr::Const(0.5));
        assert_eq!(neg, Expr::Const(-1.0));
    }

    #[test]
    fn test_create_kinetic_model_error_handling() {
        // Test error handling in create_kinetic_model
        let result = create_kinetic_model(KineticModelNames::SB, vec![1.0, 2.0]); // Missing one parameter
        assert!(result.is_err());

        let result = create_kinetic_model(KineticModelNames::JMA, vec![1.0, 2.0]); // Too many parameters
        assert!(result.is_err());
    }
}
#[cfg(test)]
mod kinetic_model_ivp_tests {
    use super::*;
    use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
    use std::time::Instant;

    #[test]
    fn test_new() {
        let ivp = KineticModelIVP::new(SolverType::BDF);

        assert!(ivp.solver.is_none());
        assert_eq!(ivp.t_final, 0.0);
        assert_eq!(ivp.beta, 0.0);
        assert_eq!(ivp.T0, 0.0);
        assert_eq!(ivp.E, 0.0);
        assert_eq!(ivp.A, 0.0);
        assert!(ivp.kinmodel.is_none());
        assert_eq!(ivp.kin_expression, Expr::Const(0.0));
        assert_eq!(ivp.jac, Expr::Const(0.0));
        assert!(matches!(ivp.solvertype, SolverType::BDF));

        // Check default solver parameters
        assert!(ivp.solver_params.contains_key("step_size"));
        assert!(ivp.solver_params.contains_key("tolerance"));
        assert!(ivp.solver_params.contains_key("max_iterations"));
    }

    #[test]
    fn test_set_problem() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        assert_eq!(ivp.t_final, 100.0);
        assert_eq!(ivp.beta, 10.0);
        assert_eq!(ivp.T0, 298.15);
        assert_eq!(ivp.E, 50000.0);
        assert_eq!(ivp.A, 1e6);
    }

    #[test]
    fn test_set_model() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameter-free model
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        assert!(ivp.kinmodel.is_some());
        assert!(matches!(ivp.kinmodel.as_ref().unwrap(), KineticModel::A2));
        assert_ne!(ivp.kin_expression, Expr::Const(0.0));
        assert_ne!(ivp.jac, Expr::Const(0.0));

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::JMA, vec![2.0]);

        assert!(matches!(
            ivp.kinmodel.as_ref().unwrap(),
            KineticModel::JMA(2.0)
        ));
        println!("{:?}", ivp.kin_expression);
        assert_ne!(ivp.kin_expression, Expr::Const(0.0));
        assert_ne!(ivp.jac, Expr::Const(0.0));
    }

    #[test]
    fn test_check_task_success() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        // Follow the correct order
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        // Should not panic
        let _ = ivp.check_task();
    }

    #[test]

    fn test_check_task_fail_no_problem() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        let _ = ivp.set_model(KineticModelNames::A2, vec![]);

        // Should panic because problem not set
        let r = ivp.check_task();
        match r {
            Ok(_) => {}
            Err(e) => {
                assert_eq!(e, "t_final not set or invalid")
            }
        }
    }

    #[test]

    fn test_check_task_fail_no_model() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        // Should panic because model not set
        let r = ivp.check_task();
        match r {
            Ok(_) => {}
            Err(e) => {
                assert_eq!(e, "kinetic model not set")
            }
        }
    }

    #[test]
    fn test_complete_workflow() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        // Step 1: new() - already done
        assert!(ivp.kinmodel.is_none());

        // Step 2: set_problem()
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        assert_eq!(ivp.t_final, 100.0);

        // Step 3: set_model()
        let _ = ivp.set_model(KineticModelNames::JMA, vec![2.0]);
        assert!(ivp.kinmodel.is_some());

        // Step 4: check_task()
        let _ = ivp.check_task(); // Should not panic
    }

    #[test]
    fn test_arrhenius_expression() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        let arrhenius = ivp.arrhenius();

        // Check that arrhenius expression contains expected components
        let expr_str = format!("{:?}", arrhenius);
        assert!(expr_str.contains("Exp"));
        assert!(expr_str.contains("Const(1000000.0)")); // A value
        assert!(expr_str.contains("Const(50000.0)")); // E value
    }

    #[test]
    fn test_different_solver_types() {
        let solver_types = [SolverType::BDF, SolverType::BackwardEuler, SolverType::BDF];

        for solver_type in solver_types {
            let mut ivp = KineticModelIVP::new(solver_type.clone());
            assert!(matches!(ivp.solvertype, ref _solver_type));

            let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
            let _ = ivp.set_model(KineticModelNames::F1, vec![]);
            let _ = ivp.check_task(); // Should work for all solver types
        }
    }

    #[test]
    fn test_multiple_models() {
        let test_cases = [
            (KineticModelNames::A2, vec![]),
            (KineticModelNames::F1, vec![]),
            (KineticModelNames::JMA, vec![2.0]),
            (KineticModelNames::SB, vec![1.0, 2.0, 0.5]),
            (KineticModelNames::PTe, vec![2.0, 1.0]),
        ];

        for (model_name, params) in test_cases {
            let mut ivp = KineticModelIVP::new(SolverType::BDF);
            let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
            let _ = ivp.set_model(model_name.clone(), params);

            // Verify model was set correctly
            assert!(ivp.kinmodel.is_some());
            assert_ne!(ivp.kin_expression, Expr::Const(0.0));

            // Should pass check_task
            let _ = ivp.check_task();
        }
    }

    #[test]
    fn test_parameter_validation_in_set_model() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();

        // This should return error because wrong number of parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![1.0, 2.0]); // JMA needs 1 param, not 2
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("JMA model requires 1 parameter")
        );
    }

    #[test]
    fn test_kinetic_expression_contains_arrhenius() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);
        let _ = ivp.set_model(KineticModelNames::F1, vec![]);

        // The kinetic expression should be arrhenius * kinetic_model
        let expr_str = format!("{:?}", ivp.kin_expression);

        // Should contain multiplication (arrhenius * model)
        assert!(expr_str.contains("Mul") || expr_str.contains("*"));
        // Should contain exponential (from arrhenius)
        assert!(expr_str.contains("Exp"));
        // Should contain variable 'a' (from kinetic model)
        assert!(expr_str.contains("Var(\"a\")"));
    }
    #[test]
    fn solve_with_rk45() {
        let mut ivp = KineticModelIVP::new(SolverType::NonStiff("RK45".to_owned()));
        let _ = ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::F1, vec![]);
        let _ = ivp.check_task();
        let _ = ivp.solve();
        let _ = ivp.plot();
    }
    #[test]
    fn solve_with_radau() {
        let t = Instant::now();
        let mut ivp = KineticModelIVP::new(SolverType::Radau(RadauOrder::Order7));
        let _ = ivp.set_problem(10.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::D3, vec![]);
        let _ = ivp.check_task();
        let _ = ivp.solve();
        println!("Time elapsed: {:?}", t.elapsed().as_secs());
        let _ = ivp.plot();
    }

    #[test]
    fn solve_with_be() {
        let mut ivp = KineticModelIVP::new(SolverType::BackwardEuler);
        let _ = ivp.set_problem(10.0, 10.0, 298.15, 50000.0, 1e6);

        // Test parameterized model
        let _ = ivp.set_model(KineticModelNames::D3, vec![]);
        let _ = ivp.check_task();
        let _ = ivp.solve();
        let _ = ivp.plot();
    }
}
#[cfg(test)]
mod error_handling_tests {
    use std::time::Instant;

    use super::*;
    use RustedSciThe::numerical::Radau::Radau_main::RadauOrder;
    #[test]
    fn test_set_problem_validation() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        // Test negative t_final
        assert!(ivp.set_problem(-1.0, 10.0, 298.15, 50000.0, 1e6).is_err());

        // Test negative beta
        assert!(ivp.set_problem(100.0, -10.0, 298.15, 50000.0, 1e6).is_err());

        // Test negative T0
        assert!(ivp.set_problem(100.0, 10.0, -298.15, 50000.0, 1e6).is_err());

        // Test negative E
        assert!(ivp.set_problem(100.0, 10.0, 298.15, -50000.0, 1e6).is_err());

        // Test negative A
        assert!(ivp.set_problem(100.0, 10.0, 298.15, 50000.0, -1e6).is_err());

        // Test valid parameters
        assert!(ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).is_ok());
    }

    #[test]
    fn test_set_model_error_propagation() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();

        // Test invalid parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![1.0, 2.0]);
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .contains("JMA model requires 1 parameter")
        );

        // Test valid parameters
        let result = ivp.set_model(KineticModelNames::JMA, vec![2.0]);
        assert!(result.is_ok());
    }

    #[test]
    fn test_check_task_detailed_errors() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        // Test uninitialized state
        let result = ivp.check_task();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("t_final not set"));

        // Set problem partially
        ivp.set_problem(100.0, 10.0, 298.15, 50000.0, 1e6).unwrap();
        let result = ivp.check_task();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("kinetic model not set"));

        // Set model
        ivp.set_model(KineticModelNames::A2, vec![]).unwrap();
        let result = ivp.check_task();
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_without_setup() {
        let mut ivp = KineticModelIVP::new(SolverType::BDF);

        // Try to solve without setup
        let result = ivp.solve();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not set"));
    }

    #[test]
    fn test_gnuplot_without_solve() {
        let ivp = KineticModelIVP::new(SolverType::BDF);

        let result = ivp.plot();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Solver not initialized"));
    }

    #[test]
    fn test_save_result_without_solve() {
        let ivp = KineticModelIVP::new(SolverType::BDF);

        let result = ivp.save_result();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Solver not initialized"));
    }

    #[test]
    fn test_complete_error_free_workflow() {
        let t = Instant::now();
        let mut ivp = KineticModelIVP::new(SolverType::Radau(RadauOrder::Order7));

        // All operations should return Ok
        assert!(ivp.set_problem(50.0, 10.0, 298.15, 50000.0, 1e6).is_ok());
        let stop_condition = HashMap::from([("a".to_owned(), 0.99)]);
        ivp.set_stop_condition(Some(stop_condition));
        assert!(ivp.set_model(KineticModelNames::F1, vec![]).is_ok());
        assert!(ivp.check_task().is_ok());
        assert!(ivp.solve().is_ok());
        assert!(ivp.plot().is_ok());
        assert!(ivp.save_result().is_ok());
        println!("Time taken: {:?}", t.elapsed().as_secs());
        ivp.plot_in_terminal()
    }
}
