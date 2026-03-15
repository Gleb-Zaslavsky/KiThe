//! KineticMethod (trait)
//!      │
//!      ▼
//!IsoconversionalSolver (dispatcher)
//!      │
//!      ├── IntegralIsoconversionalSolver
//!      │        ├── OFW
//!      │        ├── KAS
//!      │        └── Starink
//!      │
//!      ├── DifferentialIsoconversionalSolver
//!      │        └── Friedman
//!      │
//!      └── VariationalIsoconversionalSolver
//!               └── Vyazovkin
//! A WRAPPER FOR ALL ISOCONVERSIONAL METHODS
use crate::Kinetics::experimental_kinetics::kinetic_methods::Friedman::{
    DifferentialFriedmanSolver, FriedmanIntegralSolver,
};
use crate::Kinetics::experimental_kinetics::kinetic_methods::Vyazovkin::VyazovkinSolver;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::IntegralIsoConfig;
use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::{
    IntIsoconversionalMethod, IntegralIsoconversionalSolver, IsoconversionalResult,
};

use crate::Kinetics::experimental_kinetics::kinetic_methods::{
    ConversionGrid, ConversionGridBuilder, KineticDataView, KineticMethod, KineticRequirements,
    TGADomainError,
};
#[derive(Clone, Debug)]
pub enum IsoconversionalMethod {
    OFW,
    KAS,
    Starink,

    FriedmanDifferential,
    FriedmanIntegral,

    Vyazovkin,
}
impl IsoconversionalMethod {
    pub fn into_int_isomethod(&self) -> Result<IntIsoconversionalMethod, TGADomainError> {
        match self {
            IsoconversionalMethod::KAS => Ok(IntIsoconversionalMethod::KAS),
            IsoconversionalMethod::OFW => Ok(IntIsoconversionalMethod::OFW),
            IsoconversionalMethod::Starink => Ok(IntIsoconversionalMethod::Starink),

            IsoconversionalMethod::FriedmanDifferential
            | IsoconversionalMethod::FriedmanIntegral
            | IsoconversionalMethod::Vyazovkin => Err(TGADomainError::InvalidOperation(
                "kinetic method mismatch".into(),
            )),
        }
    }
}
#[derive(Clone, Debug)]
pub struct IsoconversionalSolver {
    pub method: IsoconversionalMethod,

    pub integral_config: IntegralIsoConfig,

    pub vyazovkin_config: VyazovkinSolver,
}

impl IsoconversionalSolver {
    pub fn solve(&self, grid: &ConversionGrid) -> Result<IsoconversionalResult, TGADomainError> {
        match self.method {
            IsoconversionalMethod::OFW
            | IsoconversionalMethod::KAS
            | IsoconversionalMethod::Starink => IntegralIsoconversionalSolver {
                config: self.integral_config.clone(),

                method: self.method.into_int_isomethod()?,
            }
            .solve(grid),

            IsoconversionalMethod::FriedmanIntegral => FriedmanIntegralSolver.solve(&grid),

            IsoconversionalMethod::FriedmanDifferential => DifferentialFriedmanSolver.solve(&grid),

            IsoconversionalMethod::Vyazovkin => self.vyazovkin_config.solve(grid),
        }
    }
}

pub struct IsoconversionalKineticMethod {
    pub method: IsoconversionalMethod,
}

impl KineticMethod for IsoconversionalKineticMethod {
    type Output = IsoconversionalResult;

    fn name(&self) -> &'static str {
        match self.method {
            IsoconversionalMethod::OFW => "Ozawa-Flynn-Wall",

            IsoconversionalMethod::KAS => "Kissinger-Akahira-Sunose",

            IsoconversionalMethod::Starink => "Starink",

            IsoconversionalMethod::FriedmanIntegral => "FriedmanIntegral",
            IsoconversionalMethod::FriedmanDifferential => "FriedmanDifferential",

            IsoconversionalMethod::Vyazovkin => "Vyazovkin",
        }
    }

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError> {
        let compute_dt = matches!(self.method, IsoconversionalMethod::Vyazovkin);

        let grid = ConversionGridBuilder::default()
            .compute_dt(compute_dt)
            .build(data)?;

        let solver = IsoconversionalSolver {
            method: self.method.clone(),

            integral_config: IntegralIsoConfig::default(),

            vyazovkin_config: VyazovkinSolver::default(),
        };

        solver.solve(&grid)
    }

    fn requirements(&self) -> KineticRequirements {
        KineticRequirements {
            min_experiments: 3,

            needs_conversion: true,

            needs_conversion_rate: matches!(self.method, IsoconversionalMethod::FriedmanIntegral),

            needs_temperature: true,

            needs_heating_rate: true,
        }
    }
}
