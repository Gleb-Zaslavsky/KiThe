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
    TGADomainError, check_requirements,
};
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::ColumnNature;
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum IsoconversionalMethod {
    OFW,
    KAS,
    Starink,

    FriedmanDifferential,
    FriedmanIntegral,

    Vyazovkin,
}

impl Default for IsoconversionalMethod {
    fn default() -> Self {
        IsoconversionalMethod::OFW
    }
}

impl IsoconversionalMethod {
    /// All methods available for selection in the UI.
    pub fn all() -> &'static [IsoconversionalMethod] {
        &[
            IsoconversionalMethod::OFW,
            IsoconversionalMethod::KAS,
            IsoconversionalMethod::Starink,
            IsoconversionalMethod::FriedmanDifferential,
            IsoconversionalMethod::FriedmanIntegral,
            IsoconversionalMethod::Vyazovkin,
        ]
    }

    /// A human-friendly name for UI display.
    pub fn display_name(&self) -> &'static str {
        match self {
            IsoconversionalMethod::OFW => "Ozawa-Flynn-Wall (OFW)",
            IsoconversionalMethod::KAS => "Kissinger-Akahira-Sunose (KAS)",
            IsoconversionalMethod::Starink => "Starink",
            IsoconversionalMethod::FriedmanDifferential => "Friedman (Differential)",
            IsoconversionalMethod::FriedmanIntegral => "Friedman (Integral)",
            IsoconversionalMethod::Vyazovkin => "Vyazovkin",
        }
    }

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

    pub fn requirements(&self) -> KineticRequirements {
        match self {
            IsoconversionalMethod::OFW
            | IsoconversionalMethod::KAS
            | IsoconversionalMethod::Starink
            | IsoconversionalMethod::Vyazovkin => KineticRequirements {
                min_experiments: 3,
                needs_conversion: true,
                needs_conversion_rate: false,
                needs_temperature: false,
                needs_heating_rate: true,
            },
            IsoconversionalMethod::FriedmanDifferential => KineticRequirements {
                min_experiments: 3,
                needs_conversion: true,
                needs_conversion_rate: true,
                needs_temperature: false,
                needs_heating_rate: false,
            },
            IsoconversionalMethod::FriedmanIntegral => KineticRequirements {
                min_experiments: 3,
                needs_conversion: true,
                needs_conversion_rate: false,
                needs_temperature: true,
                needs_heating_rate: false,
            },
        }
    }

    pub fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        match self {
            IsoconversionalMethod::FriedmanIntegral => {
                vec![ColumnNature::Conversion, ColumnNature::Time]
            }
            IsoconversionalMethod::OFW
            | IsoconversionalMethod::KAS
            | IsoconversionalMethod::Starink
            | IsoconversionalMethod::FriedmanDifferential
            | IsoconversionalMethod::Vyazovkin => vec![
                ColumnNature::Conversion,
                ColumnNature::Temperature,
                ColumnNature::Time,
                ColumnNature::ConversionRate,
            ],
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

        let grid = match self.method {
            IsoconversionalMethod::FriedmanIntegral => {
                ConversionGridBuilder::default().build_isothermal(data)?
            }
            _ => ConversionGridBuilder::default()
                .compute_dt(compute_dt)
                .build_nonisothermal(data)?,
        };

        // Use method-specific constants for OFW/KAS/Starink.
        let integral_config = match self.method {
            IsoconversionalMethod::OFW
            | IsoconversionalMethod::KAS
            | IsoconversionalMethod::Starink => {
                let int_method = self.method.into_int_isomethod()?;
                IntegralIsoconversionalSolver::new(int_method).config
            }
            _ => IntegralIsoConfig::default(),
        };

        let solver = IsoconversionalSolver {
            method: self.method.clone(),

            integral_config,

            vyazovkin_config: VyazovkinSolver::default(),
        };

        solver.solve(&grid)
    }

    fn requirements(&self) -> KineticRequirements {
        self.method.requirements()
    }

    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        self.method.required_columns_by_nature()
    }

    fn check_input(&self, data: &KineticDataView) -> Result<(), TGADomainError> {
        let req = self.method.requirements();
        check_requirements(data, &req)
    }
}
