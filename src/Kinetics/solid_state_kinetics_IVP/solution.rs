use super::models::KineticModel;
use nalgebra::{DMatrix, DVector};

/// Published scalar conversion trajectory and its derived thermochemical axes.
#[derive(Debug, Clone, PartialEq)]
pub struct IVPSolution {
    pub time: Vec<f64>,
    pub temperature: Vec<f64>,
    pub conversion: Vec<f64>,
    pub conversion_rate: Vec<f64>,
}

impl IVPSolution {
    pub(crate) fn from_trajectory(
        time: DVector<f64>,
        states: DMatrix<f64>,
        initial_temperature: f64,
        heating_rate: f64,
        activation_energy: f64,
        pre_exponential_factor: f64,
        model: &KineticModel,
    ) -> Result<Self, String> {
        let conversion: Vec<f64> = if states.ncols() == 1 {
            states.column(0).iter().copied().collect()
        } else if states.nrows() == 1 {
            states.row(0).iter().copied().collect()
        } else {
            return Err("conversion result is not 1D".to_owned());
        };
        let time: Vec<f64> = time.iter().copied().collect();
        let temperature: Vec<f64> = time
            .iter()
            .map(|t| initial_temperature + heating_rate * t)
            .collect();
        let rate_model = model.get_fn();
        let gas_constant = 8.31446261815324_f64;
        let conversion_rate = time
            .iter()
            .zip(conversion.iter())
            .map(|(t, a)| {
                let temperature = initial_temperature + heating_rate * t;
                let arrhenius = pre_exponential_factor
                    * (-activation_energy / (gas_constant * temperature)).exp();
                arrhenius * rate_model(*a)
            })
            .collect();

        Ok(Self {
            time,
            temperature,
            conversion,
            conversion_rate,
        })
    }
}
