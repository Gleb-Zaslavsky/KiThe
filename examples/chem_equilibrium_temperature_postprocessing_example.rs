//! Guide: postprocess a temperature sweep into a smoother report.
//!
//! This example keeps the solver separate from export/preview logic and shows
//! how to resample a solved sweep for plotting or logs.

use KiThe::Thermodynamics::ChemEquilibrium::equilibrium_temperature_postprocessing::{
    TemperatureInterpolationPolicy, TemperatureInterpolationSpace, TemperaturePostprocessingPolicy,
    TemperatureResamplingGrid, postprocess_temperature_series,
};

fn main() {
    let rows = vec![
        (1000.0, vec![0.80, 0.15, 0.05]),
        (1300.0, vec![0.72, 0.20, 0.08]),
        (1600.0, vec![0.60, 0.26, 0.14]),
    ];

    let policy = TemperaturePostprocessingPolicy {
        grid: TemperatureResamplingGrid::Uniform { points: 8 },
        interpolation: TemperatureInterpolationPolicy {
            space: TemperatureInterpolationSpace::Log,
            clamp: true,
        },
    };

    let report = postprocess_temperature_series(
        vec!["CO".to_string(), "CO2".to_string(), "O2".to_string()],
        &rows,
        &policy,
    )
    .expect("temperature postprocessing failed");

    println!("{}", report.render_table());
}
