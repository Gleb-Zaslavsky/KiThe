use crate::Kinetics::experimental_kinetics::testing_mod::{
    AdvancedTGAConfig, ExperimentMode, KineticModel, NoiseConfig, NoiseKind, SimulatedTGA,
    VirtualTGA,
};
use tabled::{Table, Tabled, settings::Style};

#[derive(Clone, Copy, Debug)]
pub(super) enum ThermalScenario {
    Isothermal,
    NonIsothermal,
}

impl ThermalScenario {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Isothermal => "isothermal",
            Self::NonIsothermal => "non-isothermal",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) enum NoiseProfile {
    Gaussian,
    Uniform,
    Drift,
    Impulse,
    Flicker,
}

impl NoiseProfile {
    pub(super) const ALL: [Self; 5] = [
        Self::Gaussian,
        Self::Uniform,
        Self::Drift,
        Self::Impulse,
        Self::Flicker,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Gaussian => "Gaussian(sigma=0.12mg)",
            Self::Uniform => "Uniform(amplitude=0.20mg)",
            Self::Drift => "Drift(slope=1e-4mg/point)",
            Self::Impulse => "Impulse(p=0.04, amplitude=0.8mg)",
            Self::Flicker => "Flicker(amplitude=0.05mg)",
        }
    }

    fn config(self) -> NoiseConfig {
        NoiseConfig {
            kind: match self {
                Self::Gaussian => NoiseKind::Gaussian { sigma: 0.12 },
                Self::Uniform => NoiseKind::Uniform { amplitude: 0.20 },
                Self::Drift => NoiseKind::Drift { slope: 1e-4 },
                Self::Impulse => NoiseKind::Impulse {
                    probability: 0.04,
                    amplitude: 0.8,
                },
                Self::Flicker => NoiseKind::Flicker { amplitude: 0.05 },
            },
        }
    }
}

pub(super) fn config(
    scenario: ThermalScenario,
    noise: NoiseProfile,
    seed: u64,
) -> AdvancedTGAConfig {
    AdvancedTGAConfig {
        n_points: 360,
        dt: 2.0,
        experiment_mode: match scenario {
            ThermalScenario::Isothermal => ExperimentMode::Isothermal {
                temperatures: vec![500.0, 530.0, 560.0],
            },
            ThermalScenario::NonIsothermal => ExperimentMode::NonIsothermal {
                t0: 500.0,
                heating_rates: vec![2.0, 5.0, 10.0],
            },
        },
        kinetic_model: KineticModel::ArrheniusSingle {
            m0: 100.0,
            k0: 2.0e5,
            e: 85_000.0,
            r: 8.314,
        },
        mass_noise: Some(noise.config()),
        temp_noise: None,
        spikes: None,
        seed,
    }
}

pub(super) fn noisy_and_truth(
    scenario: ThermalScenario,
    noise: NoiseProfile,
    seed: u64,
) -> (Vec<SimulatedTGA>, Vec<SimulatedTGA>) {
    let noisy_config = config(scenario, noise, seed);
    let noisy = VirtualTGA::generate_series(&noisy_config).unwrap();
    let mut truth_config = noisy_config;
    truth_config.mass_noise = None;
    let truth = VirtualTGA::generate_series(&truth_config).unwrap();
    assert_eq!(noisy.len(), truth.len());
    (noisy, truth)
}

#[derive(Tabled)]
struct GeneratorRow {
    #[tabled(rename = "Thermal case")]
    scenario: String,
    #[tabled(rename = "Noise profile")]
    noise: String,
    #[tabled(rename = "Seed")]
    seed: u64,
    #[tabled(rename = "Runs")]
    runs: usize,
    #[tabled(rename = "Points/run")]
    points: usize,
    #[tabled(rename = "RMSE noisy vs clean (mg)")]
    noise_rmse: String,
    #[tabled(rename = "Repeatable")]
    repeatable: bool,
}

#[test]
fn generator_stage_emits_deterministic_thermal_and_noise_matrix() {
    let mut rows = Vec::new();
    for scenario in [ThermalScenario::Isothermal, ThermalScenario::NonIsothermal] {
        for noise in NoiseProfile::ALL {
            let seed = 71_000 + noise as u64 * 101 + scenario as u64;
            let (runs, truth) = noisy_and_truth(scenario, noise, seed);
            let repeated = VirtualTGA::generate_series(&config(scenario, noise, seed)).unwrap();
            assert_eq!(runs.len(), 3);
            assert_eq!(runs.len(), truth.len());
            assert_eq!(runs[0].virtual_tga.mass, repeated[0].virtual_tga.mass);

            let run = &runs[0].virtual_tga;
            let clean = &truth[0].virtual_tga;
            let noise_rmse = crate::gui::experimental_kinetics_gui::story_tests::support::rmse(
                &run.mass,
                &clean.mass,
            );
            assert!(noise_rmse.is_finite() && noise_rmse > 0.0);
            rows.push(GeneratorRow {
                scenario: scenario.label().into(),
                noise: noise.label().into(),
                seed,
                runs: runs.len(),
                points: run.time.len(),
                noise_rmse: format!("{noise_rmse:.5}"),
                repeatable: true,
            });
        }
    }
    println!(
        "\n=== STORY MATRIX / 1. GENERATOR ===\n{}",
        Table::new(rows).with(Style::rounded())
    );
}
