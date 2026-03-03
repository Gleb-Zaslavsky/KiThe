use polars::prelude::{DataFrame, LazyFrame, PolarsResult};

use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnMeta, ColumnOrigin, History, TGADataset, TGADomainError, TGASchema, Unit,
};
use polars::prelude::*;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use rand_distr::{Distribution, Normal, Uniform};
use std::collections::HashMap;
use std::io::Write;
use tempfile::NamedTempFile;

impl TGADataset {
    pub fn create_from_synthetic_data(synthetic: &VirtualTGA) -> Result<Self, TGADomainError> {
        let n = synthetic.time.len();
        if synthetic.temperature.len() != n || synthetic.mass.len() != n {
            return Err(TGADomainError::InvalidOperation(
                "Synthetic vectors must have equal lengths".to_string(),
            ));
        }
        let height = synthetic.time.len();
        let time = Column::new("time".into(), synthetic.time.as_slice());
        let mass = Column::new("mass".into(), synthetic.mass.as_slice());
        let temperature = Column::new("temperature".into(), synthetic.temperature.as_slice());
        let frame =
            DataFrame::new(height, vec![time.into(), mass.into(), temperature.into()])?.lazy();

        let mut columns = HashMap::new();
        columns.insert(
            "time".to_string(),
            ColumnMeta {
                name: "time".to_string(),
                unit: Unit::Second,
                origin: ColumnOrigin::Raw,
            },
        );
        columns.insert(
            "mass".to_string(),
            ColumnMeta {
                name: "mass".to_string(),
                unit: Unit::Milligram,
                origin: ColumnOrigin::Raw,
            },
        );
        columns.insert(
            "temperature".to_string(),
            ColumnMeta {
                name: "temperature".to_string(),
                unit: Unit::Celsius,
                origin: ColumnOrigin::Raw,
            },
        );

        Ok(Self {
            frame,
            schema: TGASchema {
                columns,
                time: Some("time".to_string()),
                temperature: Some("temperature".to_string()),
                mass: Some("mass".to_string()),
                alpha: None,
                dm_dt: None,
                eta: None,
                deta_dt: None,
                dalpha_dt: None,
                dT_dt: None,
            },
            oneframeplot: None,
            history_of_operations: History {
                vector_of_operations: Vec::new(),
                columns_has_changed: Vec::new(),
            },
        })
    }
}
//==========================================================================================================
//  COMPLICATED MODEL WITH VARIOUS KINETICS AND NOSES
//==========================================================================================================
#[derive(Clone, Debug)]
pub struct SimulatedTGA {
    pub virtual_tga: VirtualTGA,
    pub meta: ExperimentMeta,
}

#[derive(Clone, Debug)]
pub struct VirtualTGA {
    pub time: Vec<f64>,
    pub temperature: Vec<f64>,
    pub mass: Vec<f64>,
}

// NOISE CONFIGURATION
/// Типы шума для генерации синтетических данных
#[derive(Clone, Debug)]
pub enum NoiseKind {
    Gaussian { sigma: f64 },
    Uniform { amplitude: f64 },
    Drift { slope: f64 },
}

#[derive(Clone, Debug)]
pub struct NoiseConfig {
    pub kind: NoiseKind,
}

//==========================================================================================================
// kinetics config
#[derive(Clone, Debug)]
pub enum KineticModel {
    /// Single-step Arrhenius
    ArrheniusSingle { m0: f64, k0: f64, e: f64, r: f64 },

    /// Two parallel Arrhenius reactions
    ArrheniusTwoComponent {
        m01: f64,
        k01: f64,
        e1: f64,

        m02: f64,
        k02: f64,
        e2: f64,

        r: f64,
    },
}
//=============================================================================================

#[derive(Clone, Debug)]
pub enum ExperimentMode {
    Isothermal {
        temperatures: Vec<f64>, // K
    },
    NonIsothermal {
        t0: f64,                 // start temperature
        heating_rates: Vec<f64>, // K/min
    },
}
//================================================================================
/// CONFIG FOR COMPLEX SYNTHETIC DATA
#[derive(Clone, Debug)]

pub struct AdvancedTGAConfig {
    pub n_points: usize,
    pub dt: f64, // seconds

    pub experiment_mode: ExperimentMode,

    pub kinetic_model: KineticModel,

    pub mass_noise: Option<NoiseConfig>,
    pub temp_noise: Option<NoiseConfig>,
    pub spikes: Option<SpikeModel>,

    pub seed: u64,
}
impl VirtualTGA {
    fn generate_one_isothermal(
        cfg: &AdvancedTGAConfig,
        temperature: f64,
        index: u64,
    ) -> SimulatedTGA {
        let mut rng = StdRng::seed_from_u64(cfg.seed + index);

        let time: Vec<f64> = (0..cfg.n_points).map(|i| i as f64 * cfg.dt).collect();

        let mut temp_vec = vec![temperature; cfg.n_points];

        if let Some(noise) = &cfg.temp_noise {
            Self::apply_noise(&mut temp_vec, noise, &mut rng);
        }

        let mut mass = Self::generate_mass(&time, &temp_vec, &cfg.kinetic_model, cfg.dt);

        if let Some(noise) = &cfg.mass_noise {
            Self::apply_noise(&mut mass, noise, &mut rng);
        }

        if let Some(spike) = &cfg.spikes {
            Self::apply_spikes(&mut mass, spike, &mut rng);
        }

        SimulatedTGA {
            virtual_tga: VirtualTGA {
                time,
                temperature: temp_vec,
                mass,
            },
            meta: ExperimentMeta {
                id: format!("T = {:.1} K", temperature),
                heating_rate: None,
                isothermal_temperature: Some(temperature),
                comment: None,
            },
        }
    }

    fn generate_one_ramp(
        cfg: &AdvancedTGAConfig,
        t0: f64,
        beta_min: f64, // K/min
        index: u64,
    ) -> SimulatedTGA {
        let mut rng = StdRng::seed_from_u64(cfg.seed + index);

        let beta = beta_min / 60.0; // convert to K/s

        let time: Vec<f64> = (0..cfg.n_points).map(|i| i as f64 * cfg.dt).collect();

        let mut temperature: Vec<f64> = time.iter().map(|t| t0 + beta * t).collect();

        if let Some(noise) = &cfg.temp_noise {
            Self::apply_noise(&mut temperature, noise, &mut rng);
        }

        let mut mass = Self::generate_mass(&time, &temperature, &cfg.kinetic_model, cfg.dt);

        if let Some(noise) = &cfg.mass_noise {
            Self::apply_noise(&mut mass, noise, &mut rng);
        }

        if let Some(spike) = &cfg.spikes {
            Self::apply_spikes(&mut mass, spike, &mut rng);
        }

        SimulatedTGA {
            virtual_tga: VirtualTGA {
                time,
                temperature,
                mass,
            },
            meta: ExperimentMeta {
                id: format!("β = {:.1} K/min", beta_min),
                heating_rate: Some(beta_min),
                isothermal_temperature: None,
                comment: None,
            },
        }
    }

    fn generate_mass(time: &[f64], temperature: &[f64], model: &KineticModel, dt: f64) -> Vec<f64> {
        match model {
            KineticModel::ArrheniusSingle { m0, k0, e, r } => {
                let mut m = vec![*m0];
                for i in 0..time.len() - 1 {
                    let t = temperature[i];
                    let k = k0 * (-e / (r * t)).exp();
                    let next = m[i] * (-k * dt).exp();
                    m.push(next);
                }
                m
            }

            KineticModel::ArrheniusTwoComponent {
                m01,
                k01,
                e1,
                m02,
                k02,
                e2,
                r,
            } => {
                let mut m1 = *m01;
                let mut m2 = *m02;

                let mut total = vec![m1 + m2];

                for i in 0..time.len() - 1 {
                    let t = temperature[i];

                    let k1 = k01 * (-e1 / (r * t)).exp();
                    let k2 = k02 * (-e2 / (r * t)).exp();

                    m1 *= (-k1 * dt).exp();
                    m2 *= (-k2 * dt).exp();

                    total.push(m1 + m2);
                }

                total
            }
        }
    }

    fn apply_noise(data: &mut [f64], cfg: &NoiseConfig, rng: &mut StdRng) {
        match &cfg.kind {
            NoiseKind::Gaussian { sigma } => {
                let dist = Normal::new(0.0, *sigma).unwrap();
                for v in data {
                    *v += dist.sample(rng);
                }
            }

            NoiseKind::Uniform { amplitude } => {
                let dist = Uniform::new(-amplitude, *amplitude).unwrap();
                for v in data {
                    *v += dist.sample(rng);
                }
            }

            NoiseKind::Drift { slope } => {
                for (i, v) in data.iter_mut().enumerate() {
                    *v += slope * i as f64;
                }
            }
        }
    }

    fn apply_spikes(data: &mut [f64], spike: &SpikeModel, rng: &mut StdRng) {
        if spike.probability <= 0.0 {
            return;
        }

        let spike_dist = Normal::new(0.0, spike.amplitude).unwrap();

        for v in data.iter_mut() {
            let r: f64 = rng.random();
            if r < spike.probability {
                *v += spike_dist.sample(rng);
            }
        }
    }

    pub fn generate_series(cfg: &AdvancedTGAConfig) -> Vec<SimulatedTGA> {
        match &cfg.experiment_mode {
            ExperimentMode::Isothermal { temperatures } => temperatures
                .iter()
                .enumerate()
                .map(|(i, &temp)| Self::generate_one_isothermal(cfg, temp, i as u64))
                .collect(),

            ExperimentMode::NonIsothermal { t0, heating_rates } => heating_rates
                .iter()
                .enumerate()
                .map(|(i, &beta_min)| Self::generate_one_ramp(cfg, *t0, beta_min, i as u64))
                .collect(),
        }
    }
}
//=======================================================================================================
// EASY EXPONENTIAL DECAY WITH NOISE AND OPTIONAL SPIKES
//=======================================================================================================
#[derive(Clone, Debug)]
pub struct NoiseModel {
    pub sigma: f64, // стандартное отклонение
}

#[derive(Clone, Debug)]
pub struct SpikeModel {
    pub probability: f64, // вероятность выброса на точку
    pub amplitude: f64,   // характерная амплитуда
}

#[derive(Clone, Debug)]
pub struct VirtualTGAConfig {
    pub n_points: usize,
    pub dt: f64,

    pub temperature: f64,
    pub temp_noise: NoiseModel,

    pub m0: f64,
    pub k: f64,
    pub mass_noise: NoiseModel,
    pub spikes: Option<SpikeModel>,

    pub seed: u64,
}
//===========================================================================================
impl VirtualTGA {
    pub fn generate(cfg: &VirtualTGAConfig) -> Self {
        let mut rng = StdRng::seed_from_u64(cfg.seed);

        let time: Vec<f64> = (0..cfg.n_points).map(|i| i as f64 * cfg.dt).collect();

        let temp_noise = Normal::new(0.0, cfg.temp_noise.sigma).unwrap();
        let mass_noise = Normal::new(0.0, cfg.mass_noise.sigma).unwrap();

        let temperature: Vec<f64> = time
            .iter()
            .map(|_| cfg.temperature + temp_noise.sample(&mut rng))
            .collect();

        let mut mass: Vec<f64> = time
            .iter()
            .map(|&t| {
                let clean = cfg.m0 * (-cfg.k * t).exp();
                clean + mass_noise.sample(&mut rng)
            })
            .collect();

        if let Some(spike) = &cfg.spikes {
            let spike_noise = Normal::new(0.0, spike.amplitude).unwrap();
            for m in &mut mass {
                let r: f64 = rng.random();
                if r < spike.probability {
                    *m += spike_noise.sample(&mut rng);
                }
            }
        }

        Self {
            time,
            temperature,
            mass,
        }
    }

    pub fn write_txt(&self) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();

        file.write_all(b"time mass temperature\n").unwrap();

        for i in 0..self.time.len() {
            let line = format!(
                "{} {} {}\n",
                self.time[i], self.mass[i], self.temperature[i],
            );
            file.write_all(line.as_bytes()).unwrap();
        }

        file
    }
}
//================================================================================================================
/// wrapper builder for VirtualTGA
pub struct VirtualTGABuilder {
    cfg: VirtualTGAConfig,
}

impl VirtualTGABuilder {
    pub fn new(n_points: usize, dt: f64) -> Self {
        Self {
            cfg: VirtualTGAConfig {
                n_points,
                dt,
                temperature: 300.0,
                temp_noise: NoiseModel { sigma: 0.0 },
                m0: 1.0,
                k: 1e-3,
                mass_noise: NoiseModel { sigma: 0.0 },
                spikes: None,
                seed: 42,
            },
        }
    }

    pub fn temperature(mut self, t: f64) -> Self {
        self.cfg.temperature = t;
        self
    }

    pub fn mass_kinetics(mut self, m0: f64, k: f64) -> Self {
        self.cfg.m0 = m0;
        self.cfg.k = k;
        self
    }

    pub fn mass_noise(mut self, sigma: f64) -> Self {
        self.cfg.mass_noise = NoiseModel { sigma };
        self
    }

    pub fn spikes(mut self, probability: f64, amplitude: f64) -> Self {
        self.cfg.spikes = Some(SpikeModel {
            probability,
            amplitude,
        });
        self
    }

    pub fn seed(mut self, seed: u64) -> Self {
        self.cfg.seed = seed;
        self
    }

    pub fn build(self) -> VirtualTGA {
        VirtualTGA::generate(&self.cfg)
    }
}
//================================================================================================================
//  PIPELINR INVARIANT TEST
pub struct PipelineInvariantTest {
    pub df: DataFrame,
}

impl PipelineInvariantTest {
    pub fn from_lazy(frame: &LazyFrame) -> PolarsResult<Self> {
        Ok(Self {
            df: frame.clone().collect()?,
        })
    }

    pub fn no_nulls(&self, cols: &[&str]) {
        for &c in cols {
            let n = self.df.column(c).unwrap().null_count();
            assert_eq!(n, 0, "Column '{}' has {} nulls", c, n);
        }
    }

    pub fn same_length(&self, cols: &[&str]) {
        let len = self.df.height();
        for &c in cols {
            assert_eq!(
                self.df.column(c).unwrap().len(),
                len,
                "Column '{}' length mismatch",
                c
            );
        }
    }

    pub fn monotonic_increasing(&self, col: &str) {
        let s = self.df.column(col).unwrap().f64().unwrap();
        let mut prev = None;
        for v in s.into_no_null_iter() {
            if let Some(p) = prev {
                assert!(v > p, "Column '{}' not monotonic", col);
            }
            prev = Some(v);
        }
    }

    pub fn all_no_nulls(&self) {
        for s in self.df.columns() {
            let n = s.null_count();
            assert_eq!(n, 0, "Column '{}' has {} nulls", s.name(), n);
        }
    }

    pub fn min_length(&self, min: usize) {
        let h = self.df.height();
        assert!(h >= min, "DataFrame too short: {}, expected >= {}", h, min);
    }

    pub fn same_length_all(&self) {
        let len = self.df.height();
        for s in self.df.columns() {
            assert_eq!(s.len(), len, "Column '{}' length mismatch", s.name());
        }
    }
}

//===================================================================================
#[cfg(test)]
mod tests3 {
    use super::*;
    const R: f64 = 8.314;
    fn base_config_single() -> AdvancedTGAConfig {
        AdvancedTGAConfig {
            n_points: 1000,
            dt: 0.1,
            experiment_mode: ExperimentMode::NonIsothermal {
                t0: 300.0,
                heating_rates: vec![5.0],
            },
            kinetic_model: KineticModel::ArrheniusSingle {
                m0: 100.0,
                k0: 1e5,
                e: 80000.0,
                r: R,
            },
            mass_noise: None,
            temp_noise: None,
            spikes: None,

            seed: 42,
        }
    }

    #[test]
    fn mass_monotonic_without_noise() {
        let cfg = base_config_single();
        let series = VirtualTGA::generate_series(&cfg);
        let exp = &series[0];

        for i in 1..exp.virtual_tga.mass.len() {
            assert!(
                exp.virtual_tga.mass[i] <= exp.virtual_tga.mass[i - 1],
                "Mass is not monotonically decreasing"
            );
        }
    }

    #[test]
    fn temperature_linear_ramp() {
        let cfg = base_config_single();
        let series = VirtualTGA::generate_series(&cfg);
        let exp = &series[0];

        for i in 1..exp.virtual_tga.temperature.len() {
            assert!(
                exp.virtual_tga.temperature[i] > exp.virtual_tga.temperature[i - 1],
                "Temperature is not increasing"
            );
        }
    }
    #[test]
    fn reproducible_by_seed() {
        let cfg = base_config_single();

        let s1 = VirtualTGA::generate_series(&cfg);
        let s2 = VirtualTGA::generate_series(&cfg);

        assert_eq!(s1[0].virtual_tga.mass, s2[0].virtual_tga.mass);
        assert_eq!(s1[0].virtual_tga.temperature, s2[0].virtual_tga.temperature);
    }

    #[test]
    fn spikes_modify_signal() {
        let mut cfg = base_config_single();
        cfg.spikes = Some(SpikeModel {
            probability: 0.5,
            amplitude: 50.0,
        });

        let series = VirtualTGA::generate_series(&cfg);
        let exp = &series[0];

        let mut spike_detected = false;

        for i in 1..exp.virtual_tga.mass.len() {
            if (exp.virtual_tga.mass[i] - exp.virtual_tga.mass[i - 1]).abs() > 10.0 {
                spike_detected = true;
                break;
            }
        }

        assert!(spike_detected, "No spikes detected");
    }

    #[test]
    fn two_component_decay_shape() {
        let cfg = AdvancedTGAConfig {
            n_points: 1000,
            dt: 0.1,
            experiment_mode: ExperimentMode::NonIsothermal {
                t0: 300.0,
                heating_rates: vec![5.0],
            },
            kinetic_model: KineticModel::ArrheniusTwoComponent {
                m01: 50.0,
                k01: 1e6,
                e1: 90000.0,
                m02: 50.0,
                k02: 1e3,
                e2: 60000.0,
                r: R,
            },
            mass_noise: None,
            temp_noise: None,
            spikes: None,

            seed: 7,
        };

        let series = VirtualTGA::generate_series(&cfg);
        let exp = &series[0];

        let initial_drop = exp.virtual_tga.mass[10] - exp.virtual_tga.mass[0];
        let later_drop = exp.virtual_tga.mass[900] - exp.virtual_tga.mass[890];

        assert!(
            initial_drop.abs() > later_drop.abs(),
            "Two-stage behavior not observed"
        );
    }

    #[test]
    fn gaussian_noise_changes_variance() {
        let mut cfg = base_config_single();
        cfg.mass_noise = Some(NoiseConfig {
            kind: NoiseKind::Gaussian { sigma: 5.0 },
        });

        let noisy = VirtualTGA::generate_series(&cfg);

        cfg.mass_noise = None;
        let clean = VirtualTGA::generate_series(&cfg);

        let var_noisy: f64 = variance(&noisy[0].virtual_tga.mass);
        let var_clean: f64 = variance(&clean[0].virtual_tga.mass);

        assert!(var_noisy > var_clean);
    }

    fn variance(data: &[f64]) -> f64 {
        let mean = data.iter().sum::<f64>() / data.len() as f64;
        data.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / data.len() as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cut_interval_keeps_column_lengths_equal() {
        let cfg = VirtualTGAConfig {
            n_points: 10_000,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.2 },
            m0: 10.0,
            k: 1e-4,
            mass_noise: NoiseModel { sigma: 1e-3 },
            spikes: Some(SpikeModel {
                probability: 1e-4,
                amplitude: 0.1,
            }),
            seed: 42,
        };
    }

    #[test]
    fn virtual_tga_generate_creates_correct_number_of_points() {
        let cfg = VirtualTGAConfig {
            n_points: 100,
            dt: 0.1,
            temperature: 500.0,
            temp_noise: NoiseModel { sigma: 0.1 },
            m0: 5.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 1e-4 },
            spikes: None,
            seed: 123,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        println!("time {:?}", virtual_tga.time);

        assert_eq!(virtual_tga.time.len(), cfg.n_points);
        assert_eq!(virtual_tga.temperature.len(), cfg.n_points);
        assert_eq!(virtual_tga.mass.len(), cfg.n_points);
    }

    #[test]
    fn virtual_tga_write_txt_creates_file_with_data() {
        let cfg = VirtualTGAConfig {
            n_points: 50,
            dt: 0.2,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 0.05 },
            m0: 2.0,
            k: 2e-3,
            mass_noise: NoiseModel { sigma: 5e-5 },
            spikes: Some(SpikeModel {
                probability: 0.01,
                amplitude: 0.02,
            }),
            seed: 456,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let file = virtual_tga.write_txt();

        // Check that the file exists
        assert!(file.path().exists());

        // Read the content
        let content = std::fs::read_to_string(file.path()).unwrap();

        // Check header
        assert!(content.starts_with("time mass temperature\n"));

        // Check number of lines (header + data lines)
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), cfg.n_points + 1);

        // Check that data lines have three values each
        for line in &lines[1..] {
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(parts.len(), 3);
            // Optionally, check that they are parseable as f64
            for part in parts {
                assert!(part.parse::<f64>().is_ok());
            }
        }
    }

    #[test]
    fn virtual_tga_write_txt_data_matches_struct() {
        let cfg = VirtualTGAConfig {
            n_points: 10,
            dt: 1.0,
            temperature: 300.0,
            temp_noise: NoiseModel { sigma: 0.0 }, // No noise for exact match
            m0: 1.0,
            k: 0.1,
            mass_noise: NoiseModel { sigma: 0.0 },
            spikes: None,
            seed: 789,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let file = virtual_tga.write_txt();

        let content = std::fs::read_to_string(file.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        for (i, line) in lines.iter().skip(1).enumerate() {
            let parts: Vec<f64> = line
                .split_whitespace()
                .map(|s| s.parse().unwrap())
                .collect();
            assert_eq!(parts[0], virtual_tga.time[i]);
            assert_eq!(parts[1], virtual_tga.mass[i]);
            assert_eq!(parts[2], virtual_tga.temperature[i]);
        }
    }

    #[test]
    fn virtual_tga_generate_with_spikes_includes_spikes() {
        let cfg = VirtualTGAConfig {
            n_points: 1000,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.0 },
            m0: 10.0,
            k: 1e-4,
            mass_noise: NoiseModel { sigma: 0.0 },
            spikes: Some(SpikeModel {
                probability: 0.1,
                amplitude: 1.0,
            }),
            seed: 999,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);

        // Since spikes are added, mass should not be exactly the exponential decay
        let expected_mass: Vec<f64> = virtual_tga
            .time
            .iter()
            .map(|&t| cfg.m0 * (-cfg.k * t).exp())
            .collect();

        let mut has_spike = false;
        for (actual, expected) in virtual_tga.mass.iter().zip(expected_mass.iter()) {
            if (actual - expected).abs() > 1e-6 {
                // Allow small floating point differences
                has_spike = true;
                break;
            }
        }
        assert!(has_spike, "Expected at least one spike in the data");
    }
}

#[cfg(test)]
mod tests2 {
    use super::*;
    use crate::Kinetics::experimental_kinetics::exp_engine_api::ViewRange;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::Unit;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::{
        ds_from_csv, make_csv,
    };
    #[test]
    fn sample_columns_respects_max_points() {
        let csv = make_csv(100_000, 42);
        let ds = ds_from_csv(&csv);

        let series = ds.sample_columns("time", &["mass"], None, 2000).unwrap();

        assert_eq!(series.len(), 1);
        assert!(series[0].x.len() <= 2000);
        assert_eq!(series[0].x.len(), series[0].y.len());
    }
    #[test]
    fn sample_columns_respects_max_points2() {
        let csv = make_csv(100_000, 123);
        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap();

        let out = ds.sample_columns("time", &["mass"], None, 2000).unwrap();

        assert_eq!(out.len(), 1);

        let s = &out[0];
        assert!(s.x.len() <= 2000);
        assert_eq!(s.x.len(), s.y.len());
    }
    #[test]
    fn sample_columns_respects_range() {
        let csv = make_csv(10_000, 7);
        let ds = ds_from_csv(&csv);

        let series = ds
            .sample_columns(
                "time",
                &["mass"],
                Some(ViewRange {
                    t_min: 10.0,
                    t_max: 20.0,
                }),
                500,
            )
            .unwrap();

        let s = &series[0];
        assert!(!s.x.is_empty());

        assert!(s.x.iter().all(|&t| t >= 10.0 && t <= 20.0));
    }

    #[test]
    fn sample_columns_respects_view_range2() {
        let csv = make_csv(50_000, 999);
        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap();

        let out = ds
            .sample_columns(
                "time",
                &["mass"],
                Some(ViewRange {
                    t_min: 100.0,
                    t_max: 200.0,
                }),
                1000,
            )
            .unwrap();

        let s = &out[0];
        assert!(!s.x.is_empty());

        for &t in &s.x {
            assert!(t >= 100.0 && t <= 200.0);
        }
    }
}

#[cfg(test)]
mod synthetic_dataset_tests {
    use super::*;

    #[test]
    fn create_from_synthetic_data_sets_schema_fields() {
        let virtual_tga = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            mass: vec![10.0, 9.5, 9.0],
            temperature: vec![300.0, 301.0, 302.0],
        };

        let ds = TGADataset::create_from_synthetic_data(&virtual_tga).unwrap();

        assert_eq!(ds.schema.time.as_deref(), Some("time"));
        assert_eq!(ds.schema.mass.as_deref(), Some("mass"));
        assert_eq!(ds.schema.temperature.as_deref(), Some("temperature"));
        assert_eq!(ds.schema.columns.get("time").unwrap().unit, Unit::Second);
        assert_eq!(ds.schema.columns.get("mass").unwrap().unit, Unit::Milligram);
        assert_eq!(
            ds.schema.columns.get("temperature").unwrap().unit,
            Unit::Celsius
        );
        assert_eq!(
            ds.schema.columns.get("time").unwrap().origin,
            ColumnOrigin::Raw
        );
        assert_eq!(
            ds.schema.columns.get("mass").unwrap().origin,
            ColumnOrigin::Raw
        );
        assert_eq!(
            ds.schema.columns.get("temperature").unwrap().origin,
            ColumnOrigin::Raw
        );
    }

    #[test]
    fn create_from_synthetic_data_puts_time_mass_temperature_in_correct_columns() {
        let virtual_tga = VirtualTGA {
            time: vec![0.0, 0.5, 1.0, 1.5],
            mass: vec![5.0, 4.8, 4.6, 4.5],
            temperature: vec![295.0, 300.0, 305.0, 310.0],
        };

        let ds = TGADataset::create_from_synthetic_data(&virtual_tga).unwrap();
        let df = ds.frame.collect().unwrap();

        let time: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let temperature: Vec<f64> = df
            .column("temperature")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        assert_eq!(time, virtual_tga.time);
        assert_eq!(mass, virtual_tga.mass);
        assert_eq!(temperature, virtual_tga.temperature);
    }

    #[test]
    fn create_from_synthetic_data_rejects_mismatched_lengths() {
        let virtual_tga = VirtualTGA {
            time: vec![0.0, 1.0, 2.0],
            mass: vec![10.0, 9.0],
            temperature: vec![300.0, 301.0, 302.0],
        };

        let err = TGADataset::create_from_synthetic_data(&virtual_tga).unwrap_err();
        match err {
            TGADomainError::InvalidOperation(msg) => {
                assert!(msg.contains("equal lengths"));
            }
            _ => panic!("Unexpected error variant"),
        }
    }
}
