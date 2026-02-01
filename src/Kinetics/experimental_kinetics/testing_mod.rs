use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use std::io::Write;
use tempfile::NamedTempFile;

#[derive(Clone, Debug)]
pub struct NoiseModel {
    pub sigma: f64, // стандартное отклонение
}

#[derive(Clone, Debug)]
pub struct SpikeModel {
    pub probability: f64, // вероятность выброса на точку
    pub amplitude: f64,   // характерная амплитуда
}

pub struct VirtualTGA {
    pub time: Vec<f64>,
    pub temperature: Vec<f64>,
    pub mass: Vec<f64>,
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
                if rng.random::<f64>() < spike.probability {
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

        // дальше: normalize_txt → TGADataset → cut_interval → assert lengths
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
