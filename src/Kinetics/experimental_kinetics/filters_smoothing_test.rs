#[cfg(test)]
mod tests {

    use crate::Kinetics::experimental_kinetics::exp_kinetics_column_manipulation::column_stats;
    use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{TGADataset, Unit};
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::{
        ds_from_csv, make_csv,
    };
    use crate::Kinetics::experimental_kinetics::splines::SplineKind;
    use crate::Kinetics::experimental_kinetics::testing_mod::{
        NoiseModel, PipelineInvariantTest, SpikeModel, VirtualTGA, VirtualTGAConfig,
    };
    use polars::prelude::*;
    //================================================================================================
    // ROLLING MEAN
    #[test]
    fn rolling_mean_smoothes_signal() {
        let csv = make_csv(100_000, 9876);
        let ds = ds_from_csv(&csv);
        // Use a relatively large window to ensure smoothing effect
        let ds_sm = ds.clone().rolling_mean("mass", 125);
        let df0 = ds.frame.collect().unwrap();
        let df_smooth = ds_sm.frame.collect().unwrap();

        let m0: Vec<f64> = df0
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let m_smooth: Vec<f64> = df_smooth
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        // lengths equal
        assert_eq!(m0.len(), m_smooth.len());
        // variance should not increase (rough heuristic)
        let mean0 = m0.iter().sum::<f64>() / m0.len() as f64;
        let var0 = m0.iter().map(|v| (v - mean0).powi(2)).sum::<f64>() / m0.len() as f64;
        let mean_smooth = m_smooth.iter().sum::<f64>() / m_smooth.len() as f64;
        let var_smooth = m_smooth
            .iter()
            .map(|v| (v - mean_smooth).powi(2))
            .sum::<f64>()
            / m_smooth.len() as f64;
        println!(
            "mean0 = {}, var0 = {}, mean1 = {}, var1 = {}",
            mean0, var0, mean_smooth, var_smooth
        );
        assert!(var_smooth <= var0 * 1.05); // allow small numerical tolerance
    }
    //================================================================================================
    //              HAMPEL/MAD
    #[test]
    fn hampel_filter_replace_with_median() {
        use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
        let csv = make_csv(1000, 777);
        let mut ds = ds_from_csv(&csv);

        // Add artificial outliers
        let df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        mass[100] = mass[100] * 10.0; // outlier
        mass[500] = mass[500] * 0.1; // outlier

        ds.frame = ds.frame.with_column(Series::new("mass".into(), mass).lit());

        let ds_filtered = ds
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap();
        let df_filtered = ds_filtered.frame.collect().unwrap();

        assert_eq!(df_filtered.height(), df.height());
        assert!(df_filtered.column("mass").is_ok());
    }

    #[test]
    fn hampel_filter_drop_outliers() {
        use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::HampelStrategy;
        let csv = make_csv(200, 888);
        let mut ds = ds_from_csv(&csv);

        let df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        mass[50] = mass[50] * 20.0; // strong outlier

        ds.frame = ds.frame.with_column(Series::new("mass".into(), mass).lit());

        let ds_filtered = ds
            .hampel_filter("mass", 7, 2.0, HampelStrategy::Drop)
            .unwrap();
        let df_filtered = ds_filtered.frame.collect().unwrap();

        assert!(df_filtered.height() < df.height());
    }

    #[test]
    fn tga_hampel_replace_keeps_length_and_damps_outliers() {
        let csv = make_csv(1000, 777);
        let mut ds = ds_from_csv(&csv);

        // eager modification (SAFE)
        let mut df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        let original = mass[100];
        mass[100] *= 10.0;

        df.replace("mass", Series::new("mass".into(), mass))
            .unwrap();
        ds.frame = df.clone().lazy();

        let ds2 = ds
            .hampel_filter_null_safe("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap();

        let df2 = ds2.frame.collect().unwrap();

        assert_eq!(df2.height(), df.height());

        let filtered = df2.column("mass").unwrap().f64().unwrap().get(100).unwrap();

        assert!(filtered < original * 2.0);
    }

    #[test]
    fn tga_hampel_drop_reduces_rows_consistently() {
        let csv = make_csv(300, 888);
        let mut ds = ds_from_csv(&csv);

        let mut df = ds.frame.clone().collect().unwrap();
        let mut mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        mass[50] *= 20.0;

        df.replace("mass", Series::new("mass".into(), mass))
            .unwrap();
        ds.frame = df.clone().lazy();

        let ds2 = ds
            .hampel_filter_null_safe("mass", 7, 2.0, HampelStrategy::Drop)
            .unwrap();

        let df2 = ds2.frame.collect().unwrap();

        assert!(df2.height() < df.height());

        // column alignment invariant
        assert_eq!(
            df2.column("time").unwrap().len(),
            df2.column("mass").unwrap().len()
        );
    }
    //==========================================================================
    // SAVITSKY-GOLAY
    #[test]
    fn tga_sg_filter_smooths_mass() {
        let csv = make_csv(500, 42);
        let ds = ds_from_csv(&csv);

        let ds2 = ds.sg_filter_column("mass", 11, 3, 0, 1.0).unwrap();

        let df = ds2.frame.collect().unwrap();

        let mass = df.column("mass").unwrap().f64().unwrap();

        let roughness: f64 = mass
            .into_no_null_iter()
            .collect::<Vec<_>>()
            .windows(2)
            .map(|w| (w[1] - w[0]).abs())
            .sum();

        assert!(roughness.is_finite());
    }
    //==============================================================================
    // FULL PIPELINE

    #[test]
    fn pipeline_smoke_test() {
        let csv = make_csv(1000, 123);
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .trim_edges(5, 5)
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .assert_no_nulls("hampel")
            .unwrap()
            .rolling_mean("mass", 9)
            .trim_null_edges()
            .assert_no_nulls("rolling_mean")
            .unwrap()
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .derive_mass_rate("dm_dt")
            .unwrap();

        let df = ds2.frame.collect().unwrap();

        assert!(df.column("dm_dt").is_ok());
        assert!(df.height() > 10);
    }

    #[test]
    fn golden_tga_pipeline_produces_physical_dm_dt() {
        let cfg = VirtualTGAConfig {
            n_points: 10_0,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.2 },
            m0: 10.0,
            k: 0.01,
            mass_noise: NoiseModel { sigma: 0.02 },
            spikes: Some(SpikeModel {
                probability: 0.01,
                amplitude: 1.0,
            }),
            seed: 2024,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();
        let ds = ds_from_csv(&csv);

        let ds2 = ds
            .trim_edges(10, 10)
            .hampel_filter("mass", 15, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .assert_no_nulls("hampel")
            .unwrap()
            .sg_filter_column("mass", 15, 3, 0, 1.0)
            .unwrap()
            .check_nulls_for_operation("SG")
            .derive_mass_rate("dm_dt")
            .unwrap();

        let df = ds2.frame.collect().unwrap();
        let dm_dt = df.column("dm_dt").unwrap().f64().unwrap();
        let dm_dt_vec: Vec<f64> = dm_dt.clone().into_no_null_iter().collect();
        println!("{:?}", dm_dt_vec);
        // физический инвариант: dm/dt < 0
        let neg_fraction =
            dm_dt.into_no_null_iter().filter(|v| *v < 0.0).count() as f64 / df.height() as f64;

        assert!(neg_fraction > 0.9);
    }

    #[test]
    fn golden_pipeline_smoke_test() {
        let cfg = VirtualTGAConfig {
            n_points: 5000,
            dt: 1.0,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 2.0 },
            m0: 100.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 0.5 },
            spikes: Some(SpikeModel {
                probability: 0.002,
                amplitude: 20.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();
        let ds = ds_from_csv(&csv);

        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            // калибровка
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            // отсечь начальный мусор
            .trim_edges(50, 0)
            // единицы
            .celsius_to_kelvin()
            .seconds_to_hours()
            // выбросы
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            // сглаживание
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            // синхронизация
            .trim_null_edges()
            // степень превращения
            .dimensionless_mass(0.0, 10., "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            // финальная чистка
            .trim_null_edges();

        let df = ds.frame.collect().unwrap();

        // --- Инварианты ---
        assert!(df.height() > 1000);
        for col in df.get_columns() {
            assert_eq!(col.null_count(), 0);
        }
    }

    #[test]
    fn golden_pipeline_test2() {
        // --- Virtual experiment ---

        let cfg = VirtualTGAConfig {
            n_points: 5000,
            dt: 1.0,
            temperature: 400.0,
            temp_noise: NoiseModel { sigma: 2.0 },
            m0: 100.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 0.5 },
            spikes: Some(SpikeModel {
                probability: 0.002,
                amplitude: 20.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);

        // --- Full pipeline ---
        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            // calibration
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            // trim startup
            .trim_edges(50, 0)
            // units
            .celsius_to_kelvin()
            .seconds_to_hours()
            // outliers
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            // smoothing
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            // synchronization after smoothing
            .trim_null_edges()
            // kinetics
            .dimensionless_mass(0.0, 10.0, "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            // final cleanup
            .trim_null_edges();
        assert!(ds.get_dalpha_dt().is_ok());
        assert!(ds.get_alpha().is_ok());
        // --- Invariants ---
        let inv = PipelineInvariantTest::from_lazy(&ds.frame).unwrap();

        inv.same_length(&["time", "mass", "temperature", "alpha", "dalpha_dt"]);
        inv.no_nulls(&["time", "mass", "temperature", "alpha", "dalpha_dt"]);
        inv.monotonic_increasing("time");

        assert!(inv.df.height() > 1000);
    }

    #[test]
    fn spline_smooth_removes_spikes() {
        let cfg = VirtualTGAConfig {
            n_points: 50_000,
            dt: 1.0,
            temperature: 500.0,
            temp_noise: NoiseModel { sigma: 0.5 },
            m0: 10.0,
            k: 15e-5,
            mass_noise: NoiseModel { sigma: 0.05 },
            spikes: Some(SpikeModel {
                probability: 0.02,
                amplitude: 2.0,
            }),
            seed: 42,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();

        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = ds_from_csv(&csv);

        let ds = ds
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap()
            .trim_edges(10, 10)
            .hampel_filter_null_safe("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .trim_null_edges()
            .assert_no_nulls("hampel")
            .unwrap()
            .spline_resample_columns("time", "new_time", &["mass"], 300, SplineKind::Cosine)
            .unwrap();
        let inv = PipelineInvariantTest::from_lazy(&ds.frame).unwrap();
        let mass = ds.get_mass().unwrap();
        println!("mass {:?}", mass);
        let df = ds.frame.collect().unwrap();
        let stats = column_stats(&df, "mass").unwrap();
        println!("columstats {:?}", stats);
        for col in df.get_columns() {
            assert_eq!(col.null_count(), 0);
        }
        // --- Invariants ---
    }
}
