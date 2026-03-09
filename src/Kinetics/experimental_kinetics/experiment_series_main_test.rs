#[cfg(test)]
mod tests {
    use crate::Kinetics::experimental_kinetics::exp_engine_api::{
        GoldenPipelineConfig, SavGolConfig, SplineConfig,
    };
    use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::*;
    use crate::Kinetics::experimental_kinetics::experiment_series_main::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::TGADataset;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::*;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset_test::tests::{
        ds_from_csv, make_csv,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::*;
    #[test]
    fn golden_pipeline_experiment() {
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

        let exp = TGAExperiment::new(ds)
            .with_heating_rate(10.0)
            .with_comment("synthetic noisy test");

        let exp = exp // <TGADataset as Clone>::clone(&exp)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::MilliVolt)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            .calibrate_mass(1.0, 0.1, "mass")
            .unwrap()
            .trim_edges(50, 0)
            .celsius_to_kelvin()
            .seconds_to_hours()
            .hampel_filter("mass", 11, 3.0, HampelStrategy::ReplaceWithMedian)
            .unwrap()
            .sg_filter_column("mass", 11, 3, 0, 1.0)
            .unwrap()
            .rolling_mean("mass", 9)
            .trim_null_edges()
            .dimensionless_mass(0.0, 10.0, "alpha")
            .unwrap()
            .derive_mass_rate("dalpha_dt")
            .unwrap()
            .trim_null_edges();

        let df = exp.dataset.frame.collect().unwrap();

        assert!(df.height() > 1000);
        for col in df.columns() {
            assert_eq!(col.null_count(), 0);
        }
    }

    #[test]
    fn experiment_wrapper_mean_functions() {
        let csv = make_csv(1000, 456);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds);
        let _ = exp.mean_on_interval("mass", "time", 5.0, 10.0).unwrap();
        let _ = exp
            .mean_on_interval_on_own_range("mass", 5.0, 10.0)
            .unwrap();
        let _ = exp.mean_on_column("mass").unwrap();
    }

    #[test]
    fn series_wrapper_mean_functions() {
        let csv = make_csv(500, 789);
        let ds = ds_from_csv(&csv);
        let mut series = TGASeries::new();
        series.push(TGAExperiment::new(ds.clone()).with_id("a"));
        let _ = series
            .mean_on_interval("a", "mass", "time", 0.0, 100.0)
            .unwrap();
        let _ = series
            .mean_on_interval_on_own_range("a", "mass", 0.0, 100.0)
            .unwrap();
        let _ = series.mean_on_column("a", "mass").unwrap();
    }

    fn golden_cfg(save_to_new_experiment: bool, del_old_experiment: bool) -> GoldenPipelineConfig {
        GoldenPipelineConfig {
            k: 1.0,
            b: 0.0,
            time_cut_before: None,
            time_cut_after: None,
            sav_gol_config: SavGolConfig {
                window_size: 11,
                poly_degree: 3,
                deriv: 0,
                delta: 1.0,
            },
            averaging_time: 1.0,
            spline_config: SplineConfig {
                n_points: 300,
                n_internal_points: 200,
                degree: 3,
            },
            save_to_new_experiment,
            del_old_experiment,
        }
    }

    #[test]
    fn series_apply_golden_pipeline_keeps_source_when_configured() {
        let csv = make_csv(1200, 555001);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds).with_id("exp1");
        let mut series = TGASeries::new();
        series.push(exp);

        series
            .apply_golden_pipeline("exp1", golden_cfg(true, false))
            .unwrap();

        assert_eq!(series.len(), 2);
        assert!(series.index_by_id("exp1").is_ok());
        assert!(series.index_by_id("proceeded_exp1").is_ok());
    }

    #[test]
    fn series_apply_golden_pipeline_drops_source_when_configured() {
        let csv = make_csv(1200, 555002);
        let ds = ds_from_csv(&csv);
        let exp = TGAExperiment::new(ds).with_id("exp1");
        let mut series = TGASeries::new();
        series.push(exp);

        series
            .apply_golden_pipeline("exp1", golden_cfg(true, true))
            .unwrap();

        assert_eq!(series.len(), 1);
        assert!(series.index_by_id("exp1").is_err());
        assert!(series.index_by_id("proceeded_exp1").is_ok());
    }
}
