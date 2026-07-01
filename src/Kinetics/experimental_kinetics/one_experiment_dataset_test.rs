#[cfg(test)]
pub mod tests {
    use crate::Kinetics::experimental_kinetics::column_provenance::ColumnProvenance;
    use crate::Kinetics::experimental_kinetics::exp_engine_api::XY;
    use crate::Kinetics::experimental_kinetics::one_experiment_dataset::*;

    // use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::*;
    use crate::Kinetics::experimental_kinetics::testing_mod::*;
    use polars::prelude::*;

    pub fn make_csv(n_points: usize, seed: u64) -> tempfile::NamedTempFile {
        let cfg = VirtualTGAConfig {
            n_points,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.1 },
            m0: 10.0,
            k: 1e-4,
            mass_noise: NoiseModel { sigma: 1e-3 },
            spikes: None,
            seed,
        };

        let virtual_tga = VirtualTGA::generate(&cfg);
        let txt = virtual_tga.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();
        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();
        csv
    }

    pub fn ds_from_csv(csv: &tempfile::NamedTempFile) -> TGADataset {
        TGADataset::from_csv(csv.path().to_str().unwrap(), "time", "temperature", "mass").unwrap()
    }

    #[test]
    fn normalize_txt_to_csv_writes_header_and_rows() {
        let csv = make_csv(1_000, 42);
        let content = std::fs::read_to_string(csv.path()).unwrap();
        let mut lines = content.lines();
        assert_eq!(lines.next().unwrap(), "time,mass,temperature");
        // header + data lines
        assert_eq!(content.lines().count(), 1_000 + 1);
        // check parseability of a data line (skip header)
        let parts: Vec<&str> = lines.next().unwrap().split(',').collect();
        assert_eq!(parts.len(), 3);
        for p in parts {
            p.parse::<f64>().unwrap();
        }
    }

    #[test]
    fn from_csv_sets_schema_and_can_collect() {
        let csv = make_csv(500, 1);
        let ds = ds_from_csv(&csv);
        assert_eq!(ds.schema.time.clone().unwrap(), "time");
        assert_eq!(ds.schema.temperature.clone().unwrap(), "temperature");
        assert_eq!(ds.schema.mass.clone().unwrap(), "mass");
        assert_eq!(ds.schema.columns.get("time").unwrap().unit, Unit::Second);
        assert_eq!(
            ds.schema.columns.get("temperature").unwrap().unit,
            Unit::Celsius
        );
        assert_eq!(ds.schema.columns.get("mass").unwrap().unit, Unit::Milligram);

        let df = ds.frame.collect().unwrap();
        assert_eq!(df.height(), 500);
        assert!(df.column("time").is_ok());
        assert!(df.column("mass").is_ok());
        assert!(df.column("temperature").is_ok());
    }

    #[test]
    fn sort_by_bound_time_keeps_row_alignment() {
        let csv = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            csv.path(),
            "time,mass,temperature\n0.0,10.0,100.0\n2.0,20.0,200.0\n1.0,15.0,150.0\n",
        )
        .unwrap();

        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap()
            .sort_by_bound_time()
            .unwrap();

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
        let temp: Vec<f64> = df
            .column("temperature")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        assert_eq!(time, vec![0.0, 1.0, 2.0]);
        assert_eq!(mass, vec![10.0, 15.0, 20.0]);
        assert_eq!(temp, vec![100.0, 150.0, 200.0]);
    }

    #[test]
    fn time_monotonicity_report_detects_decrease_and_preserves_old_contract() {
        let csv = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            csv.path(),
            "time,mass,temperature\n0.0,10.0,100.0\n2.0,20.0,200.0\n1.0,15.0,150.0\n",
        )
        .unwrap();

        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap();

        let report = ds.time_monotonicity_report().unwrap();
        assert_eq!(report.time_column, "time");
        assert_eq!(report.len, 3);
        assert!(report.sortable_without_data_loss);
        assert!(!report.is_strictly_increasing);
        assert!(report.null_indices.is_empty());
        assert!(report.nan_indices.is_empty());
        assert!(report.infinite_indices.is_empty());
        assert_eq!(report.decreasing_pairs.len(), 1);
        assert_eq!(report.decreasing_pairs[0].index, 2);
        assert_eq!(report.decreasing_pairs[0].previous, 2.0);
        assert_eq!(report.decreasing_pairs[0].current, 1.0);
        assert_eq!(ds.monotony_of_time_check().unwrap(), vec![1.0]);
    }

    #[test]
    fn sort_by_bound_time_rejects_invalid_time_values() {
        let csv = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(
            csv.path(),
            "time,mass,temperature\n0.0,10.0,100.0\nNaN,20.0,200.0\n1.0,15.0,150.0\n",
        )
        .unwrap();

        let ds = ds_from_csv(&csv)
            .bind_time("time", Unit::Second)
            .unwrap()
            .bind_mass("mass", Unit::Milligram)
            .unwrap()
            .bind_temperature("temperature", Unit::Celsius)
            .unwrap();

        let err = ds.sort_by_bound_time().unwrap_err();
        match err {
            TGADomainError::InvalidOperation(msg) => {
                assert!(msg.contains("NaNs"));
                assert!(msg.contains("time"));
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

    #[test]
    fn cut_interval_preserves_alignment() {
        let csv = make_csv(1_000, 42);
        let ds = ds_from_csv(&csv);

        let ds2 = ds.cut_time_interval(10.0, 20.0);
        let df = ds2.frame.collect().unwrap();

        let n = df.height();
        assert!(n > 0);
        assert_eq!(df.column("time").unwrap().len(), n);
        assert_eq!(df.column("mass").unwrap().len(), n);
        assert_eq!(df.column("temperature").unwrap().len(), n);
    }

    #[test]
    fn trim_range_limits_time_inclusive() {
        let csv = make_csv(2_000, 7);
        let ds = ds_from_csv(&csv);
        let ds_trim = ds.trim_range("time", 5.0, 50.0);
        let df = ds_trim.frame.collect().unwrap();
        assert!(df.height() > 0);
        let time: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let min_t = time.iter().copied().fold(f64::INFINITY, f64::min);
        let max_t = time.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(min_t >= 5.0 - 1e-12);
        assert!(max_t <= 50.0 + 1e-12);
    }

    #[test]
    fn trim_range_inverse_removes_time_interval_inclusive() {
        let csv = make_csv(2_000, 70);
        let ds = ds_from_csv(&csv);
        let ds_trim = ds.trim_range_inverse("time", 5.0, 50.0);
        let df = ds_trim.frame.collect().unwrap();
        assert!(df.height() > 0);
        let time: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert!(time.iter().all(|&t| t < 5.0 || t > 50.0));
    }

    #[test]
    fn cut_range_inverse_x_or_y_removes_selected_axis_interval() {
        let csv = make_csv(2_000, 71);
        let ds = ds_from_csv(&csv)
            .set_oneframeplot_x("time")
            .unwrap()
            .set_oneframeplot_y("mass")
            .unwrap()
            .cut_range_inverse_x_or_y(XY::X, 5.0, 50.0)
            .unwrap();

        let x = ds.get_x_as_vec().unwrap();
        let y = ds.get_y_as_vec().unwrap();
        assert_eq!(x.len(), y.len());
        assert!(x.iter().all(|&t| t < 5.0 || t > 50.0));
    }

    #[test]
    fn scale_columns_multiplies_mass() {
        let csv = make_csv(500, 9);
        let ds_base = ds_from_csv(&csv);
        let ds_scaled = ds_from_csv(&csv).scale_columns(&["mass"], 2.0);

        let df_base = ds_base.frame.collect().unwrap();
        let df_scaled = ds_scaled.frame.collect().unwrap();

        let m0: Vec<f64> = df_base
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let m1: Vec<f64> = df_scaled
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert_eq!(m0.len(), m1.len());
        for (a, b) in m0.iter().zip(m1.iter()) {
            assert!((b - 2.0 * a).abs() < 1e-9);
        }
    }

    #[test]
    fn offset_column_shifts_time() {
        let csv = make_csv(300, 10);
        let ds_base = ds_from_csv(&csv);
        let ds_off = ds_from_csv(&csv).offset_column("time", 5.0);

        let df0 = df_collect_col(&ds_base, "time");
        let df1 = df_collect_col(&ds_off, "time");
        assert_eq!(df0.len(), df1.len());
        for (a, b) in df0.iter().zip(df1.iter()) {
            assert!((b - (a + 5.0)).abs() < 1e-12);
        }
    }

    #[test]
    fn celsius_to_kelvin_updates_schema_and_values() {
        let csv = make_csv(400, 11);
        let ds_base = ds_from_csv(&csv);
        let ds_kelvin = ds_from_csv(&csv).celsius_to_kelvin();
        let temp_col = ds_kelvin.schema.temperature.clone().unwrap();
        assert_eq!(
            ds_kelvin.schema.columns.get(&temp_col).unwrap().unit,
            Unit::Kelvin
        );
        assert_eq!(
            ds_kelvin.schema.columns.get(&temp_col).unwrap().origin,
            ColumnOrigin::PolarsDerived
        );

        let t0 = df_collect_col(&ds_base, "temperature");
        let t1 = df_collect_col(&ds_kelvin, "temperature");
        assert_eq!(t0.len(), t1.len());
        for (a, b) in t0.iter().zip(t1.iter()) {
            assert!((b - (a + 273.15)).abs() < 1e-9);
        }
    }

    #[test]
    fn seconds_to_hours_updates_schema_and_values() {
        let csv = make_csv(400, 12);
        let ds_base = ds_from_csv(&csv);
        let ds_hr = ds_from_csv(&csv).seconds_to_hours();
        let time_col = ds_hr.schema.time.clone().unwrap();
        assert_eq!(
            ds_hr.schema.columns.get(&time_col).unwrap().unit,
            Unit::Hour
        );
        assert_eq!(
            ds_hr.schema.columns.get(&time_col).unwrap().origin,
            ColumnOrigin::PolarsDerived
        );

        let t0 = df_collect_col(&ds_base, "time");
        let t1 = df_collect_col(&ds_hr, "time");
        for (a, b) in t0.iter().zip(t1.iter()) {
            assert!((b - (a / 3600.0)).abs() < 1e-12);
        }
    }

    #[test]
    fn seconds_to_hours_is_idempotent_when_time_is_already_in_hours() {
        let csv = make_csv(400, 12);
        let ds_hr = ds_from_csv(&csv).seconds_to_hours();
        let once = df_collect_col(&ds_hr, "time");
        let twice = ds_hr.clone().seconds_to_hours();
        let after = df_collect_col(&twice, "time");

        assert_eq!(twice.schema.columns.get("time").unwrap().unit, Unit::Hour);
        assert_eq!(once, after);
    }

    #[test]
    fn calibrate_mass_from_voltage_applies_linear_transform_and_updates_meta() {
        let csv = make_csv(350, 13);
        let ds_base = ds_from_csv(&csv);
        let ds_cal = ds_from_csv(&csv).calibrate_mass_from_voltage(3.0, 5.0);
        let mass_col = ds_cal.schema.mass.clone().unwrap();
        assert_eq!(
            ds_cal.schema.columns.get(&mass_col).unwrap().unit,
            Unit::Milligram
        );
        assert_eq!(
            ds_cal.schema.columns.get(&mass_col).unwrap().origin,
            ColumnOrigin::PolarsDerived
        );

        let m0 = df_collect_col(&ds_base, "mass");
        let m1 = df_collect_col(&ds_cal, "mass");
        for (a, b) in m0.iter().zip(m1.iter()) {
            assert!((b - (3.0 * a + 5.0)).abs() < 1e-9);
        }
    }

    #[test]
    fn materialize_returns_expected_vectors() {
        let csv = make_csv(123, 21);
        let ds = ds_from_csv(&csv);
        let nb = ds.materialize("time", &["mass", "temperature"]).unwrap();
        assert_eq!(nb.grid.len(), 123);
        assert!(nb.columns.contains_key("mass"));
        assert!(nb.columns.contains_key("temperature"));
        assert_eq!(nb.columns["mass"].len(), 123);
        assert_eq!(nb.columns["temperature"].len(), 123);

        // spot check equality with collected DataFrame values
        let df = ds_from_csv(&csv).frame.collect().unwrap();
        let time: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert!((nb.grid[0] - time[0]).abs() < 1e-12);
        assert!((nb.grid[122] - time[122]).abs() < 1e-12);
    }

    #[test]
    fn add_numeric_column_inserts_column_and_is_retrievable() {
        let csv = make_csv(200, 30);
        let ds = ds_from_csv(&csv);
        // build a deterministic vector: index as f64
        let df0 = ds.frame.clone().collect().unwrap();
        let n = df0.height();
        let data: Arc<Vec<f64>> = Arc::new((0..n).map(|i| i as f64).collect());

        let ds2 = ds
            .add_numeric_column(
                "row_index",
                Unit::Dimensionless,
                data,
                ColumnNature::Unknown,
            )
            .unwrap();

        // column is present and materializable
        let nb = ds2.materialize("time", &["row_index"]).unwrap();
        assert_eq!(nb.grid.len(), n);
        let idx = nb.columns.get("row_index").unwrap();
        assert_eq!(idx.len(), n);
        assert!((idx[0] - 0.0).abs() < 1e-12);
        assert!((idx[n - 1] - (n as f64 - 1.0)).abs() < 1e-12);

        // schema updated
        let meta = ds2.schema.columns.get("row_index").unwrap();
        assert_eq!(meta.name, "row_index");
        assert_eq!(meta.unit, Unit::Dimensionless);
        assert_eq!(meta.origin, ColumnOrigin::NumericDerived);
    }

    #[test]
    fn cut_before_time_keeps_min_time() {
        let csv = make_csv(600, 77);
        let ds = ds_from_csv(&csv);
        let t_start = 12.3;
        let df = ds.cut_before_time(t_start).frame.collect().unwrap();
        let t: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let min_t = t.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(min_t >= t_start - 1e-12);
    }

    #[test]
    fn filter_by_mask_applies_predicate() {
        let csv = make_csv(500, 88);
        let ds = ds_from_csv(&csv);
        let th = 15.0; // seconds
        let ds2 = ds.filter_by_mask(col("time").lt(lit(th)));
        let df = ds2.frame.collect().unwrap();
        let t: Vec<f64> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert!(t.iter().all(|&x| x < th + 1e-12));
        assert!(df.height() > 0);
    }

    #[test]
    fn with_column_expr_adds_and_computes_column() {
        let csv = make_csv(250, 99);
        let ds = ds_from_csv(&csv);
        let meta = ColumnMeta::new(
            "temp_plus10",
            Unit::Celsius,
            ColumnOrigin::PolarsDerived,
            ColumnNature::Temperature,
            ColumnProvenance::manual(
                "temp_plus10",
                "with_column_expr",
                Some("temperature + 10".to_string()),
            ),
        );
        let ds2 = ds
            .with_column_expr(meta.clone(), col("temperature") + lit(10.0))
            .unwrap();
        // schema contains the derived column
        assert!(ds2.schema.columns.contains_key("temp_plus10"));

        let df = ds2.frame.collect().unwrap();
        let t: Vec<f64> = df
            .column("temperature")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let t10: Vec<f64> = df
            .column("temp_plus10")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (a, b) in t.iter().zip(t10.iter()) {
            assert!((b - (a + 10.0)).abs() < 1e-12);
        }
    }

    // helper to collect a single f64 column from a lazy dataset
    fn df_collect_col(ds: &TGADataset, name: &str) -> Vec<f64> {
        ds.frame
            .clone()
            .collect()
            .unwrap()
            .column(name)
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect()
    }

    #[test]
    fn test_dataset_correctness() {
        let cfg = VirtualTGAConfig {
            n_points: 2000,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.1 },
            m0: 10.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 1e-4 },
            spikes: None,
            seed: 1,
        };

        let v = VirtualTGA::generate(&cfg);
        let txt = v.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();

        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = TGADataset::from_csv(csv.path().to_str().unwrap(), "time", "temperature", "mass")
            .unwrap();
        let schema = &ds.schema;
        println!("schema {:?}", schema);
        let temp = schema.temperature.clone().unwrap();
        let mass = schema.mass.clone().unwrap();
        println!("names {}, {}", temp, mass);
        let termo_meta = ds.schema.columns.get(&temp).unwrap();
        let mass_meta = ds.schema.columns.get(&mass).unwrap();
        println!("{:?}, {:?}", termo_meta, mass_meta);
    }

    #[test]
    fn mass_rate_is_negative_on_average() {
        let cfg = VirtualTGAConfig {
            n_points: 2000,
            dt: 0.1,
            temperature: 600.0,
            temp_noise: NoiseModel { sigma: 0.1 },
            m0: 10.0,
            k: 1e-3,
            mass_noise: NoiseModel { sigma: 1e-4 },
            spikes: None,
            seed: 1,
        };

        let v = VirtualTGA::generate(&cfg);
        let txt = v.write_txt();
        let csv = tempfile::NamedTempFile::new().unwrap();

        TGADataset::normalize_txt_to_csv(txt.path(), csv.path()).unwrap();

        let ds = TGADataset::from_csv(csv.path().to_str().unwrap(), "time", "temperature", "mass")
            .unwrap();

        let ds = ds.derive_mass_rate("dm_dt").unwrap();
        let df = ds.frame.collect().unwrap();

        let mean = df.column("dm_dt").unwrap().f64().unwrap().mean().unwrap();

        assert!(mean < 0.0);
    }

    #[test]
    fn mean_on_interval_matches_manual_average() {
        let csv = make_csv(1000, 123);
        let ds = ds_from_csv(&csv);
        let df = ds.frame.clone().collect().unwrap();
        let from = 5.0;
        let to = 10.0;
        let m = ds.mean_on_interval("mass", "time", from, to).unwrap();

        let mask: Vec<bool> = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(|t| t >= from && t <= to)
            .collect();
        let mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let mut sum = 0.0;
        let mut cnt = 0usize;
        for (v, keep) in mass.into_iter().zip(mask.into_iter()) {
            if keep {
                sum += v;
                cnt += 1;
            }
        }
        let avg = sum / cnt as f64;
        assert!((m - avg).abs() < 1e-9);
    }

    #[test]
    fn mean_on_interval_on_own_range_matches_manual() {
        let csv = make_csv(1000, 123);
        let ds = ds_from_csv(&csv);
        let df = ds.frame.clone().collect().unwrap();
        let from = 5.0;
        let to = 10.0;
        let m = ds.mean_on_interval_on_own_range("mass", from, to).unwrap();

        let mask: Vec<bool> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(|v| v >= from && v <= to)
            .collect();
        let mass: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let mut sum = 0.0;
        let mut cnt = 0usize;
        for (v, keep) in mass.into_iter().zip(mask.into_iter()) {
            if keep {
                sum += v;
                cnt += 1;
            }
        }
        let avg = sum / cnt as f64;
        assert!((m - avg).abs() < 1e-9);
    }

    #[test]
    fn mean_on_column_matches_manual() {
        let csv = make_csv(1000, 123);
        let ds = ds_from_csv(&csv);
        let df = ds.frame.clone().collect().unwrap();
        let m = ds.mean_on_column("mass").unwrap();
        let manual = df.column("mass").unwrap().f64().unwrap().mean().unwrap();
        assert_eq!(m, manual);
    }

    #[test]
    fn exp_and_ln_columns_added_and_correct() {
        let csv = make_csv(200, 321);
        let ds = ds_from_csv(&csv);

        // shift mass positive if needed (generator should already be >0)
        let ds_ln = ds.clone().ln_column("mass").unwrap();
        let ds_exp = ds.clone().exp_column("mass").unwrap();

        // ln_mass exists and equals ln(mass)
        let df_ln = ds_ln.frame.collect().unwrap();
        let mass: Vec<f64> = df_ln
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let ln_mass: Vec<f64> = df_ln
            .column("ln_mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (m, l) in mass.iter().zip(ln_mass.iter()) {
            assert!((l - m.ln()).abs() < 1e-9);
        }

        // exp_mass exists and equals exp(mass)
        let df_exp = ds_exp.frame.collect().unwrap();
        let mass2: Vec<f64> = df_exp
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let exp_mass: Vec<f64> = df_exp
            .column("exp_mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (m, e) in mass2.iter().zip(exp_mass.iter()) {
            assert!((e - m.exp()).abs() < 1e-9);
        }
    }

    #[test]
    fn dimensionless_and_conversion_columns() {
        let csv = make_csv(800, 456);
        let ds = ds_from_csv(&csv);
        // Choose averaging window on early time where mass ~ m0
        let ds_dim = ds
            .derive_dimensionless_mass(0.0, 5.0)
            .unwrap()
            .derive_conversion()
            .unwrap();

        // check schema has both
        assert!(ds_dim.schema.columns.contains_key("alpha"));
        assert!(ds_dim.schema.columns.contains_key("eta"));

        let df = ds_dim.frame.collect().unwrap();
        let m_dim: Vec<f64> = df
            .column("alpha")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let conv: Vec<f64> = df
            .column("eta")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (z, c) in m_dim.iter().zip(conv.iter()) {
            assert!((c - (1.0 - z)).abs() < 1e-12);
            assert!((*z + *c - 1.0).abs() < 1e-12);
        }
        assert_eq!(
            ds_dim.schema.columns.get("alpha").unwrap().unit,
            Unit::Dimensionless
        );
        assert_eq!(
            ds_dim.schema.columns.get("eta").unwrap().unit,
            Unit::Dimensionless
        );
    }

    #[test]
    fn dimensionless_mass_rejects_empty_reference_window() {
        let csv = make_csv(200, 33);
        let ds = ds_from_csv(&csv);
        let err = ds.derive_dimensionless_mass(0.0, 0.0).unwrap_err();
        assert!(matches!(err, TGADomainError::InvalidConversionRange));
    }

    #[test]
    fn conversion_rejects_empty_reference_window() {
        let csv = make_csv(200, 34);
        let ds = ds_from_csv(&csv);
        let ds_dim = ds.derive_dimensionless_mass(0.0, 20.0).unwrap();
        let err = ds_dim.conversion(0.0, 0.0, "conv").unwrap_err();
        assert!(matches!(err, TGADomainError::InvalidConversionRange));
    }

    #[test]
    fn derive_temperature_rate_basic_properties() {
        let csv = make_csv(600, 789);
        let ds = ds_from_csv(&csv);
        let ds = ds.derive_temperature_rate("dT_dt").unwrap();
        assert_eq!(
            ds.schema.columns.get("dT_dt").unwrap().unit,
            Unit::CelsiusPerSecond
        );
        let df = ds.frame.collect().unwrap();
        let dTdt = df.column("dT_dt").unwrap().f64().unwrap();
        let dTdt_mean_polars = dTdt.mean().unwrap();
        assert!(dTdt.mean().unwrap().abs() < 1e-2);
        // let dTdt_vec: Vec<f64> = dTdt.into_no_null_iter().collect();
        // let dTdt_mean_custom: f64 = dTdt_vec.iter().sum::<f64>() / dTdt.len() as f64;
        println!("dTdt mean ={}", dTdt_mean_polars);
        //  assert!( (dTdt_mean_custom - dTdt_mean_polars).abs() <1e-4);
        // temperature should not increase because there is isothermic regime
    }

    #[test]
    fn derive_conversion_rate_uses_time_unit_in_metadata() {
        let csv = make_csv(600, 790);
        let ds_alpha = ds_from_csv(&csv)
            .seconds_to_hours()
            .derive_dimensionless_mass(0.0, 5.0)
            .unwrap()
            .derive_dalpha_dt("dalpha_dt")
            .unwrap();
        assert_eq!(
            ds_alpha.schema.columns.get("dalpha_dt").unwrap().unit,
            Unit::PerHour
        );
        assert_eq!(
            ds_alpha.schema.columns.get("dalpha_dt").unwrap().nature,
            ColumnNature::DimensionlessMassRate
        );

        let ds_eta = ds_from_csv(&csv)
            .seconds_to_hours()
            .derive_dimensionless_mass(0.0, 5.0)
            .unwrap()
            .derive_conversion()
            .unwrap()
            .derive_deta_dt("deta_dt")
            .unwrap();
        assert_eq!(
            ds_eta.schema.columns.get("deta_dt").unwrap().unit,
            Unit::PerHour
        );
        assert_eq!(
            ds_eta.schema.columns.get("deta_dt").unwrap().nature,
            ColumnNature::ConversionRate
        );

        let dalpha_dt: Vec<f64> = ds_alpha
            .frame
            .clone()
            .collect()
            .unwrap()
            .column("dalpha_dt")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();
        let deta_dt: Vec<f64> = ds_eta
            .frame
            .clone()
            .collect()
            .unwrap()
            .column("deta_dt")
            .unwrap()
            .f64()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();

        let max_sign_gap = dalpha_dt
            .iter()
            .zip(deta_dt.iter())
            .filter(|(a, b)| a.is_finite() && b.is_finite())
            .map(|(a, b)| (a + b).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            max_sign_gap < 1e-10,
            "deta_dt should be the negation of dalpha_dt, got max gap {}",
            max_sign_gap
        );
    }

    #[test]
    fn to_csv_with_units_writes_expected_header_and_data() {
        let csv = make_csv(200, 4242);
        let ds = ds_from_csv(&csv);

        let out = tempfile::NamedTempFile::new().unwrap();
        ds.to_csv_with_units(out.path()).unwrap();

        let content = std::fs::read_to_string(out.path()).unwrap();
        let mut lines = content.lines();
        // header line 1: meta tag
        assert_eq!(lines.next().unwrap(), "# KiThe TGA Dataset");
        // header line 2: units for each column
        let meta_line = lines.next().unwrap();
        assert!(meta_line.starts_with("# "));
        // should contain tokens like "time [s]", "mass [mg]", "temperature [C]" (order may vary)
        assert!(meta_line.contains("time [s]"));
        assert!(meta_line.contains("mass [mg]"));
        assert!(meta_line.contains("temperature [C]") || meta_line.contains("temperature [K]"));
        // then CSV header row
        let header_row = lines.next().unwrap();
        let cols: Vec<&str> = header_row.split(',').collect();
        assert!(cols.contains(&"time"));
        assert!(cols.contains(&"mass"));
        assert!(cols.contains(&"temperature"));
        // and at least one data row
        assert!(lines.next().is_some());
    }

    #[test]
    fn roundtrip_to_from_csv_with_units_preserves_schema_and_values() {
        let csv = make_csv(300, 5151);
        // modify schema to ensure units change appears in meta
        let ds = ds_from_csv(&csv).celsius_to_kelvin();

        let out = tempfile::NamedTempFile::new().unwrap();
        ds.to_csv_with_units(out.path()).unwrap();

        let ds2 = TGADataset::from_csv_with_units(out.path()).unwrap();
        // println!("ds2 {:?}", ds2);
        // compare schema units for primary columns
        let time_col = ds.schema.time.clone().unwrap();
        let temp_col = ds.schema.temperature.clone().unwrap();
        let mass_col = ds.schema.mass.clone().unwrap();
        assert_eq!(
            ds2.schema.columns.get(&time_col).unwrap().unit,
            ds.schema.columns.get(&time_col).unwrap().unit
        );
        assert_eq!(
            ds2.schema.columns.get(&temp_col).unwrap().unit,
            ds.schema.columns.get(&temp_col).unwrap().unit
        );
        assert_eq!(
            ds2.schema.columns.get(&mass_col).unwrap().unit,
            ds.schema.columns.get(&mass_col).unwrap().unit
        );

        // compare few values across columns
        let df1 = ds.frame.clone().collect().unwrap();
        let df2 = ds2.frame.clone().collect().unwrap();
        assert_eq!(df1.get_column_names(), df2.get_column_names());
        for name in df1.get_column_names() {
            let v1: Vec<f64> = df1
                .column(name)
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .take(5)
                .collect();
            let v2: Vec<f64> = df2
                .column(name)
                .unwrap()
                .f64()
                .unwrap()
                .into_no_null_iter()
                .take(5)
                .collect();
            for (a, b) in v1.iter().zip(v2.iter()) {
                assert!((a - b).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn from_csv_universal_reads_units_format_via_metadata_header() {
        let csv = make_csv(180, 7171);
        let ds = ds_from_csv(&csv).celsius_to_kelvin();

        let out = tempfile::NamedTempFile::new().unwrap();
        ds.to_csv_with_units(out.path()).unwrap();

        let ds2 = TGADataset::from_csv_universal(out.path()).unwrap();

        let time_col = ds.schema.time.clone().unwrap();
        let temp_col = ds.schema.temperature.clone().unwrap();
        let mass_col = ds.schema.mass.clone().unwrap();

        assert_eq!(
            ds2.schema.columns.get(&time_col).unwrap().unit,
            ds.schema.columns.get(&time_col).unwrap().unit
        );
        assert_eq!(
            ds2.schema.columns.get(&temp_col).unwrap().unit,
            ds.schema.columns.get(&temp_col).unwrap().unit
        );
        assert_eq!(
            ds2.schema.columns.get(&mass_col).unwrap().unit,
            ds.schema.columns.get(&mass_col).unwrap().unit
        );
    }

    #[test]
    fn from_csv_universal_reads_raw_csv_without_metadata_header() {
        let csv = make_csv(160, 8181);

        let ds_universal = TGADataset::from_csv_universal(csv.path()).unwrap();
        let ds_raw = TGADataset::from_csv_raw(csv.path()).unwrap();

        assert!(ds_universal.schema.time.is_none());
        assert!(ds_universal.schema.temperature.is_none());
        assert!(ds_universal.schema.mass.is_none());

        let df_u = ds_universal.frame.collect().unwrap();
        let df_r = ds_raw.frame.collect().unwrap();
        assert_eq!(df_u.shape(), df_r.shape());
        assert_eq!(df_u.get_column_names(), df_r.get_column_names());
    }

    #[test]
    fn rename_column_updates_schema_and_references() {
        let csv = make_csv(100, 111);
        let ds = ds_from_csv(&csv).rename_column("time", "t").unwrap();
        assert_eq!(ds.schema.time.unwrap(), "t");
        assert!(ds.schema.columns.contains_key("t"));
        assert!(!ds.schema.columns.contains_key("time"));
        let df = ds.frame.collect().unwrap();
        assert!(df.column("t").is_ok());
    }

    #[test]
    fn cut_mass_and_temperature_intervals() {
        let csv = make_csv(1000, 222);
        let ds = ds_from_csv(&csv);
        let df_orig = ds.frame.clone().collect().unwrap();
        let mass_vals: Vec<f64> = df_orig
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let temp_vals: Vec<f64> = df_orig
            .column("temperature")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        let mass_min = mass_vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let mass_max = mass_vals.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        let temp_min = temp_vals.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let temp_max = temp_vals.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

        let mass_cut = (mass_min + mass_max) / 2.0;
        let temp_cut = (temp_min + temp_max) / 2.0;

        let ds_mass_cut = ds.clone().cut_mass_interval(mass_cut - 0.1, mass_cut + 0.1);
        let ds_temp_cut = ds.cut_temperature_interval(temp_cut - 1.0, temp_cut + 1.0);

        assert!(ds_mass_cut.frame.collect().unwrap().height() < df_orig.height());
        assert!(ds_temp_cut.frame.collect().unwrap().height() < df_orig.height());
    }

    #[test]
    fn trim_edges_removes_specified_rows() {
        let csv = make_csv(100, 333);
        let ds = ds_from_csv(&csv);
        let orig_height = ds.frame.clone().collect().unwrap().height();
        let ds_trimmed = ds.trim_edges(5, 10);
        let new_height = ds_trimmed.frame.collect().unwrap().height();
        assert_eq!(new_height, orig_height - 15);
    }

    #[test]
    fn binding_rejects_physically_invalid_units() {
        let csv = make_csv(120, 777);

        let time_err = ds_from_csv(&csv)
            .bind_time("time", Unit::Milligram)
            .unwrap_err();
        assert!(matches!(
            time_err,
            TGADomainError::InvalidUnitForBinding {
                field: "time",
                unit: Unit::Milligram
            }
        ));

        let temp_err = ds_from_csv(&csv)
            .bind_temperature("temperature", Unit::Hour)
            .unwrap_err();
        assert!(matches!(
            temp_err,
            TGADomainError::InvalidUnitForBinding {
                field: "temperature",
                unit: Unit::Hour
            }
        ));

        let mass_err = ds_from_csv(&csv)
            .bind_mass("mass", Unit::Celsius)
            .unwrap_err();
        assert!(matches!(
            mass_err,
            TGADomainError::InvalidUnitForBinding {
                field: "mass",
                unit: Unit::Celsius
            }
        ));

        let generic_err = ds_from_csv(&csv)
            .bind_column(ColumnRole::Time, "time", Unit::Gram)
            .unwrap_err();
        assert!(matches!(
            generic_err,
            TGADomainError::InvalidUnitForBinding {
                field: "time",
                unit: Unit::Gram
            }
        ));
    }

    #[test]
    fn bind_time_accepts_minute_and_roundtrips_units() {
        let csv = make_csv(180, 778);
        let ds = ds_from_csv(&csv).bind_time("time", Unit::Minute).unwrap();
        assert_eq!(ds.schema.columns.get("time").unwrap().unit, Unit::Minute);

        let out = tempfile::NamedTempFile::new().unwrap();
        ds.to_csv_with_units(out.path()).unwrap();

        let content = std::fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("time [min]"));

        let ds2 = TGADataset::from_csv_with_units(out.path()).unwrap();
        assert_eq!(ds2.schema.columns.get("time").unwrap().unit, Unit::Minute);
    }

    #[test]
    fn csv_roundtrip_distinguishes_unknown_and_dimensionless() {
        let csv = make_csv(90, 779);
        let ds = ds_from_csv(&csv)
            .add_numeric_column(
                "aux_unknown",
                Unit::Unknown,
                vec![1.0; 90].into(),
                ColumnNature::Unknown,
            )
            .unwrap()
            .add_numeric_column(
                "aux_dimensionless",
                Unit::Dimensionless,
                vec![2.0; 90].into(),
                ColumnNature::Unknown,
            )
            .unwrap();

        let out = tempfile::NamedTempFile::new().unwrap();
        ds.to_csv_with_units(out.path()).unwrap();

        let content = std::fs::read_to_string(out.path()).unwrap();
        assert!(content.contains("aux_unknown [Unknown]"));
        assert!(content.contains("aux_dimensionless [Dimensionless]"));

        let ds2 = TGADataset::from_csv_with_units(out.path()).unwrap();
        assert_eq!(
            ds2.schema.columns.get("aux_unknown").unwrap().unit,
            Unit::Unknown
        );
        assert_eq!(
            ds2.schema.columns.get("aux_dimensionless").unwrap().unit,
            Unit::Dimensionless
        );
    }

    #[test]
    fn duplicate_column_name_is_rejected_for_creation() {
        let csv = make_csv(60, 782);
        let ds = ds_from_csv(&csv);

        let err = ds
            .clone()
            .add_numeric_column(
                "time",
                Unit::Dimensionless,
                vec![0.0; 60].into(),
                ColumnNature::Unknown,
            )
            .unwrap_err();
        assert!(matches!(err, PolarsError::ComputeError(_)));

        let err = ds
            .clone()
            .with_column_expr(
                ColumnMeta::new(
                    "mass",
                    Unit::Milligram,
                    ColumnOrigin::PolarsDerived,
                    ColumnNature::Mass,
                    ColumnProvenance::manual(
                        "mass",
                        "with_column_expr",
                        Some("attempted rename collision".to_string()),
                    ),
                ),
                col("mass") + lit(1.0),
            )
            .unwrap_err();
        assert!(matches!(
            err,
            TGADomainError::ColumnAlreadyExists(name) if name == "mass"
        ));
    }

    #[test]
    fn rename_column_rejects_existing_target_name() {
        let csv = make_csv(60, 783);
        let ds = ds_from_csv(&csv);

        let err = ds.rename_column("mass", "time").unwrap_err();
        assert!(matches!(
            err,
            TGADomainError::ColumnAlreadyExists(name) if name == "time"
        ));
    }

    #[test]
    fn imported_dataset_writes_history_feed() {
        let csv = make_csv(40, 780);
        let ds = ds_from_csv(&csv);
        let feed = ds.history_of_operations.feed_text();

        assert!(feed.contains("from_csv"));
        assert!(feed.contains("Imported CSV"));
        assert!(!feed.trim().is_empty());
    }

    #[test]
    fn undo_last_restores_dropped_column_and_logs_reverse() {
        let csv = make_csv(60, 781);
        let ds = ds_from_csv(&csv);
        let mut ds = ds.drop_column("temperature").unwrap();

        assert!(!ds.schema.columns.contains_key("temperature"));
        assert!(ds.can_undo());

        ds.undo_last().unwrap();

        assert!(ds.schema.columns.contains_key("temperature"));
        let feed = ds.history_of_operations.feed_text();
        assert!(feed.contains("drop_column"));
        assert!(feed.contains("undo_last"));
    }

    #[test]
    fn calibrate_mass_creates_new_column_and_updates_schema() {
        let csv = make_csv(200, 444);
        let ds = ds_from_csv(&csv);
        let ds_cal = ds.calibrate_mass(2.5, 1.0, "calibrated_mass").unwrap();
        assert_eq!(ds_cal.schema.mass.unwrap(), "calibrated_mass");
        assert!(ds_cal.schema.columns.contains_key("calibrated_mass"));
        assert_eq!(
            ds_cal.schema.columns.get("calibrated_mass").unwrap().unit,
            Unit::Milligram
        );
    }

    #[test]
    fn unary_column_op_applies_function() {
        let csv = make_csv(50, 555);
        let ds = ds_from_csv(&csv);
        let op = UnaryOp {
            func: Box::new(|x| x * 2.0),
            output_unit: Unit::Milligram,
            domain_check: None,
        };
        let ds_op = ds.unary_column_op("mass", "double_mass", op).unwrap();
        assert!(ds_op.schema.columns.contains_key("double_mass"));

        let df = ds_op.frame.collect().unwrap();
        let orig: Vec<f64> = df
            .column("mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let doubled: Vec<f64> = df
            .column("double_mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (a, b) in orig.iter().zip(doubled.iter()) {
            assert!((b - 2.0 * a).abs() < 1e-12);
        }
    }

    #[test]
    fn dimensionless_mass_and_conversion_functions() {
        let csv = make_csv(500, 666);
        let ds = ds_from_csv(&csv);
        let ds_dim = ds.dimensionless_mass(0.0, 2.0, "dim_mass").unwrap();
        assert!(ds_dim.schema.columns.contains_key("dim_mass"));
        assert_eq!(
            ds_dim.schema.columns.get("dim_mass").unwrap().unit,
            Unit::Dimensionless
        );

        let ds_conv = ds_dim.conversion(0.0, 2.0, "conv").unwrap();
        assert!(ds_conv.schema.columns.contains_key("conv"));

        let df = ds_conv.frame.collect().unwrap();
        let dim: Vec<f64> = df
            .column("dim_mass")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let conv: Vec<f64> = df
            .column("conv")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        for (d, c) in dim.iter().zip(conv.iter()) {
            assert!((c - (1.0 - d)).abs() < 1e-12);
        }
    }

    #[test]
    fn explicit_source_column_is_used_for_dimensionless_mass_and_conversion() {
        let csv = make_csv(500, 667);
        let ds = ds_from_csv(&csv);

        let mut shifted_meta = ds.schema.columns.get("mass").unwrap().clone();
        shifted_meta.name = "mass_shifted".to_string();
        shifted_meta.origin = ColumnOrigin::PolarsDerived;
        shifted_meta.nature = ColumnNature::Mass;

        let ds = ds
            .with_column_expr(shifted_meta, col("mass") + lit(1.0))
            .unwrap();

        let alpha_raw = ds
            .clone()
            .dimensionless_mass(0.0, 5.0, "alpha_raw")
            .unwrap();
        let alpha_shifted = ds
            .clone()
            .dimensionless_mass_from_column("mass_shifted", 0.0, 5.0, "alpha_shifted")
            .unwrap();
        let eta_raw = ds.clone().conversion(0.0, 5.0, "eta_raw").unwrap();
        let eta_shifted = ds
            .clone()
            .conversion_from_column("mass_shifted", 0.0, 5.0, "eta_shifted")
            .unwrap();

        let raw_df = alpha_raw.frame.collect().unwrap();
        let shifted_df = alpha_shifted.frame.collect().unwrap();
        let raw_alpha = raw_df.column("alpha_raw").unwrap().f64().unwrap();
        let shifted_alpha = shifted_df.column("alpha_shifted").unwrap().f64().unwrap();
        let sample = 30;
        assert!(
            (raw_alpha.get(sample).unwrap() - shifted_alpha.get(sample).unwrap()).abs() > 1e-6,
            "explicit mass source should change alpha values"
        );

        let raw_eta_df = eta_raw.frame.collect().unwrap();
        let shifted_eta_df = eta_shifted.frame.collect().unwrap();
        let raw_eta = raw_eta_df.column("eta_raw").unwrap().f64().unwrap();
        let shifted_eta = shifted_eta_df.column("eta_shifted").unwrap().f64().unwrap();
        assert!(
            (raw_eta.get(sample).unwrap() - shifted_eta.get(sample).unwrap()).abs() > 1e-6,
            "explicit mass source should change conversion values"
        );
    }

    #[test]
    fn explicit_source_column_is_used_for_mass_and_temperature_rates() {
        let csv = make_csv(500, 668);
        let ds = ds_from_csv(&csv);

        let mut mass_meta = ds.schema.columns.get("mass").unwrap().clone();
        mass_meta.name = "mass_scaled".to_string();
        mass_meta.origin = ColumnOrigin::PolarsDerived;
        mass_meta.nature = ColumnNature::Mass;

        let mut temp_meta = ds.schema.columns.get("temperature").unwrap().clone();
        temp_meta.name = "temperature_scaled".to_string();
        temp_meta.origin = ColumnOrigin::PolarsDerived;
        temp_meta.nature = ColumnNature::Temperature;

        let ds = ds
            .with_column_expr(mass_meta, col("mass") * lit(2.0))
            .unwrap()
            .with_column_expr(temp_meta, col("temperature") * lit(2.0))
            .unwrap();

        let mass_rate_raw = ds.clone().derive_mass_rate("dm_dt_raw").unwrap();
        let mass_rate_scaled = ds
            .clone()
            .derive_mass_rate_from_column("mass_scaled", "dm_dt_scaled")
            .unwrap();
        let temp_rate_raw = ds.clone().derive_temperature_rate("dT_dt_raw").unwrap();
        let temp_rate_scaled = ds
            .clone()
            .derive_temperature_rate_from_column("temperature_scaled", "dT_dt_scaled")
            .unwrap();

        let raw_mass_df = mass_rate_raw.frame.collect().unwrap();
        let scaled_mass_df = mass_rate_scaled.frame.collect().unwrap();
        let raw_mass = raw_mass_df.column("dm_dt_raw").unwrap().f64().unwrap();
        let scaled_mass = scaled_mass_df
            .column("dm_dt_scaled")
            .unwrap()
            .f64()
            .unwrap();

        let sample = 30;
        assert!(
            (scaled_mass.get(sample).unwrap() - raw_mass.get(sample).unwrap() * 2.0).abs() < 1e-6,
            "mass rate should follow the explicit source column"
        );

        let raw_temp_df = temp_rate_raw.frame.collect().unwrap();
        let scaled_temp_df = temp_rate_scaled.frame.collect().unwrap();
        let raw_temp = raw_temp_df.column("dT_dt_raw").unwrap().f64().unwrap();
        let scaled_temp = scaled_temp_df
            .column("dT_dt_scaled")
            .unwrap()
            .f64()
            .unwrap();
        assert!(
            (scaled_temp.get(sample).unwrap() - raw_temp.get(sample).unwrap() * 2.0).abs() < 1e-6,
            "temperature rate should follow the explicit source column"
        );
    }

    #[test]
    fn move_time_to_zero_shifts_time_start() {
        let csv = make_csv(123, 2024);
        let ds = ds_from_csv(&csv);
        let df_before = ds.frame.clone().collect().unwrap();
        let t_before: Vec<f64> = df_before
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        let t0 = t_before[0];
        let ds_zeroed = ds.move_time_to_zero().unwrap();
        let df_after = ds_zeroed.frame.collect().unwrap();
        let t_after: Vec<f64> = df_after
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();
        assert_eq!(t_after.len(), t_before.len());
        // first value becomes zero
        assert!(t_after[0].abs() < 1e-12);
        // all values are shifted by -t0
        for (a, b) in t_before.iter().zip(t_after.iter()) {
            assert!((*b - (a - t0)).abs() < 1e-12);
        }
    }

    #[test]
    fn check_nulls_prints_null_counts() {
        let csv = make_csv(100, 999);
        let ds = ds_from_csv(&csv);

        // Add a column with some nulls by creating a derived column
        let ds_with_nulls = ds
            .with_column_expr(
                ColumnMeta::new(
                    "test_nulls",
                    Unit::Dimensionless,
                    ColumnOrigin::PolarsDerived,
                    ColumnNature::Unknown,
                    ColumnProvenance::manual(
                        "test_nulls",
                        "with_column_expr",
                        Some("null injection test".to_string()),
                    ),
                ),
                when(col("time").lt(lit(5.0)))
                    .then(lit(1.0))
                    .otherwise(lit(NULL)),
            )
            .unwrap();

        // This should print null counts for all columns
        let _ds_checked = ds_with_nulls.check_nulls();

        // Test passes if no panic occurs
        assert!(true);
    }

    #[test]
    fn drop_column_removes_column_and_clears_binding() {
        let csv = make_csv(50, 2025);
        // Drop a non-role column: create a derived column then drop it
        let ds = ds_from_csv(&csv);
        let ds2 = ds
            .with_column_expr(
                ColumnMeta::new(
                    "aux",
                    Unit::Dimensionless,
                    ColumnOrigin::PolarsDerived,
                    ColumnNature::Unknown,
                    ColumnProvenance::manual(
                        "aux",
                        "with_column_expr",
                        Some("auxiliary derived column".to_string()),
                    ),
                ),
                col("time") * lit(0.0),
            )
            .unwrap();
        let df2 = ds2.frame.clone().collect().unwrap();
        assert!(df2.column("aux").is_ok());
        let ds3 = ds2.drop_column("aux").unwrap();
        let df3 = ds3.frame.clone().collect().unwrap();
        assert!(df3.column("aux").is_err());
        assert!(!ds3.schema.columns.contains_key("aux"));

        // Drop a role column: time
        let ds4 = ds3.drop_column("time").unwrap();
        let df4 = ds4.frame.clone().collect().unwrap();
        assert!(df4.column("time").is_err());
        assert!(ds4.schema.time.is_none());

        // Ensure other role bindings remain
        assert!(ds4.schema.temperature.is_some());
        assert!(ds4.schema.mass.is_some());
    }
}
