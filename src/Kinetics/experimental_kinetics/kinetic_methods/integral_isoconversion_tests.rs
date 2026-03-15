//========================================================================================
// TESTS
//========================================================================================

#[cfg(test)]
pub mod tests {
    use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::integral_isoconversion::*;
    use crate::Kinetics::experimental_kinetics::kinetic_methods::tests::{
        build_view_from_cfg, build_view_from_cfg_exact_m0,
    };
    use crate::Kinetics::experimental_kinetics::kinetic_methods::*;
    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::{
        base_advanced_config_isothermal, base_advanced_config_non_isothermal,
    };

    use ndarray::{Array1, Array2};
    use polars::prelude::NamedFromOwned;
    use rayon::prelude::*;
    use std::time::Instant;
    const R: f64 = 8.314462618;

    /// simulate first-order non-isothermal TGA and produce ConversionGrid
    pub fn simulate_tga_first_order(
        t0: f64,
        k0: f64,
        ea: f64,
        heating_rates: &[f64],
        dt: f64,
        n_steps: usize,
        alpha_grid: Vec<f64>,
    ) -> ConversionGrid {
        let n_exp = heating_rates.len();
        let n_alpha = alpha_grid.len();

        let mut temperature = Array2::<f64>::zeros((n_exp, n_alpha));
        let mut time = Array2::<f64>::zeros((n_exp, n_alpha));
        let mut rate = Array2::<f64>::zeros((n_exp, n_alpha));

        for (exp_i, &beta) in heating_rates.iter().enumerate() {
            // heating_rates are provided in K/min; convert to K/s for time integration
            let beta = beta / 60.0;
            let mut t_vec = Vec::with_capacity(n_steps);
            let mut temp_vec = Vec::with_capacity(n_steps);
            let mut alpha_vec = Vec::with_capacity(n_steps);
            let mut rate_vec = Vec::with_capacity(n_steps);

            let mut alpha = 0.0;

            for step in 0..n_steps {
                let t = step as f64 * dt;

                let temp = t0 + beta * t;

                let k = k0 * (-ea / (R * temp)).exp();

                let dalpha = k * (1.0 - alpha);

                alpha += dalpha * dt;

                if alpha > 0.999 {
                    break;
                }

                t_vec.push(t);
                temp_vec.push(temp);
                alpha_vec.push(alpha);
                rate_vec.push(dalpha);
            }

            for (j, &target_alpha) in alpha_grid.iter().enumerate() {
                let mut i = 0;

                while i + 1 < alpha_vec.len() && alpha_vec[i + 1] < target_alpha {
                    i += 1;
                }

                if i + 1 >= alpha_vec.len() {
                    continue;
                }

                let a0 = alpha_vec[i];
                let a1 = alpha_vec[i + 1];

                let w = (target_alpha - a0) / (a1 - a0);

                let t = t_vec[i] + w * (t_vec[i + 1] - t_vec[i]);
                let temp = temp_vec[i] + w * (temp_vec[i + 1] - temp_vec[i]);
                let r = rate_vec[i] + w * (rate_vec[i + 1] - rate_vec[i]);

                time[[exp_i, j]] = t;
                temperature[[exp_i, j]] = temp;
                rate[[exp_i, j]] = r;
            }
        }

        let inv_temperature = temperature.mapv(|t| 1.0 / t);

        let meta = heating_rates
            .iter()
            .map(|&b| ExperimentMeta {
                heating_rate: Some(b),
                ..ExperimentMeta::default()
            })
            .collect();

        ConversionGrid {
            eta: Array1::from(alpha_grid),
            temperature,
            inv_temperature,
            time,
            conversion_rate: rate,
            dt: None,
            meta,
        }
    }
    /// simulate first-order non-isothermal TGA and produce ConversionGrid
    pub fn simulate_tga_first_order2(
        t0: f64,
        k0: f64,
        ea: f64,
        heating_rates: &[f64],
        dt: f64,
        n_steps: usize,
        alpha_grid: Option<Vec<f64>>,
    ) -> ConversionGrid {
        let n_exp = heating_rates.len();

        #[derive(Clone)]
        struct ExpTraj {
            t: Vec<f64>,
            temp: Vec<f64>,
            alpha: Vec<f64>,
            rate: Vec<f64>,
        }

        #[derive(Clone)]
        struct ExpRow {
            time: Vec<f64>,
            temp: Vec<f64>,
            rate: Vec<f64>,
        }

        fn is_sorted_non_decreasing(values: &[f64]) -> bool {
            values.windows(2).all(|w| w[0] <= w[1])
        }

        // First, integrate each experiment to collect alpha/time/temperature/rate trajectories
        let per_exp: Vec<ExpTraj> = heating_rates
            .par_iter()
            .map(|&beta_min| {
                // convert K/min to K/s
                let beta = beta_min / 60.0;
                let mut t_vec = Vec::with_capacity(n_steps);
                let mut temp_vec = Vec::with_capacity(n_steps);
                let mut alpha_vec = Vec::with_capacity(n_steps);
                let mut rate_vec = Vec::with_capacity(n_steps);

                let mut alpha = 0.0;
                for step in 0..n_steps {
                    let t = step as f64 * dt;
                    let temp = t0 + beta * t;
                    let k = k0 * (-ea / (R * temp)).exp();
                    let dalpha = k * (1.0 - alpha);
                    alpha += dalpha * dt;

                    t_vec.push(t);
                    temp_vec.push(temp);
                    alpha_vec.push(alpha);
                    rate_vec.push(dalpha);

                    if alpha >= 0.999 {
                        break;
                    }
                }

                ExpTraj {
                    t: t_vec,
                    temp: temp_vec,
                    alpha: alpha_vec,
                    rate: rate_vec,
                }
            })
            .collect();

        // Determine alpha targets
        let alpha_values: Vec<f64> = match alpha_grid {
            Some(v) => v,
            None => {
                let n_alpha = n_steps.max(2);
                (0..n_alpha)
                    .map(|i| (i as f64) / (n_alpha as f64 - 1.0))
                    .collect()
            }
        };

        let n_alpha = alpha_values.len();
        let alpha_sorted = is_sorted_non_decreasing(&alpha_values);

        // Fill by interpolation for each experiment (parallel)
        let rows: Vec<ExpRow> = (0..n_exp)
            .into_par_iter()
            .map(|exp_i| {
                let t_vec = &per_exp[exp_i].t;
                let temp_vec = &per_exp[exp_i].temp;
                let alpha_vec = &per_exp[exp_i].alpha;
                let rate_vec = &per_exp[exp_i].rate;

                let mut time_row = vec![0.0f64; n_alpha];
                let mut temp_row = vec![0.0f64; n_alpha];
                let mut rate_row = vec![0.0f64; n_alpha];

                let len = alpha_vec.len();
                if len == 0 {
                    return ExpRow {
                        time: time_row,
                        temp: temp_row,
                        rate: rate_row,
                    };
                }

                if alpha_sorted {
                    let mut i = 0usize;
                    for (j, &target_alpha) in alpha_values.iter().enumerate() {
                        while i + 1 < len && alpha_vec[i + 1] < target_alpha {
                            i += 1;
                        }
                        if i + 1 >= len {
                            let last = len - 1;
                            time_row[j] = t_vec[last];
                            temp_row[j] = temp_vec[last];
                            rate_row[j] = rate_vec[last];
                            continue;
                        }

                        let a0 = alpha_vec[i];
                        let a1 = alpha_vec[i + 1];
                        let denom = a1 - a0;

                        if denom.abs() < 1e-14 || target_alpha <= a0 {
                            time_row[j] = t_vec[i];
                            temp_row[j] = temp_vec[i];
                            rate_row[j] = rate_vec[i];
                            continue;
                        }

                        let w = (target_alpha - a0) / denom;
                        time_row[j] = t_vec[i] + w * (t_vec[i + 1] - t_vec[i]);
                        temp_row[j] = temp_vec[i] + w * (temp_vec[i + 1] - temp_vec[i]);
                        rate_row[j] = rate_vec[i] + w * (rate_vec[i + 1] - rate_vec[i]);
                    }
                } else {
                    use std::cmp::Ordering;
                    for (j, &target_alpha) in alpha_values.iter().enumerate() {
                        let idx = match alpha_vec.binary_search_by(|a| {
                            a.partial_cmp(&target_alpha).unwrap_or(Ordering::Less)
                        }) {
                            Ok(i) => i,
                            Err(i) => i,
                        };
                        if idx == 0 {
                            time_row[j] = t_vec[0];
                            temp_row[j] = temp_vec[0];
                            rate_row[j] = rate_vec[0];
                        } else if idx >= len {
                            let last = len - 1;
                            time_row[j] = t_vec[last];
                            temp_row[j] = temp_vec[last];
                            rate_row[j] = rate_vec[last];
                        } else {
                            let i = idx - 1;
                            let a0 = alpha_vec[i];
                            let a1 = alpha_vec[idx];
                            let denom = a1 - a0;
                            if denom.abs() < 1e-14 {
                                time_row[j] = t_vec[i];
                                temp_row[j] = temp_vec[i];
                                rate_row[j] = rate_vec[i];
                            } else {
                                let w = (target_alpha - a0) / denom;
                                time_row[j] = t_vec[i] + w * (t_vec[idx] - t_vec[i]);
                                temp_row[j] = temp_vec[i] + w * (temp_vec[idx] - temp_vec[i]);
                                rate_row[j] = rate_vec[i] + w * (rate_vec[idx] - rate_vec[i]);
                            }
                        }
                    }
                }

                ExpRow {
                    time: time_row,
                    temp: temp_row,
                    rate: rate_row,
                }
            })
            .collect();

        let mut temperature_data = Vec::with_capacity(n_exp * n_alpha);
        let mut time_data = Vec::with_capacity(n_exp * n_alpha);
        let mut rate_data = Vec::with_capacity(n_exp * n_alpha);

        for row in &rows {
            time_data.extend_from_slice(&row.time);
            temperature_data.extend_from_slice(&row.temp);
            rate_data.extend_from_slice(&row.rate);
        }

        let temperature = Array2::from_shape_vec((n_exp, n_alpha), temperature_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));
        let time = Array2::from_shape_vec((n_exp, n_alpha), time_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));
        let rate = Array2::from_shape_vec((n_exp, n_alpha), rate_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));

        let inv_temperature = temperature.mapv(|t| 1.0 / t);

        let meta = heating_rates
            .iter()
            .map(|&b| ExperimentMeta {
                heating_rate: Some(b),
                ..ExperimentMeta::default()
            })
            .collect();

        ConversionGrid {
            eta: Array1::from(alpha_values),
            temperature,
            inv_temperature,
            time,
            conversion_rate: rate,
            dt: None,
            meta,
        }
    }

    /// Simulate first-order isothermal TGA at several constant temperatures and produce ConversionGrid.
    /// Each experiment runs at a fixed temperature T[i]; the "experiments" axis corresponds to
    /// different isothermal temperatures rather than different heating rates.
    pub fn simulate_tga_first_order_isothermal(
        t0: f64,
        k0: f64,
        ea: f64,
        T: &[f64],
        dt: f64,
        n_steps: usize,
        alpha_grid: Option<Vec<f64>>,
    ) -> ConversionGrid {
        let n_exp = T.len();

        #[derive(Clone)]
        struct ExpTraj {
            t: Vec<f64>,
            alpha: Vec<f64>,
            rate: Vec<f64>,
        }

        let per_exp: Vec<ExpTraj> = T
            .par_iter()
            .map(|&temp| {
                let k = k0 * (-ea / (R * temp)).exp();
                let mut t_vec = Vec::with_capacity(n_steps);
                let mut alpha_vec = Vec::with_capacity(n_steps);
                let mut rate_vec = Vec::with_capacity(n_steps);
                let mut alpha = 0.0f64;
                for step in 0..n_steps {
                    let t = step as f64 * dt;
                    let dalpha = k * (1.0 - alpha);
                    alpha += dalpha * dt;
                    t_vec.push(t);
                    alpha_vec.push(alpha);
                    rate_vec.push(dalpha);
                    if alpha >= 0.999 {
                        break;
                    }
                }
                ExpTraj {
                    t: t_vec,
                    alpha: alpha_vec,
                    rate: rate_vec,
                }
            })
            .collect();

        fn is_sorted_non_decreasing(v: &[f64]) -> bool {
            v.windows(2).all(|w| w[0] <= w[1])
        }

        let alpha_values: Vec<f64> = match alpha_grid {
            Some(v) => v,
            None => {
                let n = n_steps.max(2);
                (0..n).map(|i| i as f64 / (n as f64 - 1.0)).collect()
            }
        };
        let n_alpha = alpha_values.len();
        let alpha_sorted = is_sorted_non_decreasing(&alpha_values);

        let mut temperature_data = vec![0.0f64; n_exp * n_alpha];
        let mut time_data = vec![0.0f64; n_exp * n_alpha];
        let mut rate_data = vec![0.0f64; n_exp * n_alpha];

        for (exp_i, traj) in per_exp.iter().enumerate() {
            let temp = T[exp_i];
            let len = traj.alpha.len();
            for (j, &target_alpha) in alpha_values.iter().enumerate() {
                let idx = if alpha_sorted {
                    let mut i = 0usize;
                    while i + 1 < len && traj.alpha[i + 1] < target_alpha {
                        i += 1;
                    }
                    i
                } else {
                    use std::cmp::Ordering;
                    match traj.alpha.binary_search_by(|a| {
                        a.partial_cmp(&target_alpha).unwrap_or(Ordering::Less)
                    }) {
                        Ok(i) | Err(i) => i.saturating_sub(1),
                    }
                };

                let base = exp_i * n_alpha + j;
                temperature_data[base] = temp;

                if idx + 1 >= len {
                    time_data[base] = traj.t[len - 1];
                    rate_data[base] = traj.rate[len - 1];
                } else {
                    let a0 = traj.alpha[idx];
                    let a1 = traj.alpha[idx + 1];
                    let denom = a1 - a0;
                    if denom.abs() < 1e-14 || target_alpha <= a0 {
                        time_data[base] = traj.t[idx];
                        rate_data[base] = traj.rate[idx];
                    } else {
                        let w = (target_alpha - a0) / denom;
                        time_data[base] = traj.t[idx] + w * (traj.t[idx + 1] - traj.t[idx]);
                        rate_data[base] =
                            traj.rate[idx] + w * (traj.rate[idx + 1] - traj.rate[idx]);
                    }
                }
            }
        }

        let temperature = Array2::from_shape_vec((n_exp, n_alpha), temperature_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));
        let time = Array2::from_shape_vec((n_exp, n_alpha), time_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));
        let rate = Array2::from_shape_vec((n_exp, n_alpha), rate_data)
            .unwrap_or_else(|_| Array2::<f64>::zeros((n_exp, n_alpha)));
        let inv_temperature = temperature.mapv(|t| 1.0 / t);

        let meta = T
            .iter()
            .map(|&t_iso| ExperimentMeta {
                isothermal_temperature: Some(t_iso),
                ..ExperimentMeta::default()
            })
            .collect();

        ConversionGrid {
            eta: Array1::from(alpha_values),
            temperature,
            inv_temperature,
            time,
            conversion_rate: rate,
            dt: None,
            meta,
        }
    }

    pub fn simulate_tga_first_order_with_dt(
        t0: f64,
        k0: f64,
        ea: f64,
        heating_rates: &[f64],
        dt: f64,
        n_steps: usize,
        alpha_grid: Option<Vec<f64>>,
    ) -> ConversionGrid {
        let mut grid =
            simulate_tga_first_order2(t0, k0, ea, heating_rates, dt, n_steps, alpha_grid);

        let time = grid.time.clone();
        let dt = ConversionGridBuilder::compute_dt_matrix(&time);
        grid.dt = Some(dt);
        grid
    }

    #[test]
    fn solver_basic_shapes_are_consistent() {
        let eta = Array1::from(vec![0.2, 0.5, 0.8]);
        let temperature =
            Array2::from_shape_vec((2, 3), vec![350.0, 360.0, 370.0, 400.0, 410.0, 420.0]).unwrap();
        let inv_temperature = temperature.mapv(|t| 1.0 / t);
        let time = Array2::from_shape_vec((2, 3), vec![0.0, 1.0, 2.0, 0.0, 1.5, 3.0]).unwrap();
        let conversion_rate =
            Array2::from_shape_vec((2, 3), vec![0.05, 0.06, 0.07, 0.08, 0.09, 0.1]).unwrap();

        let meta = vec![
            ExperimentMeta {
                heating_rate: Some(5.0),
                ..ExperimentMeta::default()
            },
            ExperimentMeta {
                heating_rate: Some(10.0),
                ..ExperimentMeta::default()
            },
        ];

        let grid = ConversionGrid {
            eta,
            temperature,
            inv_temperature,
            time,
            conversion_rate,
            dt: None,
            meta,
        };

        let solver = IntegralIsoconversionalSolver::ofw();
        let result = solver.solve(&grid).unwrap();

        assert_eq!(result.layers.len(), 3);
        for layer in &result.layers {
            assert!(layer.ea.is_finite());
            assert!(layer.regression.r2.is_finite());
        }
    }
    //============================================================================================
    //         TESTS WITH FULL PIPELINE ON MOCK DATA
    //============================================================================================

    #[test]
    fn ofw_compute_with_mock_non_isothermal_data() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();

        let result = OFW {}.compute(&view).unwrap();

        assert_eq!(result.layers.len(), 50); // default grid segments in ConversionGridBuilder
        for layer in &result.layers {
            let eta = layer.eta;
            if eta > 0.05 && eta < 0.95 {
                assert!(layer.ea.is_finite());
                assert!((0.0..=1.0).contains(&layer.regression.r2));
            }
        }
    }

    #[test]

    fn ofw_compute_with_mock_non_isothermal_data2() {
        let now = Instant::now();
        let cfg = base_advanced_config_non_isothermal(
            420.0,
            1e6,
            100_000.0,
            1.0, // in seconds
            30_000,
            vec![0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], //K/min
        );
        let now0 = Instant::now();
        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();
        println!(
            "kinetic data view created in {}",
            now0.elapsed().as_millis()
        );
        for exp in &view.experiments {
            let heating_rate = exp.meta.heating_rate;
            let conversion: Vec<f64> = exp
                .conversion
                .clone()
                .iter()
                .map(|x| *x)
                .step_by(1000)
                .collect();
            let temp: Vec<f64> = exp
                .temperature
                .clone()
                .iter()
                .map(|x| *x)
                .step_by(1000)
                .collect();
            println!(
                "\n =================heating rate {:?}============================= \n coversion {:?}, \n T ={:?} ",
                heating_rate, conversion, temp
            );
        }

        let result = OFW {}.compute(&view).unwrap();
        println!("elapsed: {}", now.elapsed().as_millis());
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }
    #[test]

    fn ofw_compute_panics_without_heating_rates() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, vec![600.0, 700.0], 0.1, 10_000);
        let view = build_view_from_cfg(&cfg).unwrap();

        let _ = OFW {}.compute(&view);
    }
    //=======================================================================================================
    //  TESTS WITH SIMULATED GRIDS
    //====================================================================================================
    /*
    #[test]
    fn mock_grid_kas_recovers_activation_energy() {
        let e_true = 80_000.0;
        let grid = simulate_tga_first_order(
            500.0,
            1e5,
            e_true,
            &[2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            0.5,
            100_00,
            (5..90).map(|x| x as f64 / 100.0).collect(),
        );
        grid.report();
        let solver = IntegralIsoconversionalSolver::kas();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn mock_grid_starink_recovers_activation_energy() {
        let e_true = 120_000.0;
        let grid = simulate_tga_first_order(
            600.0,
            1e8,
            e_true,
            &[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0],
            0.3,
            150_00,
            (5..90).map(|x| x as f64 / 100.0).collect(),
        );
        grid.report();
        let solver = IntegralIsoconversionalSolver::starink();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn mock_grid_ofw_basic_test() {
        let grid = simulate_tga_first_order(
            550.0,
            5e6,
            95_000.0,
            &[0.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            1.0,
            80_00,
            (5..90).map(|x| x as f64 / 100.0).collect(),
        );
        grid.report();
        let solver = IntegralIsoconversionalSolver::ofw();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn test_kas_solver() {
        let grid = simulate_tga_first_order(
            300.0,
            1e13,
            120_000.0,
            &[0.5, 1.0, 2.0, 3.0, 4.0],
            0.1,
            20_00,
            (5..90).map(|x| x as f64 / 100.0).collect(),
        );
        grid.report();
        let solver = IntegralIsoconversionalSolver::kas();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }
    */
    //=========================================================================================
    #[test]
    fn mock_grid_kas_recovers_activation_energ2() {
        let e_true = 85_000.0;
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            e_true,
            &[2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            0.5,
            500_00,
            None,
        );
        let solver = IntegralIsoconversionalSolver::kas();
        grid.report();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn mock_grid_starink_recovers_activation_energy2() {
        let e_true = 150_000.0;
        let grid = simulate_tga_first_order2(
            600.0,
            1e7,
            e_true,
            &[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0],
            1.0,
            175000,
            None,
        );
        let solver = IntegralIsoconversionalSolver::starink();
        grid.report();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn mock_grid_ofw_basic_test2() {
        let grid = simulate_tga_first_order2(
            550.0,
            5e7,
            150_000.0,
            &[0.5, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            1.0,
            8000_0,
            None,
        );
        let T = grid.temperature.row(0).clone();
        println!("T={:?}, {}", T, T.len());
        let solver = IntegralIsoconversionalSolver::ofw();
        grid.report();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.05, 0.95, 100, Some(0.99));
    }

    #[test]
    fn test_kas_solver2() {
        let now = Instant::now();
        let grid = simulate_tga_first_order2(
            400.0,
            1e12,
            130_000.0,
            &[0.5, 1.0, 2.0, 3.0, 4.0],
            0.1,
            200_000,
            None,
        );
        println!(
            "Synthetic Grid generated at {} milliseconds",
            now.elapsed().as_millis()
        );
        let now = Instant::now();
        let solver = IntegralIsoconversionalSolver::kas();
        grid.report();
        let result = solver.solve(&grid).unwrap();
        result.pretty_print_and_assert(0.10, 0.95, 100, Some(0.99));
        println!(
            "problem solved at {} milliseconds",
            now.elapsed().as_millis()
        );

        println!("test_kas_solver2 passed")
    }

    //============================================================================================
    // SIMULATOR SANITY CHECKS
    //=============================================================================================
    #[test]
    fn compare_grid_with_simulation_first_order() {
        let heating_rates = vec![0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5];
        let cfg = base_advanced_config_non_isothermal(
            520.0,
            1e5,
            76_500.0,
            0.1,
            10_000,
            heating_rates.clone(),
        );

        let view = build_view_from_cfg_exact_m0(&cfg).unwrap();

        let grid = ConversionGridBuilder::new()
            //.eta_range(0.05, 0.95)
            .auto_range()
            .segments(50)
            .interpolation(GridInterpolation::Linear)
            .build(&view)
            .unwrap();

        let sim_grid = simulate_tga_first_order2(
            520.0,
            1e5,
            76_500.0,
            &heating_rates,
            cfg.dt,
            cfg.n_points,
            Some(grid.eta.to_vec()),
        );

        let mut max_abs_diff = 0.0f64;
        let mut max_rel_diff = 0.0f64;
        for ((a, b), _idx) in grid
            .temperature
            .iter()
            .zip(sim_grid.temperature.iter())
            .zip(0..)
        {
            let diff = (a - b).abs();
            max_abs_diff = max_abs_diff.max(diff);
            if b.abs() > 1e-9 {
                max_rel_diff = max_rel_diff.max(diff / b.abs());
            }
        }

        println!(
            "max_abs_temp_diff={}, max_rel_temp_diff={}",
            max_abs_diff, max_rel_diff
        );

        assert!(
            max_abs_diff < 1.0,
            "Temperature grid deviates too much from simulation"
        );
    }

    // Basic tests showing simulate_tga_first_order works with alpha_grid = None and Some
    #[test]
    fn simulate_uses_n_steps_when_no_alpha_grid() {
        let heating_rates = [1.0, 2.0, 3.0];
        let n_steps = 101;
        let dt = 0.1;
        let grid =
            simulate_tga_first_order2(500.0, 1e5, 80_000.0, &heating_rates, dt, n_steps, None);
        grid.report();
        //   assert_eq!(grid.temperature.len(), n_steps);
        let temp = grid.temperature.row(0);
        println!("temp={:?}", temp);
        assert_eq!(grid.eta.len(), n_steps);
        // Check equidistribution: first ~0.0, last ~1.0, constant step
        let tol = 1e-12;
        let first = grid.eta[0];
        let last = grid.eta[n_steps - 1];
        assert!((first - 0.0).abs() < tol);
        assert!((last - 1.0).abs() < tol);
        let expected_step = 1.0 / (n_steps as f64 - 1.0);
        for w in grid.eta.windows(2) {
            let step = w[1] - w[0];
            assert!((step - expected_step).abs() < 1e-9);
        }
    }

    #[test]
    fn simulate_respects_provided_alpha_grid() {
        let heating_rates = [1.0, 2.0, 3.0];
        // Explicit alpha grid (0.05..0.89 step 0.01 -> 85 points)
        let alpha_grid: Vec<f64> = (5..90).map(|x| x as f64 / 100.0).collect();
        let n_steps = 10_000; // should not affect eta when alpha_grid is provided
        let grid = simulate_tga_first_order2(
            500.0,
            1e5,
            80_000.0,
            &heating_rates,
            0.1,
            n_steps,
            Some(alpha_grid.clone()),
        );
        assert_eq!(grid.eta.len(), alpha_grid.len());
        for (a, b) in grid.eta.iter().zip(alpha_grid.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn simulate_tga_first_order2_temperature_monotone_with_alpha() {
        let heating_rates = [1.0];
        let alpha_grid: Vec<f64> = (0..=100).map(|x| x as f64 / 100.0).collect();
        let grid = simulate_tga_first_order2(
            450.0,
            1e5,
            90_000.0,
            &heating_rates,
            0.2,
            50_000,
            Some(alpha_grid),
        );
        let temp = grid.temperature.row(0);
        for w in temp.windows(2) {
            assert!(w[1] >= w[0]);
        }
    }

    #[test]
    fn simulate_tga_first_order2_handles_unsorted_alpha_grid() {
        let heating_rates = [2.0];
        let alpha_grid = vec![0.2, 0.1, 0.9, 0.5, 0.8];
        let grid = simulate_tga_first_order2(
            500.0,
            1e6,
            110_000.0,
            &heating_rates,
            0.3,
            40_000,
            Some(alpha_grid.clone()),
        );
        assert_eq!(grid.eta.to_vec(), alpha_grid);
        let temp = grid.temperature.row(0);
        let time = grid.time.row(0);
        for &v in temp.iter() {
            assert!(v.is_finite());
        }
        for &v in time.iter() {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn simulate_tga_first_order2_no_alpha_grid_time_is_non_decreasing() {
        let heating_rates = [0.5, 1.0];
        let grid =
            simulate_tga_first_order2(520.0, 1e5, 85_000.0, &heating_rates, 0.2, 5_000, None);
        for exp_i in 0..heating_rates.len() {
            let time = grid.time.row(exp_i);
            for w in time.windows(2) {
                assert!(w[1] >= w[0]);
            }
        }
    }
}
