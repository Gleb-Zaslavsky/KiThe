pub mod Friedman;
pub mod Vyazovkin;
pub mod integral_isoconversion;
mod integral_isoconversion_tests;
pub mod isoconversion;
pub mod kinetic_regression;
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::Kinetics::experimental_kinetics::experiment_series2::UnitedDataset;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, TGADomainError,
};
use ndarray::{Array1, Array2};
use splines::{Interpolation, Key, Spline};

pub struct KineticDataView {
    pub experiments: Vec<ExperimentData>,
}
impl KineticDataView {
    pub fn from_united_dataset(ds: &UnitedDataset) -> Result<Self, TGADomainError> {
        let mut experiments = Vec::new();

        for meta in &ds.meta {
            let time = ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Time)?;

            let temperature =
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Temperature)?;

            let conversion =
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Conversion)?;

            let conversion_rate =
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::ConversionRate)?;

            experiments.push(ExperimentData {
                meta: meta.clone(),
                time,
                temperature,
                conversion,
                conversion_rate,
                mass: None,
                mass_rate: None,
            });
        }

        Ok(Self { experiments })
    }
}
/// materialized data of one experiment
pub struct ExperimentData {
    pub meta: ExperimentMeta,

    pub time: Vec<f64>,
    pub temperature: Vec<f64>,
    // pub temperature_rate: Option<Vec<f64>>,
    pub conversion: Vec<f64>,
    pub conversion_rate: Vec<f64>,

    pub mass: Option<Vec<f64>>,
    pub mass_rate: Option<Vec<f64>>,
    // pub relative_mass:  Option<Vec<f64>>,
    //  pub relative_mass_rate:  Option<Vec<f64>>,
}
/// common trait for kinetic methods
pub trait KineticMethod {
    type Output;

    fn name(&self) -> &'static str;

    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError>;
    // method-specific checks
    fn check_input(&self, _data: &KineticDataView) -> Result<(), TGADomainError> {
        Ok(())
    }

    fn requirements(&self) -> KineticRequirements;
}

impl<O> dyn KineticMethod<Output = O> {
    pub fn run<M: KineticMethod<Output = O>>(
        method: &M,
        data: &KineticDataView,
    ) -> Result<O, TGADomainError> {
        method.check_input(data)?;

        method.compute(data)
    }
}

pub fn require_min_experiments(data: &KineticDataView, n: usize) -> Result<(), TGADomainError> {
    if data.experiments.len() < n {
        return Err(TGADomainError::InvalidOperation(format!(
            "Method requires at least {} experiments",
            n
        )));
    }

    Ok(())
}

pub fn require_heating_rates(data: &KineticDataView) -> Result<(), TGADomainError> {
    for exp in &data.experiments {
        if exp.meta.heating_rate.is_none() {
            return Err(TGADomainError::InvalidOperation(
                "Heating rate missing".into(),
            ));
        }
    }

    Ok(())
}

pub fn require_isothermal(data: &KineticDataView) -> Result<(), TGADomainError> {
    for exp in &data.experiments {
        if exp.meta.isothermal_temperature.is_none() {
            return Err(TGADomainError::InvalidOperation(
                "Method requires isothermal experiments".into(),
            ));
        }
    }

    Ok(())
}

pub fn require_conversion_rate(data: &KineticDataView) -> Result<(), TGADomainError> {
    for exp in &data.experiments {
        if exp.conversion_rate.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Conversion rate missing".into(),
            ));
        }
    }

    Ok(())
}

pub fn require_conversion(data: &KineticDataView) -> Result<(), TGADomainError> {
    for exp in &data.experiments {
        if exp.conversion.is_empty() {
            return Err(TGADomainError::InvalidOperation(
                "Conversion  missing".into(),
            ));
        }
    }

    Ok(())
}

#[derive(Default)]
pub struct KineticRequirements {
    pub min_experiments: usize,
    needs_conversion: bool,

    needs_conversion_rate: bool,

    needs_temperature: bool,

    needs_heating_rate: bool,
}
// declarative checks
pub fn check_requirements(
    data: &KineticDataView,
    req: &KineticRequirements,
) -> Result<(), TGADomainError> {
    require_min_experiments(data, req.min_experiments)?;

    if req.needs_heating_rate {
        require_heating_rates(data)?;
    }

    if req.needs_temperature {
        require_isothermal(data)?;
    }

    if req.needs_conversion_rate {
        require_conversion_rate(data)?;
    }
    if req.needs_conversion {
        require_conversion(data)?
    }
    Ok(())
}

pub fn run_method_with_requirements<M: KineticMethod>(
    method: &M,
    data: &KineticDataView,
) -> Result<M::Output, TGADomainError> {
    let req = method.requirements();

    check_requirements(data, &req)?;

    method.check_input(data)?;

    method.compute(data)
}
//===============================================================================================================
// CONVERSION GRID
//===========================================================================================================
#[derive(Clone, Copy, Debug)]
pub enum GridInterpolation {
    Linear,
    Spline,
}

pub struct ConversionGrid {
    /// conversion grid
    pub eta: Array1<f64>,

    /// temperature[experiment, eta]
    pub temperature: Array2<f64>,
    /// 1/T cached
    pub inv_temperature: Array2<f64>,
    /// time[experiment, eta]
    pub time: Array2<f64>,
    /// rate eta/dt[experiment, eta]
    pub conversion_rate: Array2<f64>,
    pub dt: Option<Array2<f64>>,
    pub meta: Vec<ExperimentMeta>,
}
impl ConversionGrid {
    pub fn report(&self) {
        use tabled::{Table, Tabled, settings::Style};

        println!("\n=== CONVERSION GRID REPORT ===");
        println!("\n1) Eta (conversion) range:");
        println!(
            "   Min: {:.6}",
            self.eta.iter().cloned().fold(f64::INFINITY, f64::min)
        );
        println!(
            "   Max: {:.6}",
            self.eta.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        );
        println!("   Grid points: {}", self.eta.len());

        println!("\n2) Maximum eta values for each experiment:");

        #[derive(Tabled)]
        struct MaxEtaRow {
            #[tabled(rename = "Experiment")]
            exp_id: String,
            #[tabled(rename = "Time (s)")]
            time: String,
            #[tabled(rename = "Temperature (K)")]
            temperature: String,
        }

        let mut rows = Vec::new();
        for (i, meta) in self.meta.iter().enumerate() {
            let max_eta_idx = self.eta.len() - 1;
            let time_val = self.time[[i, max_eta_idx]];
            let temp_val = self.temperature[[i, max_eta_idx]];

            rows.push(MaxEtaRow {
                exp_id: meta.id.clone(),
                time: format!("{:.4}", time_val),
                temperature: format!("{:.2}", temp_val),
            });
        }

        let table = Table::new(rows).with(Style::rounded()).to_string();
        println!("{}", table);

        println!("\n3) Zero detection:");
        let temp_zeros = self.temperature.iter().filter(|&&v| v == 0.0).count();
        let time_zeros = self.time.iter().filter(|&&v| v == 0.0).count();

        if temp_zeros > 0 {
            println!(
                "   WARNING: {} zero values in temperature array",
                temp_zeros
            );
            for i in 0..self.temperature.nrows() {
                let row_zeros = self
                    .temperature
                    .row(i)
                    .iter()
                    .filter(|&&v| v == 0.0)
                    .count();
                if row_zeros > 0 {
                    println!(
                        "     - Row {} (exp {}): {} zeros",
                        i, self.meta[i].id, row_zeros
                    );
                }
            }
        } else {
            println!("   Temperature: No zeros detected");
        }

        if time_zeros > 0 {
            println!("   WARNING: {} zero values in time array", time_zeros);
            for i in 0..self.time.nrows() {
                let row_zeros = self.time.row(i).iter().filter(|&&v| v == 0.0).count();
                if row_zeros > 0 {
                    println!(
                        "     - Row {} (exp {}): {} zeros",
                        i, self.meta[i].id, row_zeros
                    );
                }
            }
        } else {
            println!("   Time: No zeros detected");
        }

        println!("\n4) Eta values close to 1.0:");
        let threshold = 0.99;
        let close_to_one = self.eta.iter().filter(|&&v| v >= threshold).count();
        let fraction = close_to_one as f64 / self.eta.len() as f64;
        println!(
            "   Elements >= {}: {} / {} ({:.2}%)",
            threshold,
            close_to_one,
            self.eta.len(),
            fraction * 100.0
        );

        println!("\n5) Additional information:");
        println!("   Number of experiments: {}", self.meta.len());
        println!(
            "   Grid dimensions: {} x {}",
            self.temperature.nrows(),
            self.temperature.ncols()
        );
        println!("   dt matrix computed: {}", self.dt.is_some());

        let temp_min = self
            .temperature
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let temp_max = self
            .temperature
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        println!("   Temperature range: {:.2} - {:.2} K", temp_min, temp_max);

        let time_min = self.time.iter().cloned().fold(f64::INFINITY, f64::min);
        let time_max = self.time.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("   Time range: {:.4} - {:.4} s", time_min, time_max);

        let rate_min = self
            .conversion_rate
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let rate_max = self
            .conversion_rate
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "   Conversion rate range: {:.6} - {:.6} s⁻¹",
            rate_min, rate_max
        );

        println!("\n=== END OF REPORT ===\n");
    }
}

// syntax
// let grid = ConversionGridBuilder::new()
//     .eta_range(0.05, 0.9)
//     .segments(200)
//     .interpolation(GridInterpolation::Linear)
//     .build(&view)?;

pub enum AlphaRangeMode {
    Auto,

    Manual { min: f64, max: f64 },
}
// delta t required for Vyazovkin
pub struct GridExtras {
    pub dt_matrix: bool,
}
impl Default for GridExtras {
    fn default() -> Self {
        Self { dt_matrix: false }
    }
}
pub struct ConversionGridBuilder {
    eta_mode: AlphaRangeMode,
    segments: usize,
    interpolation: GridInterpolation,
    /// safety margins
    margin_left: f64,
    margin_right: f64,
    extras: GridExtras,
}

impl Default for ConversionGridBuilder {
    fn default() -> Self {
        Self {
            eta_mode: AlphaRangeMode::Auto,

            segments: 50,

            interpolation: GridInterpolation::Linear,
            margin_left: 0.0,
            margin_right: 0.0,
            extras: GridExtras::default(),
        }
    }
}

impl ConversionGridBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn eta_range(mut self, min: f64, max: f64) -> Self {
        self.eta_mode = AlphaRangeMode::Manual { min, max };
        self
    }

    pub fn auto_range(mut self) -> Self {
        self.eta_mode = AlphaRangeMode::Auto;
        self
    }

    pub fn segments(mut self, n: usize) -> Self {
        self.segments = n;
        self
    }

    pub fn interpolation(mut self, mode: GridInterpolation) -> Self {
        self.interpolation = mode;
        self
    }
    pub fn safety_margin(mut self, left: f64, right: f64) -> Self {
        self.margin_left = left;
        self.margin_right = right;

        self
    }
    pub fn with_dt_matrix(mut self) -> Self {
        self.extras.dt_matrix = true;

        self
    }

    pub fn compute_dt(self, compute: bool) -> Self {
        if compute {
            return self.with_dt_matrix();
        } else {
            return self;
        }
    }
    // for Vyazovkin method
    pub fn compute_dt_matrix(time: &Array2<f64>) -> Array2<f64> {
        let (n_exp, n_alpha) = time.dim();

        let mut dt = Array2::<f64>::zeros((n_exp, n_alpha));

        for i in 0..n_exp {
            let row = time.row(i);

            for j in 0..n_alpha - 1 {
                dt[[i, j]] = row[j + 1] - row[j];
            }

            // последнюю колонку можно повторить
            dt[[i, n_alpha - 1]] = dt[[i, n_alpha - 2]];
        }

        dt
    }

    fn calc_auto_eta_range(&self, data: &KineticDataView) -> Result<(f64, f64), TGADomainError> {
        let mut global_min = f64::NEG_INFINITY;
        let mut global_max = f64::INFINITY;

        for eta in data.experiments.iter().map(|e| &e.conversion) {
            let local_min = eta.iter().cloned().fold(f64::INFINITY, f64::min);

            let local_max = eta.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            global_min = global_min.max(local_min);
            global_max = global_max.min(local_max);
        }

        if global_max <= global_min {
            return Err(TGADomainError::InvalidConversionRange);
        }

        Ok((global_min, global_max))
    }

    pub fn build_withcopy(self, data: &KineticDataView) -> Result<ConversionGrid, TGADomainError> {
        let (mut eta_min, mut eta_max) = match self.eta_mode {
            AlphaRangeMode::Auto => {
                let range = self.calc_auto_eta_range(data)?;
                (range.0, range.1)
            }
            AlphaRangeMode::Manual { min, max } => (min, max),
        };
        eta_min += self.margin_left;
        eta_max -= self.margin_right;

        if eta_max <= eta_min {
            return Err(TGADomainError::InvalidConversionRange);
        }

        let eta_grid = Self::build_eta_grid(eta_min, eta_max, self.segments);

        let n_exp = data.experiments.len();
        let n_eta = eta_grid.len();

        let mut temperature = Vec::with_capacity(n_exp * n_eta);
        let mut time = Vec::with_capacity(n_exp * n_eta);
        let mut rate = Vec::with_capacity(n_exp * n_eta);

        for exp in &data.experiments {
            let (t_vec, temp_vec, rate_vec) =
                Self::interpolate_experiment(exp, &eta_grid, self.interpolation)?;

            if t_vec.len() != n_eta || temp_vec.len() != n_eta || rate_vec.len() != n_eta {
                return Err(TGADomainError::InvalidOperation(
                    "Interpolation output has unexpected length".into(),
                ));
            }

            time.extend(t_vec);
            temperature.extend(temp_vec);
            rate.extend(rate_vec);
        }

        let temperature = Array2::from_shape_vec((n_exp, n_eta), temperature).map_err(|e| {
            TGADomainError::InvalidOperation(format!("Invalid temperature grid: {}", e))
        })?;
        let time = Array2::from_shape_vec((n_exp, n_eta), time)
            .map_err(|e| TGADomainError::InvalidOperation(format!("Invalid time grid: {}", e)))?;
        let conversion_rate = Array2::from_shape_vec((n_exp, n_eta), rate)
            .map_err(|e| TGADomainError::InvalidOperation(format!("Invalid rate grid: {}", e)))?;

        let inv_temperature = temperature.mapv(|t| 1.0 / t);
        let dt = if self.extras.dt_matrix {
            Some(Self::compute_dt_matrix(&time))
        } else {
            None
        };

        Ok(ConversionGrid {
            eta: Array1::from(eta_grid),
            temperature,
            inv_temperature,
            time,
            conversion_rate,
            dt,
            meta: data.experiments.iter().map(|e| e.meta.clone()).collect(),
        })
    }

    pub fn build(self, data: &KineticDataView) -> Result<ConversionGrid, TGADomainError> {
        let (mut eta_min, mut eta_max) = match self.eta_mode {
            AlphaRangeMode::Auto => self.calc_auto_eta_range(data)?,
            AlphaRangeMode::Manual { min, max } => (min, max),
        };

        eta_min += self.margin_left;
        eta_max -= self.margin_right;

        let eta_grid = Self::build_eta_grid(eta_min, eta_max, self.segments);

        let n_exp = data.experiments.len();
        let n_eta = eta_grid.len();

        // allocate arrays once
        let mut temperature = Array2::<f64>::zeros((n_exp, n_eta));
        let mut time = Array2::<f64>::zeros((n_exp, n_eta));
        let mut rate = Array2::<f64>::zeros((n_exp, n_eta));

        for (i, exp) in data.experiments.iter().enumerate() {
            let mut t_row = time.row_mut(i);
            let mut T_row = temperature.row_mut(i);
            let mut r_row = rate.row_mut(i);

            Self::interpolate_experiment_into(
                exp,
                &eta_grid,
                self.interpolation,
                &mut t_row,
                &mut T_row,
                &mut r_row,
            )?;
        }

        let inv_temperature = temperature.mapv(|t| 1.0 / t);

        let dt = if self.extras.dt_matrix {
            Some(Self::compute_dt_matrix(&time))
        } else {
            None
        };

        Ok(ConversionGrid {
            eta: Array1::from(eta_grid),
            temperature,
            inv_temperature,
            time,
            conversion_rate: rate,
            dt,
            meta: data.experiments.iter().map(|e| e.meta.clone()).collect(),
        })
    }

    fn build_eta_grid(min: f64, max: f64, n: usize) -> Vec<f64> {
        let step = (max - min) / n as f64;

        (0..n).map(|i| min + step * i as f64).collect()
    }

    fn interpolate_experiment(
        exp: &ExperimentData,
        eta_grid: &[f64],
        method: GridInterpolation,
    ) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), TGADomainError> {
        match method {
            GridInterpolation::Linear => interpolate_linear(exp, eta_grid),

            GridInterpolation::Spline => interpolate_spline(exp, eta_grid),
        }
    }

    pub fn interpolate_experiment_into(
        exp: &ExperimentData,
        eta_grid: &[f64],
        method: GridInterpolation,
        t_out: &mut ndarray::ArrayViewMut1<f64>,
        T_out: &mut ndarray::ArrayViewMut1<f64>,
        r_out: &mut ndarray::ArrayViewMut1<f64>,
    ) -> Result<(), TGADomainError> {
        match method {
            GridInterpolation::Linear => {
                interpolate_linear_into(exp, eta_grid, t_out, T_out, r_out)
            }

            GridInterpolation::Spline => {
                interpolate_spline_into(exp, eta_grid, t_out, T_out, r_out)
            }
        }
    }

    pub fn build_isothermal(
        self,
        data: &KineticDataView,
    ) -> Result<ConversionGrid, TGADomainError> {
        let (mut eta_min, mut eta_max) = match self.eta_mode {
            AlphaRangeMode::Auto => self.calc_auto_eta_range(data)?,
            AlphaRangeMode::Manual { min, max } => (min, max),
        };

        eta_min += self.margin_left;
        eta_max -= self.margin_right;

        let eta_grid = Self::build_eta_grid(eta_min, eta_max, self.segments);

        let n_exp = data.experiments.len();
        let n_eta = eta_grid.len();

        let mut temperature = Array2::<f64>::zeros((n_exp, n_eta));
        let mut time = Array2::<f64>::zeros((n_exp, n_eta));

        for (i, exp) in data.experiments.iter().enumerate() {
            let mut t_row = time.row_mut(i);
            let mut T_row = temperature.row_mut(i);

            Self::interpolate_time_only(exp, &eta_grid, self.interpolation, &mut t_row)?;

            let t_iso = exp.meta.isothermal_temperature.ok_or_else(|| {
                TGADomainError::InvalidOperation("Isothermal temperature missing".into())
            })?;

            for j in 0..n_eta {
                T_row[j] = t_iso;
            }
        }

        let inv_temperature = temperature.mapv(|t| 1.0 / t);

        Ok(ConversionGrid {
            eta: Array1::from(eta_grid),

            temperature,
            inv_temperature,

            time,

            conversion_rate: Array2::zeros((n_exp, n_eta)),

            dt: None,

            meta: data.experiments.iter().map(|e| e.meta.clone()).collect(),
        })
    }

    pub fn interpolate_time_only(
        exp: &ExperimentData,
        eta_grid: &[f64],
        method: GridInterpolation,
        t_out: &mut ndarray::ArrayViewMut1<f64>,
    ) -> Result<(), TGADomainError> {
        match method {
            GridInterpolation::Linear => interpolate_time_linear_only(exp, eta_grid, t_out),

            GridInterpolation::Spline => interpolate_time_spline_only(exp, eta_grid, t_out),
        }
    }
}

fn interpolate_linear(
    exp: &ExperimentData,
    eta_grid: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), TGADomainError> {
    let eta = &exp.conversion;
    let time = &exp.time;
    let temp = &exp.temperature;
    let rate = &exp.conversion_rate;

    let mut out_t = Vec::with_capacity(eta_grid.len());
    let mut out_T = Vec::with_capacity(eta_grid.len());
    let mut out_r = Vec::with_capacity(eta_grid.len());

    for &a in eta_grid {
        let idx = eta
            .windows(2)
            .position(|w| w[0] <= a && w[1] >= a)
            .ok_or_else(|| TGADomainError::InvalidOperation("eta out of bounds".into()))?;

        let a0 = eta[idx];
        let a1 = eta[idx + 1];

        let w = (a - a0) / (a1 - a0);

        let t = time[idx] + w * (time[idx + 1] - time[idx]);
        let T = temp[idx] + w * (temp[idx + 1] - temp[idx]);
        let r = rate[idx] + w * (rate[idx + 1] - rate[idx]);

        out_t.push(t);
        out_T.push(T);
        out_r.push(r);
    }

    Ok((out_t, out_T, out_r))
}

fn interpolate_spline(
    exp: &ExperimentData,
    eta_grid: &[f64],
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), TGADomainError> {
    let eta = &exp.conversion;

    let interp = if eta.len() < 4 {
        Interpolation::Linear
    } else {
        Interpolation::CatmullRom
    };

    let time_keys: Vec<_> = eta
        .iter()
        .zip(&exp.time)
        .map(|(&a, &t)| Key::new(a, t, interp))
        .collect();

    let temp_keys: Vec<_> = eta
        .iter()
        .zip(&exp.temperature)
        .map(|(&a, &t)| Key::new(a, t, interp))
        .collect();

    let rate_keys: Vec<_> = eta
        .iter()
        .zip(&exp.conversion_rate)
        .map(|(&a, &r)| Key::new(a, r, interp))
        .collect();

    let spline_time = Spline::from_vec(time_keys);
    let spline_temp = Spline::from_vec(temp_keys);
    let spline_rate = Spline::from_vec(rate_keys);

    let mut t = Vec::new();
    let mut T = Vec::new();
    let mut r = Vec::new();

    for &a in eta_grid {
        t.push(spline_time.sample(a).ok_or_else(|| {
            TGADomainError::InvalidOperation("spline sample out of bounds".into())
        })?);
        T.push(spline_temp.sample(a).ok_or_else(|| {
            TGADomainError::InvalidOperation("spline sample out of bounds".into())
        })?);
        r.push(spline_rate.sample(a).ok_or_else(|| {
            TGADomainError::InvalidOperation("spline sample out of bounds".into())
        })?);
    }

    Ok((t, T, r))
}

pub fn interpolate_linear_into(
    exp: &ExperimentData,
    eta_grid: &[f64],
    t_out: &mut ndarray::ArrayViewMut1<f64>,
    T_out: &mut ndarray::ArrayViewMut1<f64>,
    r_out: &mut ndarray::ArrayViewMut1<f64>,
) -> Result<(), TGADomainError> {
    let eta = &exp.conversion;
    let time = &exp.time;
    let temp = &exp.temperature;
    let rate = &exp.conversion_rate;
    if !eta.windows(2).all(|w| w[0] <= w[1]) {
        return Err(TGADomainError::InvalidOperation(
            "Conversion must be  monotonic".into(),
        ));
    }
    /*
     for (j, &a) in eta_grid.iter().enumerate() {
        let idx = eta
            .windows(2)
            .position(|w| w[0] <= a && w[1] >= a)
            .ok_or_else(|| TGADomainError::InvalidOperation("eta out of bounds".into()))?;

        let a0 = eta[idx];
        let a1 = eta[idx + 1];

        let w = (a - a0) / (a1 - a0);

        t_out[j] = time[idx] + w * (time[idx + 1] - time[idx]);
        T_out[j] = temp[idx] + w * (temp[idx + 1] - temp[idx]);
        r_out[j] = rate[idx] + w * (rate[idx + 1] - rate[idx]);
    }

     */
    let mut i = 0;
    for (j, &a) in eta_grid.iter().enumerate() {
        while i + 1 < eta.len() && eta[i + 1] < a {
            i += 1;
        }

        if i + 1 >= eta.len() {
            return Err(TGADomainError::InvalidOperation(
                "eta grid outside experiment range".into(),
            ));
        }

        let a0 = eta[i];
        let a1 = eta[i + 1];

        let w = (a - a0) / (a1 - a0);

        t_out[j] = time[i] + w * (time[i + 1] - time[i]);
        T_out[j] = temp[i] + w * (temp[i + 1] - temp[i]);
        r_out[j] = rate[i] + w * (rate[i + 1] - rate[i]);
    }

    Ok(())
}

fn interpolate_spline_into(
    exp: &ExperimentData,
    eta_grid: &[f64],
    t_out: &mut ndarray::ArrayViewMut1<f64>,
    T_out: &mut ndarray::ArrayViewMut1<f64>,
    r_out: &mut ndarray::ArrayViewMut1<f64>,
) -> Result<(), TGADomainError> {
    let eta = &exp.conversion;
    if !eta.windows(2).all(|w| w[0] <= w[1]) {
        return Err(TGADomainError::InvalidOperation(
            "Conversion must be monotonic".into(),
        ));
    }

    let time = &exp.time;
    let temp = &exp.temperature;
    let rate = &exp.conversion_rate;

    if eta.len() < 3 {
        return Err(TGADomainError::InvalidOperation(
            "Not enough points for spline interpolation".into(),
        ));
    }

    // --- build spline keys

    let keys_t: Vec<_> = eta
        .iter()
        .zip(time)
        .map(|(&a, &t)| Key::new(a, t, Interpolation::CatmullRom))
        .collect();

    let keys_T: Vec<_> = eta
        .iter()
        .zip(temp)
        .map(|(&a, &T)| Key::new(a, T, Interpolation::CatmullRom))
        .collect();

    let keys_r: Vec<_> = eta
        .iter()
        .zip(rate)
        .map(|(&a, &r)| Key::new(a, r, Interpolation::CatmullRom))
        .collect();

    let spline_t = Spline::from_vec(keys_t);
    let spline_T = Spline::from_vec(keys_T);
    let spline_r = Spline::from_vec(keys_r);

    // --- evaluate

    for (j, &a) in eta_grid.iter().enumerate() {
        t_out[j] = spline_t
            .sample(a)
            .ok_or_else(|| TGADomainError::InvalidOperation("Spline t failed".into()))?;

        T_out[j] = spline_T
            .sample(a)
            .ok_or_else(|| TGADomainError::InvalidOperation("Spline T failed".into()))?;

        r_out[j] = spline_r
            .sample(a)
            .ok_or_else(|| TGADomainError::InvalidOperation("Spline rate failed".into()))?;
    }

    Ok(())
}

pub fn interpolate_time_linear_only(
    exp: &ExperimentData,
    eta_grid: &[f64],
    t_out: &mut ndarray::ArrayViewMut1<f64>,
) -> Result<(), TGADomainError> {
    let eta = &exp.conversion;
    let time = &exp.time;

    if !eta.windows(2).all(|w| w[0] <= w[1]) {
        return Err(TGADomainError::InvalidOperation(
            "Conversion must be monotonic".into(),
        ));
    }

    let mut i = 0;

    for (j, &a) in eta_grid.iter().enumerate() {
        while i + 1 < eta.len() && eta[i + 1] < a {
            i += 1;
        }

        if i + 1 >= eta.len() {
            return Err(TGADomainError::InvalidOperation(
                "eta grid outside experiment range".into(),
            ));
        }

        let a0 = eta[i];
        let a1 = eta[i + 1];

        let w = (a - a0) / (a1 - a0);

        t_out[j] = time[i] + w * (time[i + 1] - time[i]);
    }

    Ok(())
}

fn interpolate_time_spline_only(
    exp: &ExperimentData,
    eta_grid: &[f64],
    t_out: &mut ndarray::ArrayViewMut1<f64>,
) -> Result<(), TGADomainError> {
    let eta = &exp.conversion;
    let time = &exp.time;

    if eta.len() < 3 {
        return Err(TGADomainError::InvalidOperation(
            "Not enough points for spline interpolation".into(),
        ));
    }

    let keys: Vec<_> = eta
        .iter()
        .zip(time)
        .map(|(&a, &t)| Key::new(a, t, Interpolation::CatmullRom))
        .collect();

    let spline = Spline::from_vec(keys);

    for (j, &a) in eta_grid.iter().enumerate() {
        t_out[j] = spline
            .sample(a)
            .ok_or_else(|| TGADomainError::InvalidOperation("Spline time failed".into()))?;
    }

    Ok(())
}

pub enum KineticMethodKind {
    OFW,
    KAS,
    Starink,
    Friedman,
}
//=====================================================================================================
// TESTS
//=====================================================================================================

#[cfg(test)]
pub mod tests {
    use super::*;

    use crate::Kinetics::experimental_kinetics::testing_mod::tests_afvanced_config::{
        base_advanced_config_isothermal, base_advanced_config_non_isothermal,
        build_series_from_cfg, build_series_from_cfg_with_m0,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::{AdvancedTGAConfig, KineticModel};
    use approx::assert_relative_eq;

    fn mock_experiment() -> ExperimentData {
        ExperimentData {
            meta: ExperimentMeta::default(),
            time: vec![0.0, 5.0, 10.0],
            temperature: vec![300.0, 350.0, 400.0],
            conversion: vec![0.0, 0.5, 1.0],
            conversion_rate: vec![0.1, 0.2, 0.3],
            mass: None,
            mass_rate: None,
        }
    }

    #[test]
    fn interpolate_linear_matches_expected() {
        let exp = mock_experiment();
        let eta_grid = vec![0.25, 0.75];

        let (t, temp, rate) = interpolate_linear(&exp, &eta_grid).unwrap();

        assert_relative_eq!(t[0], 2.5, epsilon = 1e-12);
        assert_relative_eq!(t[1], 7.5, epsilon = 1e-12);
        assert_relative_eq!(temp[0], 325.0, epsilon = 1e-12);
        assert_relative_eq!(temp[1], 375.0, epsilon = 1e-12);
        assert_relative_eq!(rate[0], 0.15, epsilon = 1e-12);
        assert_relative_eq!(rate[1], 0.25, epsilon = 1e-12);
    }

    #[test]
    fn interpolate_spline_is_consistent_for_linear_data() {
        let exp = mock_experiment();
        let eta_grid = vec![0.25, 0.75];

        let (t, temp, rate) = interpolate_spline(&exp, &eta_grid).unwrap();

        assert_relative_eq!(t[0], 2.5, epsilon = 1e-8);
        assert_relative_eq!(t[1], 7.5, epsilon = 1e-8);
        assert_relative_eq!(temp[0], 325.0, epsilon = 1e-8);
        assert_relative_eq!(temp[1], 375.0, epsilon = 1e-8);
        assert_relative_eq!(rate[0], 0.15, epsilon = 1e-8);
        assert_relative_eq!(rate[1], 0.25, epsilon = 1e-8);
    }

    pub fn build_view_from_cfg(cfg: &AdvancedTGAConfig) -> Result<KineticDataView, TGADomainError> {
        let series = build_series_from_cfg(cfg, -1.0, 0.0, 1e-4)?;

        let united = series
            .concat_into_vertical_stack(
                None,
                vec![
                    ColumnNature::Time,
                    ColumnNature::Conversion,
                    ColumnNature::Temperature,
                    ColumnNature::ConversionRate,
                ],
            )
            .unwrap();

        KineticDataView::from_united_dataset(&united)
    }

    pub fn build_view_from_cfg_exact_m0(
        cfg: &AdvancedTGAConfig,
    ) -> Result<KineticDataView, TGADomainError> {
        let m0_raw = match cfg.kinetic_model {
            KineticModel::ArrheniusSingle { m0, .. } => m0,
            KineticModel::ArrheniusTwoComponent { m01, m02, .. } => m01 + m02,
        };
        let k = -1.0;
        let b = 0.0;
        let m0_calibrated = (k * m0_raw + b).abs();

        let series = build_series_from_cfg_with_m0(cfg, k, b, m0_calibrated)?;

        let united = series
            .concat_into_vertical_stack(
                None,
                vec![
                    ColumnNature::Time,
                    ColumnNature::Conversion,
                    ColumnNature::Temperature,
                    ColumnNature::ConversionRate,
                ],
            )
            .unwrap();

        KineticDataView::from_united_dataset(&united)
    }

    fn assert_basic_view_sanity(view: &KineticDataView, expected_experiments: usize) {
        assert_eq!(view.experiments.len(), expected_experiments);
        for exp in &view.experiments {
            assert!(!exp.time.is_empty());
            assert_eq!(exp.time.len(), exp.temperature.len());
            assert_eq!(exp.time.len(), exp.conversion.len());
            assert_eq!(exp.time.len(), exp.conversion_rate.len());
            //    println!("mass {:?}", &exp.mass.unwrap()[0..10].unwrap());
            println!("conversion {:?}", &exp.conversion[0..10]);
            assert!(exp.conversion.iter().all(|&v| v >= -2e-2 && v <= 1.0));

            for w in exp.time.windows(2) {
                assert!(w[1] > w[0]);
            }
        }
    }

    #[test]
    fn test_with_virtual_tga() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 3);

        let grid = ConversionGridBuilder::new()
            .eta_range(0.0, 1.0)
            .segments(200)
            .interpolation(GridInterpolation::Linear)
            .build(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 200);
        assert_eq!(grid.temperature.dim(), (3, 200));
        assert_eq!(grid.time.dim(), (3, 200));
        assert_eq!(grid.conversion_rate.dim(), (3, 200));
        assert_eq!(grid.meta.len(), 3);
        assert_relative_eq!(grid.eta[0], 0.0, epsilon = 1e-12);
        assert!(grid.eta[199] < 1.0);
    }

    #[test]
    fn test_with_virtual_tga_with_auto_range() {
        let cfg = base_advanced_config_non_isothermal(
            700.0,
            1e5,
            80_000.0,
            0.1,
            10_000,
            vec![0.5, 3.0, 5.0],
        );
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 3);

        let grid = ConversionGridBuilder::new()
            .auto_range()
            .segments(200)
            .interpolation(GridInterpolation::Linear)
            .build(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 200);
        assert_eq!(grid.temperature.dim(), (3, 200));
        assert_eq!(grid.time.dim(), (3, 200));
        assert_eq!(grid.conversion_rate.dim(), (3, 200));
        assert_eq!(grid.meta.len(), 3);
        assert_relative_eq!(grid.eta[0], 0.0, epsilon = 1e-1);
        assert!(grid.eta[199] < 1.0);
    }
    #[test]
    fn test_with_virtual_tga_spline_grid_isothermal() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, vec![600.0, 700.0], 0.1, 10_000);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);
        for exp in &view.experiments {
            assert!(exp.meta.isothermal_temperature.is_some());
            assert!(exp.meta.heating_rate.is_none());
        }

        let grid = ConversionGridBuilder::new()
            .eta_range(0.1, 0.9)
            .segments(120)
            .interpolation(GridInterpolation::Spline)
            .build(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 120);
        assert_eq!(grid.temperature.dim(), (2, 120));
        assert_eq!(grid.time.dim(), (2, 120));
        assert_eq!(grid.conversion_rate.dim(), (2, 120));
        assert_relative_eq!(grid.eta[0], 0.1, epsilon = 1e-12);
        assert!(grid.eta[119] < 0.9);
    }

    #[test]
    fn test_with_virtual_tga_line_grid_isothermal() {
        let cfg = base_advanced_config_isothermal(1e5, 80_000.0, vec![600.0, 700.0], 0.1, 10_000);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);
        for exp in &view.experiments {
            assert!(exp.meta.isothermal_temperature.is_some());
            assert!(exp.meta.heating_rate.is_none());
        }

        let grid = ConversionGridBuilder::new()
            .eta_range(0.1, 0.9)
            .segments(120)
            .interpolation(GridInterpolation::Linear)
            .build(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 120);
        assert_eq!(grid.temperature.dim(), (2, 120));
        assert_eq!(grid.time.dim(), (2, 120));
        assert_eq!(grid.conversion_rate.dim(), (2, 120));
        assert_relative_eq!(grid.eta[0], 0.1, epsilon = 1e-12);
        assert!(grid.eta[119] < 0.9);
    }
    #[test]
    fn test_with_virtual_tga_narrow_eta_grid() {
        let cfg =
            base_advanced_config_non_isothermal(500.0, 1e5, 80_000.0, 0.1, 10_000, vec![3.5, 4.5]);
        let view = build_view_from_cfg(&cfg).unwrap();

        assert_basic_view_sanity(&view, 2);

        let grid = ConversionGridBuilder::new()
            .eta_range(0.15, 0.75)
            .segments(60)
            .interpolation(GridInterpolation::Linear)
            .build(&view)
            .unwrap();

        assert_eq!(grid.eta.len(), 60);
        assert_eq!(grid.temperature.dim(), (2, 60));
        assert_eq!(grid.time.dim(), (2, 60));
        assert_eq!(grid.conversion_rate.dim(), (2, 60));
        assert_relative_eq!(grid.eta[0], 0.15, epsilon = 1e-12);
        assert!(grid.eta[59] < 0.75);
        for w in grid.eta.windows(2) {
            assert!(w[1] > w[0]);
        }
    }

    #[test]
    fn test_calc_auto_eta_range() {
        // Helper to create a simple experiment with given conversion range
        fn make_experiment(conversion: Vec<f64>) -> ExperimentData {
            let n = conversion.len();
            ExperimentData {
                meta: ExperimentMeta::default(),
                time: (0..n).map(|i| i as f64).collect(),
                temperature: vec![300.0; n],
                conversion,
                conversion_rate: vec![0.1; n],
                mass: None,
                mass_rate: None,
            }
        }

        // Single experiment
        let view = KineticDataView {
            experiments: vec![make_experiment(vec![0.1, 0.2, 0.3, 0.4])],
        };
        let builder = ConversionGridBuilder::new();
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.1, epsilon = 1e-12);
        assert_relative_eq!(max, 0.4, epsilon = 1e-12);

        // Two experiments with overlapping ranges
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.2, 0.5, 0.8]),
                make_experiment(vec![0.3, 0.6, 0.9]),
            ],
        };
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.3, epsilon = 1e-12); // intersection min = max of mins (0.2,0.3) = 0.3
        assert_relative_eq!(max, 0.8, epsilon = 1e-12); // intersection max = min of maxs (0.8,0.9) = 0.8

        // Three experiments, one with narrower range
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.0, 0.5, 1.0]),
                make_experiment(vec![0.2, 0.3, 0.4]),
                make_experiment(vec![0.1, 0.6, 0.7]),
            ],
        };
        let (min, max) = builder.calc_auto_eta_range(&view).unwrap();
        assert_relative_eq!(min, 0.2, epsilon = 1e-12); // max of mins: max(0.0,0.2,0.1) = 0.2
        assert_relative_eq!(max, 0.4, epsilon = 1e-12);

        // Non-overlapping ranges -> error
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.1, 0.2]),
                make_experiment(vec![0.5, 0.6]),
            ],
        };
        let result = builder.calc_auto_eta_range(&view);
        assert!(matches!(
            result,
            Err(TGADomainError::InvalidConversionRange)
        ));

        // Touching ranges (max of first equals min of second) -> intersection is a point, should error
        let view = KineticDataView {
            experiments: vec![
                make_experiment(vec![0.1, 0.3]),
                make_experiment(vec![0.3, 0.5]),
            ],
        };
        let result = builder.calc_auto_eta_range(&view);
        assert!(matches!(
            result,
            Err(TGADomainError::InvalidConversionRange)
        ));

        // Empty conversion vector (should not happen in practice, but test edge case)
        // The fold will produce INFINITY and NEG_INFINITY, leading to global_min = NEG_INFINITY, global_max = INFINITY? Actually local_min = INFINITY, local_max = NEG_INFINITY.
        // Then global_min = max(NEG_INFINITY, INFINITY) = INFINITY? Wait, global_min starts as NEG_INFINITY, then max with INFINITY yields INFINITY.
        // global_max starts as INFINITY, then min with NEG_INFINITY yields NEG_INFINITY.
        // At the end global_max <= global_min? INFINITY <= NEG_INFINITY? false. But we'll skip this test as it's unrealistic.
    }
}
