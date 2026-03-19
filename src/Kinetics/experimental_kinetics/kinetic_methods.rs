//! # Kinetic methods — shared infrastructure
//!
//! This module is the **foundation layer** for all isoconversional kinetic
//! analysis in KiThe.  It defines the data types, the trait, the validation
//! helpers, and the grid-building machinery that every concrete method
//! (OFW, KAS, Starink, Friedman, Vyazovkin) builds on.
//!
//! ## Architecture overview
//!
//! ```text
//! TGASeries / UnitedDataset
//!        │
//!        ▼
//! KineticDataView          ← one ExperimentData per experiment
//!        │
//!        ▼
//! ConversionGridBuilder    ← builder pattern; configures α-range,
//!   ├─ build()             │  segment count, interpolation, extras
//!   ├─ build_isothermal()  │
//!   └─ build_universal()   │
//!        │
//!        ▼
//! ConversionGrid           ← ndarray matrices [n_exp × n_α]
//!        │
//!        ▼
//! KineticMethod::compute   ← trait implemented by every solver
//!        │
//!        ▼
//! IsoconversionalResult    ← Vec<IsoLayerResult>, one per α point
//! ```
//!
//! ## Sub-modules
//! | Module | Contents |
//! |--------|----------|
//! | `integral_isoconversion` | OFW, KAS, Starink solvers + `VyazovkinMethod` adapter |
//! | `Friedman` | Differential and integral Friedman solvers |
//! | `Vyazovkin` | Non-linear Vyazovkin solver |
//! | `isoconversion` | Unified dispatcher `IsoconversionalSolver` |
//! | `kinetic_regression` | OLS linear regression used by all linear methods |
//!
//! ## Main data structures
//! - [`KineticDataView`] — a thin, materialised view over a series of
//!   experiments; the universal input to every `KineticMethod`.
//! - [`ExperimentData`] — one experiment's time-series vectors plus metadata.
//! - [`KineticRequirements`] — declarative checklist (min experiments,
//!   needs heating rate, etc.) returned by each method.
//! - [`ConversionGrid`] — the central ndarray structure: all experiments
//!   resampled onto a common α grid, with `1/T` pre-cached.
//! - [`ConversionGridBuilder`] — fluent builder that interpolates each
//!   experiment onto the α grid and assembles `ConversionGrid`.
//!
//! ## Non-trivial techniques
//!
//! ### α-grid intersection for `auto_range`
//! When no explicit α range is given, `calc_auto_eta_range` computes the
//! *intersection* of all experiments' conversion ranges:
//! `η_min = max(local_mins)`, `η_max = min(local_maxs)`.  This guarantees
//! every experiment has data at every grid point, preventing out-of-bounds
//! interpolation without requiring the caller to know the individual ranges.
//!
//! ### Single-pass linear interpolation with a sliding pointer
//! `interpolate_linear_into` uses a single forward pointer `i` that only
//! advances, exploiting the fact that both the source conversion vector and
//! the target α grid are monotonically increasing.  This reduces the
//! per-experiment interpolation from O(n·m) (naïve search) to O(n + m).
//!
//! ### Pre-cached `1/T` matrix
//! `ConversionGrid::inv_temperature` stores `1/T` for every `(experiment, α)`
//! cell at build time.  All solvers read this field directly, avoiding
//! repeated division in the inner regression loop.
//!
//! ### `build_isothermal` fills temperature from metadata
//! For isothermal experiments the measured temperature column is noisy and
//! irrelevant to the isoconversional regression.  `build_isothermal` fills
//! each row of the temperature matrix with the single constant value stored
//! in `ExperimentMeta::isothermal_temperature`, giving a perfectly clean
//! `1/T` column for the solver.
pub mod Criado;
pub mod Criado2;
pub mod Friedman;
pub mod Kissinger;
pub mod Vyazovkin;
pub mod combined;
pub mod integral_isoconversion;
mod integral_isoconversion_tests;
pub mod is_this_a_sublimation;
pub mod isoconversion;
pub mod kinetic_regression;
use crate::Kinetics::experimental_kinetics::experiment_series_main::ExperimentMeta;
use crate::Kinetics::experimental_kinetics::experiment_series2::UnitedDataset;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    ColumnNature, TGADomainError,
};
use ndarray::{Array1, Array2};
use splines::{Interpolation, Key, Spline};

/// Materialised, read-only view over a series of TGA experiments.
///
/// This is the universal input type for every `KineticMethod`.  It is
/// constructed either from a `UnitedDataset` (production path) or assembled
/// directly in tests.  Each element of `experiments` corresponds to one
/// experimental run (one heating rate or one isothermal temperature).
#[derive(Debug, Clone)]
pub struct KineticDataView {
    /// One entry per experiment, ordered as they appear in the source dataset.
    pub experiments: Vec<ExperimentData>,
}
impl KineticDataView {
    /// Constructs a `KineticDataView` by materialising all required columns
    /// (time, temperature, conversion, conversion rate) from a `UnitedDataset`.
    ///
    /// Returns an error if any column is missing for any experiment.
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
    pub fn from_united_dataset_by_nature(
        ds: &UnitedDataset,
        nature: Vec<ColumnNature>,
    ) -> Result<Self, TGADomainError> {
        let mut experiments = Vec::new();

        for meta in &ds.meta {
            let time = if nature.contains(&ColumnNature::Time) {
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Time)?
            } else {
                Vec::new()
            };

            let temperature = if nature.contains(&ColumnNature::Temperature) {
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Temperature)?
            } else {
                Vec::new()
            };

            let conversion = if nature.contains(&ColumnNature::Conversion) {
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Conversion)?
            } else {
                Vec::new()
            };

            let conversion_rate = if nature.contains(&ColumnNature::ConversionRate) {
                ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::ConversionRate)?
            } else {
                Vec::new()
            };

            let mass = if nature.contains(&ColumnNature::Mass) {
                Some(ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::Mass)?)
            } else {
                None
            };

            let mass_rate = if nature.contains(&ColumnNature::MassRate) {
                Some(ds.materialize_vertical_column_by_nature(&meta.id, ColumnNature::MassRate)?)
            } else {
                None
            };
            experiments.push(ExperimentData {
                meta: meta.clone(),
                time,
                temperature,
                conversion,
                conversion_rate,
                mass,
                mass_rate,
            });
        }

        Ok(Self { experiments })
    }
}
/// Materialised time-series data for a single TGA experiment.
///
/// All vectors are guaranteed to have the same length after construction.
/// `mass` and `mass_rate` are optional because most kinetic methods work
/// exclusively with the derived `conversion` and `conversion_rate` columns.
#[derive(Debug, Clone)]
pub struct ExperimentData {
    /// Experiment-level metadata (id, heating rate, isothermal temperature).
    pub meta: ExperimentMeta,
    /// Time axis (seconds).
    pub time: Vec<f64>,
    /// Temperature axis (Kelvin).
    pub temperature: Vec<f64>,
    /// Conversion α ∈ [0, 1].
    pub conversion: Vec<f64>,
    /// dα/dt (s⁻¹).
    pub conversion_rate: Vec<f64>,
    /// Raw sample mass (optional; not used by isoconversional solvers).
    pub mass: Option<Vec<f64>>,
    /// dm/dt (optional).
    pub mass_rate: Option<Vec<f64>>,
}
/// Common trait implemented by every kinetic method in this crate.
///
/// A method receives a [`KineticDataView`] and returns its `Output` type
/// (typically [`IsoconversionalResult`]).  The default `check_input`
/// implementation is a no-op; override it for method-specific validation
/// that goes beyond the declarative [`KineticRequirements`] checks.
pub trait KineticMethod {
    /// The result type produced by this method.
    type Output;

    /// Human-readable name of the method (e.g. `"KAS"`).
    fn name(&self) -> &'static str;

    /// Core computation: build the grid and run the solver.
    fn compute(&self, data: &KineticDataView) -> Result<Self::Output, TGADomainError>;

    /// Optional method-specific input validation, called before `compute`.
    /// The default implementation accepts any input.
    fn check_input(&self, _data: &KineticDataView) -> Result<(), TGADomainError> {
        Ok(())
    }

    /// Returns the declarative requirements for this method so that
    /// `run_method_with_requirements` can validate the view automatically.
    fn requirements(&self) -> KineticRequirements;

    fn required_columns_by_nature(&self) -> Vec<ColumnNature> {
        Vec::new()
    }
}

impl<O> dyn KineticMethod<Output = O> {
    /// Convenience runner: validates input via `check_input`, then calls `compute`.
    pub fn run<M: KineticMethod<Output = O>>(
        method: &M,
        data: &KineticDataView,
    ) -> Result<O, TGADomainError> {
        method.check_input(data)?;

        method.compute(data)
    }
}

/// Returns an error if the view contains fewer than `n` experiments.
pub fn require_min_experiments(data: &KineticDataView, n: usize) -> Result<(), TGADomainError> {
    if data.experiments.len() < n {
        return Err(TGADomainError::InvalidOperation(format!(
            "Method requires at least {} experiments",
            n
        )));
    }

    Ok(())
}

/// Returns an error if any experiment is missing a heating rate.
/// Required by non-isothermal methods (OFW, KAS, Starink, Differential Friedman).
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

/// Returns an error if any experiment is missing an isothermal temperature.
/// Required by isothermal methods (Integral Friedman).
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

/// Returns an error if any experiment has an empty `conversion_rate` vector.
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

/// Returns an error if any experiment has an empty `conversion` vector.
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

/// Declarative checklist of what a kinetic method needs from the input view.
///
/// Returned by `KineticMethod::requirements`; consumed by
/// `check_requirements` / `run_method_with_requirements` to validate the
/// view before calling `compute`.
#[derive(Default)]
pub struct KineticRequirements {
    /// Minimum number of experiments required (typically 3).
    pub min_experiments: usize,
    /// Whether every experiment must have a non-empty `conversion` vector.
    needs_conversion: bool,
    /// Whether every experiment must have a non-empty `conversion_rate` vector.
    needs_conversion_rate: bool,
    /// Whether every experiment must have a temperature column
    /// (used to gate `require_isothermal` for isothermal methods).
    needs_temperature: bool,
    /// Whether every experiment must carry a heating rate in its metadata.
    needs_heating_rate: bool,
}
/// Runs all declarative checks encoded in `req` against `data`.
///
/// Called automatically by `run_method_with_requirements`; can also be
/// invoked manually when composing custom pipelines.
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

/// Preferred entry point for running a kinetic method.
///
/// Calls `check_requirements`, then `check_input`, then `compute` in order,
/// returning the first error encountered.
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
/// Selects the interpolation strategy used when resampling experiments onto
/// the common α grid.
#[derive(Clone, Copy, Debug)]
pub enum GridInterpolation {
    /// Piecewise linear interpolation — fast, always stable, preferred default.
    Linear,
    /// Catmull–Rom spline interpolation — smoother but requires ≥ 3 points;
    /// falls back to linear automatically when fewer points are available.
    Spline,
}

/// Central data structure for isoconversional analysis.
///
/// All experiments are resampled onto a common α grid so that every solver
/// can work with aligned ndarray columns instead of ragged per-experiment
/// vectors.  Shape of all 2-D arrays is `[n_experiments × n_alpha_points]`.
///
/// Build via [`ConversionGridBuilder`]; do not construct directly.
pub struct ConversionGrid {
    /// The shared α (conversion) grid, shape `[n_alpha]`.
    pub eta: Array1<f64>,
    /// Temperature at each `(experiment, α)` cell (K).
    pub temperature: Array2<f64>,
    /// `1/T` pre-cached to avoid repeated division in solver inner loops.
    pub inv_temperature: Array2<f64>,
    /// Time at which experiment `i` reached conversion `α[j]` (seconds).
    pub time: Array2<f64>,
    /// dα/dt at each `(experiment, α)` cell (s⁻¹).
    pub conversion_rate: Array2<f64>,
    /// Time-step increments Δt between consecutive α layers, required by the
    /// Vyazovkin method.  `None` unless built with `with_dt_matrix()`.
    pub dt: Option<Array2<f64>>,
    /// Per-experiment metadata (id, heating rate, isothermal temperature).
    pub meta: Vec<ExperimentMeta>,
}
impl ConversionGrid {
    /// Prints a diagnostic report to stdout: α range, per-experiment max-α
    /// time and temperature, zero-value warnings, and array dimensions.
    /// Useful for sanity-checking grids built from real or synthetic data.
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

/// Distinguishes isothermal from non-isothermal experiment series.
/// Used by `ConversionGridBuilder::build_universal` to select the correct
/// grid-building path.
pub enum ExperimentKind {
    /// All experiments run at a fixed temperature; heating rate is absent.
    Isothermal,
    /// Experiments run with a linear temperature ramp; heating rate is present.
    NonIsothermal,
}
/// Controls how the α range of the grid is determined.
pub enum AlphaRangeMode {
    /// Automatically compute the intersection of all experiments' conversion
    /// ranges so every experiment has data at every grid point.
    Auto,
    /// Use an explicit `[min, max]` range supplied by the caller.
    Manual { min: f64, max: f64 },
}
/// Optional extras computed during grid construction.
// delta t required for Vyazovkin
pub struct GridExtras {
    /// If `true`, the `dt` matrix (time-step increments between α layers) is
    /// computed and stored in `ConversionGrid::dt`.  Required by Vyazovkin.
    pub dt_matrix: bool,
}
impl Default for GridExtras {
    fn default() -> Self {
        Self { dt_matrix: false }
    }
}
/// Fluent builder for [`ConversionGrid`].
///
/// Typical usage:
/// ```rust,ignore
/// let grid = ConversionGridBuilder::new()
///     .eta_range(0.05, 0.95)
///     .segments(100)
///     .interpolation(GridInterpolation::Linear)
///     .build(&view)?;
/// ```
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
    /// Creates a builder with default settings: auto α range, 50 segments,
    /// linear interpolation, no safety margins, no dt matrix.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets an explicit `[min, max]` α range.
    pub fn eta_range(mut self, min: f64, max: f64) -> Self {
        self.eta_mode = AlphaRangeMode::Manual { min, max };
        self
    }

    /// Switches to automatic α range (intersection of all experiments).
    pub fn auto_range(mut self) -> Self {
        self.eta_mode = AlphaRangeMode::Auto;
        self
    }

    /// Sets the number of uniformly-spaced α grid points.
    pub fn segments(mut self, n: usize) -> Self {
        self.segments = n;
        self
    }

    /// Sets the interpolation strategy used when resampling experiments.
    pub fn interpolation(mut self, mode: GridInterpolation) -> Self {
        self.interpolation = mode;
        self
    }

    /// Adds safety margins that shrink the α range inward by `left` on the
    /// low end and `right` on the high end, avoiding edge artefacts.
    pub fn safety_margin(mut self, left: f64, right: f64) -> Self {
        self.margin_left = left;
        self.margin_right = right;
        self
    }

    /// Enables computation of the Δt matrix required by the Vyazovkin method.
    pub fn with_dt_matrix(mut self) -> Self {
        self.extras.dt_matrix = true;
        self
    }

    /// Conditionally enables the Δt matrix; equivalent to calling
    /// `with_dt_matrix()` when `compute` is `true`.
    pub fn compute_dt(self, compute: bool) -> Self {
        if compute {
            return self.with_dt_matrix();
        } else {
            return self;
        }
    }

    /// Computes the Δt matrix from a time array: `dt[i, j] = time[i, j+1] - time[i, j]`.
    /// The last column is filled by repeating the second-to-last value.
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

    /// Computes the intersection of all experiments' conversion ranges.
    ///
    /// Returns `(global_min, global_max)` where `global_min = max(local_mins)`
    /// and `global_max = min(local_maxs)`.  Errors if the intersection is empty.
    pub fn calc_auto_eta_range(
        &self,
        data: &KineticDataView,
    ) -> Result<(f64, f64), TGADomainError> {
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

    /// Allocating variant of the non-isothermal grid builder (Vec-based).
    /// Prefer `build_nonisothermal` for better memory locality.
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

    /// Builds a `ConversionGrid` for **non-isothermal** experiments.
    ///
    /// Allocates the output arrays once and fills them row-by-row using
    /// `interpolate_experiment_into`, avoiding intermediate `Vec` allocations.
    /// This is the standard path for OFW, KAS, Starink, Differential Friedman,
    /// and Vyazovkin.
    pub fn build_nonisothermal(
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

    /// Builds a uniformly-spaced α grid of `n` points in `[min, max)`.
    fn build_eta_grid(min: f64, max: f64, n: usize) -> Vec<f64> {
        let step = (max - min) / n as f64;

        (0..n).map(|i| min + step * i as f64).collect()
    }

    /// Dispatches to the correct interpolation function, returning
    /// `(time, temperature, rate)` vectors of length `eta_grid.len()`.
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

    /// In-place variant of `interpolate_experiment`: writes directly into
    /// pre-allocated ndarray row views, avoiding extra allocations.
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

    /// Builds a `ConversionGrid` for **isothermal** experiments.
    ///
    /// Only the time column is interpolated from the data; the temperature
    /// matrix is filled with the constant `isothermal_temperature` value from
    /// each experiment's metadata.  `conversion_rate` is left as zeros because
    /// isothermal methods (Integral Friedman) do not use it.
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

    /// Dispatches to `interpolate_time_linear_only` or
    /// `interpolate_time_spline_only` depending on `method`.
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

    /// Unified entry point that selects `build_nonisothermal` or
    /// `build_isothermal` based on `kind`.
    pub fn build_universal(
        self,
        data: &KineticDataView,
        kind: ExperimentKind,
    ) -> Result<ConversionGrid, TGADomainError> {
        match kind {
            ExperimentKind::NonIsothermal => self.build_nonisothermal(data),

            ExperimentKind::Isothermal => self.build_isothermal(data),
        }
    }
}

/// Linear interpolation of time, temperature, and rate onto `eta_grid`.
///
/// Allocating variant; used by `build_withcopy`.
pub fn interpolate_linear(
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

/// Catmull–Rom spline interpolation of time, temperature, and rate onto
/// `eta_grid`.  Falls back to linear when fewer than 4 source points exist.
///
/// Allocating variant; used by `build_withcopy`.
pub fn interpolate_spline(
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

/// In-place linear interpolation of time, temperature, and rate.
///
/// Uses a single forward pointer that advances monotonically, giving O(n + m)
/// complexity instead of the O(n·m) naïve search.  Requires the source
/// conversion vector to be monotonically non-decreasing.
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

/// In-place Catmull–Rom spline interpolation of time, temperature, and rate.
/// Requires ≥ 3 source points and a monotonically non-decreasing conversion vector.
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

/// In-place linear interpolation of the time column only.
/// Used by `build_isothermal` where temperature is taken from metadata.
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

/// In-place Catmull–Rom spline interpolation of the time column only.
/// Requires ≥ 3 source points.
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

/// Convenience enum for selecting a kinetic method by name in high-level APIs.
pub enum KineticMethodKind {
    OFW,
    KAS,
    Starink,
    Friedman,
}
