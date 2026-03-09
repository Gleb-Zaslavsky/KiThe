use crate::Kinetics::experimental_kinetics::LSQSplines::SolverKind;
use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    OneFramePlot, TGADataset, TGADomainError,
};
use crate::Kinetics::experimental_kinetics::splines::SplineKind;
use log::info;
use std::time::Instant;
#[derive(Clone, Debug)]
pub struct ViewRange {
    pub t_min: f64,
    pub t_max: f64,
}

#[derive(Clone, Debug)]
pub struct PlotSeries {
    pub name_x: String,
    pub name_y: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

/// Selects which axis column from `oneframeplot` is used.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XY {
    X,
    Y,
}
impl Default for XY {
    fn default() -> Self {
        XY::X
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Ranges {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}
impl Default for Ranges {
    fn default() -> Self {
        Self {
            x_min: 0.0,
            x_max: 0.0,
            y_min: 0.0,
            y_max: 0.0,
        }
    }
}
impl TGADataset {
    //==================================================================================
    // PLOT SETTTERS
    //=======================================================================
    /// Resolves and returns the configured `(x, y)` column names from `oneframeplot`.
    /// Получение имен колонок x и y из oneframeplot
    fn oneframeplot_xy_names(&self) -> Result<(String, String), TGADomainError> {
        let oneframe = self.oneframeplot.as_ref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function  oneframeplot_xy_names:oneframeplot is not configured".to_string(),
            )
        })?;

        let x_col = oneframe.x.clone().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function  oneframeplot_xy_names:oneframeplot.x is not configured".to_string(),
            )
        })?;
        let y_col = oneframe.y.clone().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function  oneframeplot_xy_names:oneframeplot.y is not configured".to_string(),
            )
        })?;

        Ok((x_col, y_col))
    }
    pub fn get_oneframeplot_x_name(&self) -> Result<String, TGADomainError> {
        self.oneframeplot_axis_name(XY::X)
    }
    pub fn get_oneframeplot_y_name(&self) -> Result<String, TGADomainError> {
        self.oneframeplot_axis_name(XY::Y)
    }
    /// Resolves and returns one configured axis column name from `oneframeplot`.
    /// Получение имени колонки для заданной оси из oneframeplot
    pub fn oneframeplot_axis_name(&self, axis: XY) -> Result<String, TGADomainError> {
        let (x, y) = self.oneframeplot_xy_names()?;
        Ok(match axis {
            XY::X => x,
            XY::Y => y,
        })
    }

    /// Sets the x-axis column name in `oneframeplot`.
    ///
    /// The provided column must exist in the dataset schema.
    /// If `oneframeplot` is not initialized yet, it is created.
    /// Установка имени колонки для оси x в oneframeplot
    pub fn set_oneframeplot_x(mut self, x_col: &str) -> Result<Self, TGADomainError> {
        println!("\n setting x");
        if !self.schema.columns.contains_key(x_col) {
            return Err(TGADomainError::ColumnNotFound(x_col.into()));
        }

        match self.oneframeplot.as_mut() {
            Some(plot) => plot.x = Some(x_col.to_string()),
            None => {
                self.oneframeplot = Some(OneFramePlot {
                    x: Some(x_col.to_string()),
                    y: None,
                });
            }
        }

        Ok(self)
    }

    /// Sets the y-axis column name in `oneframeplot`.
    ///
    /// The provided column must exist in the dataset schema.
    /// If `oneframeplot` is not initialized yet, it is created.
    /// Установка имени колонки для оси y в oneframeplot
    pub fn set_oneframeplot_y(mut self, y_col: &str) -> Result<Self, TGADomainError> {
        println!("\n setting y");
        if !self.schema.columns.contains_key(y_col) {
            return Err(TGADomainError::ColumnNotFound(y_col.into()));
        }

        match self.oneframeplot.as_mut() {
            Some(plot) => plot.y = Some(y_col.to_string()),
            None => {
                self.oneframeplot = Some(OneFramePlot {
                    x: None,
                    y: Some(y_col.to_string()),
                });
            }
        }

        Ok(self)
    }
    //==================================================================
    //  GETTERS
    //==================================================================
    /// Получение значений для заданной оси в виде вектора
    pub fn get_axis_as_vec(&self, axis: XY) -> Result<Vec<f64>, TGADomainError> {
        let selected = self.oneframeplot_axis_name(axis)?;
        println!("\n column: {} ", selected);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(&selected).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// returns full x column as Vec
    /// Получение значений оси x в виде вектора
    pub fn get_x_as_vec(&self) -> Result<Vec<f64>, TGADomainError> {
        let x = XY::X;
        let x_vec = self.get_axis_as_vec(x)?;
        Ok(x_vec)
    }
    /// returns full y column as Vec
    /// Получение значений оси y в виде вектора
    pub fn get_y_as_vec(&self) -> Result<Vec<f64>, TGADomainError> {
        let x_vec = self.get_axis_as_vec(XY::Y)?;
        Ok(x_vec)
    }
    /// returns struct
    ///  PlotSeries {
    /// pub name_x: String,
    /// pub name_y: String,
    ///  pub x: Vec<f64>,
    ///  pub y: Vec<f64>,
    /// Получение данных для построения графика
    pub fn get_plotseries(&self) -> Result<PlotSeries, TGADomainError> {
        let x = self.get_x_as_vec()?;
        let y = self.get_y_as_vec()?;
        let oneframe = self.oneframeplot.as_ref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function   get_plotseries:oneframeplot is not configured".to_string(),
            )
        })?;

        let x_col = oneframe.x.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function   get_plotseries:oneframeplot.x is not configured".to_string(),
            )
        })?;
        let y_col = oneframe.y.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function   get_plotseries:oneframeplot.y is not configured".to_string(),
            )
        })?;
        let series = PlotSeries {
            name_x: x_col.to_string(),
            name_y: y_col.to_string(),
            x,
            y,
        };
        Ok(series)
    }
    //==================================================================
    // PLOT DATA MANIPULATIONS
    //===================================================================
    /// Samples one plot series using x/y names stored in `oneframeplot`.
    ///
    /// This method applies the same range filtering and downsampling logic
    /// as `sample_columns`, but always returns a single `PlotSeries`.
    /// Выборка данных для одной серии графика
    pub fn sample_oneframeplot(
        &self,
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<PlotSeries, TGADomainError> {
        let oneframe = self.oneframeplot.as_ref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function sample_oneframeplot: oneframeplot is not configured".to_string(),
            )
        })?;

        let x_col = oneframe.x.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function sample_oneframeplot:oneframeplot.x is not configured".to_string(),
            )
        })?;
        let y_col = oneframe.y.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function sample_oneframeplot: oneframeplot.y is not configured".to_string(),
            )
        })?;

        let mut sampled = self.sample_columns(x_col, &[y_col], range, max_points)?;
        sampled.pop().ok_or_else(|| {
            TGADomainError::InvalidOperation("no data points in selected range".to_string())
        })
    }

    /// Spline-resamples the configured one-frame plot pair.
    ///
    /// Internally forwards:
    /// - `oneframeplot.x` as `time_col`
    /// - `oneframeplot.y` as `columns = &[y]`
    ///
    /// The resampled time column is written under `new_time_col`.
    /// Сплайн-ресэмплинг данных для графика
    pub fn spline_resample_oneframeplot(
        self,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
    ) -> Result<Self, TGADomainError> {
        let oneframe = self.oneframeplot.clone().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot is not configured".to_string(),
            )
        })?;

        let x_col = oneframe.x.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot.x is not configured"
                    .to_string(),
            )
        })?;
        let y_col = oneframe.y.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot.y is not configured"
                    .to_string(),
            )
        })?;

        self.spline_resample_columns(x_col, new_time_col, &[y_col], n_points, kind)
    }
    pub fn spline_resample_oneframeplot_as(
        self,
        new_time_col: &str,
        n_points: usize,
        kind: SplineKind,
        out_col: Option<&str>,
    ) -> Result<Self, TGADomainError> {
        let oneframe = self.oneframeplot.clone().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot is not configured".to_string(),
            )
        })?;

        let x_col = oneframe.x.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot.x is not configured"
                    .to_string(),
            )
        })?;
        let y_col = oneframe.y.as_deref().ok_or_else(|| {
            TGADomainError::InvalidOperation(
                "function spline_resample_oneframeplot: oneframeplot.y is not configured"
                    .to_string(),
            )
        })?;
        self.spline_resample_columns_as(x_col, new_time_col, &[y_col], &[out_col], n_points, kind)
    }
    /// Applies a generic operation to one selected axis column (`X` or `Y`)
    /// configured in `oneframeplot`.
    ///
    /// The closure receives `(dataset, selected_column_name)` and must return
    /// the updated dataset.
    /// Применение операции к одной оси графика
    pub fn with_x_or_y<F>(self, axis: XY, op: F) -> Result<Self, TGADomainError>
    where
        F: FnOnce(Self, &str) -> Result<Self, TGADomainError>,
    {
        let selected = self.oneframeplot_axis_name(axis)?;
        op(self, &selected)
    }

    /// Applies a generic operation that engages both configured axis columns
    /// from `oneframeplot`.
    ///
    /// The closure receives `(dataset, x_column_name, y_column_name)` and must
    /// return the updated dataset.
    /// Применение операции к обеим осям графика
    pub fn with_x_and_y<F>(self, op: F) -> Result<Self, TGADomainError>
    where
        F: FnOnce(Self, &str, &str) -> Result<Self, TGADomainError>,
    {
        let (x_col, y_col) = self.oneframeplot_xy_names()?;
        op(self, &x_col, &y_col)
    }

    /// Cuts all rows before `start_value` on the selected axis (`X` or `Y`).
    ///
    /// Because filtering is row-wise, points in the paired axis are removed
    /// at the same indices, so x/y lengths remain equal.
    /// Обрезка данных до заданного значения по оси
    pub fn cut_before_x_or_y(self, axis: XY, start_value: f64) -> Result<Self, TGADomainError> {
        self.with_x_or_y(axis, |ds, selected_col| {
            Ok(ds.trim_range(selected_col, start_value, f64::INFINITY))
        })
    }

    /// Cuts all rows after `end_value` on the selected axis (`X` or `Y`).
    ///
    /// Because filtering is row-wise, points in the paired axis are removed
    /// at the same indices, so x/y lengths remain equal.
    /// Обрезка данных после заданного значения по оси
    pub fn cut_after_x_or_y(self, axis: XY, end_value: f64) -> Result<Self, TGADomainError> {
        self.with_x_or_y(axis, |ds, selected_col| {
            Ok(ds.trim_range(selected_col, f64::NEG_INFINITY, end_value))
        })
    }

    /// Cuts rows to the closed interval `[from, to]` on the selected axis (`X` or `Y`).
    ///
    /// Because filtering is row-wise, points in the paired axis are removed
    /// at the same indices, so x/y lengths remain equal.
    /// Обрезка данных в заданном диапазоне по оси
    pub fn cut_range_x_or_y(self, axis: XY, from: f64, to: f64) -> Result<Self, TGADomainError> {
        self.with_x_or_y(axis, |ds, selected_col| {
            Ok(ds.trim_range(selected_col, from, to))
        })
    }

    /// Removes rows inside the closed interval `[from, to]` on the selected axis (`X` or `Y`).
    ///
    /// Because filtering is row-wise, points in the paired axis are removed
    /// at the same indices, so x/y lengths remain equal.
    pub fn cut_range_inverse_x_or_y(
        self,
        axis: XY,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        self.with_x_or_y(axis, |ds, selected_col| {
            Ok(ds.trim_range_inverse(selected_col, from, to))
        })
    }

    /// Returns the smallest Euclidean distance from `point` to the nearest
    /// data point formed by configured `oneframeplot` `(x, y)` columns.
    ///
    /// Rows where either x or y is null are ignored.
    /// Вычисление минимального расстояния до точки данных
    pub fn min_distance_to_oneframeplot_point(
        &self,
        point: (f64, f64),
    ) -> Result<f64, TGADomainError> {
        let (x_col, y_col) = self.oneframeplot_xy_names()?;
        println!("\n column: {} ", x_col);
        println!("\n column: {} ", y_col);
        let df = self.frame.clone().collect()?;

        let x = df.column(&x_col)?.f64()?;
        let y = df.column(&y_col)?.f64()?;

        let mut best: Option<f64> = None;
        for (ox, oy) in x.into_iter().zip(y.into_iter()) {
            if let (Some(px), Some(py)) = (ox, oy) {
                let dx = px - point.0;
                let dy = py - point.1;
                let d = (dx * dx + dy * dy).sqrt();
                best = Some(match best {
                    Some(cur) => cur.min(d),
                    None => d,
                });
            }
        }

        best.ok_or_else(|| {
            TGADomainError::InvalidOperation(format!(
                "No valid (x, y) points found in columns '{}' and '{}'",
                x_col, y_col
            ))
        })
    }

    /// Получение диапазонов значений для осей графика
    pub fn plot_xy_ranges(&self) -> Result<Ranges, TGADomainError> {
        let (x_col, y_col) = self.oneframeplot_xy_names()?;
        println!("\n column: {} ", x_col);
        println!("\n column: {} ", y_col);
        let df = self.frame.clone().collect()?;

        let x = df.column(&x_col)?.f64()?;
        let y = df.column(&y_col)?.f64()?;

        let (x_min, x_max) = x
            .into_iter()
            .flatten()
            .fold(None, |acc: Option<(f64, f64)>, v| {
                Some(match acc {
                    Some((mn, mx)) => (mn.min(v), mx.max(v)),
                    None => (v, v),
                })
            })
            .ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "x column '{}' has no valid values",
                    x_col
                ))
            })?;
        let (y_min, y_max) = y
            .into_iter()
            .flatten()
            .fold(None, |acc: Option<(f64, f64)>, v| {
                Some(match acc {
                    Some((mn, mx)) => (mn.min(v), mx.max(v)),
                    None => (v, v),
                })
            })
            .ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "y column '{}' has no valid values",
                    y_col
                ))
            })?;

        Ok(Ranges {
            x_min,
            x_max,
            y_min,
            y_max,
        })
    }

    pub fn offset_y_column_in_range_by_x_reference(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (x_col, y_col) = self.oneframeplot_xy_names()?;
        let r = self.offset_column_in_range_by_reference(&y_col, &x_col, offset, from, to);
        Ok(r)
    }
    pub fn offset_y_column_in_its_range(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (_, y_col) = self.oneframeplot_xy_names()?;
        let r = self.offset_column_in_its_range(&y_col, offset, from, to);
        Ok(r)
    }

    pub fn offset_x_column_in_its_range(
        self,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (x_col, _) = self.oneframeplot_xy_names()?;
        let r = self.offset_column_in_its_range(&x_col, offset, from, to);
        Ok(r)
    }
    pub fn scale_y_column_in_range_by_x_reference(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (x_col, y_col) = self.oneframeplot_xy_names()?;
        let r = self.scale_column_in_range_by_reference(&y_col, &x_col, s, from, to);
        Ok(r)
    }
    pub fn scale_y_column_in_its_range(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (_, y_col) = self.oneframeplot_xy_names()?;
        let r = self.scale_column_in_its_range(&y_col, s, from, to);
        Ok(r)
    }

    pub fn scale_x_column_in_its_range(
        self,
        s: f64,
        from: f64,
        to: f64,
    ) -> Result<Self, TGADomainError> {
        let (x_col, _) = self.oneframeplot_xy_names()?;
        let r = self.scale_column_in_its_range(&x_col, s, from, to);
        Ok(r)
    }
    //===================================================================================================================
    // MISC

    /// Samples a single column with evenly distributed points.
    ///
    /// # Arguments
    /// * `col_name` - Name of the column to sample
    /// * `range` - Optional range (min, max) to filter values. If None, samples entire column
    /// * `n_points` - Number of sample points to return
    ///
    /// # Returns
    /// Vector of f64 values with approximately evenly distributed sample points
    ///
    /// # Behavior
    /// 1. If `range` is None: returns evenly distributed points from the entire column
    /// 2. If `range` is Some((min, max)): finds points within the range and samples them evenly
    ///
    /// Выборка данных из одной колонки с равномерным распределением точек
    pub fn sample_column(
        &self,
        col_name: &str,
        range: Option<(f64, f64)>,
        n_points: usize,
    ) -> Result<Vec<f64>, TGADomainError> {
        let df = self.frame.clone().collect()?;
        let col = df.column(col_name)?.f64()?;

        let indices: Vec<usize> = col
            .into_iter()
            .enumerate()
            .filter_map(|(i, val)| {
                val.and_then(|v| {
                    if range.map_or(true, |(min, max)| v >= min && v <= max) {
                        Some(i)
                    } else {
                        None
                    }
                })
            })
            .collect();

        if indices.is_empty() {
            return Ok(vec![]);
        }

        let step = (indices.len() as f64 / n_points as f64).ceil().max(1.0) as usize;
        let sampled: Result<Vec<f64>, TGADomainError> = indices
            .into_iter()
            .step_by(step)
            .map(|i| {
                col.get(i).ok_or(TGADomainError::InvalidOperation(format!(
                    "no column number {}",
                    i
                )))
            })
            .collect();

        Ok(sampled?)
    }

    pub fn sample_columns(
        &self,
        time_col: &str,
        value_cols: &[&str],
        range: Option<ViewRange>,
        max_points: usize,
    ) -> Result<Vec<PlotSeries>, TGADomainError> {
        println!("\n column: {} ", time_col);
        for &col in value_cols {
            println!("\n column: {} ", col);
        }
        let df = self.frame.clone().collect()?;
        let row_count = df.height();

        let mut out = Vec::new();

        for &col in value_cols {
            let time = df.column(time_col)?.f64()?;
            let s = df.column(col)?.f64()?;

            let mut points: Vec<(f64, f64)> = (0..row_count)
                .filter_map(|i| {
                    let t = time.get(i)?;
                    if !range
                        .as_ref()
                        .map_or(true, |r| t >= r.t_min && t <= r.t_max)
                    {
                        return None;
                    }

                    let y = s.get(i)?;
                    Some((t, y))
                })
                .collect();

            if points.is_empty() {
                continue;
            }

            let step = (points.len() as f64 / max_points as f64).ceil().max(1.0) as usize;
            points = points.into_iter().step_by(step).collect();
            let (x, y): (Vec<f64>, Vec<f64>) = points.into_iter().unzip();

            out.push(PlotSeries {
                name_x: time_col.to_string(),
                name_y: col.to_string(),
                x,
                y,
            });
        }

        Ok(out)
    }
    pub fn get_time_col(&self) -> Result<String, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?;

        Ok(time_col.to_string())
    }

    /// Получение имени колонки массы
    pub fn get_mass_col(&self) -> Result<String, TGADomainError> {
        let mass_col = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?;

        Ok(mass_col.to_string())
    }

    /// Получение имени колонки температуры
    pub fn get_temperature_col(&self) -> Result<String, TGADomainError> {
        let temp_col = self
            .schema
            .temperature
            .as_ref()
            .ok_or(TGADomainError::TemperatureNotBound)?;

        Ok(temp_col.to_string())
    }

    //===========================================================================
    //GOLDEN PIPELINE
    //============================================================================
    // the following methods are intended to be used in a chained manner for convenient data manipulation

    pub fn apply_golden_pipeline(
        self,
        config: GoldenPipelineConfig,
    ) -> Result<(Self, Vec<String>), TGADomainError> {
        let vector_of_non_monotonic_points = self.monotony_of_time_check()?;
        if vector_of_non_monotonic_points.len() > 0 {
            return Err(TGADomainError::TimeNonMonotonic(format!(
                "Non-monotonic points found at indices: {:?}",
                vector_of_non_monotonic_points
            )));
        }
        let k = config.k;
        let b = config.b;
        let time_cut_before = config.time_cut_before;
        let time_cut_after = config.time_cut_after;
        let averaging_time = config.averaging_time;
        // sav gol parameters
        let window_size = config.sav_gol_config.window_size;
        let poly_degree = config.sav_gol_config.poly_degree;
        let deriv = config.sav_gol_config.deriv;
        let delta = config.sav_gol_config.delta;
        // spines params
        let n_points = config.spline_config.n_points;
        let n_internal_points = config.spline_config.n_internal_points;
        let degree = config.spline_config.degree;

        let initial_list_columns = self.list_of_columns();

        let time_col = self.get_time_col()?;
        let temp_col = self.get_temperature_col()?;
        let mass_col = self.get_mass_col()?;

        let ds = self.calibrate_mass(k, b, &mass_col)?; // calibration
        // trim startup
        let ds = if let Some(time_cut_before) = time_cut_before {
            ds.cut_before_time(time_cut_before)
        } else {
            ds
        };
        let ds = if let Some(time_cut_after) = time_cut_after {
            ds.cut_after_time(time_cut_after)
        } else {
            ds
        };

        let ds = ds
            .celsius_to_kelvin() // units
            .seconds_to_hours()
            .move_time_to_zero()? // time moved to zero
            .sg_filter_column(&mass_col, window_size, poly_degree, deriv, delta)? // smoothing
            .dimensionless_mass(0.0, averaging_time / 3600.0, "alpha")?
            .conversion(0.0, averaging_time, "eta")?
            .lsq_spline_resample_columns_as(
                &time_col,
                "time_splined",
                &[&mass_col, &temp_col, "alpha", "eta"],
                &[None; 4],
                n_points,
                degree,
                n_internal_points,
                SolverKind::Banded,
            )?
            //.trim_null_edges()
            .derive_temperature_rate("dT/dt")?
            .derive_dalpha_dt("dalpha_dt")?
            .derive_deta_dt("deta_dt")?;
        //.trim_null_edges();
        let modifed_list_of_columns = ds.list_of_columns();
        let start = Instant::now();
        let new_columns = modifed_list_of_columns
            .into_iter()
            .filter(|col| !initial_list_columns.contains(col))
            .collect::<Vec<_>>();
        info!("collcted new columns in {} ms", start.elapsed().as_millis());
        println!("\n new columns added by golden pipeline: {:?}", new_columns);
        Ok((ds, new_columns))
    }
}
#[derive(Debug, Clone)]
pub struct SavGolConfig {
    pub window_size: usize,
    pub poly_degree: usize,
    pub deriv: usize,
    pub delta: f64,
}

impl Default for SavGolConfig {
    fn default() -> Self {
        Self {
            window_size: 11,
            poly_degree: 3,
            deriv: 0,
            delta: 1.0,
        }
    }
}
#[derive(Debug, Clone)]
pub struct SplineConfig {
    pub n_points: usize,
    pub n_internal_points: usize,
    pub degree: usize,
}
impl Default for SplineConfig {
    fn default() -> Self {
        Self {
            n_points: 300,
            n_internal_points: 10_000,
            degree: 3,
        }
    }
}
#[derive(Debug, Clone)]
pub struct GoldenPipelineConfig {
    pub k: f64,
    pub b: f64,
    pub time_cut_before: Option<f64>,
    pub time_cut_after: Option<f64>,
    pub sav_gol_config: SavGolConfig,
    pub averaging_time: f64,
    pub spline_config: SplineConfig,
    pub save_to_new_experiment: bool,
    pub del_old_experiment: bool,
}
impl Default for GoldenPipelineConfig {
    fn default() -> Self {
        Self {
            k: 1.0,
            b: 0.0,
            time_cut_before: None,
            time_cut_after: None,
            sav_gol_config: SavGolConfig::default(),
            averaging_time: 1000.0,
            spline_config: SplineConfig::default(),
            save_to_new_experiment: true,
            del_old_experiment: true,
        }
    }
}
