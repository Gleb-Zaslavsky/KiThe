/*
Basic Transformations:

✅ with_column_expr - logs new column addition

✅ filter_rows - logs row filtering (irreversible)

✅ cut_interval - logs interval cutting (irreversible)

✅ trim_edges - logs edge trimming (irreversible)

✅ trim_null_edges - logs null edge trimming (irreversible)

✅ trim_null_edges_for_columns - logs selective null trimming (irreversible)

✅ rename_column - logs column renaming (reversible)

✅ drop_column - logs column deletion (irreversible)

Unit Transformations:

✅ celsius_to_kelvin - logs temperature unit conversion (reversible)

✅ seconds_to_hours - logs time unit conversion (reversible)

Algebraic Transformations:

✅ scale_columns - logs multi-column scaling (reversible)

✅ scale_column - logs single column scaling (reversible)

✅ scale_column_in_range_by_reference - logs conditional scaling (reversible)

✅ offset_column - logs column offset (reversible)

✅ offset_column_in_range_by_reference - logs conditional offset (reversible)

Calibration:

✅ calibrate_mass_from_voltage - logs mass calibration (reversible)

✅ calibrate_mass - logs mass calibration to new column (reversible)

✅ unary_column_op - logs generic unary operations (reversible)

Mathematical Operations:

✅ exp_column - logs exponential transformation (reversible)

✅ ln_column - logs logarithmic transformation (reversible)

Dimensionless Transformations:

✅ derive_dimensionless_mass - logs alpha calculation (reversible)

✅ dimensionless_mass - logs alpha calculation to custom column (reversible)

✅ derive_conversion - logs eta calculation (reversible)

✅ conversion - logs eta calculation to custom column (reversible)

*/

use crate::Kinetics::experimental_kinetics::one_experiment_dataset::{
    AffectedColumns, ColumnMeta, ColumnNature, ColumnOrigin, ColumnTypes, TGADataset,
    TGADomainError, UnaryOp, Unit,
};
use log::info;
use polars::error::PolarsResult;
use polars::prelude::DataType;
use polars::prelude::Expr;
use polars::prelude::*;
use std::time::Instant;
use tabled::builder::Builder;
use tabled::settings::Style;

impl TGADataset {
    //================================================================
    // GETTERS

    /// Получение значений времени
    pub fn get_time(&self) -> Result<Vec<f64>, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?;
        println!("\n column: {} ", time_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(time_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }

    /// Получение значений массы
    pub fn get_mass(&self) -> Result<Vec<f64>, TGADomainError> {
        let mass_col = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?;
        println!("\n column: {} ", mass_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(mass_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }

    /// Получение значений температуры
    pub fn get_temperature(&self) -> Result<Vec<f64>, TGADomainError> {
        let temp_col = self
            .schema
            .temperature
            .as_ref()
            .ok_or(TGADomainError::TemperatureNotBound)?;
        println!("\n column: {} ", temp_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(temp_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }

    /// Получение значений указанной колонки
    pub fn get_column(&self, name: &str) -> Result<Vec<f64>, TGADomainError> {
        println!("\n column: {} ", name);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(name).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// Получение значений производной температуры
    pub fn get_dT_dt(&self) -> Result<Vec<f64>, TGADomainError> {
        let dT_dt_col = self
            .schema
            .dT_dt
            .as_ref()
            .ok_or(TGADomainError::DtDtNotFound)?;
        println!("\n column: {} ", dT_dt_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(dT_dt_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// Получение значений производной массы
    pub fn get_dm_dt(&self) -> Result<Vec<f64>, TGADomainError> {
        let dm_dt_col = self
            .schema
            .dm_dt
            .as_ref()
            .ok_or(TGADomainError::DmDtNotFound)?;
        println!("\n column: {} ", dm_dt_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(dm_dt_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// Получение значений безразмерной массы
    pub fn get_alpha(&self) -> Result<Vec<f64>, TGADomainError> {
        let alpha_col = self
            .schema
            .alpha
            .as_ref()
            .ok_or(TGADomainError::AlphaNotFound)?;
        println!("\n column: {} ", alpha_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(alpha_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// Получение значений производной безразмерной массы
    pub fn get_dalpha_dt(&self) -> Result<Vec<f64>, TGADomainError> {
        let dalpha_dt_col = self
            .schema
            .dm_dt
            .as_ref()
            .ok_or(TGADomainError::DalphaDtNotFound)?;
        println!("\n column: {} ", dalpha_dt_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df
            .column(dalpha_dt_col)
            .map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }

    /// Получение значений степени превращения
    pub fn get_eta(&self) -> Result<Vec<f64>, TGADomainError> {
        let eta_col = self
            .schema
            .eta
            .as_ref()
            .ok_or(TGADomainError::EtaNotFound)?;
        println!("\n column: {} ", eta_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df.column(eta_col).map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }
    /// Получение значений производной степени превращения
    pub fn get_deta_dt(&self) -> Result<Vec<f64>, TGADomainError> {
        let deta_dt_col = self
            .schema
            .dm_dt
            .as_ref()
            .ok_or(TGADomainError::DetaDtNotFound)?;
        println!("\n column: {} ", deta_dt_col);
        let df = self
            .frame
            .clone()
            .collect()
            .map_err(TGADomainError::PolarsError)?;
        let series = df
            .column(deta_dt_col)
            .map_err(TGADomainError::PolarsError)?;
        let f64_series = series.f64().map_err(TGADomainError::PolarsError)?;
        Ok(f64_series.into_no_null_iter().collect())
    }

    //================================================================
    // BASIC TRANSFORMATIONS: COLUMN CUT AND FILTER
    /// Добавление новой колонки с заданным выражением
    pub fn with_column_expr(mut self, meta: ColumnMeta, expr: Expr) -> Self {
        println!("\n column: {} ", meta.name);
        let col_name = meta.name.clone();
        self.frame = self.frame.with_column(expr.clone().alias(&meta.name));
        self.schema.columns.insert(meta.name.clone(), meta);

        self.log_operation(
            "with_column_expr",
            AffectedColumns::Specific(vec![col_name.clone()]),
            Some(expr),
            format!("Added column {}", col_name),
            true,
        );

        self
    }
    /// Works on all columns at once
    ///✔ Zero copy
    /// ✔ Lazy
    /// Фильтрация строк по заданному предикату
    pub fn filter_rows(mut self, predicate: Expr) -> Self {
        let pred_clone = predicate.clone();
        let frame = self.frame.clone().filter(predicate);

        self.log_operation(
            "filter_rows",
            AffectedColumns::All,
            Some(pred_clone),
            "Filtered rows by predicate".to_string(),
            false,
        );

        self.frame = frame;
        self
    }
    /// Фильтрация строк по маске
    pub fn filter_by_mask(self, mask: Expr) -> Self {
        self.filter_rows(mask)
    }
    /// Handy helper: cut interval
    /// Обрезка данных вне заданного интервала
    pub fn cut_interval(mut self, colmn: &str, from: f64, to: f64) -> Self {
        println!("\n column: {} ", colmn);
        let pred = (col(colmn).gt(lit(from)).and(col(colmn).lt(lit(to)))).not();
        self = self.filter_rows(pred);

        self.log_operation(
            "cut_interval",
            AffectedColumns::Specific(vec![colmn.to_string()]),
            None,
            format!("Cut interval ({}, {}) from column {}", from, to, colmn),
            false,
        );

        self
    }
    /// Обрезка данных по времени вне заданного интервала
    pub fn cut_time_interval(self, from: f64, to: f64) -> Self {
        let time_col = self.schema.time.as_ref().unwrap().clone();
        println!("\n column: {} ", time_col);
        self.cut_interval(&time_col, from, to)
    }
    /// Обрезка данных по температуре вне заданного интервала
    pub fn cut_temperature_interval(self, from: f64, to: f64) -> Self {
        let temp_col = self.schema.temperature.as_ref().unwrap().clone();
        println!("\n column: {} ", temp_col);
        self.cut_interval(&temp_col, from, to)
    }
    /// Обрезка данных по массе вне заданного интервала
    pub fn cut_mass_interval(self, from: f64, to: f64) -> Self {
        let mass_col = self.schema.mass.as_ref().unwrap().clone();
        println!("\n column: {} ", mass_col);
        self.cut_interval(&mass_col, from, to)
    }
    /// Обрезка данных до времени t_end
    /// Обрезка данных до времени t_start
    pub fn cut_before_time(self, t_start: f64) -> Self {
        let time_col = self.schema.time.as_ref().unwrap().clone();
        println!("\n column: {} ", time_col);

        self.filter_rows(col(&time_col).gt_eq(lit(t_start)))
    }

    pub fn cut_after_time(self, t_after: f64) -> Self {
        let time_col = self.schema.time.as_ref().unwrap().clone();
        println!("\n column: {} ", time_col);

        self.filter_rows(col(&time_col).lt_eq(lit(t_after)))
    }
    /// trim edges
    /// Обрезка данных в заданном диапазоне
    pub fn trim_range(self, column: &str, from: f64, to: f64) -> Self {
        println!("\n column: {} ", column);
        self.filter_rows(
            // from ≤ x ≤ to
            col(column).gt_eq(lit(from)).and(col(column).lt_eq(lit(to))),
        )
    }
    /// Removes rows inside the closed interval `[from, to]` and keeps the rest.
    pub fn trim_range_inverse(self, column: &str, from: f64, to: f64) -> Self {
        println!("\n column: {} ", column);
        self.filter_rows(col(column).lt(lit(from)).or(col(column).gt(lit(to))))
    }
    // When differentiation is used the last and the first
    // points of the rate column is nill
    /// Обрезка краев данных
    pub fn trim_edges(mut self, left: usize, right: usize) -> Self {
        let df = self.frame.clone().collect().unwrap();
        let total = df.height();
        let length = total.saturating_sub(left + right);
        let sliced_df = df.slice(left as i64, length);
        let frame = sliced_df.lazy();

        self.log_operation(
            "trim_edges",
            AffectedColumns::All,
            None,
            format!("Trimmed {} rows from left, {} from right", left, right),
            false,
        );

        self.frame = frame;
        self.oneframeplot = None;
        self
    }

    /// Обрезка краев с нулевыми значениями
    pub fn trim_null_edges(mut self) -> Self {
        let start = Instant::now();
        let df = self.frame.clone().collect().unwrap();

        let mut left = 0usize;
        let mut right = 0usize;

        for col in df.columns() {
            let s = col;
            let n = s.len();
            let mut l = 0;
            while l < n && s.is_null().get(l).unwrap_or(false) {
                l += 1;
            }
            let mut r = 0;
            while r < n && s.is_null().get(n - 1 - r).unwrap_or(false) {
                r += 1;
            }
            left = left.max(l);
            right = right.max(r);
        }

        self.frame = self
            .frame
            .slice(left as i64, (df.height() - left - right) as u32);

        self.log_operation(
            "trim_null_edges",
            AffectedColumns::All,
            None,
            format!("Trimmed {} null rows from left, {} from right", left, right),
            false,
        );
        info!("null trim finished in {} ms", start.elapsed().as_millis());
        self
    }

    /// Обрезка краев с нулевыми значениями для указанных колонок
    pub fn trim_null_edges_for_columns(mut self, cols_to_trim: &[&str]) -> Self {
        let df = self.frame.clone().collect().unwrap();

        let mut left = 0usize;
        let mut right = 0usize;

        for &col_name in cols_to_trim {
            if let Ok(col) = df.column(col_name) {
                let s = col;
                let n = s.len();
                let mut l = 0;
                while l < n && s.is_null().get(l).unwrap_or(false) {
                    l += 1;
                }
                let mut r = 0;
                while r < n && s.is_null().get(n - 1 - r).unwrap_or(false) {
                    r += 1;
                }
                left = left.max(l);
                right = right.max(r);
            }
        }

        self.frame = self
            .frame
            .slice(left as i64, (df.height() - left - right) as u32);

        self.log_operation(
            "trim_null_edges_for_columns",
            AffectedColumns::Specific(cols_to_trim.iter().map(|s| s.to_string()).collect()),
            None,
            format!(
                "Trimmed {} null rows from left, {} from right for {} columns",
                left,
                right,
                cols_to_trim.len()
            ),
            false,
        );

        self
    }
    /// Переименование колонки
    pub fn rename_column(mut self, old: &str, new: &str) -> Result<Self, TGADomainError> {
        println!("\n column: {} ", old);
        println!("\n column: {} ", new);
        self.frame = self.frame.rename([old], [new], false);

        if self.schema.time == Some(old.to_string()) {
            self.schema.time = Some(new.to_string());
        }
        if self.schema.temperature == Some(old.to_string()) {
            self.schema.temperature = Some(new.to_string());
        }
        if self.schema.mass == Some(old.to_string()) {
            self.schema.mass = Some(new.to_string());
        }

        if let Some(meta) = self.schema.columns.remove(old) {
            let mut meta = meta;
            meta.name = new.into();
            self.schema.columns.insert(new.into(), meta);
        }

        self.log_operation(
            "rename_column",
            AffectedColumns::Specific(vec![new.to_string()]),
            None,
            format!("Renamed column {} to {}", old, new),
            true,
        );

        Ok(self)
    }

    /// Drop a column by name from the dataset. This removes the column from the frame and metadata.
    /// If the column is currently bound as time, temperature, or mass, the corresponding binding is set to None.
    /// Удаление колонки
    pub fn drop_column(mut self, col_name: &str) -> Result<Self, TGADomainError> {
        println!("\n column: {} ", col_name);
        let exists_in_meta = self.schema.columns.contains_key(col_name);
        let exists_in_frame = self
            .frame
            .clone()
            .collect_schema()?
            .iter_names()
            .any(|n| n == col_name);
        if !exists_in_meta && !exists_in_frame {
            return Err(TGADomainError::ColumnNotFound(col_name.into()));
        }

        let schema = self.frame.clone().collect_schema()?;
        let remaining_cols: Vec<Expr> = schema
            .iter_names()
            .filter(|&name| name != col_name)
            .map(|name| col(name.as_str()))
            .collect();
        self.frame = self.frame.select(remaining_cols);

        self.schema.columns.remove(col_name);

        if self.schema.time.as_deref() == Some(col_name) {
            self.schema.time = None;
        }
        if self.schema.temperature.as_deref() == Some(col_name) {
            self.schema.temperature = None;
        }
        if self.schema.mass.as_deref() == Some(col_name) {
            self.schema.mass = None;
        }

        self.log_operation(
            "drop_column",
            AffectedColumns::Specific(vec![col_name.to_string()]),
            None,
            format!("Dropped column {}", col_name),
            false,
        );

        Ok(self)
    }
    //======================================================================
    // UNIT TRANSFORMATIONS
    /// Это не просто offset, это семантическая операция, меняющая единицы.
    pub fn celsius_to_kelvin(mut self) -> Self {
        let col_name = self.schema.temperature.as_ref().unwrap().clone();
        println!("\n column: {} ", col_name);
        let expr = (col(&col_name) + lit(273.15)).alias(&col_name);
        self.frame = self.frame.with_column(expr.clone());

        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Kelvin;
            meta.origin = ColumnOrigin::PolarsDerived;
        }

        self.log_operation(
            "celsius_to_kelvin",
            AffectedColumns::Semantic(vec![ColumnTypes::Temperature]),
            Some(expr),
            format!("Converted {} from Celsius to Kelvin", col_name),
            true,
        );

        self
    }
    /// Делим время на 3600 и меняем единицы на Hour
    pub fn seconds_to_hours(mut self) -> Self {
        let col_name = self.schema.time.as_ref().unwrap().clone();
        println!("\n column: {} ", col_name);
        let expr = (col(&col_name) / lit(3600.0)).alias(&col_name);
        self.frame = self.frame.with_column(expr.clone());

        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Hour;
            meta.origin = ColumnOrigin::PolarsDerived;
        }

        self.log_operation(
            "seconds_to_hours",
            AffectedColumns::Semantic(vec![ColumnTypes::Time]),
            Some(expr),
            format!("Converted {} from seconds to hours", col_name),
            true,
        );

        self
    }
    //=======================================================================
    // ARBITRARY ALGEBRAIC TRANSFORMATIONS ON COLUMN DATA
    /// Scale specified columns by factor
    /// Unit handling is not implemented yet
    pub fn scale_columns(mut self, cols: &[&str], factor: f64) -> Self {
        for &c in cols {
            println!("\n column: {} ", c);
        }
        let exprs: Vec<Expr> = cols
            .iter()
            .map(|&c| (col(c) * lit(factor)).alias(c))
            .collect();

        self.frame = self.frame.with_columns(exprs);

        for &c in cols {
            if let Some(meta) = self.schema.columns.get_mut(c) {
                meta.origin = ColumnOrigin::PolarsDerived;
            }
        }

        self.log_operation(
            "scale_columns",
            AffectedColumns::Specific(cols.iter().map(|s| s.to_string()).collect()),
            None,
            format!("Scaled {} columns by factor {}", cols.len(), factor),
            true,
        );

        self
    }
    pub fn scale_column(mut self, colmn: &str, factor: f64) -> Self {
        let expr = (col(colmn) * lit(factor)).alias(colmn);
        self.frame = self.frame.with_column(expr.clone());

        self.log_operation(
            "scale_column",
            AffectedColumns::Specific(vec![colmn.to_string()]),
            Some(expr),
            format!("Scaled column {} by factor {}", colmn, factor),
            true,
        );

        self
    }

    /// Scale values in `target_col` only where `reference_col` is inside closed interval `[from, to]`.
    pub fn scale_column_in_range_by_reference(
        mut self,
        target_col: &str,
        reference_col: &str,
        factor: f64,
        from: f64,
        to: f64,
    ) -> Self {
        println!("\n column: {} ", target_col);
        println!("\n column: {} ", reference_col);
        let expr = when(
            col(reference_col)
                .gt_eq(lit(from))
                .and(col(reference_col).lt_eq(lit(to))),
        )
        .then(col(target_col) * lit(factor))
        .otherwise(col(target_col))
        .alias(target_col);
        self.frame = self.frame.with_column(expr.clone());

        self.log_operation(
            "scale_column_in_range_by_reference",
            AffectedColumns::Specific(vec![target_col.to_string()]),
            Some(expr),
            format!(
                "Scaled {} by {} in range [{}, {}] of {}",
                target_col, factor, from, to, reference_col
            ),
            true,
        );

        self
    }

    /// Scale values in `colmn` only where the same column is inside closed interval `[from, to]`.
    pub fn scale_column_in_its_range(self, colmn: &str, factor: f64, from: f64, to: f64) -> Self {
        self.scale_column_in_range_by_reference(colmn, colmn, factor, from, to)
    }

    //
    /// Сдвиг времени к нулю
    pub fn move_time_to_zero(mut self) -> Result<Self, TGADomainError> {
        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        // Materialize the very first time value t0
        let df = self.frame.clone().select([col(&time)]).limit(1).collect()?;
        let t0 = df.column(&time)?.f64()?.get(0).ok_or_else(|| {
            TGADomainError::InvalidOperation("Time column is empty; cannot offset".into())
        })?;

        // Offset the time column by -t0
        self = self.offset_column(&time, -t0);
        Ok(self)
    }
    /// прибавление константы ко всем значениям в колонке
    /// Смещение значений в колонке на константу
    pub fn offset_column(mut self, colmn: &str, offset: f64) -> Self {
        println!("\n column: {} ", colmn);
        let expr = (col(colmn) + lit(offset)).alias(colmn);
        self.frame = self.frame.with_column(expr.clone());

        self.log_operation(
            "offset_column",
            AffectedColumns::Specific(vec![colmn.to_string()]),
            Some(expr),
            format!("Offset column {} by {}", colmn, offset),
            true,
        );

        self
    }

    /// Offset values in `target_col` only where `reference_col` is inside closed interval `[from, to]`.
    pub fn offset_column_in_range_by_reference(
        mut self,
        target_col: &str,
        reference_col: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Self {
        println!("\n column: {} ", target_col);
        println!("\n column: {} ", reference_col);
        let expr = when(
            col(reference_col)
                .gt_eq(lit(from))
                .and(col(reference_col).lt_eq(lit(to))),
        )
        .then(col(target_col) + lit(offset))
        .otherwise(col(target_col))
        .alias(target_col);
        self.frame = self.frame.with_column(expr.clone());

        self.log_operation(
            "offset_column_in_range_by_reference",
            AffectedColumns::Specific(vec![target_col.to_string()]),
            Some(expr),
            format!(
                "Offset {} by {} in range [{}, {}] of {}",
                target_col, offset, from, to, reference_col
            ),
            true,
        );

        self
    }

    /// Offset values in `colmn` only where the same column is inside closed interval `[from, to]`.
    pub fn offset_column_in_its_range(self, colmn: &str, offset: f64, from: f64, to: f64) -> Self {
        self.offset_column_in_range_by_reference(colmn, colmn, offset, from, to)
    }

    /// калибровочная прямая, от милливольт прибора к мг массы
    pub fn calibrate_mass_from_voltage(mut self, k: f64, b: f64) -> Self {
        let col_name = self.schema.mass.as_ref().unwrap().clone();
        println!("\n column: {} ", col_name);

        let expr = (col(&col_name) * lit(k) + lit(b)).alias(&col_name);
        self.frame = self.frame.with_column(expr.clone());

        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Milligram;
            meta.origin = ColumnOrigin::PolarsDerived;
        }

        self.log_operation(
            "calibrate_mass_from_voltage",
            AffectedColumns::Semantic(vec![ColumnTypes::Mass]),
            Some(expr),
            format!("Calibrated mass from voltage: {}*x + {}", k, b),
            true,
        );

        self
    }

    /// Калибровка массы с заданными параметрами
    pub fn calibrate_mass(mut self, k: f64, b: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let src = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
        println!("\n column: {} ", src);
        let k = k.clone();
        let b = b.clone();
        let func = Box::new(move |u| k.clone() * u + b.clone());
        self = self.unary_column_op(
            &src,
            new_col,
            UnaryOp {
                func: func,
                output_unit: Unit::Milligram,
                domain_check: None,
            },
        )?;
        self.schema.mass = Some(new_col.into());

        self.log_operation(
            "calibrate_mass",
            AffectedColumns::Specific(vec![new_col.to_string()]),
            None,
            format!("Calibrated mass to new column {}: {}*x + {}", new_col, k, b),
            true,
        );

        Ok(self)
    }

    /// Унарная операция над колонкой
    pub fn unary_column_op(
        mut self,
        src: &str,
        dst: &str,
        op: UnaryOp,
    ) -> Result<Self, TGADomainError> {
        println!("\n column: {} ", src);
        if let Some(check) = op.domain_check {
            check(&self, src)?;
        }

        let dt = |_schema: &Schema, field: &Field| {
            Ok(Field::new(field.name().clone(), DataType::Float64))
        };

        let out_unit = op.output_unit;
        let func = op.func;

        self.frame = self.frame.with_column(
            col(src)
                .map(
                    move |s| {
                        let f = &func;
                        let out_series = s.f64()?.apply_values(|x| f(x)).into_series();
                        Ok(out_series.into())
                    },
                    dt,
                )
                .alias(dst),
        );

        self.schema.columns.insert(
            dst.to_string(),
            ColumnMeta {
                name: dst.to_string(),
                unit: out_unit,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Unknown,
            },
        );

        self.log_operation(
            "unary_column_op",
            AffectedColumns::Specific(vec![dst.to_string()]),
            None,
            format!(
                "Applied unary operation: {} -> {} (unit: {:?})",
                src, dst, out_unit
            ),
            true,
        );

        Ok(self)
    }

    /// Экспонента колонки
    pub fn exp_column(mut self, col_name: &str) -> Result<Self, TGADomainError> {
        println!("\n column: {} ", col_name);
        let df = self
            .frame
            .clone()
            .select([col(col_name).min().alias("min")])
            .collect()?;

        let min_val = df.column("min").unwrap().f64().unwrap().get(0).unwrap();

        if min_val <= 0.0 {
            return Err(TGADomainError::NonPositiveValues {
                column: col_name.to_string(),
            });
        }
        let dt = |_schema: &Schema, field: &Field| {
            Ok(Field::new(field.name().clone(), DataType::Float64))
        };
        let new_name = format!("exp_{}", col_name);
        self.frame = self.frame.with_column(
            col(col_name)
                .map(
                    |s| Ok(s.f64()?.apply_values(|x| x.exp()).into_series().into()),
                    dt,
                )
                .alias(&new_name),
        );

        self.schema.columns.insert(
            new_name.clone(),
            ColumnMeta {
                name: new_name.clone(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Unknown,
            },
        );

        self.log_operation(
            "exp_column",
            AffectedColumns::Specific(vec![new_name.clone()]),
            None,
            format!("Computed exp({}) -> {}", col_name, new_name),
            true,
        );

        Ok(self)
    }

    /// Логарифмическая операция
    pub fn ln_op() -> Result<UnaryOp, TGADomainError> {
        let uo = UnaryOp {
            func: Box::new(|x| x.ln()),
            output_unit: Unit::Dimensionless,
            domain_check: Some(|ds, column| {
                let min = ds
                    .frame
                    .clone()
                    .select([col(column).min()])
                    .collect()?
                    .column(column)?
                    .f64()?
                    .get(0)
                    .unwrap();

                if min <= 0.0 {
                    Err(TGADomainError::NonPositiveValues {
                        column: column.to_string(),
                    })
                } else {
                    Ok(())
                }
            }),
        };
        Ok(uo)
    }
    /// логарифмирование колонок с проверкой на неотрицательность
    pub fn ln_column(mut self, col_name: &str) -> Result<Self, TGADomainError> {
        println!("\n column: {} ", col_name);
        let df = self
            .frame
            .clone()
            .select([col(col_name).min().alias("min")])
            .collect()
            .unwrap();

        let min_val = df.column("min").unwrap().f64().unwrap().get(0).unwrap();

        if min_val <= 0.0 {
            return Err(TGADomainError::NonPositiveValues {
                column: col_name.to_string(),
            });
        }
        let dt = |_schema: &Schema, field: &Field| {
            Ok(Field::new(field.name().clone(), DataType::Float64))
        };
        let new_name = format!("ln_{}", col_name);
        self.frame = self.frame.with_column(
            col(col_name)
                .map(
                    |s| Ok(s.f64()?.apply_values(|x| x.ln()).into_series().into()),
                    dt,
                )
                .alias(&new_name),
        );

        self.schema.columns.insert(
            new_name.clone(),
            ColumnMeta {
                name: new_name.clone(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Unknown,
            },
        );

        self.log_operation(
            "ln_column",
            AffectedColumns::Specific(vec![new_name.clone()]),
            None,
            format!("Computed ln({}) -> {}", col_name, new_name),
            true,
        );

        Ok(self)
    }
    /// среднее по интервалу
    /// Вычисление среднего значения по интервалу
    pub fn mean_on_interval(
        &self,
        value_col: &str,
        time_col: &str,
        from: f64,
        to: f64,
    ) -> PolarsResult<f64> {
        println!("\n column: {} ", value_col);
        println!("\n column: {} ", time_col);
        let df = self
            .frame
            .clone()
            .filter(
                col(time_col)
                    .gt_eq(lit(from))
                    .and(col(time_col).lt_eq(lit(to))),
            )
            .select([col(value_col).mean()])
            .collect()?;

        Ok(df.column(value_col)?.f64()?.get(0).unwrap())
    }

    /// Compute mean of `col` restricted to its own value range `[from, to]`.
    /// This is just a thin wrapper around `mean_on_interval` using the same
    /// column for values and the range filter.
    pub fn mean_on_interval_on_own_range(
        &self,
        colmn: &str,
        from: f64,
        to: f64,
    ) -> PolarsResult<f64> {
        let df = self
            .frame
            .clone()
            .filter(col(colmn).gt_eq(lit(from)).and(col(colmn).lt_eq(lit(to))))
            .select([col(colmn).mean()])
            .collect()?;

        Ok(df.column(colmn)?.f64()?.get(0).unwrap())
    }

    /// Compute the average of all entries in the specified column.
    pub fn mean_on_column(&self, col_name: &str) -> PolarsResult<f64> {
        let df = self
            .frame
            .clone()
            .select([col(col_name).mean()])
            .collect()?;
        Ok(df.column(col_name)?.f64()?.get(0).unwrap())
    }
    //======================================================================
    // DIMENSIONLESS
    /// dimensionless mass
    /// Получение безразмерной массы
    pub fn derive_dimensionless_mass(mut self, from: f64, to: f64) -> Result<Self, TGADomainError> {
        let mass_col = self.schema.mass.as_ref().unwrap().clone();
        println!("\n column: {} ", mass_col);
        let time_col = self.schema.time.as_ref().unwrap().clone();
        println!("\n column: {} ", time_col);

        let m0 = self.mean_on_interval(&mass_col, &time_col, from, to)?;

        let new_name = "alpha".to_string();

        let expr = (col(&mass_col) / lit(m0)).alias(&new_name);
        self.frame = self.frame.with_column(expr.clone());

        self.schema.columns.insert(
            new_name.clone(),
            ColumnMeta {
                name: new_name.clone(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::DimensionlessMass,
            },
        );

        self.schema.alpha = Some("alpha".to_string());

        self.log_operation(
            "derive_dimensionless_mass",
            AffectedColumns::Specific(vec!["alpha".to_string()]),
            Some(expr),
            format!(
                "Computed dimensionless mass alpha from {} (m0={:.4})",
                mass_col, m0
            ),
            true,
        );

        Ok(self)
    }

    /// Безразмерная масса
    pub fn dimensionless_mass(
        mut self,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        let mass = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
        println!("\n column: {} ", mass);

        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();
        println!("\n column: {} ", time);

        let m0 = self
            .frame
            .clone()
            .filter(col(&time).gt_eq(lit(from)).and(col(&time).lt_eq(lit(to))))
            .select([col(&mass).mean()])
            .collect()?
            .column(&mass)?
            .f64()?
            .get(0)
            .unwrap();

        if m0 <= 0.0 {
            return Err(TGADomainError::InvalidReferenceMass);
        }

        let expr = (col(&mass) / lit(m0)).alias(new_col);
        self.frame = self.frame.with_column(expr.clone());

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta {
                name: new_col.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::DimensionlessMass,
            },
        );
        self.schema.alpha = Some(new_col.to_string());

        self.log_operation(
            "dimensionless_mass",
            AffectedColumns::Specific(vec![new_col.to_string()]),
            Some(expr),
            format!(
                "Computed dimensionless mass {} from {} (m0={:.4})",
                new_col, mass, m0
            ),
            true,
        );
        info!(
            "dimensionless_mass completed in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(self)
    }
    /// conversion = 1 - dimensionless_mass
    /// Получение конверсии
    pub fn derive_conversion(mut self) -> Self {
        let src = "alpha";
        println!("\n column: {} ", src);
        let dst = "eta";

        let expr = (lit(1.0) - col(src)).alias(dst);
        self.frame = self.frame.with_column(expr.clone());

        self.schema.columns.insert(
            dst.into(),
            ColumnMeta {
                name: dst.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Conversion,
            },
        );

        self.log_operation(
            "derive_conversion",
            AffectedColumns::Specific(vec!["eta".to_string()]),
            Some(expr),
            "Computed conversion eta = 1 - alpha".to_string(),
            true,
        );

        self
    }

    /// Конверсия
    pub fn conversion(mut self, from: f64, to: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let start = Instant::now();
        let mass = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
        println!("\n column: {} ", mass);

        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();
        println!("\n column: {} ", time);

        let m0 = self
            .frame
            .clone()
            .filter(col(&time).gt_eq(lit(from)).and(col(&time).lt_eq(lit(to))))
            .select([col(&mass).mean()])
            .collect()?
            .column(&mass)?
            .f64()?
            .get(0)
            .ok_or(TGADomainError::InvalidReferenceMass)?;

        if m0 <= 0.0 {
            return Err(TGADomainError::InvalidReferenceMass);
        }

        let expr = (lit(1.0) - col(&mass) / lit(m0)).alias(new_col);
        self.frame = self.frame.with_column(expr.clone());

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta {
                name: new_col.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
                nature: ColumnNature::Conversion,
            },
        );
        self.schema.eta = Some(new_col.to_string());

        self.log_operation(
            "conversion",
            AffectedColumns::Specific(vec![new_col.to_string()]),
            Some(expr),
            format!(
                "Computed conversion {} from {} (m0={:.4})",
                new_col, mass, m0
            ),
            true,
        );
        info!(
            "conversion completed in {:?} ms",
            start.elapsed().as_millis()
        );
        Ok(self)
    }
    //================================================================================================
    //      DIAGNOSTICS AND TESTING
    /// Check for null values in all columns and print information
    /// Проверка на наличие null значений
    pub fn check_nulls(self) -> Self {
        let df = self.frame.clone().collect().unwrap();
        for col_name in df.get_column_names() {
            let null_count = df.column(col_name).unwrap().null_count();
            println!("Column '{}': {} null values", col_name, null_count);
        }
        self
    }

    /// Проверка на наличие null значений для операции
    pub fn check_nulls_for_operation(self, operation: &str) -> Self {
        let df = self.frame.clone().collect().unwrap();
        for col_name in df.get_column_names() {
            let null_count = df.column(col_name).unwrap().null_count();
            if null_count > 0 {
                println!(
                    "Opaeration {}. Column '{}': {} null values",
                    operation, col_name, null_count
                )
            };
        }
        self
    }

    /// Проверка на наличие null значений для операции (заимствованная версия)
    pub fn check_nulls_for_operation_borrowed(&self, operation: &str) {
        let df = self.frame.clone().collect().unwrap();
        for col_name in df.get_column_names() {
            let null_count = df.column(col_name).unwrap().null_count();
            if null_count > 0 {
                println!(
                    "Opaeration {}. Column '{}': {} null values",
                    operation, col_name, null_count
                )
            };
        }
    }
    /// Проверка отсутствия null значений
    /// Checks that time values are strictly increasing: t[i] > t[i-1].
    /// Returns each offending time value t[i] where monotonicity fails.
    pub fn monotony_of_time_check(&self) -> Result<Vec<f64>, TGADomainError> {
        let time_col = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?;
        let df = self.frame.clone().collect()?;
        let series = df.column(time_col)?;
        let time = series.f64()?;

        let mut failures = Vec::new();
        let mut prev: Option<f64> = None;

        for (idx, maybe_t) in time.into_iter().enumerate() {
            let t = maybe_t.ok_or_else(|| {
                TGADomainError::InvalidOperation(format!(
                    "Time monotony check failed: null at index {} in column '{}'",
                    idx, time_col
                ))
            })?;

            if let Some(prev_t) = prev {
                if t <= prev_t {
                    failures.push(t);
                }
            }
            prev = Some(t);
        }

        Ok(failures)
    }

    pub fn assert_no_nulls(self, operation: &str) -> Result<Self, TGADomainError> {
        let df = self.frame.clone().collect()?;

        let mut offenders = Vec::new();

        for col in df.get_column_names() {
            let n = df.column(col)?.null_count();
            if n > 0 {
                offenders.push((col.to_string(), n));
            }
        }

        if !offenders.is_empty() {
            return Err(TGADomainError::InvalidOperation(format!(
                "Operation '{}': nulls detected {:?}",
                operation, offenders
            )));
        }

        Ok(self)
    }

    /// Вывод статистики по колонкам
    pub fn print_column_stats(&self, cols: &[&str]) {
        let df = self.frame.clone().collect().unwrap();

        for &col in cols {
            let s = df.column(col).unwrap().f64().unwrap();
            let values: Vec<f64> = s.into_no_null_iter().collect();

            let mean = values.iter().copied().sum::<f64>() / values.len() as f64;
            let std = (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                / values.len() as f64)
                .sqrt();

            let diffs: Vec<f64> = values.windows(2).map(|w| w[1] - w[0]).collect();
            let noise = diffs.iter().map(|v| v.powi(2)).sum::<f64>().sqrt() / diffs.len() as f64;

            println!(
                "Column '{}': n={}, mean={:.4}, std={:.4}, diff_noise={:.4}",
                col,
                values.len(),
                mean,
                std,
                noise
            );
        }
    }

    /// Проверка наличия обязательных колонок
    pub fn validate_required_columns(&self, cols: &[&str]) -> Result<(), TGADomainError> {
        let schema = &self.schema;

        for &col in cols {
            if !schema.columns.contains_key(col) {
                return Err(TGADomainError::InvalidOperation(format!(
                    "Required column '{}' is missing",
                    col
                )));
            }
        }
        Ok(())
    }

    /// Samples all dataframe columns into `n` evenly spaced layers.
    ///
    /// Outer vector: sampled layers (rows),
    /// inner vector: values for each column at that sampled index.
    /// Null values are represented as `NaN`.
    /// Выборка всех колонок в n равномерно распределенных слоев
    pub fn sample_all_columns_even_layers(
        &self,
        n: usize,
    ) -> Result<Vec<Vec<f64>>, TGADomainError> {
        let df = self.frame.clone().collect()?;
        let width = df.width();
        let height = df.height();

        if n == 0 || width == 0 || height == 0 {
            return Ok(Vec::new());
        }

        let sample_count = n.min(height);
        let indices: Vec<usize> = if sample_count == 1 {
            vec![0]
        } else {
            (0..sample_count)
                .map(|i| i * (height - 1) / (sample_count - 1))
                .collect()
        };

        let mut columns_sampled: Vec<Vec<f64>> = Vec::with_capacity(width);

        for col in df.columns() {
            let col_f64 = col.cast(&DataType::Float64)?;
            let chunked = col_f64.f64()?;
            let sampled = indices
                .iter()
                .map(|&idx| chunked.get(idx).unwrap_or(f64::NAN))
                .collect::<Vec<f64>>();
            columns_sampled.push(sampled);
        }

        let mut layers: Vec<Vec<f64>> = vec![Vec::with_capacity(width); sample_count];
        for col_values in &columns_sampled {
            for (row_i, value) in col_values.iter().enumerate() {
                layers[row_i].push(*value);
            }
        }

        Ok(layers)
    }

    /// Builds a table from `sample_all_columns_even_layers` using `tabled`.
    ///
    /// `headers` should correspond to dataframe columns by position.
    /// If a header entry is `None` or missing, the dataframe column name is used.
    /// Построение таблицы из выборки колонок
    pub fn sample_all_columns_even_layers_table(
        &self,
        n: usize,
        headers: Vec<Option<String>>,
    ) -> Result<(), TGADomainError> {
        let df = self.frame.clone().collect()?;
        let layers = self.sample_all_columns_even_layers(n)?;

        if df.width() == 0 {
            return Ok(());
        }

        let mut builder = Builder::default();
        let header_row: Vec<String> = if headers.is_empty() {
            let header_row: Vec<String> = self.list_of_columns();
            header_row
        } else {
            let header_row: Vec<String> = df
                .columns()
                .iter()
                .enumerate()
                .map(|(i, col)| {
                    headers
                        .get(i)
                        .and_then(|h| h.clone())
                        .unwrap_or_else(|| col.name().to_string())
                })
                .collect();
            header_row
        };
        builder.push_record(header_row);

        for row in layers {
            builder.push_record(row.into_iter().map(|v| format!("{:.6}", v)));
        }

        let mut table = builder.build();
        table.with(Style::rounded());
        println!("{}", table);
        Ok(())
        //Ok(table.to_string())
    }
}

#[derive(Debug, Clone)]
pub struct ColumnStats {
    pub mean: f64,
    pub std: f64,
    pub mad: f64,
    pub min: f64,
    pub max: f64,
}
use crate::Kinetics::experimental_kinetics::exp_kinetics_smooth_filter::extract_f64_column;
/// Статистика по колонке
pub fn column_stats(df: &DataFrame, col: &str) -> Result<ColumnStats, TGADomainError> {
    let s = df
        .column(col)?
        .f64()
        .map_err(|_| TGADomainError::InvalidColumnType(col.into()))?;

    let values: Vec<f64> = s.into_no_null_iter().collect();
    if values.is_empty() {
        return Err(TGADomainError::EmptyColumn(col.into()));
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let std = (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64).sqrt();

    let mut sorted = values.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = sorted[sorted.len() / 2];
    let mad = sorted.iter().map(|v| (v - med).abs()).sum::<f64>() / sorted.len() as f64;

    Ok(ColumnStats {
        mean,
        std,
        mad,
        min: *sorted.first().unwrap(),
        max: *sorted.last().unwrap(),
    })
}

#[derive(Debug)]
pub struct StatsDiff {
    pub before: ColumnStats,
    pub after: ColumnStats,
}

impl StatsDiff {
    pub fn std_ratio(&self) -> f64 {
        self.after.std / self.before.std
    }

    pub fn mad_ratio(&self) -> f64 {
        self.after.mad / self.before.mad
    }

    pub fn improvement(&self) -> bool {
        self.std_ratio() < 1.0 && self.mad_ratio() < 1.0
    }
}

/// Сравнение статистики колонок
pub fn diff_column_stats(
    before: &DataFrame,
    after: &DataFrame,
    col: &str,
) -> Result<StatsDiff, TGADomainError> {
    Ok(StatsDiff {
        before: column_stats(before, col)?,
        after: column_stats(after, col)?,
    })
}

#[derive(Debug, Clone, Copy)]
pub struct OperationGuarantees {
    pub preserves_length: bool,
    pub produces_nulls: bool,
    pub preserves_monotonicity: bool,
}

use std::collections::HashMap;
/// Гарантии операций
pub fn operation_guarantees() -> HashMap<&'static str, OperationGuarantees> {
    use OperationGuarantees as G;

    HashMap::from([
        (
            "trim_edges",
            G {
                preserves_length: false,
                produces_nulls: false,
                preserves_monotonicity: true,
            },
        ),
        (
            "hampel_filter",
            G {
                preserves_length: true,
                produces_nulls: false,
                preserves_monotonicity: true,
            },
        ),
        (
            "rolling_mean",
            G {
                preserves_length: true,
                produces_nulls: true,
                preserves_monotonicity: true,
            },
        ),
        (
            "sg_filter",
            G {
                preserves_length: true,
                produces_nulls: false,
                preserves_monotonicity: true,
            },
        ),
        (
            "spline_smooth",
            G {
                preserves_length: true,
                produces_nulls: false,
                preserves_monotonicity: true,
            },
        ),
        (
            "derive_rate",
            G {
                preserves_length: false,
                produces_nulls: true,
                preserves_monotonicity: false,
            },
        ),
    ])
}
