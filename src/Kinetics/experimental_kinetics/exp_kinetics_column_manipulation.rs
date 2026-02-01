use crate::Kinetics::experimental_kinetics::exp_kinetics_main::{
    ColumnMeta, ColumnOrigin, TGADataset, TGADomainError, UnaryOp, Unit,
};
use polars::error::PolarsResult;
use polars::prelude::DataType;
use polars::prelude::Expr;
use polars::prelude::*;
impl TGADataset {
    //================================================================
    // BASIC TRANSFORMATIONS: COLUMN CUT AND FILTER
    pub fn with_column_expr(mut self, meta: ColumnMeta, expr: Expr) -> Self {
        self.frame = self.frame.with_column(expr.alias(&meta.name));
        self.schema.columns.insert(meta.name.clone(), meta);
        self
    }
    /// Works on all columns at once
    ///✔ Zero copy
    /// ✔ Lazy
    pub fn filter_rows(self, predicate: Expr) -> Self {
        let frame = self.frame.filter(predicate);
        Self {
            frame,
            schema: self.schema,
        }
    }
    pub fn filter_by_mask(self, mask: Expr) -> Self {
        self.filter_rows(mask)
    }
    /// Handy helper: cut interval
    pub fn cut_interval(self, colmn: &str, from: f64, to: f64) -> Self {
        self.filter_rows(
            // NOT (from < x < to)
            // Логика:
            // col > from AND col < to — строки внутри интервал
            //.not() — оставить всё, кроме этого интервала
            (col(colmn).gt(lit(from)).and(col(colmn).lt(lit(to)))).not(),
        )
    }
    pub fn cut_time_interval(self, from: f64, to: f64) -> Self {
        let time_col = self.schema.time.as_ref().unwrap().clone();
        self.cut_interval(&time_col, from, to)
    }
    pub fn cut_temperature_interval(self, from: f64, to: f64) -> Self {
        let temp_col = self.schema.temperature.as_ref().unwrap().clone();
        self.cut_interval(&temp_col, from, to)
    }
    pub fn cut_mass_interval(self, from: f64, to: f64) -> Self {
        let mass_col = self.schema.mass.as_ref().unwrap().clone();
        self.cut_interval(&mass_col, from, to)
    }
    /// Обрезка данных до времени t_end
    pub fn cut_before_time(self, t_start: f64) -> Self {
        let time_col = self.schema.time.as_ref().unwrap().clone();

        self.filter_rows(col(&time_col).gt_eq(lit(t_start)))
    }
    /// trim edges
    pub fn trim_range(self, column: &str, from: f64, to: f64) -> Self {
        self.filter_rows(
            // from ≤ x ≤ to
            col(column).gt_eq(lit(from)).and(col(column).lt_eq(lit(to))),
        )
    }
    // When differentiation is used the last and the first
    // points of the rate column is nill
    pub fn trim_edges(self, left: usize, right: usize) -> Self {
        let df = self.frame.collect().unwrap();
        let total = df.height();
        let length = total.saturating_sub(left + right);
        let sliced_df = df.slice(left as i64, length);
        let frame = sliced_df.lazy();

        Self {
            frame,
            schema: self.schema,
        }
    }
    pub fn rename_column(mut self, old: &str, new: &str) -> Result<Self, TGADomainError> {
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

        Ok(self)
    }

    //======================================================================
    // UNIT TRANSFORMATIONS
    /// Это не просто offset, это семантическая операция, меняющая единицы.
    pub fn celsius_to_kelvin(mut self) -> Self {
        let col_name = self.schema.temperature.as_ref().unwrap().clone();

        self.frame = self
            .frame
            .with_column((col(&col_name) + lit(273.15)).alias(&col_name));
        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Kelvin;
            meta.origin = ColumnOrigin::PolarsDerived;
        }
        self
    }
    /// Делим время на 3600 и меняем единицы на Hour
    pub fn seconds_to_hours(mut self) -> Self {
        let col_name = self.schema.time.as_ref().unwrap().clone();
        self.frame = self
            .frame
            .with_column((col(&col_name) / lit(3600.0)).alias(&col_name));

        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Hour;
            meta.origin = ColumnOrigin::PolarsDerived;
        }

        self
    }
    //=======================================================================
    // ARBITRARY ALGEBRAIC TRANSFORMATIONS ON COLUMN DATA
    /// Scale specified columns by factor
    /// Unit handling is not implemented yet
    pub fn scale_columns(mut self, cols: &[&str], factor: f64) -> Self {
        let exprs: Vec<Expr> = cols
            .iter()
            .map(|&c| (col(c) * lit(factor)).alias(c))
            .collect();

        self.frame = self.frame.with_columns(exprs);

        for &c in cols {
            if let Some(meta) = self.schema.columns.get_mut(c) {
                meta.origin = ColumnOrigin::PolarsDerived;
                // unit пока оставляем, но позже тут будет умножение единиц
            }
        }
        self
    }
    /// прибавление константы ко всем значениям в колонке
    pub fn offset_column(mut self, colmn: &str, offset: f64) -> Self {
        self.frame = self
            .frame
            .with_column((col(colmn) + lit(offset)).alias(colmn));
        self
    }

    /// калибровочная прямая, от милливольт прибора к мг массы
    pub fn calibrate_mass_from_voltage(mut self, k: f64, b: f64) -> Self {
        let col_name = self.schema.mass.as_ref().unwrap().clone();

        self.frame = self
            .frame
            .with_column((col(&col_name) * lit(k) + lit(b)).alias(&col_name));

        if let Some(meta) = self.schema.columns.get_mut(&col_name) {
            meta.unit = Unit::Milligram;
            meta.origin = ColumnOrigin::PolarsDerived;
        }

        self
    }

    pub fn calibrate_mass(mut self, k: f64, b: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let src = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();
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
        // обновляем роль mass
        self.schema.mass = Some(new_col.into());
        // The ColumnMeta is already inserted by unary_column_op, assuming it does that.

        Ok(self)
    }

    pub fn unary_column_op(
        mut self,
        src: &str,
        dst: &str,
        op: UnaryOp,
    ) -> Result<Self, TGADomainError> {
        // 1. Domain checks (if any)
        if let Some(check) = op.domain_check {
            check(&self, src)?;
        }

        // 2. Output type callback for polars map
        let dt = |_schema: &Schema, field: &Field| {
            Ok(Field::new(field.name().clone(), DataType::Float64))
        };

        // Extract fields before moving op
        let out_unit = op.output_unit;
        let func = op.func; // Box<dyn Fn(f64) -> f64 + Send + Sync>

        // 3. Polars Expr using the boxed function
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

        // 4. Schema update
        self.schema.columns.insert(
            dst.to_string(),
            ColumnMeta {
                name: dst.to_string(),
                unit: out_unit,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        Ok(self)
    }

    pub fn exp_column(mut self, col_name: &str) -> Result<Self, TGADomainError> {
        // Проверка: есть ли x <= 0
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
                name: new_name,
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );
        Ok(self)
    }

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
        // Проверка: есть ли x <= 0
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
                name: new_name,
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );
        Ok(self)
    }
    /// среднее по интервалу
    pub fn mean_on_interval(
        &self,
        value_col: &str,
        time_col: &str,
        from: f64,
        to: f64,
    ) -> PolarsResult<f64> {
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
    //======================================================================
    // DIMENSIONLESS
    /// dimensionless mass
    pub fn derive_dimensionless_mass(mut self, from: f64, to: f64) -> PolarsResult<Self> {
        let mass_col = self.schema.mass.as_ref().unwrap().clone();
        let time_col = self.schema.time.as_ref().unwrap().clone();

        let m0 = self.mean_on_interval(&mass_col, &time_col, from, to)?;

        let new_name = "dimensionless_mass".to_string();

        self.frame = self
            .frame
            .with_column((col(&mass_col) / lit(m0)).alias(&new_name));

        self.schema.columns.insert(
            new_name.clone(),
            ColumnMeta {
                name: new_name,
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        Ok(self)
    }

    pub fn dimensionless_mass(
        mut self,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let mass = self
            .schema
            .mass
            .as_ref()
            .ok_or(TGADomainError::MassNotBound)?
            .clone();

        let time = self
            .schema
            .time
            .as_ref()
            .ok_or(TGADomainError::TimeNotBound)?
            .clone();

        // вычисляем m0
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

        self.frame = self
            .frame
            .with_column((col(&mass) / lit(m0)).alias(new_col));

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta {
                name: new_col.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        Ok(self)
    }
    /// conversion = 1 - dimensionless_mass
    pub fn derive_conversion(mut self) -> Self {
        let src = "dimensionless_mass";
        let dst = "conversion";

        self.frame = self.frame.with_column((lit(1.0) - col(src)).alias(dst));

        self.schema.columns.insert(
            dst.into(),
            ColumnMeta {
                name: dst.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        self
    }

    pub fn conversion(mut self, src: &str, new_col: &str) -> Result<Self, TGADomainError> {
        let meta = self
            .schema
            .columns
            .get(src)
            .ok_or(TGADomainError::ColumnNotFound(src.into()))?;

        if meta.unit != Unit::Dimensionless {
            return Err(TGADomainError::InvalidUnitForConversion);
        }

        self.frame = self.frame.with_column((lit(1.0) - col(src)).alias(new_col));

        self.schema.columns.insert(
            new_col.into(),
            ColumnMeta {
                name: new_col.into(),
                unit: Unit::Dimensionless,
                origin: ColumnOrigin::PolarsDerived,
            },
        );

        Ok(self)
    }
}
