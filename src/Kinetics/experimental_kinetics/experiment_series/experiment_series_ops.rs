//! Операции преобразования для `TGAExperiment`.
//!
//! Здесь собраны тонкие обёртки над `TGADataset`, которые меняют
//! содержимое одного эксперимента, но сохраняют его `meta`.
//! Это место для scale/offset/unary/mass/conversion helpers.

use super::*;

impl TGAExperiment {
    //======================================================================
    // ARBITRARY ALGEBRAIC TRANSFORMATIONS ON COLUMN DATA
    //======================================================================
    pub fn scale_columns(self, cols: &[&str], factor: f64) -> Self {
        let dataset = self.dataset.scale_columns(cols, factor);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn scale_column(self, col: &str, factor: f64) -> Self {
        let dataset = self.dataset.scale_column(col, factor);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn scale_column_in_its_range(self, colmn: &str, s: f64, from: f64, to: f64) -> Self {
        let dataset = self.dataset.scale_column_in_its_range(colmn, s, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn scale_column_in_range_by_reference(
        self,
        target_col: &str,
        reference_col: &str,
        s: f64,
        from: f64,
        to: f64,
    ) -> Self {
        let dataset =
            self.dataset
                .scale_column_in_range_by_reference(target_col, reference_col, s, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn offset_column(self, colmn: &str, offset: f64) -> Self {
        let dataset = self.dataset.offset_column(colmn, offset);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn offset_column_in_its_range(self, colmn: &str, offset: f64, from: f64, to: f64) -> Self {
        let dataset = self
            .dataset
            .offset_column_in_its_range(colmn, offset, from, to);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn offset_column_in_range_by_reference(
        self,
        target_col: &str,
        reference_col: &str,
        offset: f64,
        from: f64,
        to: f64,
    ) -> Self {
        let dataset = self.dataset.offset_column_in_range_by_reference(
            target_col,
            reference_col,
            offset,
            from,
            to,
        );
        Self {
            dataset,
            meta: self.meta,
        }
    }

    pub fn calibrate_mass_from_voltage(self, k: f64, b: f64) -> Self {
        let dataset = self.dataset.calibrate_mass_from_voltage(k, b);
        Self {
            dataset,
            meta: self.meta,
        }
    }

    //======================================================================
    // UNARY / DERIVED COLUMNS
    //======================================================================
    pub fn unary_column_op(
        self,
        src: &str,
        dst: &str,
        op: UnaryOp,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.unary_column_op(src, dst, op)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn exp_column(self, col_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.exp_column(col_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn ln_column(self, col_name: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.ln_column(col_name)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    //======================================================================
    // DIMENSIONLESS TRANSFORMS
    //======================================================================
    pub fn derive_dimensionless_mass(self, from: f64, to: f64) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_dimensionless_mass(from, to)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn derive_mass_rate_from_column(
        self,
        source_col: &str,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .derive_mass_rate_from_column(source_col, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn dimensionless_mass(
        self,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.dimensionless_mass(from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn dimensionless_mass_from_column(
        self,
        source_col: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .dimensionless_mass_from_column(source_col, from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn conversion(self, from: f64, to: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.conversion(from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn conversion_from_column(
        self,
        source_col: &str,
        from: f64,
        to: f64,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .conversion_from_column(source_col, from, to, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn conversion_with_m0(self, m0: f64, new_col: &str) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.conversion_with_m0(m0, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn derive_temperature_rate_from_column(
        self,
        source_col: &str,
        new_col: &str,
    ) -> Result<Self, TGADomainError> {
        let dataset = self
            .dataset
            .derive_temperature_rate_from_column(source_col, new_col)?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }

    pub fn derive_conversion(self) -> Result<Self, TGADomainError> {
        let dataset = self.dataset.derive_conversion()?;
        Ok(Self {
            dataset,
            meta: self.meta,
        })
    }
}
