use crate::gui::experimental_kinetics_gui::model::TGAGUIError;
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};

const DEFAULT_N_POINTS: usize = 1000;
const DEFAULT_SYMBOLIC_EXPRESSION: &str = "mass";

#[derive(Debug, Clone, Copy, Deserialize)]
pub struct CalibrationLine {
    k: f64,
    b: f64,
}

impl CalibrationLine {
    pub fn new(k: f64, b: f64) -> Self {
        Self { k, b }
    }

    pub fn k(&self) -> f64 {
        self.k
    }

    pub fn b(&self) -> f64 {
        self.b
    }

    pub fn set_k(&mut self, k: f64) {
        self.k = k;
    }

    pub fn set_b(&mut self, b: f64) {
        self.b = b;
    }
}

impl Default for CalibrationLine {
    fn default() -> Self {
        Self { k: 1.0, b: 0.0 }
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct Settings {
    calibration_line: Option<CalibrationLine>,
    n_points: Option<usize>,
    symbolic_expression: Option<String>,
}

#[derive(Debug, Deserialize)]
struct SettingsFile {
    calibration_line: Option<CalibrationLine>,
    n_points: Option<usize>,
    symbolic_expression: Option<String>,
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            calibration_line: Some(CalibrationLine::default()),
            n_points: Some(DEFAULT_N_POINTS),
            symbolic_expression: Some(DEFAULT_SYMBOLIC_EXPRESSION.to_string()),
        }
    }
}

impl Settings {
    pub fn new() -> Result<Self, TGAGUIError> {
        Self::load_from_default_locations()
    }

    pub fn calibration_line(&self) -> Option<CalibrationLine> {
        self.calibration_line
    }

    pub fn n_points(&self) -> Option<usize> {
        self.n_points
    }

    pub fn symbolic_expression(&self) -> Option<&str> {
        self.symbolic_expression.as_deref()
    }

    pub fn set_calibration_line(&mut self, calibration_line: Option<CalibrationLine>) {
        self.calibration_line = calibration_line;
    }

    pub fn set_calibration_coeffs(&mut self, k: f64, b: f64) {
        self.calibration_line = Some(CalibrationLine::new(k, b));
    }

    pub fn set_n_points(&mut self, n_points: Option<usize>) {
        self.n_points = n_points;
    }

    pub fn set_symbolic_expression<S: Into<String>>(&mut self, symbolic_expression: Option<S>) {
        self.symbolic_expression = symbolic_expression.map(Into::into);
    }

    pub fn load_from_path(path: &Path) -> Result<Self, TGAGUIError> {
        let raw = fs::read_to_string(path).map_err(|e| {
            TGAGUIError::SettingsErrors(format!(
                "failed to read settings file '{}': {e}",
                path.display()
            ))
        })?;

        let parsed: SettingsFile = serde_json::from_str(&raw).map_err(|e| {
            TGAGUIError::SettingsErrors(format!(
                "failed to parse settings json '{}': {e}",
                path.display()
            ))
        })?;

        let mut settings = Self::default();
        if let Some(calibration_line) = parsed.calibration_line {
            settings.calibration_line = Some(calibration_line);
        }
        if let Some(n_points) = parsed.n_points {
            settings.n_points = Some(n_points);
        }
        if let Some(symbolic_expression) = parsed.symbolic_expression {
            settings.symbolic_expression = Some(symbolic_expression);
        }
        Ok(settings)
    }

    pub fn load_from_default_locations() -> Result<Self, TGAGUIError> {
        for candidate in Self::default_config_candidates() {
            if candidate.is_file() {
                return Self::load_from_path(&candidate);
            }
        }
        Ok(Self::default())
    }

    fn default_config_candidates() -> Vec<PathBuf> {
        let names = [
            "tgasettings",
            "tgasettings.json",
            "tga_settings",
            "tga_settings.json",
            "tga-settings",
            "tga-settings.json",
            "TGASettings",
            "TGASettings.json",
        ];

        names.into_iter().map(PathBuf::from).collect()
    }
}
