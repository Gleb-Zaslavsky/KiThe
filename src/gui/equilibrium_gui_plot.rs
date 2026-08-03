//! Shared plot-domain adapter for chemical-equilibrium results.
//!
//! Both the embedded egui/Plotters window and KiThePlot consume this value.
//! They therefore cannot silently disagree about labels, temperature order,
//! or the selected physical result basis.

use crate::gui::equilibrium_gui_model::GuiResultBasis;
use crate::gui::equilibrium_gui_result::EquilibriumGuiResultSnapshot;
use crate::gui::gui_plot::PlotWindow;
use RustedSciThe::numerical::optimization::inter_n_extrapolate::{InterpolationSpace, Pchip};
use eframe::egui;
use kithe_plot::controller::PlotController;
use kithe_plot::model::DataSource;
use kithe_plot::view::PlotEditorView;
use nalgebra::{DMatrix, DVector};

/// Unit carried by one equilibrium plot family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumGuiPlotUnit {
    Moles,
    MoleFraction,
    Unknown,
}

/// Columnar equilibrium data shared by both plot implementations.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumGuiPlotData {
    x_label: String,
    x_values: Vec<f64>,
    columns: Vec<(String, Vec<f64>)>,
    column_masks: Vec<Vec<bool>>,
    y_unit: EquilibriumGuiPlotUnit,
}

fn split_plot_columns(
    items: Vec<((String, Vec<f64>), Vec<bool>)>,
    y_unit: EquilibriumGuiPlotUnit,
) -> (
    Vec<(String, Vec<f64>)>,
    Vec<Vec<bool>>,
    EquilibriumGuiPlotUnit,
) {
    let (columns, masks) = items.into_iter().unzip();
    (columns, masks, y_unit)
}

impl EquilibriumGuiPlotData {
    /// Builds one aligned column source from an immutable result snapshot.
    pub fn from_snapshot(
        snapshot: &EquilibriumGuiResultSnapshot,
        basis: GuiResultBasis,
    ) -> Result<Self, String> {
        let x_values = snapshot.temperatures();
        if x_values.is_empty() {
            return Err("equilibrium result has no accepted points".into());
        }
        let (columns, column_masks, y_unit) = match basis {
            GuiResultBasis::ComponentMoles => split_plot_columns(
                snapshot
                    .component_labels()
                    .iter()
                    .enumerate()
                    .map(|(index, label)| {
                        let values = snapshot
                            .component_mole_series(index)
                            .ok_or_else(|| format!("missing mole series for '{label}'"))?;
                        let masks = snapshot
                            .points()
                            .iter()
                            .map(|point| point.component_active()[index])
                            .collect();
                        Ok(((label.clone(), values), masks))
                    })
                    .collect::<Result<Vec<((String, Vec<f64>), Vec<bool>)>, String>>()?,
                EquilibriumGuiPlotUnit::Moles,
            ),
            GuiResultBasis::MoleFractions => split_plot_columns(
                snapshot
                    .component_labels()
                    .iter()
                    .enumerate()
                    .map(|(index, label)| {
                        let values = snapshot
                            .component_fraction_series(index)
                            .ok_or_else(|| format!("missing mole-fraction series for '{label}'"))?;
                        let masks = snapshot
                            .points()
                            .iter()
                            .map(|point| point.component_active()[index])
                            .collect();
                        Ok(((label.clone(), values), masks))
                    })
                    .collect::<Result<Vec<((String, Vec<f64>), Vec<bool>)>, String>>()?,
                EquilibriumGuiPlotUnit::MoleFraction,
            ),
            GuiResultBasis::PhaseTotals => split_plot_columns(
                snapshot
                    .phase_labels()
                    .iter()
                    .enumerate()
                    .map(|(index, label)| {
                        let values = snapshot
                            .phase_total_series(index)
                            .ok_or_else(|| format!("missing phase-total series for '{label}'"))?;
                        let masks = snapshot
                            .points()
                            .iter()
                            .map(|point| point.phase_active()[index])
                            .collect();
                        Ok(((label.clone(), values), masks))
                    })
                    .collect::<Result<Vec<((String, Vec<f64>), Vec<bool>)>, String>>()?,
                EquilibriumGuiPlotUnit::Moles,
            ),
        };
        Self::from_columns_with_metadata("temperature_k", x_values, columns, column_masks, y_unit)
    }

    /// Creates an aligned column source for adapter tests and future derived
    /// result views. Every y-column must have exactly the x length.
    pub fn from_columns(
        x_label: impl Into<String>,
        x_values: Vec<f64>,
        columns: Vec<(String, Vec<f64>)>,
    ) -> Result<Self, String> {
        let masks = columns
            .iter()
            .map(|(_, values)| vec![true; values.len()])
            .collect::<Vec<_>>();
        Self::from_columns_with_metadata(
            x_label,
            x_values,
            columns,
            masks,
            EquilibriumGuiPlotUnit::Unknown,
        )
    }

    fn from_columns_with_metadata(
        x_label: impl Into<String>,
        x_values: Vec<f64>,
        columns: Vec<(String, Vec<f64>)>,
        column_masks: Vec<Vec<bool>>,
        y_unit: EquilibriumGuiPlotUnit,
    ) -> Result<Self, String> {
        if x_values.is_empty() {
            return Err("plot x-axis must contain at least one point".into());
        }
        if x_values.iter().any(|value| !value.is_finite()) {
            return Err("plot x-axis must contain finite values".into());
        }
        if columns.is_empty() {
            return Err("plot source must contain at least one y-column".into());
        }
        if column_masks.len() != columns.len() {
            return Err("plot column masks must match the column count".into());
        }
        for (label, values) in &columns {
            if label.trim().is_empty() {
                return Err("plot column labels must not be empty".into());
            }
            if values.len() != x_values.len() {
                return Err(format!(
                    "plot column '{label}' has {} values but x-axis has {} points",
                    values.len(),
                    x_values.len()
                ));
            }
            if values.iter().any(|value| !value.is_finite()) {
                return Err(format!("plot column '{label}' contains a non-finite value"));
            }
        }
        for (index, mask) in column_masks.iter().enumerate() {
            if mask.len() != x_values.len() {
                return Err(format!(
                    "plot column mask {index} has {} values but x-axis has {} points",
                    mask.len(),
                    x_values.len()
                ));
            }
        }
        Ok(Self {
            x_label: x_label.into(),
            x_values,
            columns,
            column_masks,
            y_unit,
        })
    }

    pub fn x_label(&self) -> &str {
        &self.x_label
    }

    pub fn x_values(&self) -> &[f64] {
        &self.x_values
    }

    pub fn column_labels(&self) -> impl Iterator<Item = &str> {
        self.columns.iter().map(|(label, _)| label.as_str())
    }

    pub fn column_values(&self, label: &str) -> Option<&[f64]> {
        self.columns
            .iter()
            .find(|(name, _)| name == label)
            .map(|(_, values)| values.as_slice())
    }

    pub fn column_mask(&self, label: &str) -> Option<&[bool]> {
        self.columns
            .iter()
            .position(|(name, _)| name == label)
            .map(|index| self.column_masks[index].as_slice())
    }

    pub fn y_unit(&self) -> EquilibriumGuiPlotUnit {
        self.y_unit
    }

    /// Builds a display-only log10 projection shared by both plot adapters.
    ///
    /// The accepted snapshot and its zeros are untouched. A non-positive
    /// series is rejected instead of being converted to `-inf` or silently
    /// discarded, because neither adapter has a common missing-value axis
    /// contract yet.
    pub fn with_log10_y_scale(&self) -> Result<Self, String> {
        let mut columns = Vec::with_capacity(self.columns.len());
        for (label, values) in &self.columns {
            if values.iter().any(|value| *value <= 0.0) {
                return Err(format!(
                    "log10 display requires strictly positive values in '{label}'"
                ));
            }
            columns.push((
                format!("log10({label})"),
                values.iter().map(|value| value.log10()).collect(),
            ));
        }
        Self::from_columns_with_metadata(
            self.x_label.clone(),
            self.x_values.clone(),
            columns,
            self.column_masks.clone(),
            self.y_unit,
        )
    }

    /// Returns a display-only view containing only the requested series.
    ///
    /// Visibility is deliberately applied after a result has been accepted.
    /// It changes neither the immutable solver snapshot nor the exact grid,
    /// and therefore cannot turn a display choice into a stale calculation.
    pub fn retain_columns<I, S>(&self, visible_labels: I) -> Result<Self, String>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        let visible = visible_labels
            .into_iter()
            .map(|label| label.as_ref().to_string())
            .collect::<std::collections::HashSet<_>>();
        let mut columns = Vec::new();
        let mut masks = Vec::new();
        for (index, (label, values)) in self.columns.iter().enumerate() {
            if visible.contains(label) {
                columns.push((label.clone(), values.clone()));
                masks.push(self.column_masks[index].clone());
            }
        }
        if columns.is_empty() {
            return Err("plot visibility leaves no selected series".into());
        }
        Self::from_columns_with_metadata(
            self.x_label.clone(),
            self.x_values.clone(),
            columns,
            masks,
            self.y_unit,
        )
    }

    /// Creates a display-only PCHIP source without changing the exact solver
    /// grid kept by the accepted result snapshot.
    pub fn resample_pchip(
        &self,
        output_points: usize,
        space: InterpolationSpace,
        clamp: bool,
    ) -> Result<Self, String> {
        if output_points < 2 {
            return Err("PCHIP display grid requires at least two points".into());
        }
        if self.x_values.len() < 2 {
            return Err("PCHIP display resampling needs at least two solved points".into());
        }
        let descending = if self.x_values[1] > self.x_values[0] {
            false
        } else if self.x_values[1] < self.x_values[0] {
            true
        } else {
            return Err("PCHIP temperature grid must be strictly monotonic".into());
        };
        if !self.x_values.windows(2).all(|window| {
            if descending {
                window[1] < window[0]
            } else {
                window[1] > window[0]
            }
        }) {
            return Err("PCHIP temperature grid must be strictly monotonic".into());
        }

        let x_ascending = if descending {
            self.x_values.iter().rev().copied().collect::<Vec<_>>()
        } else {
            self.x_values.clone()
        };
        let target_ascending = plot_linspace(
            x_ascending[0],
            *x_ascending.last().expect("validated x grid is non-empty"),
            output_points,
        );
        let x_values = if descending {
            target_ascending.iter().rev().copied().collect()
        } else {
            target_ascending.clone()
        };

        let mut columns = Vec::with_capacity(self.columns.len());
        for (label, values) in &self.columns {
            let y_ascending = if descending {
                values.iter().rev().copied().collect::<Vec<_>>()
            } else {
                values.clone()
            };
            if matches!(space, InterpolationSpace::Log)
                && y_ascending.iter().any(|value| *value <= 0.0)
            {
                return Err(format!(
                    "log-space PCHIP requires strictly positive values in '{label}'"
                ));
            }
            let interpolator = Pchip::new(&x_ascending, &y_ascending, space);
            let mut resampled = target_ascending
                .iter()
                .map(|x| interpolator.eval(*x, clamp))
                .collect::<Vec<_>>();
            if descending {
                resampled.reverse();
            }
            columns.push((label.clone(), resampled));
        }
        // Lifecycle masks are categorical, not interpolated physical values.
        // Carry the nearest solved state onto the display grid instead of
        // turning an active/inactive transition into a fractional boolean.
        let mut masks = Vec::with_capacity(self.column_masks.len());
        for mask in &self.column_masks {
            let mask_ascending = if descending {
                mask.iter().rev().copied().collect::<Vec<_>>()
            } else {
                mask.clone()
            };
            let mut resampled = target_ascending
                .iter()
                .map(|target| {
                    let mut nearest = 0;
                    let mut nearest_distance = f64::INFINITY;
                    for (index, value) in x_ascending.iter().enumerate() {
                        let distance = (*value - *target).abs();
                        if distance < nearest_distance {
                            nearest = index;
                            nearest_distance = distance;
                        }
                    }
                    mask_ascending[nearest]
                })
                .collect::<Vec<_>>();
            if descending {
                resampled.reverse();
            }
            masks.push(resampled);
        }
        Self::from_columns_with_metadata(
            self.x_label.clone(),
            x_values,
            columns,
            masks,
            self.y_unit,
        )
    }

    /// Converts to the existing embedded plot window without changing data.
    pub fn to_embedded_plot_window(&self) -> PlotWindow {
        let mut matrix_data = Vec::with_capacity(self.x_values.len() * self.columns.len());
        for row in 0..self.x_values.len() {
            for (_, values) in &self.columns {
                matrix_data.push(values[row]);
            }
        }
        PlotWindow::new(
            self.x_label.clone(),
            self.columns
                .iter()
                .map(|(label, _)| label.clone())
                .collect(),
            DVector::from_vec(self.x_values.clone()),
            DMatrix::from_row_slice(self.x_values.len(), self.columns.len(), &matrix_data),
        )
    }
}

impl DataSource for EquilibriumGuiPlotData {
    fn column(&self, name: &str) -> Option<Vec<f64>> {
        if name == self.x_label {
            Some(self.x_values.clone())
        } else {
            self.column_values(name).map(ToOwned::to_owned)
        }
    }

    fn column_names(&self) -> Vec<String> {
        std::iter::once(self.x_label.clone())
            .chain(self.columns.iter().map(|(label, _)| label.clone()))
            .collect()
    }

    fn len(&self) -> usize {
        self.x_values.len()
    }
}

fn plot_linspace(start: f64, end: f64, count: usize) -> Vec<f64> {
    let denominator = (count - 1) as f64;
    (0..count)
        .map(|index| start + (end - start) * index as f64 / denominator)
        .collect()
}

/// Owned KiThePlot editor state for equilibrium result data.
///
/// The controller receives a clone of the validated columnar adapter. It can
/// therefore outlive the solver worker and cannot observe mutable document or
/// solver state after a result has been accepted.
#[derive(Default)]
pub struct EquilibriumGuiKiThePlotWindow {
    pub open: bool,
    controller: Option<PlotController>,
    view: Option<PlotEditorView>,
}

impl EquilibriumGuiKiThePlotWindow {
    const VIEWPORT_ID: &'static str = "equilibrium_kithe_plot_viewport";

    pub fn open_from_data(&mut self, data: &EquilibriumGuiPlotData) -> Result<(), String> {
        if data.column_names().len() < 2 || data.len() == 0 {
            return Err("equilibrium plot needs an x-axis and at least one series".into());
        }
        let mut controller = PlotController::new();
        controller
            .load_from_data_source(data)
            .map_err(|error| format!("failed to load equilibrium plot data: {error}"))?;
        self.controller = Some(controller);
        self.view = Some(PlotEditorView::new());
        self.open = true;
        Ok(())
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        if !self.open {
            return;
        }
        let Some(controller) = self.controller.as_mut() else {
            self.open = false;
            return;
        };
        let Some(view) = self.view.as_mut() else {
            self.open = false;
            return;
        };
        let viewport_id = egui::ViewportId::from_hash_of(Self::VIEWPORT_ID);
        let builder = egui::ViewportBuilder::default()
            .with_title("Chemical equilibrium plots")
            .with_inner_size([1200.0, 860.0]);
        let open = &mut self.open;
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            if ctx.input(|input| input.viewport().close_requested()) {
                *open = false;
                return;
            }
            match class {
                egui::ViewportClass::EmbeddedWindow => {
                    egui::Window::new("Chemical equilibrium plots")
                        .open(open)
                        .default_size([1200.0, 860.0])
                        .show(ctx, |ui| {
                            let actions = view.draw(ui, controller);
                            for action in actions {
                                let _ = controller.dispatch(action);
                            }
                        });
                }
                _ => {
                    let actions = view.draw(ctx, controller);
                    for action in actions {
                        let _ = controller.dispatch(action);
                    }
                }
            }
        });
        if !self.open {
            self.controller = None;
            self.view = None;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kithe_plot::model::DataSource;

    #[test]
    fn embedded_and_kithe_plot_adapters_share_the_same_columns() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![
                ("gas::H2O".into(), vec![0.9, 0.7]),
                ("liquid::H2O".into(), vec![0.1, 0.3]),
            ],
        )
        .expect("aligned data builds");
        assert_eq!(
            data.column_names(),
            vec!["temperature_k", "gas::H2O", "liquid::H2O"]
        );
        assert_eq!(data.column("gas::H2O"), Some(vec![0.9, 0.7]));
        assert_eq!(data.len(), 2);

        let window = data.to_embedded_plot_window();
        assert_eq!(window.arg, "temperature_k");
        assert_eq!(window.values, vec!["gas::H2O", "liquid::H2O"]);
        assert_eq!(window.t_result.as_slice(), &[300.0, 600.0]);
        assert_eq!(window.y_result[(0, 0)], 0.9);
        assert_eq!(window.y_result[(1, 0)], 0.7);
        assert_eq!(window.y_result[(0, 1)], 0.1);
        assert_eq!(window.y_result[(1, 1)], 0.3);
    }

    #[test]
    fn adapter_rejects_misaligned_or_nonfinite_series() {
        let length_error = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("H2O".into(), vec![1.0])],
        )
        .expect_err("misaligned series must fail");
        assert!(length_error.contains("has 1 values"));

        let finite_error = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, f64::NAN],
            vec![("H2O".into(), vec![1.0, 0.0])],
        )
        .expect_err("nonfinite axis must fail");
        assert!(finite_error.contains("finite"));
    }

    #[test]
    fn kithe_plot_adapter_accepts_the_same_validated_columns() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("gas::H2O".into(), vec![1.0, 0.5])],
        )
        .expect("aligned data builds");
        let mut window = EquilibriumGuiKiThePlotWindow::default();
        window
            .open_from_data(&data)
            .expect("KiThePlot accepts aligned columns");
        assert!(window.open);
    }

    #[test]
    fn pchip_resampling_preserves_descending_direction_and_clamped_bounds() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![1000.0, 800.0, 600.0],
            vec![("gas::fuel".into(), vec![0.9, 0.5, 0.1])],
        )
        .expect("aligned data builds");
        let resampled = data
            .resample_pchip(9, InterpolationSpace::Linear, true)
            .expect("descending PCHIP succeeds");
        assert_eq!(resampled.len(), 9);
        assert!(resampled.x_values()[0] > resampled.x_values()[8]);
        let values = resampled.column_values("gas::fuel").expect("series exists");
        assert!(values.iter().all(|value| (0.1..=0.9).contains(value)));
    }

    #[test]
    fn log_pchip_rejects_zero_display_values_instead_of_panicking() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("gas::trace".into(), vec![1.0, 0.0])],
        )
        .expect("aligned data builds");
        let error = data
            .resample_pchip(5, InterpolationSpace::Log, false)
            .expect_err("zero is invalid in log space");
        assert!(error.contains("strictly positive"));
    }

    #[test]
    fn log10_display_projection_is_shared_and_does_not_mutate_source() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("gas::fuel".into(), vec![1.0, 0.1])],
        )
        .expect("aligned data builds");
        let projected = data
            .with_log10_y_scale()
            .expect("positive values log safely");
        assert_eq!(data.column_values("gas::fuel"), Some(&[1.0, 0.1][..]));
        assert_eq!(
            projected.column_names(),
            vec!["temperature_k", "log10(gas::fuel)"]
        );
        assert_eq!(projected.column("log10(gas::fuel)"), Some(vec![0.0, -1.0]));
        let embedded = projected.to_embedded_plot_window();
        assert_eq!(embedded.values, vec!["log10(gas::fuel)"]);
        assert_eq!(embedded.y_result[(1, 0)], -1.0);
    }

    #[test]
    fn log10_display_projection_rejects_zero_without_rewriting_it() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("gas::trace".into(), vec![1.0, 0.0])],
        )
        .expect("zero is valid in the physical result");
        let error = data
            .with_log10_y_scale()
            .expect_err("zero cannot be represented on a log10 display");
        assert!(error.contains("strictly positive"));
        assert_eq!(data.column_values("gas::trace"), Some(&[1.0, 0.0][..]));
    }

    #[test]
    fn pchip_is_rejected_for_a_single_point() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0],
            vec![("gas::H2O".into(), vec![1.0])],
        )
        .expect("single point is a valid plot source");
        let error = data
            .resample_pchip(10, InterpolationSpace::Linear, false)
            .expect_err("single point must not become a fake curve");
        assert!(error.contains("at least two solved points"));
    }

    #[test]
    fn point_plot_preserves_one_exact_solver_point() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![900.0],
            vec![("gas::H2O".into(), vec![1.0])],
        )
        .expect("a point result is a valid plot source");
        let window = data.to_embedded_plot_window();
        assert_eq!(data.len(), 1);
        assert_eq!(data.x_values(), &[900.0]);
        assert_eq!(window.t_result.as_slice(), &[900.0]);
        assert_eq!(window.y_result.shape(), (1, 1));
    }

    #[test]
    fn display_filter_keeps_grid_and_values_of_selected_series() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![
                ("gas::H2O".into(), vec![0.8, 0.4]),
                ("gas::CO2".into(), vec![0.2, 0.6]),
            ],
        )
        .expect("aligned data builds");
        let filtered = data
            .retain_columns(["gas::CO2"])
            .expect("one visible series remains");
        assert_eq!(filtered.x_values(), data.x_values());
        assert_eq!(
            filtered.column_labels().collect::<Vec<_>>(),
            vec!["gas::CO2"]
        );
        assert_eq!(
            filtered.column_values("gas::CO2"),
            Some([0.2, 0.6].as_slice())
        );
        assert!(data.column_values("gas::H2O").is_some());
    }

    #[test]
    fn display_filter_rejects_an_empty_visible_set() {
        let data = EquilibriumGuiPlotData::from_columns(
            "temperature_k",
            vec![300.0, 600.0],
            vec![("gas::H2O".into(), vec![0.8, 0.4])],
        )
        .expect("aligned data builds");
        let error = data
            .retain_columns(std::iter::empty::<&str>())
            .expect_err("empty display selection must be visible to the user");
        assert!(error.contains("no selected series"));
    }

    #[test]
    fn metadata_survives_filtering_and_resampling_without_interpolating_masks() {
        let data = EquilibriumGuiPlotData::from_columns_with_metadata(
            "temperature_k",
            vec![300.0, 600.0, 900.0],
            vec![("liquid::H2O".into(), vec![1.0, 0.5, 0.0])],
            vec![vec![true, false, false]],
            EquilibriumGuiPlotUnit::Moles,
        )
        .expect("valid metadata builds");
        assert_eq!(data.y_unit(), EquilibriumGuiPlotUnit::Moles);
        assert_eq!(
            data.column_mask("liquid::H2O"),
            Some([true, false, false].as_slice())
        );

        let filtered = data
            .retain_columns(["liquid::H2O"])
            .expect("selected series remains");
        assert_eq!(filtered.y_unit(), EquilibriumGuiPlotUnit::Moles);
        assert_eq!(
            filtered.column_mask("liquid::H2O"),
            Some([true, false, false].as_slice())
        );

        let resampled = filtered
            .resample_pchip(5, InterpolationSpace::Linear, true)
            .expect("display resampling succeeds");
        assert_eq!(resampled.y_unit(), EquilibriumGuiPlotUnit::Moles);
        assert_eq!(resampled.column_mask("liquid::H2O").unwrap().len(), 5);
        assert_eq!(resampled.column_mask("liquid::H2O").unwrap()[0], true);
        assert_eq!(resampled.column_mask("liquid::H2O").unwrap()[4], false);
    }
}
