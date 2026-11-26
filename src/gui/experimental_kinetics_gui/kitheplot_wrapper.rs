use crate::Kinetics::experimental_kinetics::experiment_series2::SampledColumns;
use crate::gui::experimental_kinetics_gui::model::{PlotModel, TGAGUIError};
use eframe::egui;
use kithe_plot::controller::PlotController;
use kithe_plot::model::DataSource;
use kithe_plot::view::PlotEditorView;

// This window prepares publication/export figures. The live workflow plots are
// owned by PlotModel and rendered separately by the interactive egui plot view.

impl DataSource for SampledColumns {
    fn column(&self, name: &str) -> Option<Vec<f64>> {
        SampledColumns::column(self, name)
    }

    fn column_names(&self) -> Vec<String> {
        SampledColumns::column_names(self)
    }

    fn len(&self) -> usize {
        SampledColumns::len(self)
    }
}

#[derive(Default)]
pub struct KiThePlotWindowState {
    pub open: bool,
    controller: Option<PlotController>,
    view: Option<PlotEditorView>,
    font_setup_pending: bool,
    reposition_on_open: bool,
}

impl KiThePlotWindowState {
    pub(crate) fn close_window_on_layer(&mut self, layer: Option<egui::LayerId>) -> bool {
        if self.open
            && layer.is_some_and(|layer| layer.id == egui::Id::unique(Some("KiThe Plot Redactor")))
        {
            self.open = false;
            return true;
        }
        false
    }

    pub fn open_from_model(
        &mut self,
        ctx: &egui::Context,
        model: &mut PlotModel,
    ) -> Result<(), TGAGUIError> {
        if self.open {
            let layer = egui::LayerId::new(
                egui::Order::Foreground,
                egui::Id::unique(Some("KiThe Plot Redactor")),
            );
            ctx.move_to_top(layer);
            ctx.request_repaint();
            return Ok(());
        }

        let source = model.create_data_for_plot_window()?;
        if source.column_names().len() < 2 || source.len() == 0 {
            return Err(TGAGUIError::BindingError(
                "Not enough sampled columns to open plot redactor".to_string(),
            ));
        }

        let mut controller = PlotController::new();
        controller.load_from_data_source(&source).map_err(|err| {
            TGAGUIError::BindingError(format!("Failed to load data source: {err}"))
        })?;

        self.controller = Some(controller);
        self.view = Some(PlotEditorView::new());
        install_plot_editor_font_fallbacks(ctx);
        self.font_setup_pending = true;
        self.reposition_on_open = true;
        self.open = true;
        Ok(())
    }

    pub fn show(&mut self, ctx: &egui::Context) {
        if !self.open {
            return;
        }

        if self.font_setup_pending {
            self.font_setup_pending = false;
            ctx.request_repaint();
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

        let available = ctx.content_rect();
        let size = egui::vec2(
            available.width().min(1200.0).max(1.0),
            available.height().min(860.0).max(1.0),
        );
        let position = egui::pos2(
            (available.center().x - size.x * 0.5).max(available.left()),
            (available.center().y - size.y * 0.5).max(available.top()),
        );

        let mut window = egui::Window::new("KiThe Plot Redactor")
            .open(&mut self.open)
            .order(egui::Order::Foreground)
            .resizable(true)
            .default_pos(position)
            .default_size(size)
            .min_size([size.x.min(480.0), size.y.min(360.0)])
            .constrain_to(available);
        if self.reposition_on_open {
            window = window.current_pos(position);
            self.reposition_on_open = false;
        }
        window.show(ctx, |ui| {
            let actions = view.draw(ui, controller);
            for action in actions {
                let _ = controller.dispatch(action);
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
    use crate::Kinetics::experimental_kinetics::experiment_series_main::{
        ExperimentMeta, ExperimentThermalMode,
    };
    use crate::Kinetics::experimental_kinetics::testing_mod::VirtualTGA;

    #[test]
    fn reopening_editor_preserves_the_existing_runtime_session() {
        let data = VirtualTGA {
            time: (0..40).map(|index| index as f64).collect(),
            temperature: vec![500.0; 40],
            mass: (0..40).map(|index| 100.0 - index as f64).collect(),
        };
        let mut model = PlotModel::default();
        model
            .series
            .create_from_synthetic_data(
                &data,
                ExperimentMeta {
                    id: "editor_session".into(),
                    isothermal_temperature: Some(500.0),
                    thermal_mode: ExperimentThermalMode::Isothermal,
                    ..ExperimentMeta::default()
                },
            )
            .unwrap();
        model.set_x("editor_session", "time").unwrap();
        model.set_y("editor_session", "mass").unwrap();
        model.create_points_for_curve("editor_session").unwrap();

        let ctx = egui::Context::default();
        let mut editor = KiThePlotWindowState::default();
        ctx.begin_pass(Default::default());
        editor.open_from_model(&ctx, &mut model).unwrap();
        let mut output = ctx.end_pass();
        output.textures_delta.clear();
        let controller_before = editor.controller.as_ref().unwrap() as *const _;
        ctx.begin_pass(Default::default());
        editor.open_from_model(&ctx, &mut model).unwrap();
        let mut output = ctx.end_pass();
        output.textures_delta.clear();
        let controller_after = editor.controller.as_ref().unwrap() as *const _;

        assert_eq!(controller_before, controller_after);
        assert!(editor.open);
    }
}

fn install_plot_editor_font_fallbacks(ctx: &egui::Context) {
    let mut fonts = ctx.fonts(|fonts| fonts.definitions().clone());
    let proportional_fallbacks = fonts
        .families
        .get(&egui::FontFamily::Proportional)
        .cloned()
        .unwrap_or_default();
    for family in ["KiThePlot Sans", "KiThePlot Serif", "KiThePlot Mono"] {
        fonts.families.insert(
            egui::FontFamily::Name(family.into()),
            proportional_fallbacks.clone(),
        );
    }
    ctx.set_fonts(fonts);
}

impl PlotModel {
    pub fn create_data_for_plot_window(&mut self) -> Result<SampledColumns, TGAGUIError> {
        let (_, sampled) = self.column_samples_for_all_experiment_for_plotting(500)?;
        Ok(sampled)
    }
}
