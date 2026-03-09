use crate::Kinetics::experimental_kinetics::experiment_series2::SampledColumns;
use crate::gui::experimental_kinetics_gui::model::{PlotModel, TGAGUIError};
use eframe::egui;
use kithe_plot::controller::PlotController;
use kithe_plot::model::DataSource;
use kithe_plot::view::PlotEditorView;

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
}

impl KiThePlotWindowState {
    const VIEWPORT_ID: &'static str = "kithe_plot_redactor_viewport";

    pub fn open_from_model(&mut self, model: &mut PlotModel) -> Result<(), TGAGUIError> {
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
            .with_title("KiThe Plot Redactor")
            .with_inner_size([1200.0, 860.0]);
        let open = &mut self.open;

        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            if ctx.input(|i| i.viewport().close_requested()) {
                *open = false;
                return;
            }

            match class {
                egui::ViewportClass::Embedded => {
                    egui::Window::new("KiThe Plot Redactor")
                        .open(open)
                        .default_size([1200.0, 860.0])
                        .show(ctx, |_ui| {
                            let actions = view.draw(ctx, controller);
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

impl PlotModel {
    pub fn create_data_for_plot_window(&mut self) -> Result<SampledColumns, TGAGUIError> {
        let (_, sampled) = self.column_samples_for_all_experiment_for_plotting(500)?;
        Ok(sampled)
    }
}
