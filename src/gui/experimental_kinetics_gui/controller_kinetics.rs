use crate::gui::experimental_kinetics_gui::model::PlotModel;

// Group: Direct Problem (stubs)
pub struct DirectProblem;

impl DirectProblem {
    pub fn show(ui: &mut egui::Ui, _model: &mut PlotModel) {
        if ui.button("Solve IVP").clicked() {
            println!("Stub: Solve IVP");
        }

        if ui.button("Solve BVP").clicked() {
            println!("Stub: Solve BVP");
        }
    }
}
