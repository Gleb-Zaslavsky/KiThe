use eframe::egui;
use crate::Thermodynamics::thermo_lib_api::ThermoData;

pub struct AllLibsGui {
    thermo_data: ThermoData,
    show_differences: bool,
}

impl Default for AllLibsGui {
    fn default() -> Self {
        Self::new()
    }
}

impl AllLibsGui {
    pub fn new() -> Self {
        Self {
            thermo_data: ThermoData::new(),
            show_differences: false,
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Library Review")
            .open(open)
            .default_size([1200.0, 1800.0])
            .show(ctx, |ui| {
                self.ui(ui);
            });
    }

    fn ui(&mut self, ui: &mut egui::Ui) {
        ui.heading("Library Review");
        ui.separator();
        
        ui.horizontal(|ui| {
            // Left panel with both scroll areas
            ui.vertical(|ui| {
                ui.set_width(200.0);
                ui.set_min_height(1000.0);
                
                ui.strong("Substances addresses");
                egui::ScrollArea::vertical()
                    .id_salt("all_keys_subs")
                    .max_height(400.0)
                    .show(ui, |ui| {
                        for (library, substance) in &self.thermo_data.VecOfSubsAdresses {
                            if !self.show_differences || !self.exists_in_lib_thermo_data(library, substance) {
                                ui.colored_label(egui::Color32::BLACK, format!("{} - {}", substance, library));
                            }
                        }
                    });
                 ui.separator();
                ui.add_space(10.0);
                 ui.separator();
                ui.strong("Substance-Library Pairs in the database");
                egui::ScrollArea::vertical()
                    .id_salt("substance_libs")
                    .max_height(400.0)
                    .show(ui, |ui| {
                        for (library, substances) in &self.thermo_data.LibThermoData {
                            for substance in substances.keys() {
                                if !self.show_differences || !self.exists_in_vec_addresses(library, substance) {
                                    ui.colored_label(egui::Color32::BLACK, format!("{} - {}", substance, library));
                                }
                            }
                        }
                    });
                     ui.separator();
            });
            
            ui.separator();
        });
        
        ui.separator();
        if ui.button("Compare databases").clicked() {
            self.show_differences = !self.show_differences;
        }
        
        if self.show_differences {
            ui.label("Showing only differences between databases");
        } else {
            ui.label("Showing all data");
        }
    }
    
    fn exists_in_lib_thermo_data(&self, library: &str, substance: &str) -> bool {
        self.thermo_data.LibThermoData
            .get(library)
            .map_or(false, |substances| substances.contains_key(substance))
    }
    
    fn exists_in_vec_addresses(&self, library: &str, substance: &str) -> bool {
        self.thermo_data.VecOfSubsAdresses
            .iter()
            .any(|(lib, sub)| lib == library && sub == substance)
    }
}