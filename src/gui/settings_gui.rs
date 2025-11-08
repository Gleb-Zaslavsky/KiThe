use crate::settings::Settings;
use eframe::egui;
use std::collections::HashMap;

#[derive(Default)]
pub struct SettingsGui {
    settings: Settings,
    selected_library: String,
    new_path: String,
    status_message: String,
    show_file_dialog: bool,
}

impl SettingsGui {
    pub fn new() -> Self {
        Self {
            settings: Settings::new(),
            selected_library: "Substance Base".to_string(),
            new_path: String::new(),
            status_message: String::new(),
            show_file_dialog: false,
        }
    }

    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("Library Settings")
            .open(open)
            .default_size([400.0, 300.0])
            .show(ctx, |ui| {
                self.ui(ui);
            });
    }

    fn ui(&mut self, ui: &mut egui::Ui) {
        ui.heading("Library Version Manager");
        ui.separator();

        // Current configuration display
        ui.label("Current Library Configuration:");
        egui::Grid::new("library_grid")
            .num_columns(2)
            .spacing([40.0, 4.0])
            .striped(true)
            .show(ui, |ui| {
                for library in self.settings.get_available_libraries() {
                    ui.label(library);
                    ui.label(self.settings.get_library_version(library).unwrap());
                    ui.end_row();
                }
            });

        ui.separator();

        // Library selection and update
        ui.horizontal(|ui| {
            ui.label("Select Library:");
            egui::ComboBox::from_label("")
                .selected_text(&self.selected_library)
                .show_ui(ui, |ui| {
                    for library in self.settings.get_available_libraries() {
                        ui.selectable_value(&mut self.selected_library, library.clone(), library);
                    }
                });
        });

        ui.horizontal(|ui| {
            ui.label("New Path:");
            ui.text_edit_singleline(&mut self.new_path);

            if ui.button("Browse...").clicked() {
                self.show_file_dialog = true;
            }
        });

        ui.horizontal(|ui| {
            if ui.button("Update Library").clicked() {
                self.update_library();
            }

            if ui.button("Reset to Defaults").clicked() {
                self.reset_to_defaults();
            }
        });

        // Status message
        if !self.status_message.is_empty() {
            ui.separator();
            ui.colored_label(
                if self.status_message.starts_with("✓") {
                    egui::Color32::GREEN
                } else {
                    egui::Color32::RED
                },
                &self.status_message,
            );
        }

        // File dialog simulation (in real implementation, use rfd or similar)
        if self.show_file_dialog {
            self.show_file_dialog_ui(ui);
        }
    }

    fn update_library(&mut self) {
        if self.new_path.is_empty() {
            self.status_message = "✗ Please enter a path".to_string();
            return;
        }

        match self
            .settings
            .set_library_version(&self.selected_library, &self.new_path)
        {
            Ok(_) => {
                self.status_message = format!("✓ Updated {} successfully", self.selected_library);
                self.new_path.clear();
            }
            Err(e) => {
                self.status_message = format!("✗ Failed to update: {}", e);
            }
        }
    }

    fn reset_to_defaults(&mut self) {
        match self.settings.reset_to_defaults() {
            Ok(_) => {
                self.status_message = "✓ Reset to defaults successfully".to_string();
                self.new_path.clear();
            }
            Err(e) => {
                self.status_message = format!("✗ Failed to reset: {}", e);
            }
        }
    }

    fn show_file_dialog_ui(&mut self, ui: &mut egui::Ui) {
        ui.separator();
        ui.label("File Dialog (Simulation):");
        ui.label("In a real implementation, use 'rfd' crate for native file dialogs");

        // Common file suggestions
        let suggestions = vec![
            "substance_base_v2.json",
            "substance_base_v3.json",
            "all_keys_substance.json",
            "all_keys_substance_v2.json",
            "Reactbase.json",
            "Reactbase_v2.json",
            "dict_reaction.json",
            "dict_reaction_v2.json",
        ];

        ui.label("Quick select:");
        for suggestion in suggestions {
            if ui.small_button(suggestion).clicked() {
                self.new_path = suggestion.to_string();
                self.show_file_dialog = false;
            }
        }

        if ui.button("Cancel").clicked() {
            self.show_file_dialog = false;
        }
    }

    // Batch update method for advanced users
    pub fn batch_update(&mut self, updates: HashMap<&str, &str>) -> Result<(), String> {
        self.settings.update_multiple_libraries(updates)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_settings_gui_new() {
        let gui = SettingsGui::new();
        assert_eq!(gui.selected_library, "Substance Base");
        assert!(gui.new_path.is_empty());
        assert!(gui.status_message.is_empty());
    }
}
