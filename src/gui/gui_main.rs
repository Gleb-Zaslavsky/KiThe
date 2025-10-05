use crate::gui::combustion;
use crate::gui::gui_main::egui::IconData;
use crate::gui::gui_main::egui::Image;
use crate::gui::kinetics;
use eframe::egui;
use prettytable::color;
pub fn gui_main() -> Result<(), eframe::Error> {
    // Try to load icon from assets, fallback to programmatic icon
    let icon = match std::fs::read("src/assets/icon.png") {
        Ok(icon_bytes) => match eframe::icon_data::from_png_bytes(&icon_bytes) {
            Ok(icon_data) => icon_data,
            Err(_) => create_programmatic_icon(),
        },
        Err(_) => create_programmatic_icon(),
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([800.0, 600.0])
            .with_title("Main menu ")
            .with_icon(icon),
        ..Default::default()
    };

    eframe::run_native(
        "Chemical Analysis Suite",
        options,
        Box::new(|_cc| Ok(Box::new(MainApp::default()))),
    )
}

fn create_programmatic_icon() -> IconData {
    let size = 32;
    let mut rgba = vec![0u8; size * size * 4];
    for y in 0..size {
        for x in 0..size {
            let idx = (y * size + x) * 4;
            let center_x = size as f32 / 2.0;
            let center_y = size as f32 / 2.0;
            let dist = ((x as f32 - center_x).powi(2) + (y as f32 - center_y).powi(2)).sqrt();
            if dist < 12.0 {
                rgba[idx] = 255;
                rgba[idx + 1] = (200.0 * (1.0 - dist / 12.0)) as u8;
                rgba[idx + 2] = 0;
                rgba[idx + 3] = 255;
            }
        }
    }
    IconData {
        rgba,
        width: 32,
        height: 32,
    }
}

#[derive(Default)]
struct MainApp {
    kinetics_open: bool,
    kinetics_app: Option<kinetics::KineticsApp>,
    combustion_open: bool,
    combustion_app: Option<combustion::CombustionApp>,
}

impl eframe::App for MainApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Set button text color and size to be more visible
        ctx.style_mut(|style| {
            style.visuals.widgets.inactive.fg_stroke.color = egui::Color32::BLACK;
            style.visuals.widgets.hovered.fg_stroke.color = egui::Color32::BLACK;
            style.visuals.widgets.active.fg_stroke.color = egui::Color32::WHITE;
            style.text_styles.insert(
                egui::TextStyle::Button,
                egui::FontId::new(16.0, egui::FontFamily::Proportional),
            );
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.vertical_centered(|ui| {
                ui.add_space(50.0);
                // Main title
                ui.heading("Heat and mass transfer, chemical reactors \ncombustion and macrokinetics suite");
                ui.add_space(30.0);
                // Menu buttons in a grid layout
                ui.horizontal(|ui| {
                    ui.add_space(50.0);
                    ui.vertical(|ui| {
                        // First row of buttons
                        ui.horizontal(|ui| {
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🧪 Kinetics NOT READY!")).clicked() {
                                self.kinetics_open = true;
                                if self.kinetics_app.is_none() {
                                    self.kinetics_app = Some(kinetics::KineticsApp::new());
                                }
                            }
                            ui.add_space(20.0);
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🌡️ Thermodynamics NOT READY!")).clicked() {
                                println!("Thermodynamics module clicked");
                            }
                        });
                        ui.add_space(20.0);
                        // Second row of buttons
                        ui.horizontal(|ui| {
                            if ui.add_sized([200.0, 80.0], egui::Button::new("⚖️ Mass Transfer NOT READY!")).clicked() {
                                println!("Mass Transfer module clicked");
                            }
                            ui.add_space(20.0);
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🔥 Gas-phase combustuion/steady state plug flow")).clicked() {
                                self.combustion_open = true;
                                if self.combustion_app.is_none() {
                                    self.combustion_app = Some(combustion::CombustionApp::new());
                                }
                            }
                        });
                        ui.add_space(20.0);
                        // Third row of buttons
                        ui.horizontal(|ui| {
                            if ui.add_sized([200.0, 80.0], egui::Button::new("📊 Data Analysis")).clicked() {
                                println!("Data Analysis module clicked");
                            }
                            ui.add_space(20.0);
                            if ui.add_sized([200.0, 80.0], egui::Button::new("⚙️ Settings")).clicked() {
                                println!("Settings module clicked");
                            }
                        });
                    });
                });
                ui.add_space(50.0);
                // Footer information
                ui.separator();
                ui.add_space(20.0);
                // Try to load and display logo
                if let Ok(logo_bytes) = std::fs::read("src/assets/logo.png") {
                    if let Ok(color_image) = eframe::icon_data::from_png_bytes(&logo_bytes) {
                        let egui_image = egui::ColorImage::from_rgba_unmultiplied(
                            [color_image.width as usize, color_image.height as usize],
                            &color_image.rgba
                        );
                        let texture = ctx.load_texture("logo", egui_image, egui::TextureOptions::default());
                        ui.add(egui::Image::new(&texture).max_width(200.0));
                        ui.add_space(10.0);
                    } else {
                        ui.label("Logo failed to load");
                    }
                } else {
                    ui.label("Logo not found");
                }
                ui.label("developed by Gleb E. Zaslavsky 2024-2025 (c)");
            }) ;
        });

        // Show kinetics window if opened
        if self.kinetics_open {
            if let Some(kinetics_app) = &mut self.kinetics_app {
                kinetics_app.show(ctx, &mut self.kinetics_open);
            }
        }

        // Show combustion window if opened
        if self.combustion_open {
            if let Some(combustion_app) = &mut self.combustion_app {
                combustion_app.show(ctx, &mut self.combustion_open);
            }
        }
    }
}
