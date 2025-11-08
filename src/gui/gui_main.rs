use crate::gui::all_libs_gui;
use crate::gui::combustion;
use crate::gui::gui_main::egui::IconData;
use crate::gui::kinetics_gui;
use crate::gui::settings_gui;
use crate::gui::thermochemistry_gui;
use crate::gui::transport_gui;
use eframe::CreationContext;
use eframe::egui;
use egui::ColorImage;
use egui::TextureHandle;
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
    /*
    eframe::run_native(
        "Chemical Analysis Suite",
        options,
        Box::new(|cc| Ok(Box::new(MainApp::new(cc)))),
    )
    */
    eframe::run_native(
        "Chemical Analysis Suite",
        options,
        Box::new(|cc: &eframe::CreationContext<'_>| Ok(Box::new(MainApp::new(cc)))),
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
    kinetics_app: Option<kinetics_gui::KineticsApp>,
    combustion_open: bool,
    combustion_app: Option<combustion::CombustionApp>,
    thermochemistry_open: bool,
    thermochemistry_app: Option<thermochemistry_gui::ThermochemistryApp>,
    transport_open: bool,
    transport_app: Option<transport_gui::TransportApp>,
    settings_open: bool,
    settings_app: Option<settings_gui::SettingsGui>,
    all_libs_open: bool,
    all_libs_app: Option<all_libs_gui::AllLibsGui>,
    logo_texture: Option<egui::TextureHandle>,
}
impl MainApp {
    pub fn new(cc: &CreationContext<'_>) -> Self {
        let ctx = &cc.egui_ctx;

        Self::setup_visuals(ctx);
        Self::setup_fonts(ctx);
        let logo_texture = Self::load_logo_texture(ctx);

        Self {
            kinetics_open: false,
            kinetics_app: None,
            combustion_open: false,
            combustion_app: None,
            thermochemistry_open: false,
            thermochemistry_app: None,
            transport_open: false,
            transport_app: None,
            settings_open: false,
            settings_app: None,
            all_libs_open: false,
            all_libs_app: None,
            logo_texture,
        }
    }

    fn setup_visuals(ctx: &egui::Context) {
        use egui::Visuals;
        ctx.set_visuals(Visuals::dark()); // or .light()
    }

    fn setup_fonts(ctx: &egui::Context) {
        use egui::FontDefinitions;
        let mut fonts = FontDefinitions::default();
        // You can load custom fonts here if desired.
        ctx.set_fonts(fonts);
    }

    fn load_logo_texture(ctx: &egui::Context) -> Option<TextureHandle> {
        if let Ok(logo_bytes) = std::fs::read("src/assets/logo.png") {
            if let Ok(color_image) = eframe::icon_data::from_png_bytes(&logo_bytes) {
                let egui_image = ColorImage::from_rgba_unmultiplied(
                    [color_image.width as usize, color_image.height as usize],
                    &color_image.rgba,
                );
                return Some(ctx.load_texture(
                    "kithe_logo",
                    egui_image,
                    egui::TextureOptions::default(),
                ));
            }
        }
        None
    }
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
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🧪 Kinetics!")).clicked() {
                                self.kinetics_open = true;
                                if self.kinetics_app.is_none() {
                                    self.kinetics_app = Some(kinetics_gui::KineticsApp::new());
                                }
                            }
                            ui.add_space(20.0);
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🌡️ Thermodynamical properties!")).clicked() {
                                self.thermochemistry_open = true;
                                if self.thermochemistry_app.is_none() {
                                    self.thermochemistry_app = Some(thermochemistry_gui::ThermochemistryApp::new());
                                }
                            }
                        });
                        ui.add_space(20.0);
                        // Second row of buttons
                        ui.horizontal(|ui| {
                            if ui.add_sized([200.0, 80.0], egui::Button::new("🚚 Transport properties!")).clicked() {
                                self.transport_open = true;
                                if self.transport_app.is_none() {
                                    self.transport_app = Some(transport_gui::TransportApp::new());
                                }
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
                            if ui.add_sized([200.0, 80.0], egui::Button::new("📚 Library Review")).clicked() {
                                self.all_libs_open = true;
                                if self.all_libs_app.is_none() {
                                    self.all_libs_app = Some(all_libs_gui::AllLibsGui::new());
                                }
                            }
                            ui.add_space(20.0);
                            if ui.add_sized([200.0, 80.0], egui::Button::new("⚙️ Settings")).clicked() {
                                self.settings_open = true;
                                if self.settings_app.is_none() {
                                    self.settings_app = Some(settings_gui::SettingsGui::new());
                                }
                            }
                        });
                    });
                });
                ui.add_space(50.0);
                // Footer information
                ui.separator();
                ui.add_space(20.0);
                // Load and display logo (cached)
                if self.logo_texture.is_none() {
                    if let Ok(logo_bytes) = std::fs::read("src/assets/logo.png") {
                        if let Ok(color_image) = eframe::icon_data::from_png_bytes(&logo_bytes) {
                            let egui_image = egui::ColorImage::from_rgba_unmultiplied(
                                [color_image.width as usize, color_image.height as usize],
                                &color_image.rgba
                            );
                            self.logo_texture = Some(ctx.load_texture("logo", egui_image, egui::TextureOptions::default()));
                        }
                    }
                }
                if let Some(texture) = &self.logo_texture {
                    ui.add(egui::Image::new(texture).max_width(200.0));
                    ui.add_space(10.0);
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

        // Show thermochemistry window if opened
        if self.thermochemistry_open {
            if let Some(thermochemistry_app) = &mut self.thermochemistry_app {
                thermochemistry_app.show(ctx, &mut self.thermochemistry_open);
            }
        }

        // Show transport window if opened
        if self.transport_open {
            if let Some(transport_app) = &mut self.transport_app {
                transport_app.show(ctx, &mut self.transport_open);
            }
        }

        // Show settings window if opened
        if self.settings_open {
            if let Some(settings_app) = &mut self.settings_app {
                settings_app.show(ctx, &mut self.settings_open);
            }
        }

        // Show all libs window if opened
        if self.all_libs_open {
            if let Some(all_libs_app) = &mut self.all_libs_app {
                all_libs_app.show(ctx, &mut self.all_libs_open);
            }
        }
    }
}
