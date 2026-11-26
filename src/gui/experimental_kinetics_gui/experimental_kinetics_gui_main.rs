//! Демонстрационное приложение интерактивного графика с архитектурой MVC
//!
//! ## Общие сведения об архитектуре MVC
//!
//! MVC (Model-View-Controller) - это архитектурный паттерн проектирования, который разделяет
//! приложение на три основных компонента:
//!
//! 1. **Model (Модель)** - отвечает за данные и за научно-инженерную либо бизнес-логику приложения
//! 2. **View (Представление)** - отвечает за отображение данных пользователю
//! 3. **Controller (Контроллер)** - обрабатывает взаимодействие пользователя с приложением
//!
//! ## Реализация MVC в данном приложении
//!
//! ### Model (Модель) - файл model.rs
//! В данном приложении модель представлена структурой PlotModel, которая содержит:
//! - Данные графиков (вектор GraphData)
//! - Состояние выделения (прямоугольник SelectionRect)
//! - Параметры отображения (диапазоны по осям X и Y)
//! - Состояние взаимодействия (флаги панорамирования, выделения)
//!
//! Модель отвечает за:
//! - Хранение всех данных графиков
//! - Вычисление точек графиков
//! - Управление состоянием выделения
//! - Реализацию операций масштабирования и панорамирования
//! - Поиск ближайших графиков
//!
//! ### View (Представление) - файл view.rs
//! Представление отвечает за визуализацию данных и взаимодействие с пользователем:
//! - Отрисовка графиков с использованием библиотеки egui_plot
//! - Отображение прямоугольника выделения
//! - Форматирование осей и сетки
//! - Отображение информационной панели с инструкциями
//! - Обработка визуальных эффектов (выделение, подсветка)
//!
//! ### Controller (Контроллер) - файл controller.rs
//! Контроллер обрабатывает все взаимодействия пользователя с приложением:
//! - Обработка перетаскивания мыши (панорамирование и выделение)
//! - Обработка кликов мыши (выбор графиков)
//! - Обработка прокрутки колеса мыши (масштабирование)
//! - Обработка клавиш стрелок (панорамирование)
//! - Обработка двойного клика (сброс вида)
//! - Обработка кнопок интерфейса (сброс, очистка и т.д.)
//! - Обработка элементов управления параметрами графиков
//!
//! ## Почему такая структура?
//!
//! Такое разделение позволяет:
//! 1. **Упростить сопровождение кода** - каждая часть отвечает за свою функциональность
//! 2. **Повысить переиспользуемость** - компоненты можно использовать в других приложениях
//! 3. **Облегчить тестирование** - каждую часть можно тестировать отдельно
//! 4. **Упростить модификацию** - можно изменять один компонент, не затрагивая другие
//!
//! Например, если мы захотим изменить способ отображения графиков,
//! нам нужно будет изменить только View, не трогая Model и Controller.
//! Аналогично, если мы захотим добавить новые способы взаимодействия с графиком,
//! мы можем это сделать в Controller, не изменяя Model и View.

/// Импортируем контроллер графика
use crate::gui::experimental_kinetics_gui::controller::PlotController;
/// Импорт верхних выпадающих меню
use crate::gui::experimental_kinetics_gui::controller_buttons_and_panels::{
    NewExperimentDialogState, QuickActionPanelState, TopDropDownMenues, WrightPanelControllers,
};
use crate::gui::experimental_kinetics_gui::controller_filters::Mathematics;
use crate::gui::experimental_kinetics_gui::controller_kinetics::DirectProblemDialogState;
use crate::gui::experimental_kinetics_gui::controller_methods::KineticMethodsWindowState;
use crate::gui::experimental_kinetics_gui::controller_table::{
    ColumnManagerState, show_column_manager_window,
};
use crate::gui::experimental_kinetics_gui::kitheplot_wrapper::KiThePlotWindowState;
/// Импортируем модель данных графика
use crate::gui::experimental_kinetics_gui::model::PlotModel;
/// Импортируем TestOptions для управления настройками
use crate::gui::experimental_kinetics_gui::test_options::TestOptions;
/// Импортируем представление графика
use crate::gui::experimental_kinetics_gui::view::PlotView;
/// Импортируем необходимые компоненты из библиотеки egui
use eframe::egui;

/// Основная структура приложения
///
/// Содержит модель данных графика
pub struct PlotApp {
    /// Модель данных графика
    pub model: PlotModel,
    pub quick_actions_state: QuickActionPanelState,
    pub new_experiment_dialog: NewExperimentDialogState,
    pub mathematics: Mathematics,
    pub test_options: TestOptions,
    pub kinetic_methods_state: KineticMethodsWindowState,
    pub column_manager_state: ColumnManagerState,
    pub kithe_plot_window: KiThePlotWindowState,
    pub direct_problem_state: DirectProblemDialogState,
    pub(crate) keyboard_input_focused_last_frame: bool,
}

impl PlotApp {
    fn handle_escape(&mut self, ctx: &egui::Context) {
        if self.keyboard_input_focused_last_frame
            || ctx.egui_wants_keyboard_input()
            || egui::Popup::is_any_open(ctx)
        {
            return;
        }
        let escape =
            ctx.input_mut(|input| input.consume_key(egui::Modifiers::NONE, egui::Key::Escape));
        if !escape {
            return;
        }

        let top_layer = ctx.top_layer_id();
        let closed = self.new_experiment_dialog.close_window_on_layer(top_layer)
            || self.column_manager_state.close_window_on_layer(top_layer)
            || self.mathematics.close_window_on_layer(top_layer)
            || self.kinetic_methods_state.close_window_on_layer(top_layer)
            || self.test_options.close_window_on_layer(top_layer)
            || self.direct_problem_state.close_window_on_layer(top_layer)
            || self.kithe_plot_window.close_window_on_layer(top_layer);
        if closed {
            ctx.request_repaint();
        }
    }

    /// Создает новое приложение с моделью данных по умолчанию
    pub fn new() -> Self {
        Self {
            model: PlotModel::new(),
            quick_actions_state: QuickActionPanelState::default(),
            new_experiment_dialog: NewExperimentDialogState::new(),
            mathematics: Mathematics::new(),
            test_options: TestOptions::new(),
            kinetic_methods_state: KineticMethodsWindowState::default(),
            column_manager_state: ColumnManagerState::new(),
            kithe_plot_window: KiThePlotWindowState::default(),
            direct_problem_state: DirectProblemDialogState::new(),
            keyboard_input_focused_last_frame: false,
        }
    }

    /// Отображает окно приложения
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("📊 Experimental Kinetics")
            .open(open)
            .resizable(true)
            .default_size([1000.0, 700.0])
            .min_size([600.0, 400.0])
            .show(ctx, |ui| {
                self.render_ui(ui);
            });
        self.kithe_plot_window.show(ctx);
    }

    /// Отрисовка пользовательского интерфейса
    fn render_ui(&mut self, ui: &mut egui::Ui) {
        self.handle_escape(ui.ctx());
        let document_generation = self.new_experiment_dialog.document_generation();
        // Top drop-down menus (File Manager, Math, Kinetic Methods, ...)
        TopDropDownMenues::top_menus(
            ui,
            &mut self.model,
            &mut self.mathematics,
            &mut self.new_experiment_dialog,
            &mut self.test_options,
            &mut self.kithe_plot_window,
            &mut self.kinetic_methods_state,
            &mut self.direct_problem_state,
        );
        self.new_experiment_dialog
            .show_new_experiment_dialogue(ui.ctx(), &mut self.model);
        self.new_experiment_dialog
            .show_manage_plot_dialog(ui.ctx(), &mut self.model);
        self.new_experiment_dialog
            .show_save_series_dialog(ui.ctx(), &mut self.model);
        self.new_experiment_dialog
            .show_unsaved_changes_dialog(ui.ctx(), &mut self.model);
        self.test_options.show_load_testing_data_confirmation(
            ui.ctx(),
            &mut self.model,
            &mut self.new_experiment_dialog,
        );
        if self.new_experiment_dialog.document_generation() != document_generation {
            self.column_manager_state.clear_pending_confirmations();
            self.kinetic_methods_state.fit_sessions.clear();
            self.kinetic_methods_state.fit_pending_batch = None;
            self.kinetic_methods_state.fit_model_result = None;
            self.kinetic_methods_state.fit_model_status.clear();
            self.kinetic_methods_state.fit_model_experiment_id.clear();
            self.kinetic_methods_state.fit_model_x_col.clear();
            self.kinetic_methods_state.fit_model_y_col.clear();
            self.kinetic_methods_state
                .fit_reconstructed_experiment_id
                .clear();
            self.kinetic_methods_state
                .isoconversional_selected_experiments
                .clear();
        }
        let workspace_cleared = show_column_manager_window(
            ui.ctx(),
            &mut self.model,
            &mut self.column_manager_state,
            &mut self.kinetic_methods_state,
            &mut self.new_experiment_dialog.column_manager_open,
        );
        if workspace_cleared {
            self.quick_actions_state = QuickActionPanelState::default();
            self.new_experiment_dialog = NewExperimentDialogState::new();
            self.mathematics = Mathematics::new();
            self.test_options = TestOptions::new();
            self.kinetic_methods_state = KineticMethodsWindowState::default();
            self.column_manager_state = ColumnManagerState::default();
            self.kithe_plot_window = KiThePlotWindowState::default();
            self.direct_problem_state = DirectProblemDialogState::new();
            self.keyboard_input_focused_last_frame = false;
        }

        // Show settings window if open
        if self.test_options.is_settings_window_open() {
            self.test_options
                .show_settings_window(ui.ctx(), &mut self.model);
        }
        // Show help window if open (managed inside TestOptions)
        self.test_options.show_help_window_ui(ui.ctx());
        // Show synthetic data window if requested
        self.test_options
            .show_synthetic_data_window(ui.ctx(), &mut self.model);
        self.test_options
            .show_logs_window(ui.ctx(), &mut self.model);
        // Создаем две колонки для разделения графика и элементов управления
        ui.columns(2, |columns| {
            // Левая колонка: Область отображения графика
            columns[0].vertical(|ui| {
                // Заголовок графика
                ui.heading("Data Plot");
                // Разделитель
                ui.separator();

                // Отображаем график с помощью представления
                let plot_response = PlotView::render_plot(ui, &mut self.model, |plot_ui, model| {
                    // Обрабатываем взаимодействия с графиком через контроллер
                    PlotController::handle_plot_interactions(plot_ui, model);
                });

                // Обрабатываем контекстное меню графика
                plot_response.response.context_menu(|ui| {
                    // Кнопка сброса вида
                    if ui.button("Reset View").clicked() {
                        self.model.reset_view();
                        ui.ctx().request_repaint();
                    }
                    // Кнопка очистки всех выделений
                    if ui.button("Clear All Selections").clicked() {
                        self.model.clear_selection();
                        self.model.clear_selection_rect();
                    }
                });
            });

            // Правая колонка: Элементы управления и информация
            columns[1].vertical(|ui| {
                PlotView::render_info_panel(ui, &self.model);
                ui.separator();
                let _ = WrightPanelControllers::quick_action_panel(
                    ui,
                    &mut self.model,
                    &mut self.quick_actions_state,
                    &mut self.new_experiment_dialog,
                );
                WrightPanelControllers::handle_ui_interactions(
                    ui,
                    &mut self.model,
                    &mut self.new_experiment_dialog,
                );
            });
        });

        self.direct_problem_state.show(ui.ctx(), &mut self.model);
        self.keyboard_input_focused_last_frame = ui.ctx().egui_wants_keyboard_input();
    }
}
