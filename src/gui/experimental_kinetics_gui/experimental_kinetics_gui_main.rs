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
use crate::gui::experimental_kinetics_gui::controller_table::{
    ColumnManagerState, show_column_manager_window,
};
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
    pub column_manager_state: ColumnManagerState,
}

impl PlotApp {
    /// Создает новое приложение с моделью данных по умолчанию
    pub fn new() -> Self {
        Self {
            model: PlotModel::new(),
            quick_actions_state: QuickActionPanelState::default(),
            new_experiment_dialog: NewExperimentDialogState::new(),
            mathematics: Mathematics::new(),
            test_options: TestOptions::new(),
            column_manager_state: ColumnManagerState::new(),
        }
    }

    /// Отображает окно приложения
    pub fn show(&mut self, ctx: &egui::Context, open: &mut bool) {
        egui::Window::new("📊 Experimental Kinetics")
            .open(open)
            .default_size([1000.0, 700.0])
            .min_size([600.0, 400.0])
            .show(ctx, |ui| {
                self.render_ui(ui);
            });
    }

    /// Отрисовка пользовательского интерфейса
    fn render_ui(&mut self, ui: &mut egui::Ui) {
        // Top drop-down menus (File Manager, Math, Kinetic Methods, ...)
        TopDropDownMenues::top_menus(
            ui,
            &mut self.model,
            &mut self.mathematics,
            &mut self.new_experiment_dialog,
            &mut self.test_options,
        );
        self.new_experiment_dialog
            .show_new_experiment_dialogue(ui.ctx(), &mut self.model);
        self.new_experiment_dialog
            .show_manage_plot_dialog(ui.ctx(), &mut self.model);
        self.new_experiment_dialog
            .show_save_series_dialog(ui.ctx(), &mut self.model);
        show_column_manager_window(
            ui.ctx(),
            &mut self.model,
            &mut self.column_manager_state,
            &mut self.new_experiment_dialog.column_manager_open,
        );

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
                WrightPanelControllers::quick_action_panel(
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
    }
}
