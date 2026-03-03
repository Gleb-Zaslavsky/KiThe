use crate::gui::experimental_kinetics_gui::model::{PlotCurve, PlotModel};
use egui::{Color32, Ui};
use egui_plot::{Line, Plot, PlotBounds, PlotPoints, PlotUi};

/// Представление (View) в архитектуре MVC
///
/// Отвечает за отображение графиков и пользовательского интерфейса.
/// Содержит методы для рендеринга графиков, информационной панели
/// и вспомогательных функций для форматирования осей.
pub struct PlotView;

impl PlotView {
    /// Отображает основную область графика
    ///
    /// Создает и настраивает виджет графика, отображает все графики,
    /// прямоугольник выделения и обрабатывает пользовательские взаимодействия.
    ///
    /// # Параметры
    /// * `ui` - изменяемая ссылка на пользовательский интерфейс egui
    /// * `model` - изменяемая ссылка на модель данных графика
    /// * `on_plot_ui` - функция обратного вызова для обработки взаимодействий с графиком
    ///
    /// # Возвращает
    /// Результат отображения графика (PlotResponse)
    pub fn render_plot<F>(
        ui: &mut Ui,
        model: &mut PlotModel,
        mut on_plot_ui: F,
    ) -> egui_plot::PlotResponse<()>
    where
        F: FnMut(&mut PlotUi, &mut PlotModel),
    {
        // Отключаем встроенные функции перетаскивания/масштабирования,
        // так как основное приложение управляет видом через модель
        let plot = Plot::new("interactive_plot")
            .legend(egui_plot::Legend::default()) // Отображаем легенду графиков
            .allow_drag(false) // Отключаем встроенное перетаскивание
            .allow_zoom(true) // Разрешаем масштабирование
            .allow_boxed_zoom(true) // Разрешаем масштабирование прямоугольником
            .x_grid_spacer(Self::aligned_x_grid_spacer) // Используем собственный форматтер сетки
            .show_axes(true) // Отображаем оси
            .show_grid(true) // Отображаем сетку
            .height(400.0); // Устанавливаем высоту графика

        plot.show(ui, |plot_ui| {
            {
                // Устанавливаем границы графика из модели, чтобы метки осей
                // соответствовали ручному панорамированию/масштабированию
                let bounds = PlotBounds::from_min_max(
                    [
                        model.interaction.view_range.0,
                        model.interaction.view_y_range.0,
                    ], // Минимальные значения [x, y]
                    [
                        model.interaction.view_range.1,
                        model.interaction.view_y_range.1,
                    ], // Максимальные значения [x, y]
                );
                // Применяем границы, если они корректны
                if bounds.is_valid() {
                    plot_ui.set_plot_bounds(bounds);
                }

                // Render each curve
                for curve in &model.plots {
                    Self::render_curve(plot_ui, curve);
                }

                // Отображаем прямоугольник выделения (если есть активное выделение)
                if let Some(rect) = &model.interaction.selection_rect {
                    // Проверяем, активен ли прямоугольник выделения
                    if rect.is_active() {
                        // Нормализуем прямоугольник для корректного отображения
                        let normalized = rect.normalized();
                        // Создаем точки для отображения прямоугольника
                        let rect_points = PlotPoints::from(vec![
                            [normalized.start[0], normalized.start[1]], // Левый нижний угол
                            [normalized.end[0], normalized.start[1]],   // Правый нижний угол
                            [normalized.end[0], normalized.end[1]],     // Правый верхний угол
                            [normalized.start[0], normalized.end[1]],   // Левый верхний угол
                            [normalized.start[0], normalized.start[1]], // Замыкаем контур
                        ]);
                        // Создаем линию для отображения прямоугольника выделения
                        let line = egui_plot::Line::new("", rect_points)
                            .color(Color32::RED) // Красный цвет
                            .width(1.0) // Толщина линии
                            .fill(0.1); // Прозрачная заливка
                        plot_ui.line(line);
                    }
                }

                // Обеспечиваем правильный диапазон по Y, добавляя почти прозрачный
                // прямоугольник, охватывающий текущие диапазоны по X и Y.
                // Это гарантирует, что egui_plot включит границы по Y в отображение.
                let (y_min, y_max) = model.interaction.view_y_range;
                // Проверяем корректность границ по Y
                if y_min.is_finite() && y_max.is_finite() && y_max > y_min {
                    let x0 = model.interaction.view_range.0; // Левая граница по X
                    let x1 = model.interaction.view_range.1; // Правая граница по X
                    // Создаем точки для вспомогательного прямоугольника
                    let rect_points = PlotPoints::from(vec![
                        [x0, y_min], // Левый нижний угол
                        [x1, y_min], // Правый нижний угол
                        [x1, y_max], // Правый верхний угол
                        [x0, y_max], // Левый верхний угол
                        [x0, y_min], // Замыкаем контур
                    ]);
                    // Создаем почти прозрачную линию черного цвета с минимальным альфа-каналом,
                    // чтобы рендерер по-прежнему учитывал точки при вычислении диапазонов осей
                    let helper_line = Line::new("", rect_points)
                        .color(Color32::from_rgba_unmultiplied(0, 0, 0, 1)) // Почти прозрачный черный
                        .width(0.5); // Тонкая линия
                    plot_ui.line(helper_line);
                }
            }

            // Вызываем функцию обратного вызова для обработки пользовательских взаимодействий
            on_plot_ui(plot_ui, model);
        })
    }

    /// ОСНОВНАЯ ФФУНКЦИЯ РИСОВАНИЯ ОТДЕЛЬНЫХ ГРАФИКОВ
    /// Отображает отдельный график
    ///
    /// Рисует график функции с учетом его состояния (выбран/подсвечен)
    /// и параметров отображения из модели.
    ///
    /// # Параметры
    /// * `plot_ui` - изменяемая ссылка на интерфейс графика
    /// * `curve` - ссылка на данные графика для отображения
    /// * `model` - ссылка на модель данных графика
    pub fn render_curve(plot_ui: &mut PlotUi, curve: &PlotCurve) {
        if !curve.is_shown() {
            return;
        }

        let plot_points = PlotPoints::from_iter(curve.points.clone());
        let base_color = Color32::from_rgb(curve.color[0], curve.color[1], curve.color[2]);
        let color = if curve.selected {
            Color32::GOLD
        } else {
            base_color
        };

        let line = Line::new("", plot_points)
            .color(color)
            .name(curve.get_label())
            .width(2.0);
        plot_ui.line(line);

        if curve.highlighted {
            let highlighted_pts = PlotPoints::from_iter(curve.points.clone());
            let highlight_line = Line::new("", highlighted_pts)
                .color(Color32::RED)
                .width(3.0);
            plot_ui.line(highlight_line);
        }
    }

    /// Отображает информационную панель
    ///
    /// Показывает информацию о выбранном графике, активном выделении,
    /// текущем диапазоне отображения и инструкции по использованию.
    ///
    /// # Параметры
    /// * `ui` - изменяемая ссылка на пользовательский интерфейс egui
    /// * `model` - ссылка на модель данных графика
    pub fn render_info_panel(ui: &mut Ui, model: &PlotModel) {
        // Заголовок информационной панели
        ui.heading("TGA Plot Information");
        ui.separator();
        // Инструкции по использованию
        ui.label("📖 Instructions:");
        ui.label("• Click on a curve to select it");
        ui.label("• Drag to select area");
        ui.label("• Scroll to zoom, drag to pan");
        ui.label("• Double-click to reset view");
        ui.separator();
        // Информация о выбранном графике
        if let Some(index) = model.get_selected_curve_index() {
            // Получаем выбранный график
            let curve = &model.plots[index];
            // Отображаем имя выбранного графика
            ui.label(format!("📌 Selected: {}", curve.get_name()));
            ui.label(format!("📐 Plot: {} vs {}", curve.y_name, curve.x_name));
        } else {
            // Если график не выбран, показываем соответствующее сообщение
            ui.label("No curve selected");
        }

        if let Some(rect) = &model.interaction.selection_rect {
            let (x_min, x_max, y_min, y_max) = rect.bounds();
            ui.label(format!(
                "📐 Selection: x=[{:.2}, {:.2}], y=[{:.2}, {:.2}]",
                x_min, x_max, y_min, y_max
            ));

            let highlighted: Vec<String> = model
                .plots
                .iter()
                .filter(|c| c.highlighted)
                .map(|c| c.get_name())
                .collect();

            if !highlighted.is_empty() {
                ui.label(format!("🔴 Highlighted: {}", highlighted.join(", ")));
            }
        }

        ui.label(format!(
            "👁️ View: x=[{:.2}, {:.2}], y=[{:.2}, {:.2}]",
            model.interaction.view_range.0,
            model.interaction.view_range.1,
            model.interaction.view_y_range.0,
            model.interaction.view_y_range.1
        ));
        ui.label(format!("{}", model.message));
    }

    /// Форматирует метки сетки по оси X
    ///
    /// Создает равномерно распределенные метки сетки по оси X
    /// с "красивыми" значениями (1, 2, 5, 10 и т.д.)
    ///
    /// # Параметры
    /// * `input` - входные данные для форматирования сетки
    ///
    /// # Возвращает
    /// Вектор меток сетки
    fn aligned_x_grid_spacer(input: egui_plot::GridInput) -> Vec<egui_plot::GridMark> {
        // Получаем границы диапазона
        let (min, max) = input.bounds;
        // Проверяем корректность границ
        if !min.is_finite() || !max.is_finite() || max <= min {
            return Vec::new();
        }

        // Вычисляем шаг сетки
        let range = max - min;
        let step = Self::nice_step(range);
        // Проверяем корректность шага
        if !step.is_finite() || step <= 0.0 {
            return Vec::new();
        }

        // Создаем вектор меток
        let mut marks = Vec::new();
        // Вычисляем начальное значение
        let start = (min / step).floor() * step;
        let mut value = start;
        // Ограничитель для предотвращения бесконечного цикла
        let mut guard = 0;
        // Генерируем метки в пределах диапазона
        while value <= max && guard < 10_000 {
            // Добавляем метку, если она в пределах диапазона
            if value >= min {
                marks.push(egui_plot::GridMark {
                    value,           // Значение метки
                    step_size: step, // Размер шага
                });
            }
            // Переходим к следующей метке
            value += step;
            guard += 1;
        }

        // Если метки не созданы, добавляем одну в середине диапазона
        if marks.is_empty() {
            let mid = (min + max) * 0.5;
            marks.push(egui_plot::GridMark {
                value: mid,
                step_size: range,
            });
        }

        marks
    }

    /// Вычисляет "красивый" шаг для сетки
    ///
    /// Определяет оптимальный шаг сетки, который будет "красивым"
    /// (1, 2, 5, 10 и т.д.) в зависимости от диапазона.
    ///
    /// # Параметры
    /// * `range` - диапазон значений
    ///
    /// # Возвращает
    /// Оптимальный шаг сетки
    fn nice_step(range: f64) -> f64 {
        // Проверяем корректность диапазона
        if !range.is_finite() || range <= 0.0 {
            return 0.0;
        }
        // Вычисляем целевой шаг (примерно 6 делений на диапазон)
        let target = range / 6.0;
        if target <= 0.0 {
            return 0.0;
        }
        // Вычисляем степень 10 для целевого шага
        let pow10 = 10.0_f64.powf(target.abs().log10().floor());
        // Нормализуем шаг
        let mut step = target / pow10;
        // Округляем до ближайшего "красивого" значения
        step = if step <= 1.0 {
            1.0
        } else if step <= 2.0 {
            2.0
        } else if step <= 5.0 {
            5.0
        } else {
            10.0
        };
        // Возвращаем окончательный шаг
        step * pow10
    }
}
