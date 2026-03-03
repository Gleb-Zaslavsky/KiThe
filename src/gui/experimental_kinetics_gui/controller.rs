use crate::gui::experimental_kinetics_gui::model::PlotModel;
use egui_plot::PlotUi;

/// Контроллер для обработки взаимодействий пользователя с графиком
///
/// Этот структурный тип не содержит данных, а служит как контейнер
/// для методов обработки взаимодействий с графиком и пользовательским интерфейсом.
pub struct PlotController;

impl PlotController {
    /// Обрабатывает взаимодействия с областью графика (мышь, клавиатура)
    ///
    /// Эта функция отвечает за:
    /// - Обработку перетаскивания (панорамирование или выделение области)
    /// - Обработку кликов мыши (выбор графика)
    /// - Обработку прокрутки колеса мыши (масштабирование)
    /// - Обработку клавиш стрелок (панорамирование)
    /// - Обработку двойного клика (сброс вида)
    ///
    /// # Параметры
    /// * `plot_ui` - ссылка на пользовательский интерфейс графика
    /// * `model` - изменяемая ссылка на модель данных графика
    pub fn handle_plot_interactions(plot_ui: &mut PlotUi, model: &mut PlotModel) {
        // Получаем ответ от пользовательского интерфейса графика (структура Response отвечает за обработку взаимодействий)
        let response = plot_ui.response();

        // Получаем позицию указателя мыши в координатах графика
        if let Some(pointer_pos) = plot_ui.pointer_coordinate() {
            // Преобразуем координаты в массив f64 для дальнейшей работы
            let mouse_pos = [pointer_pos.x as f64, pointer_pos.y as f64];

            // Обрабатываем начало перетаскивания (определяем тип: панорамирование или выделение)
            if response.drag_started() {
                // Получаем состояние ввода для определения нажатой кнопки мыши
                let input = response.ctx.input(|i| i.clone());
                // Проверяем, нажата ли левая кнопка мыши
                let primary_down = input.pointer.button_down(egui::PointerButton::Primary);
                if primary_down {
                    // Левая кнопка мыши - начинаем панорамирование
                    model.interaction.start_pan(mouse_pos);
                } else {
                    // Другая кнопка мыши (например, правая) - начинаем выделение области
                    model.start_selection(mouse_pos);
                }
            }

            // Обрабатываем процесс перетаскивания
            if response.dragged() {
                // Если в данный момент происходит панорамирование
                if model.interaction.is_panning {
                    // Обновляем панорамирование с новой позицией мыши
                    model.update_pan(mouse_pos);
                    // Запрашиваем перерисовку для отображения изменений
                    response.ctx.request_repaint();
                // Если в данный момент происходит выделение области
                } else if model.interaction.is_selecting {
                    // Обновляем выделение с новой позицией мыши
                    model.update_selection(mouse_pos);
                    // Запрашиваем перерисовку для отображения изменений
                    response.ctx.request_repaint();
                }
            }

            // Обрабатываем завершение перетаскивания
            if response.drag_stopped() {
                // Если завершается панорамирование
                if model.interaction.is_panning {
                    model.interaction.end_pan();
                // Если завершается выделение области
                } else if model.interaction.is_selecting {
                    model.end_selection();
                }
                // Запрашиваем перерисовку для отображения изменений
                response.ctx.request_repaint();
            }

            // Обрабатываем клик мыши (выбор графика)
            if response.clicked()
                && !model.interaction.is_selecting
                && !model.interaction.is_panning
            {
                // Ищем ближайший график к точке клика
                if let Some(curve_index) = model.find_nearest_curve(mouse_pos) {
                    // Выбираем найденный график
                    model.select_curve(curve_index);
                    // Очищаем выделение области
                    model.clear_selection_rect();
                } else if model.is_only_one_shown() {
                    if let Some(curve_index) = model.plots.iter().position(|curve| curve.is_shown())
                    {
                        model.select_curve(curve_index);
                        model.clear_selection_rect();
                    } else {
                        model.clear_selection();
                    }
                } else {
                    // Если график не найден, очищаем выбор
                    model.clear_selection();
                }
            }
        }

        // Обрабатываем масштабирование (прокрутка колеса мыши)
        let scroll_delta = response.ctx.input(|i| i.raw_scroll_delta);
        // Проверяем, было ли движение колеса мыши по вертикали
        if scroll_delta.y != 0.0 {
            // Определяем направление прокрутки (вверх или вниз)
            let dir = if scroll_delta.y > 0.0 { 1.0 } else { -1.0 };
            // Получаем величину прокрутки
            let magnitude = scroll_delta.y.abs();
            // Базовый коэффициент масштабирования (определяет скорость масштабирования по мере прокрутки)
            let base = 1.02_f64;
            // Вычисляем коэффициент масштабирования с учетом направления и величины
            let zoom_factor = base.powf(dir * (magnitude as f64 * 0.5));
            // Получаем позицию указателя мыши для масштабирования относительно точки
            if let Some(pointer_pos) = plot_ui.pointer_coordinate() {
                // Применяем масштабирование относительно позиции указателя
                model.zoom(zoom_factor, [pointer_pos.x as f64, pointer_pos.y as f64]);
                // Запрашиваем перерисовку для отображения изменений
                response.ctx.request_repaint();
            }
        }

        // Обрабатываем панорамирование с помощью клавиш стрелок
        let input = response.ctx.input(|i| i.clone());
        // Если нажата стрелка влево, сдвигаем вид влево
        if input.key_down(egui::Key::ArrowLeft) {
            model.interaction.pan(-0.1);
        }
        // Если нажата стрелка вправо, сдвигаем вид вправо
        if input.key_down(egui::Key::ArrowRight) {
            model.interaction.pan(0.1);
        }

        // Сброс вида при двойном клике
        if response.double_clicked() {
            model.reset_view();
            response.ctx.request_repaint();
        }
    }
}
