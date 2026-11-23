# GUI: тестирование workflow, UX и клавиатурное управление

Цель: GUI является основным способом работы с движком экспериментальной кинетики, поэтому
проверяем не только наличие виджетов, но и завершённые пользовательские маршруты: от загрузки
данных до преобразования, анализа и сохранения результата. Тесты должны запускать настоящее
egui-приложение через `egui_kittest` и проверять как UI-взаимодействия, так и состояние модели.

## 1. Сценарный тестовый каркас

- [ ] Развить существующие `with_plot_app_state_harness` и тестовые генераторы в лёгкий GUI
  scenario harness: общие операции запуска кадра, открытия меню/окна, ввода текста и проверки
  сообщения/состояния модели.
- [ ] Сохранить разделение по уровням: smoke-тесты рендеринга; тесты одного действия; сквозные
  сценарии нескольких шагов. Для каждого workflow использовать устойчивые accessibility labels,
  не привязываясь к внутреннему порядку виджетов.
- [x] Проверить доступность выпадающих списков: `egui_kittest` открывает ComboBox по подписи и
  выбирает пункт popup по доступному имени (`Model:` → `Linear`), без координат.
- [ ] Добавить компактные helpers для создания приложения с пустой серией, одним экспериментом,
  несколькими экспериментами и несколькими кривыми одного опыта. Использовать реальные fixtures
  и детерминированные синтетические данные там, где важны численные результаты.
- [ ] Для операций, меняющих данные, проверять не только клик: результат в серии/dataset,
  binding и единицы, список колонок, кривые, выбранную кривую, системные сообщения, History/undo.
- [ ] Добавить сценарии ошибок и отмены: некорректный ввод или неудачная операция не должны
  частично менять dataset, историю, графики либо выбор пользователя.

## 2. Основные пользовательские маршруты

- [x] **Импорт и создание опыта:** открыть импорт, выбрать файл, задать/проверить метаданные,
  связать time/temperature/mass, создать кривую и увидеть точки на графике.
- [x] Проверить автоматический binding по распространённым заголовкам (`t`, `T`, `m`) и создание
  кривых после импорта; пустой ID по умолчанию заполняется именем файла.
- [x] Проверить ручное исправление binding на пользовательских заголовках (`elapsed`, `probe`,
  `signal`) и успешный импорт через GUI.
- [x] **Серия опытов:** импортировать два опыта через GUI, проверить разные ID и соответствие
  каждой кривой своему experiment/X/Y; повторный ID отклоняется без изменения числа опытов и графиков.
- [x] **Предобработка массы и температуры:** для `T vs t` и `m vs t` GUI-сценарии проверяют
  создание `sg_T`/`sg_m`, переключение активного источника и появление новой кривой; обработка `T`
  также подтверждает неизменность времени и создание `dT/dt`, сценарий массы подтверждает, что
  температурная кривая не меняется. Проверку `dm/dt` выполнять только на данных с физической
  единицей массы: имеющийся fixture `m` помечен как mV и правильно отклоняется rate-методом.
- [x] **Кинетическая обработка:** GUI-сценарий через Quick Actions выполняет калибровку массы
  mV → mg → conversion → conversion rate; проверяет появление `conversion`/`deta_dt`, сохранность
  выбранной кривой и единицы `Dimensionless`/`PerSecond`. Предобработка source → smoothing/filter
  отдельно проверена для массы и температуры выше.
- [x] **Кинетический метод:** из трёх синтетических опытов подготовить conversion/rate, открыть
  Isoconversional через меню и выполнить OFW кнопкой Calculate; проверить статус, опыт результата,
  кривую `eta vs Ea` и отказ на пустой серии без частичных данных.
- [x] **Direct Problem / IVP:** через меню открыть окно, решить детерминированную модель F1, проверить
  создание `ivp_result`, колонки `time`/`conversion` и сразу видимую кривую без Manage Plot.
- [x] **Фиттинг:** кнопкой Fit выполнить линейный fit выбранных `t`/`T`; проверить метрики,
  fitted-колонку, provenance и видимую fitted-кривую; модель выбирается в UI через ComboBox.
  [x] Таблица формул сделана сворачиваемой и закрыта по умолчанию, чтобы основная кнопка Fit
  оставалась в видимой части окна. GUI end-to-end проверяет выбор модели через ComboBox.
  [x] Начальные X/Y привязаны к осям выбранной кривой; если она не подходит, fallback детерминирован
  по имени колонок, а не по порядку `HashMap`.
- [x] **Сохранение и повторная загрузка:** пройти Save As CSV через GUI-диалог, загрузить результат
  обратно и сверить число опытов и колонки; binding и преобразованные единицы также проверены.
- [x] **Явный контракт экспорта:** Typed CSV сохраняет всю серию либо один выбранный опыт вместе
  с ExperimentMeta, bindings, единицами, nature/origin и provenance. Plain TXT является обычной
  табличной выгрузкой одного выбранного опыта без метаданных и не обещает round-trip.
- [x] В диалоге экспорта показывать область сохранения и список колонок с выбором всех колонок
  по умолчанию. Для серии колонки выбираются отдельно внутри каждого опыта; нет предположения,
  что одинаковые заголовки в разных опытах означают одинаковые данные. Сохранённые bindings
  очищаются только если связанная колонка исключена. Не разрешать сохранить пустой набор колонок.
- [x] Развести реализации CSV и TXT: `Save As CSV` использует typed-series serializer, `Save As TXT`
  пишет настоящий tab-delimited plain text одного опыта. Сохранение одного опыта в Typed CSV должно
  оставаться полноценным форматом серии из одного опыта, чтобы не терять ExperimentMeta.
- [x] Сохранить в shortcut Save последнюю фактическую конфигурацию экспорта (format, scope, columns,
  path), чтобы повторное нажатие Ctrl/Cmd+S не меняло формат и область незаметно.
- [x] Покрыть доменные serializers отборами опытов/колонок, метаданные и bindings после round-trip,
  длины опытов и заголовки plain TXT. Добавить egui_kittest сценарии выбора scope/columns, typed CSV
  для серии и одного опыта, plain TXT, пустого выбора и повторного Save.
- [ ] **Plot и Column Manager:** проверить несколько живых графиков, выбор/скрытие/удаление кривой,
  Column Manager, samples, provenance и quality details; проверить, что преобразования и обрезка
  обновляют все затронутые графики.
  [x] Проверены синхронизация Undo после in-place преобразования и после удаления выходной
  колонки, сохранение кривой при переименовании оси и обновление индексов после удаления опыта.
  [x] GUI-сценарий Column Manager открывает окно, видит числовые samples, выбирает колонку по
  accessible role+label и проверяет provenance/history. Для сглаженной колонки также проверяется
  таблица quality-метрик: RMSE, normalized RMSE, correlation и оба roughness ratio с пояснениями.
- [x] **Undo:** проверить восстановление provenance/истории и сообщение при пустом undo stack.
  Добавлены кнопка Undo на панели действий и `Ctrl/Cmd+Z`; GUI-сценарий проверяет восстановление данных и provenance,
  удаление выходной колонки, обновление графика, запись undo в History и отказ без мутации при пустом
  стеке снимков.

## 3. Покрытие основных операций

- [ ] Сгруппировать действия меню Math и Quick Actions и проверить каждое доступное пользователю
  действие хотя бы на одном успешном и одном граничном/ошибочном сценарии.
  [x] Успешный UI-путь проверен для Hampel, rolling mean, spline, LSQ spline и LOWESS;
  каждый сценарий проверяет новую колонку и появление соответствующей кривой. SG проверяется отдельно.
- [x] Проверить cut selected, cut before time, move time to zero и преобразования единиц на
  нескольких кривых одного опыта; графики перестраиваются согласованно.
- [x] Проверить идемпотентные команды единиц: повторное нажатие не меняет данные, а сообщает,
  что преобразование уже выполнено.
- [x] Проверить арифметику на выбранной кривой при наличии нескольких X/Y-кривых одного опыта:
  затрагиваются именно выбранные колонки, не первые попавшиеся и не общий default plot.
- [ ] Для фильтров и fitting проверить параметры на границе допустимости и ошибочные значения.
  [x] SG без выбранной кривой показывает ошибку и не меняет серию/графики.
  [ ] Добавить граничные/ошибочные сценарии для остальных фильтров и fitting.
  [ ] Уточнить и тестировать выходные имена: одинаковое имя source/output означает in-place,
  а занятое имя другой колонки требует отдельного контракта.
- [x] Проверить selection/zoom/reset/delete plot и сохранение выбора после перерисовки кривых.
  [ ] Доработать `Clear Selected`: в текущем Harness-кейсе виджет находится, но click не доходит
  до обработчика; существующий тест отдельно очищает модель и не является end-to-end проверкой кнопки.

## 4. UX-аудит и устранение ложных affordances

- [ ] Провести полную инвентаризацию всех видимых пунктов меню и кнопок: каждое действие должно
  работать, быть явно недоступным с причиной либо быть удалено.
  [x] Удалены `Run Tests` (заглушка) и дублирующий мёртвый пункт Golden Pipeline; Criado явно
  помечен как недоступный до научной валидации; Kissinger больше не содержит кнопок с `todo!()`;
  `Save As Image` переименован в `Plot Editor`, так как он открывает редактор. egui_kittest
  закрепляет отсутствие/недоступность этих ложных команд.
- [x] Удалить из File Manager заглушки `Append data to experiment`, `Export Data`, `Manage Data` и
  `Save As Excel`; egui_kittest проверяет, что они не отображаются, а рабочие Import/Manage Plot/
  Save CSV/TXT остаются доступны.
- [x] Согласовать единый маршрут импорта: `Import Data` и `New experiment` не должны выглядеть как
  два разных действия, если фактически запускают один процесс.
- [ ] Для ошибок показывать короткое пользовательское объяснение рядом с действием или в системных
  сообщениях; сохранять технический контекст в логах. Проверить ошибки чтения, binding, фильтра,
  kinetic method, fit и экспорта.
  [x] Проверен отказ SG без выбранной кривой: сообщение отображается, состояние не мутируется.
- [ ] Проверить, что после каждой мутации UI отражает актуальные данные сразу, без закрытия и
  повторного открытия графика/менеджера.
- [x] Проверить фокус клавиатуры и текстовые поля: ввод в поле импорта не запускает Save/Undo;
  Esc в TextEdit не закрывает окно.
- [x] Показывать доступные shortcuts в подписях меню/tooltip и не назначать конфликтующие
  сочетания на системные или специфичные для Plot взаимодействия (`O`, `S`, `Shift+S`, `Z`, `Esc`).

## 5. Стандартные сочетания клавиш

- [x] Реализовать `Ctrl/Cmd+O` для открытия импорта данных (не обозначать его как открытие проекта,
  пока отдельный формат проекта не существует).
- [x] Реализовать `Ctrl/Cmd+S` для сохранения серии в CSV: при первом сохранении открывается диалог,
  после успешного сохранения повторная команда записывает в тот же путь. `Ctrl/Cmd+Shift+S` открывает
  Save As.
- [x] Реализовать `Ctrl/Cmd+Z` для undo выбранного опыта; определить и реализовать redo отдельно
  только при наличии соответствующей доменной поддержки, не имитировать redo повторным запуском.
- [x] `Esc` закрывает только окно на верхнем egui-слое; при фокусе ввода клавиша остаётся штатному
  поведению egui и не закрывает окно. Поведение проверено для двух открытых окон и TextEdit.
- [ ] Определить, нужен ли отдельный shortcut для очистки выбора; не перехватывать `Delete`, пока
  не определено единое ожидаемое поведение для графика и таблиц.
- [ ] Расширить проверки сочетаний в `egui_kittest` через реальные key events: команда срабатывает один раз,
  обрабатывает Ctrl и Command в соответствии с платформой, не срабатывает в неподходящем окне или при
  открытом модальном диалоге.
  [x] Реальные key events проверяют импорт, Save/Save As, Undo и то, что команды не перехватываются
  текстовым полем; `Esc` закрывает только верхнее окно и не закрывает окно при редактировании.
  Shortcut не запускается при открытом menu popup. В приложении сейчас нет настоящих modal dialogs;
  отдельная политика модальных окон понадобится, если они появятся.
- [x] Убедиться, что shortcut Save и пункт меню Save вызывают один обработчик;
  не дублировать бизнес-логику в обработчике клавиатуры.

## 6. Критерии готовности

- [ ] Для каждого основного workflow выше есть хотя бы один `egui_kittest` сценарий, который проходит
  через пользовательские виджеты и проверяет конечное состояние модели.
- [ ] Все доступные действия меню имеют тест успешного пути; критические мутации имеют проверку
  отказа без частичного изменения состояния.
- [ ] Shortcut-тесты покрывают действие, фокус текстового поля, модальные окна и platform modifiers.
  [x] Импорт/Save/Save As/Undo проверены реальными key events и не перехватываются при фокусе в поле;
  Esc закрывает только верхнее окно. Отдельный тест модального диалога и явное сравнение Ctrl/Command
  на разных платформах ещё не добавлены.
- [x] Запустить GUI-тесты экспериментальной кинетики последовательно:
  `cargo test gui::experimental_kinetics_gui:: -- --test-threads=1` — 73 теста прошли для lib target
  и 73 для bin test target.

## 7. Column Manager: достоверность данных и UX

### Согласованный UX-рефакторинг окна

- [x] В обычном режиме выбирать один активный опыт; таблица содержит только его колонки.
- [x] Клик по заголовку выбирает колонку для инспектора; checkbox используется только для массового
  удаления и отключён в режиме сравнения.
- [x] Переименование выполнять для строки инспектора; поле имени и Rename находятся рядом,
  checkbox для выбора цели переименования не используется.
- [x] Массовое удаление колонок отключено без выбора, показывает список имён и требует отдельного
  подтверждения до первой мутации.
- [x] Удаление опыта выполняется отдельной командой для явно показанного активного ID, без связи
  с выделением колонок, и имеет отдельное подтверждение.
- [x] Comparison mode read-only: отключены Rename/Delete/Metadata/Sort/Time check.
- [x] Сложить редкие действия (Metadata, проверка и сортировка времени) в More; основные команды
  просмотра и удаления оставить на виду.
- [x] Разделить инспектор на Overview / Processing / Quality / History. Overview показывает имя,
  единицу, ColumnNature, schema binding, активный источник расчёта и состояние графика.
- [x] egui_kittest проходит по вкладкам Overview → Processing → Quality → History и проверяет
  отображение метаданных, quality-метрик и истории операций.
- [x] Сделать Samples вертикально регулируемой областью с собственным scroll; Inspector остаётся
  отдельной прокручиваемой областью и имеет компактные вкладки.
- [x] Сообщения успеха/предупреждения/ошибки/информации используют явный тип и разные цвета;
  уровень включён в доступный текст сообщения, цвет не является единственным индикатором.
  Unit-тест проверяет типизированное сообщение и доступное представление.
- [x] Основные команды, вторичные действия, подтверждения, выборы колонок и поля metadata имеют
  tooltip с нужными параметрами и последствиями действия; egui_kittest проверяет tooltip по hover.
- [x] Проверить полный набор действий Column Manager и пояснить цель checkbox отдельно от клика
  по названию колонки.

- [ ] Использовать Column Manager как инспектор данных и веток обработки, а не как spreadsheet:
  не добавлять редактирование произвольных ячеек, формулы, pivot tables или каскадное изменение
  исходных данных.
- [x] Перестроить предпросмотр так, чтобы значения разных колонок одного опыта относились
  к одной строке исходного DataFrame. Выбирать общие индексы строк, показывать исходный индекс и
  null показывать явно, а не пропускать независимо для каждой колонки. Для выбора строки теперь
  используется один frame collect; тест проверяет разные позиции null и сохранение alignment.
- [ ] Не смешивать строки разных экспериментов в одной таблице по умолчанию. Сначала выбирать один
  эксперимент; сравнение нескольких экспериментов сделать отдельным явным режимом с указанием,
  что временные сетки могут различаться.
- [x] В режиме сравнения дать выбрать несколько опытов и полный набор колонок каждого опыта,
  а не ограничивать сценарий одной парой X/Y. Колонки разных опытов могут иметь разные имена и
  единицы; неявно сопоставлять их по позиции нельзя. Общие имена/роли можно предлагать как подсказку,
  но финальный выбор остаётся за пользователем.
- [x] Сравнительный preview семплирует одинаковое число строк от каждого выбранного опыта, если
  данных достаточно; показывать для каждого опыта свои исходные индексы и значения всех выбранных
  колонок. Не утверждать, что строки разных опытов соответствуют одному времени. Интерполяция на
  общую сетку, если понадобится, должна быть отдельной явно выбранной операцией.
- [ ] Развести просмотр колонки и выбор цели операции: Details и группа checkbox должны явно
  показывать, какая колонка будет переименована/удалена. Проверить групповые действия для нескольких
  экспериментов и обеспечить детерминированный порядок либо единый предварительно проверенный batch.
- [x] Удаление эксперимента требует выбора колонок ровно из одного опыта; неоднозначный выбор
  отклоняется вместо выбора первого ID из HashMap. Атомарность групповых удалений и undo остаются.
- [x] Показывать числа с адаптивным форматом, чтобы малые производные не отображались как `0.0000`.
- [ ] Добавить копирование полных значений с заголовками и единицами, а не форматированного preview.
- [ ] Различать успешные сообщения, предупреждения и ошибки визуально и по доступности; не красить
  все сообщения как ошибки.
- [ ] В инспекторе отдельно показывать классификацию ColumnNature, активную привязку схемы,
  выбранный GUI-источник расчёта и наличие колонки на графике. Выбор для просмотра не должен менять
  binding или источник вычислений.
- [ ] Организовать детали выбранной колонки в компактные разделы Overview / Processing / Quality /
  History; сделать samples прокручиваемой и регулируемой областью, чтобы Quality/History оставались
  доступны в небольшом окне.
- [ ] Добавить диагностику качества данных: число строк, null/NaN/inf, диапазон, а для времени также
  дубликаты и нарушения порядка. Тяжёлые проверки запускать по запросу и кэшировать до изменения данных.
- [ ] Добавить сравнение нескольких обработанных колонок: общий интервал, явная QualityReference,
  метрики изменения сигнала отдельно от точности относительно эталона; не ранжировать фильтры
  автоматически.
- [ ] Добавить команды показа выбранной колонки на графике и показа выбранного интервала, используя
  существующий plot workflow и явно показывая выбранные X/Y.
- [ ] Для реконструированных данных показывать preview общей X-сетки, диапазона, числа точек,
  моделей, параметров и R². Отмена evaluation не должна удалять сохранённые fit-модели.
- [ ] Перед удалением предупреждать о роли колонки в schema, графиках и активных источниках;
  не удалять потомков provenance каскадно.
- [ ] Рассмотреть поиск, скрытие/закрепление колонок, копирование с заголовками и единицами и экспорт
  выбранных колонок. Явно различать экспорт полного набора и видимого sample.
- [ ] Не выполнять тяжёлую материализацию/диагностику каждый кадр egui; получать согласованный sample
  одним запросом на эксперимент и обновлять только после изменений данных.
- [ ] Добавить egui_kittest-регрессии на соответствие sample исходным строкам при разных null,
  формат малых значений, выбор цели операции, неоднозначное удаление эксперимента и доступность
  Quality/History при ограниченной высоте окна.
- [x] Назначить `F1` открытие Column Manager без модификаторов; проверено, что обработчик уважает
  существующий общий gate клавиатуры: не срабатывает при текстовом фокусе или открытом popup.
  egui_kittest проверяет открытие окна по F1.
- [x] Назначить `F2` / `F3` / `F4` для открытия меню File Manager / Math / Kinetic Methods;
  меню открываются реальными key events в egui_kittest и остаются заблокированы при открытом popup
  или фокусе ввода. Клавиши видны в tooltip-ах меню и Column Manager.
- [ ] Проверить поведение `F1` при повторном нажатии и после появления в приложении модальных окон;
  сейчас модальных окон в этом workflow нет, а повторное нажатие идемпотентно открывает то же окно.
# Production-Readiness Story-Test Gates

These are end-to-end user-story gates, not isolated widget tests. Drive the real
egui application with `egui_kittest` wherever native file dialogs permit it;
use deterministic temporary CSV fixtures and the existing synthetic-data
generators. Assert both visible UI feedback and resulting dataset/series state.
A gate must report a real failure when the workflow is broken; do not weaken its
assertions to make an incomplete workflow appear successful.

- [x] **Continue one experiment from a second file.** Import both CSV segments
  through the GUI, first checking that an out-of-order continuation is rejected
  without changing the source, then appending a reset-time continuation. Verify
  row alignment across `t`, `T`, and `m`, the boundary and total row count,
  bindings/units, refreshed plots, and preservation of the source experiment.
- [x] **Reconstruct a missing-data interval.** Generate noisy synthetic `t`,
  `T`, and `m` data with a contiguous interval absent from acquisition (no
  placeholder zeros or null rows). Through the GUI, fit the available data on
  both sides as one fit set, inspect fit quality, evaluate temperature and mass
  on one shared X-grid, and commit the reconstruction as a new experiment
  without mutating the source. Verify gap boundaries, row alignment, units,
  fit report, and that a deliberately unsuitable fit is rejected.
- [x] **Three-experiment isothermal kinetics workflow.** Generate three
  deterministic noisy isothermal experiments with `t`, `T`, and `m` headers.
  Import and process them through the GUI: bind seconds and millivolts,
  calibrate mass, convert time to hours, derive conversion, verify smoothing
  reduces roughness, differentiate, and clean derivative edge rows. Export
  only final columns (`t`, `T`, `sg_mass`, `eta`, `deta_dt`) through the normal
  CSV action, re-import them, and run Friedman Integral. Also fit the three
  isothermal rates and recover the fixture's activation energy and
  pre-exponential factor within the test tolerances. Negative metadata coverage
  verifies that unsuitable inputs are rejected.
- [x] **Non-isothermal kinetics workflow.** Process a deterministic noisy
  three-heating-rate series through GUI import, binding/calibration, smoothing,
  conversion, and rate preparation. Export only final columns, re-import, and
  run Differential Friedman. Compare the recovered activation energy against
  the synthetic truth within the documented test tolerance; a missing-heating-
  rate metadata case verifies an informative failure.
- [x] **Quantitative generator/filter/method matrix.** Split the synthetic
  story coverage into deterministic noise generation, GUI-handler filter
  comparisons (SG, Hampel/MAD, rolling mean, LOWESS, linear resampling, and
  LSQ spline over Gaussian/uniform/drift/impulse/flicker profiles), and GUI
  kinetic-method runs. Print seeds, objective filter metrics, truth RMSE,
  method R2 where available, Ea/A estimates, and parameter errors as `tabled`
  console output. Treat this as an exploratory scientific report: do not
  require every filter to improve every noise profile, and do not silently
  accept or normalize systematic method bias. Results and interpretation are
  recorded in [`story_tests.md`](story_tests.md).
- [x] Run the full Golden Pipeline from its GUI dialog across isothermal and
  non-isothermal fixtures, five deterministic noise profiles, and SG windows
  11/31/51; print all five `FilterQualityReport` metrics and paired-truth RMSE
  before/after resampling. The test enlarges the egui_kittest viewport, clicks
  the actual front-panel Golden Pipeline button and dialog Calculate widget,
  then verifies the resulting experiment.
- [x] Show the SG mass-filter quality report in the Golden Pipeline dialog on
  success: both roughness ratios, RMSE, normalized RMSE, correlation, compared
  columns, and a note that these metrics are not truth accuracy. Keep the
  dialog open so the user can inspect the report.
- [ ] Add a distinct quantitative diagnostic for failures and a visible
  explanation that LSQ resampling changes the grid, so its quality cannot be
  compared directly by the current same-grid report.

For every gate, record fixture seed/configuration, selected GUI actions, expected
scientific parameters, numerical tolerances, final exported/imported schema, and
failure diagnostics. Run the successful baseline deterministically; keep
negative cases as explicit assertions about safe failure and unchanged state.
Numerical filter/fit reports and the observed fixture results are documented in
[`story_tests.md`](story_tests.md); simulated workflow tests print `[STORY]`
lines with quality metrics and recovered parameters.

## Isothermal experiments without a measured temperature column

- [x] Make thermal mode explicit in New Experiment / Import Data. Allow an
  isothermal import with a manually entered temperature and no T binding; keep
  measured-temperature binding optional for that mode. Preserve current
  non-isothermal import and validation behavior.
- [x] Support Append to Experiment for isothermal targets without a T column;
  require an incoming temperature column only when the target's data schema
  actually contains one. Preserve mode and temperature metadata on commit.
- [x] Replace the unconditional `Temperature source` presentation with a
  measured-temperature source control shown only when a usable column exists.
  When metadata defines an isothermal temperature, show it as the calculation
  temperature instead of an empty source selector.
- [x] When a measured temperature column exists, provide an explicit action to
  calculate its mean over the experiment and adopt that value as the isothermal
  temperature. Convert the selected unit to Kelvin, retain the original column,
  report the valid sample count and result, and do not infer isothermal mode
  merely because the column appears flat.
- [x] Make Golden Pipeline's GUI preflight and result explanation mode-aware.
  Isothermal data without a T column must complete mass/conversion processing
  using metadata; do not fabricate measured T or dT/dt columns. Keep measured T
  available when present, and leave the non-isothermal path unchanged.
- [x] Label kinetic methods by thermal compatibility and prevent running
  non-isothermal methods on isothermal experiments. Ensure isothermal methods
  accept metadata temperature without requiring a bound T column. All failures
  must reach the main system-message label.
- [x] Add egui_kittest routes for manual-temperature import of t/m-only data,
  explicit sensor-column averaging, Golden Pipeline with and without measured
  temperature, Append without a T column, compatible isothermal method
  execution, and clear rejection of an incompatible non-isothermal method.
  Assert visible status and resulting dataset/schema/metadata, not only clicks.
- [x] Final validation for this refactor: run the complete
  `Kinetics::experimental_kinetics` domain test module and the complete
  `gui::experimental_kinetics_gui` test module, including story tests.

## GUI Production Hardening Follow-up

- [x] P1: Make import thermal-mode controls mutually consistent. Clearing or
  hiding the sensor-mean option must clear its state; validate the active mode,
  not stale hidden fields. Test mode switching both directions and verify
  resulting metadata and required bindings.
- [x] P1: Make Column Manager metadata edits atomic. Parse and validate the full
  draft first (including finite positive heating rate/setpoint and incompatible
  simultaneous values), then commit once; on invalid input preserve every field.
  Add egui_kittest coverage for partial-update regressions and valid updates.
- [x] P1: Separate Golden Pipeline data commit from plot refresh. Preserve all
  unrelated experiment curves; refresh only the processed experiment for
  in-place mode; for save-as-new add only its curve. If plotting fails, keep
  processed data and History, remove partial/stale curves for that experiment,
  and publish a visible post-commit warning. Test with two live experiments and
  injected plot failure.
- [ ] P2: Make source selectors communicate role-compatible columns and units;
  retain access to unusual sources through an explicit advanced choice. Before
  transformations, show source -> output, units, and whether the column is
  created or overwritten. Test selection and operation-target consistency.
- [ ] P2: Expand GUI guarantees beyond dialog-local `last_error`: assert the
  main status label describes failure or partial success and the available
  recovery action. Add selected-column type/unit and large-series Column
  Manager latency benchmarks before optimizing its sample cache.
- [x] Final verification: focused regressions, complete
  `gui::experimental_kinetics_gui` tests, and `git diff --check`.

## Reaction-Network Direct-Problem Series

- [x] Add `Complex Direct Problem` alongside `Simple One-Step Kinetics`.
- [x] Support isothermal temperature lists and linear heating-rate lists in
  K/min, with shared time bounds and independent runs from the same initial
  concentrations.
- [x] Convert accepted runs to separate experiments with explicit simulated
  provenance, temperature/concentration units, and run-condition metadata.
- [x] Preflight experiment IDs and publish the batch only after every solve and
  conversion succeeds; create visible concentration-vs-time curves afterward.
- [x] Cover isothermal GUI publication and heating-rate parsing with tests.
- [ ] Deferred: add native mid-solve cancellation when `StiffAPI` supports
  combining parameter rebinding and `CancellationToken`.
- [ ] Deferred: support condition-specific time intervals if a concrete
  workflow requires equal final temperatures at different heating rates.
