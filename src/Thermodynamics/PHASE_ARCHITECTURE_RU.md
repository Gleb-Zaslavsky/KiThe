# Архитектура фазовой подсистемы

Этот документ описывает подсистему \`phase_*\` в термодинамическом слое. Он
охватывает границу между поиском термодинамических данных (\`SubsData\`) и
верхнеуровневым кодом равновесия или solver. Он не определяет алгоритм
равновесия и не выбирает неявно физическую модель активности.

## Основные понятия

| Термин | Значение в данном коде |
| --- | --- |
| **Transactional operation** | Многофазная операция работает с копиями payload и временными результатами, а публикует их только после успеха всех фаз и проверок. |
| **Cache** | Восстанавливаемое производное состояние: численные Gibbs/entropy, closures, символические выражения и контекст вычисления. Это не исходные данные библиотек. |
| **Snapshot** | Типизированный read-only снимок cache, согласованный с revision и содержащий layout и условия валидности. |
| **Cache bundle** | \`PhaseThermoCacheBundle\`: cache одной фазы, сгруппированный по семейству свойства и представлению. |
| **Layout** | \`SystemLayout\`, единственный канонический порядок solver-компонентов с учетом фазы. |
| **Resolution** | Ошибочный по природе переход от декларативного \`PhaseSpec\` к выбранным библиотечным записям, layout и provenance. |

## Назначение и границы

Подсистема находится между двумя уровнями:

\`\`\`text
Thermodynamic libraries / SubsData
    -> поиск, parsing, coefficients, transport and thermo records

phase_* subsystem
    -> typed phase declarations, deterministic component order,
       pure evaluation, transactional transitions, cache validity

Chemical-equilibrium and reactor consumers
    -> equations, equilibrium constraints, numerical solvers
\`\`\`

Ее задача - предоставить верхнему уровню согласованные данные о фазах: какие
компоненты существуют, в каком порядке они идут в solver, какие
термодинамические функции им соответствуют и при каких условиях получен
результат. Она не подменяет собой физическую модель или решатель равновесия.

## Структуры и зависимости

\`\`\`text
PhaseSpec --------------------------> ResolvedPhaseSystem
  |                                      |
  | PhaseId, components,                 | resolved SubsData payloads
  | physical state, model                | SystemLayout + provenance
  v                                      v
phase_domain.rs                     PhaseSystem
                                         |
              +--------------------------+--------------------------+
              |                          |                          |
              v                          v                          v
        phase_layout.rs            phase_evaluation.rs       phase_cache.rs
        canonical order            typed pure requests       bundles and snapshots
              |                          |                          |
              +--------------------------+--------------------------+
                                         |
                                         v
                    phase_operations.rs + phase_property_builders.rs
                    staged fallible work, no state publication
\`\`\`

\`PhaseSystem\` - единственный orchestration-модуль. Он владеет фазовыми
\`SubsData\`, cache bundles, revision, resolved layout и контекстом численных
вычислений. Вынесенные модули не получают права самостоятельно публиковать
состояние: они либо описывают данные, либо строят временный результат, либо
дают read-only представление.

Это защищает от частичной публикации. Если одна фаза уже получила новые
coefficients или closures, а обработка следующей фазы завершилась ошибкой,
внешний код не должен увидеть смешанное состояние.

## Жизненный цикл данных

\`\`\`text
1. Declare: PhaseSpec + lookup policy
        |
2. Resolve: ThermoRepository / SubsData lookup
        |
        +--> ResolvedPhaseSystem + provenance + SystemLayout
        |
3. Install: PhaseSystem installs payloads and invalidates derived state
        |
4. Evaluate/build: numeric request OR symbolic/closure staged construction
        |
5. Validate: phases, component order, finite values, conditions, revision
        |
6. Commit: one PhaseSystem publication and revision update
        |
7. Read: borrowed view or typed snapshot aligned to current layout
\`\`\`

\`revision\` и cache context - часть контракта, а не отладочная метка. После
изменения payload, layout или symbolic expression старый результат не должен
выглядеть актуальным.

## Почему компонент привязан к фазе

\`PhaseComponentId\` включает и фазу, и вещество. \`H2O(g)\` и \`H2O(l)\` могут
иметь одинаковый элементный состав, но являются разными solver-неизвестными:
для них могут отличаться стандартные состояния, химические потенциалы и модель
смешения. Хранение только по имени вещества привело бы к ошибочному слиянию
или перезаписи компонента.

\`SystemLayout\` задает устойчивый порядок таких компонентов. \`HashMap\` остается
формой для lookup и совместимости с library data, но не определяет порядок
solver-вектора. Это предотвращает недетерминированные векторы, перестановку
коэффициентов и ошибки матриц.

## Transactional mutation

\`\`\`text
visible PhaseSystem state
        |
        v
clone every phase payload
        |
        v
run a fallible operation for every phase
        |
   error? ---- yes ---> discard the whole stage
        |
        no
        v
validate results, request shape, finite values, and layout
        |
        v
one commit + revision update + cache invalidation/publication
\`\`\`

То же относится к numeric cache publication: сначала проверяются полный request,
покрытие фаз и компонентов, конечность значений и условия; лишь затем
заменяется cache family и context. Поэтому ошибка не может оставить новую
Gibbs-cache и старый context или наоборот.

## Cache, bundle и snapshot

\`\`\`text
PhaseThermoCacheBundle for one phase
|
+-- Gibbs property cache
|   +-- numeric: substance -> f64
|   +-- functions: substance -> closure
|   \`-- symbolic: substance -> Expr
|
\`-- Entropy property cache
    +-- numeric: substance -> f64
    +-- functions: substance -> closure
    \`-- symbolic: substance -> Expr
\`\`\`

Число, closure и \`Expr\` не считаются независимыми глобальными картами. Это
разные представления одного свойства одной фазы. Bundle удерживает данную
связь в структуре данных и уменьшает риск обновить \`dG_sym\`, забыв
инвалидировать согласованный численный \`dG\`.

Snapshot добавляет к cache условия валидности: revision, \`T\`, \`P\` и
нормализованный состав для numeric evaluation. Он также может проецировать
численные значения в плотный вектор строго в порядке \`SystemLayout\`, что важно
для solver hot paths.

## Почему pure evaluation отделен от publication

\`phase_evaluation.rs\` описывает typed request и валидирует его до запуска
калькулятора. \`phase_property_builders.rs\` строит symbolic expressions и
closures на копиях payload. Эти шаги могут вернуть \`SubsDataResult\`, но не
меняют \`PhaseSystem\`.

Это дает три свойства:

1. Ошибка одной фазы не портит другую.
2. Вычисление можно тестировать без проверки скрытого cache state.
3. Только одно место отвечает за revision, invalidation и publication policy.

Эта граница введена не ради абстракции. Она оставляет mutable \`SubsData\` там,
где он пока нужен для получения coefficient data, но не распространяет его
неявные мутации на верхние слои.

## Public API и legacy boundary

Внешний код работает через \`OnePhase\`, \`PhaseOrSolution\`, typed specs и
read-only views. \`PhaseSystem\` crate-private, потому что прямое изменение его
maps обходило бы invalidation и layout revision.

Legacy closures типа \`Fn(...) -> f64\` остаются совместимым представлением, но
не способны вернуть typed error во время вызова. Поэтому новые numeric paths
должны предпочитать fallible typed request/evaluation API. Полная миграция
старых equilibrium consumers отложена до работы над \`ChemEquilibrium\`.

## Ответственность модулей

| Модуль | Ответственность | Чем он не должен владеть |
| --- | --- | --- |
| \`phase_domain.rs\` | Объявления фаз, physical state/model, resolved records, provenance | cache mutation или evaluation policy |
| \`phase_layout.rs\` | Канонический упорядоченный phase/component layout | database lifecycle или cache publication |
| \`phase_evaluation.rs\` | Typed numeric request validation и pure evaluation | владение cache |
| \`phase_operations.rs\` | Transactional staging над копиями \`SubsData\` payload | commit или revision policy |
| \`phase_property_builders.rs\` | Staged symbolic и closure construction | публикация видимого state |
| \`phase_cache.rs\` | Типы cache bundle, views, snapshots, dense projections | mutable lifecycle ownership |
| \`phase_system.rs\` | Commit, invalidation, revisioning, cache lifecycle | второй public facade |
| facade modules | Ergonomic и compatibility-facing API | независимое engine state |

## Правило для будущих изменений

Для новой физической модели, свойства или solver adapter сначала нужно решить,
является ли он декларативными domain data, pure evaluation, staged
transformation или \`PhaseSystem\` lifecycle transition. Только последняя
категория может менять видимое phase state или публиковать cache. Это сохраняет
подсистему детерминированной, тестируемой и устойчивой к частичным обновлениям.

