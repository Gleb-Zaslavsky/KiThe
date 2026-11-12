# Запрос в RST: произвольная начальная сетка для Damped и CorrectionNewton

> Исторический запрос, подготовленный до RST 0.4.35. Explicit mesh и
> node-aligned guess теперь поддерживаются. Старые `UnsupportedRoute` ниже
> описывают прежнюю версию; текущие численные проверки находятся в
> `rst_solver_selection.rs`, а production-facade проверки — в
> `production_mesh_matrix.rs` и Story MF-05. Не использовать старую команду
> rejection-теста как текущий контракт.

## Цель

Нужно дать `BVP_Damp::Damped` и `BVP_CorrectionNewton` возможность принимать
явную начальную строго возрастающую сетку и seed, заданный в её узлах. Это
нужно не только для удобства ввода: на одной и той же задаче кластеризованная
сетка проходит через BVP_sci/SciPy, тогда как равномерная сетка исчерпывает
лимит узлов. Сейчас оба указанных маршрута отвергают кластеризованную сетку
в `BVP_api` до запуска решателя.

Запрос относится только к численному BVP API. Ниже приведена автономная
система `u'(z)=f(z,u)`; физическая интерпретация для реализации не требуется.

## Standalone system

Состояние `u=[u0,u1,u2,u3,u4,u5]`, `z in [0,1]`. Параметры:

```text
delta = 1e-4
beta = 0.05
gamma = 0.01
lambda = 0.1
D = 1e-4
Pe = 12017.324603855763
```

Определить

```text
S(u) = delta * exp(u0 / (1 + beta*u0)) * u2 / (1 + beta*u0)
```

и решить:

```text
u0' = u1 / lambda
u1' = Pe*u1 - lambda*Pe^2*S(u)
u2' = u3 / D
u3' = Pe*u3 + D*gamma*Pe^2*S(u)
u4' = u5 / D
u5' = Pe*u5 - D*gamma*Pe^2*S(u)
```

Boundary conditions:

```text
u0(0) = 0
u1(1) = 0
u2(0) = 1
u3(1) = 0
u4(0) = 0
u5(1) = 0
```

## Shared inputs and control results

Воспроизводитель и данные находятся в KiThe:

```text
src/SimpleReactorsBVP/literature_mf1963/rst_solver_selection.rs
src/SimpleReactorsBVP/literature_mf1963/six_state_parity.rs
src/SimpleReactorsBVP/literature_mf1963/six_state_parity/prepare_seed.py
src/SimpleReactorsBVP/literature_mf1963/six_state_parity/solve_scipy.py
src/SimpleReactorsBVP/literature_mf1963/fixtures/mf1963_mf_low_six_state_seed.csv
src/SimpleReactorsBVP/literature_mf1963/fixtures/mf1963_mf_low_six_state_scipy_solution.csv
src/SimpleReactorsBVP/literature_mf1963/fixtures/mf1963_mf_low_uniform_six_state_seed.csv
```

CSV seed содержит 2049 узлов и 17-значные значения `z,u0,...,u5`. Сетка
кластеризована около узкого переходного участка и построена из одного
замороженного reference-профиля. Отдельный файл задаёт равномерную mesh control
с теми же 2049 узлами. Оба генератора используют только числа/таблицы и
SciPy; код не зависит от типов KiThe runtime.

На кластеризованной mesh, с `tol=2e-7` и `max_nodes=6000`:

```text
SciPy solve_bvp 1.18.1: success, 3817 final nodes, 9 iterations
RST BVP_api::Sci: success, 3985 final nodes, 25 Newton iterations, 7 refinements
max |RST - SciPy| on the common 2049 input coordinates:
  [u0,u1,u2,u3,u4,u5] =
  [4.596e-11, 2.856e-8, 4.597e-13, 2.856e-13, 4.597e-13, 2.856e-13]
```

На равномерном seed, с теми же уравнениями, допуском и лимитом узлов, оба
маршрута достигают `MaxNodesExceeded`. Подробности и целевые команды записаны в
`src/SimpleReactorsBVP/STORY_TESTS.md`.

В текущем RST 0.4.34 точный общий seed/problem отвергается при создании обоих
solver-ов:

```text
BVP_Damp Damped: UnsupportedRoute {
  reason: "this route accepts uniform interval endpoints/count, not an arbitrary mesh"
}
BVP_CorrectionNewton: UnsupportedRoute {
  reason: "this route accepts uniform interval endpoints/count, not an arbitrary mesh"
}
```

Regression-команда для этих typed capability результатов:

```powershell
cargo test --lib --no-default-features mf1963_explicit_mesh -- --nocapture --test-threads=1
```

Positive parity gates и генерация SciPy fixtures:

```powershell
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/prepare_seed.py
uv run --with numpy --with scipy python src/SimpleReactorsBVP/literature_mf1963/six_state_parity/solve_scipy.py
cargo test --lib --no-default-features mf1963_six_state_rst_sci_ -- --nocapture --test-threads=1
```

## Requested API behavior

1. Allow a caller to provide a validated explicit mesh and a node-aligned
   initial profile. Do not silently replace or resample the supplied mesh to a
   uniform interval/count grid. Keep the existing uniform convenience route.
2. For `CorrectionNewton`, use each local interval width `h_i=z[i+1]-z[i]` in
   the fixed-mesh trapezoid residual and its Jacobian. This remains a fixed-mesh
   solve; no refinement behavior is requested for this solver.
3. For `Damped`, initialize its Newton state from the supplied mesh/profile and
   use actual local interval widths in residual/Jacobian evaluation. If adaptive
   refinement is enabled, refine the existing nonuniform mesh locally and
   preserve unaffected nodes; do not first flatten it to a uniform grid.
4. Make the node/profile shape convention explicit in the public API. For an
   explicit mesh, callers should be able to provide one state vector per mesh
   node; endpoint elimination/reduction, if required internally, belongs in
   the adapter and must preserve the coordinate-to-value pairing.
5. Validate finite strictly increasing nodes, exact endpoint/domain agreement,
   finite profile values and boundary compatibility with typed errors.
6. Add RST tests for this fixture: BVP_sci positive control; Damped and
   CorrectionNewton on the same nonuniform mesh/seed; uniform-mesh regression;
   local-step residual/Jacobian checks; and at least one Damped refinement pass
   that proves the seed mesh is not discarded. Preserve all current uniform
   mesh behavior.

The key experiment after implementation is an apples-to-apples run of
`Damped`, `CorrectionNewton`, and `BVP_sci` from the same 2049-node clustered
MF-LOW seed, followed by the same-size uniform control. This isolates solver
family behavior from the already demonstrated mesh-distribution effect.
