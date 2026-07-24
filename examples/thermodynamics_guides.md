# Thermodynamics Guides

This file collects the stable user-facing thermodynamics workflows that are
covered by tests and examples.

## Supported Library Formats

- `NASA`
- `NIST`
- `CEA`
- `TRANSPORT`

Each format carries its own schema and property expectations. The supported
lookup path is intentionally typed: the caller selects a library family and
then asks for the properties that family can supply.

### Required Fields At A Glance

- `NASA`: polynomial coefficients, composition, and a model tag.
- `NIST`: interval bounds plus interval coefficients for the selected fit.
- `CEA`: calculator-specific thermochemical or transport payloads.
- `TRANSPORT`: transport-property coefficients plus the phase/context data
  used to evaluate them.

## Units And Contract

- Temperature is expected in kelvin.
- Pressure is expected in pascal unless the API call explicitly documents a
  different unit.
- Molar mass is expected in `kg/mol` or `g/mol` depending on the calculator
  entry point.
- Transport calculations require the usual prerequisite state: pressure,
  temperature, molar mass, and composition-derived inputs.
- Finite physical inputs are required. NaN, infinities, zero where invalid,
  and negative values are rejected before derived calculations run.
- Temperature windows are library-specific, so callers should treat the
  documented fit range as part of the contract rather than a hint.

## Typical Workflows

### Local-library lookup

Use the canonical lookup API to resolve a substance from local catalog data,
then inspect the typed search report.

See:
- `examples/thermodynamics_guide_canonical_accessors.rs`
- `src/Thermodynamics/User_substances_tests2.rs`

```rust,ignore
let mut subs_data = SubsData::new();
subs_data.set_substances(vec!["H2O".to_string(), "CO2".to_string()]);
subs_data.set_library_priority_id(LibraryId::NasaGas, LibraryPriority::Priority);
subs_data.search_substances()?;

let report = subs_data.search_summary_report();
println!("{}", report.render_table());
```

### Explicit mixed-source lookup

Explicit instructions can route one substance to one library while the rest
use the normal priority search order. The canonical search report preserves the
source and the property category separately.

```rust,ignore
subs_data.set_explicit_search_instruction("H2O".to_string(), LibraryId::NasaGas);
subs_data.set_library_priority_id(LibraryId::AramcoTransport, LibraryPriority::Priority);
subs_data.search_substances()?;
```

### Numeric values

After lookup, numeric caches publish values such as `Cp`, `dH`, `dS`,
`Lambda`, and `Visc` through explicit readiness accessors.

```rust,ignore
let cp_state = subs_data.therm_value_state("H2O", DataType::Cp);
assert!(cp_state.is_ready());
```

### Closures

Temperature-dependent closures are exposed as read-only cached callables. The
readiness state is explicit, so callers do not need to inspect placeholder
values.

```rust,ignore
let lambda_fun_state = subs_data.transport_function_state("CO", DataType::Lambda_fun);
if let Some(lambda_fun) = lambda_fun_state.value() {
    println!("{}", lambda_fun(400.0));
}
```

### Symbolic expressions

Symbolic results are published as `Expr` values and can be pretty-printed or
lambdified when needed.

```rust,ignore
let sym_state = subs_data.transport_symbolic_state("CO", DataType::Lambda_sym);
if let Some(expr) = sym_state.value() {
    println!("{}", expr.pretty_print());
}
```

### Diffusion coefficients

Binary diffusion calculations are built from the transport side and exposed as
numeric, closure, and symbolic forms.

```rust,ignore
subs_data.initialize_diffusion_data()?;
subs_data.calculate_all_pairs()?;
subs_data.calculate_all_closures()?;
subs_data.calculate_all_symbolic()?;
```

## Batch Semantics

- Fail-fast batch methods publish results only after the full loop succeeds.
- Best-effort paths remain diagnostic-only and should not replace a valid cache
  with a partial one.

## Repository Loading

- Repository loading is explicit and can be pointed at temporary fixture files
  in tests.
- Reloading is not hidden behind a process-global side effect.
- Offline NIST behaviour is part of the default contract; live network tests
  should remain feature-gated or explicitly opted in.
