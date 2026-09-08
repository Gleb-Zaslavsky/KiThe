# Frozen External Reference Evidence

This test-only subtree owns the I5 external-reference layer:

## Regression policy

Every frozen-reference test, including an ignored release story, must retain
at least one executable assertion. Printed tables are evidence for human
review, never the sole failure detector. Assertions are chosen by evidence
kind: exact source identity and topology are strict; conservation,
complementarity, and independently derived oracles use numerical contracts;
external source/model differences use the fixture-specific regression envelope
recorded in `STORY_TESTS.md`. Wall time, iteration counts, and a backend's
optional diagnostic outcome remain characterization unless a test explicitly
declares them as part of its contract.

- `mod.rs` defines the strict read-only loader, provenance schema, typed row
  contracts, and structured errors;
- `loader_tests.rs` covers generic parser, metadata, schema, and read-only
  behavior;
- `iapws_tests.rs` owns the first external IAPWS water/liquid case, direct
  gas-to-liquid and liquid-to-gas lifecycle regressions, and the diagnostic
  comparison of independent, gas-side TPD, and condensed-side TPD roots;
- `iapws_water_saturation.rs` owns typed row evidence, report rendering, and
  the source/model-specific characterization and two-sided symmetry contracts
  for that IAPWS case;
- `iapws_ice_sublimation.rs` owns the corresponding typed three-way evidence
  for ice-Ih sublimation. Its validity domain and phase semantics stay
  separate from the liquid-vapor case despite the shared `(T, p)` shape;
- `iapws_ice_tests.rs` resolves the same direct local `H2O(g)/H2O(s)` fixture,
  proves the low-pressure trace path remains finite, and records the distinct
  deposition and sublimation lifecycle endpoints;
- `janaf_boudouard_thermochemistry.rs` owns the reaction-thermochemistry I5
  contract for `2 CO(g) <=> CO2(g) + C(gr)`. It intentionally precedes all
  Boudouard phase-boundary and lifecycle evidence;
- `janaf_boudouard_boundary.rs` and its companion tests own the separate
  Boudouard `P,T` boundary layer: a frozen-JANAF analytical 50/50-gas oracle,
  independent I1/I2 and canonical TPD roots, then local phase lifecycle. The
  external comparison is characterization-only because current local records
  expose their standard-state pressure as `Undeclared`;
- `data/` keeps reviewed frozen rows and provenance beside the tests that read
  them.
- `nasa_cea_h2_o2_hp.rs` owns the first full multicomponent I5 equilibrium
  characterization. It freezes NASA CEA's declared eleven-component
  `H2/O2`, `P,H` universe, reconstructs the input mixture only from local
  KiThe record data, and compares the general production route by semantic
  component identity. Its first executable layer solves the exact nine-gas
  subset because CEA reports both condensed-water rows as zero, while retaining
  those rows as typed external exclusions. This does not claim that local TPD
  rejected liquid water or ice: current condensed records do not cover the HP
  temperature. It is intentionally not a pure-phase or reaction-extent
  validator; published CEA values remain `CharacterizationOnly` until an
  explicit source/model acceptance study exists.
- `argonne_stanjan_chon.rs` owns the first general fixed-`P,T` CHON I5 case:
  a frozen 16-row Argonne/STANJAN source table and a deliberately exact
  15-species local NASA-gas universe. The published zero `C5H12` row remains
  visible as `ExternallyZeroReactantNotSolved`; the local feed is explicitly
  proven element-equivalent instead of requiring a local pentane record.
  This isolates the production 15-component, 11-reaction Gibbs formulation
  from pure-phase and scalar-reaction validation. External agreement is also
  `CharacterizationOnly` until the source/model comparison policy is reviewed.
- `nasa_tp1907_chon_graphite.rs` owns the first general heterogeneous I5
  composition case: a NASA TP-1907 Table 11.3E CHON+Ar gas universe plus an
  initially inactive pure `C(gr)` phase. Gas and condensed external amounts
  have separate types and explicitly use system-total normalisation. The
  fixture runs the normal bounded active-set solver; it does not embed a
  special carbon algorithm. Its local input follows the TP-1906 dry-air basis,
  including atmospheric `CO2`; F/A and chemical ER are independent rounded
  diagnostics. All four source rows and both continuation directions now
  characterize topology. A direct `P,T` result and an independent one-point
  typed range are required to match component-by-component before the
  multi-point continuation evidence is interpreted.
- `nasa_tp1906_chon_graphite_ph.rs` joins distinct TP-1906 Table 11.3E
  heterogeneous `H [J/g]` rows to the TP-1907 composition family only after
  checking their physical metadata. It derives total P,H targets from the
  exact executable inventory mass, characterizes the local P,T enthalpy
  reference before interpreting temperature error, and uses nested P,H as the
  conservative lifecycle reference route. The fixed-active monolithic H700
  regression now reaches the same canonical P,T graphite branch after the
  dimensionless reaction-row scaling defect was corrected; no fitted enthalpy
  offset or benchmark-specific numerical exception is used.
- `nist_benzene_toluene_vle.rs` and its companion tests own the first binary
  candidate *solution* story. The frozen NIST ThermoML dataset is P-x evidence
  only; an independent Raoult oracle never invents experimental vapor
  composition. The production fixture resolves `C6H6/C7H8` gas and
  `C6H6(L)/C7H8(L)` liquid records offline, then checks the canonical
  two-component ideal-solution TPD argmin, pressure-sign activation, transfer
  chemical-potential equality, and phase-qualified conservation. The optional
  NIST-versus-local-NASA P-x table is characterization only, not a source
  accuracy tolerance.
- `nist_toluene_ethylbenzene_chlorobenzene_preflight.rs` owns the capability
  gate for the first genuinely ternary, two-dimensional liquid-simplex
  candidate. The official NIST ThermoML payload has 48 `T-x` rows at four
  pressures and publishes liquid, but not vapor, composition. The current
  local repository lacks exact liquid ethylbenzene and both chlorobenzene
  states, so the gate returns `ValidationNotApplicable` and deliberately stops
  before freezing rows or implementing ternary TPD/lifecycle tests.
- `nist_ternary_vle_antoine_gauge.rs` owns a deliberately test-only `P,T`
  gauge for the same ternary identities. It freezes independent pure Antoine
  records and four genuine ThermoML `T/P/x` interior rows; the bounded scalar
  Raoult root derives, rather than freezes, the vapor simplex. It never enters
  production `SubsData` lookup and cannot serve a `P,H` calculation.
- `nist_ternary_vle_thermochemistry_inventory.md` records the separate source
  selection work for a possible explicit, test-only thermochemistry universe.
  It is intentionally not numerical fixture data: the production negative
  preflight remains authoritative until all missing state closures and their
  reference-state alignment are independently reviewed.
- `nist_ethylbenzene_pure_component.rs` and its frozen NIST WebBook data own
  the approved ethylbenzene numerical anchors, vaporization-enthalpy evidence,
  and independent Antoine pressure oracle. They intentionally do **not** yet
  manufacture a gas/liquid `G0(T)` closure: that needs an explicit,
  source-consistent reference-state review before it can join a ternary fixture.
- `nist_ternary_vle_antoine_gauge.rs` owns the separate, strictly test-only
  three-component `P,T` gauge.  Its three independent Antoine relations set
  `G0_gas=0` and derive only the gas/liquid Gibbs **differences** at explicit
  `p0=1 bar`; its data is independent of the ternary ThermoML rows.  It must
  never be exposed as production thermochemistry or used by `P,H` workflows.

The parent `ChemEquilibrium.rs` retains the stable internal module names using
explicit `#[path]` declarations. This organization separates I5 data evidence
from I1-I4 cross-validation while avoiding unnecessary import churn.
