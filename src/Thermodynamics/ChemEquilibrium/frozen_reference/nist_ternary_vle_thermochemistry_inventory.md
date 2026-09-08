# Ternary VLE Thermochemistry Source Inventory

This document is a source-selection record for the tentative
`toluene + ethylbenzene + chlorobenzene` ternary VLE benchmark.  It is not a
frozen thermochemistry dataset and contains no coefficients or fitted closure.
Its sole purpose is to prevent a future test-only fixture from silently mixing
incompatible reference states or being calibrated from the VLE data it should
later characterize.

## Preserved production contract

`nist_toluene_ethylbenzene_chlorobenzene_preflight.rs` remains the production
repository probe.  With online lookup disabled it must keep reporting
`ValidationNotApplicable`: only toluene gas/liquid and ethylbenzene gas resolve
from the exact six requested records.  Nothing listed here is eligible for
normal `SubsData` search or production JSON libraries.

## VLE evidence kept separate

The eventual characterization rows are NIST ThermoML DOI `10.1021/je020186c`:
48 ternary liquid-composition, boiling-temperature rows at 26.66, 53.33, 79.99
and 101.32 kPa.  The payload does not publish vapour composition.  It must not
be used to construct `G0(T)`, `H0(T)`, heat capacity, or a phase-pressure
correction.

## Missing exact component states

| Component | Formula / CAS | State | Candidate authoritative route | Evidence currently visible | Closure status |
| --- | --- | --- | --- | --- | --- |
| ethylbenzene | `C8H10` / 100-41-4 | liquid | NIST Chemistry WebBook (SRD 69), condensed-phase thermochemistry and phase-change tables | liquid formation enthalpy; liquid heat-capacity data is advertised; phase-change data is available | Candidate only. Select one internally consistent compilation, then record `Cp(T)`, `H0(Tref)`, `S0(Tref)`, temperature range, and standard pressure. |
| chlorobenzene | `C6H5Cl` / 108-90-7 | gas | NIST Chemistry WebBook (SRD 69), gas-phase thermochemistry tables | gas-phase thermochemistry section and phase-change data are available | Candidate only. Directly transcribe one reviewed source route; do not combine unrelated WebBook entries by default. |
| chlorobenzene | `C6H5Cl` / 108-90-7 | liquid | NIST Chemistry WebBook (SRD 69), condensed-phase thermochemistry and phase-change tables | condensed-phase section, normal boiling point, and vaporization data are available | Candidate only. The selected data must support a bounded liquid `G0(T)` closure, not merely one normal-boiling point. |

NIST's published pages establish that the candidate evidence exists; they do
not yet establish a complete, compatible six-closure thermochemistry route.
In particular, the retrieved pages do not provide a machine-readable common
standard-pressure declaration suitable for automatic mixing with local NASA
records.

## Approved ethylbenzene numerical seeds

The first reviewed NIST WebBook seed is frozen in
`data/nist_webbook/ethylbenzene_pure_component_seed.*.json`, independently of
the ternary ThermoML rows:

| Quantity | Approved value |
| --- | --- |
| `Delta_f H0_liquid(298.15 K)` | `-12.5 +/- 0.84 kJ/mol` |
| `S0_liquid(298.15 K)` | `255.01 J/mol/K` |
| `Delta_f H0_gas(298.15 K)` | `+29.8 +/- 0.84 kJ/mol` |
| formation-enthalpy difference | `42.3 kJ/mol` |
| `Delta_vap H(294.01 K)` | `42.490 kJ/mol` |
| liquid `Cp` anchors | six values at `293.31` and `298.15 K` |
| vaporization characterization | eleven values from `294.01` to `564 K` |
| Majer-Svoboda `Delta_vap H` | `A=58.32 kJ/mol`, `beta=0.2823`, `Tc=617.1 K`, valid `295..437 K` |
| Antoine `log10(P/bar)` | `A=4.07488`, `B=1419.315 K`, `C=-60.539 K`, valid `329.74..410.27 K` |

The frozen module evaluates the Antoine and Majer-Svoboda relations only as
bounded independent oracles.  It intentionally does not interpolate the six
liquid `Cp` anchors, derive a liquid closure from the ternary VLE data, or
declare the NASA/external reference-state alignment solved.

## Non-negotiable selection gates

1. Record the primary paper/table, representation, temperature domain,
   reference temperature, reference pressure, formation convention, and
   elemental reference state for every selected closure.
2. Treat the local NASA records' `Undeclared` standard pressure as genuinely
   unknown.  A mixed production/frozen fixture is prohibited until a reviewed
   reference-state-alignment adapter proves the conversion; otherwise use a
   compatible external gas/liquid pair for that substance.
3. Prefer a reviewed polynomial.  A tabulated closure or `Cp` integration is
   permitted only when `H0(Tref)` and `S0(Tref)` make its reference state
   unambiguous.  Never fit to DOI `10.1021/je020186c`.
4. Before the ternary fixture, compare each independent
   `G0_gas(T) - G0_liquid(T)` to pure-component saturation data that is not
   derived from the ternary table.
5. The eventual frozen file must be read-only, retain semantic external
   provenance, enforce its bounds through typed errors, and enter the solver
   solely through `ResolvedThermochemistry::from_functions`.

## Primary discovery links

- NIST Chemistry WebBook, ethylbenzene (CAS 100-41-4):
  <https://webbook.nist.gov/cgi/cbook.cgi?ID=C100414&Mask=1EFF&Units=SI>
- NIST Chemistry WebBook, chlorobenzene (CAS 108-90-7):
  <https://webbook.nist.gov/cgi/cbook.cgi?ID=C108907&Mask=1D&Units=SI>
- NIST ThermoML ternary VLE payload:
  <https://trc.nist.gov/ThermoML/10.1021/je020186c.html>
