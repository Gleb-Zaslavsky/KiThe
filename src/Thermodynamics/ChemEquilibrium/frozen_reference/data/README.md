# Frozen Reference Data

## NIST ThermoML benzene/toluene VLE

`nist_thermoml/benzene_toluene_vle_353_15k.*.json` is a reviewed five-row
P-x subset from the official ThermoML JSON for DOI `10.1016/j.fluid.2009.03.012`.
It preserves the published liquid benzene composition and total pressure at
353.15 K, including the two pure endpoints. The source is P-x evidence only:
no experimental vapor composition is synthesized or stored here.

This directory stores small, reviewed, read-only external reference datasets
used by the chemical-equilibrium validation suite. It is colocated with the
test-only frozen-reference loader and stories so data, schema, and provenance
remain one reviewable unit.

Each dataset keeps numerical rows separate from provenance metadata. Future
authoritative data belongs under source-specific directories such as
`iapws/`, `janaf/`, `atct/`, or `nasa_cea/`. Regression tests must never
download, refresh, or rewrite these files.

The `synthetic/` fixture verifies the loader and validator only. It is not I5
physical evidence and must not be cited as agreement with an external source.

`iapws/` currently contains pressure-boundary data. `janaf/` contains primary
species thermochemistry for the Boudouard reaction; it deliberately does not
store a pre-derived reaction constant so a review can trace every diagnostic
quantity back to the external CO/CO2 columns.

`nasa_cea/` contains complete published-equilibrium characterizations. Its
first `H2/O2` HP case stores the exact CEA-declared product universe and the
rounded published state by species identity. It does not contain a KiThe
solution, regenerated CEA values, or a hidden species-selection policy. The
current executable comparison declares its nine-gas subset in Rust and keeps
the two zero condensed rows visible as excluded external evidence; the frozen
JSON itself remains the complete eleven-row source transcription.

`argonne_stanjan/` contains a reviewed fixed-`P,T` CHON composition table from
Argonne's STANJAN comparison. The original pentane/methane/air input basis,
element totals, and zero-pentane output row remain in JSON. The executable
15-species local universe and its element-equivalent feed are case logic in
the accompanying Rust adapter, never generated reference data.

`nasa_tp1907/` contains heterogeneous composition rows from NASA TP-1907 Table
11.3E. The files keep printed gas and condensed components separately because
their source normalisation is a physical contract, not a presentation detail.
The reviewed rows currently support the graphite topology story and
characterization; they do not contain generated KiThe values or a hidden
full-55-species NASA universe.

`nasa_tp1906/` stores the separate Table 11.3E heterogeneous equilibrium
specific-enthalpy rows in their published `J/g` unit. They are not a molar
thermochemistry table and do not duplicate TP-1907 composition. The P,H
fixture joins both sources by explicit physical-family fields and converts
`J/g` to a total target only from its reviewed local closed-inventory mass.

## Dataset Convention

Every reviewed dataset uses two files with one shared `dataset_id` and a
positive `dataset_format_version` owned by KiThe. This format version is
separate from the cited source release:

- `<name>.metadata.json` contains source identity, citation, transcription,
  column meanings, units, and optional source uncertainty or precision;
- `<name>.rows.json` contains only the dataset identity and typed numerical
  rows.

The metadata `data_file` must name the rows file exactly. A Rust row type owns
the expected columns, units, and physical row validation. Comparison tolerance
belongs to the individual validation story, not to this loader.

## Manual Review Workflow

1. Select an authoritative publication or maintained reference correlation and
   record its stable identifier, release/version, and exact table identity.
2. Transcribe a small reviewable subset without converting units silently.
   Declare the stored unit for every column, format version, source precision,
   and record any conversion or correlation evaluation in the transcription
   statement.
3. Have a second review compare every frozen value with the cited source.
4. Add a typed row schema, structural tests, and a source-justified comparison
   policy. The comparison failure must print dataset, source, version, and row.
5. Commit changes as an ordinary reviewed code/data change.

There is deliberately no refresh command. Do not generate these files from
KiThe output, fetch them during tests, or overwrite them after a comparison.
