# DBhandlers Story Tests

These stories verify record-level correctness and diagnostics. Audit mode is
non-transactional and read-only; production calculator and `SubsData` workflows
remain fail-fast and transactional.

## Story DBH-01: Full NASA gas library audit

**Name:** Every NASA gas record is inspected and the audit is persisted

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features story_nasa_gas_library_audit_scans_and_persists_every_record -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Can the public DBhandlers audit scan a large library
record-by-record, preserve any failures without aborting, report stage timings,
and persist the complete report?

**Fixture:** The complete bundled `NASA_gas` library, using the already-open
shared `ThermoData` repository snapshot.

**Acceptance:** Every library key has exactly one audit result; valid and failed
counts sum to the catalog size; telemetry accounts for all selected records;
the serialized JSON report can be read back with every row intact. The story
does not presume that every future data revision is defect-free.

**Observed result (2026-10-26, debug):** 748 records; 748 valid, 0 failed;
payload deserialization 2.621 ms, coefficient parsing 0.407 ms, record lookup
0.531 ms, total 3.887 ms; 1 passed, 0 failed.

**Conclusion and limits:** This establishes the full-library NASA gas audit
path and its opt-in record-stage telemetry. It is not a runtime property
accuracy gate and does not cover NIST, NASA condensed, CEA, or transport
catalogs yet.

## Story DBH-02: Audit continuation at a damaged record

**Name:** A malformed middle payload does not hide later records

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features record_audit_keeps_location_and_continues_after_middle_payload_failure -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does an invalid payload retain its library, record
key, failure stage, typed category, and raw source while the next valid record
still succeeds?

**Acceptance:** The fixture produces two successful rows around one JSON-shape
failure. The failure is identified as `NASA_gas` / `broken` at
`PayloadDeserialize`, and telemetry reports 3 requested, 3 found, 2 valid, and
1 failed.

**Observed result (2026-10-26, debug):** 1 passed, 0 failed.

**Conclusion and limits:** Audit mode intentionally retains partial outcomes.
This behavior must not be used to publish a partially validated production
selection.

## Story DBH-03: Caller-owned repository lifecycle

**Name:** Closing one database view leaves other views and the shared cache valid

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features test_open_close_releases_only_the_callers_repository_view -- --nocapture --test-threads=1
```

**Question / hypothesis:** Can one caller open a repository, release its view,
and later reopen the shared catalog without invalidating another caller's
already-open snapshot?

**Acceptance:** `close()` clears the caller's local view; a second view remains
usable; a subsequent `open()` reuses the process-wide immutable repository.

**Observed result (2026-10-26, debug):** 1 passed, 0 failed.

**Conclusion and limits:** `ThermoData::open/close` manage caller-owned handles.
The process-wide default cache intentionally remains loaded for reuse; use
`try_new_from_paths` for isolated repository snapshots.

## Story DBH-04: Bundled catalog family audit

**Name:** Full-catalog audit accounts for every bundled record family

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features story_remaining_supported_catalog_families_are_auditable_without_aborting -- --include-ignored --nocapture --test-threads=1
```

**Question / hypothesis:** Can the non-transactional audit process every row
in the remaining bundled libraries, persist each family report, and continue
through malformed records without losing their locations?

**Observed result (2026-10-26, debug):**

| Library | Rows | Valid | Failed | Total ms | Decode ms | Coefficient parse ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `nuig_thermo` | 2852 | 2852 | 0 | 8.816 | 4.621 | 1.453 |
| `NASA_cond` | 382 | 378 | 4 | 1.677 | 1.065 | 0.160 |
| `CEA` | 107 | 107 | 0 | 0.414 | 0.263 | 0.000 |
| `Aramco_transport` | 1383 | 1383 | 0 | 3.280 | 2.056 | 0.000 |

The four `NASA_cond` coefficient-parse failures are `Cr(cr)`, `Fe(a)`,
`Ni(cr)`, and `ZnSO4(a)`, each reporting `Invalid temperature range in
coefficient data`. The audit persisted all 382 outcomes and completed the
subsequent library scans. NIST is not part of this bundled catalog snapshot.

**Conclusion and limits:** The public audit covers all bundled library
families and exposes pre-existing condensed-record defects precisely. This is
diagnostic evidence, not authorization to alter source data; each record needs
thermochemical provenance review before repair.

## Story DBH-05: Production record failure is contextual and transactional

**Name:** A malformed selected row reports its source and publishes no partial search state

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features production_record_parse_error_names_library_and_key_and_rolls_back -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does the fail-fast production path identify the
canonical library, exact selected record, and parsing stage while preserving
all pre-existing search results after a failure?

**Acceptance:** A malformed `NASA_gas` / `CH4` payload returns the typed
`DatabaseRecordLoadFailed` variant with stage `thermochemistry payload parse`
and retains the original calculator error as its source. A pre-existing
sentinel search result survives; no result for the malformed request escapes.

**Observed result (2026-10-26, debug):** 1 passed, 0 failed.

**Conclusion and limits:** The primary `search_substances` record-parse
boundary is contextual and transactional. This does not yet claim equivalent
coverage for every coefficient-publication or NIST retrieval route.

## Story DBH-06: Repository-open telemetry survives partial failure

**Name:** Explicit repository loading reports per-file work and the failing stage

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features explicit_repository_open_reports_file_timings_and_partial_json_failure -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does opt-in loading telemetry identify each file's
bytes and open/read/deserialize duration, then retain measurements already
collected when a later catalog file contains malformed JSON?

**Acceptance:** A successful isolated snapshot reports all three files,
schema validation, index materialization, and total duration. A malformed
second file returns a typed JSON error with its `Deserialize` failure stage
while preserving earlier file observations and total elapsed time.

**Observed result (2026-10-26, debug):** 1 passed, 0 failed.

**Conclusion and limits:** This verifies the isolated fresh-snapshot path.
For process-shared `open_with_telemetry`, a cache hit correctly reports the
cache outcome without pretending file I/O happened again; inspect a fresh
snapshot when actual load-stage timings are required.

## Story DBH-07: Shared cache telemetry distinguishes reuse from I/O

**Name:** Reopening a shared repository reports a cache hit without fictitious file timings

**Command:**

```powershell
cargo test -j 1 --lib --no-default-features shared_repository_open_telemetry_distinguishes_cache_reuse -- --nocapture --test-threads=1
```

**Question / hypothesis:** Does repeated `open_with_telemetry()` distinguish
the initial shared-cache outcome from the subsequent reuse, and does the
second open avoid reporting file work that did not occur?

**Acceptance:** Both caller-owned views are open; the second telemetry result
is `Hit`, has a total duration, and has no per-file entries. Both handles can
close without clearing the process-wide cache.

**Observed result (2026-10-26, debug):** 1 passed, 0 failed.
