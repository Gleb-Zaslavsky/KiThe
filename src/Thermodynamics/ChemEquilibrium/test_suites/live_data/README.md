# Live-Data Regression Suite

This directory owns production-evidence tests that resolve bundled local
thermochemistry rather than using synthetic coefficient fixtures.

`equilibrium_live_data_tests.rs` is deliberately kept intact in this first
filesystem pass. Its established test names, ignored release commands, and
shared repository/JSON-integrity helpers remain stable. It covers real lookup
and provenance, fixed `P,T`, `P,H`, range continuation, backend matrices,
phase transitions, Jacobian checks, timings, and library immutability.

It is not a cross-validation-only suite: K_eq and phase-boundary checks are
only independent evidence inside a much wider canonical production path.

Future decomposition should extract a local `support.rs` and split the file
into catalog/lookup, P,T, P,H, range/backend, and phase-transition modules.
That work should happen only with a conscious decision about retaining or
renaming the current test IDs.
