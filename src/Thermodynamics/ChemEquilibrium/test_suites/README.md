# Equilibrium Test Suites

This directory groups large test-only evidence suites by test-data source and
contract family. Production modules should not depend on it.

`live_data/` contains offline tests that resolve the bundled thermochemistry
repository. Additional suites belong here only when they have a distinct
fixture or evidence boundary; small unit tests should stay beside their owner
module.
