# Native reference-mapping regression (2026-09-14)

The published `5537fb5f34ab543b5a672596c3ad8d46d0ae4f06` source was built on the
physical M4 Pro Mac mini and `swift test --filter SingleCellReferenceTests` ran
all three tests successfully. The compressed log and execution record retain the
host, command, source revision and hashes.

This is a native software regression for the explicit reference-mapping novelty
gates. It does not measure held-out annotation accuracy, calibrate thresholds,
validate biological identities or establish prediction of a biological outcome.
