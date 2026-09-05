# Portable posterior checks

This directory hosts focused Linux checks for the native Swift inference and assay mathematics. The runner compiles selected production source files; it is not a full Apple package build or a Metal test.

The artifact identity adapter in the portable runner uses system OpenSSL for SHA-256 and Foundation JSON settings equivalent to the production interface. That adapter does not qualify CryptoKit, production artifact validation, or filesystem persistence. The random-number generator is extracted from the production source without editing its implementation.

Initial verification in this development pass: the existing posterior model, numerical operations, and tempered sampler compiled with Swift 6.2.1, optimization enabled, Swift 6 language mode, and warnings as errors. A 256-particle correlated Gaussian example reached beta=1 in eight stages, with means 0.35438664601606923 and 0.6525812425082524 for reference centers 0.35 and 0.65. This single run does not establish posterior calibration. Expanded checks and their exact source identities are recorded separately when executed.
