# Portable posterior and pharmacology checks

From the repository root, with Swift 6, Python 3 and system OpenSSL development libraries on Linux:

```sh
bash Tools/Posterior/run_portable_checks.sh /tmp/numivivo-research-checks
```

Choose a new output directory. The runner compiles 15 selected production numerical source files into an optimized, testable Swift module, compiles the two test drivers separately, and retains reports, logs and source identities. It does not build the complete Apple package, compile Metal, execute a GPU, or establish biological validity.

The retained development invocation passed **94 assertions**: 72 inference/assay checks and 22 finite-pool/design/benchmark checks. Reports are under `Results/`. The final invocation reused the already successfully compiled optimized numerical library after checking that every production source and support adapter was unchanged; updated test drivers were separately compiled. This is not represented as a fresh full-tree build. The manifest records this boundary and the library digest.

The suite exercises stable censored Gaussian probabilities against high-precision references, correlated Gaussian densities against a separate covariance solve, repeated-seed correlated posterior sampling, exact stage restart, failed/omitted evaluator output, calibration/held-out separation, inferred residual measurement noise, correlated/censored prediction, finite-pool conservation, exact reaction limits, and information-gain limiting cases. There are also 24 fixed-seed prior-predictive trials with an analytic truncated-normal posterior. Their 90% intervals covered the generating parameter in 21 of 24 trials. This small descriptive result, using dependent SMC particles, is not a general simulation-based calibration certificate.

The finite-pool end state differed from the independent Radau reference by at most 5.931103484153953e-14 M in the retained eight-state fixture. A second independent DOP853 run agrees with Radau within 1e-16 M for that fixture. These are scoped numerical comparisons, not universal error bounds.

## Harness boundary

`PortableSupport.swift` supplies only the artifact interfaces needed to compile the numerical sources. SHA-256 delegates to system OpenSSL on Linux and CryptoKit where available. JSON settings match the public production convention. The adapter does not qualify production fingerprint decoding, CryptoKit implementation, signatures, artifact storage, or filesystem persistence. It is not compiled into the production package.

The default runner extracts the final public `VivoSplitMix64` declaration verbatim from `VivoAdaptiveEnsembleOptimizer.swift`; it does not compile the optimizer. A provided source excerpt can be used with `--rng-fragment` for a selected-source checkout, and that choice is recorded. The retained run used a verbatim excerpt. Modifying a source invalidates its retained source-hash match; an old report must not be promoted to a claim about changed code.

Numerical source compilation uses `-O -whole-module-optimization -swift-version 6 -warnings-as-errors -enable-testing`. Test drivers use `-Onone -parse-as-library -swift-version 6 -warnings-as-errors`. The Linux verification used Swift 6.2.1. The Darwin branch of the runner has not been executed here.

## Independent reference generation

Committed reference JSON allows ordinary checks without SciPy, NumPy or mpmath. Regeneration is explicit and overwrites the corresponding reference files:

```sh
python3 Tools/Posterior/generate_reference_data.py
python3 Tools/Posterior/generate_domain_references.py
```

The retained references used mpmath 1.3.0 at 100 decimal digits, NumPy 2.3.5 and SciPy 1.17.0. These are optional validation tools, not production dependencies. The first generator also regenerates the explicitly synthetic assay-fitting example. Neither reference generator supplies clinical measurements.

See `Examples/target-engagement/PORTABLE_RESEARCH.md` for application commands and `Documentation/Design/PORTABLE_RESEARCH_INCREMENT.md` for model boundaries.
