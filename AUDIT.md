# NumiVivo audit status

The execution-critical audit found substantive defects; previous architecture descriptions must not be read as verification results.

[Execution and artifact audit, 4 September 2026](Documentation/Audit/2026-09-04_EXECUTION_AUDIT.md) documents the repairs, compatibility consequences, remaining source gaps, and Apple-silicon verification gates. [Machine-readable findings](Documentation/Audit/findings.json) distinguish implemented fixes, rejected unsupported behavior, open source work, and environment-dependent checks.

The executed portable checks passed 55 native pack/source-semantic assertions, 20 real-filesystem assertions and 12 coupling-iteration assertions. The native harness ran with AddressSanitizer and UndefinedBehaviorSanitizer. Seven selected Swift files passed syntax parsing; this is not Apple SDK typechecking. No full package build, Metal compilation, GPU execution, numerical qualification or benchmark was performed.

```sh
bash Tools/Audit/run_portable_checks.sh
```

[Results](Documentation/Audit/portable-results.json), [raw log](Documentation/Audit/portable-results.log), and [checked source identities](Documentation/Audit/checked-source-identities.json) describe exactly what was exercised. The source manifest is a partial audit workspace inventory, not a complete list of repository files.

The next acceptance target is a correctness-qualified execution release. Close the remaining typed-IR/parser, stochastic retry and global-allocation work, then complete the ordered Apple gates before expanding performance or biological claims.
