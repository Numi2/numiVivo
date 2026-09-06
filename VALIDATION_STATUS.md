# Validation status

Validated source revision before this note: `60e417b78544639b4803304d671fca439eacdb69` on Apple M4 Pro, macOS 26.6, Xcode 26.6.

PASS:

- Release product build and test build with testable module access.
- Full Swift suite: 12 tests in 2 suites passed.
- Real Metal target-engagement test and production CLI: 12,800 committed steps, FP64 comparison passed, exit 0.
- Metal pipeline compilation, populated target ABI records, MD empty/populated tables, neighbor-list off/on, checkpoint/restart, minimization, NVE, NVT, NPT and PME smoke qualification.
- Native reaction qualification: 36 checks and 36 independent PySCF comparisons.
- Correlated FCI/CASCI C-PCM density/reaction-field closure, cache reuse, and equilibrium full-CI nuclear checks: 37 checks passed.
- Mapped H3 connectivity: four independent displacement/step trials, twelve pair comparisons, explicit mapped distinct endpoints.
- Production reaction CLI/store checks: 21 checks passed.

FAIL:

- No exercised Apple-native build, test, Metal, MD, reaction, CLI, or cache path remains failing at this revision.

BLOCKED ON INPUT:

- Acrylamide-methanethiolate reproduction: mapped precomplex, transition-state, product geometries and required paper inputs are not present.
- BTK reproduction: prepared BTK snapshots and associated mapped structures are not present.

NOT SCIENTIFICALLY QUALIFIED:

- No kinetic rate is exported or qualified. The target result uses declared synthetic/assumed parameters, and the reaction harmonic barrier remains explicitly non-kinetic.
- The reaction fixtures are finite-basis conformance cases, not the paper reaction. Correlated solvent checks cover fixed-orbital FCI/CASCI C-PCM; they do not establish global ECC-DMET/PCM self-consistency.

TOOLING LIMITATION:

- The requested debug AddressSanitizer link form was rejected by this Xcode linker (`unknown options: -sanitize=address`); the optimized Metal fault was instead localized with LLDB and fixed at the typed-buffer ABI write boundary.
