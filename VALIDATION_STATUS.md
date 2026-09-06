# Validation status

## Current verified implementation

Code revision: `0a1cae6be73c7a6e12d01d9ba8bee5f2bb4eade9`.
The subsequent audit/status commit changes documentation only.

The continuation preserves the Mac mini baseline at
`a7f4017a5573d20b40b8be9738f4f8c26cd29801`. The typed Metal ABI writes,
mapped-connectivity implementation, shared nuclear budget, populated MD tests
and deepest-existing-parent output-alias checks have not been reverted.
The preserved source files and exact new source hashes are recorded in
`Documentation/Audit/barrier-convergence-observations.json`.

PASS:

- Complete release product and all test targets built from the exact committed
  source with `swift build -c release --build-tests --jobs 3 -Xswiftc -enable-testing`.
- Full Swift suite: **17/17 tests** in three suites. This includes the existing
  eight Apple execution tests and four mapped-connectivity tests, plus five new
  barrier-convergence tests containing 41 explicit checks.
- New production barrier campaigns: **100/100 plain-frame** and **106/106
  ensemble-frame** independent PySCF/CLI/cache checks. Both actual template
  commands, input equality, all energies, physical overlaps, negative scientific
  assessments and cache reconstruction were exercised.
- Real Metal pipeline compilation and execution remain passing. The separate
  product workflow passed all four isolated GPU gates and the production target
  CLI. Target execution committed 12,800 steps; maximum published fraction
  difference from FP64 was `2.0485708886353038e-5` and balance error was
  `2.0734297468405494e-5`.
- The existing nuclear reaction qualification workflow passed again on the same
  code commit, including the complete-product build and native reaction,
  thermochemistry, reference and correlated-equilibrium-solvent checks.
- No tracked source was rewritten during validation. Binary, source archive,
  input, reference, result and report SHA-256 identities were retained.

CI environment for this continuation: Apple Swift 6.1.2, arm64 macOS target,
Darwin 24.6.0; Metal reports **Apple Paravirtual device**. This is actual Metal
backend execution on the hosted runner, not a new run on the user's M4 Pro.
The independent oracle uses pinned PySCF 2.8.0, NumPy 1.26.4, SciPy 1.13.1 and
h5py 3.11.0. Production numerical execution does not depend on Python.

Successful workflows at the code revision above:

- Barrier convergence conformance: `34013510686`.
- Complete native product contracts: `34013510672`.
- Nuclear reaction qualification: `34013510661`.

## Scientific result of the new campaign

Nine explicitly mapped H3/6-31G scan geometries were evaluated with complete FCI
and nested three-, four-, five- and six-orbital CAS spaces. These are supplied
scan geometries, not optimized stationary points or the paper's reaction.

The default accuracy threshold is 0.001 Hartree for barrier/reaction differences,
relative profile error and successive-level change. At least two consecutive
**genuinely reduced** levels must satisfy it. Full-space agreement does not count
as a second reduced level.

For the five-orbital model, the plain transported core-orbital frame gave maximum
barrier error `0.1870260274894946` Hartree and relative-profile error
`0.18702602748949015` Hartree. A single uniformly weighted path-density
natural-orbital frame reduced those to `0.0008422109054371241` and
`0.0009336877697174728` Hartree respectively, without changing the criteria or
fitting the orbital frame to energy differences.

However, the adjacent four-orbital ensemble model still failed the target and
the four-to-five-orbital profile change was `0.020967628349120027` Hartree.
Both campaigns therefore correctly return **`reducedAccuracyNotEstablished`**
with no accepted reduced-level identifier. This is a successfully computed
negative scientific assessment, not a software-test failure.

Maximum native-versus-PySCF discrepancy across the tested CAS level energies was
`1.9895196601282805e-13` Hartree; across the reported barrier/reaction differences
it was `2.191580250610059e-13` Hartree. Agreement between implementations does
not establish reduced-model accuracy outside the stated reference comparison.

The ensemble policy requires full-reference CI density matrices. It is a
reference-assisted development tool, not evidence of a production speedup or
reference-free prediction. Separate ECC ladder checks retain reduced-bath
particle-closure failures rather than substituting CAS results or removing
unfavorable levels.

## Previous Mac mini validation retained

The user's Apple M4 Pro/macOS 26.6/Xcode 26.6 run at implementation revision
`60e417b78544639b4803304d671fca439eacdb69` (documented by `a7f4017...`) passed
12 tests, the real target GPU crash regression, populated/empty MD tables,
neighbor-list modes, checkpoint/restart, minimization, NVE/NVT/NPT/PME smoke
checks and production MD protocol commands. That run also reported 36 native
reaction checks, 36 independent PySCF comparisons, 21 CLI/store checks and 37
correlated-solvent checks. The new CI observations above are a separate record,
not an assertion that the current commit was rerun on that same Mac mini.

## Outstanding scientific requirements

- A stable reduced-space reaction-barrier hierarchy is **not established by the
  new H3/6-31G campaign**, despite the improved five-orbital result.
- Acrylamide-methanethiolate exact reproduction still lacks mapped author-model
  precomplex, transition-state/product structures and required numerical inputs.
- BTK exact reproduction still lacks prepared mapped snapshots and associated
  author-model settings. Independently generated examples must be labeled as
  such, not represented as supplied paper data.
- Fixed-orbital FCI/CASCI equilibrium C-PCM does not establish a global,
  physically representable ECC-DMET/PCM density functional for overlapping
  fragments.
- No kinetic rate is scientifically qualified or exported by these campaigns.
  Target simulations still retain assumed-parameter and uncertainty flags.
  Fixed-geometry energy differences do not replace saddle characterization,
  endpoint connectivity, thermal/standard-state and dynamical qualification.

No exercised build/test/CLI/Metal path failed in the current verified workflows.
This is not a claim that every possible workload, GPU kernel combination or
chemical approximation is qualified.

## Sanitizer status

The earlier explicit `-Xlinker -sanitize=address` invocation failed because that
compiler-driver flag was sent to the linker. That failure does not establish
lack of AddressSanitizer support. The corrected SwiftPM form is
`swift test -c debug --sanitize=address`; release testing may add
`--sanitize=address -Xswiftc -enable-testing`. These corrected sanitizer runs
were not performed in this continuation and are not included in PASS above.
