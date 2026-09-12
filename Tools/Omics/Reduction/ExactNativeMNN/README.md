# Exact CPU matcher integrated into the native Swift owner

`VivoMNNIntegration` now calls an allocation-free C++ exact Manhattan matcher.
The actual Swift owner produces **byte-identical scores, anchors and reports**
on all **82,567 original Hagai, Kang and Ding cells**. Plans, method identifiers,
work accounting, tie rules, panorama assembly and Gaussian options stay the same.
This is an implementation replacement, not a newly promoted biological method.

| Actual Swift owner | Cells | Seconds | Scores / anchors / report |
| --- | ---: | ---: | --- |
| Hagai | 13,863 | 2.851 | Byte-identical |
| Kang | 24,673 | 21.816 | Byte-identical |
| Ding | 44,031 | 47.551 | Byte-identical |

The entire driver takes 72.89 seconds with maximum RSS 133,906,432 bytes. These
are scoped execution measurements, not an end-to-end application, scverse or
Metal speedup. The unchanged biological evidence includes four unavailable Kang
classifier strata and partial Ding labels. No fresh biological scoring is needed
to establish byte equivalence, and none is claimed by this implementation check.
Independent validation, HIRISA preservation and million-cell MNN remain open.

## Ownership and numerical contract

[OmicsMNN.cpp](../../../../Sources/NumiVivoCore/OmicsMNN.cpp) evaluates each
cross-level pair once and updates the two existing k-wide output heaps. Ties use
(level, original source row), matching the original grouped-source-index rule.
The source-row traversal changes, but selected neighbor membership and all
published downstream arrays remain exact on the complete benchmarks. There is
no internal heap allocation or N-by-N distance array. The Swift owner retains
its current resident PCA and metadata bounds; this change is not out-of-core.

The C interface validates capacities, levels, finite inputs and the complete
pair-count × dimensions budget before matching. It handles finite-input distance
overflow separately and polls cancellation during validation and matching.
Partial outputs are invalid on failure. Swift forwards task cancellation and
requires the returned scalar-work count to equal its admitted count. The scoped
H5AD library/CLI and file-expression test linker include the new object; SwiftPM
includes the new C++ source through the existing target directory.

## Checks and reproduction

- The actual scoped library and product CLI compile. Full package/app builds and
  fresh CLI bundle publication were not run in this qualification.
- Complete owner execution checks prior source artifacts against the retained
  original manifest, then compares every output byte and report on all cohorts.
- ASan/UBSan tests pass for interleaved levels, exact ties, budget/capacity errors,
  immediate and mid-matching callback cancellation, invalid levels, nonfinite
  values and finite-input overflow.
- A real Swift task cancellation request returns `CancellationError` in about
  0.060 seconds. This is an owner task-boundary check, not precise tracing of
  which instruction received that request.

`./test.sh NEW_DIRECTORY` compiles and runs the native sanitizer tests.
`python3 verify.py` checks the archive, current owning source hashes, CLI identity,
complete output equivalence and retained failure evidence. It is evidence
verification, not a fresh owner run. The archive retains executed Swift/C++ test
sources, logs, changed owner sources and build source hashes. Full binaries,
libraries and score/anchor files remain remotely retained and hash-bound in
[manifest.json](manifest.json). [verification.json](verification.json) identifies
the source base commit and exact changed files.

The initial cohort driver requested an oversized score-reader tile and failed
before fitting. Its source, log and binary identity are retained. The corrected
driver reads 2,048-row tiles; no product admission limit was relaxed. Existing
compiler warnings remain in the archived build log.

For a fresh owner run, restore the recorded source base and changed files, run
`Tools/Omics/H5AD/build.sh OUT --with-cli`, and compile the archived `Main.swift`
against that testing-enabled module plus `OmicsHNSW.o`, `OmicsGaussian.o` and
`OmicsMNN.o`. Supply the original cohort root and a new output directory. Restore
input identities from the retained prior manifest; do not relabel historical
receipts with a new executable identity.
