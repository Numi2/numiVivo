# Single-cell count foundation: source and execution scope

The first increment is `ca5e4351d033594c5e5a8ae0bfb22a0e96962b7e`, adding
native sparse count import, quality metrics, separate normalization, replicate-
aware pseudobulk and 94 executed portable checks. The General platform workflows
run `34280325148` completed its full Apple build for that commit. Its remaining
regression/CLI steps were still running when inspected; no complete workflow
pass was inferred from the build.

The following integration increment adds strict manifest decoding, original-byte
snapshots, safe relative-file loading, full campaign composition, common artifact
storage, public commands, a common workflow adapter, a synthetic example and a
dedicated Apple qualification workflow. Local execution now passes **120 checks**
with Swift 6.2.1 on x86_64 Linux. The compiled sources include the real existing
safe filesystem implementation, whose Git blob matches
`442a8a246b46e39a37eb5d0057cc4fdc7569ee97`; no substitute reader was tested.

The [machine-readable record](2026-09-08_SINGLECELL_FOUNDATION.json) binds the
local source bytes and command. Native checks cover exact counts above 2^53,
overflow, malformed CSR/MEX, independent dense-reference comparisons, sample and
replicate semantics, missing annotations, unknown options, altered reports,
source changes, symlinks and FIFOs. A CRLF defect was found and repaired before
the first commit; text-record admission was strengthened in the second.

The CLI/artifact wrappers have only been syntax-parsed in the local Linux
workspace. `Single-cell count workflows` is the separate Apple build and public-
CLI gate; the presence of that workflow is not a successful execution result.
Its checker exercises fresh/repeated runs, no-clobber output, corrupted and
correctly hashed false artifacts, source snapshots, generic workflow execution,
export and validated cache reuse. Retain the exact run/commit when reporting its
outcome.

No real experimental dataset, differential-expression method, full single-cell
analysis package, GPU speed claim or whole-product release is qualified by this
increment. The previously reported molecular/prepared-reaction release failure
has not been diagnosed or closed by this work. The next biological milestone
needs real-data interoperability and externally evaluated analysis, using the
count and identity contracts already implemented here.
