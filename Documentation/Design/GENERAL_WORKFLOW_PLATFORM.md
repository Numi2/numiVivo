# General native workflow platform

The workflow recipe is a reusable interface for molecular structures, electronic
calculations, finite-drug and exposure-driven kinetics, and bounded Metal MD
segments. It is independent of any particular publication, protein or dataset.
It composes existing native engines; it does not create another numerical runtime,
artifact store, or scientific evidence classification.

## Lifecycle

`workflow-catalog` exposes the registered operations, typed input/output ports,
backend and output-validation scope. `workflow-plan` checks the entire graph,
operation versions, configuration shapes, port kinds, dependency cycles and
resource admission before numerical work or store writes. Invalid downstream
nodes cannot leave partially executed upstream jobs during static planning.

`workflow-run` executes ready independent branches in deterministic admission
batches. Each task reserves its declared numerical, input and output byte limits.
Both concurrent-task and Metal-task limits apply, and all node reservations must
fit the aggregate cap. These are **declared admission reservations**, not measured
RSS, a universal GPU allocator, or a global machine-wide limit. The numerical
engines retain their own budgets; two separately invoked recipes have separate
admission limits. Existing MD protocols remain the path for adaptive NPT and long
trajectory archives. Recipe MD nodes provide fixed-cell NVE/NVT segments with
explicit particle, mesh and step limits and conservative buffer preflight.

All execution and per-output validation goes through `VivoChemistryWorkflow`.
Identical concurrent tasks share its existing in-flight execution. A completed
stage is content-addressed before it can feed a child. Rerunning the same recipe
and store revalidates prior task outputs and executes missing stages. A changed
input, method, configuration or implementation changes its task identity; merely
changing unrelated branches does not invalidate unaffected tasks.

Numerical/resource failures are recorded. Descendants are `blocked`, while
independent branches continue. The partial run report retains every requested
node and marks unavailable exports explicitly. Exit 1 means an executed workflow
had failures; exit 65 means invalid input or a top-level operation error. A failure
is not automatically retried with a different method or a relaxed tolerance.

Cancellation is cooperative at admission/task boundaries. Already admitted work
may finish and leave individually validated cache entries. No completed recipe
report is emitted by a cancelled invocation. The core's in-flight task sharing is
not a hard real-time GPU cancellation mechanism.

## Recipe structure

The schema is `numivivo.org/workflow-recipe/v1`. A recipe contains named artifacts,
registered operation nodes, typed dependencies, named exports and a scheduling
policy. Inputs are either inline JSON or existing content fingerprints. No shell
commands, arbitrary libraries, file-system discovery or automatic downloading are
part of the recipe. Raw molecular text can be represented as a JSON string with
kind `vivo.structure-source-text` and passed to `vivo.platform.structure-import`.
Multi-record SMILES/SDF/MOL2 input is rejected rather than silently truncated.

Each node declares operation/version, input bindings, its existing native
configuration payload, and `VivoChemistryResourceContract`. Input kinds are checked
against producer kinds, not guessed from filenames. Atom order and all arrays
inside scientific payloads remain unchanged. Sorting named nodes or exports does
not change the semantic recipe identity.

The structure-to-electronic-system adapter requires an explicit conformer and
alpha/beta electron populations. It converts nanometres to Bohr and preserves
atom indices; it does not add missing hydrogens, infer protonation or turn a
periodic structure into an isolated cluster.

## Examples

```bash
swift build -c release --jobs 3
BIN=.build/release/numivivo
$BIN workflow-catalog
$BIN workflow-template molecular-analysis --output recipe.json
$BIN workflow-plan recipe.json
$BIN workflow-run recipe.json --store ./artifacts --output report.json
```

The molecular example constructs a mapped H2 electronic system, shares AO
integrals and one RHF Hamiltonian, then evaluates FCI and MP2 in independent
branches. It is an integration example, not an assertion about general chemical
accuracy. The fixed-size basis is not a recommendation for research energetics.

`workflow-template md-segments` supplies two linked harmonic NVE segments. The
second restores the first's checkpoint, including accepted-step and random
namespace identity; it does not thermalize or reset the clock. A failed segment
cannot publish a successful downstream checkpoint.

`workflow-template md-electronic-analysis` connects those segments to a mapped
snapshot, Gaussian integrals, RHF and independent FCI/MP2 calculations. The
snapshot conversion checks the original structure/system/checkpoint identities,
reconstructs atom-to-particle correspondence, excludes non-atomic virtual
particles, retains nanometre coordinates and publishes the explicit mapping.
It assigns a new conformer/structure fingerprint instead of reusing the starting
geometry identity. Periodic cells remain periodic: the isolated electronic-system
adapter rejects them rather than silently omitting periodic or MM interactions.
This example is a finite harmonic-system integration test, not an equilibrium
ensemble or validated MD force field for subsequent quantum energetics.

`Examples/workflows/target-panel.json` runs the existing synthetic exposure
fixture and a half-exposure variant using the FP64 target reference. Assumed
parameters and uncertainty flags remain in the numerical outputs. The recipe
compares conditional model behavior, not clinical predictions or fitted rates.

## Artifact and export operations

`workflow-import input.json --kind KIND --store STORE` imports bounded canonical
JSON. `workflow-export SHA256 --kind KIND` verifies content and unwraps a task
output. `workflow-verify REPORT_SHA256` binds the saved recipe/plan to the current
implementation, re-executes the workflow through verified cache entries, and
compares task identities, outputs and recorded failures. Verification may execute
previously failed stages; changed results are a new run, not the old certificate.

With `--output`, a sibling `.receipt.json` contains the run-report artifact hash.
Output and receipt paths are rejected when they alias input files or any location
inside the artifact store, including tested symlink/deep-new-parent cases. Existing
files need `--force`. Writes use the existing bounded pinned-root document I/O.
A report and its convenience receipt are separate atomic files, not a multi-file
transaction; the store remains the durable authority if file export fails.

The command-line implementation identity hashes the actual executable and records
OS, architecture, and the default Metal device identity. Library users must
provide an equally meaningful implementation fingerprint themselves. Digests do
not authenticate authorship or qualify the underlying scientific model.

## Validation boundaries

Pure new adapters deterministically reconstruct outputs. Existing chemistry
operations retain their existing numerical validators. MD cache validation checks
checkpoint identity, current numerical contract, exact clocks, source binding,
shapes and finite observables; it does **not** rerun the GPU trajectory. The catalog
records this difference explicitly. No backend is silently changed.

`PlatformWorkflowTests` covers planning, type errors, cycles, admission, actual
parallelism, duplicate-task sharing, failures/blocked descendants, restart,
corrupted inputs, and actual molecular and Metal engines. `PlatformSnapshotTests`
checks mapped MD-to-electronic execution, non-identity particle ordering,
periodicity rejection and multi-record import protection. The production CLI
check adds output aliasing, no-clobber, corrupt cached results, manifest tampering,
verified export/resume and the complete nine-stage MD-to-electronic workflow.

The dedicated read-only Apple workflow builds the complete unchanged product and
all test targets. It executes the focused platform, snapshot, real-Metal and
selected-CI suites before the actual CLI checks, and rejects a filter that did not
execute each named suite. The existing domain conformance workflows separately
retain the long full-suite chemistry/connectivity regressions. Logs and precise
source/binary identities are retained for failures too. This specification does
not itself establish a passing CI result; consult measured reports.

General platform development continues with wider registered adapters, streamed
ensemble reduction, improved long-workflow cancellation, resource measurements,
and cross-scale data transformations with explicit units and evidence. Paper
reproduction is an optional application, not the platform's acceptance criterion.
