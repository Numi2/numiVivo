# Native scale-aware mutual-neighbor integration

The file-backed PCA workflow now offers an explicit mutual-nearest-neighbor
(MNN) method. It corrects donor or batch effects in global median-PCA-norm units,
retains original counts and PCA, and publishes the actual anchor pairs used to
assemble corrected coordinates. Existing ridge plans remain the default.

This implements the scaled reference investigated in
[the preceding experiments](INTEGRATION_ALTERNATIVES.md). That scale was selected
after inspecting earlier failures. Evaluation on these same cohorts establishes
native implementation and measured development behavior; independent unseen-study
validation remains open. The earlier ridge and unscaled-reference failures remain
part of the evidence.

## Use and ownership

Write a plan such as:

```json
{
  "schemaVersion": 1,
  "inputKind": "fitted",
  "mnn": {
    "covariate": "donor",
    "neighbors": 20,
    "sigma": 15,
    "minimumAlignment": 0.1,
    "maximumWork": 100000000000,
    "maximumResidentBytes": 536870912
  }
}
```

```sh
numivivo singlecell-pca-integrate /original-pca --plan mnn.json --output /new-mnn
numivivo singlecell-pca-integrate-verify /new-mnn
```

`inputKind: "query"` also accepts a verified query-PCA bundle. This runs a new
transductive integration across its cells; it does not learn an unseen-donor
transform. Recursive integrated input is rejected. Use `covariate: "batch"` only
when that identity is known in source metadata. Downstream PCA-neighbor plans
must explicitly select `inputKind: "integrated"`. Both exact and HNSW graph
construction, graph clustering and embedding consume the verified result.

`VivoMNNIntegration.swift` owns numerical correction and resource admission.
`VivoPCAMNNIntegration.swift` owns file publication and reconstruction.
`VivoPCAIntegration` selects the method and retains the original parent receipt.
This option is available through the PCA-bundle route; it is not added to the
older monolithic resident-analysis plan.

Schema version 1 and encoded legacy plans/receipts remain compatible. Omitting
both method fields continues to select ridge. Supplying both method objects is
an error. Swift callers must handle the now-optional `plan.integration`,
`receipt.memberships` and `receipt.assignmentScores` fields. MNN receipts instead
contain `anchors`; no synthetic membership or assignment files are fabricated.
As before, parent and child receipts are bound to the executable identity; old
binary receipts cannot be relabeled as fresh executions.

## Numerical contract

The independently implemented native algorithm follows the pinned
[Scanorama 1.7.4](https://github.com/brianhie/scanorama) panorama assembly, whose
MIT notice is retained in `Documentation/ThirdParty/Scanorama-LICENSE.txt`.

1. Require unique cell identities, known selected covariates, and a connected
   condition/covariate design. Conditions participate only in admission. Source
   cell-type labels and RNA programs do not enter matching or correction.
2. Match exact Manhattan neighbors on row-normalized original PCA. For each
   level, searches cover the unions of preceding and following levels. Ties use
   literal level order and original source-row order within each level.
3. Retain mutual pairs. Pair alignment is the greater endpoint coverage fraction;
   select values strictly above `minimumAlignment`. Sort descending score, then
   descending level indices. Preserve panorama merge multiplicity and the
   reference's visit-count rule.
4. Divide all coordinates by one positive scalar, the median original row norm.
   Apply normalized Gaussian-weighted anchor biases with gamma `0.5 * sigma`.
   Snapshot source coordinates and biases before each correction step. Multiply
   the final result back into original PCA units.

Matching uses normalized directions, while the Gaussian correction uses the
single global scale. These serve different purposes. Uniform input rescaling is
covered by a focused invariant test. Zero-direction rows remain explicit; an
all-zero median norm is rejected. Exact kernel underflow leaves the affected row
uncorrected at that step and increments `zeroWeightCells`. Unconnected panoramas
remain separate. These states are reported rather than asserted to be successful
biological integration.

`anchors.bin` contains 16-byte little-endian records: first original row (`u32`),
second original row (`u32`), and normalized Manhattan distance (`f64`). Records
are grouped by level pair and sorted by original row indices. The report records
pair offsets, coverage, assembly order, step multiplicities, panoramas, kernel
weight extrema, zero-weight counts and RMS correction magnitude. `scores.bin`
uses the existing row/column/Double record layout. Verification rebuilds all
outputs from the original parent and compares their receipts and hashes.

## Resource limits

Counts remain sparse and file-backed. MNN's latent arrays are resident; matching
uses two k-wide heaps per cell, without an N-by-N distance matrix. This is not
out-of-core, million-cell or Metal qualification.

Admission permits 2–1,000,000 cells, 1–64 components, 2–128 levels, and 1–100
neighbors, with at least k cells per level. Those are validation bounds, not
qualified cohort sizes. Exact matching remains quadratic in cross-level cell
pairs. `maximumWork` counts cross-pair distance scalar terms plus each executed
kernel's cells × anchors × components; it is not an instruction counter.
The default is 100 billion terms, with explicit ceilings up to 10 trillion.

The base latent estimate is
`n * (3*d*8 + 2*k*16 + 80 + k*40)` bytes. Each assembly additionally admits two
anchor-coordinate buffers, `anchors * d * 16`. Metadata, JSON encoding and
source-PCA verification memory are outside this estimate. The default latent
limit is 512 MiB, with explicit ceilings up to 4 GiB. OS peak RSS is measured
separately in the native command logs.

Ding initially reached the default total-work limit during assembly. That failed
attempt is retained. Only its explicit work ceiling was raised to 1 trillion;
source, PCA, numerical parameters and biological margins were unchanged. The
runner resumes from verified prior artifacts and retains every command attempt.

## Verification and measured scope

The production executable SHA-256 is
`12d037d4993a98bca43a4d5394ad8f9b7d3f9216f086f79148875cfe512eff5a`.
The source manifest records the base commit and every production source hash.
Builds, tests and full native runs used the physical M4 Pro Mac mini; the copied
identical executable also passed local lifecycle checks.

- Release regression: 91 tests across 28 selected single-cell, H5AD, multiassay,
  perturbation and composition suites. Release product and scoped H5AD builds
  passed. The initial test-fixture compilation failure is retained separately.
- MNN lifecycle: 45 CLI commands, including 11 expected rejections. Fitted/query,
  donor/batch, exact/HNSW downstream graphs, clustering and embedding pass.
  Rehashed score, anchor and report tampering is rejected; staging is removed.
  Independent fixture reconstruction has maximum coordinate error 9.66e-15.
- Legacy lifecycle: 60 commands, including 16 expected rejections. Fixed-mode
  matrix bits match the frozen independent executable oracle; adaptive behavior
  and downstream receipt checks pass.
- Scoped invariants cover deterministic equal-distance ties, kernel underflow,
  early and late resource limits, nonfinite options and duplicate identities.

Full-cohort metrics and command measurements are retained in the
[evidence archive](evidence/2026-09-10-native-mnn/manifest.json).
The exact original 20-component PCA bits and metadata must
match earlier frozen inputs before evaluation. All native anchors and assembly
orders are compared with independently invoked Scanpy/Scanorama, and every
corrected coordinate is checked. Original 30-neighbor evaluators and acceptance
margins are reused without dropping failed or unavailable strata.

| Cohort | Cells | Original mixing excess | MNN mixing excess | Original type balanced accuracy | MNN type balanced accuracy | Measured preservation |
|---|---:|---:|---:|---:|---:|---|
| Hagai | 13,863 | 0.487372 | 0.281448 | not applicable | not applicable | all applicable gates pass |
| Kang | 24,673 | 0.057125 | 0.034965 | 0.910778 | 0.907235 | measured margins pass; four classifier strata unavailable |
| Ding | 44,031 | 0.704231 | 0.436751 | 0.660902 | 0.654672 | measured margins pass; partial source labels |

Kang NK recall is 0.936217 versus original 0.939737. Ding pDC recall is
0.328914 versus 0.369736, and cytotoxic T-cell recall is 0.630250 versus
0.673474. These decreases remain within the original 0.05 margin. Ding's
within-stratum T-receptor Spearman is 0.173048 versus 0.201722. Passing a
tolerance does not mean unchanged biological signal, and the margin headroom
for some Ding types remains small.

Every native anchor and assembly order matches the pinned independent reference.
Maximum absolute coordinate differences are 1.97e-13 (Hagai), 3.75e-11 (Kang)
and 4.57e-13 (Ding); anchor-distance errors are at most 2.23e-15. No step in
these three full runs has zero kernel weights. The numerical reference uses
NumPy/scikit-learn/Scanorama, separately from native Swift scalar reductions.
All global and stratified 30-neighbor indices also match the frozen reference
exactly on all three cohorts; both sets of arrays remain retained.

| Native publication command | Wall seconds | Peak RSS bytes | Matching plus kernel scalar terms |
|---|---:|---:|---:|
| Hagai MNN | 18.10 | 229,244,928 | 6,197,613,240 |
| Kang MNN | 36.67 | 227,459,072 | 46,417,554,220 |
| Ding MNN, raised admission | 80.12 | 333,496,320 | 104,540,599,940 |

These are individual CPU publication measurements, including parent verification
and artifact writing; they are not repeated timing estimates or a GPU/scverse
speedup comparison. Full reconstruction took 18.18, 36.86 and 81.08 seconds,
respectively. The native runner records 12 commands, including the retained Ding
budget failure and Baron confounding rejection. A separate 13th full-source
command rejects Ding's unreported donor identity without publishing output.

Kang's four unavailable condition-classifier strata prevent a complete pass.
Ding's assigned source labels cover 29,411 of 44,031 UMI-method cells; all cells
remain in correction and neighborhood/program evaluation. Source labels are
evaluation evidence, not independent experimental truth. Hagai has one source
type and does not establish rare-type preservation. This method remains opt-in
pending validation on an untouched independent cohort.

## Reproduce the retained checks

From this directory, verify all compressed and uncompressed evidence hashes:

```sh
python archive_adaptive_integration.py --verify --out evidence/2026-09-10-native-mnn
```

The verified archive contains 1,385 logical files in 335 unique objects
(215,375,112 raw logical bytes; 60,249,687 compressed object bytes), plus 12
human-readable summaries. Initial compilation and resource-admission failures
remain available alongside the successful runs.

The manifest maps each logical file to a deterministic gzip object. Decompress
that object to restore the logical path; verify its raw SHA-256 and byte count.
The original H5AD files, executable, HDF5 library and frozen legacy executable
remain external with identities recorded in `externalInputs`. The prior
Scanorama reference artifacts remain in the separately published alternatives
archive, whose manifest and all consumed reference files were rechecked.

With those exact inputs restored, `run_native_mnn.py` performs source hashing,
fresh native PCA, integration, full replay and exact non-H5AD artifact transfer.
The retained original and resumed cohort specifications distinguish the failed
Ding budget from the successful explicit ceiling. `--resume` rechecks completed
artifacts and continues an interrupted run without overwriting its logs.
`check_full_mnn.py` compares each native cohort with the frozen reference and
reruns the original biological evaluators. `check_mnn_integration.py` exercises
the fixture lifecycle using Scanorama 1.7.4; the preceding reference report
records its isolated dependency installation. `check_file_integration.py` retains
the unchanged legacy oracle checks. All scripts expose their required paths
through `--help`.
