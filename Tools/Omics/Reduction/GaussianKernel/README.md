# Tiled all-anchor Gaussian MNN correction

The native MNN integration route now offers `"kernel": "tiledGaussian"`.
It retains the original exact matcher and complete Gaussian anchor sum while
accelerating its evaluation with bounded CPU tiles. Every original Hagai, Kang
and Ding cell participates. The original scalar option remains the default;
omitting `kernel` preserves historical plan encoding and scalar replay.

This follows the [failed local-kernel trial](../LocalMNN/README.md). No anchors
are omitted or selected by distance, and no bandwidth or biological margin was
changed. Complete measured biological metrics and saved neighbor/prediction
arrays match the original native MNN results exactly. Kang's four unavailable
classifier strata and Ding's partial source annotations remain unavailable.
These are inspected-cohort development results, not independent validation.

## Contract and use

Add the option to the existing `mnn` object:

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
    "maximumResidentBytes": 536870912,
    "kernel": "tiledGaussian"
  }
}
```

Use the existing commands and parent-validation path:

```sh
numivivo singlecell-pca-integrate /original-pca --plan mnn.json --output /new-mnn
numivivo singlecell-pca-integrate-verify /new-mnn
```

The implementation processes up to 32 queries and 256 anchor records per tile.
Squared distances use direct coordinate differences in the original dimension
order, without contraction/reassociation. Accelerate vector exponentials and
FP64 BLAS form the complete weighted bias sum. Anchor duplicates, source snapshots,
normalization, alignment order and panorama rules are unchanged. Floating-point
reduction order differs from the scalar implementation, so exact scalar score
bytes are not promised. The report names the selected implementation explicitly.

For total weight below `1e-280`, the implementation recomputes that row with the
original scalar exp and summation order. This preserves subnormal/underflow
behavior; it does not rescale weights or change the mathematical model. Every
extra distance term is charged to `maximumWork`, and insufficient remaining work
rejects the result. Cancellation is checked at every anchor tile, including the
fallback. A zero total preserves the original row.

The explicit extra buffers use 116,992 bytes at 20 components and at most 229,632
bytes at 64 components. This is added to latent memory admission; it excludes
Accelerate's internal workspace, allocator overhead and source metadata. Queries,
anchor coordinates and biases remain resident. Matching and the full kernel
retain their original worst-case work growth: this is **not million-cell MNN,
out-of-core integration or Metal qualification**.

## Complete-cohort evidence

The [original protocol](PROTOCOL.md), [subnormal repair](PROTOCOL_UNDERFLOW.md)
[compiler isolation](PROTOCOL_BRIDGE.md)
and [scalar owner repair](PROTOCOL_SCALAR.md) retain the implementation decisions
and their order. The initial complete trial passes, but an expanded exponent-740
case exposed a `4.5777e-5` normalized delta error in BLAS subnormal accumulation.
The original source, failed check and outputs are retained. The scalar fallback
repairs this case; all 31 expanded independent cases pass, with maximum delta
error `4.45e-16`, alongside nine rejection paths and local ASan/UBSan checks.

| Complete cohort | Cells | Production scalar (s) | Current scalar (s) | Tiled / replay (s) | Tiled / production speed ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Hagai | 13,863 | 4.505 | 3.304 | 1.854 / 1.852 | 2.43x |
| Kang | 24,673 | 30.662 | 24.807 | 12.676 / 12.675 | 2.42x |
| Ding | 44,031 | 67.914 | 54.435 | 28.241 / 28.257 | 2.40x |

These final full-cohort timings measure the
native numerical owner, excluding input decoding and artifact publication; they
are not matched end-to-end comparisons with Scanpy/Scanorama. The historical
HIRISA clustering baseline remains active. All intermediate timings are retained,
including the scalar compiler regression and unsuccessful Swift helper attempts.
The published scalar phase validates immutable anchors once and uses native C++
with the original arithmetic order;
its complete original score bytes, report and plan encoding remain mandatory.

Every original anchor byte and assembly decision matches. Maximum final score
errors versus the original scalar output are `1.51e-14` (Hagai), `1.30e-13` (Kang)
and `2.34e-14` (Ding), below the declared tolerance. Repaired scores match the
initial trial byte for byte, permitting its completed biological evaluation to
be reused through verified full-payload identity. Complete biological metrics,
4,954,020 neighbor entries and all Ding prediction arrays match the original
native results. No outcome-dependent subset replaces the original cohorts.

Five native Swift tests include legacy encoding, reference/scale preservation,
replay, memory admission and budgeted zero-weight fallback. The complete Hagai
public API publication and verification path also reconstructs its original
count/PCA parent and reproduces the qualified scores, anchors and report.
Final Hagai publication took 9.284 seconds and
verification 9.214 seconds. All three cohorts
also reproduce the same full outputs without an external BLAS thread cap.
Execution and compiler identities are retained separately from the original
implementation tag required by its historical parent receipt. Qualification of
this parent-reuse route does not relabel old receipts with a new executable hash.

## Reproduction and remaining work

The [evidence archive](../evidence/2026-09-11-gaussian-kernel) contains complete
numerical, biological and lifecycle reports, corrected scores, neighbor/prediction
arrays, failures and source. Its restoration map avoids duplicating unchanged
original anchors/scalar coordinates and identical replay outputs. Original cohort
inputs and executable binaries remain externally hash-bound.

Build the scoped library with `Tools/Omics/H5AD/build.sh`, then compile `Main.swift`
and `LifecycleMain.swift` against it, linking `OmicsHNSW.o`, `OmicsGaussian.o`,
Accelerate and libc++. `run_cohorts.py` executes all original cohorts and checks
full replay. `check_kernel.py` independently checks the C bridge. `evaluate.py`
uses the unchanged original diagnostic owners and pinned evaluation inputs;
`compare_metrics.py` compares every retained biological value and array.
Commands, input identities and the exact executed scripts are archived.

This option improves evaluation of the qualified all-anchor model; it does not
resolve HIRISA's separate annotation-retention failures or prove general
biological prediction. Next scaling work must preserve or explicitly bound
kernel error and retain all original biological gates. Independent-study
prediction still requires authoritative experimental roles and measured outcomes.

After archive verification, [duplicate-payload consolidation](../evidence/2026-09-11-gaussian-kernel-storage.json)
removed 90 identical intermediate score/anchor files and recovered 478,367,744
bytes of observed Mac mini free space. Original cohorts, the complete selected
native-phase outputs, all runtimes and the active baseline remain intact. The
receipt records every removed path, retained counterpart and hash; use it with
the archive restoration map before replaying an older intermediate artifact.
