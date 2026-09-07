# Accepted sampling export — 8 September 2026

The adaptive sampler now publishes valid checkpoint references, resumes a verified completed prefix, and exports an exact selected state into the existing chemistry workflow. A real two-replica Metal fixture reached a seven-stage electronic workflow through the public CLI. Its sampling criteria remain explicitly unmet; this is a continuation and integration result, not equilibrium or chemical-accuracy qualification.

The [interface guide](../Design/MOLECULAR_SAMPLING_EXPORT.md) documents selection, typed outputs, limits, receipt verification and the remaining scientific boundaries. The [compact observations](2026-09-08-sampling-export-observations.json) bind the source archives, logs, native prefixes, executable/resource hashes and raw evidence manifest.

## Verified source and checks

Tested source is `3e9440632897e2dd58520ae8bd3674595e3d8fdb`. It passed a complete release build including all test targets in **212.08 seconds** on the physical Apple M4 Pro Mac mini, with macOS 26.6 (25G72), Xcode 26.6 (17F113) and Swift 6.3.3.

| Gate | Result |
| --- | --- |
| Focused host and native regressions | **116 executed test functions in 18 suites**, 5.654 seconds. |
| New sampling/export/electronic CLI campaign | **87 checks**, real processes and retained native input. |
| Existing general workflow CLI campaign | **61 checks**, including discovery of the new commands. |
| Target-panel workflow | **2 fresh tasks**, both successful. |
| Retained large-archive CLI campaign | **6 fresh invocations**, including the expected explicit-limit exit 65. |

Swift Testing's footer reports 117 functions including the one explicit opt-in 100,001-chunk generation test that was skipped. The six archive commands inspect the already retained 100,001-to-100,002 fixture; this campaign does not regenerate it or change its earlier qualification. The 18-suite gate includes the new reader/export/mapping/output tests, preparation, workflow cancellation, MD protocol and archive recovery, snapshot/QM/MM interfaces, and the existing compact PME/thermal/NPT and numerical-identity/resource tests. This is not a rerun of every historical suite.

| Executable | SHA-256 |
| --- | --- |
| CLI | `a245e33d4f83085d8cd90883d2f2cd26fee8ff54a1f514e8aec471c6f3dbcf3f` |
| Test executable | `8719f7d1c810183d5690b7e5b5f23e1d8a927429fc6ebfae68179b463e64898a` |

All 23 shader-resource hashes match the preceding sampling candidate. The remote checkout remained clean and both executable hashes were unchanged after testing. Native windows were coordinated with the separate Human workload; these elapsed times are run records, not performance measurements. Publication adds documentation and compact evidence to the tested implementation.

## Defect found and repaired

The first reader candidate, `32831ef4a56f7fbe955829bf74c25df0f50ea2fb`, built successfully, but one of its 15 host test functions failed. The production runner generated `molecular-sampling/<hash>/checkpoint`; the existing store rejects reference names containing slashes. This prevented the runner from publishing its initial durable cursor.

The owner now generates `molecular-sampling-<hash>-checkpoint`. Store reference validation was preserved. Source `62c971e2f937be3053e8c247f867aa1311034335` then passed all 15 reader functions and the new real native continuation fixture. Its complete build took 210.71 seconds; the host suite took 0.925 seconds and the native fixture 0.244 seconds. The original failing log remains retained rather than being replaced by the passing result.

The final integration adds a bounded reference-target read before allocation, including the preliminary named-reference lookup. Shared workflow receipt validation similarly bounds the receipt and base64 output envelopes before reading them. Tests distinguish a size-limit rejection from reading/hashing a corrupted oversized object, and verify that neither immutable verification nor a cache lookup executes replacement work.

## Actual native prefix and exact continuation

The fixture uses the explicitly synthetic two-particle harmonic model in `MolecularSamplingNativeBridgeTests`: masses 1 Da, bond length 0.125 nm, force constant 1000 in the native bonded convention, zero charge/LJ, Langevin NVT at 300 K and a 1/1024 ps timestep. Those coefficients are routing inputs, not a validated H2 force field.

Two replicas use seeds 17 and 29, initial clocks 0.125 and 0.375 ps, four equilibration steps, four blocks of eight production steps and one sample every four steps. The deliberately unmet requirement of 64 retained frames per replica prevents this short run from claiming convergence. The final receipt records `budgetExhausted` and false declared-observable criteria.

| Replica | Seed | Accepted step | Final clock, ps | Retained production frames |
| --- | --- | --- | --- | --- |
| 0 | 17 | 36 | 0.16015625 | 8 |
| 1 | 29 | 36 | 0.41015625 | 8 |

The test selects the immutable block-1 cursor actually published during the initial native run, then continues the original request. The resulting complete receipt, diagnostic identities, checkpoints and trajectory manifests must match the original block-4 result exactly. Terminal resume must produce the same receipt without advancement. A block-zero cursor remains inspectable but has no selectable production state.

The aggregate block-4 cursor is `d29abc6c0f44ac83267865e6f3f49ffc0ed13a4785493c29ebb5f4808ab1e488`. The independently retained runs on sources `62c971e` and `3e94406` have matching request, earlier cursor, final cursor, per-replica checkpoint/manifest and terminal receipt identities. This records these two executions on the same hardware and numerical profile, not general cross-device determinism. The final fixture passed in 0.199 seconds within the combined gate.

## Export and downstream route

The CLI checker copies the small retained native store and leaves the original fixture byte-identical. It selects both replicas, verifies exact checkpoint bytes and accepted clocks/velocities/cells, checks partial status and explicit convergence policy, exercises read/publication limits, and rejects malformed inputs and aliased/existing output paths. Its oldest-payload corruption case demonstrates the documented difference between final-payload and all-payload verification.

Nine typed workflow outputs retain the original checkpoint, source structure, classical system, replica configuration, all-particle snapshot, mapped structure/frame, atom mapping and provenance. They use the existing chemistry-output envelopes. The original `md-checkpoint` artifact is not assigned a conflicting kind. Host cases separately cover valid noncanonical checkpoint formatting: its original byte hash remains distinct from the canonical decoded-state hash used by the existing mapper.

Fresh verification reconstructs the source selection and exact outputs, then checks a specific immutable workflow receipt. Tests delete or repoint the mutable task cache and require successful verification without writes. A stronger all-payload verification reports its effective scope while preserving the original export receipt's final-payload scope. Task cancellation tests require no successful publication; CLI documentation separately identifies cooperative cancellation and unchanged process signal behavior.

Output destinations undergo collective validation. The CLI then captures directory authority before asynchronous publication, and later writes through those retained handles. Deterministic parent-rename and ancestor-to-symlink tests establish that replacement paths cannot redirect the prepared writes. Existing leaves remain no-clobber. Optional files precede the final receipt, but multiple file publications are not one transaction.

The checker replaces the two MD producer nodes in the existing `md-electronic-analysis` template with the exported source/system/checkpoint envelopes. All seven remaining stages execute successfully with the template's explicit H2 STO-3G and one-alpha/one-beta-electron settings. The downstream mapping exactly matches the export mapping. All seven stages reuse verified task outputs on the repeated run, and the stored workflow report verifies under the same CLI implementation identity.

This route performs electronic analysis of one accepted geometry. It does not reweight the classical samples into a QM ensemble, qualify a stationary point or reaction path, derive a barrier/rate, or complete the preparation-through-observable research campaign.

## Reproduction and retained evidence

Raw evidence resides at `/Users/n/numivivo-sampling-reader-evidence-20260908`, with a local mirror under `/Users/home/numivivo-sampling-reader-evidence-20260908`. The source-1/2/3 archives, corresponding build logs and executable/resource hash records retain the failed reader candidate, first passing native prefix and final integrated candidate. Frozen executable bundles include matching shader resources beside the binaries and inside the test bundles.

`tests-3.log`, `native-3`, `sampling-cli-3`, `workflow-cli-3`, `archive-cli-3` and `target-panel-3.json` retain the final outcomes. The raw manifest covers **1,144 files**, with SHA-256 `3451458fdf95b511f29cb4febb98d0e59c2a80e54d19804d4435954eaf4d3652`. The compact observations are an extraction record; the successful native log and source-bound raw artifacts establish the execution. Independent review checked source archives against their Git blobs, original native objects, exported envelope/payload identities, the downstream route and the raw manifest. The additional `independent-review-3.json` is recorded separately in the compact observations and is not counted in that frozen raw manifest.

To reproduce the focused prefix/export route in a fresh output directory:

```sh
swift build -c release --build-tests --jobs 3 -Xswiftc -enable-testing
export NUMIVIVO_TEST_ARTIFACTS=/absolute/new/native
export NUMIVIVO_TEST_SOURCE_COMMIT="$(git rev-parse HEAD)"
MTL_DEBUG_LAYER=1 swift test -c release --skip-build --no-parallel \
  --filter 'MolecularSampling.*Tests|KineticsPreparedOutputTests|PlatformSnapshot.*Tests'
python3 Tools/Platform/check_sampling_export_cli.py --binary .build/release/numivivo \
  --fixture /absolute/new/native/molecular-sampling-prefix-UUID/sampling-prefix-receipt.json \
  --out /absolute/new/sampling-export-cli
```

Use the unique fixture directory printed by the native test. The full recorded 18-suite selection is in the compact observations and retained `gate-3.json`. CI requires the new suites, checks their executed-pass records and runs the CLI campaign from the retained native fixture. YAML, shell blocks and embedded/checker Python were syntax-checked; no hosted CI execution is claimed here.

The [completion roadmap](../COMPLETION_ROADMAP.md) remains active. A mapped prepared-system-through-kinetic-observable scientific campaign, broader interacting ensembles, demonstrated useful scale, remaining runtime semantics and release qualification are still separate acceptance requirements.
