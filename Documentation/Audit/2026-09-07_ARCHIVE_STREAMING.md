# Long MD archive recovery and inspection — 7 September 2026

A real 100,001-chunk archive passed streaming validation, resumed through the production writer, and extended to 100,002 chunks while preserving its immutable prefixes and references. Fresh CLI processes verified the retained archives. The former implicit 100,000-chunk limit no longer prevents streaming inspection or restart.

## Source and environment

- Initial archive implementation: `d9c9de49c069dcd0635acf6fbb43d006e4cd1ec0`.
- Integrated source: `b2572343f6c0a6587d13ce126068e327c79e2a77`, containing the archive implementation rebased onto nuclear/surrogate qualification revision `59e18685820f49c1a990cca11289a88abd3b9b55`, the CLI checker, and the native CI gate.
- Retained comparison source: `c8fcddc6429684c3e9a0769d7bf76152e3daff33`.
- Physical Apple M4 Pro Mac mini, reached through `ssh macmini`; macOS 26.6 build 25G72, Xcode 26.6, Swift 6.3.3. The CLI checker also ran under system Python 3.9.6.
- Raw logs, source archives, binary hashes and the retained fixture are outside Git under `/Users/n/numivivo-archive-evidence-20260907`. Compact evidence is mirrored on the coordinating MacBook. The fixture is reproducible from the committed test; its roughly 400,000 files are not committed or uploaded as CI evidence.

The archive implementation, test and checker are unchanged by integration with `59e1868`. Results below identify whether they belong to the initial or integrated source. Build and execution evidence are distinct from hosted CI status and scientific trajectory accuracy.

## Behavior and integrity

The reader traverses one predecessor at a time, checks its exact descending ordinal and frame boundary, and validates hashes, identities, shapes, byte lengths and ranges. It retains no complete link array or visited set. Read limits reach the rooted file reader before allocation. Corrupt huge claimed counts cannot trigger eager allocation.

Validation has three scopes: links alone, links plus the final restart payload, and links plus every payload. An internal immutable result binds successful restart validation to its exact manifest and rooted store actor, allowing the runner to hand it to the writer without a second traversal. Cancellation is cooperative between bounded reads; incomplete walks cannot return successful validation.

The materializing `index(maximumChunks: 100_000)` API retains its explicit array bound. Streaming validation, writer resume and CLI inspection have no unrelated implicit chunk ceiling. A supplied `maximumChunks` still rejects excessive work before traversal. Restart does not claim to reread every older coordinate payload; exhaustive verification is an explicit scope.

## Native campaign

The initial source compiled all release and test targets in 216.90 s. The routine storage suite executed 14 test functions and explicitly skipped the separately enabled long campaign. Structural corruption, bounded reads, validation scopes, rooted-store authority, cancellation and forged-count rejection passed.

With `NUMIVIVO_LONG_ARCHIVE_TESTS=1`, the named long campaign passed in 122.405 s. It constructs actual hash-addressed objects with production formats; test-only fixture writes avoid hundreds of thousands of durability synchronizations. Validation, resume and extension use production code. The successful receipt is emitted only after the required checks finish.

| Prefix | Chunks / frames | Manifest SHA-256 | Verified payload bytes |
| --- | ---: | --- | ---: |
| Benchmark | 10,000 | `d04d67d2a1cc030fbb236914a08f08a82c6d1f1170a0ff724d63b5d594d29c22` | 2,680,000 |
| Original | 100,001 | `969425ec1b83308cbdcc80be0b47726bb3cab4c77e82bdedff3bcf7b1b5a6519` | 26,800,268 |
| Extended | 100,002 | `be1ae04b398b947d62fad8484f6669a3b2d22d49da496050f8da33b376595462` | All payloads checked by CLI; the Swift receipt records restart scope. |

The campaign checks both old references, original manifest preservation, and rejection by the bounded materializing API. Its whole-process peak RSS was 70,123,520 bytes, including fixture generation and the test runner; this is not the reader's isolated allocation cost.

The [retained campaign receipt](2026-09-07-archive-validation-receipt.json) identifies every prefix and preservation check. The [CLI comparison record](2026-09-07-archive-cli-comparison.json) retains full commands, binary hashes and process observations.

## Fresh-process measurements

The initial candidate and retained comparison binary read the same one-particle, 268-byte-per-chunk fixture. Every invocation was a separate process measured with macOS `/usr/bin/time -l`; RSS is in bytes. The old binary received an explicit raised limit for the larger successful comparison. Separate invocations confirmed that its default limit rejects the valid archive, while the candidate rejects only when that limit is explicitly supplied. All nine CLI checks passed.

| Operation | 10,000 chunks: maximum RSS | 100,001 chunks: maximum RSS | Elapsed seconds, 10,000 / 100,001 |
| --- | ---: | ---: | ---: |
| Candidate, links | 10,387,456 | 10,354,688 | 1.78 / 10.85 |
| Candidate, all payloads | 10,436,608 | 10,452,992 | 3.00 / 30.13 |
| Retained baseline, all payloads | 12,713,984 | 30,703,616 | 5.27 / 44.37 |

The extended 100,002-chunk candidate passed exhaustive inspection at 10,502,144 bytes maximum RSS and 27.04 s. Candidate binary SHA-256: `89e72aa155a5de0ea47f67765c68294c8f810df5c3ed953a27263098c1dbd8ad`. Comparison binary SHA-256: `d250e2c5dce3b7388b4e9e28c74674dbebf68dd622ec8846195f38ee8e713fe7`.

These observations support memory bounded independently of archive length for this fixture, consistent with the owning traversal. They do not bound RSS for the largest allowed coordinate chunk. Filesystem cache state and unrelated CPU work were uncontrolled; elapsed times do not establish a general speedup or throughput qualification.

## Integrated-source verification

The combined `b257234` source passed a complete release/all-test-target build in 214.50 s. All six fresh archive CLI invocations then passed again against the retained fixture. The [integrated CLI record](2026-09-07-archive-cli-integrated.json) binds those results to binary SHA-256 `6451d466cba85a509aafde470c6c72ac6afb93e6d2a9cc01b5a1b2889601aab3`. Maximum RSS for exhaustive 10,000-, 100,001- and 100,002-chunk inspection was 10,600,448, 10,567,680 and 10,534,912 bytes respectively.

The final integrated native run executed **84 test functions across 12 suites**, with the opt-in long campaign explicitly skipped. Swift Testing reports 85 tests including that skip; the separately enabled long campaign above is its execution evidence. The run passed in 7.691 s with `MTL_DEBUG_LAYER=1`, external artifact retention and this exact filter:

```sh
swift test -c release --skip-build --no-parallel \
  --filter 'MDProtocolResumeTests|MDTrajectoryArchiveValidationTests|WorkflowCancellationTests|HybridRuntimeTests|AppleExecutionTests|PlatformWorkflowTests|PlatformSnapshotTests|MDQMMMInterfaceTests|PlatformQMMMWorkflowTests|RingPolymerSamplingTests|ReactiveSurrogateTests|PrecisionPlatformTests'
```

The same binary passed **57 production workflow CLI checks** under system Python 3.9.6 and completed both nodes of `Examples/workflows/target-panel.json`. The checkout remained clean and its binary hash unchanged after all gates. The final publication adds documentation and compact evidence only. No full-suite run of this combined revision is claimed: its complete build, affected native suites and CLI routes are the integrated evidence.

## Reproduction and ongoing gate

See the exact generation and CLI commands in [the archive contract](../Design/MD_TRAJECTORY_ARCHIVE.md). `Tools/Platform/check_md_archive_cli.py` retains each command, exit code, stdout, stderr, receipt identity and binary hash. It checks 10,000- and 100,001-chunk link/full scopes, exhaustive verification of the extension, and explicit work-limit rejection.

The platform workflow now requires the separately enabled long campaign and six fresh CLI checks. It retains logs and the compact success receipt, keeping the large fixture outside uploaded evidence. YAML and extracted shell syntax were checked locally; this record does not claim a hosted CI run.

The [completion roadmap](../COMPLETION_ROADMAP.md) remains active. This milestone qualifies archive structure, integrity and continuation for the stated cases. It does not establish physical trajectory accuracy, the complete MD force/ensemble matrix, biological scale, or product completion.
