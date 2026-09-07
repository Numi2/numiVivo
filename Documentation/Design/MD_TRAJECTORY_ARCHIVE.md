# MD trajectory archive v1

Trajectories reuse VivoArtifactStore. Coordinates are bounded binary chunk objects; immutable link objects form a chronological index; a manifest identifies a flushed prefix. No separate storage directory convention, mutable multi-gigabyte JSON array, or duplicate provenance authority is introduced.

## Binary contract

All numbers are little endian. Header: ASCII `NVMDTRJ1` (8), UInt32 version (1), UInt32 flags (bit 0 = velocities), UInt32 particle count, UInt32 frame count, UInt64 first frame ordinal, system fingerprint ASCII hex (64), configuration fingerprint ASCII hex (64): 160 bytes.

Each frame: UInt64 accepted step, Float64 time in ps, UInt32 cell-present flag, UInt32 zero reserved, nine Float64 lattice-vector components in nm, then 3*particleCount Float32 positions, then optional 3*particleCount Float32 velocities in nm/ps. Absent cells have nine zero scalars. Metadata is 96 bytes/frame; coordinate payload is 12 bytes/particle, or 24 with velocities. The codec rejects values that cannot round-trip exactly through FP32; it does not silently quantize arbitrary high-precision input. These files are trajectory samples, not restart checkpoints.

Chunks are limited to 64 MiB and 4096 frames. The writer targets 8 MiB by default, with at least one frame if it fits the hard limit. Buffered objects and encoding copies increase peak host memory beyond the on-disk chunk size, but memory does not scale with the entire run. Metadata reads are limited to 64 KiB and payload reads to the validated link's exact byte count through `VivoArtifactStore.data(for:maximumBytes:)`. The rooted file reader checks the opened file size before allocating Data, independently of a descriptor's claimed size. Codec shape checks still precede coordinate-array allocation. Trajectory times and cells retain FP64.

## Publication and recovery

`append` accepts newer states of one system/configuration only. An explicit actor reservation prevents overlapping append/snapshot/finish operations across artifact-store awaits. `snapshot` flushes buffered frames and persists a restartable manifest. `finish` flushes and seals the segment. The writer advances its persisted tail only after the binary payload and link have both been stored. Failed writes may leave unreachable immutable objects, but cannot invalidate previously published prefixes.

A manifest contains counts and first/last ranges, not an ever-growing list. Each flush adds one binary chunk and one fixed-size link. Restart from a nonsealed manifest forks from its immutable prefix; another writer cannot rewrite that prefix. Sealed segments require a new segment.

`VivoMDTrajectoryArchiveReader.validate(scope:maximumChunks:)` uses one pull-based backward traversal. It retains only the current predecessor, exact remaining chunk ordinal and remaining frame ordinal. Every link is hash-verified, its system/configuration and wire-layout metadata checked, and its range compared with the newer link. The walk must terminate at precisely zero remaining chunks/frames with the declared first and last ranges. Revisited hashes identify the same bytes and thus the same chunk ordinal; they cannot match the strictly decreasing exact expected ordinal a second time. This preserves cycle rejection without an ever-growing visited set.

The validation scopes are explicit:

- `index` checks every link and structural boundary, with no coordinate payload reads.
- `restart` checks every link and the final coordinate payload.
- `allPayloads` checks every link and every coordinate payload. `verify()` uses this scope.

`readChunk` verifies the content hash and checks header, flags, byte lengths, frame values and the indexed range. Payload validation retains one wire chunk and its decoded FP64 arrays at a time. Memory is bounded independently of total trajectory length; this is not a claim that process RSS is limited to the 64 MiB wire size. Restart does not claim to re-read older coordinate payloads. Use `allPayloads` or `verify()` for exhaustive archive integrity.

Streaming validation and writer resume have no default 100,000-chunk cutoff. The finite manifest count bounds the traversal; callers can additionally supply `maximumChunks` to bound total work. The old `index(maximumChunks: 100_000)` remains an explicitly materializing convenience API returning chronological fingerprints. It retains its array bound and uses the same structural walker. Neither resume nor exhaustive verification calls that array API.

Cancellation is checked before and after each awaited archive read and before returning a completed validation result. The pull walker has no background producer or buffered stream. Cancelled or incomplete walks return no successful validation. A validation result has an inaccessible initializer, no decoding conformance, an immutable scope and manifest identity, and an internal binding to the concrete rooted-store actor. Protocol resume hands the active stage's `restart` result directly to the writer, avoiding a duplicate index walk. This immediate internal handoff cannot substitute a different store or manifest; the captured actor owns its pinned directory authority even if a pathname changes. The result records integrity observed during that traversal, not a promise against later external file changes. Subsequent archive reads continue verifying their hashes.

`md-trajectory-inspect` now performs one streaming validation and returns a UInt64 chunk count. `--verify` chooses `allPayloads`; otherwise it chooses `index`. `--maximum-chunks N` remains an optional explicit work limit, with no implicit array-allocation limit when omitted.

Sampled steps may be nonconsecutive. A trajectory is not an event log and does not prove that unsampled states were accurate. SHA-256 establishes byte integrity, not scientific validation or authorship. Cell vectors may change between NPT frames.

`MDTrajectoryArchiveValidationTests` covers bounded reads, structural corruption, scope distinctions, cancellation and materialization limits with actual rooted storage. The separate long fixture is enabled with `NUMIVIVO_LONG_ARCHIVE_TESTS=1`; it constructs more than 100,000 real links and payloads and checks streaming validation and continuation. Test-source presence does not prove that campaign ran. Record source revision, command and results from the target host before claiming long-archive qualification. Peak-memory measurements and numerical/statistical trajectory qualification remain separate evidence.

The [7 September native archive audit](../Audit/2026-09-07_ARCHIVE_STREAMING.md) records the executed 100,001-chunk campaign, extension to 100,002 chunks, and fresh-process CLI memory observations with exact source and binary identities.

After building the complete release test targets, retain a fixture and verify it through fresh CLI processes:

```sh
NUMIVIVO_LONG_ARCHIVE_TESTS=1 NUMIVIVO_TEST_ARTIFACTS=/absolute/results/native \
  swift test -c release --skip-build --no-parallel \
  --filter archiveBeyondOneHundredThousandChunksResumesAndPreservesItsPrefix
# Use the NUMIVIVO_ARCHIVE_RECEIPT path emitted by that successful campaign.
python3 Tools/Platform/check_md_archive_cli.py --binary .build/release/numivivo \
  --receipt /absolute/path/to/archive-validation-receipt.json --out /absolute/results/cli
```

The CLI checker verifies index/full-payload scopes, the extended prefix and explicit work-limit rejection, retaining every response and `/usr/bin/time -l` observation. On macOS its maximum resident set is reported in bytes. `--baseline-binary` optionally compares a retained older binary against the same 10,000- and 100,001-chunk prefixes with explicit limits. Each command runs in a fresh process, excluding fixture generation from its RSS; filesystem cache state and unrelated host workloads are not controlled, so timings alone do not establish a throughput claim.

Setting `NUMIVIVO_TEST_ARTIFACTS` to an external directory retains the long campaign's unique archive root there; otherwise its temporary fixture is removed. A successful run writes `archive-validation-receipt.json` inside that root with the 10,000-chunk prefix, original 100,001-chunk prefix and extended manifest identities, counts, validation scope and checked payload bytes, and reference-preservation result. The retained prefixes support separate fresh-process `md-trajectory-inspect --store <root> --manifest <hash> --verify` measurements with `/usr/bin/time -l`, excluding fixture generation from the measured CLI process. The test does not produce a successful receipt when its required checks fail; retained incomplete roots remain available for diagnosis.
