# MD trajectory archive v1

Trajectories reuse VivoArtifactStore. Coordinates are bounded binary chunk objects; immutable link objects form a chronological index; a manifest identifies a flushed prefix. No separate storage directory convention, mutable multi-gigabyte JSON array, or duplicate provenance authority is introduced.

## Binary contract

All numbers are little endian. Header: ASCII `NVMDTRJ1` (8), UInt32 version (1), UInt32 flags (bit 0 = velocities), UInt32 particle count, UInt32 frame count, UInt64 first frame ordinal, system fingerprint ASCII hex (64), configuration fingerprint ASCII hex (64): 160 bytes.

Each frame: UInt64 accepted step, Float64 time in ps, UInt32 cell-present flag, UInt32 zero reserved, nine Float64 lattice-vector components in nm, then 3*particleCount Float32 positions, then optional 3*particleCount Float32 velocities in nm/ps. Absent cells have nine zero scalars. Metadata is 96 bytes/frame; coordinate payload is 12 bytes/particle, or 24 with velocities. The codec rejects values that cannot round-trip exactly through FP32; it does not silently quantize arbitrary high-precision input. These files are trajectory samples, not restart checkpoints.

Chunks are limited to 64 MiB and 4096 frames. The writer targets 8 MiB by default, with at least one frame if it fits the hard limit. Buffered objects and encoding copies increase peak host memory beyond the on-disk chunk size, but memory does not scale with the entire run. Raw and decoded sizes are checked before coordinate-array allocation. Trajectory times and cells retain FP64.

## Publication and recovery

`append` accepts newer states of one system/configuration only. An explicit actor reservation prevents overlapping append/snapshot/finish operations across artifact-store awaits. `snapshot` flushes buffered frames and persists a restartable manifest. `finish` flushes and seals the segment. The writer advances its persisted tail only after the binary payload and link have both been stored. Failed writes may leave unreachable immutable objects, but cannot invalidate previously published prefixes.

A manifest contains counts and first/last ranges, not an ever-growing list. Each flush adds one binary chunk and one fixed-size link. Restart from a nonsealed manifest forks from its immutable prefix; another writer cannot rewrite that prefix. Sealed segments require a new segment. Reader `index` validates bounded link traversal, identity, ordinal continuity, range ordering and cycle absence. `readChunk` verifies the content hash and checks header, flags, byte lengths, frame values and the indexed range. `verify` streams through every chunk. Restart validates the link sequence and final payload; it does not claim to re-read every older coordinate payload. Use `verify` for an exhaustive archive-integrity check.

Sampled steps may be nonconsecutive. A trajectory is not an event log and does not prove that unsampled states were accurate. SHA-256 establishes byte integrity, not scientific validation or authorship. Cell vectors may change between NPT frames.

Source implementation only. No Apple build, GPU run, statistical qualification or performance benchmark was executed in this development pass.
