import Foundation

public struct VivoMDTrajectoryManifest: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/md-trajectory-manifest/v1"
    public let schema: String
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let particleCount: UInt32
    public let includeVelocities: Bool
    public let frameCount: UInt64
    public let chunkCount: UInt64
    public let tail: VivoFingerprint?
    public let firstStep: UInt64?
    public let lastStep: UInt64?
    public let firstTimePS: Double?
    public let lastTimePS: Double?
    public let sealed: Bool

    public func validate() throws {
        _ = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: particleCount, includeVelocities: includeVelocities)
        guard schema == Self.schemaID, chunkCount <= frameCount else { throw bad("trajectory manifest schema/count mismatch") }
        if frameCount == 0 {
            guard chunkCount == 0, tail == nil, firstStep == nil, lastStep == nil,
                  firstTimePS == nil, lastTimePS == nil else { throw bad("empty trajectory has nonempty metadata") }
        } else {
            guard chunkCount > 0, tail != nil, let firstStep, let lastStep, firstStep <= lastStep,
                  let firstTimePS, let lastTimePS, firstTimePS.isFinite, lastTimePS.isFinite,
                  firstTimePS >= 0, firstTimePS <= lastTimePS,
                  (frameCount == 1 || (firstStep < lastStep && firstTimePS < lastTimePS)) else {
                throw bad("trajectory manifest range is invalid")
            }
        }
    }
}

public struct VivoMDTrajectoryChunkLink: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/md-trajectory-link/v1"
    public let schema: String
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let particleCount: UInt32
    public let includeVelocities: Bool
    public let chunkOrdinal: UInt64
    public let firstFrameOrdinal: UInt64
    public let frameCount: UInt32
    public let previous: VivoFingerprint?
    public let payload: VivoFingerprint
    public let payloadBytes: UInt64
    public let firstStep: UInt64
    public let lastStep: UInt64
    public let firstTimePS: Double
    public let lastTimePS: Double
}

public enum VivoMDTrajectoryValidationScope: String, Sendable, Equatable {
    /// Every index link and its structural relationship; no coordinate payloads.
    case index
    /// Every index link and the newest coordinate payload, as required by resume.
    case restart
    /// Every index link and every coordinate payload.
    case allPayloads
}

/// Evidence of one completed traversal, not a promise against later external
/// filesystem changes. The private reader binding retains the actual rooted
/// store authority, not merely a path that could name a different directory.
/// There is intentionally no public initializer or decoding conformance.
public struct VivoMDTrajectoryValidation: Sendable {
    public let manifestFingerprint: VivoFingerprint
    public let scope: VivoMDTrajectoryValidationScope
    public let indexedChunks: UInt64
    public let verifiedPayloads: UInt64
    public let verifiedPayloadBytes: UInt64
    fileprivate let reader: VivoMDTrajectoryArchiveReader

    fileprivate init(reader: VivoMDTrajectoryArchiveReader, scope: VivoMDTrajectoryValidationScope,
                     indexedChunks: UInt64, verifiedPayloads: UInt64, verifiedPayloadBytes: UInt64) {
        self.reader = reader; manifestFingerprint = reader.manifestFingerprint; self.scope = scope
        self.indexedChunks = indexedChunks; self.verifiedPayloads = verifiedPayloads
        self.verifiedPayloadBytes = verifiedPayloadBytes
    }
}

/// Validation retains one link and at most one decoded chunk at a time. The
/// materializing index() remains a separate explicitly bounded convenience API.
public struct VivoMDTrajectoryArchiveReader: Sendable {
    public let store: VivoArtifactStore
    public let manifest: VivoMDTrajectoryManifest
    public let manifestFingerprint: VivoFingerprint

    private init(store: VivoArtifactStore, manifest: VivoMDTrajectoryManifest, manifestFingerprint: VivoFingerprint) {
        self.store = store; self.manifest = manifest; self.manifestFingerprint = manifestFingerprint
    }

    public static func open(store: VivoArtifactStore, manifest fingerprint: VivoFingerprint) async throws -> Self {
        try Task.checkCancellation()
        let data = try await store.data(for: fingerprint, maximumBytes: 64 * 1024)
        try Task.checkCancellation()
        let value = try VivoCanonicalJSON.decode(VivoMDTrajectoryManifest.self, from: data)
        try value.validate()
        try Task.checkCancellation()
        return .init(store: store, manifest: value, manifestFingerprint: fingerprint)
    }

    /// Returns chronological link identities. This API deliberately allocates
    /// O(chunkCount) output; streaming validation and resume never call it.
    public func index(maximumChunks: UInt64 = 100_000) async throws -> [VivoFingerprint] {
        try Task.checkCancellation()
        guard manifest.chunkCount <= UInt64(Int.max) else { throw bad("trajectory index cannot fit in an array") }
        var traversal = try VivoMDTrajectoryIndexTraversal(reader: self, maximumChunks: maximumChunks)
        var result: [VivoFingerprint] = []
        while let entry = try await traversal.next() { result.append(entry.fingerprint) }
        try Task.checkCancellation()
        return result.reversed()
    }

    /// Traverses the finite manifest count without a length-dependent allocation.
    /// An optional caller work limit rejects before any link read. The default
    /// has no unrelated chunk-count ceiling; cancellation remains cooperative at
    /// each metadata/payload read. Only a complete traversal returns evidence.
    public func validate(scope: VivoMDTrajectoryValidationScope,
                         maximumChunks: UInt64? = nil) async throws -> VivoMDTrajectoryValidation {
        try Task.checkCancellation()
        var traversal = try VivoMDTrajectoryIndexTraversal(reader: self, maximumChunks: maximumChunks)
        var indexed: UInt64 = 0, payloads: UInt64 = 0, bytes: UInt64 = 0
        while let entry = try await traversal.next() {
            if scope == .allPayloads || (scope == .restart && indexed == 0) {
                _ = try await readChunk(entry.link)
                let total = bytes.addingReportingOverflow(entry.link.payloadBytes)
                guard !total.overflow else { throw bad("verified trajectory byte count overflow") }
                bytes = total.partialValue; payloads += 1
            }
            indexed += 1
        }
        try Task.checkCancellation()
        return .init(reader: self, scope: scope, indexedChunks: indexed,
                     verifiedPayloads: payloads, verifiedPayloadBytes: bytes)
    }

    public func readLink(_ fingerprint: VivoFingerprint) async throws -> VivoMDTrajectoryChunkLink {
        try Task.checkCancellation()
        let data = try await store.data(for: fingerprint, maximumBytes: 64 * 1024)
        try Task.checkCancellation()
        let value = try VivoCanonicalJSON.decode(VivoMDTrajectoryChunkLink.self, from: data)
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: manifest.particleCount,
                                                                includeVelocities: manifest.includeVelocities)
        guard value.schema == VivoMDTrajectoryChunkLink.schemaID,
              value.systemFingerprint == manifest.systemFingerprint,
              value.configurationFingerprint == manifest.configurationFingerprint,
              value.particleCount == manifest.particleCount, value.includeVelocities == manifest.includeVelocities,
              value.frameCount > 0, value.frameCount <= UInt32(VivoMDTrajectoryChunkCodec.maximumFramesPerChunk),
              UInt64(value.frameCount) <= UInt64.max - value.firstFrameOrdinal,
              value.payloadBytes == UInt64(VivoMDTrajectoryChunkCodec.headerBytes) + UInt64(value.frameCount) * UInt64(stride),
              value.payloadBytes <= UInt64(VivoMDTrajectoryChunkCodec.maximumChunkBytes),
              value.firstStep <= value.lastStep, value.firstTimePS.isFinite, value.lastTimePS.isFinite,
              value.firstTimePS >= 0, value.firstTimePS <= value.lastTimePS,
              (value.chunkOrdinal == 0) == (value.previous == nil) else { throw bad("invalid trajectory link contract") }
        try Task.checkCancellation()
        return value
    }

    public func readChunk(_ linkFingerprint: VivoFingerprint) async throws -> [VivoMDArchiveFrame] {
        try await readChunk(readLink(linkFingerprint))
    }

    private func readChunk(_ link: VivoMDTrajectoryChunkLink) async throws -> [VivoMDArchiveFrame] {
        try Task.checkCancellation()
        // readLink() already bounded this field to the codec's maximum. Apply
        // its exact byte count to the opened file before allocating any Data.
        let data = try await store.data(for: link.payload, maximumBytes: Int(link.payloadBytes))
        try Task.checkCancellation()
        guard UInt64(data.count) == link.payloadBytes else { throw bad("trajectory payload byte count mismatch") }
        let frames = try VivoMDTrajectoryChunkCodec.decode(data,
            systemFingerprint: manifest.systemFingerprint, configurationFingerprint: manifest.configurationFingerprint,
            particleCount: manifest.particleCount, includeVelocities: manifest.includeVelocities,
            firstFrameOrdinal: link.firstFrameOrdinal)
        guard frames.count == Int(link.frameCount), frames.first?.stepIndex == link.firstStep,
              frames.last?.stepIndex == link.lastStep, frames.first?.timePS == link.firstTimePS,
              frames.last?.timePS == link.lastTimePS else { throw bad("trajectory payload range disagrees with index") }
        try Task.checkCancellation()
        return frames
    }

    /// Exhaustive payload verification is streaming. Callers may explicitly cap
    /// total work without changing the default ability to verify a long archive.
    public func verify() async throws {
        _ = try await validate(scope: .allPayloads)
    }
    public func verify(maximumChunks: UInt64) async throws {
        _ = try await validate(scope: .allPayloads, maximumChunks: maximumChunks)
    }
}

/// A pull walker shared by index(), validation and resume. No detached producer,
/// buffered sequence or growing visited set is needed. Revisited hashes have
/// identical verified bytes and therefore the same chunkOrdinal; they cannot
/// match the strictly decreasing EXACT expected ordinal a second time.
struct VivoMDTrajectoryIndexTraversal: Sendable {
    struct Entry: Sendable {
        let fingerprint: VivoFingerprint
        let link: VivoMDTrajectoryChunkLink
    }
    private let reader: VivoMDTrajectoryArchiveReader
    private var nextFingerprint: VivoFingerprint?
    private var chunksRemaining: UInt64
    private var endOrdinal: UInt64
    private var newer: VivoMDTrajectoryChunkLink?

    init(reader: VivoMDTrajectoryArchiveReader, maximumChunks: UInt64? = nil) throws {
        try reader.manifest.validate()
        if let maximumChunks, reader.manifest.chunkCount > maximumChunks {
            throw bad("trajectory traversal exceeds caller's chunk limit")
        }
        self.reader = reader; nextFingerprint = reader.manifest.tail
        chunksRemaining = reader.manifest.chunkCount; endOrdinal = reader.manifest.frameCount
    }

    mutating func next() async throws -> Entry? {
        try Task.checkCancellation()
        guard let fingerprint = nextFingerprint else {
            guard chunksRemaining == 0, endOrdinal == 0 else { throw bad("trajectory prefix is truncated") }
            if let first = newer {
                guard first.firstStep == reader.manifest.firstStep, first.firstTimePS == reader.manifest.firstTimePS else {
                    throw bad("trajectory first range disagrees with manifest")
                }
            }
            return nil
        }
        guard chunksRemaining > 0 else { throw bad("trajectory link cycle or excess links") }
        let link = try await reader.readLink(fingerprint)
        try Task.checkCancellation()
        guard link.chunkOrdinal == chunksRemaining - 1,
              UInt64(link.frameCount) <= endOrdinal,
              link.firstFrameOrdinal == endOrdinal - UInt64(link.frameCount) else {
            throw bad("trajectory link ordinal gap")
        }
        if let newer {
            guard link.lastStep < newer.firstStep, link.lastTimePS < newer.firstTimePS else {
                throw bad("overlapping trajectory chunk ranges")
            }
        } else {
            guard link.lastStep == reader.manifest.lastStep, link.lastTimePS == reader.manifest.lastTimePS else {
                throw bad("trajectory tail range disagrees with manifest")
            }
        }
        nextFingerprint = link.previous; chunksRemaining -= 1
        endOrdinal = link.firstFrameOrdinal; newer = link
        return .init(fingerprint: fingerprint, link: link)
    }
}

public actor VivoMDTrajectoryArchiveWriter {
    private let store: VivoArtifactStore
    private let system: VivoFingerprint
    private let configuration: VivoFingerprint
    private let particles: UInt32
    private let velocities: Bool
    private let framesPerChunk: Int
    private var buffered: [VivoMDArchiveFrame] = []
    private var persistedFrames: UInt64 = 0
    private var chunkCount: UInt64 = 0
    private var tail: VivoFingerprint?
    private var firstStep: UInt64?
    private var lastStep: UInt64?
    private var firstTime: Double?
    private var lastTime: Double?
    private var sealed = false
    private var busy = false

    public init(store: VivoArtifactStore, systemFingerprint: VivoFingerprint,
                configurationFingerprint: VivoFingerprint, particleCount: UInt32,
                includeVelocities: Bool = false, targetChunkBytes: Int = 8 * 1024 * 1024) throws {
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: particleCount, includeVelocities: includeVelocities)
        guard targetChunkBytes > VivoMDTrajectoryChunkCodec.headerBytes,
              targetChunkBytes <= VivoMDTrajectoryChunkCodec.maximumChunkBytes else { throw bad("invalid trajectory chunk budget") }
        self.store = store; system = systemFingerprint; configuration = configurationFingerprint
        particles = particleCount; velocities = includeVelocities
        framesPerChunk = min(VivoMDTrajectoryChunkCodec.maximumFramesPerChunk,
                             max(1, (targetChunkBytes - VivoMDTrajectoryChunkCodec.headerBytes) / stride))
    }

    private init(store: VivoArtifactStore, manifest: VivoMDTrajectoryManifest, framesPerChunk: Int) {
        self.store = store; system = manifest.systemFingerprint; configuration = manifest.configurationFingerprint
        particles = manifest.particleCount; velocities = manifest.includeVelocities
        self.framesPerChunk = framesPerChunk; persistedFrames = manifest.frameCount
        chunkCount = manifest.chunkCount; tail = manifest.tail
        firstStep = manifest.firstStep; lastStep = manifest.lastStep
        firstTime = manifest.firstTimePS; lastTime = manifest.lastTimePS
    }

    public static func resume(store: VivoArtifactStore, manifest fingerprint: VivoFingerprint,
                              targetChunkBytes: Int = 8 * 1024 * 1024) async throws -> VivoMDTrajectoryArchiveWriter {
        let reader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: fingerprint)
        guard !reader.manifest.sealed else { throw bad("sealed trajectory cannot be extended; create a new segment") }
        let validated = try await reader.validate(scope: .restart)
        return try resume(validated: validated, targetChunkBytes: targetChunkBytes)
    }

    /// Immediate in-module handoff of an already verified prefix. Its authority
    /// is the captured store actor and manifest, never caller-substituted fields.
    /// This avoids a second full traversal in protocol resume without exporting
    /// a public API that claims old evidence is fresh filesystem verification.
    static func resume(validated: VivoMDTrajectoryValidation,
                       targetChunkBytes: Int = 8 * 1024 * 1024) throws -> VivoMDTrajectoryArchiveWriter {
        try Task.checkCancellation()
        let reader = validated.reader
        guard validated.scope == .restart || validated.scope == .allPayloads else {
            throw bad("trajectory resume requires verification of its final payload")
        }
        guard !reader.manifest.sealed else { throw bad("sealed trajectory cannot be extended; create a new segment") }
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: reader.manifest.particleCount,
                                                                includeVelocities: reader.manifest.includeVelocities)
        guard targetChunkBytes > VivoMDTrajectoryChunkCodec.headerBytes,
              targetChunkBytes <= VivoMDTrajectoryChunkCodec.maximumChunkBytes else { throw bad("invalid trajectory chunk budget") }
        let capacity = min(VivoMDTrajectoryChunkCodec.maximumFramesPerChunk,
                           max(1, (targetChunkBytes - VivoMDTrajectoryChunkCodec.headerBytes) / stride))
        try Task.checkCancellation()
        return .init(store: reader.store, manifest: reader.manifest, framesPerChunk: capacity)
    }

    public func append(_ snapshot: VivoMDStateSnapshot) async throws {
        try reserve(); defer { busy = false }
        guard snapshot.systemFingerprint == system, snapshot.configurationFingerprint == configuration else {
            throw bad("trajectory sample model/configuration mismatch")
        }
        let frame = VivoMDArchiveFrame(snapshot: snapshot, includeVelocities: velocities)
        try VivoMDTrajectoryChunkCodec.validate(frame, particleCount: particles, includeVelocities: velocities)
        guard persistedFrames <= UInt64.max - UInt64(buffered.count) - 1 else { throw bad("trajectory frame index exhausted") }
        if let lastStep, let lastTime {
            guard frame.stepIndex > lastStep, frame.timePS > lastTime else { throw bad("trajectory sample is not a newer accepted state") }
        }
        if buffered.count == framesPerChunk { try await flushReserved() }
        buffered.append(frame)
        if firstStep == nil { firstStep = frame.stepIndex; firstTime = frame.timePS }
        lastStep = frame.stepIndex; lastTime = frame.timePS
    }

    /// Flushes all buffered frames before producing a restartable prefix. Prior
    /// manifests remain valid if a later append or storage operation fails.
    public func snapshot() async throws -> VivoStoredArtifact {
        try reserve(); defer { busy = false }
        return try await publish(seal: false)
    }
    public func finish() async throws -> VivoStoredArtifact {
        try reserve(); defer { busy = false }
        let result = try await publish(seal: true)
        sealed = true
        return result
    }
    private func publish(seal: Bool) async throws -> VivoStoredArtifact {
        try await flushReserved()
        let manifest = VivoMDTrajectoryManifest(schema: VivoMDTrajectoryManifest.schemaID,
            systemFingerprint: system, configurationFingerprint: configuration, particleCount: particles,
            includeVelocities: velocities, frameCount: persistedFrames, chunkCount: chunkCount, tail: tail,
            firstStep: firstStep, lastStep: lastStep, firstTimePS: firstTime, lastTimePS: lastTime, sealed: seal)
        try manifest.validate()
        return try await store.put(data: VivoCanonicalJSON.encode(manifest), kind: "md-trajectory-manifest",
                                    mediaType: "application/vnd.numivivo.md-trajectory-manifest+json")
    }
    private func flushReserved() async throws {
        guard let first = buffered.first, let last = buffered.last else { return }
        guard chunkCount < UInt64.max else { throw bad("trajectory chunk index exhausted") }
        let data = try VivoMDTrajectoryChunkCodec.encode(buffered, systemFingerprint: system,
            configurationFingerprint: configuration, particleCount: particles, includeVelocities: velocities,
            firstFrameOrdinal: persistedFrames)
        let payload = try await store.put(data: data, kind: "md-trajectory-chunk",
                                          mediaType: VivoMDTrajectoryChunkCodec.mediaType)
        let link = VivoMDTrajectoryChunkLink(schema: VivoMDTrajectoryChunkLink.schemaID,
            systemFingerprint: system, configurationFingerprint: configuration, particleCount: particles,
            includeVelocities: velocities, chunkOrdinal: chunkCount, firstFrameOrdinal: persistedFrames,
            frameCount: UInt32(buffered.count), previous: tail, payload: payload.fingerprint,
            payloadBytes: payload.byteCount, firstStep: first.stepIndex, lastStep: last.stepIndex,
            firstTimePS: first.timePS, lastTimePS: last.timePS)
        let storedLink = try await store.put(data: VivoCanonicalJSON.encode(link), kind: "md-trajectory-link",
                                             mediaType: "application/vnd.numivivo.md-trajectory-link+json")
        // Publication occurs only after both immutable objects were persisted.
        persistedFrames += UInt64(buffered.count); chunkCount += 1; tail = storedLink.fingerprint
        buffered.removeAll(keepingCapacity: true)
    }
    private func reserve() throws {
        guard !busy, !sealed else { throw bad("trajectory writer is busy or sealed") }
        busy = true
    }
}

private func bad(_ message: String) -> VivoArtifactValidationError { .invalid(message) }
