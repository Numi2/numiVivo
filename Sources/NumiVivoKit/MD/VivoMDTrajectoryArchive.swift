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

/// One bounded chunk of decoded frames at a time. The persistent index is a
/// content-addressed linked sequence, so every flush writes O(1) metadata rather
/// than rewriting a growing JSON array. A manifest is an immutable prefix view.
public struct VivoMDTrajectoryArchiveReader: Sendable {
    public let store: VivoArtifactStore
    public let manifest: VivoMDTrajectoryManifest
    public let manifestFingerprint: VivoFingerprint

    public static func open(store: VivoArtifactStore, manifest fingerprint: VivoFingerprint) async throws -> Self {
        let data = try await store.data(for: fingerprint)
        guard data.count <= 64 * 1024 else { throw bad("oversized trajectory manifest") }
        let value = try VivoCanonicalJSON.decode(VivoMDTrajectoryManifest.self, from: data)
        try value.validate()
        return .init(store: store, manifest: value, manifestFingerprint: fingerprint)
    }

    /// Returns chronological link identities, without loading coordinate data.
    /// Bound before allocation/traversal. Cycles, missing links, skipped ordinals,
    /// reordered frames and mismatched model/configuration identities are errors.
    public func index(maximumChunks: UInt64 = 100_000) async throws -> [VivoFingerprint] {
        guard manifest.chunkCount <= maximumChunks, manifest.chunkCount <= UInt64(Int.max) else {
            throw bad("trajectory index exceeds caller's chunk limit")
        }
        var result: [VivoFingerprint] = []
        var visited = Set<VivoFingerprint>()
        var next = manifest.tail
        var chunksRemaining = manifest.chunkCount
        var endOrdinal = manifest.frameCount
        var newer: VivoMDTrajectoryChunkLink?
        while let fingerprint = next {
            guard chunksRemaining > 0, visited.insert(fingerprint).inserted else { throw bad("trajectory link cycle or excess links") }
            let link = try await readLink(fingerprint)
            guard link.chunkOrdinal == chunksRemaining - 1,
                  UInt64(link.frameCount) <= endOrdinal,
                  link.firstFrameOrdinal == endOrdinal - UInt64(link.frameCount) else {
                throw bad("trajectory link ordinal gap")
            }
            if let newer {
                guard link.lastStep < newer.firstStep, link.lastTimePS < newer.firstTimePS else { throw bad("overlapping trajectory chunk ranges") }
            } else {
                guard link.lastStep == manifest.lastStep, link.lastTimePS == manifest.lastTimePS else { throw bad("trajectory tail range disagrees with manifest") }
            }
            result.append(fingerprint)
            next = link.previous; chunksRemaining -= 1; endOrdinal = link.firstFrameOrdinal; newer = link
        }
        guard chunksRemaining == 0, endOrdinal == 0 else { throw bad("trajectory prefix is truncated") }
        if let first = newer {
            guard first.firstStep == manifest.firstStep, first.firstTimePS == manifest.firstTimePS else { throw bad("trajectory first range disagrees with manifest") }
        }
        return result.reversed()
    }

    public func readLink(_ fingerprint: VivoFingerprint) async throws -> VivoMDTrajectoryChunkLink {
        let data = try await store.data(for: fingerprint)
        guard data.count <= 64 * 1024 else { throw bad("oversized trajectory link") }
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
        return value
    }

    public func readChunk(_ linkFingerprint: VivoFingerprint) async throws -> [VivoMDArchiveFrame] {
        let link = try await readLink(linkFingerprint)
        let data = try await store.data(for: link.payload)
        guard UInt64(data.count) == link.payloadBytes else { throw bad("trajectory payload byte count mismatch") }
        let frames = try VivoMDTrajectoryChunkCodec.decode(data,
            systemFingerprint: manifest.systemFingerprint, configurationFingerprint: manifest.configurationFingerprint,
            particleCount: manifest.particleCount, includeVelocities: manifest.includeVelocities,
            firstFrameOrdinal: link.firstFrameOrdinal)
        guard frames.count == Int(link.frameCount), frames.first?.stepIndex == link.firstStep,
              frames.last?.stepIndex == link.lastStep, frames.first?.timePS == link.firstTimePS,
              frames.last?.timePS == link.lastTimePS else { throw bad("trajectory payload range disagrees with index") }
        return frames
    }

    public func verify(maximumChunks: UInt64 = 100_000) async throws {
        for fingerprint in try await index(maximumChunks: maximumChunks) {
            try Task.checkCancellation()
            _ = try await readChunk(fingerprint)
        }
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
        _ = try await reader.index()
        if let tail = reader.manifest.tail { _ = try await reader.readChunk(tail) }
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: reader.manifest.particleCount,
                                                                includeVelocities: reader.manifest.includeVelocities)
        guard targetChunkBytes > VivoMDTrajectoryChunkCodec.headerBytes,
              targetChunkBytes <= VivoMDTrajectoryChunkCodec.maximumChunkBytes else { throw bad("invalid trajectory chunk budget") }
        let capacity = min(VivoMDTrajectoryChunkCodec.maximumFramesPerChunk,
                           max(1, (targetChunkBytes - VivoMDTrajectoryChunkCodec.headerBytes) / stride))
        return .init(store: store, manifest: reader.manifest, framesPerChunk: capacity)
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
