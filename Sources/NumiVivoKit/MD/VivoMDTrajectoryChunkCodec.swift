import Foundation

/// Storage frames are not checkpoints: velocity may be omitted, but coordinate
/// quantization is never implicit. FP32 coordinates round-trip bit-for-bit.
public struct VivoMDArchiveFrame: Codable, Sendable, Equatable {
    public let stepIndex: UInt64
    public let timePS: Double
    public let positionsNM: [VivoVector3D]
    public let velocitiesNMPerPS: [VivoVector3D]?
    public let periodicCell: VivoPeriodicCell?

    public init(snapshot: VivoMDStateSnapshot, includeVelocities: Bool) {
        stepIndex = snapshot.stepIndex
        timePS = snapshot.timePS
        positionsNM = snapshot.positionsNM
        velocitiesNMPerPS = includeVelocities ? snapshot.velocitiesNMPerPS : nil
        periodicCell = snapshot.periodicCell
    }

    public init(stepIndex: UInt64, timePS: Double, positionsNM: [VivoVector3D],
                velocitiesNMPerPS: [VivoVector3D]?, periodicCell: VivoPeriodicCell?) {
        self.stepIndex = stepIndex; self.timePS = timePS
        self.positionsNM = positionsNM; self.velocitiesNMPerPS = velocitiesNMPerPS
        self.periodicCell = periodicCell
    }
}

/// Portable little-endian trajectory chunks. Header = 160 bytes; frame metadata
/// = 96 bytes; positions = 12 bytes/particle; optional velocities add 12 bytes.
/// Cells and times retain FP64. Integrity is provided by VivoArtifactStore's
/// SHA-256 object identity rather than a second, incompatible persistence layer.
public enum VivoMDTrajectoryChunkCodec {
    public static let mediaType = "application/vnd.numivivo.md-trajectory-chunk"
    public static let maximumChunkBytes = 64 * 1024 * 1024
    public static let maximumFramesPerChunk = 4096
    public static let headerBytes = 160
    private static let magic = Data("NVMDTRJ1".utf8)

    public static func frameBytes(particleCount: UInt32, includeVelocities: Bool) throws -> Int {
        guard particleCount > 0 else { throw invalid("empty trajectory particle set") }
        let bytes = UInt64(particleCount) * (includeVelocities ? 24 : 12) + 96
        guard bytes <= UInt64(maximumChunkBytes - headerBytes) else {
            throw invalid("one trajectory frame exceeds the 64 MiB chunk limit")
        }
        return Int(bytes)
    }

    public static func validate(_ frame: VivoMDArchiveFrame, particleCount: UInt32,
                                includeVelocities: Bool) throws {
        guard frame.positionsNM.count == Int(particleCount), frame.timePS.isFinite,
              frame.timePS >= 0, frame.periodicCell?.isValid != false,
              (frame.velocitiesNMPerPS != nil) == includeVelocities,
              frame.velocitiesNMPerPS.map({ $0.count == Int(particleCount) }) ?? true else {
            throw invalid("trajectory frame shape, cell, velocity policy or time mismatch")
        }
        for value in frame.positionsNM { try validateFP32(value) }
        for value in frame.velocitiesNMPerPS ?? [] { try validateFP32(value) }
    }

    public static func encode(_ frames: [VivoMDArchiveFrame], systemFingerprint: VivoFingerprint,
                              configurationFingerprint: VivoFingerprint, particleCount: UInt32,
                              includeVelocities: Bool, firstFrameOrdinal: UInt64) throws -> Data {
        let stride = try frameBytes(particleCount: particleCount, includeVelocities: includeVelocities)
        guard !frames.isEmpty, frames.count <= maximumFramesPerChunk,
              frames.count <= (maximumChunkBytes - headerBytes) / stride,
              UInt64(frames.count) <= UInt64.max - firstFrameOrdinal else {
            throw invalid("trajectory chunk frame count or ordinal exceeds bounds")
        }
        var data = Data(capacity: headerBytes + frames.count * stride)
        data.append(magic)
        append(UInt32(1), to: &data)
        append(UInt32(includeVelocities ? 1 : 0), to: &data)
        append(particleCount, to: &data)
        append(UInt32(frames.count), to: &data)
        append(firstFrameOrdinal, to: &data)
        data.append(contentsOf: systemFingerprint.hex.utf8)
        data.append(contentsOf: configurationFingerprint.hex.utf8)
        guard data.count == headerBytes else { throw invalid("fingerprint header width mismatch") }
        var previous: VivoMDArchiveFrame?
        for frame in frames {
            try validate(frame, particleCount: particleCount, includeVelocities: includeVelocities)
            if let previous {
                guard frame.stepIndex > previous.stepIndex, frame.timePS > previous.timePS else {
                    throw invalid("trajectory frames must increase in accepted step and physical time")
                }
            }
            append(frame.stepIndex, to: &data)
            append(frame.timePS.bitPattern, to: &data)
            append(UInt32(frame.periodicCell == nil ? 0 : 1), to: &data)
            append(UInt32(0), to: &data)
            let cell = frame.periodicCell
            for vector in [cell?.a ?? .zero, cell?.b ?? .zero, cell?.c ?? .zero] {
                for component in [vector.x, vector.y, vector.z] { append(component.bitPattern, to: &data) }
            }
            for vector in frame.positionsNM { append(vector, to: &data) }
            for vector in frame.velocitiesNMPerPS ?? [] { append(vector, to: &data) }
            previous = frame
        }
        return data
    }

    public static func decode(_ data: Data, systemFingerprint: VivoFingerprint,
                              configurationFingerprint: VivoFingerprint, particleCount: UInt32,
                              includeVelocities: Bool, firstFrameOrdinal: UInt64) throws -> [VivoMDArchiveFrame] {
        guard data.count >= headerBytes, data.count <= maximumChunkBytes,
              data.prefix(8) == magic else { throw invalid("invalid trajectory magic or chunk length") }
        var reader = Reader(data: data, offset: 8)
        let version: UInt32 = try reader.word()
        let flags: UInt32 = try reader.word()
        let particles: UInt32 = try reader.word()
        let frameCount: UInt32 = try reader.word()
        let ordinal: UInt64 = try reader.word()
        let system = try reader.bytes(64)
        let configuration = try reader.bytes(64)
        guard version == 1, flags == (includeVelocities ? 1 : 0), particles == particleCount,
              ordinal == firstFrameOrdinal,
              system == Data(systemFingerprint.hex.utf8), configuration == Data(configurationFingerprint.hex.utf8),
              frameCount > 0, frameCount <= UInt32(maximumFramesPerChunk),
              UInt64(frameCount) <= UInt64.max - ordinal else {
            throw invalid("trajectory header identity, policy, version or shape mismatch")
        }
        let stride = try frameBytes(particleCount: particles, includeVelocities: includeVelocities)
        guard UInt64(data.count) == UInt64(headerBytes) + UInt64(frameCount) * UInt64(stride) else {
            throw invalid("trajectory payload is truncated or has trailing data")
        }
        var frames: [VivoMDArchiveFrame] = []
        frames.reserveCapacity(Int(frameCount))
        for _ in 0..<frameCount {
            let step: UInt64 = try reader.word()
            let timeBits: UInt64 = try reader.word()
            let cellFlag: UInt32 = try reader.word()
            let reserved: UInt32 = try reader.word()
            guard cellFlag <= 1, reserved == 0 else { throw invalid("unknown trajectory frame flags") }
            var vectors: [VivoVector3D] = []
            for _ in 0..<3 {
                let x: UInt64 = try reader.word(), y: UInt64 = try reader.word(), z: UInt64 = try reader.word()
                vectors.append(.init(Double(bitPattern: x), Double(bitPattern: y), Double(bitPattern: z)))
            }
            if cellFlag == 0, !vectors.allSatisfy({ $0 == .zero }) { throw invalid("absent cell has nonzero payload") }
            let cell: VivoPeriodicCell? = cellFlag == 0 ? nil : .init(a: vectors[0], b: vectors[1], c: vectors[2])
            var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
            positions.reserveCapacity(Int(particles))
            if includeVelocities { velocities.reserveCapacity(Int(particles)) }
            for _ in 0..<particles { positions.append(try reader.vector()) }
            if includeVelocities { for _ in 0..<particles { velocities.append(try reader.vector()) } }
            let frame = VivoMDArchiveFrame(stepIndex: step, timePS: Double(bitPattern: timeBits),
                                            positionsNM: positions, velocitiesNMPerPS: includeVelocities ? velocities : nil,
                                            periodicCell: cell)
            try validate(frame, particleCount: particles, includeVelocities: includeVelocities)
            if let last = frames.last, frame.stepIndex <= last.stepIndex || frame.timePS <= last.timePS {
                throw invalid("nonmonotone trajectory payload")
            }
            frames.append(frame)
        }
        return frames
    }

    private static func validateFP32(_ value: VivoVector3D) throws {
        guard value.isFinite else { throw invalid("nonfinite trajectory coordinate") }
        for x in [value.x, value.y, value.z] {
            guard Double(Float(x)) == x else {
                throw invalid("trajectory codec accepts exact FP32 runtime values; explicit conversion is required for higher precision input")
            }
        }
    }
    private static func append<T: FixedWidthInteger>(_ value: T, to data: inout Data) {
        var little = value.littleEndian
        withUnsafeBytes(of: &little) { data.append(contentsOf: $0) }
    }
    private static func append(_ vector: VivoVector3D, to data: inout Data) {
        for value in [vector.x, vector.y, vector.z] { append(Float(value).bitPattern, to: &data) }
    }
    private static func invalid(_ message: String) -> VivoArtifactValidationError { .invalid(message) }
    private struct Reader {
        let data: Data
        var offset: Int
        mutating func word<T: FixedWidthInteger>() throws -> T {
            let width = MemoryLayout<T>.size
            guard offset <= data.count, width <= data.count - offset else { throw invalid("truncated trajectory scalar") }
            let value: T = data.withUnsafeBytes { raw in
                T(littleEndian: raw.loadUnaligned(fromByteOffset: offset, as: T.self))
            }
            offset += width
            return value
        }
        mutating func bytes(_ count: Int) throws -> Data {
            guard count >= 0, offset <= data.count, count <= data.count - offset else { throw invalid("truncated trajectory bytes") }
            let result = data.subdata(in: offset..<(offset + count))
            offset += count
            return result
        }
        mutating func vector() throws -> VivoVector3D {
            let x: UInt32 = try word(), y: UInt32 = try word(), z: UInt32 = try word()
            return .init(Double(Float(bitPattern: x)), Double(Float(bitPattern: y)), Double(Float(bitPattern: z)))
        }
    }
}
