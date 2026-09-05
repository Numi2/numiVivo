import Foundation
@preconcurrency import Metal

public struct VivoMolecularCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-checkpoint/v2"

    public let schema: String
    public let programFingerprint: VivoFingerprint
    public let sourceProgramFingerprint: VivoFingerprint
    public let fidelity: VivoFidelity
    public let seed: VivoRuntimeSeed
    public let stepIndex: UInt32
    public let absoluteTimeSeconds: Double
    public let speciesCount: UInt32
    public let laneCount: UInt32
    public let parameterCount: UInt32
    public let parameterEnvironmentCount: UInt32
    public let temporalStateCount: UInt32
    public let stateFP32LE: Data
    public let parametersFP32LE: Data
    public let temporalStateFP32LE: Data
    public let transportRecordLE: Data
    public let velocityFP32LE: Data
    public let volumeFractionFP32LE: Data

    public init(
        programFingerprint: VivoFingerprint,
        sourceProgramFingerprint: VivoFingerprint,
        fidelity: VivoFidelity,
        seed: VivoRuntimeSeed,
        stepIndex: UInt32,
        absoluteTimeSeconds: Double,
        speciesCount: UInt32,
        laneCount: UInt32,
        parameterCount: UInt32,
        parameterEnvironmentCount: UInt32,
        temporalStateCount: UInt32,
        stateFP32LE: Data,
        parametersFP32LE: Data,
        temporalStateFP32LE: Data,
        transportRecordLE: Data,
        velocityFP32LE: Data,
        volumeFractionFP32LE: Data
    ) {
        self.schema = Self.schema
        self.programFingerprint = programFingerprint
        self.sourceProgramFingerprint = sourceProgramFingerprint
        self.fidelity = fidelity
        self.seed = seed
        self.stepIndex = stepIndex
        self.absoluteTimeSeconds = absoluteTimeSeconds
        self.speciesCount = speciesCount
        self.laneCount = laneCount
        self.parameterCount = parameterCount
        self.parameterEnvironmentCount = parameterEnvironmentCount
        self.temporalStateCount = temporalStateCount
        self.stateFP32LE = stateFP32LE
        self.parametersFP32LE = parametersFP32LE
        self.temporalStateFP32LE = temporalStateFP32LE
        self.transportRecordLE = transportRecordLE
        self.velocityFP32LE = velocityFP32LE
        self.volumeFractionFP32LE = volumeFractionFP32LE
    }

    public func validate() throws {
        guard schema == Self.schema,
              absoluteTimeSeconds.isFinite,
              absoluteTimeSeconds >= 0,
              speciesCount > 0,
              laneCount > 0,
              parameterEnvironmentCount > 0 else {
            throw VivoArtifactValidationError.invalid(
                "molecular checkpoint header is invalid"
            )
        }

        try requireFP32(
            stateFP32LE,
            count: try product(UInt64(speciesCount), UInt64(laneCount), label: "state"),
            label: "state"
        )
        try requireFP32(
            parametersFP32LE,
            count: try product(UInt64(parameterCount), UInt64(parameterEnvironmentCount), label: "parameters"),
            label: "parameters"
        )
        let temporalScalars = try product(
            try product(UInt64(temporalStateCount), UInt64(laneCount), label: "temporal state"),
            2,
            label: "temporal float2 state"
        )
        try requireFP32(
            temporalStateFP32LE,
            count: max(temporalScalars, 2),
            label: "temporal state"
        )
        let expectedTransport = try product(
            UInt64(max(speciesCount, 1)),
            UInt64(VivoTransportRecordLE.stride),
            label: "transport"
        )
        guard expectedTransport <= UInt64(Int.max),
              transportRecordLE.count == Int(expectedTransport) else {
            throw VivoArtifactValidationError.invalid(
                "molecular checkpoint transport byte count is invalid"
            )
        }
        try requireFP32(
            velocityFP32LE,
            count: try product(UInt64(laneCount), 4, label: "velocity"),
            label: "velocity"
        )
        try requireFP32(
            volumeFractionFP32LE,
            count: UInt64(laneCount),
            label: "volume fractions"
        )

        let state = try VivoLittleEndianFP32.decode(stateFP32LE)
        let parameters = try VivoLittleEndianFP32.decode(parametersFP32LE)
        let temporal = try VivoLittleEndianFP32.decode(temporalStateFP32LE)
        let velocity = try VivoLittleEndianFP32.decode(velocityFP32LE)
        let volume = try VivoLittleEndianFP32.decode(volumeFractionFP32LE)
        guard state.allSatisfy(\.isFinite),
              parameters.allSatisfy(\.isFinite),
              temporal.allSatisfy(\.isFinite),
              velocity.allSatisfy(\.isFinite),
              volume.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }) else {
            throw VivoArtifactValidationError.invalid(
                "molecular checkpoint contains non-finite or invalid scalar data"
            )
        }
        _ = try VivoTransportRecordLE.decode(transportRecordLE)
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    private func requireFP32(
        _ data: Data,
        count: UInt64,
        label: String
    ) throws {
        let bytes = try product(count, 4, label: "\(label) bytes")
        guard bytes <= UInt64(Int.max), data.count == Int(bytes) else {
            throw VivoArtifactValidationError.invalid(
                "molecular checkpoint \(label) byte count does not match its shape"
            )
        }
    }

    private func product(
        _ lhs: UInt64,
        _ rhs: UInt64,
        label: String
    ) throws -> UInt64 {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else {
            throw VivoArtifactValidationError.invalid(
                "molecular checkpoint \(label) size overflow"
            )
        }
        return result.partialValue
    }
}

struct VivoMolecularArenaCheckpointData: Sendable {
    let state: [Float]
    let parameters: [Float]
    let temporalState: [Float]
    let transport: [VivoSpeciesTransportABI]
    let velocity: [Float]
    let volumeFractions: [Float]
}

extension VivoMetalArena {
    func captureCheckpoint(
        commandQueue: MTLCommandQueue
    ) throws -> VivoMolecularArenaCheckpointData {
        let sources: [(MTLBuffer, Int)] = [
            (currentState, currentState.length),
            (parameterBuffer, parameterBuffer.length),
            (temporalCurrent, temporalCurrent.length),
            (transport, transport.length),
            (velocity, velocity.length),
            (volumeFraction, volumeFraction.length)
        ]
        let total = try sources.reduce(into: 0) { partial, item in
            let addition = partial.addingReportingOverflow(item.1)
            guard !addition.overflow else {
                throw VivoRuntimeError.allocationFailed(
                    "molecular checkpoint staging size overflow"
                )
            }
            partial = addition.partialValue
        }
        guard let staging = device.makeBuffer(
            length: max(total, 1),
            options: .storageModeShared
        ),
        let command = commandQueue.makeCommandBuffer(),
        let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.allocationFailed(
                "could not allocate molecular checkpoint readback"
            )
        }
        staging.label = "NumiVivo.Checkpoint.Readback"
        command.label = "NumiVivo.Checkpoint.Capture"
        var offset = 0
        for (source, length) in sources {
            blit.copy(
                from: source,
                sourceOffset: 0,
                to: staging,
                destinationOffset: offset,
                size: length
            )
            offset += length
        }
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error {
            throw VivoRuntimeError.commandFailed(
                "molecular checkpoint capture: \(error)"
            )
        }

        var cursor = 0
        func floats(length: Int) -> [Float] {
            defer { cursor += length }
            let pointer = staging.contents()
                .advanced(by: cursor)
                .assumingMemoryBound(to: Float.self)
            return Array(
                UnsafeBufferPointer(start: pointer, count: length / 4)
            )
        }
        let state = floats(length: currentState.length)
        let parameters = floats(length: parameterBuffer.length)
        let temporalState = floats(length: temporalCurrent.length)

        let transportLength = transport.length
        let transportPointer = staging.contents()
            .advanced(by: cursor)
            .assumingMemoryBound(to: VivoSpeciesTransportABI.self)
        let transportRecords = Array(
            UnsafeBufferPointer(
                start: transportPointer,
                count: transportLength / MemoryLayout<VivoSpeciesTransportABI>.stride
            )
        )
        cursor += transportLength
        let velocity = floats(length: self.velocity.length)
        let volumeFractions = floats(length: volumeFraction.length)

        return VivoMolecularArenaCheckpointData(
            state: state,
            parameters: parameters,
            temporalState: temporalState,
            transport: transportRecords,
            velocity: velocity,
            volumeFractions: volumeFractions
        )
    }

    func restoreCheckpoint(
        _ data: VivoMolecularArenaCheckpointData,
        commandQueue: MTLCommandQueue
    ) throws {
        guard data.state.count * 4 == currentState.length,
              data.parameters.count * 4 == parameterBuffer.length,
              data.temporalState.count * 4 == temporalCurrent.length,
              data.transport.count * MemoryLayout<VivoSpeciesTransportABI>.stride == transport.length,
              data.velocity.count * 4 == velocity.length,
              data.volumeFractions.count * 4 == volumeFraction.length,
              data.state.allSatisfy(\.isFinite),
              data.parameters.allSatisfy(\.isFinite),
              data.temporalState.allSatisfy(\.isFinite),
              data.velocity.allSatisfy(\.isFinite),
              data.volumeFractions.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }) else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular checkpoint does not match the allocated runtime arena"
            )
        }

        let stateStage = try checkpointStaging(data.state, label: "State")
        let parameterStage = try checkpointStaging(data.parameters, label: "Parameters")
        let temporalStage = try checkpointStaging(data.temporalState, label: "Temporal")
        let transportStage = try checkpointStaging(data.transport, label: "Transport")
        let velocityStage = try checkpointStaging(data.velocity, label: "Velocity")
        let volumeStage = try checkpointStaging(data.volumeFractions, label: "Volume")

        guard let command = commandQueue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Checkpoint.Restore"
        for destination in [currentState, baseState, stageState, candidateState] {
            blit.copy(
                from: stateStage,
                sourceOffset: 0,
                to: destination,
                destinationOffset: 0,
                size: currentState.length
            )
        }
        blit.fill(buffer: derivativeK1, range: 0..<derivativeK1.length, value: 0)
        blit.copy(from: parameterStage, sourceOffset: 0, to: parameterBuffer, destinationOffset: 0, size: parameterBuffer.length)
        blit.copy(from: temporalStage, sourceOffset: 0, to: temporalCurrent, destinationOffset: 0, size: temporalCurrent.length)
        blit.copy(from: temporalStage, sourceOffset: 0, to: temporalCandidate, destinationOffset: 0, size: temporalCandidate.length)
        blit.copy(from: transportStage, sourceOffset: 0, to: transport, destinationOffset: 0, size: transport.length)
        blit.copy(from: velocityStage, sourceOffset: 0, to: velocity, destinationOffset: 0, size: velocity.length)
        blit.copy(from: volumeStage, sourceOffset: 0, to: volumeFraction, destinationOffset: 0, size: volumeFraction.length)
        blit.fill(buffer: reactionEvents, range: 0..<reactionEvents.length, value: 0)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        if let error = command.error {
            throw VivoRuntimeError.commandFailed(
                "molecular checkpoint restore: \(error)"
            )
        }
    }

    private func checkpointStaging<T>(
        _ values: [T],
        label: String
    ) throws -> MTLBuffer {
        let byteCount = values.count * MemoryLayout<T>.stride
        let buffer: MTLBuffer? = values.withUnsafeBytes { bytes in
            guard let base = bytes.baseAddress else { return nil }
            return device.makeBuffer(
                bytes: base,
                length: max(byteCount, 1),
                options: .storageModeShared
            )
        }
        guard let buffer else {
            throw VivoRuntimeError.allocationFailed(
                "could not allocate checkpoint staging buffer \(label)"
            )
        }
        buffer.label = "NumiVivo.Checkpoint.Stage.\(label)"
        return buffer
    }
}
