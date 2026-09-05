import Foundation
@preconcurrency import Metal

public struct VivoNumiVivoChannelTransform: Codable, Sendable, Equatable {
    public var scale: Float
    public var offset: Float
    public var minimum: Float?
    public var maximum: Float?

    public init(
        scale: Float = 1,
        offset: Float = 0,
        minimum: Float? = nil,
        maximum: Float? = nil
    ) {
        self.scale = scale
        self.offset = offset
        self.minimum = minimum
        self.maximum = maximum
    }

    public func validate(subject: String) throws {
        guard scale.isFinite,
              offset.isFinite,
              minimum?.isFinite ?? true,
              maximum?.isFinite ?? true,
              minimum == nil || maximum == nil || minimum! <= maximum! else {
            throw VivoNumiLabCouplingError.invalidChannel(
                "binding transform for \(subject) is invalid"
            )
        }
    }

    public func apply(_ value: Float, subject: String) throws -> Float {
        guard value.isFinite else {
            throw VivoNumiLabCouplingError.invalidOutput(
                "binding \(subject) received a non-finite value"
            )
        }
        var transformed = value * scale + offset
        guard transformed.isFinite else {
            throw VivoNumiLabCouplingError.invalidOutput(
                "binding \(subject) produced a non-finite transformed value"
            )
        }
        if let minimum { transformed = max(transformed, minimum) }
        if let maximum { transformed = min(transformed, maximum) }
        return transformed
    }
}

public struct VivoNumiVivoInputBinding: Codable, Sendable, Equatable {
    public let channelID: String
    public let unit: String
    public let speciesIdentifier: String
    public let laneIndices: [UInt32]
    public let mode: VivoCouplingMode
    public let transform: VivoNumiVivoChannelTransform

    public init(
        channelID: String,
        unit: String,
        speciesIdentifier: String,
        laneIndices: [UInt32],
        mode: VivoCouplingMode = .replace,
        transform: VivoNumiVivoChannelTransform = .init()
    ) {
        self.channelID = channelID
        self.unit = unit
        self.speciesIdentifier = speciesIdentifier
        self.laneIndices = laneIndices
        self.mode = mode
        self.transform = transform
    }
}

public struct VivoNumiVivoOutputBinding: Codable, Sendable, Equatable {
    public let channelID: String
    public let unit: String
    public let speciesIdentifier: String
    public let laneIndices: [UInt32]
    public let transform: VivoNumiVivoChannelTransform

    public init(
        channelID: String,
        unit: String,
        speciesIdentifier: String,
        laneIndices: [UInt32],
        transform: VivoNumiVivoChannelTransform = .init()
    ) {
        self.channelID = channelID
        self.unit = unit
        self.speciesIdentifier = speciesIdentifier
        self.laneIndices = laneIndices
        self.transform = transform
    }
}

public struct VivoNumiVivoParticipantConfiguration: Codable, Sendable, Equatable {
    public var inputs: [VivoNumiVivoInputBinding]
    public var outputs: [VivoNumiVivoOutputBinding]

    public init(
        inputs: [VivoNumiVivoInputBinding],
        outputs: [VivoNumiVivoOutputBinding]
    ) {
        self.inputs = inputs
        self.outputs = outputs
    }
}

public actor VivoNumiVivoCoupledParticipant: VivoCoupledParticipant {
    private enum Phase {
        case prepared
        case staged
    }

    private struct PreparedInputBinding: Sendable {
        let source: VivoNumiVivoInputBinding
        let speciesIndex: UInt32
    }

    private struct PreparedOutputBinding: Sendable {
        let source: VivoNumiVivoOutputBinding
        let speciesIndex: UInt32
        let publicationOffset: Int
    }

    private struct Pending: Sendable {
        let transactionID: String
        let runtimeTransactionID: UUID
        let prepared: VivoPreparedMolecularStep
        var phase: Phase
    }

    private struct StateDigestMaterial: Codable {
        let transactionID: String
        let programFingerprint: VivoFingerprint
        let stepIndex: UInt32
        let candidateTimeStep: Float
        let runtimeStatus: VivoRuntimeStatus
        let outputs: [VivoCouplingSample]
    }

    public nonisolated let participantID: VivoNumiLabParticipant = .numiVivo
    public nonisolated let programFingerprint: VivoFingerprint
    public nonisolated let sourceProgramFingerprint: VivoFingerprint
    public nonisolated let capabilities: VivoMetalCapabilities

    private let runtime: VivoTransactionalMolecularRuntime
    private let minimumTimeStep: Float
    private let inputBindings: [PreparedInputBinding]
    private let outputBindings: [PreparedOutputBinding]
    private var pending: Pending?
    private var releaseFailure: String?

    public static func make(
        pack: VivoProgramPack,
        runtimeConfiguration: VivoRuntimeConfiguration,
        participantConfiguration: VivoNumiVivoParticipantConfiguration,
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoNumiVivoCoupledParticipant {
        let metadata = try pack.speciesMetadata()
        let byIdentifier = Dictionary(
            uniqueKeysWithValues: metadata.map { ($0.identifier, $0) }
        )
        var inputDestinations = Set<UInt64>()
        var preparedInputs: [PreparedInputBinding] = []
        var preparedOutputs: [PreparedOutputBinding] = []
        var outputOffset = 0

        for binding in participantConfiguration.inputs {
            try binding.transform.validate(subject: binding.channelID)
            guard !binding.channelID.isEmpty,
                  !binding.unit.isEmpty,
                  !binding.laneIndices.isEmpty,
                  Set(binding.laneIndices).count == binding.laneIndices.count,
                  binding.laneIndices.allSatisfy({ $0 < runtimeConfiguration.laneCount }),
                  let species = byIdentifier[binding.speciesIdentifier],
                  species.isExternallyOwned || species.isInput else {
                throw VivoNumiLabCouplingError.invalidChannel(
                    "input binding \(binding.channelID) is unresolved, duplicated, out of range, or targets internally owned state"
                )
            }
            for lane in binding.laneIndices {
                let key = UInt64(species.index) << 32 | UInt64(lane)
                guard inputDestinations.insert(key).inserted else {
                    throw VivoNumiLabCouplingError.invalidChannel(
                        "multiple input bindings write species \(species.identifier) lane \(lane)"
                    )
                }
            }
            preparedInputs.append(.init(source: binding, speciesIndex: species.index))
        }

        for binding in participantConfiguration.outputs {
            try binding.transform.validate(subject: binding.channelID)
            guard !binding.channelID.isEmpty,
                  !binding.unit.isEmpty,
                  !binding.laneIndices.isEmpty,
                  Set(binding.laneIndices).count == binding.laneIndices.count,
                  binding.laneIndices.allSatisfy({ $0 < runtimeConfiguration.laneCount }),
                  let species = byIdentifier[binding.speciesIdentifier] else {
                throw VivoNumiLabCouplingError.invalidChannel(
                    "output binding \(binding.channelID) is unresolved, duplicated, or out of range"
                )
            }
            preparedOutputs.append(.init(
                source: binding,
                speciesIndex: species.index,
                publicationOffset: outputOffset
            ))
            let addition = outputOffset.addingReportingOverflow(binding.laneIndices.count)
            guard !addition.overflow else {
                throw VivoNumiLabCouplingError.invalidChannel(
                    "output publication count overflow"
                )
            }
            outputOffset = addition.partialValue
        }

        guard Set(participantConfiguration.inputs.map(\.channelID)).count == participantConfiguration.inputs.count,
              Set(participantConfiguration.outputs.map(\.channelID)).count == participantConfiguration.outputs.count else {
            throw VivoNumiLabCouplingError.invalidChannel(
                "NumiVivo participant bindings contain duplicate channel ids"
            )
        }

        let runtime = try await VivoTransactionalMolecularRuntime.make(
            pack: pack,
            configuration: runtimeConfiguration,
            device: requestedDevice
        )
        return VivoNumiVivoCoupledParticipant(
            runtime: runtime,
            minimumTimeStep: runtimeConfiguration.minimumTimeStep,
            inputBindings: preparedInputs,
            outputBindings: preparedOutputs,
            programFingerprint: pack.header.contentFingerprint,
            sourceProgramFingerprint: pack.header.sourceFingerprint,
            capabilities: runtime.capabilities
        )
    }

    private init(
        runtime: VivoTransactionalMolecularRuntime,
        minimumTimeStep: Float,
        inputBindings: [PreparedInputBinding],
        outputBindings: [PreparedOutputBinding],
        programFingerprint: VivoFingerprint,
        sourceProgramFingerprint: VivoFingerprint,
        capabilities: VivoMetalCapabilities
    ) {
        self.runtime = runtime
        self.minimumTimeStep = minimumTimeStep
        self.inputBindings = inputBindings
        self.outputBindings = outputBindings
        self.programFingerprint = programFingerprint
        self.sourceProgramFingerprint = sourceProgramFingerprint
        self.capabilities = capabilities
    }

    public func prepare(
        _ transaction: VivoCouplingTransaction
    ) async throws -> VivoPreparedCouplingState {
        if let releaseFailure {
            throw VivoNumiLabCouplingError.stagingFailure(
                "NumiVivo participant is faulted after release failure: \(releaseFailure)"
            )
        }
        guard pending == nil else {
            throw VivoNumiLabCouplingError.invalidTransaction(
                "NumiVivo participant already owns a prepared transaction"
            )
        }
        try transaction.clock.validate()
        guard transaction.clock.deltaTime <= Double(Float.greatestFiniteMagnitude),
              transaction.clock.deltaTime >= Double(minimumTimeStep) else {
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .needsSmallerStep,
                maximumAcceptedStep: Double(minimumTimeStep),
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["requested coupling step is outside NumiVivo FP32/runtime bounds"]
            )
        }

        let samples = Dictionary(
            uniqueKeysWithValues: transaction.inputs.map { ($0.channelID, $0) }
        )
        var coupling: [VivoCouplingUpdate] = []
        for binding in inputBindings {
            guard let sample = samples[binding.source.channelID] else {
                throw VivoNumiLabCouplingError.invalidTransaction(
                    "required NumiVivo input channel \(binding.source.channelID) is absent"
                )
            }
            guard sample.unit == binding.source.unit,
                  sample.values.count == binding.source.laneIndices.count else {
                throw VivoNumiLabCouplingError.invalidOutput(
                    "NumiVivo input channel \(binding.source.channelID) unit or shape does not match its binding"
                )
            }
            for (lane, rawValue) in zip(binding.source.laneIndices, sample.values) {
                coupling.append(VivoCouplingUpdate(
                    speciesIndex: binding.speciesIndex,
                    laneIndex: lane,
                    mode: binding.source.mode,
                    value: try binding.source.transform.apply(
                        rawValue,
                        subject: binding.source.channelID
                    )
                ))
            }
        }

        let publicationRequests = outputBindings.flatMap { binding in
            binding.source.laneIndices.map { lane in
                VivoPublicationRequest(
                    speciesIndex: binding.speciesIndex,
                    laneIndex: lane
                )
            }
        }
        let runtimeTransactionID = UUID()
        let prepared = try await runtime.prepareStep(
            VivoStepRequest(
                timeStep: Float(transaction.clock.deltaTime),
                coupling: coupling,
                publications: publicationRequests,
                permitAdaptiveReduction: true
            ),
            transactionID: runtimeTransactionID
        )

        if prepared.disposition == .prepared,
           prepared.candidateTimeStep < Float(transaction.clock.deltaTime) {
            try await runtime.discardPreparedStep(transactionID: runtimeTransactionID)
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .needsSmallerStep,
                maximumAcceptedStep: Double(prepared.candidateTimeStep),
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["NumiVivo internally reduced the candidate time step"]
            )
        }

        switch prepared.disposition {
        case .prepared:
            let outputs = try makeOutputs(
                transaction: transaction,
                publications: prepared.publications
            )
            pending = Pending(
                transactionID: transaction.id,
                runtimeTransactionID: runtimeTransactionID,
                prepared: prepared,
                phase: .prepared
            )
            let digest = try VivoCanonicalJSON.fingerprint(
                VivoCanonicalJSON.encode(StateDigestMaterial(
                    transactionID: transaction.id,
                    programFingerprint: programFingerprint,
                    stepIndex: prepared.stepIndex,
                    candidateTimeStep: prepared.candidateTimeStep,
                    runtimeStatus: prepared.status,
                    outputs: outputs
                ))
            )
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .prepared,
                maximumAcceptedStep: Double(prepared.candidateTimeStep),
                outputs: outputs,
                stateDigest: String(describing: digest),
                diagnostics: []
            )
        case .requiresSmallerStep:
            let proposed = max(
                Double(minimumTimeStep),
                Double(prepared.candidateTimeStep) * 0.5
            )
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .needsSmallerStep,
                maximumAcceptedStep: proposed,
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["NumiVivo candidate requires a smaller time step"]
            )
        case .rejected:
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .rejected,
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["NumiVivo rejected its molecular candidate"]
            )
        case .reversibleShutdown:
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .reversibleShutdown,
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["NumiVivo requested reversible shutdown"]
            )
        case .permanentShutdown:
            return VivoPreparedCouplingState(
                participant: .numiVivo,
                transactionID: transaction.id,
                disposition: .permanentShutdown,
                stateDigest: String(describing: programFingerprint),
                diagnostics: ["NumiVivo requested permanent shutdown"]
            )
        }
    }

    public func stageCommit(transactionID: String) async throws {
        guard var pending,
              pending.transactionID == transactionID,
              pending.prepared.canCommit,
              pending.phase == .prepared else {
            throw VivoNumiLabCouplingError.stagingFailure(
                "NumiVivo has no matching commit-eligible prepared candidate"
            )
        }
        pending.phase = .staged
        self.pending = pending
    }

    public func releaseCommit(transactionID: String) async {
        guard let pending,
              pending.transactionID == transactionID,
              pending.phase == .staged else {
            releaseFailure = "release requested without a matching staged candidate"
            return
        }
        do {
            _ = try await runtime.commitPreparedStep(
                transactionID: pending.runtimeTransactionID
            )
            self.pending = nil
        } catch {
            releaseFailure = error.localizedDescription
            self.pending = nil
        }
    }

    public func rollback(transactionID: String) async {
        guard let pending,
              pending.transactionID == transactionID else { return }
        try? await runtime.discardPreparedStep(
            transactionID: pending.runtimeTransactionID
        )
        self.pending = nil
    }

    public func snapshot() async throws -> VivoStateSnapshot {
        try await runtime.snapshot()
    }

    public func failureState() -> String? { releaseFailure }

    private func makeOutputs(
        transaction: VivoCouplingTransaction,
        publications: [Float]
    ) throws -> [VivoCouplingSample] {
        let expectedCount = outputBindings.reduce(0) {
            $0 + $1.source.laneIndices.count
        }
        guard publications.count == expectedCount else {
            throw VivoNumiLabCouplingError.invalidOutput(
                "NumiVivo publication count does not match configured output bindings"
            )
        }

        return try outputBindings.map { binding in
            let count = binding.source.laneIndices.count
            let range = binding.publicationOffset..<(binding.publicationOffset + count)
            let values = try publications[range].map {
                try binding.source.transform.apply(
                    $0,
                    subject: binding.source.channelID
                )
            }
            return VivoCouplingSample(
                channelID: binding.source.channelID,
                producer: .numiVivo,
                unit: binding.source.unit,
                values: values,
                sourceStep: UInt64(transaction.clock.step),
                sourceTime: transaction.clock.time + transaction.clock.deltaTime
            )
        }
    }
}
