import Foundation
import Metal

public struct VivoMolecularPhysiologyStepRequest: Codable, Sendable, Equatable {
    public var timeStepSeconds: Double?
    public var molecularCoupling: [VivoCouplingUpdate]
    public var molecularPublications: [VivoPublicationRequest]
    public var physiologyPreUpdates: [VivoPhysiologyStateUpdate]
    public var physiologyPostUpdates: [VivoPhysiologyStateUpdate]
    public var physiologyPublications: [VivoPhysiologyPublicationRequest]
    public var permitAdaptiveReduction: Bool

    public init(
        timeStepSeconds: Double? = nil,
        molecularCoupling: [VivoCouplingUpdate] = [],
        molecularPublications: [VivoPublicationRequest] = [],
        physiologyPreUpdates: [VivoPhysiologyStateUpdate] = [],
        physiologyPostUpdates: [VivoPhysiologyStateUpdate] = [],
        physiologyPublications: [VivoPhysiologyPublicationRequest] = [],
        permitAdaptiveReduction: Bool = true
    ) {
        self.timeStepSeconds = timeStepSeconds
        self.molecularCoupling = molecularCoupling
        self.molecularPublications = molecularPublications
        self.physiologyPreUpdates = physiologyPreUpdates
        self.physiologyPostUpdates = physiologyPostUpdates
        self.physiologyPublications = physiologyPublications
        self.permitAdaptiveReduction = permitAdaptiveReduction
    }
}

public enum VivoMolecularPhysiologyDisposition: String, Codable, Sendable {
    case committed
    case committedWithReducedStep
    case committedWithoutConvergence
    case rejected
    case reversibleShutdown
    case permanentShutdown
}

public struct VivoMolecularPhysiologyCertificate: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoMolecularPhysiologyDisposition
    public let bridgeFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let physiologyFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let timeBeforeSeconds: Double
    public let timeAfterSeconds: Double
    public let requestedTimeStepSeconds: Double
    public let acceptedTimeStepSeconds: Double?
    public let stepNegotiationCount: UInt32
    public let fixedPointIterationCount: UInt32
    public let finalResidual: Double
    public let converged: Bool
    public let molecularStatus: VivoRuntimeStatus?
    public let physiologyStatus: VivoPhysiologyRuntimeStatus?
    public let message: String

    public var committed: Bool {
        disposition == .committed ||
            disposition == .committedWithReducedStep ||
            disposition == .committedWithoutConvergence
    }
}

public struct VivoMolecularPhysiologyStepResult: Codable, Sendable, Equatable {
    public let certificate: VivoMolecularPhysiologyCertificate
    public let molecularCertificate: VivoMolecularTransactionCertificate?
    public let physiologyCertificate: VivoPhysiologyStepCertificate?
    public let molecularEvents: [VivoEvent]
    public let bridgeExposureValues: [Float]
    public let bridgeFeedbackValues: [Float]
    public let molecularPublications: [Float]
    public let physiologyPublications: [Float]
}

public actor VivoMolecularPhysiologyCoordinator {
    public nonisolated let bridge: PreparedVivoMolecularPhysiologyBridge
    public nonisolated let programPack: VivoProgramPack
    public nonisolated let physiologyModel: PreparedVivoPhysiologyModel
    public nonisolated let molecularConfiguration: VivoRuntimeConfiguration
    public nonisolated let physiologyConfiguration: VivoPhysiologyRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities

    private let molecular: VivoTransactionalMolecularRuntime
    private let physiology: VivoPhysiologyRuntime
    private var acceptedStepIndex: UInt32 = 0
    private var previousFeedbackValues: [Float]
    private var hasAcceptedFeedback = false

    public static func make(
        programPack: VivoProgramPack,
        molecularConfiguration: VivoRuntimeConfiguration,
        physiologyModel: PreparedVivoPhysiologyModel,
        physiologyConfiguration: VivoPhysiologyRuntimeConfiguration = .init(),
        bridge: PreparedVivoMolecularPhysiologyBridge,
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoMolecularPhysiologyCoordinator {
        guard bridge.schema == PreparedVivoMolecularPhysiologyBridge.schema,
              bridge.programFingerprint == programPack.header.contentFingerprint,
              bridge.physiologyFingerprint == physiologyModel.fingerprint else {
            throw VivoArtifactValidationError.incompatible(
                "coordinator inputs do not share program and physiology fingerprints"
            )
        }
        try bridge.policy.validate()
        try physiologyConfiguration.validate(for: physiologyModel)
        try validateBridgeBounds(
            bridge,
            programPack: programPack,
            physiologyModel: physiologyModel,
            laneCount: molecularConfiguration.laneCount
        )

        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        async let molecularRuntime = VivoTransactionalMolecularRuntime.make(
            pack: programPack,
            configuration: molecularConfiguration,
            device: device
        )
        async let physiologyRuntime = VivoPhysiologyRuntime.make(
            model: physiologyModel,
            configuration: physiologyConfiguration,
            device: device
        )
        let (molecular, physiology) = try await (
            molecularRuntime,
            physiologyRuntime
        )
        return VivoMolecularPhysiologyCoordinator(
            bridge: bridge,
            programPack: programPack,
            physiologyModel: physiologyModel,
            molecularConfiguration: molecularConfiguration,
            physiologyConfiguration: physiologyConfiguration,
            capabilities: VivoMetalCapabilities(device: device),
            molecular: molecular,
            physiology: physiology
        )
    }

    private init(
        bridge: PreparedVivoMolecularPhysiologyBridge,
        programPack: VivoProgramPack,
        physiologyModel: PreparedVivoPhysiologyModel,
        molecularConfiguration: VivoRuntimeConfiguration,
        physiologyConfiguration: VivoPhysiologyRuntimeConfiguration,
        capabilities: VivoMetalCapabilities,
        molecular: VivoTransactionalMolecularRuntime,
        physiology: VivoPhysiologyRuntime
    ) {
        self.bridge = bridge
        self.programPack = programPack
        self.physiologyModel = physiologyModel
        self.molecularConfiguration = molecularConfiguration
        self.physiologyConfiguration = physiologyConfiguration
        self.capabilities = capabilities
        self.molecular = molecular
        self.physiology = physiology
        self.previousFeedbackValues = [Float](
            repeating: 0,
            count: bridge.molecularToPhysiology.count
        )
    }

    public func stepIndex() -> UInt32 { acceptedStepIndex }

    public func timeSeconds() async throws -> Double {
        try await synchronizedTimes().physiology
    }

    public func molecularSnapshot() async throws -> VivoStateSnapshot {
        try await molecular.snapshot()
    }

    public func physiologySnapshot() async throws -> VivoPhysiologySnapshot {
        try await physiology.snapshot()
    }

    public func physiologyCheckpoint() async throws -> VivoPhysiologyCheckpoint {
        try await physiology.checkpoint()
    }

    public func step(
        _ request: VivoMolecularPhysiologyStepRequest = .init()
    ) async throws -> VivoMolecularPhysiologyStepResult {
        let transactionID = UUID()
        let times = try await synchronizedTimes()
        let timeBefore = times.physiology
        let minimumStep = max(
            Double(molecularConfiguration.minimumTimeStep),
            physiologyModel.minimumTimeStepSeconds
        )
        let maximumStep = min(
            Double(molecularConfiguration.maximumTimeStep),
            physiologyModel.maximumTimeStepSeconds
        )
        let requested = request.timeStepSeconds ?? min(
            Double(molecularConfiguration.timeStep),
            physiologyModel.preferredTimeStepSeconds
        )
        guard requested.isFinite,
              requested >= minimumStep,
              requested <= maximumStep,
              requested <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoRuntimeError.invalidConfiguration(
                "joint time step lies outside the molecular/physiology intersection"
            )
        }
        try validateExternalRequest(request)

        var dt = canonicalTimeStep(requested)
        var negotiation: UInt32 = 0
        var lastMolecular: VivoPreparedMolecularStep?
        var lastPhysiology: VivoPreparedPhysiologyStep?
        var lastResidual = Double.greatestFiniteMagnitude
        var lastExposure: [Float] = []
        var lastFeedback: [Float] = []

        negotiationLoop: while negotiation < bridge.policy.maximumStepNegotiations {
            var feedbackGuess = previousFeedbackValues
            if feedbackGuess.count != bridge.molecularToPhysiology.count {
                feedbackGuess = [Float](
                    repeating: 0,
                    count: bridge.molecularToPhysiology.count
                )
            }
            var hasFeedbackGuess = hasAcceptedFeedback || bridge.molecularToPhysiology.isEmpty

            for fixedPointIteration in 0..<bridge.policy.maximumFixedPointIterations {
                var physiologyPending = false
                var molecularPending = false
                do {
                    let physiologyPublications = bridge.physiologyToMolecular.map {
                        VivoPhysiologyPublicationRequest(
                            pairIndex: $0.pairIndex,
                            environmentIndex: $0.environmentIndex
                        )
                    } + request.physiologyPublications
                    let feedbackUpdates: [VivoPhysiologyStateUpdate]
                    if hasFeedbackGuess {
                        feedbackUpdates = zip(
                            bridge.molecularToPhysiology,
                            feedbackGuess
                        ).map { link, value in
                            VivoPhysiologyStateUpdate(
                                pairIndex: link.pairIndex,
                                environmentIndex: link.environmentIndex,
                                mode: link.physiologyMode,
                                value: value
                            )
                        }
                    } else {
                        feedbackUpdates = []
                    }

                    let physiologyPrepared = try await physiology.prepareStep(
                        VivoPhysiologyPrepareRequest(
                            timeStepSeconds: dt,
                            preUpdates: request.physiologyPreUpdates,
                            postUpdates: request.physiologyPostUpdates + feedbackUpdates,
                            publications: physiologyPublications,
                            permitAdaptiveReduction: request.permitAdaptiveReduction
                        ),
                        transactionID: transactionID
                    )
                    lastPhysiology = physiologyPrepared
                    physiologyPending = physiologyPrepared.canCommit

                    if !physiologyPrepared.canCommit ||
                        !sameTimeStep(physiologyPrepared.candidateTimeStep, dt) {
                        if physiologyPending {
                            try await physiology.discardPreparedStep(
                                transactionID: transactionID
                            )
                            physiologyPending = false
                        }
                        guard let reduced = reduction(
                            current: dt,
                            preparedCandidate: physiologyPrepared.canCommit
                                ? physiologyPrepared.candidateTimeStep
                                : nil,
                            runtimeSuggestion: Double(
                                physiologyPrepared.status.suggestedTimeStep
                            ),
                            permitted: request.permitAdaptiveReduction,
                            minimum: minimumStep
                        ) else {
                            return rejectedResult(
                                transactionID: transactionID,
                                disposition: .rejected,
                                timeBefore: timeBefore,
                                requested: requested,
                                negotiation: negotiation + 1,
                                iteration: fixedPointIteration + 1,
                                residual: lastResidual,
                                molecular: lastMolecular,
                                physiology: physiologyPrepared,
                                exposure: lastExposure,
                                feedback: lastFeedback,
                                message: "physiology candidate rejected the common time step"
                            )
                        }
                        dt = reduced
                        negotiation &+= 1
                        continue negotiationLoop
                    }

                    let exposureCount = bridge.physiologyToMolecular.count
                    guard physiologyPrepared.publications.count >= exposureCount else {
                        throw VivoRuntimeError.commandFailed(
                            "physiology candidate did not publish every bridge exposure"
                        )
                    }
                    let rawExposure = Array(
                        physiologyPrepared.publications.prefix(exposureCount)
                    )
                    let exposure = try zip(
                        bridge.physiologyToMolecular,
                        rawExposure
                    ).map { link, value in
                        try link.transfer.apply(value, subject: link.identifier)
                    }
                    lastExposure = exposure

                    let bridgeCoupling = zip(
                        bridge.physiologyToMolecular,
                        exposure
                    ).map { link, value in
                        VivoCouplingUpdate(
                            speciesIndex: link.molecularSpeciesIndex,
                            laneIndex: link.molecularLaneIndex,
                            mode: link.molecularMode,
                            value: value
                        )
                    }
                    try validateMolecularDestinations(
                        bridgeCoupling + request.molecularCoupling
                    )
                    let molecularPublications = bridge.molecularToPhysiology.map {
                        VivoPublicationRequest(
                            speciesIndex: $0.molecularSpeciesIndex,
                            laneIndex: $0.molecularLaneIndex
                        )
                    } + request.molecularPublications
                    let molecularPrepared = try await molecular.prepareStep(
                        VivoStepRequest(
                            timeStep: Float(dt),
                            coupling: bridgeCoupling + request.molecularCoupling,
                            publications: molecularPublications,
                            permitAdaptiveReduction: request.permitAdaptiveReduction
                        ),
                        transactionID: transactionID
                    )
                    lastMolecular = molecularPrepared
                    molecularPending = molecularPrepared.canCommit

                    if molecularPrepared.disposition == .permanentShutdown ||
                        molecularPrepared.disposition == .reversibleShutdown {
                        if physiologyPending {
                            try await physiology.discardPreparedStep(
                                transactionID: transactionID
                            )
                            physiologyPending = false
                        }
                        return rejectedResult(
                            transactionID: transactionID,
                            disposition: molecularPrepared.disposition == .permanentShutdown
                                ? .permanentShutdown
                                : .reversibleShutdown,
                            timeBefore: timeBefore,
                            requested: requested,
                            negotiation: negotiation + 1,
                            iteration: fixedPointIteration + 1,
                            residual: lastResidual,
                            molecular: molecularPrepared,
                            physiology: physiologyPrepared,
                            exposure: exposure,
                            feedback: lastFeedback,
                            message: "molecular runtime requested shutdown"
                        )
                    }

                    if !molecularPrepared.canCommit ||
                        !sameTimeStep(Double(molecularPrepared.candidateTimeStep), dt) {
                        if molecularPending {
                            try await molecular.discardPreparedStep(
                                transactionID: transactionID
                            )
                            molecularPending = false
                        }
                        if physiologyPending {
                            try await physiology.discardPreparedStep(
                                transactionID: transactionID
                            )
                            physiologyPending = false
                        }
                        guard let reduced = reduction(
                            current: dt,
                            preparedCandidate: molecularPrepared.canCommit
                                ? Double(molecularPrepared.candidateTimeStep)
                                : nil,
                            runtimeSuggestion: 0,
                            permitted: request.permitAdaptiveReduction,
                            minimum: minimumStep
                        ) else {
                            return rejectedResult(
                                transactionID: transactionID,
                                disposition: .rejected,
                                timeBefore: timeBefore,
                                requested: requested,
                                negotiation: negotiation + 1,
                                iteration: fixedPointIteration + 1,
                                residual: lastResidual,
                                molecular: molecularPrepared,
                                physiology: physiologyPrepared,
                                exposure: exposure,
                                feedback: lastFeedback,
                                message: "molecular candidate rejected the common time step"
                            )
                        }
                        dt = reduced
                        negotiation &+= 1
                        continue negotiationLoop
                    }

                    let feedbackCount = bridge.molecularToPhysiology.count
                    guard molecularPrepared.publications.count >= feedbackCount else {
                        throw VivoRuntimeError.commandFailed(
                            "molecular candidate did not publish every bridge feedback value"
                        )
                    }
                    let rawFeedback = Array(
                        molecularPrepared.publications.prefix(feedbackCount)
                    )
                    let feedback = try zip(
                        bridge.molecularToPhysiology,
                        rawFeedback
                    ).map { link, value in
                        try link.transfer.apply(value, subject: link.identifier)
                    }
                    lastFeedback = feedback
                    let residual = hasFeedbackGuess
                        ? normalizedResidual(candidate: feedback, reference: feedbackGuess)
                        : Double.greatestFiniteMagnitude
                    lastResidual = residual
                    let finalIteration = fixedPointIteration + 1 >=
                        bridge.policy.maximumFixedPointIterations
                    let converged = feedback.isEmpty ||
                        (hasFeedbackGuess && residual <= bridge.policy.convergenceTolerance)
                    let permitUnconvergedCommit = finalIteration &&
                        !bridge.policy.requireConvergence

                    if converged || permitUnconvergedCommit {
                        async let molecularCommit = molecular.commitPreparedStep(
                            transactionID: transactionID
                        )
                        async let physiologyCommit = physiology.commitPreparedStep(
                            transactionID: transactionID
                        )
                        let (molecularCertificate, physiologyCertificate) = try await (
                            molecularCommit,
                            physiologyCommit
                        )
                        molecularPending = false
                        physiologyPending = false
                        previousFeedbackValues = feedback
                        hasAcceptedFeedback = true
                        let step = acceptedStepIndex
                        acceptedStepIndex &+= 1

                        let molecularTail = Array(
                            molecularPrepared.publications.dropFirst(feedbackCount)
                        )
                        let physiologyTail = Array(
                            physiologyPrepared.publications.dropFirst(exposureCount)
                        )
                        let reduced = !sameTimeStep(dt, requested)
                        let disposition: VivoMolecularPhysiologyDisposition
                        if !converged {
                            disposition = .committedWithoutConvergence
                        } else if reduced {
                            disposition = .committedWithReducedStep
                        } else {
                            disposition = .committed
                        }
                        return VivoMolecularPhysiologyStepResult(
                            certificate: VivoMolecularPhysiologyCertificate(
                                transactionID: transactionID,
                                disposition: disposition,
                                bridgeFingerprint: bridge.fingerprint,
                                programFingerprint: bridge.programFingerprint,
                                physiologyFingerprint: bridge.physiologyFingerprint,
                                stepIndex: step,
                                timeBeforeSeconds: timeBefore,
                                timeAfterSeconds: physiologyCertificate.timeAfter,
                                requestedTimeStepSeconds: requested,
                                acceptedTimeStepSeconds: dt,
                                stepNegotiationCount: negotiation + 1,
                                fixedPointIterationCount: fixedPointIteration + 1,
                                finalResidual: residual,
                                converged: converged,
                                molecularStatus: molecularPrepared.status,
                                physiologyStatus: physiologyPrepared.status,
                                message: converged
                                    ? "joint candidate committed"
                                    : "joint candidate committed under non-convergence policy"
                            ),
                            molecularCertificate: molecularCertificate,
                            physiologyCertificate: physiologyCertificate,
                            molecularEvents: molecularPrepared.events,
                            bridgeExposureValues: exposure,
                            bridgeFeedbackValues: feedback,
                            molecularPublications: molecularTail,
                            physiologyPublications: physiologyTail
                        )
                    }

                    try await molecular.discardPreparedStep(
                        transactionID: transactionID
                    )
                    molecularPending = false
                    try await physiology.discardPreparedStep(
                        transactionID: transactionID
                    )
                    physiologyPending = false
                    if hasFeedbackGuess {
                        feedbackGuess = relaxed(
                            previous: feedbackGuess,
                            candidate: feedback,
                            factor: bridge.policy.relaxation
                        )
                    } else {
                        feedbackGuess = feedback
                        hasFeedbackGuess = true
                    }
                } catch {
                    if molecularPending {
                        try? await molecular.discardPreparedStep(
                            transactionID: transactionID
                        )
                    }
                    if physiologyPending {
                        try? await physiology.discardPreparedStep(
                            transactionID: transactionID
                        )
                    }
                    throw error
                }
            }

            return rejectedResult(
                transactionID: transactionID,
                disposition: .rejected,
                timeBefore: timeBefore,
                requested: requested,
                negotiation: negotiation + 1,
                iteration: bridge.policy.maximumFixedPointIterations,
                residual: lastResidual,
                molecular: lastMolecular,
                physiology: lastPhysiology,
                exposure: lastExposure,
                feedback: lastFeedback,
                message: "molecular-physiology fixed-point iteration did not converge"
            )
        }

        return rejectedResult(
            transactionID: transactionID,
            disposition: .rejected,
            timeBefore: timeBefore,
            requested: requested,
            negotiation: negotiation,
            iteration: 0,
            residual: lastResidual,
            molecular: lastMolecular,
            physiology: lastPhysiology,
            exposure: lastExposure,
            feedback: lastFeedback,
            message: "common-step negotiation budget was exhausted"
        )
    }

    private static func validateBridgeBounds(
        _ bridge: PreparedVivoMolecularPhysiologyBridge,
        programPack: VivoProgramPack,
        physiologyModel: PreparedVivoPhysiologyModel,
        laneCount: UInt32
    ) throws {
        let speciesCount = programPack.runtimeContract.speciesCount
        for link in bridge.physiologyToMolecular {
            guard link.pairIndex < physiologyModel.pairCount,
                  link.environmentIndex < physiologyModel.environmentCount,
                  link.molecularSpeciesIndex < speciesCount,
                  link.molecularLaneIndex < laneCount else {
                throw VivoArtifactValidationError.invalid(
                    "prepared physiology-to-molecular link \(link.identifier) is out of bounds"
                )
            }
        }
        var replacementDestinations = Set<UInt64>()
        for link in bridge.molecularToPhysiology {
            guard link.molecularSpeciesIndex < speciesCount,
                  link.molecularLaneIndex < laneCount,
                  link.pairIndex < physiologyModel.pairCount,
                  link.environmentIndex < physiologyModel.environmentCount else {
                throw VivoArtifactValidationError.invalid(
                    "prepared molecular-to-physiology link \(link.identifier) is out of bounds"
                )
            }
            if link.physiologyMode == .replace {
                let key = UInt64(link.pairIndex) << 32 | UInt64(link.environmentIndex)
                guard replacementDestinations.insert(key).inserted else {
                    throw VivoArtifactValidationError.invalid(
                        "multiple bridge replacement links target one physiology state element"
                    )
                }
            }
        }
    }

    private func synchronizedTimes() async throws -> (
        molecular: Double,
        physiology: Double
    ) {
        async let molecularTime = molecular.time()
        async let physiologyTime = physiology.time()
        let (molecularValue, physiologyValue) = await (
            molecularTime,
            physiologyTime
        )
        let molecularSeconds = Double(molecularValue)
        let tolerance = synchronizationTolerance(
            at: max(abs(molecularSeconds), abs(physiologyValue))
        )
        guard abs(molecularSeconds - physiologyValue) <= tolerance else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular and physiology clocks differ by \(abs(molecularSeconds - physiologyValue)) seconds"
            )
        }
        return (molecularSeconds, physiologyValue)
    }

    private func validateExternalRequest(
        _ request: VivoMolecularPhysiologyStepRequest
    ) throws {
        for update in request.molecularCoupling {
            guard update.speciesIndex < programPack.runtimeContract.speciesCount,
                  update.laneIndex < molecularConfiguration.laneCount,
                  update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration(
                    "external molecular coupling update is invalid"
                )
            }
        }
        for publication in request.molecularPublications {
            guard publication.speciesIndex < programPack.runtimeContract.speciesCount,
                  publication.laneIndex < molecularConfiguration.laneCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "external molecular publication is invalid"
                )
            }
        }
        for update in request.physiologyPreUpdates + request.physiologyPostUpdates {
            guard update.pairIndex < physiologyModel.pairCount,
                  update.environmentIndex < physiologyModel.environmentCount,
                  update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration(
                    "external physiology update is invalid"
                )
            }
        }
        for publication in request.physiologyPublications {
            guard publication.pairIndex < physiologyModel.pairCount,
                  publication.environmentIndex < physiologyModel.environmentCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "external physiology publication is invalid"
                )
            }
        }
    }

    private func validateMolecularDestinations(
        _ updates: [VivoCouplingUpdate]
    ) throws {
        var destinations = Set<UInt64>()
        for update in updates {
            let key = UInt64(update.speciesIndex) << 32 | UInt64(update.laneIndex)
            guard destinations.insert(key).inserted else {
                throw VivoRuntimeError.invalidConfiguration(
                    "multiple coupling updates target molecular species \(update.speciesIndex) lane \(update.laneIndex)"
                )
            }
        }
    }

    private func reduction(
        current: Double,
        preparedCandidate: Double?,
        runtimeSuggestion: Double,
        permitted: Bool,
        minimum: Double
    ) -> Double? {
        guard permitted else { return nil }
        var candidates = [current * 0.5]
        if let preparedCandidate,
           preparedCandidate.isFinite,
           preparedCandidate > 0,
           preparedCandidate < current {
            candidates.append(preparedCandidate)
        }
        if runtimeSuggestion.isFinite,
           runtimeSuggestion > 0,
           runtimeSuggestion < current {
            candidates.append(runtimeSuggestion)
        }
        guard let raw = candidates.min() else { return nil }
        let canonical = canonicalTimeStep(raw)
        guard canonical >= minimum,
              canonical < current,
              canonical > 0 else {
            return nil
        }
        return canonical
    }

    private func normalizedResidual(
        candidate: [Float],
        reference: [Float]
    ) -> Double {
        guard candidate.count == reference.count,
              !candidate.isEmpty else {
            return candidate.isEmpty && reference.isEmpty
                ? 0
                : Double.greatestFiniteMagnitude
        }
        var maximum = 0.0
        for (newValue, oldValue) in zip(candidate, reference) {
            guard newValue.isFinite, oldValue.isFinite else {
                return Double.greatestFiniteMagnitude
            }
            let new = Double(newValue)
            let old = Double(oldValue)
            let denominator = max(1, max(abs(new), abs(old)))
            maximum = max(maximum, abs(new - old) / denominator)
        }
        return maximum
    }

    private func relaxed(
        previous: [Float],
        candidate: [Float],
        factor: Double
    ) -> [Float] {
        guard previous.count == candidate.count else { return candidate }
        return zip(previous, candidate).map { old, new in
            Float((1 - factor) * Double(old) + factor * Double(new))
        }
    }

    private func canonicalTimeStep(_ value: Double) -> Double {
        Double(Float(value))
    }

    private func sameTimeStep(_ lhs: Double, _ rhs: Double) -> Bool {
        abs(lhs - rhs) <= max(
            bridge.policy.timeSynchronizationTolerance,
            max(1, max(abs(lhs), abs(rhs))) * Double(Float.ulpOfOne) * 8
        )
    }

    private func synchronizationTolerance(at time: Double) -> Double {
        max(
            bridge.policy.timeSynchronizationTolerance,
            max(abs(time), 1) * Double(Float.ulpOfOne) * 32
        )
    }

    private func rejectedResult(
        transactionID: UUID,
        disposition: VivoMolecularPhysiologyDisposition,
        timeBefore: Double,
        requested: Double,
        negotiation: UInt32,
        iteration: UInt32,
        residual: Double,
        molecular: VivoPreparedMolecularStep?,
        physiology: VivoPreparedPhysiologyStep?,
        exposure: [Float],
        feedback: [Float],
        message: String
    ) -> VivoMolecularPhysiologyStepResult {
        VivoMolecularPhysiologyStepResult(
            certificate: VivoMolecularPhysiologyCertificate(
                transactionID: transactionID,
                disposition: disposition,
                bridgeFingerprint: bridge.fingerprint,
                programFingerprint: bridge.programFingerprint,
                physiologyFingerprint: bridge.physiologyFingerprint,
                stepIndex: acceptedStepIndex,
                timeBeforeSeconds: timeBefore,
                timeAfterSeconds: timeBefore,
                requestedTimeStepSeconds: requested,
                acceptedTimeStepSeconds: nil,
                stepNegotiationCount: negotiation,
                fixedPointIterationCount: iteration,
                finalResidual: residual,
                converged: false,
                molecularStatus: molecular?.status,
                physiologyStatus: physiology?.status,
                message: message
            ),
            molecularCertificate: nil,
            physiologyCertificate: nil,
            molecularEvents: molecular?.events ?? [],
            bridgeExposureValues: exposure,
            bridgeFeedbackValues: feedback,
            molecularPublications: [],
            physiologyPublications: []
        )
    }
}
