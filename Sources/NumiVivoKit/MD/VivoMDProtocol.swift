import Foundation

public enum VivoMDProtocolStageKind: String, Codable, Sendable { case minimization, dynamics }
public enum VivoMDProtocolPhase: String, Codable, Sendable { case running, stageFinished, blocked }
public enum VivoMDProtocolDisposition: String, Codable, Sendable { case completed, rejected, cancelled, failed }

public struct VivoMDProtocolStage: Codable, Sendable, Equatable {
    public var identifier: String
    public var kind: VivoMDProtocolStageKind
    public var configuration: VivoMDConfiguration
    public var steps: UInt64
    public var minimization: VivoMDMinimizationConfiguration?
    public var requireMinimizationConvergence: Bool
    public var velocityInitialization: VivoMDVelocityInitialization
    public var thermalizationSeed: UInt64?
    public var sampleEvery: UInt64?
    public var observablesEvery: UInt64?
    public var checkpointEvery: UInt64

    public init(identifier: String, kind: VivoMDProtocolStageKind,
                configuration: VivoMDConfiguration, steps: UInt64 = 0,
                minimization: VivoMDMinimizationConfiguration? = nil,
                requireMinimizationConvergence: Bool = true,
                velocityInitialization: VivoMDVelocityInitialization = .preserve,
                thermalizationSeed: UInt64? = nil, sampleEvery: UInt64? = nil,
                observablesEvery: UInt64? = 1000, checkpointEvery: UInt64 = 10_000) {
        self.identifier = identifier; self.kind = kind; self.configuration = configuration
        self.steps = steps; self.minimization = minimization
        self.requireMinimizationConvergence = requireMinimizationConvergence
        self.velocityInitialization = velocityInitialization; self.thermalizationSeed = thermalizationSeed
        self.sampleEvery = sampleEvery; self.observablesEvery = observablesEvery
        self.checkpointEvery = checkpointEvery
    }
    public func validate() throws {
        try configuration.validate()
        guard !identifier.isEmpty, identifier.utf8.count <= 128,
              checkpointEvery > 0, sampleEvery != 0, observablesEvery != 0 else {
            throw VivoArtifactValidationError.invalid("invalid MD stage identifier or sampling schedule")
        }
        if kind == .minimization {
            guard steps == 0, let minimization, sampleEvery == nil, observablesEvery == nil,
                  velocityInitialization == .zero else {
                throw VivoArtifactValidationError.invalid("minimization requires explicit settings, zero velocities, no dynamics steps or interval samples")
            }
            try minimization.validate()
        } else {
            guard steps > 0, steps <= 10_000_000_000, minimization == nil else {
                throw VivoArtifactValidationError.invalid("dynamics requires a bounded positive step count and no minimizer")
            }
        }
        if velocityInitialization == .maxwellBoltzmann {
            guard thermalizationSeed != nil, configuration.targetTemperatureK != nil else {
                throw VivoArtifactValidationError.invalid("thermal initialization requires an explicit seed and target temperature")
            }
        } else if thermalizationSeed != nil {
            throw VivoArtifactValidationError.invalid("unused stage thermalization seed")
        }
    }
}

public struct VivoMDProtocolPlan: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/md-protocol/v1"
    public var schema: String
    public var numericalContract: String
    public var identifier: String
    public var systemFingerprint: VivoFingerprint
    public var stages: [VivoMDProtocolStage]
    public var trajectoryIncludesVelocities: Bool
    public var trajectoryChunkBytes: Int
    public var metadata: [String: String]

    public init(identifier: String, systemFingerprint: VivoFingerprint, stages: [VivoMDProtocolStage],
                trajectoryIncludesVelocities: Bool = false, trajectoryChunkBytes: Int = 8 * 1024 * 1024,
                metadata: [String: String] = [:]) {
        schema = Self.schemaID; numericalContract = VivoMDExecutionIdentity.current
        self.identifier = identifier; self.systemFingerprint = systemFingerprint; self.stages = stages
        self.trajectoryIncludesVelocities = trajectoryIncludesVelocities
        self.trajectoryChunkBytes = trajectoryChunkBytes; self.metadata = metadata
    }
    public func validate() throws {
        guard schema == Self.schemaID, numericalContract == VivoMDExecutionIdentity.current,
              !identifier.isEmpty, identifier.utf8.count <= 256, !stages.isEmpty, stages.count <= 128,
              Set(stages.map(\.identifier)).count == stages.count,
              trajectoryChunkBytes > VivoMDTrajectoryChunkCodec.headerBytes,
              trajectoryChunkBytes <= VivoMDTrajectoryChunkCodec.maximumChunkBytes,
              metadata.count <= 64, metadata.allSatisfy({ $0.key.utf8.count <= 128 && $0.value.utf8.count <= 4096 }) else {
            throw VivoArtifactValidationError.invalid("invalid MD protocol identity, stages or storage policy")
        }
        var total: UInt64 = 0
        for stage in stages {
            try stage.validate()
            let next = total.addingReportingOverflow(stage.steps)
            guard !next.overflow else { throw VivoArtifactValidationError.invalid("protocol step count overflow") }
            total = next.partialValue
        }
        guard stages[0].velocityInitialization != .preserve else {
            throw VivoArtifactValidationError.invalid("initial-state input has no velocities; first stage must explicitly initialize them")
        }
    }
    public func fingerprint() throws -> VivoFingerprint {
        try validate(); return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
    public func validate(system: VivoClassicalSystem, initialState: VivoClassicalInitialState) throws -> [VivoMDCapabilityReport] {
        try validate()
        guard systemFingerprint == (try system.fingerprint()) else {
            throw VivoArtifactValidationError.incompatible("protocol identifies another classical system")
        }
        if stages.contains(where: { $0.sampleEvery != nil }) {
            _ = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: UInt32(system.particles.count),
                                                         includeVelocities: trajectoryIncludesVelocities)
        }
        return try stages.map { try VivoMDCapabilityAnalyzer.analyze(system: system, initialState: initialState,
                                                                     configuration: $0.configuration) }
    }

    /// Illustrative preparation schedule, NOT the paper's measured protocol or
    /// an assertion of equilibration. Users set durations and convergence gates.
    public static func preparationTemplate(systemFingerprint: VivoFingerprint,
                                            nvtSteps: UInt64 = 10_000, nptSteps: UInt64 = 10_000,
                                            productionSteps: UInt64 = 100_000) -> Self {
        let minimum = VivoMDConfiguration(ensemble: .nve, thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil)
        let nvt = VivoMDConfiguration()
        let npt = VivoMDConfiguration(ensemble: .npt, barostat: .monteCarloIsotropic, targetPressureBar: 1)
        return .init(identifier: "prepare-and-sample", systemFingerprint: systemFingerprint, stages: [
            .init(identifier: "minimize", kind: .minimization, configuration: minimum,
                  minimization: .init(), velocityInitialization: .zero, observablesEvery: nil),
            .init(identifier: "nvt", kind: .dynamics, configuration: nvt, steps: nvtSteps,
                  velocityInitialization: .maxwellBoltzmann, thermalizationSeed: 0x4e5654494e495431,
                  sampleEvery: 1000, observablesEvery: 1000),
            .init(identifier: "npt", kind: .dynamics, configuration: npt, steps: nptSteps,
                  sampleEvery: 1000, observablesEvery: 1000),
            .init(identifier: "production", kind: .dynamics, configuration: npt, steps: productionSteps,
                  sampleEvery: 1000, observablesEvery: 1000)
        ], metadata: ["status": "illustrative source workflow; not a validated equilibration protocol"])
    }
}

public struct VivoMDObservationLink: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/md-observation-link/v1"
    public let schema: String
    public let stageIdentifier: String
    public let ordinal: UInt64
    public let previous: VivoFingerprint?
    public let observation: VivoMDObservables
}

public struct VivoMDProtocolStageReport: Codable, Sendable, Equatable {
    public let schema: String
    public let planFingerprint: VivoFingerprint
    public let stageIdentifier: String
    public let stageIndex: Int
    public let successful: Bool
    public let entryCheckpoint: VivoFingerprint
    public let exitCheckpoint: VivoFingerprint
    public let transition: VivoFingerprint?
    public let committedSteps: UInt64
    public let startTimePS: Double
    public let endTimePS: Double
    public let trajectoryManifest: VivoFingerprint?
    public let observationTail: VivoFingerprint?
    public let observationCount: UInt64
    public let minimization: VivoMDMinimizationCertificate?
    public let rejected: VivoMDStepCertificate?
}

/// Small immutable restart cursor. Large coordinate state and sampled output
/// are separate verified objects. A cursor never points at buffered output.
public struct VivoMDProtocolCheckpoint: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/md-protocol-checkpoint/v1"
    public var schema: String
    public var numericalContract: String
    public var runID: UUID
    public var resumedFrom: VivoFingerprint?
    public var planFingerprint: VivoFingerprint
    public var systemFingerprint: VivoFingerprint
    public var initialStateArtifact: VivoFingerprint
    public var stageIndex: Int
    public var phase: VivoMDProtocolPhase
    public var completedStepsInStage: UInt64
    public var entryCheckpoint: VivoFingerprint
    public var currentCheckpoint: VivoFingerprint
    public var transition: VivoFingerprint?
    public var trajectoryManifest: VivoFingerprint?
    public var observationTail: VivoFingerprint?
    public var observationCount: UInt64
    public var priorStageReports: [VivoFingerprint]
    public var activeStageReport: VivoFingerprint?

    public func validate(plan: VivoMDProtocolPlan) throws {
        guard schema == Self.schemaID, numericalContract == VivoMDExecutionIdentity.current,
              planFingerprint == (try plan.fingerprint()), systemFingerprint == plan.systemFingerprint,
              stageIndex >= 0, stageIndex < plan.stages.count, priorStageReports.count == stageIndex,
              completedStepsInStage <= plan.stages[stageIndex].steps,
              (observationCount == 0) == (observationTail == nil),
              (phase == .running) == (activeStageReport == nil) else {
            throw VivoArtifactValidationError.incompatible("MD protocol checkpoint violates identity, stage or progress invariants")
        }
        let stage = plan.stages[stageIndex]
        if stage.sampleEvery == nil, trajectoryManifest != nil {
            throw VivoArtifactValidationError.invalid("unexpected trajectory for an unsampled stage")
        }
        if phase == .stageFinished, stage.kind == .dynamics, completedStepsInStage != stage.steps {
            throw VivoArtifactValidationError.invalid("finished dynamics stage has incomplete step count")
        }
    }
}

public struct VivoMDProtocolRunReceipt: Codable, Sendable, Equatable {
    public let schema: String
    public let runID: UUID
    public let planFingerprint: VivoFingerprint
    public let disposition: VivoMDProtocolDisposition
    public let completedStages: Int
    public let latestDurableCheckpoint: VivoStoredArtifact?
    public let checkpointReference: String
    public let diagnostic: String?
    public let persistenceDiagnostic: String?
}
