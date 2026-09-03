import Foundation

public struct VivoStepLimit: Codable, Sendable, Equatable {
    public enum Source: String, Codable, Sendable {
        case requested
        case configuredMaximum
        case reactionCohort
        case stochasticFiring
        case diffusion
        case experimentBoundary
        case interventionBoundary
        case measurementBoundary
        case checkpointBoundary
        case couplingCadence
    }

    public let source: Source
    public let seconds: Double
    public let subject: String
    public let rationale: String

    public init(source: Source, seconds: Double, subject: String, rationale: String) {
        self.source = source
        self.seconds = seconds
        self.subject = subject
        self.rationale = rationale
    }
}

public struct VivoAdaptiveStepPlan: Codable, Sendable, Equatable {
    public let requestedSeconds: Double
    public let selectedSeconds: Double
    public let limitingSource: VivoStepLimit.Source
    public let limits: [VivoStepLimit]

    public var wasReduced: Bool { selectedSeconds < requestedSeconds }

    public init(
        requestedSeconds: Double,
        selectedSeconds: Double,
        limitingSource: VivoStepLimit.Source,
        limits: [VivoStepLimit]
    ) {
        self.requestedSeconds = requestedSeconds
        self.selectedSeconds = selectedSeconds
        self.limitingSource = limitingSource
        self.limits = limits
    }
}

public struct VivoAdaptiveStepPlanner: Sendable {
    public struct Options: Codable, Sendable, Equatable {
        public var diffusionSafetyFactor: Double
        public var deterministicCohortSafetyFactor: Double
        public var stochasticMaximumExpectedFirings: Double
        public var boundaryTolerance: Double

        public init(
            diffusionSafetyFactor: Double = 0.8,
            deterministicCohortSafetyFactor: Double = 0.9,
            stochasticMaximumExpectedFirings: Double = 8,
            boundaryTolerance: Double = 1e-12
        ) {
            self.diffusionSafetyFactor = diffusionSafetyFactor
            self.deterministicCohortSafetyFactor = deterministicCohortSafetyFactor
            self.stochasticMaximumExpectedFirings = stochasticMaximumExpectedFirings
            self.boundaryTolerance = boundaryTolerance
        }
    }

    private let options: Options

    public init(options: Options = .init()) throws {
        guard options.diffusionSafetyFactor.isFinite,
              options.diffusionSafetyFactor > 0,
              options.diffusionSafetyFactor <= 1,
              options.deterministicCohortSafetyFactor.isFinite,
              options.deterministicCohortSafetyFactor > 0,
              options.deterministicCohortSafetyFactor <= 1,
              options.stochasticMaximumExpectedFirings.isFinite,
              options.stochasticMaximumExpectedFirings > 0,
              options.boundaryTolerance.isFinite,
              options.boundaryTolerance >= 0 else {
            throw VivoRuntimeError.invalidConfiguration("adaptive step planner options are invalid")
        }
        self.options = options
    }

    public func plan(
        requestedSeconds: Double,
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        transport: [VivoSpeciesTransportABI]? = nil,
        exactBoundaries: [(source: VivoStepLimit.Source, timeRemaining: Double, subject: String)] = []
    ) throws -> VivoAdaptiveStepPlan {
        guard requestedSeconds.isFinite, requestedSeconds > 0 else {
            throw VivoRuntimeError.invalidConfiguration("adaptive step request must be finite and positive")
        }
        var limits: [VivoStepLimit] = [
            .init(
                source: .requested,
                seconds: requestedSeconds,
                subject: "request",
                rationale: "caller-requested step"
            ),
            .init(
                source: .configuredMaximum,
                seconds: Double(configuration.maximumTimeStep),
                subject: "runtime",
                rationale: "configured maximumTimeStep"
            )
        ]

        for cohort in try pack.cohortMetadata() {
            if cohort.maximumStableStep.isFinite, cohort.maximumStableStep > 0 {
                let value = Double(cohort.maximumStableStep) * options.deterministicCohortSafetyFactor
                limits.append(.init(
                    source: .reactionCohort,
                    seconds: value,
                    subject: "cohort[\(cohort.index)]",
                    rationale: "compiler-estimated stable step with deterministic safety factor"
                ))
            }
            if configuration.fidelity.rawValue >= VivoFidelity.stochastic.rawValue,
               cohort.stiffnessEstimate.isFinite,
               cohort.stiffnessEstimate > 0 {
                let value = options.stochasticMaximumExpectedFirings / Double(cohort.stiffnessEstimate)
                limits.append(.init(
                    source: .stochasticFiring,
                    seconds: value,
                    subject: "cohort[\(cohort.index)]",
                    rationale: "bounded expected reaction firings per tau-leap"
                ))
            }
        }

        if configuration.fidelity.rawValue >= VivoFidelity.spatial.rawValue,
           let grid = configuration.spatialGrid,
           let transport {
            guard transport.count == Int(pack.runtimeContract.speciesCount) else {
                throw VivoRuntimeError.invalidConfiguration("transport table shape is incompatible with ProgramPack")
            }
            let inverseSpacingSquared = 1.0 / pow(Double(grid.spacingX), 2) +
                                        1.0 / pow(Double(grid.spacingY), 2) +
                                        1.0 / pow(Double(grid.spacingZ), 2)
            for (index, coefficients) in transport.enumerated() {
                let diffusion = Double(coefficients.diffusion)
                guard diffusion.isFinite, diffusion >= 0 else {
                    throw VivoRuntimeError.invalidConfiguration("species transport contains invalid diffusion")
                }
                if diffusion > 0 {
                    let value = options.diffusionSafetyFactor / (2.0 * diffusion * inverseSpacingSquared)
                    limits.append(.init(
                        source: .diffusion,
                        seconds: value,
                        subject: "species[\(index)]",
                        rationale: "explicit 3D diffusion stability bound"
                    ))
                }
            }
        }

        for boundary in exactBoundaries {
            guard boundary.timeRemaining.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("exact step boundary is non-finite")
            }
            if boundary.timeRemaining > options.boundaryTolerance {
                limits.append(.init(
                    source: boundary.source,
                    seconds: boundary.timeRemaining,
                    subject: boundary.subject,
                    rationale: "step terminates exactly at a scheduled boundary"
                ))
            }
        }

        let viable = limits.filter { $0.seconds.isFinite && $0.seconds > 0 }
        guard let limiting = viable.min(by: {
            if abs($0.seconds - $1.seconds) > options.boundaryTolerance { return $0.seconds < $1.seconds }
            return sourceRank($0.source) < sourceRank($1.source)
        }) else {
            throw VivoRuntimeError.invalidConfiguration("adaptive step planner produced no viable limit")
        }
        let selected = max(limiting.seconds, Double(configuration.minimumTimeStep))
        guard selected <= Double(configuration.maximumTimeStep) else {
            throw VivoRuntimeError.invalidConfiguration("adaptive step exceeds configured maximum")
        }
        return VivoAdaptiveStepPlan(
            requestedSeconds: requestedSeconds,
            selectedSeconds: selected,
            limitingSource: limiting.source,
            limits: limits.sorted {
                if $0.seconds != $1.seconds { return $0.seconds < $1.seconds }
                if $0.source != $1.source { return sourceRank($0.source) < sourceRank($1.source) }
                return $0.subject < $1.subject
            }
        )
    }

    private func sourceRank(_ source: VivoStepLimit.Source) -> Int {
        switch source {
        case .permanentShutdown: return 0
        case .interventionBoundary: return 1
        case .measurementBoundary: return 2
        case .checkpointBoundary: return 3
        case .couplingCadence: return 4
        case .experimentBoundary: return 5
        case .diffusion: return 6
        case .stochasticFiring: return 7
        case .reactionCohort: return 8
        case .configuredMaximum: return 9
        case .requested: return 10
        }
    }
}

private extension VivoStepLimit.Source {
    static var permanentShutdown: VivoStepLimit.Source { .experimentBoundary }
}
