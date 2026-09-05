import Foundation

public enum VivoIntervention: Codable, Sendable, Equatable {
    case setSignal(species: String, value: VivoQuantity, lanes: VivoLaneSelection)
    case addSignal(species: String, value: VivoQuantity, lanes: VivoLaneSelection)
    case setParameter(parameter: String, value: VivoQuantity, environments: [UInt32]?)
    case setTransport(species: String, definition: VivoSpeciesTransport)
    case reversibleShutdown(reason: String)
    case permanentShutdown(reason: String)

    private enum CodingKeys: String, CodingKey {
        case kind
        case species
        case parameter
        case value
        case lanes
        case environments
        case definition
        case reason
    }

    private enum Kind: String, Codable {
        case setSignal
        case addSignal
        case setParameter
        case setTransport
        case reversibleShutdown
        case permanentShutdown
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .setSignal:
            self = .setSignal(
                species: try container.decode(String.self, forKey: .species),
                value: try container.decode(VivoQuantity.self, forKey: .value),
                lanes: try container.decode(VivoLaneSelection.self, forKey: .lanes)
            )
        case .addSignal:
            self = .addSignal(
                species: try container.decode(String.self, forKey: .species),
                value: try container.decode(VivoQuantity.self, forKey: .value),
                lanes: try container.decode(VivoLaneSelection.self, forKey: .lanes)
            )
        case .setParameter:
            self = .setParameter(
                parameter: try container.decode(String.self, forKey: .parameter),
                value: try container.decode(VivoQuantity.self, forKey: .value),
                environments: try container.decodeIfPresent([UInt32].self, forKey: .environments)
            )
        case .setTransport:
            self = .setTransport(
                species: try container.decode(String.self, forKey: .species),
                definition: try container.decode(VivoSpeciesTransport.self, forKey: .definition)
            )
        case .reversibleShutdown:
            self = .reversibleShutdown(reason: try container.decode(String.self, forKey: .reason))
        case .permanentShutdown:
            self = .permanentShutdown(reason: try container.decode(String.self, forKey: .reason))
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .setSignal(let species, let value, let lanes):
            try container.encode(Kind.setSignal, forKey: .kind)
            try container.encode(species, forKey: .species)
            try container.encode(value, forKey: .value)
            try container.encode(lanes, forKey: .lanes)
        case .addSignal(let species, let value, let lanes):
            try container.encode(Kind.addSignal, forKey: .kind)
            try container.encode(species, forKey: .species)
            try container.encode(value, forKey: .value)
            try container.encode(lanes, forKey: .lanes)
        case .setParameter(let parameter, let value, let environments):
            try container.encode(Kind.setParameter, forKey: .kind)
            try container.encode(parameter, forKey: .parameter)
            try container.encode(value, forKey: .value)
            try container.encodeIfPresent(environments, forKey: .environments)
        case .setTransport(let species, let definition):
            try container.encode(Kind.setTransport, forKey: .kind)
            try container.encode(species, forKey: .species)
            try container.encode(definition, forKey: .definition)
        case .reversibleShutdown(let reason):
            try container.encode(Kind.reversibleShutdown, forKey: .kind)
            try container.encode(reason, forKey: .reason)
        case .permanentShutdown(let reason):
            try container.encode(Kind.permanentShutdown, forKey: .kind)
            try container.encode(reason, forKey: .reason)
        }
    }
}

public struct VivoTimedIntervention: Codable, Sendable, Equatable {
    public var identifier: String
    public var time: VivoQuantity
    public var intervention: VivoIntervention
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        time: VivoQuantity,
        intervention: VivoIntervention,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.time = time
        self.intervention = intervention
        self.evidence = evidence
    }
}

public struct VivoMeasurement: Codable, Sendable, Equatable {
    public enum Storage: String, Codable, Sendable {
        case fullPrecision
        case reducedPrecision
        case summaryOnly
        case eventOnly
    }

    public var identifier: String
    public var species: String
    public var lanes: VivoLaneSelection
    public var aggregation: VivoCouplingAggregation
    public var cadence: VivoQuantity
    public var start: VivoQuantity
    public var end: VivoQuantity?
    public var outputUnit: String
    public var storage: Storage

    public init(
        identifier: String,
        species: String,
        lanes: VivoLaneSelection = .all,
        aggregation: VivoCouplingAggregation = .direct,
        cadence: VivoQuantity,
        start: VivoQuantity = .init(value: 0, unit: "s"),
        end: VivoQuantity? = nil,
        outputUnit: String,
        storage: Storage = .fullPrecision
    ) {
        self.identifier = identifier
        self.species = species
        self.lanes = lanes
        self.aggregation = aggregation
        self.cadence = cadence
        self.start = start
        self.end = end
        self.outputUnit = outputUnit
        self.storage = storage
    }
}

public struct VivoExperimentStopCondition: Codable, Sendable, Equatable {
    public enum Response: String, Codable, Sendable {
        case finishReplicate
        case rejectReplicate
        case reversibleShutdown
        case permanentShutdown
    }

    public var identifier: String
    public var monitorIdentifier: String
    public var response: Response
    public var required: Bool

    public init(
        identifier: String,
        monitorIdentifier: String,
        response: Response,
        required: Bool = true
    ) {
        self.identifier = identifier
        self.monitorIdentifier = monitorIdentifier
        self.response = response
        self.required = required
    }
}

public struct VivoExperimentPack: Codable, Sendable, Equatable {
    public var programFingerprint: VivoFingerprint
    public var hostContextFingerprint: VivoFingerprint
    public var couplingFingerprint: VivoFingerprint?
    public var fidelity: VivoFidelity
    public var duration: VivoQuantity
    public var preferredTimeStep: VivoQuantity
    public var minimumTimeStep: VivoQuantity
    public var maximumTimeStep: VivoQuantity
    public var environmentCount: UInt32
    public var replicateCount: UInt32
    public var baseSeed: VivoRuntimeSeed
    public var interventions: [VivoTimedIntervention]
    public var measurements: [VivoMeasurement]
    public var stopConditions: [VivoExperimentStopCondition]
    public var checkpointCadence: VivoQuantity?
    public var labels: [String: String]

    public init(
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint,
        couplingFingerprint: VivoFingerprint? = nil,
        fidelity: VivoFidelity,
        duration: VivoQuantity,
        preferredTimeStep: VivoQuantity,
        minimumTimeStep: VivoQuantity,
        maximumTimeStep: VivoQuantity,
        environmentCount: UInt32,
        replicateCount: UInt32,
        baseSeed: VivoRuntimeSeed,
        interventions: [VivoTimedIntervention] = [],
        measurements: [VivoMeasurement],
        stopConditions: [VivoExperimentStopCondition] = [],
        checkpointCadence: VivoQuantity? = nil,
        labels: [String: String] = [:]
    ) {
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.fidelity = fidelity
        self.duration = duration
        self.preferredTimeStep = preferredTimeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.environmentCount = environmentCount
        self.replicateCount = replicateCount
        self.baseSeed = baseSeed
        self.interventions = interventions
        self.measurements = measurements
        self.stopConditions = stopConditions
        self.checkpointCadence = checkpointCadence
        self.labels = labels
    }
}

public struct PreparedVivoIntervention: Codable, Sendable, Equatable {
    public enum Operation: Codable, Sendable, Equatable {
        case coupling([VivoCouplingUpdate])
        case parameter(index: UInt32, environments: [UInt32], value: Float)
        case transport(index: UInt32, value: VivoSpeciesTransportABI)
        case reversibleShutdown(reason: String)
        case permanentShutdown(reason: String)

        private enum CodingKeys: String, CodingKey {
            case kind
            case coupling
            case index
            case environments
            case value
            case transport
            case reason
        }

        private enum Kind: String, Codable {
            case coupling
            case parameter
            case transport
            case reversibleShutdown
            case permanentShutdown
        }

        public init(from decoder: Decoder) throws {
            let container = try decoder.container(keyedBy: CodingKeys.self)
            switch try container.decode(Kind.self, forKey: .kind) {
            case .coupling:
                self = .coupling(try container.decode([VivoCouplingUpdate].self, forKey: .coupling))
            case .parameter:
                self = .parameter(
                    index: try container.decode(UInt32.self, forKey: .index),
                    environments: try container.decode([UInt32].self, forKey: .environments),
                    value: try container.decode(Float.self, forKey: .value)
                )
            case .transport:
                self = .transport(
                    index: try container.decode(UInt32.self, forKey: .index),
                    value: try container.decode(VivoSpeciesTransportABI.self, forKey: .transport)
                )
            case .reversibleShutdown:
                self = .reversibleShutdown(reason: try container.decode(String.self, forKey: .reason))
            case .permanentShutdown:
                self = .permanentShutdown(reason: try container.decode(String.self, forKey: .reason))
            }
        }

        public func encode(to encoder: Encoder) throws {
            var container = encoder.container(keyedBy: CodingKeys.self)
            switch self {
            case .coupling(let coupling):
                try container.encode(Kind.coupling, forKey: .kind)
                try container.encode(coupling, forKey: .coupling)
            case .parameter(let index, let environments, let value):
                try container.encode(Kind.parameter, forKey: .kind)
                try container.encode(index, forKey: .index)
                try container.encode(environments, forKey: .environments)
                try container.encode(value, forKey: .value)
            case .transport(let index, let value):
                try container.encode(Kind.transport, forKey: .kind)
                try container.encode(index, forKey: .index)
                try container.encode(value, forKey: .transport)
            case .reversibleShutdown(let reason):
                try container.encode(Kind.reversibleShutdown, forKey: .kind)
                try container.encode(reason, forKey: .reason)
            case .permanentShutdown(let reason):
                try container.encode(Kind.permanentShutdown, forKey: .kind)
                try container.encode(reason, forKey: .reason)
            }
        }
    }

    public let identifier: String
    public let timeSeconds: Double
    public let operation: Operation
}

public struct PreparedVivoMeasurement: Codable, Sendable, Equatable {
    public let identifier: String
    public let speciesIndex: UInt32
    public let lanes: [UInt32]
    public let aggregation: VivoCouplingAggregation
    public let cadenceSeconds: Double
    public let startSeconds: Double
    public let endSeconds: Double
    public let transform: VivoUnitTransform
    public let outputUnit: String
    public let storage: VivoMeasurement.Storage
}

public struct PreparedVivoExperiment: Codable, Sendable, Equatable {
    public let fingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint
    public let couplingFingerprint: VivoFingerprint?
    public let fidelity: VivoFidelity
    public let durationSeconds: Double
    public let preferredTimeStepSeconds: Double
    public let minimumTimeStepSeconds: Double
    public let maximumTimeStepSeconds: Double
    public let environmentCount: UInt32
    public let replicateCount: UInt32
    public let replicateSeeds: [VivoRuntimeSeed]
    public let interventions: [PreparedVivoIntervention]
    public let measurements: [PreparedVivoMeasurement]
    public let stopConditions: [VivoExperimentStopCondition]
    public let checkpointCadenceSeconds: Double?
}

public struct VivoExperimentCompiler: Sendable {
    private let units: VivoUnitSystem

    public init(units: VivoUnitSystem = .standard) {
        self.units = units
    }

    public func compile(
        _ experiment: VivoExperimentPack,
        programPack: VivoProgramPack,
        hostContext: PreparedVivoHostContext,
        coupling: PreparedVivoCouplingPlan? = nil,
        configuration: VivoRuntimeConfiguration
    ) throws -> PreparedVivoExperiment {
        guard experiment.programFingerprint == programPack.header.contentFingerprint else {
            throw VivoArtifactValidationError.incompatible("experiment references a different ProgramPack")
        }
        guard experiment.hostContextFingerprint == hostContext.contextFingerprint else {
            throw VivoArtifactValidationError.incompatible("experiment references a different host context")
        }
        guard experiment.environmentCount == configuration.environmentCount,
              experiment.environmentCount == hostContext.environmentCount else {
            throw VivoArtifactValidationError.incompatible("experiment environmentCount disagrees with runtime or host context")
        }
        guard experiment.replicateCount > 0 else {
            throw VivoArtifactValidationError.invalid("experiment replicateCount must be positive")
        }
        if let expected = experiment.couplingFingerprint {
            guard coupling?.fingerprint == expected else {
                throw VivoArtifactValidationError.incompatible("experiment coupling fingerprint is unresolved")
            }
        }
        guard experiment.fidelity == configuration.fidelity else {
            throw VivoArtifactValidationError.incompatible("experiment fidelity disagrees with runtime configuration")
        }

        let duration = try positiveSeconds(experiment.duration, label: "duration")
        let preferred = try positiveSeconds(experiment.preferredTimeStep, label: "preferredTimeStep")
        let minimum = try positiveSeconds(experiment.minimumTimeStep, label: "minimumTimeStep")
        let maximum = try positiveSeconds(experiment.maximumTimeStep, label: "maximumTimeStep")
        guard minimum <= preferred, preferred <= maximum else {
            throw VivoArtifactValidationError.invalid("experiment time steps must satisfy minimum <= preferred <= maximum")
        }
        let checkpoint = try experiment.checkpointCadence.map {
            try positiveSeconds($0, label: "checkpointCadence")
        }

        let species = try programPack.speciesMetadata()
        let parameters = try programPack.parameterMetadata()
        let speciesIndices = try uniqueIndex(species.map(\.identifier), label: "species")
        let parameterIndices = try uniqueIndex(parameters.map(\.identifier), label: "parameter")
        let monitors = try programPack.monitorMetadata()
        let monitorIdentifiers = Set(monitors.map(\.identifier))

        var interventionIdentifiers = Set<String>()
        var preparedInterventions: [PreparedVivoIntervention] = []
        for intervention in experiment.interventions {
            guard !intervention.identifier.isEmpty,
                  interventionIdentifiers.insert(intervention.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("intervention identifiers must be non-empty and unique")
            }
            let time = try nonnegativeSeconds(intervention.time, label: "interventions.\(intervention.identifier).time")
            guard time <= duration else {
                throw VivoArtifactValidationError.invalid("intervention \(intervention.identifier) occurs after experiment duration")
            }
            preparedInterventions.append(.init(
                identifier: intervention.identifier,
                timeSeconds: time,
                operation: try compileIntervention(
                    intervention.intervention,
                    species: species,
                    speciesIndices: speciesIndices,
                    parameters: parameters,
                    parameterIndices: parameterIndices,
                    configuration: configuration
                )
            ))
        }
        preparedInterventions.sort {
            if $0.timeSeconds != $1.timeSeconds { return $0.timeSeconds < $1.timeSeconds }
            return $0.identifier < $1.identifier
        }

        guard !experiment.measurements.isEmpty else {
            throw VivoArtifactValidationError.invalid("experiment requires at least one measurement")
        }
        var measurementIdentifiers = Set<String>()
        var preparedMeasurements: [PreparedVivoMeasurement] = []
        for measurement in experiment.measurements {
            guard !measurement.identifier.isEmpty,
                  measurementIdentifiers.insert(measurement.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("measurement identifiers must be non-empty and unique")
            }
            guard let index = speciesIndices[measurement.species] else {
                throw VivoArtifactValidationError.unresolved("measurement species \(measurement.species) is absent from ProgramPack")
            }
            guard units.definition(for: measurement.outputUnit) != nil,
                  units.areCompatible(species[index].unit, measurement.outputUnit) else {
                throw VivoArtifactValidationError.incompatible(
                    "measurement \(measurement.identifier) output unit is incompatible with species unit"
                )
            }
            let start = try nonnegativeSeconds(measurement.start, label: "measurements.\(measurement.identifier).start")
            let end = try measurement.end.map {
                try nonnegativeSeconds($0, label: "measurements.\(measurement.identifier).end")
            } ?? duration
            let cadence = try positiveSeconds(measurement.cadence, label: "measurements.\(measurement.identifier).cadence")
            guard start <= end, end <= duration else {
                throw VivoArtifactValidationError.invalid("measurement \(measurement.identifier) time window is invalid")
            }
            let transform = try unitTransform(from: species[index].unit, to: measurement.outputUnit)
            preparedMeasurements.append(.init(
                identifier: measurement.identifier,
                speciesIndex: UInt32(index),
                lanes: try measurement.lanes.resolve(configuration: configuration),
                aggregation: measurement.aggregation,
                cadenceSeconds: cadence,
                startSeconds: start,
                endSeconds: end,
                transform: transform,
                outputUnit: measurement.outputUnit,
                storage: measurement.storage
            ))
        }
        preparedMeasurements.sort { $0.identifier < $1.identifier }

        var stopIdentifiers = Set<String>()
        for condition in experiment.stopConditions {
            guard !condition.identifier.isEmpty, stopIdentifiers.insert(condition.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("stop-condition identifiers must be non-empty and unique")
            }
            guard monitorIdentifiers.contains(condition.monitorIdentifier) || !condition.required else {
                throw VivoArtifactValidationError.unresolved(
                    "required stop condition \(condition.identifier) references missing monitor \(condition.monitorIdentifier)"
                )
            }
        }

        let replicateSeeds = try deriveSeeds(base: experiment.baseSeed, count: experiment.replicateCount)
        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment))
        return PreparedVivoExperiment(
            fingerprint: fingerprint,
            programFingerprint: experiment.programFingerprint,
            hostContextFingerprint: experiment.hostContextFingerprint,
            couplingFingerprint: experiment.couplingFingerprint,
            fidelity: experiment.fidelity,
            durationSeconds: duration,
            preferredTimeStepSeconds: preferred,
            minimumTimeStepSeconds: minimum,
            maximumTimeStepSeconds: maximum,
            environmentCount: experiment.environmentCount,
            replicateCount: experiment.replicateCount,
            replicateSeeds: replicateSeeds,
            interventions: preparedInterventions,
            measurements: preparedMeasurements,
            stopConditions: experiment.stopConditions,
            checkpointCadenceSeconds: checkpoint
        )
    }

    private func compileIntervention(
        _ intervention: VivoIntervention,
        species: [VivoProgramPack.SpeciesMetadata],
        speciesIndices: [String: Int],
        parameters: [VivoProgramPack.ParameterMetadata],
        parameterIndices: [String: Int],
        configuration: VivoRuntimeConfiguration
    ) throws -> PreparedVivoIntervention.Operation {
        switch intervention {
        case .setSignal(let identifier, let quantity, let selection),
             .addSignal(let identifier, let quantity, let selection):
            guard let index = speciesIndices[identifier] else {
                throw VivoArtifactValidationError.unresolved("intervention species \(identifier) is absent from ProgramPack")
            }
            guard species[index].isExternallyOwned || species[index].isInput else {
                throw VivoArtifactValidationError.incompatible("intervention species \(identifier) is internally owned")
            }
            let converted = try units.convert(quantity, to: species[index].unit)
            guard converted >= Double(species[index].minimum), converted <= Double(species[index].maximum) else {
                throw VivoArtifactValidationError.invalid("intervention value for \(identifier) is outside species bounds")
            }
            let value = try checkedFloat(converted, label: "intervention.\(identifier)")
            let mode: VivoCouplingMode
            switch intervention {
            case .setSignal: mode = .replace
            case .addSignal: mode = .add
            default: preconditionFailure("unreachable")
            }
            return .coupling(try selection.resolve(configuration: configuration).map {
                VivoCouplingUpdate(speciesIndex: UInt32(index), laneIndex: $0, mode: mode, value: value)
            })
        case .setParameter(let identifier, let quantity, let selectedEnvironments):
            guard let index = parameterIndices[identifier] else {
                throw VivoArtifactValidationError.unresolved("intervention parameter \(identifier) is absent from ProgramPack")
            }
            let metadata = parameters[index]
            let converted = try units.convert(quantity, to: metadata.unit)
            guard converted >= metadata.minimum, converted <= metadata.maximum else {
                throw VivoArtifactValidationError.invalid("intervention value for parameter \(identifier) is outside bounds")
            }
            let environments = selectedEnvironments ?? Array(0..<configuration.environmentCount)
            guard !environments.isEmpty,
                  Set(environments).count == environments.count,
                  environments.allSatisfy({ $0 < configuration.environmentCount }) else {
                throw VivoArtifactValidationError.invalid("parameter intervention environments are empty, duplicated, or out of bounds")
            }
            return .parameter(
                index: UInt32(index),
                environments: environments.sorted(),
                value: try checkedFloat(converted, label: "intervention.parameter.\(identifier)")
            )
        case .setTransport(let identifier, let definition):
            guard identifier == definition.species else {
                throw VivoArtifactValidationError.invalid("transport intervention species and definition species disagree")
            }
            guard let index = speciesIndices[identifier] else {
                throw VivoArtifactValidationError.unresolved("transport intervention species \(identifier) is absent from ProgramPack")
            }
            try definition.validate()
            return .transport(
                index: UInt32(index),
                value: .init(
                    diffusion: try checkedFloat(units.convert(definition.diffusion, to: "m2/s"), label: "transport.diffusion"),
                    membranePermeability: try definition.membranePermeability.map {
                        try checkedFloat(units.convert($0, to: "m/s"), label: "transport.membranePermeability")
                    } ?? 0,
                    decayRate: try definition.extracellularDecayRate.map {
                        try checkedFloat(units.convert($0, to: "1/s"), label: "transport.decayRate")
                    } ?? 0,
                    flags: 0
                )
            )
        case .reversibleShutdown(let reason):
            guard !reason.isEmpty else {
                throw VivoArtifactValidationError.invalid("reversible shutdown intervention requires a reason")
            }
            return .reversibleShutdown(reason: reason)
        case .permanentShutdown(let reason):
            guard !reason.isEmpty else {
                throw VivoArtifactValidationError.invalid("permanent shutdown intervention requires a reason")
            }
            return .permanentShutdown(reason: reason)
        }
    }

    private func positiveSeconds(_ quantity: VivoQuantity, label: String) throws -> Double {
        let value = try nonnegativeSeconds(quantity, label: label)
        guard value > 0 else {
            throw VivoArtifactValidationError.invalid("\(label) must be positive")
        }
        return value
    }

    private func nonnegativeSeconds(_ quantity: VivoQuantity, label: String) throws -> Double {
        try quantity.validate(label: label, nonnegative: true)
        let seconds = try units.convert(quantity, to: "s")
        guard seconds >= 0 else {
            throw VivoArtifactValidationError.invalid("\(label) cannot be negative")
        }
        return seconds
    }

    private func unitTransform(from source: String, to destination: String) throws -> VivoUnitTransform {
        let zero = try units.convert(0, from: source, to: destination)
        let one = try units.convert(1, from: source, to: destination)
        return .init(scale: one - zero, offset: zero)
    }

    private func uniqueIndex(_ identifiers: [String], label: String) throws -> [String: Int] {
        var result: [String: Int] = [:]
        for (index, identifier) in identifiers.enumerated() {
            guard result.updateValue(index, forKey: identifier) == nil else {
                throw VivoArtifactValidationError.invalid("duplicate ProgramPack \(label) \(identifier)")
            }
        }
        return result
    }

    private func checkedFloat(_ value: Double, label: String) throws -> Float {
        guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid("\(label) cannot be represented as finite FP32")
        }
        return Float(value)
    }

    private func deriveSeeds(base: VivoRuntimeSeed, count: UInt32) throws -> [VivoRuntimeSeed] {
        guard count <= UInt32(Int.max) else {
            throw VivoArtifactValidationError.invalid("replicateCount exceeds addressable memory")
        }
        return (0..<count).map { index in
            VivoRuntimeSeed(
                low: splitMix64(base.low &+ UInt64(index) &* 0x9E3779B97F4A7C15),
                high: splitMix64(base.high &+ UInt64(index) &* 0xD2B74407B1CE6E93)
            )
        }
    }

    private func splitMix64(_ source: UInt64) -> UInt64 {
        var value = source &+ 0x9E3779B97F4A7C15
        value = (value ^ (value >> 30)) &* 0xBF58476D1CE4E5B9
        value = (value ^ (value >> 27)) &* 0x94D049BB133111EB
        return value ^ (value >> 31)
    }
}

public struct VivoExperimentCampaignAxis: Codable, Sendable, Equatable {
    public var identifier: String
    public var parameter: String
    public var values: [VivoQuantity]

    public init(identifier: String, parameter: String, values: [VivoQuantity]) {
        self.identifier = identifier
        self.parameter = parameter
        self.values = values
    }
}

public struct VivoExperimentCampaignManifest: Codable, Sendable, Equatable {
    public var baseExperiment: VivoExperimentPack
    public var axes: [VivoExperimentCampaignAxis]
    public var maximumRuns: UInt64

    public init(baseExperiment: VivoExperimentPack, axes: [VivoExperimentCampaignAxis], maximumRuns: UInt64) {
        self.baseExperiment = baseExperiment
        self.axes = axes
        self.maximumRuns = maximumRuns
    }

    public func expandedExperiments() throws -> [VivoExperimentPack] {
        guard maximumRuns > 0 else {
            throw VivoArtifactValidationError.invalid("campaign maximumRuns must be positive")
        }
        var identifiers = Set<String>()
        for axis in axes {
            guard !axis.identifier.isEmpty,
                  identifiers.insert(axis.identifier).inserted,
                  !axis.parameter.isEmpty,
                  !axis.values.isEmpty else {
                throw VivoArtifactValidationError.invalid("campaign axes require unique identifiers, a parameter, and values")
            }
        }
        let total = axes.reduce(UInt64(1)) { partial, axis in
            let product = partial.multipliedReportingOverflow(by: UInt64(axis.values.count))
            return product.overflow ? UInt64.max : product.partialValue
        }
        guard total <= maximumRuns, total <= UInt64(Int.max) else {
            throw VivoArtifactValidationError.invalid("campaign expansion exceeds maximumRuns or addressable memory")
        }
        guard !axes.isEmpty else { return [baseExperiment] }

        var results: [VivoExperimentPack] = []
        results.reserveCapacity(Int(total))
        func expand(axisIndex: Int, selections: [(VivoExperimentCampaignAxis, VivoQuantity)]) {
            if axisIndex == axes.count {
                var experiment = baseExperiment
                let suffix = selections.map { "\($0.0.identifier)=\($0.1.value)\($0.1.unit)" }.joined(separator: ",")
                experiment.labels["campaign.coordinates"] = suffix
                for (offset, selection) in selections.enumerated() {
                    experiment.interventions.append(.init(
                        identifier: "campaign.\(selection.0.identifier).\(offset)",
                        time: .init(value: 0, unit: "s"),
                        intervention: .setParameter(
                            parameter: selection.0.parameter,
                            value: selection.1,
                            environments: nil
                        )
                    ))
                }
                results.append(experiment)
                return
            }
            let axis = axes[axisIndex]
            for value in axis.values {
                expand(axisIndex: axisIndex + 1, selections: selections + [(axis, value)])
            }
        }
        expand(axisIndex: 0, selections: [])
        return results
    }
}
