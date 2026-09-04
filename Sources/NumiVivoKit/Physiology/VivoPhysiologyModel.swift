import Foundation

public struct VivoPhysiologyCompartment: Codable, Sendable, Equatable {
    public var identifier: String
    public var name: String
    public var volume: VivoQuantity
    public var annotations: [String: String]
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        name: String = "",
        volume: VivoQuantity,
        annotations: [String: String] = [:],
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.name = name
        self.volume = volume
        self.annotations = annotations
        self.evidence = evidence
    }
}

public struct VivoPhysiologyConcentrationBounds: Codable, Sendable, Equatable {
    public var minimum: VivoQuantity
    public var maximum: VivoQuantity

    public init(minimum: VivoQuantity, maximum: VivoQuantity) {
        self.minimum = minimum
        self.maximum = maximum
    }
}

public struct VivoPhysiologyInitialCondition: Codable, Sendable, Equatable {
    public var compartment: String
    public var values: [VivoQuantity]
    public var evidence: [VivoEvidenceReference]

    public init(
        compartment: String,
        values: [VivoQuantity],
        evidence: [VivoEvidenceReference] = []
    ) {
        self.compartment = compartment
        self.values = values
        self.evidence = evidence
    }
}

public struct VivoPhysiologyAnalyte: Codable, Sendable, Equatable {
    public var identifier: String
    public var name: String
    public var concentrationUnit: String
    public var bounds: VivoPhysiologyConcentrationBounds
    public var initialConditions: [VivoPhysiologyInitialCondition]
    public var annotations: [String: String]
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        name: String = "",
        concentrationUnit: String,
        bounds: VivoPhysiologyConcentrationBounds,
        initialConditions: [VivoPhysiologyInitialCondition] = [],
        annotations: [String: String] = [:],
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.name = name
        self.concentrationUnit = concentrationUnit
        self.bounds = bounds
        self.initialConditions = initialConditions
        self.annotations = annotations
        self.evidence = evidence
    }
}

public struct VivoPhysiologyFlow: Codable, Sendable, Equatable {
    public var identifier: String
    public var sourceCompartment: String
    public var destinationCompartment: String
    public var rate: VivoQuantity
    public var analytes: [String]?
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        sourceCompartment: String,
        destinationCompartment: String,
        rate: VivoQuantity,
        analytes: [String]? = nil,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.sourceCompartment = sourceCompartment
        self.destinationCompartment = destinationCompartment
        self.rate = rate
        self.analytes = analytes
        self.evidence = evidence
    }
}

public struct VivoPhysiologyExchange: Codable, Sendable, Equatable {
    public var identifier: String
    public var sourceCompartment: String
    public var destinationCompartment: String
    public var conductance: VivoQuantity
    public var destinationToSourcePartition: Double
    public var analytes: [String]?
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        sourceCompartment: String,
        destinationCompartment: String,
        conductance: VivoQuantity,
        destinationToSourcePartition: Double = 1,
        analytes: [String]? = nil,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.sourceCompartment = sourceCompartment
        self.destinationCompartment = destinationCompartment
        self.conductance = conductance
        self.destinationToSourcePartition = destinationToSourcePartition
        self.analytes = analytes
        self.evidence = evidence
    }
}

public enum VivoPhysiologyClearanceLaw: Codable, Sendable, Equatable {
    case firstOrder(rate: VivoQuantity)
    case volumeFlow(clearance: VivoQuantity)
    case saturable(maximumRate: VivoQuantity, halfSaturation: VivoQuantity)

    private enum CodingKeys: String, CodingKey {
        case kind
        case rate
        case clearance
        case maximumRate
        case halfSaturation
    }

    private enum Kind: String, Codable {
        case firstOrder
        case volumeFlow
        case saturable
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .firstOrder:
            self = .firstOrder(rate: try container.decode(VivoQuantity.self, forKey: .rate))
        case .volumeFlow:
            self = .volumeFlow(clearance: try container.decode(VivoQuantity.self, forKey: .clearance))
        case .saturable:
            self = .saturable(
                maximumRate: try container.decode(VivoQuantity.self, forKey: .maximumRate),
                halfSaturation: try container.decode(VivoQuantity.self, forKey: .halfSaturation)
            )
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .firstOrder(let rate):
            try container.encode(Kind.firstOrder, forKey: .kind)
            try container.encode(rate, forKey: .rate)
        case .volumeFlow(let clearance):
            try container.encode(Kind.volumeFlow, forKey: .kind)
            try container.encode(clearance, forKey: .clearance)
        case .saturable(let maximumRate, let halfSaturation):
            try container.encode(Kind.saturable, forKey: .kind)
            try container.encode(maximumRate, forKey: .maximumRate)
            try container.encode(halfSaturation, forKey: .halfSaturation)
        }
    }
}

public struct VivoPhysiologyClearance: Codable, Sendable, Equatable {
    public var identifier: String
    public var analyte: String
    public var compartment: String
    public var law: VivoPhysiologyClearanceLaw
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        analyte: String,
        compartment: String,
        law: VivoPhysiologyClearanceLaw,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.analyte = analyte
        self.compartment = compartment
        self.law = law
        self.evidence = evidence
    }
}

public enum VivoPhysiologyDoseAction: Codable, Sendable, Equatable {
    case concentrationDelta(VivoQuantity)
    case amount(VivoQuantity)
    case concentrationInfusion(rate: VivoQuantity, duration: VivoQuantity)
    case amountInfusion(rate: VivoQuantity, duration: VivoQuantity)

    private enum CodingKeys: String, CodingKey {
        case kind
        case value
        case rate
        case duration
    }

    private enum Kind: String, Codable {
        case concentrationDelta
        case amount
        case concentrationInfusion
        case amountInfusion
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .concentrationDelta:
            self = .concentrationDelta(try container.decode(VivoQuantity.self, forKey: .value))
        case .amount:
            self = .amount(try container.decode(VivoQuantity.self, forKey: .value))
        case .concentrationInfusion:
            self = .concentrationInfusion(
                rate: try container.decode(VivoQuantity.self, forKey: .rate),
                duration: try container.decode(VivoQuantity.self, forKey: .duration)
            )
        case .amountInfusion:
            self = .amountInfusion(
                rate: try container.decode(VivoQuantity.self, forKey: .rate),
                duration: try container.decode(VivoQuantity.self, forKey: .duration)
            )
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .concentrationDelta(let value):
            try container.encode(Kind.concentrationDelta, forKey: .kind)
            try container.encode(value, forKey: .value)
        case .amount(let value):
            try container.encode(Kind.amount, forKey: .kind)
            try container.encode(value, forKey: .value)
        case .concentrationInfusion(let rate, let duration):
            try container.encode(Kind.concentrationInfusion, forKey: .kind)
            try container.encode(rate, forKey: .rate)
            try container.encode(duration, forKey: .duration)
        case .amountInfusion(let rate, let duration):
            try container.encode(Kind.amountInfusion, forKey: .kind)
            try container.encode(rate, forKey: .rate)
            try container.encode(duration, forKey: .duration)
        }
    }
}

public struct VivoPhysiologyDose: Codable, Sendable, Equatable {
    public var identifier: String
    public var time: VivoQuantity
    public var analyte: String
    public var compartment: String
    public var environments: [UInt32]?
    public var action: VivoPhysiologyDoseAction
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        time: VivoQuantity,
        analyte: String,
        compartment: String,
        environments: [UInt32]? = nil,
        action: VivoPhysiologyDoseAction,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.time = time
        self.analyte = analyte
        self.compartment = compartment
        self.environments = environments
        self.action = action
        self.evidence = evidence
    }
}

public struct VivoPhysiologyPack: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/physiology-pack/v1"

    public var schema: String
    public var programFingerprint: VivoFingerprint?
    public var hostContextFingerprint: VivoFingerprint?
    public var environmentCount: UInt32
    public var compartments: [VivoPhysiologyCompartment]
    public var analytes: [VivoPhysiologyAnalyte]
    public var flows: [VivoPhysiologyFlow]
    public var exchanges: [VivoPhysiologyExchange]
    public var clearances: [VivoPhysiologyClearance]
    public var doses: [VivoPhysiologyDose]
    public var preferredTimeStep: VivoQuantity
    public var minimumTimeStep: VivoQuantity
    public var maximumTimeStep: VivoQuantity
    public var labels: [String: String]

    public init(
        programFingerprint: VivoFingerprint? = nil,
        hostContextFingerprint: VivoFingerprint? = nil,
        environmentCount: UInt32,
        compartments: [VivoPhysiologyCompartment],
        analytes: [VivoPhysiologyAnalyte],
        flows: [VivoPhysiologyFlow] = [],
        exchanges: [VivoPhysiologyExchange] = [],
        clearances: [VivoPhysiologyClearance] = [],
        doses: [VivoPhysiologyDose] = [],
        preferredTimeStep: VivoQuantity,
        minimumTimeStep: VivoQuantity,
        maximumTimeStep: VivoQuantity,
        labels: [String: String] = [:]
    ) {
        self.schema = Self.schema
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.environmentCount = environmentCount
        self.compartments = compartments
        self.analytes = analytes
        self.flows = flows
        self.exchanges = exchanges
        self.clearances = clearances
        self.doses = doses
        self.preferredTimeStep = preferredTimeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.labels = labels
    }

    public func validateStructure() throws {
        guard schema == Self.schema,
              environmentCount > 0,
              !compartments.isEmpty,
              !analytes.isEmpty else {
            throw VivoArtifactValidationError.invalid("physiology pack schema, environment count, compartments, or analytes are invalid")
        }
        try requireUnique(compartments.map(\.identifier), label: "physiology compartment")
        try requireUnique(analytes.map(\.identifier), label: "physiology analyte")
        try requireUnique(flows.map(\.identifier), label: "physiology flow")
        try requireUnique(exchanges.map(\.identifier), label: "physiology exchange")
        try requireUnique(clearances.map(\.identifier), label: "physiology clearance")
        try requireUnique(doses.map(\.identifier), label: "physiology dose")
        try preferredTimeStep.validate(label: "physiology.preferredTimeStep", nonnegative: true)
        try minimumTimeStep.validate(label: "physiology.minimumTimeStep", nonnegative: true)
        try maximumTimeStep.validate(label: "physiology.maximumTimeStep", nonnegative: true)
    }

    private func requireUnique(_ identifiers: [String], label: String) throws {
        guard identifiers.allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              Set(identifiers).count == identifiers.count else {
            throw VivoArtifactValidationError.invalid("\(label) identifiers must be non-empty and unique")
        }
    }
}
