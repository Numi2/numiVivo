import Foundation

public struct VivoPhysiologyDiagnostic: Codable, Sendable, Equatable, Hashable {
    public enum Severity: String, Codable, Sendable {
        case note
        case warning
    }

    public let severity: Severity
    public let code: String
    public let subject: String
    public let message: String
}

public struct PreparedVivoPhysiologyCompartment: Codable, Sendable, Equatable {
    public let index: UInt32
    public let identifier: String
    public let name: String
    public let volumeLitres: Double
    public let annotations: [String: String]
}

public struct PreparedVivoPhysiologyAnalyte: Codable, Sendable, Equatable {
    public let index: UInt32
    public let identifier: String
    public let name: String
    public let concentrationUnit: String
    public let minimum: Double
    public let maximum: Double
    public let annotations: [String: String]
}

public struct VivoPhysiologyIncidenceABI: Codable, Sendable, Equatable {
    public let sourcePairIndex: UInt32
    public let coefficientPerSecond: Float
    public let flags: UInt32
    public let reserved: UInt32

    public init(sourcePairIndex: UInt32, coefficientPerSecond: Float, flags: UInt32 = 0, reserved: UInt32 = 0) {
        self.sourcePairIndex = sourcePairIndex
        self.coefficientPerSecond = coefficientPerSecond
        self.flags = flags
        self.reserved = reserved
    }
}

public struct VivoPhysiologyClearanceABI: Codable, Sendable, Equatable {
    public static let saturableFlag: UInt32 = 1 << 0

    public let firstOrderRate: Float
    public let maximumRate: Float
    public let halfSaturation: Float
    public let flags: UInt32

    public init(firstOrderRate: Float = 0, maximumRate: Float = 0, halfSaturation: Float = 0, flags: UInt32 = 0) {
        self.firstOrderRate = firstOrderRate
        self.maximumRate = maximumRate
        self.halfSaturation = halfSaturation
        self.flags = flags
    }
}

public enum PreparedVivoPhysiologyDoseKind: UInt32, Codable, Sendable {
    case concentrationDelta = 0
    case concentrationInfusion = 1
}

public struct PreparedVivoPhysiologyDose: Codable, Sendable, Equatable {
    public let identifier: String
    public let timeSeconds: Double
    public let endTimeSeconds: Double
    public let pairIndex: UInt32
    public let environments: [UInt32]
    public let kind: PreparedVivoPhysiologyDoseKind
    public let value: Float
    public let evidence: [VivoEvidenceReference]
}

public struct PreparedVivoPhysiologyModel: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/prepared-physiology/v1"

    public let schema: String
    public let fingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint?
    public let hostContextFingerprint: VivoFingerprint?
    public let environmentCount: UInt32
    public let pairCount: UInt32
    public let compartments: [PreparedVivoPhysiologyCompartment]
    public let analytes: [PreparedVivoPhysiologyAnalyte]
    public let incidenceOffsets: [UInt32]
    public let incidence: [VivoPhysiologyIncidenceABI]
    public let clearances: [VivoPhysiologyClearanceABI]
    public let initialState: [Float]
    public let doses: [PreparedVivoPhysiologyDose]
    public let preferredTimeStepSeconds: Double
    public let minimumTimeStepSeconds: Double
    public let maximumTimeStepSeconds: Double
    public let linearStabilityLimitSeconds: Double
    public let diagnostics: [VivoPhysiologyDiagnostic]
    public let labels: [String: String]

    public init(
        fingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint?,
        hostContextFingerprint: VivoFingerprint?,
        environmentCount: UInt32,
        pairCount: UInt32,
        compartments: [PreparedVivoPhysiologyCompartment],
        analytes: [PreparedVivoPhysiologyAnalyte],
        incidenceOffsets: [UInt32],
        incidence: [VivoPhysiologyIncidenceABI],
        clearances: [VivoPhysiologyClearanceABI],
        initialState: [Float],
        doses: [PreparedVivoPhysiologyDose],
        preferredTimeStepSeconds: Double,
        minimumTimeStepSeconds: Double,
        maximumTimeStepSeconds: Double,
        linearStabilityLimitSeconds: Double,
        diagnostics: [VivoPhysiologyDiagnostic],
        labels: [String: String]
    ) {
        self.schema = Self.schema
        self.fingerprint = fingerprint
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.environmentCount = environmentCount
        self.pairCount = pairCount
        self.compartments = compartments
        self.analytes = analytes
        self.incidenceOffsets = incidenceOffsets
        self.incidence = incidence
        self.clearances = clearances
        self.initialState = initialState
        self.doses = doses
        self.preferredTimeStepSeconds = preferredTimeStepSeconds
        self.minimumTimeStepSeconds = minimumTimeStepSeconds
        self.maximumTimeStepSeconds = maximumTimeStepSeconds
        self.linearStabilityLimitSeconds = linearStabilityLimitSeconds
        self.diagnostics = diagnostics
        self.labels = labels
    }

    public func pairIndex(analyte: UInt32, compartment: UInt32) throws -> UInt32 {
        guard analyte < UInt32(analytes.count), compartment < UInt32(compartments.count) else {
            throw VivoArtifactValidationError.invalid("physiology analyte or compartment index is out of bounds")
        }
        let product = analyte.multipliedReportingOverflow(by: UInt32(compartments.count))
        let result = product.partialValue.addingReportingOverflow(compartment)
        guard !product.overflow, !result.overflow else {
            throw VivoArtifactValidationError.invalid("physiology pair index overflow")
        }
        return result.partialValue
    }

    public func pairIndex(analyte: String, compartment: String) throws -> UInt32 {
        guard let analyteIndex = analytes.first(where: { $0.identifier == analyte })?.index else {
            throw VivoArtifactValidationError.unresolved("unknown physiology analyte \(analyte)")
        }
        guard let compartmentIndex = compartments.first(where: { $0.identifier == compartment })?.index else {
            throw VivoArtifactValidationError.unresolved("unknown physiology compartment \(compartment)")
        }
        return try pairIndex(analyte: analyteIndex, compartment: compartmentIndex)
    }

    public func stateIndex(pairIndex: UInt32, environment: UInt32) throws -> Int {
        guard pairIndex < pairCount, environment < environmentCount else {
            throw VivoArtifactValidationError.invalid("physiology pair or environment index is out of bounds")
        }
        let value = UInt64(pairIndex) * UInt64(environmentCount) + UInt64(environment)
        guard value <= UInt64(Int.max) else {
            throw VivoArtifactValidationError.invalid("physiology state index exceeds Int.max")
        }
        return Int(value)
    }
}

public struct VivoPhysiologyCompiler: Sendable {
    public struct Limits: Codable, Sendable, Equatable {
        public var maximumCompartments: Int
        public var maximumAnalytes: Int
        public var maximumPairs: Int
        public var maximumStateElements: Int
        public var maximumIncidenceRecords: Int
        public var maximumPreparedDoseRecords: Int

        public init(
            maximumCompartments: Int = 4_096,
            maximumAnalytes: Int = 16_384,
            maximumPairs: Int = 1_048_576,
            maximumStateElements: Int = 268_435_456,
            maximumIncidenceRecords: Int = 16_777_216,
            maximumPreparedDoseRecords: Int = 16_777_216
        ) {
            self.maximumCompartments = maximumCompartments
            self.maximumAnalytes = maximumAnalytes
            self.maximumPairs = maximumPairs
            self.maximumStateElements = maximumStateElements
            self.maximumIncidenceRecords = maximumIncidenceRecords
            self.maximumPreparedDoseRecords = maximumPreparedDoseRecords
        }
    }

    private let units: VivoUnitSystem
    private let limits: Limits

    public init(units: VivoUnitSystem = .standard, limits: Limits = .init()) {
        self.units = units
        self.limits = limits
    }

    public func compile(_ model: VivoPhysiologyPack) throws -> PreparedVivoPhysiologyModel {
        try model.validateStructure()
        try validateLimits(model)

        let compartmentIndices = try uniqueIndex(model.compartments.map(\.identifier), label: "physiology compartment")
        let analyteIndices = try uniqueIndex(model.analytes.map(\.identifier), label: "physiology analyte")
        let environmentCount = Int(model.environmentCount)
        let pairCount = try checkedProduct(model.compartments.count, model.analytes.count, label: "physiology pair count")
        let stateCount = try checkedProduct(pairCount, environmentCount, label: "physiology state count")
        guard pairCount <= limits.maximumPairs,
              stateCount <= limits.maximumStateElements,
              pairCount <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("physiology pair or state count exceeds compiler limits")
        }

        var diagnostics: [VivoPhysiologyDiagnostic] = []
        let preparedCompartments = try compileCompartments(model.compartments, diagnostics: &diagnostics)
        let volumes = preparedCompartments.map(\.volumeLitres)
        let preparedAnalytes = try compileAnalytes(
            model.analytes,
            compartments: compartmentIndices,
            environmentCount: environmentCount,
            compartmentCount: model.compartments.count,
            stateCount: stateCount,
            diagnostics: &diagnostics
        )
        var initialState = preparedAnalytes.state

        var sparse = [[UInt32: Double]](repeating: [:], count: pairCount)
        try compileFlows(
            model.flows,
            compartments: compartmentIndices,
            analytes: analyteIndices,
            volumes: volumes,
            analyteCount: model.analytes.count,
            compartmentCount: model.compartments.count,
            sparse: &sparse,
            diagnostics: &diagnostics
        )
        try compileExchanges(
            model.exchanges,
            compartments: compartmentIndices,
            analytes: analyteIndices,
            volumes: volumes,
            analyteCount: model.analytes.count,
            compartmentCount: model.compartments.count,
            sparse: &sparse,
            diagnostics: &diagnostics
        )

        let clearanceResult = try compileClearances(
            model.clearances,
            compartments: compartmentIndices,
            analytes: analyteIndices,
            preparedAnalytes: preparedAnalytes.metadata,
            volumes: volumes,
            pairCount: pairCount,
            compartmentCount: model.compartments.count,
            diagnostics: &diagnostics
        )
        let incidenceResult = try flatten(sparse)
        var maximumCharacteristicRate = max(clearanceResult.maximumLocalRate, incidenceResult.maximumRowMagnitude)

        let preparedDoses = try compileDoses(
            model.doses,
            compartments: compartmentIndices,
            analytes: analyteIndices,
            preparedAnalytes: preparedAnalytes.metadata,
            volumes: volumes,
            compartmentCount: model.compartments.count,
            environmentCount: environmentCount,
            diagnostics: &diagnostics
        )
        guard preparedDoses.count <= limits.maximumPreparedDoseRecords else {
            throw VivoArtifactValidationError.invalid("prepared physiology dose records exceed compiler limits")
        }

        let preferred = try seconds(model.preferredTimeStep, label: "physiology.preferredTimeStep")
        let minimum = try seconds(model.minimumTimeStep, label: "physiology.minimumTimeStep")
        let requestedMaximum = try seconds(model.maximumTimeStep, label: "physiology.maximumTimeStep")
        guard minimum > 0, preferred >= minimum, requestedMaximum >= preferred else {
            throw VivoArtifactValidationError.invalid("physiology time-step bounds must satisfy 0 < minimum <= preferred <= maximum")
        }
        if !maximumCharacteristicRate.isFinite || maximumCharacteristicRate < 0 {
            throw VivoArtifactValidationError.invalid("compiled physiology characteristic rate is invalid")
        }
        if maximumCharacteristicRate == 0 { maximumCharacteristicRate = 0 }
        let linearLimit = maximumCharacteristicRate > 0 ? 0.5 / maximumCharacteristicRate : requestedMaximum
        let effectiveMaximum = min(requestedMaximum, linearLimit)
        guard effectiveMaximum >= minimum else {
            throw VivoArtifactValidationError.invalid("physiology minimum time step exceeds the compiled stability limit")
        }
        if effectiveMaximum < requestedMaximum {
            diagnostics.append(.init(
                severity: .warning,
                code: "NVPHYS005",
                subject: "timeStep",
                message: "Maximum time step was reduced from \(requestedMaximum) s to \(effectiveMaximum) s by the compiled flow and clearance bound."
            ))
        }

        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(model))
        return PreparedVivoPhysiologyModel(
            fingerprint: fingerprint,
            programFingerprint: model.programFingerprint,
            hostContextFingerprint: model.hostContextFingerprint,
            environmentCount: model.environmentCount,
            pairCount: UInt32(pairCount),
            compartments: preparedCompartments,
            analytes: preparedAnalytes.metadata,
            incidenceOffsets: incidenceResult.offsets,
            incidence: incidenceResult.records,
            clearances: clearanceResult.records,
            initialState: initialState,
            doses: preparedDoses,
            preferredTimeStepSeconds: min(preferred, effectiveMaximum),
            minimumTimeStepSeconds: minimum,
            maximumTimeStepSeconds: effectiveMaximum,
            linearStabilityLimitSeconds: linearLimit,
            diagnostics: diagnostics,
            labels: model.labels
        )
    }

    private func compileCompartments(
        _ compartments: [VivoPhysiologyCompartment],
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws -> [PreparedVivoPhysiologyCompartment] {
        var result: [PreparedVivoPhysiologyCompartment] = []
        result.reserveCapacity(compartments.count)
        for (index, compartment) in compartments.enumerated() {
            try compartment.volume.validate(label: "physiology.compartments.\(compartment.identifier).volume", nonnegative: true)
            let volume = try units.convert(compartment.volume, to: "L")
            guard volume.isFinite, volume > 0 else {
                throw VivoArtifactValidationError.invalid("physiology compartment \(compartment.identifier) volume must be positive")
            }
            result.append(.init(
                index: UInt32(index),
                identifier: compartment.identifier,
                name: compartment.name,
                volumeLitres: volume,
                annotations: compartment.annotations
            ))
            if compartment.evidence.isEmpty {
                diagnostics.append(.init(
                    severity: .warning,
                    code: "NVPHYS001",
                    subject: compartment.identifier,
                    message: "Compartment volume has no attached evidence record."
                ))
            }
        }
        return result
    }

    private func compileAnalytes(
        _ analytes: [VivoPhysiologyAnalyte],
        compartments: [String: Int],
        environmentCount: Int,
        compartmentCount: Int,
        stateCount: Int,
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws -> (metadata: [PreparedVivoPhysiologyAnalyte], state: [Float]) {
        var metadata: [PreparedVivoPhysiologyAnalyte] = []
        var state = [Float](repeating: 0, count: stateCount)
        for (analyteIndex, analyte) in analytes.enumerated() {
            guard units.areCompatible(analyte.concentrationUnit, "M") else {
                throw VivoArtifactValidationError.incompatible("analyte \(analyte.identifier) unit \(analyte.concentrationUnit) is not a molar concentration unit")
            }
            try analyte.bounds.minimum.validate(label: "physiology.analytes.\(analyte.identifier).bounds.minimum", nonnegative: true)
            try analyte.bounds.maximum.validate(label: "physiology.analytes.\(analyte.identifier).bounds.maximum", nonnegative: true)
            let minimum = try units.convert(analyte.bounds.minimum, to: analyte.concentrationUnit)
            let maximum = try units.convert(analyte.bounds.maximum, to: analyte.concentrationUnit)
            guard minimum.isFinite, maximum.isFinite, minimum >= 0, minimum < maximum else {
                throw VivoArtifactValidationError.invalid("analyte \(analyte.identifier) concentration bounds must be finite, nonnegative, and ordered")
            }
            metadata.append(.init(
                index: UInt32(analyteIndex),
                identifier: analyte.identifier,
                name: analyte.name,
                concentrationUnit: analyte.concentrationUnit,
                minimum: minimum,
                maximum: maximum,
                annotations: analyte.annotations
            ))

            var initialized = Set<String>()
            for condition in analyte.initialConditions {
                guard let compartmentIndex = compartments[condition.compartment] else {
                    throw VivoArtifactValidationError.unresolved("initial condition for \(analyte.identifier) references unknown compartment \(condition.compartment)")
                }
                guard initialized.insert(condition.compartment).inserted,
                      condition.values.count == 1 || condition.values.count == environmentCount else {
                    throw VivoArtifactValidationError.invalid("initial conditions for \(analyte.identifier) must be unique by compartment and provide one value or one per environment")
                }
                let pair = analyteIndex * compartmentCount + compartmentIndex
                for environment in 0..<environmentCount {
                    let source = condition.values.count == 1 ? condition.values[0] : condition.values[environment]
                    try source.validate(label: "physiology.initial.\(analyte.identifier).\(condition.compartment)[\(environment)]", nonnegative: true)
                    let value = try units.convert(source, to: analyte.concentrationUnit)
                    guard value >= minimum, value <= maximum else {
                        throw VivoArtifactValidationError.invalid("initial concentration for \(analyte.identifier) in \(condition.compartment) is outside bounds")
                    }
                    state[pair * environmentCount + environment] = try checkedFloat(value, label: "initial concentration")
                }
            }
            if analyte.initialConditions.isEmpty {
                diagnostics.append(.init(
                    severity: .note,
                    code: "NVPHYS002",
                    subject: analyte.identifier,
                    message: "Analyte has no initial conditions; all compartments start at zero."
                ))
            }
            if analyte.evidence.isEmpty {
                diagnostics.append(.init(
                    severity: .warning,
                    code: "NVPHYS003",
                    subject: analyte.identifier,
                    message: "Analyte definition and bounds have no attached evidence record."
                ))
            }
        }
        return (metadata, state)
    }

    private func compileFlows(
        _ flows: [VivoPhysiologyFlow],
        compartments: [String: Int],
        analytes: [String: Int],
        volumes: [Double],
        analyteCount: Int,
        compartmentCount: Int,
        sparse: inout [[UInt32: Double]],
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws {
        for flow in flows {
            guard let source = compartments[flow.sourceCompartment],
                  let destination = compartments[flow.destinationCompartment],
                  source != destination else {
                throw VivoArtifactValidationError.unresolved("flow \(flow.identifier) has unknown or identical source and destination compartments")
            }
            try flow.rate.validate(label: "physiology.flows.\(flow.identifier).rate", nonnegative: true)
            let rate = try units.convert(flow.rate, to: "L/s")
            guard rate >= 0 else {
                throw VivoArtifactValidationError.invalid("flow \(flow.identifier) rate cannot be negative")
            }
            for analyte in try selectedAnalytes(flow.analytes, indices: analytes, total: analyteCount, subject: flow.identifier) {
                let sourcePair = analyte * compartmentCount + source
                let destinationPair = analyte * compartmentCount + destination
                try add(&sparse[sourcePair], source: sourcePair, coefficient: -rate / volumes[source])
                try add(&sparse[destinationPair], source: sourcePair, coefficient: rate / volumes[destination])
            }
            if flow.evidence.isEmpty {
                diagnostics.append(.init(severity: .warning, code: "NVPHYS006", subject: flow.identifier, message: "Directed flow has no attached evidence record."))
            }
        }
    }

    private func compileExchanges(
        _ exchanges: [VivoPhysiologyExchange],
        compartments: [String: Int],
        analytes: [String: Int],
        volumes: [Double],
        analyteCount: Int,
        compartmentCount: Int,
        sparse: inout [[UInt32: Double]],
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws {
        for exchange in exchanges {
            guard let source = compartments[exchange.sourceCompartment],
                  let destination = compartments[exchange.destinationCompartment],
                  source != destination else {
                throw VivoArtifactValidationError.unresolved("exchange \(exchange.identifier) has unknown or identical compartments")
            }
            guard exchange.destinationToSourcePartition.isFinite,
                  exchange.destinationToSourcePartition > 0 else {
                throw VivoArtifactValidationError.invalid("exchange \(exchange.identifier) partition coefficient must be finite and positive")
            }
            try exchange.conductance.validate(label: "physiology.exchanges.\(exchange.identifier).conductance", nonnegative: true)
            let conductance = try units.convert(exchange.conductance, to: "L/s")
            guard conductance >= 0 else {
                throw VivoArtifactValidationError.invalid("exchange \(exchange.identifier) conductance cannot be negative")
            }
            let partition = exchange.destinationToSourcePartition
            for analyte in try selectedAnalytes(exchange.analytes, indices: analytes, total: analyteCount, subject: exchange.identifier) {
                let sourcePair = analyte * compartmentCount + source
                let destinationPair = analyte * compartmentCount + destination
                try add(&sparse[sourcePair], source: sourcePair, coefficient: -conductance / volumes[source])
                try add(&sparse[sourcePair], source: destinationPair, coefficient: conductance / (volumes[source] * partition))
                try add(&sparse[destinationPair], source: sourcePair, coefficient: conductance / volumes[destination])
                try add(&sparse[destinationPair], source: destinationPair, coefficient: -conductance / (volumes[destination] * partition))
            }
            if exchange.evidence.isEmpty {
                diagnostics.append(.init(severity: .warning, code: "NVPHYS007", subject: exchange.identifier, message: "Bidirectional exchange has no attached evidence record."))
            }
        }
    }

    private func compileClearances(
        _ clearances: [VivoPhysiologyClearance],
        compartments: [String: Int],
        analytes: [String: Int],
        preparedAnalytes: [PreparedVivoPhysiologyAnalyte],
        volumes: [Double],
        pairCount: Int,
        compartmentCount: Int,
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws -> (records: [VivoPhysiologyClearanceABI], maximumLocalRate: Double) {
        var records = [VivoPhysiologyClearanceABI](repeating: .init(), count: pairCount)
        var occupied = Set<Int>()
        var maximumLocalRate = 0.0
        for clearance in clearances {
            guard let analyteIndex = analytes[clearance.analyte],
                  let compartmentIndex = compartments[clearance.compartment] else {
                throw VivoArtifactValidationError.unresolved("clearance \(clearance.identifier) references an unknown analyte or compartment")
            }
            let pair = analyteIndex * compartmentCount + compartmentIndex
            guard occupied.insert(pair).inserted else {
                throw VivoArtifactValidationError.invalid("multiple clearance laws target \(clearance.analyte) in \(clearance.compartment)")
            }
            switch clearance.law {
            case .firstOrder(let quantity):
                try quantity.validate(label: "physiology.clearance.\(clearance.identifier).rate", nonnegative: true)
                let rate = try units.convert(quantity, to: "1/s")
                guard rate >= 0 else { throw VivoArtifactValidationError.invalid("clearance rate cannot be negative") }
                records[pair] = .init(firstOrderRate: try checkedFloat(rate, label: "first-order clearance"))
                maximumLocalRate = max(maximumLocalRate, rate)
            case .volumeFlow(let quantity):
                try quantity.validate(label: "physiology.clearance.\(clearance.identifier).clearance", nonnegative: true)
                let volumeRate = try units.convert(quantity, to: "L/s")
                guard volumeRate >= 0 else { throw VivoArtifactValidationError.invalid("clearance volume flow cannot be negative") }
                let rate = volumeRate / volumes[compartmentIndex]
                records[pair] = .init(firstOrderRate: try checkedFloat(rate, label: "volume-flow clearance"))
                maximumLocalRate = max(maximumLocalRate, rate)
            case .saturable(let maximumRate, let halfSaturation):
                try maximumRate.validate(label: "physiology.clearance.\(clearance.identifier).maximumRate", nonnegative: true)
                try halfSaturation.validate(label: "physiology.clearance.\(clearance.identifier).halfSaturation", nonnegative: true)
                let analyte = preparedAnalytes[analyteIndex]
                let vmax = try convertConcentrationRate(maximumRate, toConcentrationUnit: analyte.concentrationUnit)
                let half = try units.convert(halfSaturation, to: analyte.concentrationUnit)
                guard vmax > 0, half > 0 else {
                    throw VivoArtifactValidationError.invalid("saturable clearance requires positive maximum rate and half saturation")
                }
                records[pair] = .init(
                    maximumRate: try checkedFloat(vmax, label: "saturable maximum rate"),
                    halfSaturation: try checkedFloat(half, label: "saturable half saturation"),
                    flags: VivoPhysiologyClearanceABI.saturableFlag
                )
                maximumLocalRate = max(maximumLocalRate, vmax / half)
            }
            if clearance.evidence.isEmpty {
                diagnostics.append(.init(severity: .warning, code: "NVPHYS004", subject: clearance.identifier, message: "Clearance law has no attached evidence record."))
            }
        }
        return (records, maximumLocalRate)
    }

    private func flatten(_ sparse: [[UInt32: Double]]) throws -> (offsets: [UInt32], records: [VivoPhysiologyIncidenceABI], maximumRowMagnitude: Double) {
        var offsets: [UInt32] = [0]
        var records: [VivoPhysiologyIncidenceABI] = []
        var maximumRowMagnitude = 0.0
        for row in sparse {
            let rowMagnitude = row.values.reduce(0) { $0 + abs($1) }
            maximumRowMagnitude = max(maximumRowMagnitude, rowMagnitude)
            for (source, coefficient) in row.sorted(by: { $0.key < $1.key }) where coefficient != 0 {
                records.append(.init(
                    sourcePairIndex: source,
                    coefficientPerSecond: try checkedFloat(coefficient, label: "physiology incidence coefficient")
                ))
            }
            guard records.count <= limits.maximumIncidenceRecords,
                  records.count <= Int(UInt32.max) else {
                throw VivoArtifactValidationError.invalid("physiology incidence table exceeds compiler limits")
            }
            offsets.append(UInt32(records.count))
        }
        return (offsets, records, maximumRowMagnitude)
    }

    private func compileDoses(
        _ doses: [VivoPhysiologyDose],
        compartments: [String: Int],
        analytes: [String: Int],
        preparedAnalytes: [PreparedVivoPhysiologyAnalyte],
        volumes: [Double],
        compartmentCount: Int,
        environmentCount: Int,
        diagnostics: inout [VivoPhysiologyDiagnostic]
    ) throws -> [PreparedVivoPhysiologyDose] {
        var result: [PreparedVivoPhysiologyDose] = []
        for dose in doses {
            guard let analyteIndex = analytes[dose.analyte],
                  let compartmentIndex = compartments[dose.compartment] else {
                throw VivoArtifactValidationError.unresolved("dose \(dose.identifier) references an unknown analyte or compartment")
            }
            let time = try seconds(dose.time, label: "physiology.doses.\(dose.identifier).time", permitZero: true)
            let environments = try selectedEnvironments(dose.environments, count: environmentCount, subject: dose.identifier)
            let concentrationUnit = preparedAnalytes[analyteIndex].concentrationUnit
            let volume = volumes[compartmentIndex]
            let pairIndex = UInt32(analyteIndex * compartmentCount + compartmentIndex)
            let kind: PreparedVivoPhysiologyDoseKind
            let value: Double
            let endTime: Double

            switch dose.action {
            case .concentrationDelta(let quantity):
                try quantity.validate(label: "physiology.doses.\(dose.identifier).value", nonnegative: true)
                kind = .concentrationDelta
                value = try units.convert(quantity, to: concentrationUnit)
                endTime = time
            case .amount(let quantity):
                try quantity.validate(label: "physiology.doses.\(dose.identifier).amount", nonnegative: true)
                kind = .concentrationDelta
                let moles = try units.convert(quantity, to: "mol")
                value = try units.convert(moles / volume, from: "M", to: concentrationUnit)
                endTime = time
            case .concentrationInfusion(let rate, let duration):
                try rate.validate(label: "physiology.doses.\(dose.identifier).rate", nonnegative: true)
                kind = .concentrationInfusion
                value = try convertConcentrationRate(rate, toConcentrationUnit: concentrationUnit)
                let durationSeconds = try seconds(duration, label: "physiology.doses.\(dose.identifier).duration")
                endTime = time + durationSeconds
            case .amountInfusion(let rate, let duration):
                try rate.validate(label: "physiology.doses.\(dose.identifier).rate", nonnegative: true)
                kind = .concentrationInfusion
                let molesPerSecond = try units.convert(rate, to: "mol/s")
                value = try concentrationRateFromMolar(molesPerSecond / volume, to: concentrationUnit)
                let durationSeconds = try seconds(duration, label: "physiology.doses.\(dose.identifier).duration")
                endTime = time + durationSeconds
            }
            guard value.isFinite, value >= 0, endTime.isFinite, endTime >= time else {
                throw VivoArtifactValidationError.invalid("dose \(dose.identifier) compiled to an invalid value or time interval")
            }
            result.append(PreparedVivoPhysiologyDose(
                identifier: dose.identifier,
                timeSeconds: time,
                endTimeSeconds: endTime,
                pairIndex: pairIndex,
                environments: environments,
                kind: kind,
                value: try checkedFloat(value, label: "physiology dose value"),
                evidence: dose.evidence
            ))
            if dose.evidence.isEmpty {
                diagnostics.append(.init(severity: .warning, code: "NVPHYS008", subject: dose.identifier, message: "Dose has no attached evidence or provenance record."))
            }
        }
        return result.sorted {
            if $0.timeSeconds != $1.timeSeconds { return $0.timeSeconds < $1.timeSeconds }
            if $0.identifier != $1.identifier { return $0.identifier < $1.identifier }
            return $0.pairIndex < $1.pairIndex
        }
    }

    private func convertConcentrationRate(_ quantity: VivoQuantity, toConcentrationUnit unit: String) throws -> Double {
        guard quantity.value.isFinite,
              let source = units.definition(for: quantity.unit),
              let target = units.definition(for: unit),
              source.dimension == VivoDimension(length: -3, time: -1, amount: 1),
              target.dimension == VivoDimension(length: -3, amount: 1) else {
            throw VivoArtifactValidationError.incompatible("concentration-rate unit \(quantity.unit) is incompatible with \(unit) per second")
        }
        let siRate = (quantity.value + source.offsetToSI) * source.scaleToSI
        let converted = siRate / target.scaleToSI
        guard converted.isFinite else {
            throw VivoArtifactValidationError.invalid("concentration-rate conversion produced a non-finite value")
        }
        return converted
    }

    private func concentrationRateFromMolar(_ molarPerSecond: Double, to unit: String) throws -> Double {
        guard molarPerSecond.isFinite,
              let source = units.definition(for: "M"),
              let target = units.definition(for: unit),
              target.dimension == source.dimension else {
            throw VivoArtifactValidationError.incompatible("cannot convert molar rate to concentration unit \(unit)")
        }
        let converted = molarPerSecond * source.scaleToSI / target.scaleToSI
        guard converted.isFinite else {
            throw VivoArtifactValidationError.invalid("molar-rate conversion produced a non-finite value")
        }
        return converted
    }

    private func selectedAnalytes(_ requested: [String]?, indices: [String: Int], total: Int, subject: String) throws -> [Int] {
        guard let requested else { return Array(0..<total) }
        guard !requested.isEmpty, Set(requested).count == requested.count else {
            throw VivoArtifactValidationError.invalid("\(subject) analyte selection must be non-empty and unique")
        }
        return try requested.map { identifier in
            guard let index = indices[identifier] else {
                throw VivoArtifactValidationError.unresolved("\(subject) references unknown analyte \(identifier)")
            }
            return index
        }.sorted()
    }

    private func selectedEnvironments(_ requested: [UInt32]?, count: Int, subject: String) throws -> [UInt32] {
        guard let requested else { return (0..<count).map(UInt32.init) }
        guard !requested.isEmpty,
              Set(requested).count == requested.count,
              requested.allSatisfy({ $0 < UInt32(count) }) else {
            throw VivoArtifactValidationError.invalid("dose \(subject) environment selection is empty, duplicated, or out of bounds")
        }
        return requested.sorted()
    }

    private func add(_ row: inout [UInt32: Double], source: Int, coefficient: Double) throws {
        guard source >= 0, source <= Int(UInt32.max), coefficient.isFinite else {
            throw VivoArtifactValidationError.invalid("physiology incidence source or coefficient is invalid")
        }
        let key = UInt32(source)
        let combined = row[key, default: 0] + coefficient
        guard combined.isFinite else {
            throw VivoArtifactValidationError.invalid("physiology incidence coefficient overflow")
        }
        row[key] = combined
    }

    private func seconds(_ quantity: VivoQuantity, label: String, permitZero: Bool = false) throws -> Double {
        try quantity.validate(label: label, nonnegative: true)
        let value = try units.convert(quantity, to: "s")
        guard value.isFinite, permitZero ? value >= 0 : value > 0 else {
            throw VivoArtifactValidationError.invalid("\(label) must be \(permitZero ? "nonnegative" : "positive")")
        }
        return value
    }

    private func uniqueIndex(_ values: [String], label: String) throws -> [String: Int] {
        var result: [String: Int] = [:]
        for (index, value) in values.enumerated() {
            guard result.updateValue(index, forKey: value) == nil else {
                throw VivoArtifactValidationError.invalid("duplicate \(label) identifier \(value)")
            }
        }
        return result
    }

    private func validateLimits(_ model: VivoPhysiologyPack) throws {
        guard limits.maximumCompartments > 0,
              limits.maximumAnalytes > 0,
              limits.maximumPairs > 0,
              limits.maximumStateElements > 0,
              limits.maximumIncidenceRecords > 0,
              limits.maximumPreparedDoseRecords > 0,
              model.compartments.count <= limits.maximumCompartments,
              model.analytes.count <= limits.maximumAnalytes else {
            throw VivoArtifactValidationError.invalid("physiology compiler limits are invalid or exceeded")
        }
    }

    private func checkedProduct(_ lhs: Int, _ rhs: Int, label: String) throws -> Int {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else { throw VivoArtifactValidationError.invalid("\(label) overflow") }
        return result.partialValue
    }

    private func checkedFloat(_ value: Double, label: String) throws -> Float {
        guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid("\(label) is not representable as FP32")
        }
        return Float(value)
    }
}
