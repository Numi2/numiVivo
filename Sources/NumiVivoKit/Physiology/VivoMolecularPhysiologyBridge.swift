import Foundation

public enum VivoBridgeBoundPolicy: String, Codable, Sendable {
    case reject
    case clamp
}

public enum VivoBridgeIndexMapping: String, Codable, Sendable {
    case oneToOne
    case broadcastSource
    case broadcastDestination
    case cartesian
}

public struct VivoBridgeTransfer: Codable, Sendable, Equatable {
    public var gain: Double
    public var offset: VivoQuantity?
    public var minimum: VivoQuantity?
    public var maximum: VivoQuantity?
    public var boundPolicy: VivoBridgeBoundPolicy

    public init(
        gain: Double = 1,
        offset: VivoQuantity? = nil,
        minimum: VivoQuantity? = nil,
        maximum: VivoQuantity? = nil,
        boundPolicy: VivoBridgeBoundPolicy = .reject
    ) {
        self.gain = gain
        self.offset = offset
        self.minimum = minimum
        self.maximum = maximum
        self.boundPolicy = boundPolicy
    }
}

public struct VivoPhysiologyToMolecularLink: Codable, Sendable, Equatable {
    public var identifier: String
    public var analyte: String
    public var compartment: String
    public var physiologyEnvironments: [UInt32]?
    public var molecularSpecies: String
    public var molecularLanes: [UInt32]?
    public var mapping: VivoBridgeIndexMapping
    public var molecularMode: VivoCouplingMode
    public var transfer: VivoBridgeTransfer
    public var required: Bool

    public init(
        identifier: String,
        analyte: String,
        compartment: String,
        physiologyEnvironments: [UInt32]? = nil,
        molecularSpecies: String,
        molecularLanes: [UInt32]? = nil,
        mapping: VivoBridgeIndexMapping = .oneToOne,
        molecularMode: VivoCouplingMode = .replace,
        transfer: VivoBridgeTransfer = .init(),
        required: Bool = true
    ) {
        self.identifier = identifier
        self.analyte = analyte
        self.compartment = compartment
        self.physiologyEnvironments = physiologyEnvironments
        self.molecularSpecies = molecularSpecies
        self.molecularLanes = molecularLanes
        self.mapping = mapping
        self.molecularMode = molecularMode
        self.transfer = transfer
        self.required = required
    }
}

public struct VivoMolecularToPhysiologyLink: Codable, Sendable, Equatable {
    public var identifier: String
    public var molecularSpecies: String
    public var molecularLanes: [UInt32]?
    public var analyte: String
    public var compartment: String
    public var physiologyEnvironments: [UInt32]?
    public var mapping: VivoBridgeIndexMapping
    public var physiologyMode: VivoPhysiologyUpdateMode
    public var transfer: VivoBridgeTransfer
    public var required: Bool

    public init(
        identifier: String,
        molecularSpecies: String,
        molecularLanes: [UInt32]? = nil,
        analyte: String,
        compartment: String,
        physiologyEnvironments: [UInt32]? = nil,
        mapping: VivoBridgeIndexMapping = .oneToOne,
        physiologyMode: VivoPhysiologyUpdateMode = .rate,
        transfer: VivoBridgeTransfer = .init(),
        required: Bool = true
    ) {
        self.identifier = identifier
        self.molecularSpecies = molecularSpecies
        self.molecularLanes = molecularLanes
        self.analyte = analyte
        self.compartment = compartment
        self.physiologyEnvironments = physiologyEnvironments
        self.mapping = mapping
        self.physiologyMode = physiologyMode
        self.transfer = transfer
        self.required = required
    }
}

public struct VivoMolecularPhysiologyPolicy: Codable, Sendable, Equatable {
    public var maximumStepNegotiations: UInt32
    public var maximumFixedPointIterations: UInt32
    public var convergenceTolerance: Double
    public var relaxation: Double
    public var requireConvergence: Bool
    public var timeSynchronizationTolerance: Double

    public init(
        maximumStepNegotiations: UInt32 = 12,
        maximumFixedPointIterations: UInt32 = 6,
        convergenceTolerance: Double = 1e-5,
        relaxation: Double = 0.7,
        requireConvergence: Bool = true,
        timeSynchronizationTolerance: Double = 1e-8
    ) {
        self.maximumStepNegotiations = maximumStepNegotiations
        self.maximumFixedPointIterations = maximumFixedPointIterations
        self.convergenceTolerance = convergenceTolerance
        self.relaxation = relaxation
        self.requireConvergence = requireConvergence
        self.timeSynchronizationTolerance = timeSynchronizationTolerance
    }

    public func validate() throws {
        guard maximumStepNegotiations > 0,
              maximumFixedPointIterations > 0,
              convergenceTolerance.isFinite,
              convergenceTolerance >= 0,
              relaxation.isFinite,
              relaxation > 0,
              relaxation <= 1,
              timeSynchronizationTolerance.isFinite,
              timeSynchronizationTolerance >= 0 else {
            throw VivoArtifactValidationError.invalid(
                "molecular-physiology coupling policy is invalid"
            )
        }
    }
}

public struct VivoMolecularPhysiologyCouplingPack: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-physiology-coupling/v1"

    public var schema: String
    public var programFingerprint: VivoFingerprint
    public var physiologyFingerprint: VivoFingerprint
    public var physiologyToMolecular: [VivoPhysiologyToMolecularLink]
    public var molecularToPhysiology: [VivoMolecularToPhysiologyLink]
    public var policy: VivoMolecularPhysiologyPolicy
    public var labels: [String: String]

    public init(
        programFingerprint: VivoFingerprint,
        physiologyFingerprint: VivoFingerprint,
        physiologyToMolecular: [VivoPhysiologyToMolecularLink],
        molecularToPhysiology: [VivoMolecularToPhysiologyLink] = [],
        policy: VivoMolecularPhysiologyPolicy = .init(),
        labels: [String: String] = [:]
    ) {
        self.schema = Self.schema
        self.programFingerprint = programFingerprint
        self.physiologyFingerprint = physiologyFingerprint
        self.physiologyToMolecular = physiologyToMolecular
        self.molecularToPhysiology = molecularToPhysiology
        self.policy = policy
        self.labels = labels
    }
}

public struct PreparedVivoBridgeTransfer: Codable, Sendable, Equatable {
    public let scale: Double
    public let offset: Double
    public let minimum: Double?
    public let maximum: Double?
    public let boundPolicy: VivoBridgeBoundPolicy

    public func apply(_ source: Float, subject: String) throws -> Float {
        guard source.isFinite else {
            throw VivoArtifactValidationError.invalid(
                "bridge source \(subject) is non-finite"
            )
        }
        var result = scale * Double(source) + offset
        guard result.isFinite else {
            throw VivoArtifactValidationError.invalid(
                "bridge transfer \(subject) produced a non-finite value"
            )
        }
        if let minimum, result < minimum {
            guard boundPolicy == .clamp else {
                throw VivoArtifactValidationError.invalid(
                    "bridge transfer \(subject) is below its destination minimum"
                )
            }
            result = minimum
        }
        if let maximum, result > maximum {
            guard boundPolicy == .clamp else {
                throw VivoArtifactValidationError.invalid(
                    "bridge transfer \(subject) exceeds its destination maximum"
                )
            }
            result = maximum
        }
        guard abs(result) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid(
                "bridge transfer \(subject) is not FP32 representable"
            )
        }
        return Float(result)
    }
}

public struct PreparedVivoPhysiologyToMolecularLink: Codable, Sendable, Equatable {
    public let identifier: String
    public let pairIndex: UInt32
    public let environmentIndex: UInt32
    public let molecularSpeciesIndex: UInt32
    public let molecularLaneIndex: UInt32
    public let molecularMode: VivoCouplingMode
    public let transfer: PreparedVivoBridgeTransfer
    public let required: Bool
}

public struct PreparedVivoMolecularToPhysiologyLink: Codable, Sendable, Equatable {
    public let identifier: String
    public let molecularSpeciesIndex: UInt32
    public let molecularLaneIndex: UInt32
    public let pairIndex: UInt32
    public let environmentIndex: UInt32
    public let physiologyMode: VivoPhysiologyUpdateMode
    public let transfer: PreparedVivoBridgeTransfer
    public let required: Bool
}

public struct PreparedVivoMolecularPhysiologyBridge: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/prepared-molecular-physiology-coupling/v1"

    public let schema: String
    public let fingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let physiologyFingerprint: VivoFingerprint
    public let physiologyToMolecular: [PreparedVivoPhysiologyToMolecularLink]
    public let molecularToPhysiology: [PreparedVivoMolecularToPhysiologyLink]
    public let policy: VivoMolecularPhysiologyPolicy
    public let labels: [String: String]

    public init(
        fingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint,
        physiologyFingerprint: VivoFingerprint,
        physiologyToMolecular: [PreparedVivoPhysiologyToMolecularLink],
        molecularToPhysiology: [PreparedVivoMolecularToPhysiologyLink],
        policy: VivoMolecularPhysiologyPolicy,
        labels: [String: String]
    ) {
        self.schema = Self.schema
        self.fingerprint = fingerprint
        self.programFingerprint = programFingerprint
        self.physiologyFingerprint = physiologyFingerprint
        self.physiologyToMolecular = physiologyToMolecular
        self.molecularToPhysiology = molecularToPhysiology
        self.policy = policy
        self.labels = labels
    }
}

public struct VivoMolecularPhysiologyCouplingCompiler: Sendable {
    public struct Limits: Codable, Sendable, Equatable {
        public var maximumExpandedLinks: Int

        public init(maximumExpandedLinks: Int = 1_048_576) {
            self.maximumExpandedLinks = maximumExpandedLinks
        }
    }

    private let units: VivoUnitSystem
    private let limits: Limits

    public init(
        units: VivoUnitSystem = .standard,
        limits: Limits = .init()
    ) {
        self.units = units
        self.limits = limits
    }

    public func compile(
        _ coupling: VivoMolecularPhysiologyCouplingPack,
        programPack: VivoProgramPack,
        physiology: PreparedVivoPhysiologyModel,
        molecularLaneCount: UInt32
    ) throws -> PreparedVivoMolecularPhysiologyBridge {
        guard coupling.schema == VivoMolecularPhysiologyCouplingPack.schema,
              coupling.programFingerprint == programPack.header.contentFingerprint,
              coupling.physiologyFingerprint == physiology.fingerprint,
              molecularLaneCount > 0,
              limits.maximumExpandedLinks > 0 else {
            throw VivoArtifactValidationError.incompatible(
                "molecular-physiology coupling fingerprints, schema, lane count, or limits are invalid"
            )
        }
        try coupling.policy.validate()

        let species = try programPack.speciesMetadata()
        let speciesByIdentifier = Dictionary(
            uniqueKeysWithValues: species.map { ($0.identifier, $0) }
        )
        let analyteByIdentifier = Dictionary(
            uniqueKeysWithValues: physiology.analytes.map { ($0.identifier, $0) }
        )
        var identifiers = Set<String>()
        var exposure: [PreparedVivoPhysiologyToMolecularLink] = []
        var feedback: [PreparedVivoMolecularToPhysiologyLink] = []
        var molecularDestinations = Set<UInt64>()

        for link in coupling.physiologyToMolecular {
            guard !link.identifier.isEmpty,
                  identifiers.insert(link.identifier).inserted,
                  let analyte = analyteByIdentifier[link.analyte],
                  let target = speciesByIdentifier[link.molecularSpecies],
                  target.isExternallyOwned || target.isInput else {
                throw VivoArtifactValidationError.unresolved(
                    "physiology-to-molecular link \(link.identifier) is duplicated, unresolved, or targets an internally owned molecular species"
                )
            }
            let pairIndex = try physiology.pairIndex(
                analyte: link.analyte,
                compartment: link.compartment
            )
            let mappings = try expandedMappings(
                sources: link.physiologyEnvironments,
                sourceCount: physiology.environmentCount,
                destinations: link.molecularLanes,
                destinationCount: molecularLaneCount,
                mode: link.mapping,
                subject: link.identifier
            )
            let transfer = try compileStateTransfer(
                link.transfer,
                sourceUnit: analyte.concentrationUnit,
                destinationUnit: target.unit,
                subject: link.identifier
            )
            for (environment, lane) in mappings {
                let destinationKey = UInt64(target.index) << 32 | UInt64(lane)
                guard molecularDestinations.insert(destinationKey).inserted else {
                    throw VivoArtifactValidationError.invalid(
                        "multiple bridge links write molecular species \(target.identifier) lane \(lane)"
                    )
                }
                exposure.append(.init(
                    identifier: "\(link.identifier)[\(environment):\(lane)]",
                    pairIndex: pairIndex,
                    environmentIndex: environment,
                    molecularSpeciesIndex: target.index,
                    molecularLaneIndex: lane,
                    molecularMode: link.molecularMode,
                    transfer: transfer,
                    required: link.required
                ))
            }
        }

        for link in coupling.molecularToPhysiology {
            guard !link.identifier.isEmpty,
                  identifiers.insert(link.identifier).inserted,
                  let source = speciesByIdentifier[link.molecularSpecies],
                  let analyte = analyteByIdentifier[link.analyte] else {
                throw VivoArtifactValidationError.unresolved(
                    "molecular-to-physiology link \(link.identifier) is duplicated or unresolved"
                )
            }
            let pairIndex = try physiology.pairIndex(
                analyte: link.analyte,
                compartment: link.compartment
            )
            let mappings = try expandedMappings(
                sources: link.molecularLanes,
                sourceCount: molecularLaneCount,
                destinations: link.physiologyEnvironments,
                destinationCount: physiology.environmentCount,
                mode: link.mapping,
                subject: link.identifier
            )
            let transfer: PreparedVivoBridgeTransfer
            if link.physiologyMode == .rate {
                transfer = try compileRateTransfer(
                    link.transfer,
                    sourceUnit: source.unit,
                    destinationConcentrationUnit: analyte.concentrationUnit,
                    subject: link.identifier
                )
            } else {
                transfer = try compileStateTransfer(
                    link.transfer,
                    sourceUnit: source.unit,
                    destinationUnit: analyte.concentrationUnit,
                    subject: link.identifier
                )
            }
            for (lane, environment) in mappings {
                feedback.append(.init(
                    identifier: "\(link.identifier)[\(lane):\(environment)]",
                    molecularSpeciesIndex: source.index,
                    molecularLaneIndex: lane,
                    pairIndex: pairIndex,
                    environmentIndex: environment,
                    physiologyMode: link.physiologyMode,
                    transfer: transfer,
                    required: link.required
                ))
            }
        }

        guard !exposure.isEmpty || !feedback.isEmpty else {
            throw VivoArtifactValidationError.invalid(
                "molecular-physiology bridge contains no expanded channels"
            )
        }
        guard exposure.count <= limits.maximumExpandedLinks,
              feedback.count <= limits.maximumExpandedLinks else {
            throw VivoArtifactValidationError.invalid(
                "molecular-physiology bridge exceeds expanded-link limits"
            )
        }

        let fingerprint = try VivoCanonicalJSON.fingerprint(
            VivoCanonicalJSON.encode(coupling)
        )
        return PreparedVivoMolecularPhysiologyBridge(
            fingerprint: fingerprint,
            programFingerprint: coupling.programFingerprint,
            physiologyFingerprint: coupling.physiologyFingerprint,
            physiologyToMolecular: exposure.sorted { $0.identifier < $1.identifier },
            molecularToPhysiology: feedback.sorted { $0.identifier < $1.identifier },
            policy: coupling.policy,
            labels: coupling.labels
        )
    }

    private func compileStateTransfer(
        _ transfer: VivoBridgeTransfer,
        sourceUnit: String,
        destinationUnit: String,
        subject: String
    ) throws -> PreparedVivoBridgeTransfer {
        guard transfer.gain.isFinite,
              units.areCompatible(sourceUnit, destinationUnit) else {
            throw VivoArtifactValidationError.incompatible(
                "bridge \(subject) state units \(sourceUnit) and \(destinationUnit) are incompatible"
            )
        }
        let zero = try units.convert(0, from: sourceUnit, to: destinationUnit)
        let one = try units.convert(1, from: sourceUnit, to: destinationUnit)
        return try finishTransfer(
            transfer,
            unitScale: one - zero,
            unitOffset: zero,
            destinationUnit: destinationUnit,
            subject: subject
        )
    }

    private func compileRateTransfer(
        _ transfer: VivoBridgeTransfer,
        sourceUnit: String,
        destinationConcentrationUnit: String,
        subject: String
    ) throws -> PreparedVivoBridgeTransfer {
        guard transfer.gain.isFinite,
              let source = units.definition(for: sourceUnit),
              let destination = units.definition(for: destinationConcentrationUnit),
              source.dimension == VivoDimension(length: -3, time: -1, amount: 1),
              destination.dimension == VivoDimension(length: -3, amount: 1),
              source.offsetToSI == 0,
              destination.offsetToSI == 0 else {
            throw VivoArtifactValidationError.incompatible(
                "bridge \(subject) rate unit \(sourceUnit) is incompatible with \(destinationConcentrationUnit) per second"
            )
        }
        return try finishTransfer(
            transfer,
            unitScale: source.scaleToSI / destination.scaleToSI,
            unitOffset: 0,
            destinationUnit: destinationConcentrationUnit,
            subject: subject
        )
    }

    private func finishTransfer(
        _ transfer: VivoBridgeTransfer,
        unitScale: Double,
        unitOffset: Double,
        destinationUnit: String,
        subject: String
    ) throws -> PreparedVivoBridgeTransfer {
        let offset = try transfer.offset.map {
            try units.convert($0, to: destinationUnit)
        } ?? 0
        let minimum = try transfer.minimum.map {
            try units.convert($0, to: destinationUnit)
        }
        let maximum = try transfer.maximum.map {
            try units.convert($0, to: destinationUnit)
        }
        let scale = transfer.gain * unitScale
        let combinedOffset = transfer.gain * unitOffset + offset
        guard scale.isFinite,
              combinedOffset.isFinite,
              minimum?.isFinite ?? true,
              maximum?.isFinite ?? true,
              minimum == nil || maximum == nil || minimum! <= maximum! else {
            throw VivoArtifactValidationError.invalid(
                "bridge transfer \(subject) contains invalid scale, offset, or bounds"
            )
        }
        return PreparedVivoBridgeTransfer(
            scale: scale,
            offset: combinedOffset,
            minimum: minimum,
            maximum: maximum,
            boundPolicy: transfer.boundPolicy
        )
    }

    private func expandedMappings(
        sources: [UInt32]?,
        sourceCount: UInt32,
        destinations: [UInt32]?,
        destinationCount: UInt32,
        mode: VivoBridgeIndexMapping,
        subject: String
    ) throws -> [(UInt32, UInt32)] {
        let sourceValues = sources ?? Array(0..<sourceCount)
        let destinationValues = destinations ?? Array(0..<destinationCount)
        guard !sourceValues.isEmpty,
              !destinationValues.isEmpty,
              Set(sourceValues).count == sourceValues.count,
              Set(destinationValues).count == destinationValues.count,
              sourceValues.allSatisfy({ $0 < sourceCount }),
              destinationValues.allSatisfy({ $0 < destinationCount }) else {
            throw VivoArtifactValidationError.invalid(
                "bridge \(subject) index sets must be non-empty, unique, and in range"
            )
        }

        switch mode {
        case .oneToOne:
            guard sourceValues.count == destinationValues.count else {
                throw VivoArtifactValidationError.invalid(
                    "bridge \(subject) one-to-one mapping requires equal source and destination counts"
                )
            }
            return Array(zip(sourceValues, destinationValues))
        case .broadcastSource:
            guard sourceValues.count == 1 else {
                throw VivoArtifactValidationError.invalid(
                    "bridge \(subject) broadcastSource requires exactly one source index"
                )
            }
            return destinationValues.map { (sourceValues[0], $0) }
        case .broadcastDestination:
            guard destinationValues.count == 1 else {
                throw VivoArtifactValidationError.invalid(
                    "bridge \(subject) broadcastDestination requires exactly one destination index"
                )
            }
            return sourceValues.map { ($0, destinationValues[0]) }
        case .cartesian:
            let count = sourceValues.count.multipliedReportingOverflow(
                by: destinationValues.count
            )
            guard !count.overflow,
                  count.partialValue <= limits.maximumExpandedLinks else {
                throw VivoArtifactValidationError.invalid(
                    "bridge \(subject) Cartesian mapping exceeds expanded-link limits"
                )
            }
            return sourceValues.flatMap { source in
                destinationValues.map { (source, $0) }
            }
        }
    }
}
