import Foundation

public enum VivoBridgeBoundPolicy: String, Codable, Sendable { case reject, clamp }
public enum VivoBridgeIndexMapping: String, Codable, Sendable {
    case oneToOne, broadcastSource, broadcastDestination, cartesian
}
public struct VivoBridgeTransfer: Codable, Sendable, Equatable {
    public var gain: Double
    public var offset: VivoQuantity?
    public var minimum: VivoQuantity?
    public var maximum: VivoQuantity?
    public var boundPolicy: VivoBridgeBoundPolicy
    public init(gain: Double = 1, offset: VivoQuantity? = nil, minimum: VivoQuantity? = nil,
                maximum: VivoQuantity? = nil, boundPolicy: VivoBridgeBoundPolicy = .reject) {
        self.gain = gain; self.offset = offset; self.minimum = minimum; self.maximum = maximum; self.boundPolicy = boundPolicy
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
    public init(identifier: String, analyte: String, compartment: String, physiologyEnvironments: [UInt32]? = nil,
                molecularSpecies: String, molecularLanes: [UInt32]? = nil, mapping: VivoBridgeIndexMapping = .oneToOne,
                molecularMode: VivoCouplingMode = .replace, transfer: VivoBridgeTransfer = .init(), required: Bool = true) {
        self.identifier = identifier; self.analyte = analyte; self.compartment = compartment
        self.physiologyEnvironments = physiologyEnvironments; self.molecularSpecies = molecularSpecies
        self.molecularLanes = molecularLanes; self.mapping = mapping; self.molecularMode = molecularMode
        self.transfer = transfer; self.required = required
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
    public init(identifier: String, molecularSpecies: String, molecularLanes: [UInt32]? = nil, analyte: String,
                compartment: String, physiologyEnvironments: [UInt32]? = nil, mapping: VivoBridgeIndexMapping = .oneToOne,
                physiologyMode: VivoPhysiologyUpdateMode = .rate, transfer: VivoBridgeTransfer = .init(), required: Bool = true) {
        self.identifier = identifier; self.molecularSpecies = molecularSpecies; self.molecularLanes = molecularLanes
        self.analyte = analyte; self.compartment = compartment; self.physiologyEnvironments = physiologyEnvironments
        self.mapping = mapping; self.physiologyMode = physiologyMode; self.transfer = transfer; self.required = required
    }
}
public struct VivoMolecularPhysiologyPolicy: Codable, Sendable, Equatable {
    public var maximumStepNegotiations: UInt32
    public var maximumFixedPointIterations: UInt32
    public var convergenceTolerance: Double
    public var relaxation: Double
    public var requireConvergence: Bool
    public var timeSynchronizationTolerance: Double
    public init(maximumStepNegotiations: UInt32 = 12, maximumFixedPointIterations: UInt32 = 6,
                convergenceTolerance: Double = 1e-5, relaxation: Double = 0.7, requireConvergence: Bool = true,
                timeSynchronizationTolerance: Double = 1e-8) {
        self.maximumStepNegotiations = maximumStepNegotiations; self.maximumFixedPointIterations = maximumFixedPointIterations
        self.convergenceTolerance = convergenceTolerance; self.relaxation = relaxation
        self.requireConvergence = requireConvergence; self.timeSynchronizationTolerance = timeSynchronizationTolerance
    }
    public func validate() throws {
        guard (1...1024).contains(maximumStepNegotiations), (1...1024).contains(maximumFixedPointIterations),
              convergenceTolerance.isFinite, convergenceTolerance >= 0, convergenceTolerance <= 1,
              relaxation.isFinite, relaxation > 0, relaxation <= 1,
              timeSynchronizationTolerance.isFinite, timeSynchronizationTolerance >= 0 else {
            throw VivoArtifactValidationError.invalid("invalid or unbounded molecular-physiology policy")
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
    public init(programFingerprint: VivoFingerprint, physiologyFingerprint: VivoFingerprint,
                physiologyToMolecular: [VivoPhysiologyToMolecularLink], molecularToPhysiology: [VivoMolecularToPhysiologyLink] = [],
                policy: VivoMolecularPhysiologyPolicy = .init(), labels: [String: String] = [:]) {
        schema = Self.schema; self.programFingerprint = programFingerprint; self.physiologyFingerprint = physiologyFingerprint
        self.physiologyToMolecular = physiologyToMolecular; self.molecularToPhysiology = molecularToPhysiology
        self.policy = policy; self.labels = labels
    }
}
public struct PreparedVivoBridgeTransfer: Codable, Sendable, Equatable {
    public let scale: Double
    public let offset: Double
    public let minimum: Double?
    public let maximum: Double?
    public let boundPolicy: VivoBridgeBoundPolicy
    public func validate() throws {
        guard scale.isFinite, offset.isFinite, minimum?.isFinite ?? true, maximum?.isFinite ?? true,
              minimum == nil || maximum == nil || minimum! <= maximum! else {
            throw VivoArtifactValidationError.invalid("nonfinite or inverted prepared bridge transfer")
        }
        for value in [offset, minimum, maximum].compactMap({ $0 }) {
            guard abs(value) <= Double(Float.greatestFiniteMagnitude), value == 0 || Float(value) != 0 else {
                throw VivoArtifactValidationError.invalid("bridge offset/bounds are not FP32 representable")
            }
        }
    }
    public func apply(_ source: Float, subject: String) throws -> Float {
        try validate()
        guard source.isFinite else { throw VivoArtifactValidationError.invalid("nonfinite bridge source \(subject)") }
        var value = scale * Double(source) + offset
        guard value.isFinite else { throw VivoArtifactValidationError.invalid("nonfinite bridge result \(subject)") }
        if let minimum, value < minimum {
            guard boundPolicy == .clamp else { throw VivoArtifactValidationError.invalid("bridge \(subject) is below its minimum") }
            value = minimum
        }
        if let maximum, value > maximum {
            guard boundPolicy == .clamp else { throw VivoArtifactValidationError.invalid("bridge \(subject) exceeds its maximum") }
            value = maximum
        }
        guard abs(value) <= Double(Float.greatestFiniteMagnitude), value == 0 || Float(value) != 0 else {
            throw VivoArtifactValidationError.invalid("bridge result \(subject) overflows or underflows FP32")
        }
        return Float(value)
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
    public init(fingerprint: VivoFingerprint, programFingerprint: VivoFingerprint, physiologyFingerprint: VivoFingerprint,
                physiologyToMolecular: [PreparedVivoPhysiologyToMolecularLink], molecularToPhysiology: [PreparedVivoMolecularToPhysiologyLink],
                policy: VivoMolecularPhysiologyPolicy, labels: [String: String]) {
        schema = Self.schema; self.fingerprint = fingerprint; self.programFingerprint = programFingerprint
        self.physiologyFingerprint = physiologyFingerprint; self.physiologyToMolecular = physiologyToMolecular
        self.molecularToPhysiology = molecularToPhysiology; self.policy = policy; self.labels = labels
    }
    public func validate(program: VivoProgramPack, physiology: PreparedVivoPhysiologyModel, molecularLaneCount: UInt32) throws {
        try policy.validate()
        try VivoPhysiologyModelValidator.validate(physiology)
        guard schema == Self.schema, programFingerprint == program.header.contentFingerprint,
              physiologyFingerprint == physiology.fingerprint, molecularLaneCount > 0,
              !physiologyToMolecular.isEmpty || !molecularToPhysiology.isEmpty,
              physiologyToMolecular.count + molecularToPhysiology.count <= 1_048_576 else {
            throw VivoArtifactValidationError.incompatible("bridge identity, schema or capacity mismatch")
        }
        let species = try program.speciesMetadata()
        var ids = Set<String>(), molecularDestinations = Set<UInt64>(), physiologyReplacements = Set<UInt64>()
        for link in physiologyToMolecular {
            guard !link.identifier.isEmpty, ids.insert(link.identifier).inserted,
                  link.pairIndex < physiology.pairCount, link.environmentIndex < physiology.environmentCount,
                  Int(link.molecularSpeciesIndex) < species.count, link.molecularLaneIndex < molecularLaneCount,
                  species[Int(link.molecularSpeciesIndex)].isExternallyOwned,
                  molecularDestinations.insert(UInt64(link.molecularSpeciesIndex) << 32 | UInt64(link.molecularLaneIndex)).inserted else {
                throw VivoArtifactValidationError.invalid("invalid, duplicate or internally-owned exposure link \(link.identifier)")
            }
            try link.transfer.validate()
        }
        for link in molecularToPhysiology {
            guard !link.identifier.isEmpty, ids.insert(link.identifier).inserted,
                  link.pairIndex < physiology.pairCount, link.environmentIndex < physiology.environmentCount,
                  Int(link.molecularSpeciesIndex) < species.count, link.molecularLaneIndex < molecularLaneCount else {
                throw VivoArtifactValidationError.invalid("invalid or duplicate feedback link \(link.identifier)")
            }
            if link.physiologyMode == .replace {
                guard physiologyReplacements.insert(UInt64(link.pairIndex) << 32 | UInt64(link.environmentIndex)).inserted else {
                    throw VivoArtifactValidationError.invalid("multiple replacement links target one physiology element")
                }
            }
            try link.transfer.validate()
        }
    }
}

public struct VivoMolecularPhysiologyCouplingCompiler: Sendable {
    public struct Limits: Codable, Sendable, Equatable {
        public var maximumExpandedLinks: Int
        public init(maximumExpandedLinks: Int = 1_048_576) { self.maximumExpandedLinks = maximumExpandedLinks }
    }
    private let units: VivoUnitSystem
    private let limits: Limits
    public init(units: VivoUnitSystem = .standard, limits: Limits = .init()) { self.units = units; self.limits = limits }
    public func compile(_ coupling: VivoMolecularPhysiologyCouplingPack, programPack: VivoProgramPack,
                        physiology: PreparedVivoPhysiologyModel, molecularLaneCount: UInt32) throws -> PreparedVivoMolecularPhysiologyBridge {
        guard coupling.schema == VivoMolecularPhysiologyCouplingPack.schema,
              coupling.programFingerprint == programPack.header.contentFingerprint,
              coupling.physiologyFingerprint == physiology.fingerprint, molecularLaneCount > 0,
              (1...1_048_576).contains(limits.maximumExpandedLinks) else {
            throw VivoArtifactValidationError.incompatible("bridge identity, schema or limits mismatch")
        }
        try coupling.policy.validate()
        try VivoPhysiologyModelValidator.validate(physiology)
        let species = try programPack.speciesMetadata()
        let speciesByID = Dictionary(uniqueKeysWithValues: species.enumerated().map { ($0.element.identifier, (index: UInt32($0.offset), metadata: $0.element)) })
        let analytes = Dictionary(uniqueKeysWithValues: physiology.analytes.map { ($0.identifier, $0) })
        var exposure: [PreparedVivoPhysiologyToMolecularLink] = [], feedback: [PreparedVivoMolecularToPhysiologyLink] = []
        var identifiers = Set<String>()
        for link in coupling.physiologyToMolecular {
            guard !link.identifier.isEmpty, identifiers.insert(link.identifier).inserted,
                  let source = analytes[link.analyte], let target = speciesByID[link.molecularSpecies], target.metadata.isExternallyOwned else {
                throw VivoArtifactValidationError.unresolved("invalid exposure link \(link.identifier)")
            }
            let pair = try physiology.pairIndex(analyte: link.analyte, compartment: link.compartment)
            let mappings = try mappings(sources: link.physiologyEnvironments, sourceCount: physiology.environmentCount,
                                         destinations: link.molecularLanes, destinationCount: molecularLaneCount,
                                         mode: link.mapping, budget: limits.maximumExpandedLinks - exposure.count - feedback.count)
            let transfer = try compileTransfer(link.transfer, from: source.concentrationUnit, to: target.metadata.unit, rate: false)
            for (environment, lane) in mappings {
                exposure.append(.init(identifier: "\(link.identifier)[\(environment):\(lane)]", pairIndex: pair,
                                       environmentIndex: environment, molecularSpeciesIndex: target.index, molecularLaneIndex: lane,
                                       molecularMode: link.molecularMode, transfer: transfer, required: link.required))
            }
        }
        for link in coupling.molecularToPhysiology {
            guard !link.identifier.isEmpty, identifiers.insert(link.identifier).inserted,
                  let source = speciesByID[link.molecularSpecies], let target = analytes[link.analyte] else {
                throw VivoArtifactValidationError.unresolved("invalid feedback link \(link.identifier)")
            }
            let pair = try physiology.pairIndex(analyte: link.analyte, compartment: link.compartment)
            let mappings = try mappings(sources: link.molecularLanes, sourceCount: molecularLaneCount,
                                         destinations: link.physiologyEnvironments, destinationCount: physiology.environmentCount,
                                         mode: link.mapping, budget: limits.maximumExpandedLinks - exposure.count - feedback.count)
            let transfer = try compileTransfer(link.transfer, from: source.metadata.unit, to: target.concentrationUnit, rate: link.physiologyMode == .rate)
            for (lane, environment) in mappings {
                feedback.append(.init(identifier: "\(link.identifier)[\(lane):\(environment)]", molecularSpeciesIndex: source.index,
                                       molecularLaneIndex: lane, pairIndex: pair, environmentIndex: environment,
                                       physiologyMode: link.physiologyMode, transfer: transfer, required: link.required))
            }
        }
        let compiled = PreparedVivoMolecularPhysiologyBridge(fingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(coupling)),
                                                              programFingerprint: coupling.programFingerprint, physiologyFingerprint: coupling.physiologyFingerprint,
                                                              physiologyToMolecular: exposure.sorted { $0.identifier < $1.identifier },
                                                              molecularToPhysiology: feedback.sorted { $0.identifier < $1.identifier }, policy: coupling.policy, labels: coupling.labels)
        try compiled.validate(program: programPack, physiology: physiology, molecularLaneCount: molecularLaneCount)
        return compiled
    }
    private func compileTransfer(_ input: VivoBridgeTransfer, from sourceUnit: String, to destinationUnit: String, rate: Bool) throws -> PreparedVivoBridgeTransfer {
        guard let source = units.definition(for: sourceUnit), let destination = units.definition(for: destinationUnit), input.gain.isFinite,
              source.offsetToSI == 0, destination.offsetToSI == 0,
              destination.dimension == VivoDimension(length: -3, amount: 1) else {
            throw VivoArtifactValidationError.incompatible("bridge requires explicit molar concentration units or concentration-rate source units")
        }
        let expected = rate ? VivoDimension(length: -3, time: -1, amount: 1) : destination.dimension
        guard source.dimension == expected else { throw VivoArtifactValidationError.incompatible("bridge source dimension does not match update mode") }
        func convert(_ quantity: VivoQuantity?) throws -> Double? {
            guard let quantity else { return nil }
            guard quantity.value.isFinite, let definition = units.definition(for: quantity.unit), definition.dimension == expected,
                  definition.offsetToSI == 0 else {
                throw VivoArtifactValidationError.incompatible(rate ? "rate offsets and bounds require concentration/time units" : "state offsets and bounds require concentration units")
            }
            let value = quantity.value * (definition.scaleToSI / destination.scaleToSI)
            guard value.isFinite else { throw VivoArtifactValidationError.invalid("bridge quantity conversion overflow") }
            return value
        }
        // For a rate link destination.scaleToSI is the concentration scale;
        // seconds are the runtime's time basis. Offset/min/max follow the SAME
        // rate dimension, not the destination concentration's state dimension.
        let result = try PreparedVivoBridgeTransfer(scale: input.gain * source.scaleToSI / destination.scaleToSI,
                                                     offset: convert(input.offset) ?? 0, minimum: convert(input.minimum),
                                                     maximum: convert(input.maximum), boundPolicy: input.boundPolicy)
        try result.validate()
        return result
    }
    private func mappings(sources: [UInt32]?, sourceCount: UInt32, destinations: [UInt32]?, destinationCount: UInt32,
                          mode: VivoBridgeIndexMapping, budget: Int) throws -> [(UInt32,UInt32)] {
        let aCount = sources.map { UInt64($0.count) } ?? UInt64(sourceCount)
        let bCount = destinations.map { UInt64($0.count) } ?? UInt64(destinationCount)
        guard aCount > 0, bCount > 0, budget > 0, aCount <= UInt64(budget), bCount <= UInt64(budget) else {
            throw VivoArtifactValidationError.invalid("bridge index sets exceed remaining expansion capacity")
        }
        let output: UInt64
        switch mode {
        case .oneToOne:
            guard aCount == bCount else { throw VivoArtifactValidationError.invalid("one-to-one mapping requires equal index counts") }
            output = aCount
        case .broadcastSource:
            guard aCount == 1 else { throw VivoArtifactValidationError.invalid("broadcast source requires one index") }
            output = bCount
        case .broadcastDestination:
            guard bCount == 1 else { throw VivoArtifactValidationError.invalid("broadcast destination requires one index") }
            output = aCount
        case .cartesian:
            let product = aCount.multipliedReportingOverflow(by: bCount)
            guard !product.overflow else { throw VivoArtifactValidationError.invalid("Cartesian link count overflow") }
            output = product.partialValue
        }
        guard output <= UInt64(budget) else { throw VivoArtifactValidationError.invalid("expanded links exceed remaining capacity") }
        let a = sources ?? Array(0..<sourceCount), b = destinations ?? Array(0..<destinationCount)
        guard Set(a).count == a.count, Set(b).count == b.count,
              a.allSatisfy({ $0 < sourceCount }), b.allSatisfy({ $0 < destinationCount }) else {
            throw VivoArtifactValidationError.invalid("bridge indices are duplicated or out of range")
        }
        switch mode {
        case .oneToOne: return Array(zip(a,b))
        case .broadcastSource: return b.map { (a[0], $0) }
        case .broadcastDestination: return a.map { ($0, b[0]) }
        case .cartesian: return a.flatMap { source in b.map { (source, $0) } }
        }
    }
}
