import CryptoKit
import Foundation

public struct VivoPartitionCompartment: Codable, Hashable, Sendable {
    public let id: String
    public let volumeCubicMetres: Float

    public init(id: String, volumeCubicMetres: Float) {
        self.id = id
        self.volumeCubicMetres = volumeCubicMetres
    }
}

public struct VivoPartitionAnalyte: Codable, Hashable, Sendable {
    public let id: String
    public let unit: String
    public let minimumConcentration: Float
    public let maximumConcentration: Float

    public init(
        id: String,
        unit: String = "mol/m3",
        minimumConcentration: Float = 0,
        maximumConcentration: Float = .greatestFiniteMagnitude
    ) {
        self.id = id
        self.unit = unit
        self.minimumConcentration = minimumConcentration
        self.maximumConcentration = maximumConcentration
    }
}

/// Reversible transport driven toward a declared target/source concentration ratio.
/// At equilibrium, `target concentration / source concentration` equals
/// `partitionCoefficient` after the declared unbound fractions are applied.
public struct VivoPartitionEdge: Codable, Hashable, Sendable {
    public let id: String
    public let analyte: String
    public let sourceCompartment: String
    public let targetCompartment: String
    public let partitionCoefficient: Float
    public let clearanceCubicMetresPerSecond: Float
    public let sourceUnboundFraction: Float
    public let targetUnboundFraction: Float
    public let evidenceClass: String
    public let evidenceReference: String?

    public init(
        id: String,
        analyte: String,
        sourceCompartment: String,
        targetCompartment: String,
        partitionCoefficient: Float,
        clearanceCubicMetresPerSecond: Float,
        sourceUnboundFraction: Float = 1,
        targetUnboundFraction: Float = 1,
        evidenceClass: String = "assumed",
        evidenceReference: String? = nil
    ) {
        self.id = id
        self.analyte = analyte
        self.sourceCompartment = sourceCompartment
        self.targetCompartment = targetCompartment
        self.partitionCoefficient = partitionCoefficient
        self.clearanceCubicMetresPerSecond = clearanceCubicMetresPerSecond
        self.sourceUnboundFraction = sourceUnboundFraction
        self.targetUnboundFraction = targetUnboundFraction
        self.evidenceClass = evidenceClass
        self.evidenceReference = evidenceReference
    }
}

public struct VivoPhysiologicalPartitionModel: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let name: String
    public let compartments: [VivoPartitionCompartment]
    public let analytes: [VivoPartitionAnalyte]
    public let edges: [VivoPartitionEdge]
    public let fingerprint: String

    public init(
        name: String,
        compartments: [VivoPartitionCompartment],
        analytes: [VivoPartitionAnalyte],
        edges: [VivoPartitionEdge]
    ) throws {
        let unsigned = VivoPhysiologicalPartitionModel(
            schemaVersion: 1,
            name: name,
            compartments: compartments,
            analytes: analytes,
            edges: edges,
            fingerprint: ""
        )
        try unsigned.validate()
        self = .init(
            schemaVersion: unsigned.schemaVersion,
            name: unsigned.name,
            compartments: unsigned.compartments,
            analytes: unsigned.analytes,
            edges: unsigned.edges,
            fingerprint: try Self.fingerprint(unsigned)
        )
    }

    private init(
        schemaVersion: UInt32,
        name: String,
        compartments: [VivoPartitionCompartment],
        analytes: [VivoPartitionAnalyte],
        edges: [VivoPartitionEdge],
        fingerprint: String
    ) {
        self.schemaVersion = schemaVersion
        self.name = name
        self.compartments = compartments
        self.analytes = analytes
        self.edges = edges
        self.fingerprint = fingerprint
    }

    public func validate() throws {
        guard schemaVersion == 1 else {
            throw VivoPartitionModelError.unsupportedVersion(schemaVersion)
        }
        guard !name.isEmpty, name.utf8.count <= 256 else {
            throw VivoPartitionModelError.invalidModel("name is empty or too long")
        }
        guard !compartments.isEmpty, !analytes.isEmpty else {
            throw VivoPartitionModelError.invalidModel("compartments and analytes must be non-empty")
        }
        guard compartments.count <= Int(UInt32.max), analytes.count <= Int(UInt32.max), edges.count <= Int(UInt32.max) else {
            throw VivoPartitionModelError.invalidModel("model exceeds the UInt32 runtime ABI")
        }
        guard Set(compartments.map(\.id)).count == compartments.count else {
            throw VivoPartitionModelError.invalidModel("compartment identifiers are not unique")
        }
        guard Set(analytes.map(\.id)).count == analytes.count else {
            throw VivoPartitionModelError.invalidModel("analyte identifiers are not unique")
        }
        guard Set(edges.map(\.id)).count == edges.count else {
            throw VivoPartitionModelError.invalidModel("edge identifiers are not unique")
        }

        let compartmentIDs = Set(compartments.map(\.id))
        let analyteIDs = Set(analytes.map(\.id))
        for (index, compartment) in compartments.enumerated() {
            guard !compartment.id.isEmpty,
                  compartment.volumeCubicMetres.isFinite,
                  compartment.volumeCubicMetres > 0 else {
                throw VivoPartitionModelError.invalidCompartment(index)
            }
        }
        for (index, analyte) in analytes.enumerated() {
            guard !analyte.id.isEmpty,
                  !analyte.unit.isEmpty,
                  analyte.minimumConcentration.isFinite,
                  analyte.maximumConcentration.isFinite,
                  analyte.minimumConcentration >= 0,
                  analyte.maximumConcentration >= analyte.minimumConcentration else {
                throw VivoPartitionModelError.invalidAnalyte(index)
            }
        }

        let evidenceClasses: Set<String> = [
            "observed", "derived", "calibrated", "inferred", "assumed", "hypothetical"
        ]
        var physicalPairs = Set<String>()
        for (index, edge) in edges.enumerated() {
            let finite = [
                edge.partitionCoefficient,
                edge.clearanceCubicMetresPerSecond,
                edge.sourceUnboundFraction,
                edge.targetUnboundFraction
            ].allSatisfy(\.isFinite)
            guard !edge.id.isEmpty,
                  analyteIDs.contains(edge.analyte),
                  compartmentIDs.contains(edge.sourceCompartment),
                  compartmentIDs.contains(edge.targetCompartment),
                  edge.sourceCompartment != edge.targetCompartment,
                  finite,
                  edge.partitionCoefficient > 0,
                  edge.clearanceCubicMetresPerSecond >= 0,
                  (0...1).contains(edge.sourceUnboundFraction),
                  (0...1).contains(edge.targetUnboundFraction),
                  evidenceClasses.contains(edge.evidenceClass) else {
                throw VivoPartitionModelError.invalidEdge(index)
            }
            let canonicalCompartments = [edge.sourceCompartment, edge.targetCompartment].sorted()
            let pair = "\(edge.analyte)\u{0}\(canonicalCompartments[0])\u{0}\(canonicalCompartments[1])"
            guard physicalPairs.insert(pair).inserted else {
                throw VivoPartitionModelError.invalidModel(
                    "multiple partition edges represent the same analyte-compartment pair: \(pair)"
                )
            }
        }
    }

    private static func fingerprint(_ model: VivoPhysiologicalPartitionModel) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return SHA256.hash(data: try encoder.encode(model)).map { String(format: "%02x", $0) }.joined()
    }
}

public enum VivoPartitionModelError: Error, LocalizedError, Sendable {
    case unsupportedVersion(UInt32)
    case invalidModel(String)
    case invalidCompartment(Int)
    case invalidAnalyte(Int)
    case invalidEdge(Int)
    case arithmeticOverflow

    public var errorDescription: String? {
        switch self {
        case .unsupportedVersion(let version): return "unsupported partition model version \(version)"
        case .invalidModel(let reason): return "invalid partition model: \(reason)"
        case .invalidCompartment(let index): return "invalid partition compartment \(index)"
        case .invalidAnalyte(let index): return "invalid partition analyte \(index)"
        case .invalidEdge(let index): return "invalid partition edge \(index)"
        case .arithmeticOverflow: return "partition model arithmetic overflow"
        }
    }
}

struct VivoPartitionCompartmentABI: Sendable {
    var volume: Float
    var inverseVolume: Float
    var reserved0: Float = 0
    var reserved1: Float = 0
}

struct VivoPartitionAnalyteABI: Sendable {
    var minimumConcentration: Float
    var maximumConcentration: Float
    var reserved0: UInt32 = 0
    var reserved1: UInt32 = 0
}

struct VivoPartitionEdgeABI: Sendable {
    var analyteIndex: UInt32
    var sourceCompartment: UInt32
    var targetCompartment: UInt32
    var flags: UInt32 = 0
    var partitionCoefficient: Float
    var clearance: Float
    var sourceUnboundFraction: Float
    var targetUnboundFraction: Float
}

struct VivoCompiledPartitionModel: Sendable {
    let compartments: [VivoPartitionCompartmentABI]
    let analytes: [VivoPartitionAnalyteABI]
    let edges: [VivoPartitionEdgeABI]
    let conservativeRateBound: Float

    init(model: VivoPhysiologicalPartitionModel) throws {
        let compartmentIndex = Dictionary(uniqueKeysWithValues: model.compartments.enumerated().map { ($1.id, UInt32($0)) })
        let analyteIndex = Dictionary(uniqueKeysWithValues: model.analytes.enumerated().map { ($1.id, UInt32($0)) })
        compartments = model.compartments.map {
            .init(volume: $0.volumeCubicMetres, inverseVolume: 1 / $0.volumeCubicMetres)
        }
        analytes = model.analytes.map {
            .init(minimumConcentration: $0.minimumConcentration, maximumConcentration: $0.maximumConcentration)
        }
        edges = try model.edges.map { edge in
            guard let analyte = analyteIndex[edge.analyte],
                  let source = compartmentIndex[edge.sourceCompartment],
                  let target = compartmentIndex[edge.targetCompartment] else {
                throw VivoPartitionModelError.invalidModel("identifier resolution failed")
            }
            return .init(
                analyteIndex: analyte,
                sourceCompartment: source,
                targetCompartment: target,
                partitionCoefficient: edge.partitionCoefficient,
                clearance: edge.clearanceCubicMetresPerSecond,
                sourceUnboundFraction: edge.sourceUnboundFraction,
                targetUnboundFraction: edge.targetUnboundFraction
            )
        }

        var perStateRate = [Float](repeating: 0, count: model.compartments.count * model.analytes.count)
        for edge in edges {
            let source = Int(edge.analyteIndex) * model.compartments.count + Int(edge.sourceCompartment)
            let target = Int(edge.analyteIndex) * model.compartments.count + Int(edge.targetCompartment)
            perStateRate[source] += edge.clearance * edge.sourceUnboundFraction * compartments[Int(edge.sourceCompartment)].inverseVolume
            perStateRate[target] += edge.clearance * edge.targetUnboundFraction /
                edge.partitionCoefficient * compartments[Int(edge.targetCompartment)].inverseVolume
        }
        conservativeRateBound = perStateRate.max() ?? 0
    }

    var conservativeMaximumStep: Float {
        guard conservativeRateBound.isFinite, conservativeRateBound > 0 else { return 3_600 }
        return max(1e-9, min(3_600, 0.4 / conservativeRateBound))
    }
}
