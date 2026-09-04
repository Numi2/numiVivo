import CryptoKit
import Foundation

public enum VivoPopulationBoundaryMode: UInt32, Codable, CaseIterable, Sendable {
    case noFlux = 0
    case periodic = 1
    case absorbing = 2
}

public enum VivoPopulationTransitionMode: UInt32, Codable, CaseIterable, Sendable {
    case constitutive = 0
    case activated = 1
    case repressed = 2
}

public struct VivoPopulationGrid: Codable, Equatable, Sendable {
    public let width: UInt32
    public let height: UInt32
    public let depth: UInt32
    public let spacingXMetres: Float
    public let spacingYMetres: Float
    public let spacingZMetres: Float
    public let boundary: VivoPopulationBoundaryMode

    public init(
        width: UInt32,
        height: UInt32,
        depth: UInt32,
        spacingXMetres: Float,
        spacingYMetres: Float,
        spacingZMetres: Float,
        boundary: VivoPopulationBoundaryMode = .noFlux
    ) {
        self.width = width
        self.height = height
        self.depth = depth
        self.spacingXMetres = spacingXMetres
        self.spacingYMetres = spacingYMetres
        self.spacingZMetres = spacingZMetres
        self.boundary = boundary
    }

    public var voxelCount: UInt32? {
        let xy = UInt64(width).multipliedReportingOverflow(by: UInt64(height))
        guard !xy.overflow else { return nil }
        let xyz = xy.partialValue.multipliedReportingOverflow(by: UInt64(depth))
        guard !xyz.overflow, xyz.partialValue <= UInt64(UInt32.max) else { return nil }
        return UInt32(xyz.partialValue)
    }

    public func validate() throws {
        guard width > 0, height > 0, depth > 0, voxelCount != nil else {
            throw VivoPopulationModelError.invalidGrid("dimensions are zero or overflow UInt32")
        }
        guard spacingXMetres.isFinite,
              spacingYMetres.isFinite,
              spacingZMetres.isFinite,
              spacingXMetres > 0,
              spacingYMetres > 0,
              spacingZMetres > 0 else {
            throw VivoPopulationModelError.invalidGrid("voxel spacing must be finite and positive")
        }
    }
}

public struct VivoPopulationRegulatorField: Codable, Hashable, Sendable {
    public let id: String
    public let unit: String
    public let minimum: Float
    public let maximum: Float

    public init(id: String, unit: String, minimum: Float = 0, maximum: Float = .greatestFiniteMagnitude) {
        self.id = id
        self.unit = unit
        self.minimum = minimum
        self.maximum = maximum
    }
}

public struct VivoPopulationPhenotype: Codable, Hashable, Sendable {
    public let id: String
    public let birthRatePerSecond: Float
    public let deathRatePerSecond: Float
    public let carryingCapacityPerVoxel: Float
    public let diffusionSquareMetresPerSecond: Float
    public let minimumDensity: Float
    public let maximumDensity: Float
    public let externallyMaintained: Bool

    public init(
        id: String,
        birthRatePerSecond: Float,
        deathRatePerSecond: Float,
        carryingCapacityPerVoxel: Float,
        diffusionSquareMetresPerSecond: Float = 0,
        minimumDensity: Float = 0,
        maximumDensity: Float? = nil,
        externallyMaintained: Bool = false
    ) {
        self.id = id
        self.birthRatePerSecond = birthRatePerSecond
        self.deathRatePerSecond = deathRatePerSecond
        self.carryingCapacityPerVoxel = carryingCapacityPerVoxel
        self.diffusionSquareMetresPerSecond = diffusionSquareMetresPerSecond
        self.minimumDensity = minimumDensity
        self.maximumDensity = maximumDensity ?? carryingCapacityPerVoxel * 4
        self.externallyMaintained = externallyMaintained
    }
}

public struct VivoPopulationTransition: Codable, Hashable, Sendable {
    public let id: String
    public let sourcePhenotype: String
    public let destinationPhenotype: String
    public let mode: VivoPopulationTransitionMode
    public let regulatorField: String?
    public let baseRatePerSecond: Float
    public let maximumRegulatedRatePerSecond: Float
    public let threshold: Float
    public let hillCoefficient: Float

    public init(
        id: String,
        sourcePhenotype: String,
        destinationPhenotype: String,
        mode: VivoPopulationTransitionMode = .constitutive,
        regulatorField: String? = nil,
        baseRatePerSecond: Float,
        maximumRegulatedRatePerSecond: Float = 0,
        threshold: Float = 1,
        hillCoefficient: Float = 1
    ) {
        self.id = id
        self.sourcePhenotype = sourcePhenotype
        self.destinationPhenotype = destinationPhenotype
        self.mode = mode
        self.regulatorField = regulatorField
        self.baseRatePerSecond = baseRatePerSecond
        self.maximumRegulatedRatePerSecond = maximumRegulatedRatePerSecond
        self.threshold = threshold
        self.hillCoefficient = hillCoefficient
    }
}

/// Actor phenotype changes the target phenotype's per-capita growth rate.
/// `coefficientPerDensityPerSecond` may be positive or negative.
public struct VivoPopulationInteraction: Codable, Hashable, Sendable {
    public let actorPhenotype: String
    public let targetPhenotype: String
    public let coefficientPerDensityPerSecond: Float

    public init(
        actorPhenotype: String,
        targetPhenotype: String,
        coefficientPerDensityPerSecond: Float
    ) {
        self.actorPhenotype = actorPhenotype
        self.targetPhenotype = targetPhenotype
        self.coefficientPerDensityPerSecond = coefficientPerDensityPerSecond
    }
}

public struct VivoPopulationModel: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let name: String
    public let grid: VivoPopulationGrid
    public let regulatorFields: [VivoPopulationRegulatorField]
    public let phenotypes: [VivoPopulationPhenotype]
    public let transitions: [VivoPopulationTransition]
    public let interactions: [VivoPopulationInteraction]
    public let fingerprint: String

    public init(
        name: String,
        grid: VivoPopulationGrid,
        regulatorFields: [VivoPopulationRegulatorField] = [],
        phenotypes: [VivoPopulationPhenotype],
        transitions: [VivoPopulationTransition] = [],
        interactions: [VivoPopulationInteraction] = []
    ) throws {
        let unsigned = VivoPopulationModel(
            schemaVersion: 1,
            name: name,
            grid: grid,
            regulatorFields: regulatorFields,
            phenotypes: phenotypes,
            transitions: transitions,
            interactions: interactions,
            fingerprint: ""
        )
        try unsigned.validate()
        self = .init(
            schemaVersion: unsigned.schemaVersion,
            name: unsigned.name,
            grid: unsigned.grid,
            regulatorFields: unsigned.regulatorFields,
            phenotypes: unsigned.phenotypes,
            transitions: unsigned.transitions,
            interactions: unsigned.interactions,
            fingerprint: try Self.fingerprint(unsigned)
        )
    }

    private init(
        schemaVersion: UInt32,
        name: String,
        grid: VivoPopulationGrid,
        regulatorFields: [VivoPopulationRegulatorField],
        phenotypes: [VivoPopulationPhenotype],
        transitions: [VivoPopulationTransition],
        interactions: [VivoPopulationInteraction],
        fingerprint: String
    ) {
        self.schemaVersion = schemaVersion
        self.name = name
        self.grid = grid
        self.regulatorFields = regulatorFields
        self.phenotypes = phenotypes
        self.transitions = transitions
        self.interactions = interactions
        self.fingerprint = fingerprint
    }

    public func validate() throws {
        guard schemaVersion == 1 else { throw VivoPopulationModelError.unsupportedVersion(schemaVersion) }
        guard !name.isEmpty, name.utf8.count <= 256 else {
            throw VivoPopulationModelError.invalidModel("name is empty or too long")
        }
        try grid.validate()
        guard !phenotypes.isEmpty, phenotypes.count <= Int(UInt32.max) else {
            throw VivoPopulationModelError.invalidModel("phenotype table is empty or exceeds the ABI")
        }
        guard Set(phenotypes.map(\.id)).count == phenotypes.count else {
            throw VivoPopulationModelError.invalidModel("phenotype identifiers are not unique")
        }
        guard Set(regulatorFields.map(\.id)).count == regulatorFields.count else {
            throw VivoPopulationModelError.invalidModel("regulator field identifiers are not unique")
        }
        guard Set(transitions.map(\.id)).count == transitions.count else {
            throw VivoPopulationModelError.invalidModel("transition identifiers are not unique")
        }

        let phenotypeIDs = Set(phenotypes.map(\.id))
        let fieldIDs = Set(regulatorFields.map(\.id))
        for (index, field) in regulatorFields.enumerated() {
            guard !field.id.isEmpty, !field.unit.isEmpty,
                  field.minimum.isFinite,
                  field.maximum.isFinite,
                  field.minimum <= field.maximum else {
                throw VivoPopulationModelError.invalidField(index, "identifier, unit or bounds are invalid")
            }
        }
        for (index, phenotype) in phenotypes.enumerated() {
            let finite = [
                phenotype.birthRatePerSecond,
                phenotype.deathRatePerSecond,
                phenotype.carryingCapacityPerVoxel,
                phenotype.diffusionSquareMetresPerSecond,
                phenotype.minimumDensity,
                phenotype.maximumDensity
            ].allSatisfy(\.isFinite)
            guard !phenotype.id.isEmpty, phenotype.id.utf8.count <= 256, finite else {
                throw VivoPopulationModelError.invalidPhenotype(index, "identifier or parameter is invalid")
            }
            guard phenotype.birthRatePerSecond >= 0,
                  phenotype.deathRatePerSecond >= 0,
                  phenotype.carryingCapacityPerVoxel > 0,
                  phenotype.diffusionSquareMetresPerSecond >= 0,
                  phenotype.minimumDensity >= 0,
                  phenotype.maximumDensity >= phenotype.minimumDensity else {
                throw VivoPopulationModelError.invalidPhenotype(index, "rates, capacity, diffusion or bounds are outside valid ranges")
            }
        }
        for (index, transition) in transitions.enumerated() {
            guard !transition.id.isEmpty,
                  phenotypeIDs.contains(transition.sourcePhenotype),
                  phenotypeIDs.contains(transition.destinationPhenotype),
                  transition.sourcePhenotype != transition.destinationPhenotype else {
                throw VivoPopulationModelError.invalidTransition(index, "phenotype identities are invalid")
            }
            let finite = [
                transition.baseRatePerSecond,
                transition.maximumRegulatedRatePerSecond,
                transition.threshold,
                transition.hillCoefficient
            ].allSatisfy(\.isFinite)
            guard finite,
                  transition.baseRatePerSecond >= 0,
                  transition.maximumRegulatedRatePerSecond >= 0,
                  transition.threshold > 0,
                  transition.hillCoefficient >= 1 else {
                throw VivoPopulationModelError.invalidTransition(index, "kinetic parameters are invalid")
            }
            switch transition.mode {
            case .constitutive:
                guard transition.regulatorField == nil else {
                    throw VivoPopulationModelError.invalidTransition(index, "constitutive transition must not declare a regulator")
                }
            case .activated, .repressed:
                guard let field = transition.regulatorField, fieldIDs.contains(field) else {
                    throw VivoPopulationModelError.invalidTransition(index, "regulated transition references an unknown field")
                }
            }
        }
        var interactionPairs = Set<String>()
        for (index, interaction) in interactions.enumerated() {
            guard phenotypeIDs.contains(interaction.actorPhenotype),
                  phenotypeIDs.contains(interaction.targetPhenotype),
                  interaction.coefficientPerDensityPerSecond.isFinite else {
                throw VivoPopulationModelError.invalidInteraction(index)
            }
            let key = "\(interaction.targetPhenotype)\u{0}\(interaction.actorPhenotype)"
            guard interactionPairs.insert(key).inserted else {
                throw VivoPopulationModelError.invalidModel("duplicate population interaction \(key)")
            }
        }
    }

    private static func fingerprint(_ model: VivoPopulationModel) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return SHA256.hash(data: try encoder.encode(model)).map { String(format: "%02x", $0) }.joined()
    }
}

public enum VivoPopulationModelError: Error, LocalizedError, Sendable {
    case unsupportedVersion(UInt32)
    case invalidGrid(String)
    case invalidModel(String)
    case invalidField(Int, String)
    case invalidPhenotype(Int, String)
    case invalidTransition(Int, String)
    case invalidInteraction(Int)
    case arithmeticOverflow

    public var errorDescription: String? {
        switch self {
        case .unsupportedVersion(let version): return "unsupported population model version \(version)"
        case .invalidGrid(let reason): return "invalid population grid: \(reason)"
        case .invalidModel(let reason): return "invalid population model: \(reason)"
        case .invalidField(let index, let reason): return "invalid regulator field \(index): \(reason)"
        case .invalidPhenotype(let index, let reason): return "invalid phenotype \(index): \(reason)"
        case .invalidTransition(let index, let reason): return "invalid transition \(index): \(reason)"
        case .invalidInteraction(let index): return "invalid population interaction \(index)"
        case .arithmeticOverflow: return "population model arithmetic overflow"
        }
    }
}

struct VivoPopulationPhenotypeABI {
    var birthRate: Float
    var deathRate: Float
    var carryingCapacity: Float
    var diffusionCoefficient: Float
    var minimumDensity: Float
    var maximumDensity: Float
    var flags: UInt32
    var reserved: UInt32 = 0
}

struct VivoPopulationTransitionABI {
    var sourcePhenotype: UInt32
    var destinationPhenotype: UInt32
    var regulatorField: UInt32
    var mode: UInt32
    var baseRate: Float
    var maximumRegulatedRate: Float
    var threshold: Float
    var hillCoefficient: Float
}

struct VivoPopulationCompiledModel: Sendable {
    let phenotypes: [VivoPopulationPhenotypeABI]
    let transitions: [VivoPopulationTransitionABI]
    let interactionMatrix: [Float]
    let conservativeKineticRate: Float
    let conservativeDiffusionRate: Float

    init(model: VivoPopulationModel) throws {
        let phenotypeIndex = Dictionary(uniqueKeysWithValues: model.phenotypes.enumerated().map { ($1.id, UInt32($0)) })
        let fieldIndex = Dictionary(uniqueKeysWithValues: model.regulatorFields.enumerated().map { ($1.id, UInt32($0)) })
        phenotypes = model.phenotypes.map {
            VivoPopulationPhenotypeABI(
                birthRate: $0.birthRatePerSecond,
                deathRate: $0.deathRatePerSecond,
                carryingCapacity: $0.carryingCapacityPerVoxel,
                diffusionCoefficient: $0.diffusionSquareMetresPerSecond,
                minimumDensity: $0.minimumDensity,
                maximumDensity: $0.maximumDensity,
                flags: $0.externallyMaintained ? 1 : 0
            )
        }
        transitions = try model.transitions.map { transition in
            guard let source = phenotypeIndex[transition.sourcePhenotype],
                  let destination = phenotypeIndex[transition.destinationPhenotype] else {
                throw VivoPopulationModelError.invalidModel("transition phenotype resolution failed")
            }
            let regulator: UInt32
            if let field = transition.regulatorField {
                guard let resolved = fieldIndex[field] else {
                    throw VivoPopulationModelError.invalidModel("transition field resolution failed")
                }
                regulator = resolved
            } else {
                regulator = UInt32.max
            }
            return VivoPopulationTransitionABI(
                sourcePhenotype: source,
                destinationPhenotype: destination,
                regulatorField: regulator,
                mode: transition.mode.rawValue,
                baseRate: transition.baseRatePerSecond,
                maximumRegulatedRate: transition.maximumRegulatedRatePerSecond,
                threshold: transition.threshold,
                hillCoefficient: transition.hillCoefficient
            )
        }

        let count = model.phenotypes.count
        let matrixCount = count.multipliedReportingOverflow(by: count)
        guard !matrixCount.overflow else { throw VivoPopulationModelError.arithmeticOverflow }
        var matrix = [Float](repeating: 0, count: matrixCount.partialValue)
        for interaction in model.interactions {
            guard let actor = phenotypeIndex[interaction.actorPhenotype],
                  let target = phenotypeIndex[interaction.targetPhenotype] else {
                throw VivoPopulationModelError.invalidModel("interaction phenotype resolution failed")
            }
            matrix[Int(target) * count + Int(actor)] = interaction.coefficientPerDensityPerSecond
        }
        interactionMatrix = matrix

        var maximumKinetic: Float = 0
        for (targetIndex, phenotype) in model.phenotypes.enumerated() {
            var rate = phenotype.birthRatePerSecond + phenotype.deathRatePerSecond
            for transition in model.transitions where transition.sourcePhenotype == phenotype.id {
                rate += transition.baseRatePerSecond + transition.maximumRegulatedRatePerSecond
            }
            for actorIndex in model.phenotypes.indices {
                rate += abs(matrix[targetIndex * count + actorIndex]) * model.phenotypes[actorIndex].maximumDensity
            }
            maximumKinetic = max(maximumKinetic, rate)
        }
        conservativeKineticRate = maximumKinetic

        let inverseX2 = 1 / (model.grid.spacingXMetres * model.grid.spacingXMetres)
        let inverseY2 = 1 / (model.grid.spacingYMetres * model.grid.spacingYMetres)
        let inverseZ2 = 1 / (model.grid.spacingZMetres * model.grid.spacingZMetres)
        let maximumDiffusion = model.phenotypes.map(\.diffusionSquareMetresPerSecond).max() ?? 0
        conservativeDiffusionRate = 2 * maximumDiffusion * (inverseX2 + inverseY2 + inverseZ2)
    }

    func maximumStableStep(maximumVelocity: SIMD3<Float>) -> Float {
        let advectionRate = abs(maximumVelocity.x) / 1 + abs(maximumVelocity.y) / 1 + abs(maximumVelocity.z) / 1
        let total = conservativeKineticRate + conservativeDiffusionRate + advectionRate
        guard total.isFinite, total > 0 else { return 60 }
        return max(1e-9, min(60, 0.35 / total))
    }
}
