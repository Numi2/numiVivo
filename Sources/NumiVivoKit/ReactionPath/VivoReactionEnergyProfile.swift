import Foundation

public enum VivoReactionEnergyReference: String, Codable, Sendable {
    case precomplex, isolatedReactants, product, explicitState
}
public struct VivoReactionEnergyPoint: Codable, Sendable {
    public var identifier: String
    public var coordinate: Double
    public var energyHartree: Double
    public var calculationFingerprint: VivoFingerprint
    public var methodContractFingerprint: VivoFingerprint
    public var frozenEnvironmentFingerprint: VivoFingerprint
    public var activeOrbitalIDs: [String]
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public init(identifier: String, coordinate: Double, energyHartree: Double,
                calculationFingerprint: VivoFingerprint, methodContractFingerprint: VivoFingerprint,
                frozenEnvironmentFingerprint: VivoFingerprint, activeOrbitalIDs: [String],
                alphaElectrons: Int, betaElectrons: Int) {
        self.identifier = identifier; self.coordinate = coordinate; self.energyHartree = energyHartree
        self.calculationFingerprint = calculationFingerprint; self.methodContractFingerprint = methodContractFingerprint
        self.frozenEnvironmentFingerprint = frozenEnvironmentFingerprint; self.activeOrbitalIDs = activeOrbitalIDs
        self.alphaElectrons = alphaElectrons; self.betaElectrons = betaElectrons
    }
}
public struct VivoReactionEnergyProfile: Codable, Sendable {
    public var points: [VivoReactionEnergyPoint]
    public var reference: VivoReactionEnergyReference
    public var referencePointIdentifier: String
    public var coordinateDefinition: String
    public var coordinateUnit: String
    public init(points: [VivoReactionEnergyPoint], reference: VivoReactionEnergyReference,
                referencePointIdentifier: String, coordinateDefinition: String, coordinateUnit: String) {
        self.points = points; self.reference = reference; self.referencePointIdentifier = referencePointIdentifier
        self.coordinateDefinition = coordinateDefinition; self.coordinateUnit = coordinateUnit
    }
    public func validate() throws {
        guard let first = points.first, (2...100_000).contains(points.count),
              Set(points.map(\.identifier)).count == points.count,
              points.contains(where: { $0.identifier == referencePointIdentifier }),
              !coordinateDefinition.isEmpty, !coordinateUnit.isEmpty,
              !first.activeOrbitalIDs.isEmpty, Set(first.activeOrbitalIDs).count == first.activeOrbitalIDs.count,
              (0...first.activeOrbitalIDs.count).contains(first.alphaElectrons),
              (0...first.activeOrbitalIDs.count).contains(first.betaElectrons) else {
            throw VivoQMError.invalid("reaction energy-profile identity, reference, or electron sector")
        }
        for (index,point) in points.enumerated() {
            guard !point.identifier.isEmpty, point.coordinate.isFinite, point.energyHartree.isFinite,
                  point.methodContractFingerprint == first.methodContractFingerprint,
                  point.frozenEnvironmentFingerprint == first.frozenEnvironmentFingerprint,
                  point.activeOrbitalIDs == first.activeOrbitalIDs,
                  point.alphaElectrons == first.alphaElectrons, point.betaElectrons == first.betaElectrons else {
                throw VivoQMError.invalid("reaction profile changes method, environment, active-space identity, or electron sector")
            }
            if index > 0, point.coordinate <= points[index-1].coordinate {
                throw VivoQMError.invalid("reaction coordinates must be strictly increasing in declared path order")
            }
        }
    }
    /// The maximum on a sampled scan is NOT a verified first-order saddle, an
    /// IRC-validated transition state, a Gibbs barrier, or an activation rate.
    public func sampledElectronicBarrier() throws -> VivoSampledElectronicBarrier {
        try validate()
        let referencePoint = points.first { $0.identifier == referencePointIdentifier }!
        let maximum = points.max { $0.energyHartree < $1.energyHartree }!
        let difference = maximum.energyHartree-referencePoint.energyHartree
        return .init(reference:reference,referencePointIdentifier:referencePointIdentifier,
                     maximumSampleIdentifier:maximum.identifier,
                     differenceHartree:difference,differenceKJPerMol:difference*VivoQMUnits.hartreeInKJPerMol,
                     transitionStateVerified:false,freeEnergyCorrectionsIncluded:false)
    }
}
public struct VivoSampledElectronicBarrier: Codable, Sendable {
    public let reference: VivoReactionEnergyReference
    public let referencePointIdentifier: String
    public let maximumSampleIdentifier: String
    public let differenceHartree: Double
    public let differenceKJPerMol: Double
    public let transitionStateVerified: Bool
    public let freeEnergyCorrectionsIncluded: Bool
}
