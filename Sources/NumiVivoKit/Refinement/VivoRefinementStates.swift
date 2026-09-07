import Foundation

/// The labels name adiabatic energy ranks, not persistent diabatic characters.
/// Contiguous groups may mix internally, but their complete state subspaces
/// must survive both geometry changes and active-space expansion. Singleton
/// groups enforce ordinary physical wavefunction-overlap tracking.
public struct VivoRefinementStatePolicy: Codable, Sendable, Equatable {
    public let labels: [String]
    public let groups: [[Int]]
    public let minimumIntergroupGapHartree: Double
    public init(labels: [String], groups: [[Int]], minimumIntergroupGapHartree: Double = 1e-8) {
        self.labels = labels; self.groups = groups
        self.minimumIntergroupGapHartree = minimumIntergroupGapHartree
    }
    public static let groundState = Self(labels: ["ground"], groups: [[0]])
    public func validate() throws {
        guard (1...12).contains(labels.count), Set(labels).count == labels.count,
              labels.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }),
              !groups.isEmpty, groups.allSatisfy({ !$0.isEmpty }),
              groups.flatMap({ $0 }) == Array(labels.indices),
              minimumIntergroupGapHartree.isFinite, minimumIntergroupGapHartree > 0 else {
            throw VivoChemistryError.invalid("refinement state labels, contiguous subspaces or gap threshold")
        }
    }
    func validateSpectrum(_ energies: [Double]) throws {
        try validate()
        guard energies.count == labels.count, energies.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("refinement state spectrum shape")
        }
        for i in 1..<energies.count where energies[i] < energies[i-1] - 1e-10 {
            throw VivoChemistryError.invalid("refinement states must be in adiabatic energy order")
        }
        for group in groups.dropFirst() {
            let i = group[0]
            guard energies[i] - energies[i-1] >= minimumIntergroupGapHartree else {
                throw VivoChemistryError.convergence("near-degenerate roots straddle separate declared state groups")
            }
        }
    }
    /// Singular values test the entire overlap block. Pairwise greedy matching
    /// is not sufficient when degenerate states rotate within their subspace.
    func minimumRetainedOverlapSquared(_ overlap: VivoQMMatrix, threshold: Double) throws -> Double {
        guard overlap.rows == labels.count, overlap.columns == labels.count,
              overlap.values.allSatisfy(\.isFinite), threshold.isFinite, threshold > 0, threshold <= 1 else {
            throw VivoChemistryError.invalid("refinement state-overlap shape or threshold")
        }
        var minimum = 1.0
        for group in groups {
            var block = VivoQMMatrix(group.count,group.count)
            for (i,p) in group.enumerated() { for (j,q) in group.enumerated() { block[i,j] = overlap[p,q] } }
            let eigen = try VivoQMDenseAlgebra.symmetricEigen(block.transposed.multiplied(by: block), tolerance: 1e-14)
            guard eigen.values.allSatisfy({ $0.isFinite && $0 >= -1e-9 && $0 <= 1+1e-6 }),
                  let smallest = eigen.values.first, smallest >= threshold else {
                throw VivoChemistryError.convergence("tracked electronic state subspace left the retained wavefunction space")
            }
            minimum = min(minimum,max(0,min(1,smallest)))
        }
        return minimum
    }
}
public struct VivoRefinementRootObservation: Codable, Sendable, Equatable {
    public let label: String
    public let energyHartree: Double
    public let determinantCount: Int
    public let residualNorm: Double
    public let orbitalInformation: VivoSelectiveOrbitalInformationResult?
}
public struct VivoRefinementRootShift: Codable, Sendable, Equatable {
    public let label: String
    public let forwardHartree: Double
    public let reverseHartree: Double
    public let reactionHartree: Double
    public let maximumRelativeProfileHartree: Double
    public let maximumAbsoluteEnergyHartree: Double
}
