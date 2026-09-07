import Foundation

/// Algebraic fixtures, deliberately not named as chemical reaction benchmarks.
/// Identity overlaps are physically exact here because every point is expressed
/// in the same fixed orthonormal model basis, not inferred from equal dimensions.
public enum VivoPropertyRefinementExamples {
    public static func algebraicPath(selectedCI: Bool = false) throws -> VivoPropertyDirectedSpaceRequest {
        let names = ["reactant","held-out-early","barrier-point","held-out-late","product"]
        let strengths = [0.01,0.20,0.30,0.20,0.01]
        let scalars = [0.0,0.15,0.25,0.16,-0.05]
        let points = try names.indices.map { i -> VivoECCPathPoint in
            var one = VivoQMMatrix(3,3)
            one[1,1] = 1; one[2,2] = 1.4
            one[0,1] = 0.12; one[1,0] = 0.12
            one[0,2] = strengths[i]; one[2,0] = strengths[i]
            let h = VivoEmbeddedHamiltonian(orbitalIdentifiers: ["required","common-correction","reaction-sensitive"],
                alphaElectrons: 1,betaElectrons: 0,oneElectron: one,twoElectron: [Double](repeating: 0,count: 81),
                constantEnergyHartree: scalars[i],energyReference: "fixed orthonormal algebraic basis; not a molecular reaction",
                provenance: ["fixture":"property-directed-algebraic-v1","coordinateUnits":"dimensionless model index"])
            return .init(identifier: names[i],hamiltonian: h,overlapWithPrevious: i == 0 ? nil : try .identity(3))
        }
        let solver: VivoSpaceRefinementSolver = selectedCI
            ? .selectedCI(configuration: .init(maximumDeterminants: 3,selectionBatchSize: 1,
                minimumSelectionContributionHartree: 0,pt2ToleranceHartree: 1e-10,eigenResidualTolerance: 1e-12,
                maximumDavidsonSubspace: 8,maximumExternalDeterminants: 32,fullResidualTolerance: 1e-9))
            : .directCI(configuration: .init(residualTolerance: 1e-12))
        return .init(identifier: selectedCI ? "algebraic-selected-ci-refinement" : "algebraic-direct-ci-refinement",
            points: points,transportGroups: [[0],[1],[2]],initialSpace: .init(active: [0]),mandatoryActiveOrbitals: [0],
            candidateBlocks: [.init(identifier: "common-correction",orbitals: [1]),
                              .init(identifier: "reaction-sensitive",orbitals: [2])],
            discoveryPointIdentifiers: [names[0],names[2],names[4]],confirmationPointIdentifiers: [names[1],names[3]],
            target: .init(barrierPointIdentifier: names[2],maximumBarrierShiftHartree: 0.003,
                          maximumRelativeProfileShiftHartree: 0.003),solver: solver,
            maximumRounds: 4,maximumActiveOrbitals: 3,maximumPointEvaluations: 64,
            minimumStateOverlapSquared: 0.8,collectOrbitalInformation: false)
    }
}
