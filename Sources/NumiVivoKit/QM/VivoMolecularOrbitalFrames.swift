import Foundation

/// Shared source-bound molecular-to-orbital preparation for ECC paths and
/// property refinement. Cross-geometry overlaps are integrated in the AO basis.
public enum VivoMolecularOrbitalFrames {
    public static func prepare(snapshots: [VivoMolecularPathSnapshot], basis: VivoGaussianBasis,
                               reference configuration: VivoSCFConfiguration,
                               integrals: [VivoAOIntegrals], references: [VivoHartreeFockResult],
                               rotations: [VivoQMMatrix]? = nil, orbitalPrefix: String = "path-orbital",
                               budget: VivoChemistryBudget = .init()) throws -> [VivoECCPathPoint] {
        try budget.validate(); try configuration.validate()
        guard (2...128).contains(snapshots.count), integrals.count == snapshots.count,
              references.count == snapshots.count, rotations == nil || rotations?.count == snapshots.count,
              !orbitalPrefix.isEmpty, orbitalPrefix.utf8.count <= 256,
              configuration.reference == .restricted else {
            throw VivoChemistryError.invalid("molecular orbital-frame sources, count or restricted reference")
        }
        let first = snapshots[0].system, n = integrals[0].count
        guard n > 0, n <= budget.maximumBasisFunctions else { throw VivoChemistryError.invalid("molecular orbital-frame dimension") }
        _ = try budget.elements([snapshots.count,n,n,n,n],simultaneousArrays: 16)
        var coefficients: [VivoQMMatrix] = [], points: [VivoECCPathPoint] = []
        for i in snapshots.indices {
            let system = snapshots[i].system, ao = integrals[i], hf = references[i]
            guard ao.sourceBasis == basis, ao.sourceSystem == system, ao.count == n,
                  system.nuclei.map(\.atomicNumber) == first.nuclei.map(\.atomicNumber),
                  system.nuclei.map(\.structureAtomIndex) == first.nuclei.map(\.structureAtomIndex),
                  system.alphaElectrons == first.alphaElectrons, system.betaElectrons == first.betaElectrons,
                  system.pointCharges == first.pointCharges else {
                throw VivoChemistryError.invalid("molecular path changes source, atom mapping, electron sector or frozen charges")
            }
            try VivoHartreeFock.validate(result: hf,system: system,integrals: ao,configuration: configuration,budget: budget)
            let u = try rotations?[i] ?? VivoQMMatrix.identity(n)
            guard u.rows == n, u.columns == n,
                  try u.transposed.multiplied(by: u).adding(.identity(n),scale: -1).frobeniusNorm < 1e-8 else {
                throw VivoChemistryError.invalid("molecular preparation orbital rotation")
            }
            let c = try hf.alphaCoefficients.multiplied(by: u)
            coefficients.append(c)
            let h = try VivoEmbeddedHamiltonian.fromAO(ao,coefficients: c,
                alphaElectrons: system.alphaElectrons,betaElectrons: system.betaElectrons,
                orbitalIdentifiers: (0..<n).map { "\(orbitalPrefix)-\($0)" },
                energyReference: "physical electronic Hamiltonian; scalar inherited from AO integrals once; no thermal/standard-state correction",budget: budget)
            var overlap: VivoQMMatrix?
            if i > 0 {
                let cross = try VivoGaussianIntegralEngine.crossOverlap(leftSystem: snapshots[i-1].system,leftBasis: basis,
                    rightSystem: system,rightBasis: basis,budget: budget)
                overlap = try coefficients[i-1].transposed.multiplied(by: cross).multiplied(by: c)
            }
            points.append(.init(identifier: snapshots[i].identifier,hamiltonian: h,overlapWithPrevious: overlap))
        }
        return points
    }
}
