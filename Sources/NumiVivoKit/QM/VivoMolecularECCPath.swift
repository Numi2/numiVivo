import Foundation

public struct VivoMolecularPathSnapshot: Codable, Sendable, Equatable {
    public let identifier: String
    public let coordinate: Double
    public let system: VivoElectronicSystem
    public init(identifier: String, coordinate: Double, system: VivoElectronicSystem) {
        self.identifier=identifier; self.coordinate=coordinate; self.system=system
    }
}
/// One basis recipe and ordered atom mapping for an entire path. Reactions that
/// change composition, electron population, basis or frozen external charges
/// require a different thermodynamic construction, not this profile contract.
public struct VivoMolecularECCPathRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/molecular-ecc-path/v1"
    public let schema: String
    public let identifier: String
    public let atomIdentifiers: [String]
    public let coordinateUnit: String
    public let snapshots: [VivoMolecularPathSnapshot]
    public let basis: VivoGaussianBasis
    public let reference: VivoSCFConfiguration
    public let configuration: VivoECCPathConfiguration
    public let initialSharedRotation: VivoQMMatrix?
    public let budget: VivoChemistryBudget
    public init(identifier: String, atomIdentifiers: [String], coordinateUnit: String,
                snapshots: [VivoMolecularPathSnapshot], basis: VivoGaussianBasis,
                reference: VivoSCFConfiguration = .init(), configuration: VivoECCPathConfiguration,
                initialSharedRotation: VivoQMMatrix? = nil, budget: VivoChemistryBudget = .init()) {
        schema=Self.schema; self.identifier=identifier; self.atomIdentifiers=atomIdentifiers; self.coordinateUnit=coordinateUnit
        self.snapshots=snapshots; self.basis=basis; self.reference=reference; self.configuration=configuration
        self.initialSharedRotation=initialSharedRotation; self.budget=budget
    }
    public func validate() throws {
        try budget.validate(); try reference.validate()
        guard schema==Self.schema,!identifier.isEmpty,identifier.utf8.count<=1024,
              !coordinateUnit.isEmpty,coordinateUnit.utf8.count<=128,
              (2...128).contains(snapshots.count),let first=snapshots.first,
              atomIdentifiers.count==first.system.nuclei.count,Set(atomIdentifiers).count==atomIdentifiers.count,
              atomIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count<=1024 }),
              Set(snapshots.map(\.identifier)).count==snapshots.count,
              reference.reference == .restricted else {
            throw VivoChemistryError.invalid("molecular ECC path schema, atom mapping, snapshots or restricted reference")
        }
        try basis.validate(nucleusCount:first.system.nuclei.count)
        let atoms=first.system.nuclei.map(\.atomicNumber), mapping=first.system.nuclei.map(\.structureAtomIndex)
        var previous = -Double.infinity
        for snapshot in snapshots {
            try snapshot.system.validate()
            guard !snapshot.identifier.isEmpty,snapshot.identifier.utf8.count<=1024,
                  snapshot.coordinate.isFinite,snapshot.coordinate>previous,
                  snapshot.system.nuclei.map(\.atomicNumber)==atoms,
                  snapshot.system.nuclei.map(\.structureAtomIndex)==mapping,
                  snapshot.system.alphaElectrons==first.system.alphaElectrons,
                  snapshot.system.betaElectrons==first.system.betaElectrons,
                  snapshot.system.alphaElectrons==snapshot.system.betaElectrons,
                  snapshot.system.pointCharges==first.system.pointCharges else {
                throw VivoChemistryError.invalid("path changes atom mapping, composition, charge/spin or frozen electrostatic environment")
            }
            previous=snapshot.coordinate
        }
        let n=try VivoGaussianIntegralEngine.expanded(system:first.system,basis:basis,budget:budget).count
        _=try budget.elements([snapshots.count,n,n,n,n],simultaneousArrays:16)
        if let rotation=initialSharedRotation {
            guard rotation.rows==n,rotation.columns==n,rotation.values.allSatisfy(\.isFinite),
                  try rotation.transposed.multiplied(by:rotation).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {
                throw VivoChemistryError.invalid("molecular path initial shared rotation")
            }
        }
    }
}
public enum VivoMolecularECCPath {
    /// Preparation uses the original geometry-specific AO integrals and RHF
    /// coefficients. Cross-geometry overlap is integrated natively, never
    /// fabricated from energies or assigned an identity by convention.
    public static func prepare(_ request: VivoMolecularECCPathRequest,
                               integrals: [VivoAOIntegrals], references: [VivoHartreeFockResult]) throws -> [VivoECCPathPoint] {
        try request.validate()
        guard integrals.count==request.snapshots.count,references.count==integrals.count else {
            throw VivoChemistryError.invalid("missing path AO/reference snapshot")
        }
        var points:[VivoECCPathPoint]=[]
        for i in integrals.indices {
            let ao=integrals[i], reference=references[i], system=request.snapshots[i].system
            guard ao.sourceBasis==request.basis,ao.sourceSystem==system else {
                throw VivoChemistryError.invalid("path snapshot integral source binding")
            }
            try VivoHartreeFock.validate(result:reference,system:system,integrals:ao,configuration:request.reference,budget:request.budget)
            let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:reference.alphaCoefficients,
                alphaElectrons:system.alphaElectrons,betaElectrons:system.betaElectrons,
                orbitalIdentifiers:(0..<ao.count).map { "path-orbital-\($0)" },
                energyReference:"physical electronic Hamiltonian; scalar inherited from AO integrals once; no thermal/standard-state correction",
                budget:request.budget)
            var overlap:VivoQMMatrix?
            if i>0 {
                let cross=try VivoGaussianIntegralEngine.crossOverlap(leftSystem:request.snapshots[i-1].system,leftBasis:request.basis,
                    rightSystem:system,rightBasis:request.basis,budget:request.budget)
                overlap=try references[i-1].alphaCoefficients.transposed.multiplied(by:cross).multiplied(by:reference.alphaCoefficients)
            }
            points.append(.init(identifier:request.snapshots[i].identifier,hamiltonian:h,overlapWithPrevious:overlap))
        }
        try request.configuration.validate(points:points,budget:request.budget)
        return points
    }
    public static func hydrogenStretchTemplate() throws -> VivoMolecularECCPathRequest {
        let distances=[1.2,1.4,1.6,1.8,2.0]
        let snapshots=distances.enumerated().map { i,d in
            VivoMolecularPathSnapshot(identifier:"h2-\(i)",coordinate:d,system:.init(nuclei:[
                .init(atomicNumber:1,positionBohr:.init(0,0,-d/2),structureAtomIndex:0),
                .init(atomicNumber:1,positionBohr:.init(0,0,d/2),structureAtomIndex:1)],alphaElectrons:1,betaElectrons:1))
        }
        let embedding=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[
            .init(identifier:"bond",orbitals:[0],maximumBathOrbitals:1,clusterAlphaElectrons:1,clusterBetaElectrons:1)],
            bathSelection:.init(minimumBathOrbitals:1),localityGroups:[[0,1]],maximumOrbitalIterations:150,
            orbitalGradientTolerance:1e-5,maximumOrbitalStep:0.1)
        var generator=VivoQMMatrix(2,2); generator[0,1]=0.21; generator[1,0] = -0.21
        return .init(identifier:"H2-STO-3G-shared-orbital-stretch",atomIdentifiers:["H-left","H-right"],coordinateUnit:"Bohr",
            snapshots:snapshots,basis:.hydrogenSTO3G(nucleusIndices:[0,1]),
            configuration:.init(embedding:embedding,transportGroups:[[0],[1]],expectedBathOrbitals:[1],
                pointWeights:[Double](repeating:0.2,count:5)),
            initialSharedRotation:try VivoQMDenseAlgebra.orbitalRotation(generator:generator))
    }
}
