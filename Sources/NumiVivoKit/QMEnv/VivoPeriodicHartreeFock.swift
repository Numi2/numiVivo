import Foundation

public struct VivoPeriodicHartreeFockResult: Codable, Sendable, Equatable {
    public static let convention = "localized-HF+exact-Gaussian-near+distributed-raw-quadrupole-Ewald-images; no-MM-self; no-LJ; tin-foil-neutral-v1"
    public let sourceSystem: VivoElectronicSystem
    public let sourceBasis: VivoGaussianBasis
    public let periodicCell: VivoPeriodicCell
    public let configuration: VivoPeriodicQMMMConfiguration
    public let reference: VivoHartreeFockResult
    public let isolatedEnergyHartree: Double
    public let periodicEmbeddingEnergyHartree: Double
    public let nucleusForcesHartreePerBohr: [VivoVector3D]
    public let pointChargeForcesHartreePerBohr: [VivoVector3D]
    /// Complete electronic affine-strain derivative for independent nuclear and
    /// point-charge centers. Link/site constraints must apply their chain rules.
    public let affineStrainDerivativeHartree: VivoQMMatrix
    public let energyConvention: String
}

public enum VivoPeriodicHartreeFock {
    /// Extends the existing HF Fock/DIIS authority; no alternative SCF solver.
    /// Variational density coupling, primary/self subtractions, exact-near switches,
    /// AO moment response and Pulay derivatives all differentiate the same energy.
    public static func evaluate(system: VivoElectronicSystem,basis: VivoGaussianBasis,cell: VivoPeriodicCell,
                                configuration cfg: VivoPeriodicQMMMConfiguration = .init(),
                                budget: VivoChemistryBudget = .init()) throws -> VivoPeriodicHartreeFockResult {
        try cfg.validate();try budget.validate()
        var isolated = system;isolated.pointCharges = []
        let ao = try VivoGaussianIntegralEngine.compute(system: isolated,basis: basis,budget: budget)
        let context = try VivoPeriodicQMMMContext(source: system,integrals: ao,cell: cell,configuration: cfg,budget: budget)
        var previousDensity: VivoQMMatrix?,previousEvaluation: VivoPeriodicEmbeddingEvaluation?
        func embedding(_ da: VivoQMMatrix,_ db: VivoQMMatrix) throws -> VivoPeriodicEmbeddingEvaluation {
            let density = try da.adding(db)
            if previousDensity == density,let previousEvaluation { return previousEvaluation }
            let value = try context.evaluate(density: density,derivatives: false)
            previousDensity = density;previousEvaluation = value
            return value
        }
        let reference = try VivoHartreeFock.solveWithFockBuilder(system: isolated,overlap: ao.overlap,core: ao.coreHamiltonian,
            constant: ao.constantEnergyHartree,configuration: cfg.scf,budget: budget,
            energyFunctional: { da,db,_,_ in
                let gas = VivoHartreeFock.focks(ao,da,db)
                return try VivoHartreeFock.energy(ao,da,db,gas.0,gas.1)+embedding(da,db).energy
            },buildFocks: { da,db in
                let gas = VivoHartreeFock.focks(ao,da,db),potential = try embedding(da,db).fock
                return (try gas.0.adding(potential),try gas.1.adding(potential))
            })
        let density = try reference.alphaDensity.adding(reference.betaDensity)
        let correction = try context.evaluate(density: density,derivatives: true)
        let molecular = try VivoHartreeFockNuclearDerivatives.evaluate(integrals: ao,reference: reference,budget: budget)
        let gradients = zip(molecular.nuclearGradientsHartreePerBohr,correction.nuclearGradients).map(+)
        var strain = correction.affineStrain
        for i in system.nuclei.indices {
            let gradient = molecular.nuclearGradientsHartreePerBohr[i],g = [gradient.x,gradient.y,gradient.z]
            let r = system.nuclei[i].positionBohr
            for a in 0..<3 { for b in 0..<3 { strain[a,b] += g[a]*r[b] } }
        }
        let gasFock = VivoHartreeFock.focks(ao,reference.alphaDensity,reference.betaDensity)
        let gasEnergy = VivoHartreeFock.energy(ao,reference.alphaDensity,reference.betaDensity,gasFock.0,gasFock.1)
        guard gradients.allSatisfy(\.isFinite),strain.values.allSatisfy(\.isFinite),
              abs(gasEnergy+correction.energy-reference.energyHartree) <= 10*cfg.scf.energyToleranceHartree else {
            throw VivoChemistryError.convergence("periodic HF energy or analytic derivative reconstruction")
        }
        return .init(sourceSystem: system,sourceBasis: basis,periodicCell: cell,configuration: cfg,reference: reference,
            isolatedEnergyHartree: gasEnergy,periodicEmbeddingEnergyHartree: correction.energy,
            nucleusForcesHartreePerBohr: gradients.map { $0 * -1 },pointChargeForcesHartreePerBohr: correction.chargeGradients.map { $0 * -1 },
            affineStrainDerivativeHartree: strain,energyConvention: VivoPeriodicHartreeFockResult.convention)
    }
}
