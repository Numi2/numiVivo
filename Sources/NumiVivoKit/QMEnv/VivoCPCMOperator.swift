import Foundation

/// Immutable geometry-dependent C-PCM workspace. Gaussian surface integrals and
/// nuclear potentials are compiled once and reused throughout the SCF cycle.
/// The surface solve stays matrix-free FP64; no FP32 chemistry is substituted.
public struct VivoCPCMOperator: Sendable {
    public let configuration: VivoCPCMConfiguration
    public let tesserae: [VivoCPCMTessera]
    private let system: VivoElectronicSystem
    private let overlap: VivoQMMatrix
    private let aoPotentials: VivoQMMatrix // [tessera, flattened AO pair]
    private let nuclearPotential: [Double]
    private let budget: VivoChemistryBudget

    public init(system: VivoElectronicSystem, integrals ao: VivoAOIntegrals,
                configuration cfg: VivoCPCMConfiguration = .init(),
                budget: VivoChemistryBudget = .init()) throws {
        try system.validate(); try ao.validate(budget:budget); try cfg.validate(system:system)
        guard system == ao.sourceSystem else { throw VivoChemistryError.invalid("C-PCM operator AO/source binding") }
        let surface = try VivoCPCM.cavity(system:system,configuration:cfg.cavity)
        let count = try budget.elements([surface.count,ao.count,ao.count],simultaneousArrays:2)
        guard count <= budget.maximumOperatorApplications else {
            throw VivoChemistryError.resourceLimit("C-PCM compiled AO potential work")
        }
        var matrices = VivoQMMatrix(surface.count,ao.count*ao.count)
        var potential = [Double](repeating:0,count:surface.count)
        for i in surface.indices {
            let point = surface[i].positionBohr
            let matrix = try VivoGaussianElectrostaticPotential.matrix(integrals:ao,at:point,budget:budget)
            for pair in matrix.values.indices { matrices[i,pair] = matrix.values[pair] }
            for nucleus in system.nuclei {
                let distance = vivoQMNorm(point-nucleus.positionBohr)
                guard distance.isFinite, distance > 1e-10 else { throw VivoChemistryError.invalid("C-PCM tessera at a nucleus") }
                potential[i] += Double(nucleus.atomicNumber)/distance
            }
        }
        guard potential.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("C-PCM nuclear potential overflow") }
        self.system=system; self.overlap=ao.overlap; self.configuration=cfg; self.tesserae=surface
        self.aoPotentials=matrices; self.nuclearPotential=potential; self.budget=budget
    }

    public func evaluate(totalDensity density: VivoQMMatrix) throws -> VivoCPCMResult {
        let n = overlap.rows, cfg = configuration
        guard density.rows == n, density.columns == n, density.values.allSatisfy(\.isFinite),
              try density.adding(density.transposed,scale:-1).frobeniusNorm <= 1e-9 else {
            throw VivoChemistryError.invalid("C-PCM total density shape/symmetry")
        }
        let electrons = zip(density.values,overlap.values).reduce(0.0) { $0+$1.0*$1.1 }
        guard abs(electrons-Double(system.alphaElectrons+system.betaElectrons)) <= 1e-6 else {
            throw VivoChemistryError.invalid("C-PCM AO-metric density electron count")
        }
        let densitySpectrum = try VivoQMDenseAlgebra.symmetricEigen(density)
        guard densitySpectrum.values.allSatisfy({ $0 >= -1e-9 }) else {
            throw VivoChemistryError.invalid("C-PCM density is not positive semidefinite")
        }
        let electronic = try aoPotentials.multiplied(by:VivoQMMatrix(rows:n*n,columns:1,values:density.values))
        let potential = zip(nuclearPotential,electronic.values).map { $0-$1 }
        let scale = (cfg.dielectricConstant-1)/cfg.dielectricConstant, rhs = potential.map { -scale*$0 }
        let solution = try VivoCPCM.solve(rhs:rhs,tesserae:tesserae,tolerance:cfg.solverTolerance,
            maximumIterations:cfg.maximumSolverIterations,budget:budget)
        let charges = solution.0, applied = try VivoCPCM.applySurfaceMatrix(charges,tesserae)
        // Recompute the true residual; a small recursive PCG residual alone is
        // not sufficient evidence that the original surface equation was solved.
        let residual = zip(applied,rhs).reduce(0.0) { hypot($0,$1.0-$1.1) }
        let rhsNorm = rhs.reduce(0.0) { hypot($0,$1) }, relative = rhsNorm > 0 ? residual/rhsNorm : residual
        guard relative.isFinite, relative <= 8*cfg.solverTolerance else {
            throw VivoChemistryError.convergence("C-PCM true surface-equation residual exceeds tolerance")
        }
        let energy = 0.5*zip(charges,potential).reduce(0.0) { $0+$1.0*$1.1 }
        let row = try VivoQMMatrix(rows:1,columns:charges.count,values:charges.map { -$0 })
        let reaction = try VivoQMMatrix(rows:n,columns:n,values:row.multiplied(by:aoPotentials).values)
        let soluteCharge = Double(system.nuclei.reduce(0) { $0+$1.atomicNumber })-electrons
        let expected = -scale*soluteCharge, total = charges.reduce(0,+)
        let gauss = abs(total-expected)/max(1,abs(expected))
        guard gauss.isFinite, gauss <= cfg.maximumGaussLawRelativeError,
              energy.isFinite, reaction.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("C-PCM Gauss-law or finite-value diagnostic")
        }
        return .init(dielectricConstant:cfg.dielectricConstant,dielectricScale:scale,tesserae:tesserae,
            solutePotentialHartreePerE:potential,apparentSurfaceChargesE:charges,polarizationEnergyHartree:energy,
            reactionPotentialMatrix:reaction,linearResidualNorm:residual,linearRelativeResidual:relative,
            solverIterations:solution.2,totalSurfaceChargeE:total,expectedGaussLawSurfaceChargeE:expected,
            gaussLawRelativeError:gauss,
            status:"C-PCM x=0 point-charge cavity; cached FP64 AO potentials; true PCG residual; electrostatic polarization only, not a complete Gibbs energy")
    }
}
