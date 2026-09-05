import Foundation

public struct VivoCPCMRHFConfiguration: Codable, Sendable, Equatable {
    public var scf: VivoSCFConfiguration
    public var cpcm: VivoCPCMConfiguration
    public init(scf: VivoSCFConfiguration = .init(reference:.restricted), cpcm: VivoCPCMConfiguration = .init()) {
        self.scf=scf; self.cpcm=cpcm
    }
    public func validate(system:VivoElectronicSystem) throws {
        try scf.validate(); try cpcm.validate(system:system)
        guard scf.reference == .restricted, system.alphaElectrons==system.betaElectrons else {
            throw VivoChemistryError.unsupported("current variational C-PCM SCF is restricted closed-shell")
        }
    }
}
public struct VivoCPCMRHFIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let totalFreeEnergyHartree: Double
    public let gasEnergyFunctionalHartree: Double
    public let polarizationEnergyHartree: Double
    /// Zero in the first iteration indicates no previous density comparison.
    public let densityChange: Double
    public let commutatorNorm: Double
    public let cpcmRelativeResidual: Double
}
public struct VivoCPCMRHFResult: Codable, Sendable, Equatable {
    /// Legacy field name retained: electronic energy plus electrostatic solvent
    /// polarization, NOT a complete standard-state/thermal Gibbs free energy.
    public let totalFreeEnergyHartree: Double
    public let gasEnergyFunctionalHartree: Double
    public let polarizationEnergyHartree: Double
    public let coefficients: VivoQMMatrix
    public let orbitalEnergiesHartree: [Double]
    public let alphaDensity: VivoQMMatrix
    public let betaDensity: VivoQMMatrix
    public let cpcm: VivoCPCMResult
    public let finalCommutatorNorm: Double
    public let iterations: [VivoCPCMRHFIteration]
    public let status: String
}
private struct VivoCPCMDIISRecord { let fock:VivoQMMatrix; let error:[Double] }
public enum VivoCPCMRHF {
    public static func solve(system:VivoElectronicSystem, integrals ao:VivoAOIntegrals,
                             configuration cfg:VivoCPCMRHFConfiguration = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoCPCMRHFResult {
        try system.validate(); try ao.validate(budget:budget); try cfg.validate(system:system)
        guard system==ao.sourceSystem, system.alphaElectrons<=ao.count else { throw VivoChemistryError.invalid("C-PCM RHF source binding/occupation") }
        let n=ao.count, occupied=system.alphaElectrons
        _ = try budget.elements([n,n],simultaneousArrays:40+2*cfg.scf.diisHistory)
        let overlapEigen=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        guard overlapEigen.values.first!>=cfg.scf.minimumOverlapEigenvalue else { throw VivoChemistryError.invalid("C-PCM RHF AO linear dependence") }
        var x=overlapEigen.vectors
        for i in 0..<n { for j in 0..<n { x[i,j]/=sqrt(overlapEigen.values[j]) } }
        func canonical(_ f:VivoQMMatrix) throws -> (VivoQMMatrix,[Double]) {
            let e=try VivoQMDenseAlgebra.symmetricEigen(f.congruence(x)); return (try x.multiplied(by:e.vectors),e.values)
        }
        let solvent=try VivoCPCMOperator(system:system,integrals:ao,configuration:cfg.cpcm,budget:budget)
        let (initial,_)=try canonical(ao.coreHamiltonian)
        var da=VivoHartreeFock.density(initial,occupied:occupied), db=da
        var history:[VivoCPCMDIISRecord]=[], trace:[VivoCPCMRHFIteration]=[], previous:Double?, densityChange=0.0
        for iteration in 1...cfg.scf.maximumIterations {
            let gas=VivoHartreeFock.focks(ao,da,db), totalDensity=try da.adding(db)
            let pcm=try solvent.evaluate(totalDensity:totalDensity)
            let fa=try gas.0.adding(pcm.reactionPotentialMatrix), fb=try gas.1.adding(pcm.reactionPotentialMatrix)
            let gasEnergy=VivoHartreeFock.energy(ao,da,db,gas.0,gas.1), free=gasEnergy+pcm.polarizationEnergyHartree
            let errors=try VivoHartreeFock.error(fa,da,ao.overlap,x)+VivoHartreeFock.error(fb,db,ao.overlap,x)
            let residual=errors.reduce(0.0) { hypot($0,$1) }
            guard free.isFinite, residual.isFinite else { throw VivoChemistryError.convergence("nonfinite C-PCM RHF iteration") }
            trace.append(.init(iteration:iteration,totalFreeEnergyHartree:free,gasEnergyFunctionalHartree:gasEnergy,
                polarizationEnergyHartree:pcm.polarizationEnergyHartree,densityChange:densityChange,
                commutatorNorm:residual,cpcmRelativeResidual:pcm.linearRelativeResidual))
            if let last=previous, abs(free-last)<=cfg.scf.energyToleranceHartree,
               densityChange<=cfg.scf.densityTolerance, residual<=cfg.scf.commutatorTolerance {
                let (coefficients,orbitalEnergies)=try canonical(fa)
                let finalDA=VivoHartreeFock.density(coefficients,occupied:occupied), finalDB=finalDA
                let finalGas=VivoHartreeFock.focks(ao,finalDA,finalDB)
                let finalPCM=try solvent.evaluate(totalDensity:finalDA.adding(finalDB))
                let finalFock=try finalGas.0.adding(finalPCM.reactionPotentialMatrix)
                let finalErrors=try VivoHartreeFock.error(finalFock,finalDA,ao.overlap,x)+VivoHartreeFock.error(finalFock,finalDB,ao.overlap,x)
                let finalResidual=finalErrors.reduce(0.0) { hypot($0,$1) }
                let change=try finalDA.adding(da,scale:-1).frobeniusNorm+finalDB.adding(db,scale:-1).frobeniusNorm
                let finalGasEnergy=VivoHartreeFock.energy(ao,finalDA,finalDB,finalGas.0,finalGas.1)
                let finalFree=finalGasEnergy+finalPCM.polarizationEnergyHartree
                guard finalResidual<=cfg.scf.commutatorTolerance, change<=cfg.scf.densityTolerance,
                      abs(finalFree-free)<=cfg.scf.energyToleranceHartree else {
                    throw VivoChemistryError.convergence("final C-PCM RHF canonicalization violates SCF tolerances")
                }
                return .init(totalFreeEnergyHartree:finalFree,gasEnergyFunctionalHartree:finalGasEnergy,
                    polarizationEnergyHartree:finalPCM.polarizationEnergyHartree,coefficients:coefficients,
                    orbitalEnergiesHartree:orbitalEnergies,alphaDensity:finalDA,betaDensity:finalDB,cpcm:finalPCM,
                    finalCommutatorNorm:finalResidual,iterations:trace,
                    status:"restricted HF + self-consistent cached C-PCM; electrostatic solvation energy only, not a complete Gibbs energy")
            }
            history.append(.init(fock:fa,error:errors))
            if history.count>cfg.scf.diisHistory { history.removeFirst() }
            var guess=fa
            if history.count>=2 {
                let m=history.count
                var b=VivoQMMatrix(m+1,m+1), rhs=[Double](repeating:0,count:m+1); rhs[m] = -1
                let scale=max(1e-30,history.map { $0.error.reduce(0) { $0+$1*$1 } }.max() ?? 1)
                for i in 0..<m { for k in 0..<m { b[i,k]=zip(history[i].error,history[k].error).reduce(0) { $0+$1.0*$1.1 }/scale }; b[i,m] = -1; b[m,i] = -1 }
                if let weights=try? VivoQMDenseAlgebra.solve(b,rhs:rhs), weights.allSatisfy({abs($0)<1e6}) {
                    guess=VivoQMMatrix(n,n)
                    for i in 0..<m { guess=try guess.adding(history[i].fock,scale:weights[i]) }
                } else if history.count>2 { history.removeFirst() }
            }
            let (coefficients,_)=try canonical(guess)
            var next=VivoHartreeFock.density(coefficients,occupied:occupied)
            if iteration<=3 { next=try next.scaled(1-cfg.scf.initialDamping).adding(da,scale:cfg.scf.initialDamping) }
            densityChange=2*(try next.adding(da,scale:-1).frobeniusNorm)
            previous=free; da=next; db=next
        }
        throw VivoChemistryError.convergence("C-PCM RHF did not converge")
    }
}
