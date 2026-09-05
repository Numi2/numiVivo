import Foundation

public struct VivoLDAConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var energyToleranceHartree: Double
    public var densityTolerance: Double
    public var commutatorTolerance: Double
    public var minimumOverlapEigenvalue: Double
    public var diisHistory: Int
    public var initialDamping: Double
    public var maximumIntegratedElectronError: Double
    public var grid: VivoDFTGridConfiguration
    /// Absent in older documents: nil retains the gas-phase LDA method.
    public var cpcm: VivoCPCMConfiguration?
    public init(maximumIterations: Int = 128, energyToleranceHartree: Double = 1e-9,
                densityTolerance: Double = 1e-7, commutatorTolerance: Double = 1e-7,
                minimumOverlapEigenvalue: Double = 1e-8, diisHistory: Int = 8,
                initialDamping: Double = 0.2, maximumIntegratedElectronError: Double = 5e-4,
                grid: VivoDFTGridConfiguration = .init(), cpcm: VivoCPCMConfiguration? = nil) {
        self.maximumIterations=maximumIterations; self.energyToleranceHartree=energyToleranceHartree
        self.densityTolerance=densityTolerance; self.commutatorTolerance=commutatorTolerance
        self.minimumOverlapEigenvalue=minimumOverlapEigenvalue; self.diisHistory=diisHistory
        self.initialDamping=initialDamping; self.maximumIntegratedElectronError=maximumIntegratedElectronError
        self.grid=grid; self.cpcm=cpcm
    }
    public func validate(atomCount: Int) throws {
        try grid.validate(atomCount:atomCount)
        guard (1...10000).contains(maximumIterations), (2...32).contains(diisHistory),
              [energyToleranceHartree,densityTolerance,commutatorTolerance,minimumOverlapEigenvalue,maximumIntegratedElectronError].allSatisfy({$0.isFinite && $0>0}),
              initialDamping.isFinite, (0..<1).contains(initialDamping) else { throw VivoChemistryError.invalid("LDA configuration") }
    }
}
public struct VivoLDAIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    /// Zero at iteration one means no previous density comparison; convergence
    /// also requires a previous energy and cannot be declared at that iteration.
    public let densityChange: Double
    public let commutatorNorm: Double
    public let integratedElectrons: Double
    public let exchangeCorrelationEnergyHartree: Double
    public let polarizationEnergyHartree: Double?
}
public struct VivoLDAResult: Codable, Sendable, Equatable {
    public let functional: String
    /// Electronic energy functional plus optional electrostatic polarization.
    /// This is NOT a thermally/standard-state corrected Gibbs energy.
    public let energyHartree: Double
    public let coefficients: VivoQMMatrix
    public let orbitalEnergiesHartree: [Double]
    public let totalDensity: VivoQMMatrix
    public let exchangeCorrelationEnergyHartree: Double
    public let integratedElectrons: Double
    public let integratedElectronError: Double
    public let finalCommutatorNorm: Double
    public let gridPointCount: Int
    public let iterations: [VivoLDAIteration]
    public let status: String
    public let gasEnergyFunctionalHartree: Double?
    public let cpcm: VivoCPCMResult?
}
private struct VivoLDADIISRecord { let fock: VivoQMMatrix; let error: [Double] }

/// Geometry-fixed AO grid values, reused across every density evaluation. The
/// density projection and XC matrix contractions use FP64 GEMM/Accelerate.
private struct VivoLDAWorkspace {
    let phi: VivoQMMatrix
    let weights: [Double]
    init(ao: VivoAOIntegrals, grid: [VivoDFTGridPoint], budget: VivoChemistryBudget) throws {
        _ = try budget.elements([grid.count,ao.count],simultaneousArrays:6)
        var values=VivoQMMatrix(grid.count,ao.count)
        for (i,point) in grid.enumerated() {
            for (mu,orbital) in ao.orbitals.enumerated() {
                let d=point.position-ao.sourceSystem.nuclei[orbital.nucleusIndex].positionBohr
                let polynomial=pow(d.x,Double(orbital.angular[0]))*pow(d.y,Double(orbital.angular[1]))*pow(d.z,Double(orbital.angular[2]))
                let r2=vivoQMDot(d,d)
                var radial=0.0
                for k in orbital.weights.indices { radial += orbital.weights[k]*exp(-orbital.primitiveExponents[k]*r2) }
                values[i,mu]=polynomial*radial
            }
        }
        guard values.values.allSatisfy(\.isFinite), grid.allSatisfy({$0.weightBohr3.isFinite && $0.weightBohr3>=0}) else {
            throw VivoChemistryError.invalid("LDA AO grid values/weights")
        }
        phi=values; weights=grid.map(\.weightBohr3)
    }
    func evaluate(_ density: VivoQMMatrix) throws -> (matrix:VivoQMMatrix,energy:Double,electrons:Double) {
        let projected=try phi.multiplied(by:density)
        var weighted=phi, energy=0.0, electrons=0.0
        for i in 0..<phi.rows {
            var rho=0.0
            for mu in 0..<phi.columns { rho += phi[i,mu]*projected[i,mu] }
            if rho<0 && rho > -1e-10 { rho=0 }
            guard rho.isFinite, rho>=0 else { throw VivoChemistryError.convergence("negative/nonfinite LDA density") }
            let local=try VivoPZ81LDA.evaluate(density:rho), weight=weights[i]
            electrons += weight*rho; energy += weight*rho*local.xcEnergyPerElectronHartree
            for mu in 0..<phi.columns { weighted[i,mu] *= weight*local.xcPotentialHartree }
        }
        let potential=try phi.transposed.multiplied(by:weighted)
        guard energy.isFinite, electrons.isFinite else { throw VivoChemistryError.convergence("LDA grid reduction") }
        return (potential,energy,electrons)
    }
}
public enum VivoRestrictedLDA {
    private static func density(_ c:VivoQMMatrix, occupied:Int) -> VivoQMMatrix {
        VivoHartreeFock.density(c,occupied:occupied).scaled(2)
    }
    private static func coulomb(_ ao:VivoAOIntegrals, density p:VivoQMMatrix) -> VivoQMMatrix {
        let n=ao.count
        var j=VivoQMMatrix(n,n)
        for mu in 0..<n { for nu in 0..<n { for r in 0..<n { for s in 0..<n { j[mu,nu] += p[r,s]*ao.eri(mu,nu,r,s) } } } }
        return j
    }
    public static func solve(system:VivoElectronicSystem, integrals ao:VivoAOIntegrals,
                             configuration cfg:VivoLDAConfiguration = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoLDAResult {
        try system.validate(); try ao.validate(budget:budget); try cfg.validate(atomCount:system.nuclei.count)
        guard system==ao.sourceSystem, system.alphaElectrons==system.betaElectrons, system.alphaElectrons<=ao.count else {
            throw VivoChemistryError.unsupported("current LDA authority is closed-shell spin-unpolarized only")
        }
        let n=ao.count, occupied=system.alphaElectrons, electronCount=2*occupied
        _ = try budget.elements([n,n],simultaneousArrays:30+2*cfg.diisHistory)
        let overlapEigen=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        guard overlapEigen.values.first!>=cfg.minimumOverlapEigenvalue else { throw VivoChemistryError.invalid("LDA AO linear dependence") }
        var x=overlapEigen.vectors
        for i in 0..<n { for j in 0..<n { x[i,j]/=sqrt(overlapEigen.values[j]) } }
        func canonical(_ f:VivoQMMatrix) throws -> (VivoQMMatrix,[Double]) {
            let e=try VivoQMDenseAlgebra.symmetricEigen(f.congruence(x))
            return (try x.multiplied(by:e.vectors),e.values)
        }
        let grid=try VivoDFTQuadrature.build(system:system,configuration:cfg.grid)
        let workspace=try VivoLDAWorkspace(ao:ao,grid:grid,budget:budget)
        let solvent=try cfg.cpcm.map { try VivoCPCMOperator(system:system,integrals:ao,configuration:$0,budget:budget) }
        func evaluate(_ p:VivoQMMatrix) throws -> (fock:VivoQMMatrix,energy:Double,gas:Double,xc:Double,electrons:Double,pcm:VivoCPCMResult?) {
            let j=coulomb(ao,density:p), xc=try workspace.evaluate(p), pcm=try solvent?.evaluate(totalDensity:p)
            var fock=try ao.coreHamiltonian.adding(j).adding(xc.matrix)
            if let pcm { fock=try fock.adding(pcm.reactionPotentialMatrix) }
            let gas=ao.constantEnergyHartree
                + zip(p.values,ao.coreHamiltonian.values).reduce(0.0) { $0+$1.0*$1.1 }
                + 0.5*zip(p.values,j.values).reduce(0.0) { $0+$1.0*$1.1 }+xc.energy
            let energy=gas+(pcm?.polarizationEnergyHartree ?? 0)
            guard energy.isFinite else { throw VivoChemistryError.convergence("nonfinite LDA energy") }
            return (fock,energy,gas,xc.energy,xc.electrons,pcm)
        }
        let (initial,_)=try canonical(ao.coreHamiltonian)
        var p=density(initial,occupied:occupied), previous:Double?, densityChange=0.0
        var history:[VivoLDADIISRecord]=[], trace:[VivoLDAIteration]=[]
        for iteration in 1...cfg.maximumIterations {
            let current=try evaluate(p)
            let error=try VivoHartreeFock.error(current.fock,p,ao.overlap,x)
            let residual=error.reduce(0.0) { hypot($0,$1) }
            trace.append(.init(iteration:iteration,energyHartree:current.energy,densityChange:densityChange,
                commutatorNorm:residual,integratedElectrons:current.electrons,exchangeCorrelationEnergyHartree:current.xc,
                polarizationEnergyHartree:current.pcm?.polarizationEnergyHartree))
            if let previous, abs(current.energy-previous)<=cfg.energyToleranceHartree,
               densityChange<=cfg.densityTolerance, residual<=cfg.commutatorTolerance {
                let (coefficients,eps)=try canonical(current.fock), finalDensity=density(coefficients,occupied:occupied)
                let final=try evaluate(finalDensity)
                let finalResidual=try VivoHartreeFock.error(final.fock,finalDensity,ao.overlap,x).reduce(0.0) { hypot($0,$1) }
                let electronError=abs(final.electrons-Double(electronCount))
                guard electronError<=cfg.maximumIntegratedElectronError else { throw VivoChemistryError.convergence("LDA grid electron-count error exceeds configured bound") }
                guard finalResidual<=cfg.commutatorTolerance,
                      try finalDensity.adding(p,scale:-1).frobeniusNorm<=cfg.densityTolerance,
                      abs(final.energy-current.energy)<=cfg.energyToleranceHartree else {
                    throw VivoChemistryError.convergence("final LDA canonicalization violates SCF tolerances")
                }
                return .init(functional:"spin-unpolarized Dirac exchange + Perdew-Zunger 1981 correlation",
                    energyHartree:final.energy,coefficients:coefficients,orbitalEnergiesHartree:eps,totalDensity:finalDensity,
                    exchangeCorrelationEnergyHartree:final.xc,integratedElectrons:final.electrons,integratedElectronError:electronError,
                    finalCommutatorNorm:finalResidual,gridPointCount:grid.count,iterations:trace,
                    status:"native FP64 LDA with cached AO grid; optional self-consistent C-PCM; electronic plus electrostatic polarization, not a complete Gibbs energy",
                    gasEnergyFunctionalHartree:final.gas,cpcm:final.pcm)
            }
            history.append(.init(fock:current.fock,error:error))
            if history.count>cfg.diisHistory { history.removeFirst() }
            var guess=current.fock
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
            var next=density(coefficients,occupied:occupied)
            if iteration<=3 { next=try next.scaled(1-cfg.initialDamping).adding(p,scale:cfg.initialDamping) }
            densityChange=try next.adding(p,scale:-1).frobeniusNorm
            previous=current.energy; p=next
        }
        throw VivoChemistryError.convergence("LDA KS-SCF did not converge")
    }
}
