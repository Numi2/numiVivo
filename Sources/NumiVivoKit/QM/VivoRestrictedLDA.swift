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
    public init(maximumIterations: Int = 128, energyToleranceHartree: Double = 1e-9,
                densityTolerance: Double = 1e-7, commutatorTolerance: Double = 1e-7,
                minimumOverlapEigenvalue: Double = 1e-8, diisHistory: Int = 8,
                initialDamping: Double = 0.2, maximumIntegratedElectronError: Double = 5e-4,
                grid: VivoDFTGridConfiguration = .init()) {
        self.maximumIterations=maximumIterations;self.energyToleranceHartree=energyToleranceHartree
        self.densityTolerance=densityTolerance;self.commutatorTolerance=commutatorTolerance
        self.minimumOverlapEigenvalue=minimumOverlapEigenvalue;self.diisHistory=diisHistory
        self.initialDamping=initialDamping;self.maximumIntegratedElectronError=maximumIntegratedElectronError;self.grid=grid
    }
    public func validate(atomCount: Int) throws {
        try grid.validate(atomCount: atomCount)
        guard maximumIterations>0,(2...32).contains(diisHistory),
              [energyToleranceHartree,densityTolerance,commutatorTolerance,minimumOverlapEigenvalue,maximumIntegratedElectronError].allSatisfy({$0.isFinite&&$0>0}),
              initialDamping.isFinite,(0..<1).contains(initialDamping) else { throw VivoChemistryError.invalid("LDA configuration") }
    }
}

public struct VivoLDAIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let densityChange: Double
    public let commutatorNorm: Double
    public let integratedElectrons: Double
    public let exchangeCorrelationEnergyHartree: Double
}

public struct VivoLDAResult: Codable, Sendable, Equatable {
    public let functional: String
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
}

private struct VivoLDADIISRecord { let fock: VivoQMMatrix; let error: [Double] }

public enum VivoRestrictedLDA {
    private static func density(_ c: VivoQMMatrix, occupied: Int) -> VivoQMMatrix {
        var d=VivoQMMatrix(c.rows,c.rows)
        for p in 0..<c.rows { for q in 0..<c.rows { for i in 0..<occupied { d[p,q] += 2*c[p,i]*c[q,i] } } }
        return d
    }

    private static func aoValues(_ orbitals: [VivoCartesianOrbital], system: VivoElectronicSystem,
                                 at point: SIMD3<Double>) -> [Double] {
        var values=[Double](repeating:0,count:orbitals.count)
        for mu in orbitals.indices {
            let orbital=orbitals[mu],center=system.nuclei[orbital.nucleusIndex].positionBohr,d=point-center
            let polynomial=pow(d.x,Double(orbital.angular[0]))*pow(d.y,Double(orbital.angular[1]))*pow(d.z,Double(orbital.angular[2]))
            let r2=vivoQMDot(d,d);var radial=0.0
            for i in orbital.weights.indices { radial += orbital.weights[i]*exp(-orbital.primitiveExponents[i]*r2) }
            values[mu]=polynomial*radial
        }
        return values
    }

    private static func coulomb(_ ao: VivoAOIntegrals, density p: VivoQMMatrix) -> VivoQMMatrix {
        let n=ao.count;var j=VivoQMMatrix(n,n)
        for mu in 0..<n { for nu in 0..<n { for r in 0..<n { for s in 0..<n { j[mu,nu] += p[r,s]*ao.eri(mu,nu,r,s) } } } }
        return j
    }

    private static func xc(density p: VivoQMMatrix, ao: VivoAOIntegrals,
                           grid: [VivoDFTGridPoint]) throws -> (matrix: VivoQMMatrix, energy: Double, electrons: Double) {
        let n=ao.count;var potential=VivoQMMatrix(n,n),energy=0.0,electrons=0.0
        for point in grid {
            let phi=aoValues(ao.orbitals,system:ao.sourceSystem,at:point.position);var rho=0.0
            for mu in 0..<n { for nu in 0..<n { rho += p[mu,nu]*phi[mu]*phi[nu] } }
            if rho < 0 && rho > -1e-10 { rho=0 }
            guard rho.isFinite,rho>=0 else { throw VivoChemistryError.convergence("negative/nonfinite LDA density") }
            let local=try VivoPZ81LDA.evaluate(density:rho),weight=point.weightBohr3
            electrons += weight*rho;energy += weight*rho*local.xcEnergyPerElectronHartree
            let coefficient=weight*local.xcPotentialHartree
            if coefficient != 0 {
                for mu in 0..<n { for nu in 0...mu {
                    let value=coefficient*phi[mu]*phi[nu];potential[mu,nu]+=value
                    if mu != nu { potential[nu,mu]+=value }
                } }
            }
        }
        guard energy.isFinite,electrons.isFinite else { throw VivoChemistryError.convergence("LDA grid reduction") }
        return(potential,energy,electrons)
    }

    private static func error(_ f:VivoQMMatrix,_ d:VivoQMMatrix,_ s:VivoQMMatrix,_ x:VivoQMMatrix)throws->[Double]{
        let fds=try f.multiplied(by:d).multiplied(by:s),sdf=try s.multiplied(by:d).multiplied(by:f)
        return try fds.adding(sdf,scale:-1).congruence(x).values
    }

    public static func solve(system: VivoElectronicSystem, integrals ao: VivoAOIntegrals,
                             configuration cfg: VivoLDAConfiguration = .init(),
                             budget: VivoChemistryBudget = .init()) throws -> VivoLDAResult {
        try system.validate();try ao.validate(budget:budget);try cfg.validate(atomCount:system.nuclei.count)
        guard system==ao.sourceSystem,system.alphaElectrons==system.betaElectrons,system.alphaElectrons<=ao.count else {
            throw VivoChemistryError.unsupported("current LDA authority is closed-shell spin-unpolarized only")
        }
        let n=ao.count,occupied=system.alphaElectrons,electronCount=2*occupied
        _=try budget.elements([n,n],simultaneousArrays:30+2*cfg.diisHistory)
        let overlapEigen=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        guard overlapEigen.values.first!>=cfg.minimumOverlapEigenvalue else { throw VivoChemistryError.invalid("LDA AO linear dependence") }
        var x=overlapEigen.vectors
        for i in 0..<n { for j in 0..<n { x[i,j]/=sqrt(overlapEigen.values[j]) } }
        func canonical(_ f:VivoQMMatrix)throws->(VivoQMMatrix,[Double]){
            let e=try VivoQMDenseAlgebra.symmetricEigen(f.congruence(x));return(try x.multiplied(by:e.vectors),e.values)
        }
        let grid=try VivoDFTQuadrature.build(system:system,configuration:cfg.grid)
        let (initialCoefficients,_)=try canonical(ao.coreHamiltonian)
        var densityMatrix=density(initialCoefficients,occupied:occupied),previousEnergy:Double?
        var history:[VivoLDADIISRecord]=[],trace:[VivoLDAIteration]=[],densityChange=Double.infinity
        for iteration in 1...cfg.maximumIterations {
            let j=coulomb(ao,density:densityMatrix),xcResult=try xc(density:densityMatrix,ao:ao,grid:grid)
            let fock=try ao.coreHamiltonian.adding(j).adding(xcResult.matrix)
            let energy=ao.constantEnergyHartree
                + zip(densityMatrix.values,ao.coreHamiltonian.values).reduce(0){$0+$1.0*$1.1}
                + 0.5*zip(densityMatrix.values,j.values).reduce(0){$0+$1.0*$1.1}+xcResult.energy
            let commutator=try error(fock,densityMatrix,ao.overlap,x)
            let residual=sqrt(commutator.reduce(0){$0+$1*$1})
            trace.append(.init(iteration:iteration,energyHartree:energy,densityChange:densityChange,
                               commutatorNorm:residual,integratedElectrons:xcResult.electrons,
                               exchangeCorrelationEnergyHartree:xcResult.energy))
            if let last=previousEnergy,abs(energy-last)<=cfg.energyToleranceHartree,
               densityChange<=cfg.densityTolerance,residual<=cfg.commutatorTolerance {
                let (coefficients,orbitalEnergies)=try canonical(fock)
                let finalDensity=density(coefficients,occupied:occupied),finalJ=coulomb(ao,density:finalDensity)
                let finalXC=try xc(density:finalDensity,ao:ao,grid:grid)
                let finalFock=try ao.coreHamiltonian.adding(finalJ).adding(finalXC.matrix)
                let finalEnergy=ao.constantEnergyHartree
                    + zip(finalDensity.values,ao.coreHamiltonian.values).reduce(0){$0+$1.0*$1.1}
                    + 0.5*zip(finalDensity.values,finalJ.values).reduce(0){$0+$1.0*$1.1}+finalXC.energy
                let finalResidual=sqrt(try error(finalFock,finalDensity,ao.overlap,x).reduce(0){$0+$1*$1})
                let electronError=abs(finalXC.electrons-Double(electronCount))
                guard electronError<=cfg.maximumIntegratedElectronError else {
                    throw VivoChemistryError.convergence("LDA grid integrates \(finalXC.electrons) electrons; error \(electronError) exceeds configured bound")
                }
                guard finalResidual<=cfg.commutatorTolerance,
                      try finalDensity.adding(densityMatrix,scale:-1).frobeniusNorm<=cfg.densityTolerance else {
                    throw VivoChemistryError.convergence("final LDA canonicalization violates SCF tolerances")
                }
                return .init(functional:"spin-unpolarized Dirac exchange + Perdew-Zunger 1981 correlation",
                             energyHartree:finalEnergy,coefficients:coefficients,orbitalEnergiesHartree:orbitalEnergies,
                             totalDensity:finalDensity,exchangeCorrelationEnergyHartree:finalXC.energy,
                             integratedElectrons:finalXC.electrons,integratedElectronError:electronError,
                             finalCommutatorNorm:finalResidual,gridPointCount:grid.count,iterations:trace,
                             status:"native FP64 LDA; smooth atom partition + mapped Gauss-Legendre radial/Fibonacci angular quadrature; grid diagnostics authoritative")
            }
            history.append(.init(fock:fock,error:commutator));if history.count>cfg.diisHistory{history.removeFirst()}
            var guess=fock
            if history.count>=2 {
                let m=history.count;var b=VivoQMMatrix(m+1,m+1),rhs=[Double](repeating:0,count:m+1);rhs[m] = -1
                let scale=max(1e-30,history.map{$0.error.reduce(0){$0+$1*$1}}.max() ?? 1)
                for i in 0..<m { for k in 0..<m {
                    b[i,k]=zip(history[i].error,history[k].error).reduce(0){$0+$1.0*$1.1}/scale
                };b[i,m] = -1;b[m,i] = -1 }
                if let weights=try? VivoQMDenseAlgebra.solve(b,rhs:rhs),weights.allSatisfy({abs($0)<1e6}) {
                    guess=VivoQMMatrix(n,n);for i in 0..<m{guess=try guess.adding(history[i].fock,scale:weights[i])}
                } else if history.count>2 { history.removeFirst() }
            }
            let (coefficients,_)=try canonical(guess);var next=density(coefficients,occupied:occupied)
            if iteration<=3 { next=try next.scaled(1-cfg.initialDamping).adding(densityMatrix,scale:cfg.initialDamping) }
            densityChange=try next.adding(densityMatrix,scale:-1).frobeniusNorm
            previousEnergy=energy;densityMatrix=next
        }
        throw VivoChemistryError.convergence("LDA KS-SCF did not converge")
    }
}
