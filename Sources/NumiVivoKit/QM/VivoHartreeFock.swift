import Foundation

public enum VivoHFReference: String, Codable, Sendable { case restricted, unrestricted }
public struct VivoSCFConfiguration: Codable, Sendable, Equatable {
    public var reference: VivoHFReference
    public var maximumIterations: Int
    public var energyToleranceHartree: Double
    public var densityTolerance: Double
    public var commutatorTolerance: Double
    public var minimumOverlapEigenvalue: Double
    public var diisHistory: Int
    public var initialDamping: Double
    public init(reference: VivoHFReference = .restricted, maximumIterations: Int = 128,
                energyToleranceHartree: Double = 1e-10, densityTolerance: Double = 1e-8,
                commutatorTolerance: Double = 1e-8, minimumOverlapEigenvalue: Double = 1e-8,
                diisHistory: Int = 8, initialDamping: Double = 0.2) {
        self.reference=reference; self.maximumIterations=maximumIterations
        self.energyToleranceHartree=energyToleranceHartree; self.densityTolerance=densityTolerance
        self.commutatorTolerance=commutatorTolerance; self.minimumOverlapEigenvalue=minimumOverlapEigenvalue
        self.diisHistory=diisHistory; self.initialDamping=initialDamping
    }
    public func validate() throws {
        guard (1...10000).contains(maximumIterations), (2...32).contains(diisHistory),
              [energyToleranceHartree,densityTolerance,commutatorTolerance,minimumOverlapEigenvalue].allSatisfy({$0.isFinite && $0>0}),
              initialDamping.isFinite, (0..<1).contains(initialDamping) else { throw VivoChemistryError.invalid("SCF configuration") }
    }
}
public struct VivoSCFIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let commutatorNorm: Double
    public let densityChange: Double
}
public struct VivoHartreeFockResult: Codable, Sendable, Equatable {
    public let reference: VivoHFReference
    public let alphaElectrons: Int
    public let betaElectrons: Int
    public let energyHartree: Double
    public let alphaCoefficients: VivoQMMatrix
    public let betaCoefficients: VivoQMMatrix
    public let alphaOrbitalEnergies: [Double]
    public let betaOrbitalEnergies: [Double]
    /// Spin-resolved densities: each occupied spin orbital contributes one.
    public let alphaDensity: VivoQMMatrix
    public let betaDensity: VivoQMMatrix
    public let finalCommutatorNorm: Double
    public let spinSquared: Double
    public let iterations: [VivoSCFIteration]
}
private struct VivoDIISRecord {
    let alphaFock: VivoQMMatrix
    let betaFock: VivoQMMatrix
    let error: [Double]
}
public enum VivoHartreeFock {
    static func density(_ c: VivoQMMatrix, occupied: Int) -> VivoQMMatrix {
        var d = VivoQMMatrix(c.rows,c.rows)
        for p in 0..<c.rows { for q in 0..<c.rows { for i in 0..<occupied { d[p,q] += c[p,i]*c[q,i] } } }; return d
    }
    static func focks(_ ao: VivoAOIntegrals, _ da: VivoQMMatrix, _ db: VivoQMMatrix) -> (VivoQMMatrix,VivoQMMatrix) {
        let n=ao.count; var fa=ao.coreHamiltonian,fb=fa
        for p in 0..<n { for q in 0..<n {
            var j=0.0,ka=0.0,kb=0.0
            for r in 0..<n { for s in 0..<n {
                j += (da[r,s]+db[r,s])*ao.eri(p,q,r,s)
                ka += da[r,s]*ao.eri(p,r,q,s); kb += db[r,s]*ao.eri(p,r,q,s)
            } }
            fa[p,q] += j-ka; fb[p,q] += j-kb
        } }
        return (fa,fb)
    }
    static func error(_ f: VivoQMMatrix, _ d: VivoQMMatrix, _ s: VivoQMMatrix, _ x: VivoQMMatrix) throws -> [Double] {
        let fds=try f.multiplied(by:d).multiplied(by:s)
        let sdf=try s.multiplied(by:d).multiplied(by:f)
        return try fds.adding(sdf,scale:-1).congruence(x).values
    }
    static func energy(_ ao: VivoAOIntegrals, _ da: VivoQMMatrix, _ db: VivoQMMatrix,
                       _ fa: VivoQMMatrix, _ fb: VivoQMMatrix) -> Double {
        let h=ao.coreHamiltonian
        return ao.constantEnergyHartree + 0.5 * h.values.indices.reduce(0.0) { total,i in
            total + (da.values[i]+db.values[i])*h.values[i]+da.values[i]*fa.values[i]+db.values[i]*fb.values[i]
        }
    }
    public static func solve(system: VivoElectronicSystem, integrals ao: VivoAOIntegrals,
                             configuration cfg: VivoSCFConfiguration = .init(),
                             budget: VivoChemistryBudget = .init()) throws -> VivoHartreeFockResult {
        try system.validate(); try ao.validate(budget:budget); try cfg.validate()
        let n=ao.count,na=system.alphaElectrons,nb=system.betaElectrons
        guard system==ao.sourceSystem,na<=n,nb<=n, cfg.reference != .restricted || na==nb else { throw VivoChemistryError.invalid("HF occupations or restricted spin sector") }
        _ = try budget.elements([n,n],simultaneousArrays:30+4*cfg.diisHistory)
        let overlapEigen = try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        guard let smallest=overlapEigen.values.first,smallest>=cfg.minimumOverlapEigenvalue else {
            throw VivoChemistryError.invalid("linearly dependent AO basis; no silent rank deletion")
        }
        var x=overlapEigen.vectors
        for i in 0..<n { for j in 0..<n { x[i,j] /= sqrt(overlapEigen.values[j]) } }
        func canonical(_ f: VivoQMMatrix) throws -> (VivoQMMatrix,[Double]) {
            let e=try VivoQMDenseAlgebra.symmetricEigen(f.congruence(x))
            return (try x.multiplied(by:e.vectors),e.values)
        }
        let (c0,_) = try canonical(ao.coreHamiltonian)
        var da=density(c0,occupied:na),db=density(c0,occupied:nb)
        var previousEnergy:Double?,densityChange=0.0,history:[VivoDIISRecord]=[],trace:[VivoSCFIteration]=[]
        for iteration in 1...cfg.maximumIterations {
            let (fa,fb)=focks(ao,da,db)
            let e=energy(ao,da,db,fa,fb)
            let ea=try error(fa,da,ao.overlap,x),eb=try error(fb,db,ao.overlap,x),combined=ea+eb
            let residual=sqrt(combined.reduce(0){$0+$1*$1})
            guard e.isFinite,residual.isFinite else { throw VivoChemistryError.convergence("nonfinite HF state") }
            trace.append(.init(iteration:iteration,energyHartree:e,commutatorNorm:residual,densityChange:densityChange))
            if let last=previousEnergy,abs(e-last)<=cfg.energyToleranceHartree,
               densityChange<=cfg.densityTolerance,residual<=cfg.commutatorTolerance {
                let (ca,epsa)=try canonical(fa),(cb,epsb)=try canonical(fb)
                let finalDA=density(ca,occupied:na),finalDB=density(cb,occupied:nb)
                let change=try finalDA.adding(da,scale:-1).frobeniusNorm + finalDB.adding(db,scale:-1).frobeniusNorm
                let (finalFA,finalFB)=focks(ao,finalDA,finalDB)
                let finalErrors=try error(finalFA,finalDA,ao.overlap,x)+error(finalFB,finalDB,ao.overlap,x)
                let finalResidual=sqrt(finalErrors.reduce(0){$0+$1*$1})
                let finalEnergy=energy(ao,finalDA,finalDB,finalFA,finalFB)
                if change<=cfg.densityTolerance,finalResidual<=cfg.commutatorTolerance,
                   abs(finalEnergy-e)<=cfg.energyToleranceHartree {
                    let ab=try ca.transposed.multiplied(by:ao.overlap).multiplied(by:cb)
                    var occupiedOverlap=0.0
                    for i in 0..<na { for j in 0..<nb { occupiedOverlap += ab[i,j]*ab[i,j] } }
                    let sz=0.5*Double(na-nb)
                    let spinSquared=sz*sz+0.5*Double(na+nb)-occupiedOverlap
                    return .init(reference:cfg.reference,alphaElectrons:na,betaElectrons:nb,energyHartree:finalEnergy,
                                 alphaCoefficients:ca,betaCoefficients:cb,alphaOrbitalEnergies:epsa,betaOrbitalEnergies:epsb,
                                 alphaDensity:finalDA,betaDensity:finalDB,finalCommutatorNorm:finalResidual,
                                 spinSquared:max(0,spinSquared),iterations:trace)
                }
            }
            history.append(.init(alphaFock:fa,betaFock:fb,error:combined))
            if history.count>cfg.diisHistory { history.removeFirst() }
            var guessA=fa,guessB=fb
            if history.count>=2 {
                let m=history.count
                var b=VivoQMMatrix(m+1,m+1)
                let scale=max(1e-30,history.map{$0.error.reduce(0){$0+$1*$1}}.max() ?? 1)
                for i in 0..<m { for j in 0..<m {
                    b[i,j]=zip(history[i].error,history[j].error).reduce(0){$0+$1.0*$1.1}/scale
                }; b[i,m] = -1; b[m,i] = -1 }
                var rhs=[Double](repeating:0,count:m+1); rhs[m] = -1
                if let weights=try? VivoQMDenseAlgebra.solve(b,rhs:rhs),weights.allSatisfy({abs($0)<1e6}) {
                    guessA=VivoQMMatrix(n,n); guessB=guessA
                    for i in 0..<m {
                        guessA=try guessA.adding(history[i].alphaFock,scale:weights[i])
                        guessB=try guessB.adding(history[i].betaFock,scale:weights[i])
                    }
                } else if history.count>2 { history.removeFirst() }
            }
            let (ca,_)=try canonical(guessA),(cb,_)=try canonical(guessB)
            var nextA=density(ca,occupied:na),nextB=density(cb,occupied:nb)
            if iteration<=3 {
                nextA=try nextA.scaled(1-cfg.initialDamping).adding(da,scale:cfg.initialDamping)
                nextB=try nextB.scaled(1-cfg.initialDamping).adding(db,scale:cfg.initialDamping)
            }
            densityChange=try nextA.adding(da,scale:-1).frobeniusNorm+nextB.adding(db,scale:-1).frobeniusNorm
            previousEnergy=e;da=nextA;db=nextB
        }
        throw VivoChemistryError.convergence("HF did not satisfy energy, density and commutator criteria in \(cfg.maximumIterations) iterations")
    }
}

extension VivoHartreeFock {
    public static func validate(result: VivoHartreeFockResult, system: VivoElectronicSystem, integrals ao: VivoAOIntegrals,
                                configuration cfg: VivoSCFConfiguration, budget: VivoChemistryBudget = .init()) throws {
        try system.validate();try cfg.validate();try ao.validate(budget:budget)
        let n=ao.count
        guard system==ao.sourceSystem,result.reference==cfg.reference,
              result.alphaElectrons==system.alphaElectrons,result.betaElectrons==system.betaElectrons,
              result.alphaElectrons<=n,result.betaElectrons<=n,
              result.energyHartree.isFinite,result.spinSquared.isFinite,result.spinSquared>=0,
              result.finalCommutatorNorm.isFinite,result.finalCommutatorNorm>=0,result.finalCommutatorNorm<=cfg.commutatorTolerance,
              !result.iterations.isEmpty,result.iterations.count<=cfg.maximumIterations,
              result.alphaOrbitalEnergies.count==n,result.betaOrbitalEnergies.count==n,
              (result.alphaOrbitalEnergies+result.betaOrbitalEnergies).allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("HF result identity/convergence contract") }
        for matrix in [result.alphaCoefficients,result.betaCoefficients,result.alphaDensity,result.betaDensity] {
            guard matrix.rows==n,matrix.columns==n,matrix.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("HF result matrix shape/values") }
        }
        let identity=try VivoQMMatrix.identity(n)
        for c in [result.alphaCoefficients,result.betaCoefficients] {
            guard try ao.overlap.congruence(c).adding(identity,scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("HF result orbital metric") }
        }
        let da=density(result.alphaCoefficients,occupied:system.alphaElectrons),db=density(result.betaCoefficients,occupied:system.betaElectrons)
        guard try da.adding(result.alphaDensity,scale:-1).frobeniusNorm<=cfg.densityTolerance,
              try db.adding(result.betaDensity,scale:-1).frobeniusNorm<=cfg.densityTolerance else { throw VivoChemistryError.invalid("HF result density/orbital mismatch") }
        let (fa,fb)=focks(ao,da,db),actualEnergy=energy(ao,da,db,fa,fb)
        guard abs(actualEnergy-result.energyHartree)<=10*cfg.energyToleranceHartree else { throw VivoChemistryError.invalid("HF result physical energy mismatch") }
        let eig=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
        guard eig.values.allSatisfy({$0>=cfg.minimumOverlapEigenvalue}) else { throw VivoChemistryError.invalid("HF metric rank") }
        var x=eig.vectors
        for i in 0..<n { for j in 0..<n { x[i,j] /= sqrt(eig.values[j]) } }
        let errors=try error(fa,da,ao.overlap,x)+error(fb,db,ao.overlap,x)
        guard sqrt(errors.reduce(0){$0+$1*$1})<=cfg.commutatorTolerance else { throw VivoChemistryError.invalid("HF physical commutator exceeds tolerance") }
    }
}
