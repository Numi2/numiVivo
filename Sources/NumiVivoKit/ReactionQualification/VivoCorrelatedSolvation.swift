import Foundation

public struct VivoCorrelatedSolventConfiguration: Codable, Sendable, Equatable {
    /// Nil: full CI. An explicit partition: fixed-orbital CASCI with occupied
    /// core and empty external orbitals. No overlapping-fragment density sum.
    public var partition: VivoActiveSpace?
    public var maximumIterations: Int
    public var densityTolerance: Double
    public var potentialToleranceHartree: Double
    public var energyToleranceHartree: Double
    public var ciResidualTolerance: Double
    public var damping: Double
    public init(partition: VivoActiveSpace? = nil, maximumIterations: Int = 128,
                densityTolerance: Double = 1e-8, potentialToleranceHartree: Double = 1e-8,
                energyToleranceHartree: Double = 1e-10, ciResidualTolerance: Double = 1e-11,
                damping: Double = 0) {
        self.partition=partition; self.maximumIterations=maximumIterations; self.densityTolerance=densityTolerance
        self.potentialToleranceHartree=potentialToleranceHartree; self.energyToleranceHartree=energyToleranceHartree
        self.ciResidualTolerance=ciResidualTolerance; self.damping=damping
    }
    public func validate() throws {
        guard (2...4096).contains(maximumIterations), damping.isFinite, damping>=0, damping<1,
              [densityTolerance,potentialToleranceHartree,energyToleranceHartree,ciResidualTolerance]
                .allSatisfy({$0.isFinite && $0>0}), ciResidualTolerance<=energyToleranceHartree else {
            throw VivoChemistryError.invalid("correlated solvent iterations, tolerances or damping")
        }
    }
}
public struct VivoCorrelatedSolventRequest: Codable, Sendable, Equatable {
    public let system: VivoElectronicSystem
    public let basis: VivoGaussianBasis
    public let solvent: VivoSmoothCPCMConfiguration
    public let configuration: VivoCorrelatedSolventConfiguration
    /// Complete fixed orthonormal AO coefficient frame. Nil uses canonical
    /// overlap orthogonalization; it does not infer an active-space partition.
    public let coefficients: VivoQMMatrix?
    public let budget: VivoChemistryBudget
    public func validate() throws {
        try configuration.validate();try system.validate();try basis.validate(nucleusCount:system.nuclei.count);try budget.validate()
        guard coefficients?.values.allSatisfy(\.isFinite) ?? true else {throw VivoChemistryError.invalid("nonfinite correlated solvent frame")}
    }
    public init(system: VivoElectronicSystem, basis: VivoGaussianBasis,
                solvent: VivoSmoothCPCMConfiguration = .init(), configuration: VivoCorrelatedSolventConfiguration = .init(),
                coefficients: VivoQMMatrix? = nil, budget: VivoChemistryBudget = .init()) {
        self.system=system; self.basis=basis; self.solvent=solvent; self.configuration=configuration
        self.coefficients=coefficients; self.budget=budget
    }
}
public struct VivoCorrelatedSolventIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let densityResidual: Double
    public let potentialResidualHartree: Double
    public let energyChangeHartree: Double?
}
public struct VivoCorrelatedSolventResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/correlated-equilibrium-cpcm/v1"
    public let schema: String
    public let request: VivoCorrelatedSolventRequest
    public let coefficients: VivoQMMatrix
    public let partition: VivoActiveSpace
    public let ci: VivoCIResult
    public let totalDensityAO: VivoQMMatrix
    public let inputDensityAO: VivoQMMatrix
    public let equilibriumField: VivoSmoothCPCMResult
    public let gasEnergyHartree: Double
    public let energyHartree: Double
    public let densityResidual: Double
    public let potentialResidualHartree: Double
    public let history: [VivoCorrelatedSolventIteration]
    public let method: String
}

/// Alternating CI wavefunction/surface-field equilibrium at fixed nuclei.
/// FP64 authority. A stationary solution is not a global-minimum guarantee.
/// Response CC matrices and overlapping-fragment densities are not accepted as
/// probability operators. Electronic plus electrostatic, not Gibbs energy.
public enum VivoCorrelatedSolvation {
    public static let method="variational fixed-orbital FCI/CASCI + equilibrium smooth C-PCM; correlated density feedback; electronic plus electrostatic polarization, no thermal/non-electrostatic correction"
    private struct Workspace {
        let request: VivoCorrelatedSolventRequest
        let ao: VivoAOIntegrals
        let c: VivoQMMatrix
        let gas: VivoEmbeddedHamiltonian
        let partition: VivoActiveSpace
        let pcm: VivoSmoothCPCMOperator
        init(_ r: VivoCorrelatedSolventRequest) throws {
            try r.validate()
            let ao=try VivoGaussianIntegralEngine.compute(system:r.system,basis:r.basis,budget:r.budget),n=ao.count
            let c: VivoQMMatrix
            if let specified=r.coefficients { c=specified }
            else {
                let eig=try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
                guard eig.values.first!>=1e-8 else {throw VivoChemistryError.invalid("correlated PCM AO linear dependence")}
                var orthogonal=eig.vectors
                for i in 0..<n {for j in 0..<n {orthogonal[i,j]/=sqrt(eig.values[j])}}
                c=orthogonal
            }
            guard c.rows==n,c.columns==n,c.values.allSatisfy(\.isFinite),
                  try ao.overlap.congruence(c).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {
                throw VivoChemistryError.invalid("correlated PCM requires a complete orthonormal spatial frame")
            }
            let gas=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:c,alphaElectrons:r.system.alphaElectrons,
                betaElectrons:r.system.betaElectrons,orbitalIdentifiers:(0..<n).map{"equilibrium-orbital-\($0)"},
                energyReference:"physical electronic energy; AO scalar once; equilibrium solvent polarization accounted separately",budget:r.budget)
            let partition=r.configuration.partition ?? .init(active:Array(0..<n))
            try partition.validate(for:gas,budget:r.budget)
            self.request=r; self.ao=ao; self.c=c; self.gas=gas; self.partition=partition
            pcm=try .init(system:r.system,basis:r.basis,configuration:r.solvent,budget:r.budget)
        }
        func solveCI(field: VivoQMMatrix) throws -> VivoCIResult {
            let effective=VivoEmbeddedHamiltonian(orbitalIdentifiers:gas.orbitalIdentifiers,alphaElectrons:gas.alphaElectrons,
                betaElectrons:gas.betaElectrons,oneElectron:try gas.oneElectron.adding(field),twoElectron:gas.twoElectron,
                constantEnergyHartree:gas.constantEnergyHartree,energyReference:gas.energyReference)
            let active=try effective.frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:request.budget)
            return try VivoDirectCI.solve(active,configuration:.init(residualTolerance:request.configuration.ciResidualTolerance),budget:request.budget).roots[0]
        }
        func spatialDensity(_ state: VivoCIState) throws -> VivoQMMatrix {
            try state.validate(budget:request.budget)
            let n=gas.orbitalCount,a=partition.active.count
            guard state.orbitalCount==a,state.alphaElectrons==gas.alphaElectrons-partition.doublyOccupiedCore.count,
                  state.betaElectrons==gas.betaElectrons-partition.doublyOccupiedCore.count else {throw VivoChemistryError.invalid("solvent CI state electron sector")}
            let (work,overflow)=(2*a*a).multipliedReportingOverflow(by:state.determinants.count)
            guard !overflow,work<=request.budget.maximumOperatorApplications else {throw VivoChemistryError.resourceLimit("correlated solvent one-RDM work")}
            let index=Dictionary(uniqueKeysWithValues:state.determinants.enumerated().map{($0.element,$0.offset)})
            var p=VivoQMMatrix(n,n)
            for i in partition.doublyOccupiedCore {p[i,i]=2}
            for i in 0..<a {for j in 0..<a {for spin in 0..<2 {
                p[partition.active[i],partition.active[j]] += VivoCIDensityMatrices.expectation(state,index,
                    [.init(mode:2*j+spin,creation:false),.init(mode:2*i+spin,creation:true)])
            }}}
            return p
        }
        func aoDensity(_ p: VivoQMMatrix) throws -> VivoQMMatrix {try c.multiplied(by:p).multiplied(by:c.transposed)}
    }
    private static func contraction(_ a: VivoQMMatrix,_ b: VivoQMMatrix) -> Double {
        zip(a.values,b.values).reduce(0.0){$0+$1.0*$1.1}
    }
    public static func solve(_ request: VivoCorrelatedSolventRequest) throws -> VivoCorrelatedSolventResult {
        let w=try Workspace(request),cfg=request.configuration
        var input=try w.spatialDensity(w.solveCI(field:VivoQMMatrix(w.gas.orbitalCount,w.gas.orbitalCount)).state)
        var previous:Double?,history:[VivoCorrelatedSolventIteration]=[]
        for iteration in 1...cfg.maximumIterations {
            let inputAO=try w.aoDensity(input),oldField=try w.pcm.evaluate(totalDensity:inputAO)
            let potential=try oldField.reactionPotentialMatrix.congruence(w.c)
            let ci=try w.solveCI(field:potential),p=try w.spatialDensity(ci.state),density=try w.aoDensity(p)
            let field=try w.pcm.evaluate(totalDensity:density)
            let densityResidual=try p.adding(input,scale:-1).frobeniusNorm
            let potentialResidual=try field.reactionPotentialMatrix.adding(oldField.reactionPotentialMatrix,scale:-1).congruence(w.c).frobeniusNorm
            let gas=ci.energyHartree-contraction(p,potential),energy=gas+field.polarizationEnergyHartree
            let change=previous.map{abs(energy-$0)}
            guard energy.isFinite,densityResidual.isFinite,potentialResidual.isFinite else {throw VivoChemistryError.convergence("correlated equilibrium PCM nonfinite iterate")}
            history.append(.init(iteration:iteration,energyHartree:energy,densityResidual:densityResidual,potentialResidualHartree:potentialResidual,energyChangeHartree:change))
            if let change,change<=cfg.energyToleranceHartree,densityResidual<=cfg.densityTolerance,potentialResidual<=cfg.potentialToleranceHartree {
                return .init(schema:VivoCorrelatedSolventResult.schema,request:request,coefficients:w.c,partition:w.partition,ci:ci,
                    totalDensityAO:density,inputDensityAO:inputAO,equilibriumField:field,gasEnergyHartree:gas,energyHartree:energy,
                    densityResidual:densityResidual,potentialResidualHartree:potentialResidual,history:history,method:method)
            }
            previous=energy
            input=try p.scaled(1-cfg.damping).adding(input,scale:cfg.damping)
        }
        throw VivoChemistryError.convergence("correlated-density smooth C-PCM did not reach joint CI/field stationarity")
    }
    public static func validate(_ result: VivoCorrelatedSolventResult,request: VivoCorrelatedSolventRequest) throws {
        let w=try Workspace(request),cfg=request.configuration
        guard result.schema==VivoCorrelatedSolventResult.schema,result.request==request,result.method==method,result.ci.method == .fci,
              result.partition==w.partition,result.coefficients.rows==w.c.rows,result.coefficients.columns==w.c.columns,
              try result.coefficients.adding(w.c,scale:-1).frobeniusNorm<1e-9,
              result.history.count>=2,result.history.count<=cfg.maximumIterations,
              result.densityResidual.isFinite,result.densityResidual>=0,result.densityResidual<=cfg.densityTolerance,
              result.potentialResidualHartree.isFinite,result.potentialResidualHartree>=0,result.potentialResidualHartree<=cfg.potentialToleranceHartree else {
            throw VivoChemistryError.invalid("correlated solvent result identity, frame or convergence contract")
        }
        let p=try w.spatialDensity(result.ci.state),d=try w.aoDensity(p)
        guard try d.adding(result.totalDensityAO,scale:-1).frobeniusNorm<1e-9 else {throw VivoChemistryError.invalid("correlated solvent density does not belong to the stored CI state")}
        let field=try w.pcm.evaluate(totalDensity:d),inputField=try w.pcm.evaluate(totalDensity:result.inputDensityAO)
        let v=try inputField.reactionPotentialMatrix.congruence(w.c),actual=try w.solveCI(field:v)
        let fciEnergy=actual.energyHartree
        let active=try VivoEmbeddedHamiltonian(orbitalIdentifiers:w.gas.orbitalIdentifiers,alphaElectrons:w.gas.alphaElectrons,
            betaElectrons:w.gas.betaElectrons,oneElectron:try w.gas.oneElectron.adding(v),twoElectron:w.gas.twoElectron,
            constantEnergyHartree:w.gas.constantEnergyHartree,energyReference:w.gas.energyReference)
            .frozenCore(active:w.partition.active,doublyOccupiedCore:w.partition.doublyOccupiedCore,budget:request.budget)
        let stateEnergy=try VivoCIDensityMatrices.compute(result.ci.state,budget:request.budget).energy(of:active)
        let residual=try VivoDirectCI.residualNorm(hamiltonian:active,state:result.ci.state,energyHartree:stateEnergy,budget:request.budget)
        let inputP=try w.c.transposed.multiplied(by:w.ao.overlap).multiplied(by:result.inputDensityAO)
            .multiplied(by:w.ao.overlap).multiplied(by:w.c)
        let densityResidual=try inputP.adding(p,scale:-1).frobeniusNorm
        let potentialResidual=try field.reactionPotentialMatrix.adding(inputField.reactionPotentialMatrix,scale:-1).congruence(w.c).frobeniusNorm
        guard residual<=max(1e-12,2*cfg.ciResidualTolerance),
              densityResidual<=1.01*cfg.densityTolerance,potentialResidual<=1.01*cfg.potentialToleranceHartree,
              abs(densityResidual-result.densityResidual)<1e-10,
              abs(potentialResidual-result.potentialResidualHartree)<1e-10 else {throw VivoChemistryError.invalid("reported correlated solvent residual differs from physical equations")}
        for (i,item) in result.history.enumerated() {
            guard item.iteration==i+1,[item.energyHartree,item.densityResidual,item.potentialResidualHartree].allSatisfy(\.isFinite),
                  item.densityResidual>=0,item.potentialResidualHartree>=0 else {throw VivoChemistryError.invalid("correlated solvent iteration trace")}
            if i==0 {guard item.energyChangeHartree==nil else {throw VivoChemistryError.invalid("first solvent iteration has a fabricated prior energy")}}
            else {guard let change=item.energyChangeHartree,change.isFinite,
                abs(change-abs(item.energyHartree-result.history[i-1].energyHartree))<1e-12 else {throw VivoChemistryError.invalid("solvent energy-change trace")}}
        }
        guard let last=result.history.last,last.energyChangeHartree!<=cfg.energyToleranceHartree,
              abs(last.energyHartree-result.energyHartree)<1e-12,
              last.densityResidual==result.densityResidual,last.potentialResidualHartree==result.potentialResidualHartree else {
            throw VivoChemistryError.invalid("solvent final trace does not match returned state")
        }
        let gas=stateEnergy-contraction(p,v),energy=gas+field.polarizationEnergyHartree
        // Check density feedback at the returned correlated field, not only the
        // earlier field that generated the last CI vector.
        let closure=try w.solveCI(field:field.reactionPotentialMatrix.congruence(w.c))
        let closureDensity=try w.spatialDensity(closure.state)
        guard abs(stateEnergy-fciEnergy)<max(1e-9,10*cfg.ciResidualTolerance),
              abs(stateEnergy-result.ci.energyHartree)<1e-9,
              try closureDensity.adding(p,scale:-1).frobeniusNorm<=2*cfg.densityTolerance,
              abs(gas-result.gasEnergyHartree)<1e-9,abs(energy-result.energyHartree)<1e-9,
              result.equilibriumField.configuration==request.solvent,
              result.equilibriumField.tesseraCount==field.tesseraCount,
              result.equilibriumField.method==field.method,
              result.equilibriumField.trueRelativeResidual.isFinite,
              result.equilibriumField.trueRelativeResidual<=8*request.solvent.tolerance,
              abs(result.equilibriumField.surfaceChargeDefect-field.surfaceChargeDefect)<1e-8,
              result.equilibriumField.apparentSurfaceCharges.count==field.apparentSurfaceCharges.count,
              zip(result.equilibriumField.apparentSurfaceCharges,field.apparentSurfaceCharges).reduce(0.0,{hypot($0,$1.0-$1.1)})<1e-8,
              abs(result.equilibriumField.polarizationEnergyHartree-field.polarizationEnergyHartree)<1e-9,
              try result.equilibriumField.reactionPotentialMatrix.adding(field.reactionPotentialMatrix,scale:-1).frobeniusNorm<1e-9 else {
            throw VivoChemistryError.invalid("correlated equilibrium density, energy or reaction-field reconstruction")
        }
    }
}
