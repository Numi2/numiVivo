import Foundation

/// Versioned additions to the same native many-body solver authority. The
/// original request/result coding remains backward compatible.
public enum VivoAdvancedManyBodyRequest:Codable,Sendable,Equatable {
    case directCI(configuration:VivoDavidsonConfiguration)
    case tensorCCSD(configuration:VivoTensorCCSDConfiguration)
    case multistateCASSCF(partition:VivoActiveSpace,configuration:VivoMultiStateCASSCFConfiguration)
}
public enum VivoAdvancedManyBodyResult:Codable,Sendable,Equatable {
    case directCI(result:VivoCISpectrum)
    case tensorCCSD(result:VivoTensorCCSDResult)
    case multistateCASSCF(result:VivoMultiStateCASSCFResult)
    public var converged:Bool {
        switch self { case .directCI(let r):return !r.roots.isEmpty;case .tensorCCSD(let r):return r.converged;case .multistateCASSCF(let r):return r.converged }
    }
    public var objectiveEnergyHartree:Double {
        switch self { case .directCI(let r):return r.roots.first?.energyHartree ?? .nan;case .tensorCCSD(let r):return r.energyHartree;case .multistateCASSCF(let r):return r.weightedEnergyHartree }
    }
    /// A state-average objective must not be substituted for an individual
    /// electronic state's reaction energy. Multiroot states are selected by the
    /// caller; CC response matrices are not exposed as probability densities.
    public var isStateAverage:Bool { if case .multistateCASSCF = self { return true };return false }
}
public extension VivoManyBodySolver {
    static func solve(_ h:VivoEmbeddedHamiltonian,request:VivoAdvancedManyBodyRequest,
                      budget:VivoChemistryBudget = .init()) throws -> VivoAdvancedManyBodyResult {
        switch request {
        case .directCI(let c):return .directCI(result:try VivoDirectCI.solve(h,configuration:c,budget:budget))
        case .tensorCCSD(let c):return .tensorCCSD(result:try VivoTensorCCSD.solve(h,configuration:c,budget:budget))
        case .multistateCASSCF(let p,let c):return .multistateCASSCF(result:try VivoMultiStateCASSCF.solve(h,partition:p,configuration:c,budget:budget))
        }
    }
    static func validate(_ result:VivoAdvancedManyBodyResult,hamiltonian h:VivoEmbeddedHamiltonian,
                         request:VivoAdvancedManyBodyRequest,budget:VivoChemistryBudget = .init()) throws {
        try h.validate(budget:budget)
        guard result.converged,result.objectiveEnergyHartree.isFinite else { throw VivoChemistryError.convergence("advanced many-body result is not converged") }
        func checkCI(_ ci:VivoCIResult,_ input:VivoEmbeddedHamiltonian) throws {
            try ci.state.validate(budget:budget)
            guard ci.method == .fci,ci.state.orbitalCount==input.orbitalCount,ci.state.alphaElectrons==input.alphaElectrons,
                  ci.state.betaElectrons==input.betaElectrons,ci.eigenResidualNorm.isFinite,ci.eigenResidualNorm>=0 else { throw VivoChemistryError.invalid("multiroot CI sector") }
            let energy=try VivoCIDensityMatrices.compute(ci.state,budget:budget).energy(of:input)
            guard abs(energy-ci.energyHartree)<1e-8 else { throw VivoChemistryError.invalid("multiroot CI RDM energy") }
        }
        switch (request,result) {
        case (.directCI(let cfg),.directCI(let r)):
            guard r.roots.count==cfg.roots,!r.roots.isEmpty else { throw VivoChemistryError.invalid("direct CI root count") }
            for ci in r.roots { try checkCI(ci,h);guard ci.eigenResidualNorm<=1.01*cfg.residualTolerance else { throw VivoChemistryError.invalid("direct CI requested residual") } }
        case (.tensorCCSD(let cfg),.tensorCCSD(let r)):
            let occ=(0..<h.alphaElectrons).map { 2*$0 }+(0..<h.betaElectrons).map { 2*$0+1 }
            let vir=(h.alphaElectrons..<h.orbitalCount).map { 2*$0 }+(h.betaElectrons..<h.orbitalCount).map { 2*$0+1 },o=occ.count,v=vir.count
            guard r.occupiedSpinModes==occ,r.virtualSpinModes==vir,r.singles.rows==o,r.singles.columns==v,r.doubles.count==o*o*v*v,
                  r.singles.values.allSatisfy(\.isFinite),r.doubles.allSatisfy(\.isFinite),r.residual.isFinite,r.residual>=0,
                  r.residual<=cfg.residualTolerance,r.iterations<=cfg.maximumIterations else { throw VivoChemistryError.invalid("tensor CCSD amplitude/sector contract") }
            func g(_ p:Int,_ q:Int,_ a:Int,_ b:Int)->Double {
                (p%2==a%2 && q%2==b%2 ? h.eri(p/2,a/2,q/2,b/2):0)-(p%2==b%2 && q%2==a%2 ? h.eri(p/2,b/2,q/2,a/2):0)
            }
            var energy=h.constantEnergyHartree
            for p in occ { energy+=h.oneElectron[p/2,p/2];for q in occ { energy+=0.5*g(p,q,p,q) } }
            guard abs(energy-r.referenceEnergyHartree)<1e-8 else { throw VivoChemistryError.invalid("tensor CCSD reference energy") }
            for i in 0..<o { for a in 0..<v {
                let f=(occ[i]%2==vir[a]%2 ? h.oneElectron[occ[i]/2,vir[a]/2]:0)+occ.reduce(0.0) { $0+g(occ[i],$1,vir[a],$1) }
                energy+=f*r.singles[i,a]
                for j in 0..<o { for b in 0..<v {
                    energy+=0.25*g(occ[i],occ[j],vir[a],vir[b])*(r.doubles[((i*o+j)*v+a)*v+b]+2*r.singles[i,a]*r.singles[j,b])
                } }
            } }
            guard abs(energy-r.energyHartree)<1e-8 else { throw VivoChemistryError.invalid("tensor CCSD amplitude energy") }
        case (.multistateCASSCF(let partition,let cfg),.multistateCASSCF(let r)):
            guard r.configuration==cfg,r.states.count==cfg.weights.count,r.orbitalGradient.allSatisfy(\.isFinite),
                  r.orbitalGradient.reduce(0.0,{hypot($0,$1)})<=cfg.optimization.gradientTolerance else { throw VivoChemistryError.invalid("multistate CASSCF settings/gradient") }
            try partition.validate(for:h,budget:budget)
            let active=try h.rotated(by:r.orbitalRotation,budget:budget).frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:budget)
            for ci in r.states { try checkCI(ci,active);guard ci.eigenResidualNorm<=1.01*cfg.davidson.residualTolerance else { throw VivoChemistryError.invalid("state-average CI residual") } }
            let objective=zip(cfg.weights,r.states).reduce(0.0) { $0+$1.0*$1.1.energyHartree }
            guard abs(objective-r.weightedEnergyHartree)<1e-8 else { throw VivoChemistryError.invalid("state-average objective binding") }
        default:throw VivoChemistryError.invalid("advanced result has a different requested method")
        }
    }
}
