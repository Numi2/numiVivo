import Foundation

/// One explicitly selected classical solver at the portable Hamiltonian boundary.
/// No automatic fallback changes the requested physical/numerical method.
public enum VivoManyBodySolverRequest: Codable, Sendable, Equatable {
    case mp2(minimumGapHartree: Double)
    case configurationInteraction(method: VivoCIMethod)
    case ccsd(configuration: VivoCCSDConfiguration)
    case casci(partition: VivoActiveSpace)
    case casscf(partition: VivoActiveSpace, configuration: VivoCASSCFConfiguration)
}
public struct VivoCASCIResult: Codable, Sendable, Equatable {
    public let partition: VivoActiveSpace
    public let activeHamiltonian: VivoEmbeddedHamiltonian
    public let ci: VivoCIResult
}
public enum VivoManyBodySolverResult: Codable, Sendable, Equatable {
    case mp2(result: VivoMP2Result)
    case configurationInteraction(result: VivoCIResult)
    case ccsd(result: VivoCCSDResult)
    case casci(result: VivoCASCIResult)
    case casscf(result: VivoCASSCFResult)
    public var energyHartree: Double {
        switch self {
        case .mp2(let r): return r.totalEnergyHartree
        case .configurationInteraction(let r): return r.energyHartree
        case .ccsd(let r): return r.energyHartree
        case .casci(let r): return r.ci.energyHartree
        case .casscf(let r): return r.activeCI.energyHartree
        }
    }
    public var converged: Bool {
        switch self {
        case .ccsd(let r): return r.converged
        case .casscf(let r): return r.converged
        default: return true // Other engines throw instead of publishing a failed solve.
        }
    }
    /// Available only when the solver supplies a normalized variational CI state.
    /// CC response RDMs are not density operators for entropy/QIO and are never
    /// silently converted into one through this common interface.
    public var informationState: VivoCIState? {
        switch self {
        case .configurationInteraction(let r): return r.state
        case .casci(let r): return r.ci.state
        case .casscf(let r): return r.activeCI.state
        default: return nil
        }
    }
}
public enum VivoManyBodySolver {
    public static func solve(_ h: VivoEmbeddedHamiltonian, request: VivoManyBodySolverRequest,
                             budget: VivoChemistryBudget = .init()) throws -> VivoManyBodySolverResult {
        try h.validate(budget:budget)
        switch request {
        case .mp2(let gap): return .mp2(result: try VivoRestrictedMP2.solve(h,minimumGapHartree:gap,budget:budget))
        case .configurationInteraction(let method): return .configurationInteraction(result: try VivoConfigurationInteraction.solve(h,method:method,budget:budget))
        case .ccsd(let configuration): return .ccsd(result: try VivoCoupledCluster.solve(h,configuration:configuration,budget:budget))
        case .casci(let partition):
            try partition.validate(for:h,budget:budget)
            let active=try h.frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:budget)
            return .casci(result: .init(partition:partition,activeHamiltonian:active,ci:try VivoConfigurationInteraction.solve(active,budget:budget)))
        case .casscf(let partition,let configuration):
            return .casscf(result: try VivoCASSCF.solve(h,partition:partition,configuration:configuration,budget:budget))
        }
    }
    /// Cache/output contract checks, not external chemical validation. A cache
    /// entry also requires the workflow's input, executable and content digests.
    public static func validate(_ result: VivoManyBodySolverResult, hamiltonian h: VivoEmbeddedHamiltonian,
                                request: VivoManyBodySolverRequest, budget: VivoChemistryBudget = .init()) throws {
        try h.validate(budget:budget)
        guard result.converged, result.energyHartree.isFinite else { throw VivoChemistryError.convergence("requested many-body solver did not converge") }
        func validateCI(_ ci:VivoCIResult,_ input:VivoEmbeddedHamiltonian) throws {
            try ci.state.validate(budget:budget)
            guard ci.state.orbitalCount==input.orbitalCount, ci.state.alphaElectrons==input.alphaElectrons,
                  ci.state.betaElectrons==input.betaElectrons, ci.eigenResidualNorm.isFinite,
                  ci.eigenResidualNorm>=0, ci.eigenResidualNorm<1e-8 else { throw VivoChemistryError.invalid("CI output sector/residual") }
            let e=try VivoCIDensityMatrices.compute(ci.state,budget:budget).energy(of:input)
            guard abs(e-ci.energyHartree)<1e-8 else { throw VivoChemistryError.invalid("CI output RDM energy") }
        }
        switch (request,result) {
        case (.mp2(let gap),.mp2(let r)):
            let expected=try VivoRestrictedMP2.solve(h,minimumGapHartree:gap,budget:budget)
            guard abs(r.totalEnergyHartree-expected.totalEnergyHartree)<1e-9,
                  abs(r.referenceEnergyHartree-expected.referenceEnergyHartree)<1e-9 else { throw VivoChemistryError.invalid("MP2 input/result binding") }
        case (.configurationInteraction(let method),.configurationInteraction(let ci)):
            guard ci.method==method else { throw VivoChemistryError.invalid("CI method binding") }
            try validateCI(ci,h)
        case (.ccsd(let configuration),.ccsd(let cc)):
            try configuration.validate(); try cc.normalizedRightState.validate(budget:budget)
            guard cc.configuration==configuration, cc.orbitalIdentifiers==h.orbitalIdentifiers,
                  cc.energyReference==h.energyReference, cc.normalizedRightState.orbitalCount==h.orbitalCount,
                  cc.normalizedRightState.alphaElectrons==h.alphaElectrons, cc.normalizedRightState.betaElectrons==h.betaElectrons,
                  cc.amplitudes.count==cc.excitations.count, cc.amplitudes.allSatisfy(\.isFinite),
                  cc.projectedResidualNorm.isFinite, cc.projectedResidualNorm>=0,
                  cc.projectedResidualNorm<=configuration.residualTolerance,
                  cc.referenceEnergyHartree.isFinite, cc.rightExpectationEnergyHartree.isFinite else {
                throw VivoChemistryError.invalid("CCSD input/result contract")
            }
            if configuration.solveLambda {
                let e=try VivoCoupledCluster.responseDensityMatrices(cc,budget:budget).energy(of:h)
                guard abs(e-cc.energyHartree)<1e-8 else { throw VivoChemistryError.invalid("CCSD lambda-response energy") }
            }
        case (.casci(let partition),.casci(let r)):
            try partition.validate(for:h,budget:budget)
            let expected=try h.frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:budget)
            guard r.partition==partition, r.activeHamiltonian==expected, r.ci.method == .fci else { throw VivoChemistryError.invalid("CASCI partition/Hamiltonian binding") }
            try validateCI(r.ci,expected)
        case (.casscf(let partition,let configuration),.casscf(let r)):
            try configuration.validate(); try partition.validate(for:h,budget:budget)
            let expected=try h.rotated(by:r.orbitalRotation,budget:budget).frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:budget)
            guard r.partition==partition, r.activeHamiltonian==expected, r.activeCI.method == .fci,
                  r.orbitalGradient.allSatisfy(\.isFinite), r.orbitalGradient.count==r.rotationPairs.count,
                  r.orbitalGradient.reduce(0.0,{hypot($0,$1)})<=configuration.gradientTolerance else {
                throw VivoChemistryError.invalid("CASSCF partition/Hamiltonian/gradient binding")
            }
            try validateCI(r.activeCI,expected)
        default: throw VivoChemistryError.invalid("many-body result has a different method than the request")
        }
    }
}
