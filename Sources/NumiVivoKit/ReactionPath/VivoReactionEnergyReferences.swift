import Foundation

public enum VivoReactionEnergyKind: String, Codable, Sendable { case electronic, electronicPlusZeroPoint, gibbs }
public enum VivoReactionPointRole: String, Codable, Sendable { case precomplex, transitionState, product, scanPoint, isolatedReactants }
public struct VivoReactionEnergyContext: Codable, Sendable, Equatable {
    public let method: String
    public let basisRecipeFingerprint: VivoFingerprint
    public let molecularIdentityFingerprint: VivoFingerprint
    public let frozenEnvironmentFingerprint: VivoFingerprint
    public let pathPartitionFingerprint: VivoFingerprint
    public let energyKind: VivoReactionEnergyKind
    public let energyReference: String
    public let alphaElectrons: Int
    public let betaElectrons: Int
    public let temperatureK: Double?
    public let standardState: String?
    public init(method: String, basisRecipeFingerprint: VivoFingerprint, molecularIdentityFingerprint: VivoFingerprint,
                frozenEnvironmentFingerprint: VivoFingerprint, pathPartitionFingerprint: VivoFingerprint,
                energyKind: VivoReactionEnergyKind, energyReference: String, alphaElectrons: Int, betaElectrons: Int,
                temperatureK: Double? = nil, standardState: String? = nil) {
        self.method=method;self.basisRecipeFingerprint=basisRecipeFingerprint;self.molecularIdentityFingerprint=molecularIdentityFingerprint
        self.frozenEnvironmentFingerprint=frozenEnvironmentFingerprint;self.pathPartitionFingerprint=pathPartitionFingerprint
        self.energyKind=energyKind;self.energyReference=energyReference;self.alphaElectrons=alphaElectrons;self.betaElectrons=betaElectrons
        self.temperatureK=temperatureK;self.standardState=standardState
    }
    public func validate() throws {
        guard !method.isEmpty,!energyReference.isEmpty,alphaElectrons>=0,betaElectrons>=0 else { throw VivoChemistryError.invalid("reaction energy context") }
        if energyKind == .gibbs {
            guard let temperatureK,temperatureK.isFinite,temperatureK>0,let standardState,!standardState.isEmpty else {
                throw VivoChemistryError.invalid("Gibbs energy requires temperature and standard-state convention")
            }
        } else if temperatureK != nil || standardState != nil { throw VivoChemistryError.invalid("thermal metadata attached to a non-Gibbs energy") }
    }
}
public struct VivoReactionEnergyPoint: Codable, Sendable, Equatable {
    public let identifier: String
    public let geometryFingerprint: VivoFingerprint
    public let context: VivoReactionEnergyContext
    public let role: VivoReactionPointRole
    public let energyHartree: Double
    public let imaginaryModeCount: Int?
    public let stationaryPointEvidence: VivoFingerprint?
    public init(identifier: String, geometryFingerprint: VivoFingerprint, context: VivoReactionEnergyContext,
                role: VivoReactionPointRole, energyHartree: Double, imaginaryModeCount: Int? = nil,
                stationaryPointEvidence: VivoFingerprint? = nil) {
        self.identifier=identifier;self.geometryFingerprint=geometryFingerprint;self.context=context;self.role=role
        self.energyHartree=energyHartree;self.imaginaryModeCount=imaginaryModeCount;self.stationaryPointEvidence=stationaryPointEvidence
    }
}
public struct VivoReactionBarrier: Codable, Sendable, Equatable {
    public let precomplex: VivoReactionEnergyPoint
    public let transitionState: VivoReactionEnergyPoint
    public let deltaHartree: Double
    public var deltaKJPerMol: Double { deltaHartree*VivoAtomicUnits.hartreeInKJPerMol }
    public func validate() throws {
        try precomplex.context.validate();try transitionState.context.validate()
        guard precomplex.context==transitionState.context,precomplex.role == .precomplex,
              transitionState.role == .transitionState,precomplex.imaginaryModeCount==0,transitionState.imaginaryModeCount==1,
              precomplex.stationaryPointEvidence != nil,transitionState.stationaryPointEvidence != nil,
              !precomplex.identifier.isEmpty,!transitionState.identifier.isEmpty,precomplex.identifier != transitionState.identifier,
              precomplex.energyHartree.isFinite,transitionState.energyHartree.isFinite,deltaHartree.isFinite,deltaHartree>=0,
              abs(deltaHartree-(transitionState.energyHartree-precomplex.energyHartree))<1e-12 else {
            throw VivoChemistryError.invalid("barrier mixes contexts/references, is negative, or lacks characterized stationary points")
        }
    }
}
public enum VivoReactionEnergies {
    public static let boundComplexStandardState = "pre-reactive-bound-complex"

    public static func barrier(precomplex: VivoReactionEnergyPoint, transitionState: VivoReactionEnergyPoint) throws -> VivoReactionBarrier {
        let result=VivoReactionBarrier(precomplex:precomplex,transitionState:transitionState,
                                       deltaHartree:transitionState.energyHartree-precomplex.energyHartree)
        try result.validate();return result
    }

    /// Adapts a characterized NumiVivo Gibbs barrier into the repository's single
    /// kinetic barrier authority. Electronic barriers cannot enter kinetics here.
    public static func activationBarrier(_ barrier: VivoReactionBarrier, kineticContext: VivoKineticContext,
                                         evidence: VivoKineticEvidence,
                                         conditionalStandardDeviationKJPerMol: Double? = nil) throws -> VivoActivationBarrier {
        try barrier.validate();try kineticContext.validate();try evidence.validate(origin:.calculated)
        let context=barrier.precomplex.context
        guard context.energyKind == .gibbs,let temperature=context.temperatureK,
              abs(temperature-kineticContext.temperatureK)<=1e-9*max(1,temperature),
              context.standardState == boundComplexStandardState else {
            throw VivoChemistryError.invalid("kinetics requires a Gibbs barrier at the identical temperature and pre-reactive-bound-complex standard state")
        }
        if let sd=conditionalStandardDeviationKJPerMol {
            guard sd.isFinite,sd>=0 else { throw VivoChemistryError.invalid("barrier conditional standard deviation") }
        }
        return .init(context:kineticContext,quantity:.activationGibbsFreeEnergy,referenceState:.preReactiveBoundComplex,
                     value:barrier.deltaKJPerMol,unit:.kilojoulesPerMol,
                     conditionalStandardDeviation:conditionalStandardDeviationKJPerMol,
                     method:context.method,
                     samplingDescription:"NumiVivo characterized precomplex-to-transition-state Gibbs barrier; source geometries \(barrier.precomplex.geometryFingerprint.hex), \(barrier.transitionState.geometryFingerprint.hex)",
                     origin:.calculated,evidence:evidence)
    }

    public static func unimolecularRate(_ barrier: VivoReactionBarrier, kineticContext: VivoKineticContext,
                                        barrierEvidence: VivoKineticEvidence,
                                        conditionalBarrierStandardDeviationKJPerMol: Double? = nil,
                                        transmissionProbability: Double = 1,
                                        transmissionOrigin: VivoKineticOrigin,
                                        transmissionEvidence: VivoKineticEvidence) throws -> VivoTransitionStateRateEstimate {
        let activation=try activationBarrier(barrier,kineticContext:kineticContext,evidence:barrierEvidence,
                                             conditionalStandardDeviationKJPerMol:conditionalBarrierStandardDeviationKJPerMol)
        let request=VivoTransitionStateRateRequest(barrier:activation,transmissionProbability:transmissionProbability,
                                                   transmissionOrigin:transmissionOrigin,transmissionEvidence:transmissionEvidence)
        return try VivoTransitionStateRateEstimator.estimate(request)
    }
}
