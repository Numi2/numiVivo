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
              precomplex.energyHartree.isFinite,transitionState.energyHartree.isFinite,deltaHartree.isFinite,
              abs(deltaHartree-(transitionState.energyHartree-precomplex.energyHartree))<1e-12 else {
            throw VivoChemistryError.invalid("barrier mixes contexts/references or lacks characterized stationary points")
        }
    }
}
public struct VivoUnimolecularRateEstimate: Codable, Sendable, Equatable {
    public let temperatureK: Double
    public let activationGibbsKJPerMol: Double
    public let logRatePerSecond: Double
    public let ratePerSecond: Double
    public let transmissionCoefficient: Double
    public let underflowed: Bool
    public let sourceBarrier: VivoReactionBarrier
    public let status: String
}
public enum VivoReactionEnergies {
    public static func barrier(precomplex: VivoReactionEnergyPoint, transitionState: VivoReactionEnergyPoint) throws -> VivoReactionBarrier {
        let result=VivoReactionBarrier(precomplex:precomplex,transitionState:transitionState,
                                       deltaHartree:transitionState.energyHartree-precomplex.energyHartree)
        try result.validate();return result
    }
    public static func unimolecularRate(_ barrier: VivoReactionBarrier, transmissionCoefficient: Double = 1) throws -> VivoUnimolecularRateEstimate {
        try barrier.validate()
        guard barrier.precomplex.context.energyKind == .gibbs,let temperature=barrier.precomplex.context.temperatureK,
              transmissionCoefficient.isFinite,transmissionCoefficient>0,transmissionCoefficient<=1 else {
            throw VivoChemistryError.invalid("rate conversion requires Gibbs barrier and explicit transmission coefficient in (0,1]")
        }
        let logRate=try VivoUnimolecularEyring.logRatePerSecond(activationGibbsHartree:barrier.deltaHartree,
            temperatureK:temperature,transmissionCoefficient:transmissionCoefficient)
        guard logRate.isFinite,logRate<log(Double.greatestFiniteMagnitude) else { throw VivoChemistryError.invalid("rate overflow") }
        let value=exp(logRate)
        return .init(temperatureK:temperature,activationGibbsKJPerMol:barrier.deltaKJPerMol,logRatePerSecond:logRate,
                     ratePerSecond:value,transmissionCoefficient:transmissionCoefficient,underflowed:value==0,
                     sourceBarrier:barrier,status:"derived-rate-proposal; not automatically installed in a NumiVivo kinetic model")
    }
}
